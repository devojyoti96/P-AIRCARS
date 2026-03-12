import logging
import numpy as np
import argparse
import traceback
import copy
import time
import glob
import sys
import os
from casatools import msmetadata
from dask import delayed
from paircars.utils.basic_utils import timestamp_to_mjdsec, mjdsec_to_timestamp
from paircars.utils.image_utils import (
    create_circular_mask,
    make_stokes_wsclean_imagecube,
)
from paircars.utils.imaging import (
    calc_sun_dia,
    calc_field_of_view,
    calc_npix_in_psf,
    calc_cellsize,
    calc_multiscale_scales,
    get_multiscale_bias,
    get_fft_size,
)
from paircars.utils.logger_utils import (
    SmartDefaultsHelpFormatter,
    clean_shutdown,
    create_logger,
    init_logger,
)
from paircars.utils.ms_metadata import check_datacolumn_valid, get_ms_size
from paircars.utils.mwa_utils import get_ncoarse
from paircars.utils.mwa_ploting_utils import rename_mwasolar_image
from paircars.utils.proc_manage_utils import (
    scale_worker_and_wait,
    get_local_dask_cluster,
)
from paircars.utils.resource_utils import drop_cache
from paircars.utils.udocker_utils import (
    check_udocker_container,
    initialize_wsclean_container,
    run_wsclean,
)

logging.getLogger("distributed").setLevel(logging.ERROR)
logging.getLogger("tornado.application").setLevel(logging.CRITICAL)


def perform_imaging(
    msname="",
    workdir="",
    datacolumn="CORRECTED_DATA",
    freqrange="",
    timerange="",
    imagedir="",
    imsize=1024,
    cellsize=2,
    image_freqres=-1,
    image_timeres=-1,
    pol="I",
    weight="briggs",
    robust=0.0,
    minuv=0,
    threshold=1.0,
    use_multiscale=True,
    use_solar_mask=True,
    mask_radius=40,
    savemodel=True,
    saveres=True,
    cutout_rsun=10.0,
    make_plots=True,
    logfile="imaging.log",
    ncpu=1,
    mem=1,
):
    """
    Perform spectropolarimetric snapshot imaging of a ms

    Parameters
    ----------
    msname : str
        Name of the measurement set
    workdir : str
        Work directory name
    datacolumn : str, optional
        Data column
    freqrange : str, optional
        Frequency range to image
    timerange : str, optional
        Time range to image
    imagedir : str, optional
        Image directory name (default: workdir). Images, models, residuals will be saved in directories named images. models, residuals inside imagedir
    imsize : int, optional
        Image size in pixels
    cellsize : float, optional
        Cell size in arcseconds
    image_freqres : float, optional
        Frequency resolution of image in MHz
    image_timeres : float, optional
        Time resolution of image in seconds
    pol : str, optional
        Stokes parameters to image
    weight : str, optional
        Image weighting scheme
    robust : float, optional
        Briggs weighting robustness parameter
    minuv : float, optional
        Minimum UV-lambda to be used in imaging
    threshold : float, optional
        CLEAN threshold
    use_multiscale : bool, optional
        Use multiscale or not
    use_solar_mask : bool, optional
        Use solar mask
    mask_radius : float, optional
        Mask radius in arcminute
    savemodel : bool, optional
        Save model images or not
    saveres : bool, optional
        Save residual images or not
    cutout_rsun : float, optional
        Cutout image size in solar radii from center (default: 10.0 solar radii)
    make_plots : bool, optional
        Make radio map helioprojective plots
    logfile : str, optional
        Log file name
    ncpu : int, optional
        Number of CPU threads to use
    mem : float, optional
        Memory in GB to use

    Returns
    -------
    int
        Success message
    list
        List of images [[images],[models],[residuals]]
    """
    ncpu = max(1, ncpu)
    mem = abs(mem)

    if os.path.exists(logfile):
        os.system(f"rm -rf {logfile}")
    logger, logfile = create_logger(
        os.path.basename(logfile).split(".log")[0],
        logfile,
        verbose=False,
        get_print=True,
    )
    sub_observer = None
    if os.path.exists(f"{workdir}/.jobname_password.npy") and logfile is not None:
        time.sleep(5)
        jobname, password = np.load(
            f"{workdir}/.jobname_password.npy", allow_pickle=True
        )
        if os.path.exists(logfile):
            sub_observer = init_logger(
                "remotelogger_imaging_{os.path.basename(msname).split('.ms')[0]}",
                logfile,
                jobname=jobname,
                password=password,
            )
    try:
        msname = msname.rstrip("/")
        msname = os.path.abspath(msname)
        logger.info(f"{os.path.basename(msname)} --Perform imaging...\n")

        #########
        # Imaging
        #########
        msmd = msmetadata()
        msmd.open(msname)
        msmd.meanfreq(0, unit="MHz")
        freqs = msmd.chanfreqs(0, unit="MHz")
        freqres = msmd.chanres(0, unit="MHz")[0]
        times = msmd.timesforspws(0)
        timeres = np.nanmin(np.diff(times))
        npol = msmd.ncorrforpol()[0]
        msmd.close()

        ###################################################
        # Whether calibration solutions were applied or not
        ###################################################
        if os.path.exists(f"{msname}/.nocal"):
            cal_sol = False
        else:
            cal_sol = True

        ####################################
        # Whether pol-selfcal is done or not
        ####################################
        if os.path.exists(f"{msname}/.nopolselfcal"):
            pol_selfcal = False
        else:
            pol_selfcal = True

        ###################################
        # Finding channel and time ranges
        ###################################
        if freqrange != "":
            start_chans = []
            end_chans = []
            freq_list = freqrange.split(",")
            for f in freq_list:
                start_freq = float(f.split("~")[0])
                end_freq = float(f.split("~")[-1])
                if start_freq >= np.nanmin(freqs) and end_freq <= np.nanmax(freqs):
                    start_chan = np.argmin(np.abs(start_freq - freqs))
                    end_chan = np.argmin(np.abs(end_freq - freqs))
                    start_chans.append(start_chan)
                    end_chans.append(end_chan)
        else:
            start_chans = [0]
            end_chans = [len(freqs)]
        if len(start_chans) == 0:
            print(f"Please provide valid channel range between 0 and {len(freqs)}")
            time.sleep(5)
            clean_shutdown(sub_observer)
            return 1, []

        if timerange != "":
            start_times = []
            end_times = []
            time_list = timerange.split(",")
            for timerange in time_list:
                start_time = timestamp_to_mjdsec(timerange.split("~")[0])
                end_time = timestamp_to_mjdsec(timerange.split("~")[-1])
                if start_time >= np.nanmin(times) and end_time <= np.nanmax(times):
                    start_times.append(np.argmin(np.abs(times - start_time)))
                    end_times.append(np.argmin(np.abs(times - end_time)))
        else:
            start_times = [0]
            end_times = [len(times)]
        if len(start_times) == 0:
            print(
                f"Please provide valid time range between {mjdsec_to_timestamp(times[0])} and {mjdsec_to_timestamp(times[-1])}"
            )
            time.sleep(5)
            clean_shutdown(sub_observer)
            return 1, []

        if npol < 4 and pol == "IQUV":
            pol = "I"
        prefix = workdir + "/imaging_" + os.path.basename(msname).split(".ms")[0]
        if imagedir == "":
            imagedir = workdir
        os.makedirs(imagedir, exist_ok=True)
        if weight == "briggs":
            weight += " " + str(robust)
        if threshold <= 1:
            threshold = 1.1

        wsclean_args = [
            "-quiet",
            "-scale " + str(cellsize) + "asec",
            "-size " + str(imsize) + " " + str(imsize),
            "-no-dirty",
            "-gridder wgridder",
            "-weight " + weight,
            "-name " + prefix,
            "-pol " + str(pol),
            "-niter 10000",
            "-mgain 0.85",
            "-nmiter 5",
            "-gain 0.1",
            "-minuv-l " + str(minuv),
            "-j " + str(ncpu),
            "-abs-mem " + str(round(mem, 2)),
            "-auto-threshold 1 -auto-mask " + str(threshold),
            "-no-update-model-required",
        ]
        if datacolumn != "CORRECTED_DATA" and datacolumn != "corrected":
            wsclean_args.append("-data-column " + datacolumn)

        ngrid = max(1, int(ncpu / 2))
        if ngrid > 1:
            wsclean_args.append("-parallel-gridding " + str(ngrid))

        if pol == "I":
            wsclean_args.append("-no-negative")

        ################################################
        # Creating and using a solar mask
        ################################################
        if use_solar_mask:
            fits_mask = prefix + "_solar-mask.fits"
            if not os.path.exists(fits_mask):
                logger.info(
                    f"{os.path.basename(msname)} -- Creating solar mask of radius: {mask_radius} arcmin.\n",
                )
                fits_mask = create_circular_mask(
                    msname, cellsize, imsize, mask_radius=mask_radius
                )
            if fits_mask is not None and os.path.exists(fits_mask):
                wsclean_args.append("-fits-mask " + fits_mask)
        final_list_dic = {"image": [], "model": [], "residual": []}
        for i in range(len(start_chans)):
            for j in range(len(start_times)):
                temp_wsclean_args = copy.deepcopy(wsclean_args)
                temp_wsclean_args.append(
                    f"-channel-range {start_chans[i]} {end_chans[i]}"
                )
                temp_wsclean_args.append(f"-interval {start_times[j]} {end_times[j]}")

                #####################################
                # Spectral imaging configuration
                #####################################
                if image_freqres > 0:
                    nchan = max(1, int(image_freqres / freqres))
                    chan_chunk = int((end_chans[i] - start_chans[i]) / nchan)
                    if chan_chunk > 1:
                        temp_wsclean_args.append(f"-channels-out {chan_chunk}")
                        temp_wsclean_args.append("-no-mf-weighting")

                #####################################
                # Temporal imaging configuration
                #####################################
                if image_timeres > 0:
                    ntime = max(1, int(image_timeres / timeres))
                    time_chunk = int((end_times[j] - start_times[j]) / ntime)
                    if time_chunk > 1:
                        temp_wsclean_args.append(f"-intervals-out {time_chunk}")

                ######################################
                # Multiscale configuration
                ######################################
                if use_multiscale:
                    num_pixel_in_psf = calc_npix_in_psf(weight, robust=robust)
                    chan_number = int((start_chans[i] + end_chans[i]) / 2)
                    freqMHz = freqs[chan_number]
                    sun_dia = calc_sun_dia(freqMHz)  # Sun diameter in arcmin
                    sun_rad = sun_dia / 2
                    multiscale_scales = calc_multiscale_scales(
                        msname,
                        num_pixel_in_psf,
                        chan_number=chan_number,
                        max_scale=sun_rad,
                    )
                    temp_wsclean_args.append("-multiscale")
                    temp_wsclean_args.append("-multiscale-gain 0.1")
                    temp_wsclean_args.append(
                        "-multiscale-scales "
                        + ",".join([str(s) for s in multiscale_scales])
                    )
                    mid_freq = np.nanmean(
                        freqs[int(start_chans[i]) : int(end_chans[i])]
                    )
                    scale_bias = get_multiscale_bias(mid_freq)
                    temp_wsclean_args.append(f"-multiscale-scale-bias {scale_bias}")
                    if imsize >= 1024 and 4 * max(multiscale_scales) < 512:
                        temp_wsclean_args.append("-parallel-deconvolution 512")
                elif imsize >= 1024:
                    temp_wsclean_args.append("-parallel-deconvolution 512")

                ######################################
                # Running imaging
                ######################################
                wsclean_cmd = "wsclean " + " ".join(temp_wsclean_args) + " " + msname
                logger.info(
                    f"{os.path.basename(msname)} -- WSClean command: {wsclean_cmd}\n",
                )
                msg = run_wsclean(wsclean_cmd, "paircarswsclean", verbose=False)
                if msg == 0:
                    os.system("rm -rf " + prefix + "*psf.fits")
                    ######################
                    # Making stokes cubes
                    ######################
                    pollist = [i.upper() for i in list(pol)]
                    if len(pollist) == 1:
                        imagelist = sorted(glob.glob(prefix + "*image.fits"))
                        if not savemodel:
                            os.system("rm -rf " + prefix + "*model.fits")
                        else:
                            modellist = sorted(glob.glob(prefix + "*model.fits"))
                        if not saveres:
                            os.system("rm -rf " + prefix + "*residual.fits")
                        else:
                            reslist = sorted(glob.glob(prefix + "*residual.fits"))
                    else:
                        imagelist = []
                        stokeslist = []
                        for p in pollist:
                            stokeslist.append(
                                sorted(glob.glob(prefix + "*" + p + "-image.fits"))
                            )
                        for i in range(len(stokeslist[0])):
                            wsclean_images = sorted(
                                [stokeslist[k][i] for k in range(len(pollist))]
                            )
                            image_prefix = os.path.basename(wsclean_images[0]).split(
                                "-image"
                            )[0]
                            image_cube = make_stokes_wsclean_imagecube(
                                wsclean_images,
                                image_prefix + f"_{pol}_image.fits",
                                keep_wsclean_images=False,
                            )
                            imagelist.append(image_cube)
                        del stokeslist
                        if not savemodel:
                            os.system("rm -rf " + prefix + "*model.fits")
                        else:
                            modellist = []
                            stokeslist = []
                            for p in pollist:
                                stokeslist.append(
                                    sorted(glob.glob(prefix + f"*{p}*model.fits"))
                                )
                            for i in range(len(stokeslist[0])):
                                wsclean_models = sorted(
                                    [stokeslist[k][i] for k in range(len(pollist))]
                                )
                                model_prefix = os.path.basename(
                                    wsclean_models[0]
                                ).split("-model")[0]
                                model_cube = make_stokes_wsclean_imagecube(
                                    wsclean_models,
                                    model_prefix + f"_{pol}_model.fits",
                                    keep_wsclean_images=False,
                                )
                                modellist.append(model_cube)
                            del stokeslist
                        if not saveres:
                            os.system("rm -rf " + prefix + "*residual.fits")
                        else:
                            reslist = []
                            stokeslist = []
                            for p in pollist:
                                stokeslist.append(
                                    sorted(glob.glob(prefix + f"*{p}*residual.fits"))
                                )
                            for i in range(len(stokeslist[0])):
                                wsclean_residuals = sorted(
                                    [stokeslist[k][i] for k in range(len(pollist))]
                                )
                                res_prefix = os.path.basename(
                                    wsclean_residuals[0]
                                ).split("-residual")[0]
                                residual_cube = make_stokes_wsclean_imagecube(
                                    wsclean_residuals,
                                    res_prefix + f"_{pol}_residual.fits",
                                    keep_wsclean_images=False,
                                )
                                reslist.append(residual_cube)
                            del stokeslist

                    ######################
                    # Renaming images
                    ######################
                    if len(imagelist) > 0:
                        logger.info(f"Total {len(imagelist)} images are made.")
                        logger.info("Renaming and making plots...")
                        os.makedirs(imagedir + "/images", exist_ok=True)
                        final_image_list = []
                        for imagename in imagelist:
                            renamed_image = rename_mwasolar_image(
                                imagename,
                                imagetype="image",
                                imagedir=imagedir + "/images",
                                pol=pol,
                                cutout_rsun=cutout_rsun,
                                make_plots=make_plots,
                                pol_selfcal=pol_selfcal,
                                cal_sol=cal_sol,
                            )
                            if renamed_image is not None:
                                final_image_list.append(renamed_image)
                        final_list_dic["image"] = final_image_list
                        if savemodel and len(modellist) > 0:
                            final_model_list = []
                            os.makedirs(imagedir + "/models", exist_ok=True)
                            for modelname in modellist:
                                renamed_model = rename_mwasolar_image(
                                    modelname,
                                    imagetype="model",
                                    imagedir=imagedir + "/models",
                                    pol=pol,
                                    cutout_rsun=cutout_rsun,
                                    make_plots=False,
                                    pol_selfcal=pol_selfcal,
                                    cal_sol=cal_sol,
                                )
                                if renamed_model is not None:
                                    final_model_list.append(renamed_model)
                            final_list_dic["model"] = final_model_list
                        if saveres and len(reslist) > 0:
                            final_res_list = []
                            os.makedirs(imagedir + "/residuals", exist_ok=True)
                            for resname in reslist:
                                renamed_res = rename_mwasolar_image(
                                    resname,
                                    imagetype="residual",
                                    imagedir=imagedir + "/residuals",
                                    pol=pol,
                                    cutout_rsun=cutout_rsun,
                                    make_plots=False,
                                    pol_selfcal=pol_selfcal,
                                )
                                if renamed_res is not None:
                                    final_res_list.append(renamed_res)
                            final_list_dic["residual"] = final_res_list
            if os.path.exists(f"{imagedir}/images/dask-scratch-space"):
                os.system(f"rm -rf {imagedir}/images/dask-scratch-space")
            if use_solar_mask and os.path.exists(fits_mask):
                os.system("rm -rf " + fits_mask)
            if len(final_list_dic["image"]) == 0:
                logger.info(
                    f"{os.path.basename(msname)} -- No image is made.\n",
                )
                time.sleep(5)
                clean_shutdown(sub_observer)
                return 1, final_list_dic
            else:
                logger.info(
                    f"{os.path.basename(msname)} -- Imaging is successfully done.\n",
                )
                time.sleep(5)
                clean_shutdown(sub_observer)
                return 0, final_list_dic
        else:
            if use_solar_mask and os.path.exists(fits_mask):
                os.system("rm -rf " + fits_mask)
            logger.info(
                f"{os.path.basename(msname)} -- No image is made.\n",
            )
            time.sleep(5)
            clean_shutdown(sub_observer)
            return 1, {}
    except Exception:
        logger.info(
            f"{os.path.basename(msname)} -- Error in imaging.\n",
        )
        logger.exception(traceback.print_exc())
        time.sleep(5)
        clean_shutdown(sub_observer)
        return 1, {}
    finally:
        time.sleep(5)
        drop_cache(msname)


def run_all_imaging(
    mslist,
    dask_client,
    workdir="",
    outdir="",
    freqrange="",
    timerange="",
    datacolumn="CORRECTED_DATA",
    freqres=-1,
    timeres=-1,
    weight="briggs",
    robust=0.0,
    minuv=0,
    pol="I",
    threshold=1.0,
    use_multiscale=True,
    use_solar_mask=True,
    imaging_params={},  # TODO
    savemodel=False,
    saveres=False,
    cutout_rsun=10.0,
    make_plots=True,
    cpu_frac=0.8,
    mem_frac=0.8,
    logfile="imaging.log",
):
    """
    Run spectropolarimetric snapshot imaging on a list of measurement sets

    Parameters
    ----------
    mslist : list
        Measurement set list
    dask_client : dask.client
        Dask client
    workdir : str
        Work directory
    outdir : str
        Output directory
    freqrange : str, optional
        Frequency range to image
    timerange : str, optional
        Time range
    datacolumn : str, optional
        Data column
    freqres : float, optional
        Frequency resolution of spectral chunk in MHz
    timeres : float, optional
        Time resolution of temporal chunk in seconds
    weight : str, optional
        Image weighting
    robust : float, optional
        Briggs weighting robust parameter
    minuv : float, optional
        Minimum UV-lambda to use in imaging
    pol : str, optional
        Stokes parameters to image
    threshold : float, optional
        CLEAN threshold
    use_multiscale : bool, optional
        Use multiscale or not
    use_solar_mask : bool, optional
        Use solar mask
    savemodel : bool, optional
        Save model images or not
    saveres : bool, optional
        Save residual images or not
    cutout_rsun : float, optional
        Cutout image size (width and height is : 2 times cutout_rsun)
        Default value: 10 solar radii
        Note: default FoV is 20 solar solar radii. If cutout_rsun is chosen larger than 20 solar radii, FoV will be increased accordingly.
    make_plots : bool, optional
        Make radio image helioprojective plots
    cpu_frac : float, optional
        CPU fraction to use
    mem_frac : float, optional
        Memory fraction to use

    Returns
    -------
    int
        Success message
    int
        Succeeded ms number
    int
        Failed ms number
    int
        Total images
    """
    cpu_frac = min(0.8, abs(cpu_frac))
    mem_frac = min(0.8, abs(mem_frac))

    mslist = sorted(mslist)

    if len(mslist) == 0:
        print("Please provide a valid measurement set list.")
        return 1, 0, 0
    else:
        succeed = 0
        failed = len(mslist)
        total_images = 0

    try:
        if weight == "briggs":
            weight_str = f"{weight}_{robust}"
        else:
            weight_str = weight
        if freqres == -1 and timeres == -1:
            imagedir = outdir + f"/imagedir_f_all_t_all_pol_{pol}_w_{weight_str}"
        elif freqres != -1 and timeres == -1:
            imagedir = outdir + f"/imagedir_f_{freqres}_t_all_pol_{pol}_w_{weight_str}"
        elif freqres == -1 and timeres != -1:
            imagedir = outdir + f"/imagedir_f_all_t_{timeres}_pol_{pol}_w_{weight_str}"
        else:
            imagedir = (
                outdir + f"/imagedir_f_{freqres}_t_{timeres}_pol_{pol}_w_{weight_str}"
            )
        os.makedirs(imagedir, exist_ok=True)

        ####################################
        # Filtering any corrupted ms
        #####################################
        filtered_mslist = []  # Filtering in case any ms is corrupted
        for ms in mslist:
            checkcol = check_datacolumn_valid(ms)
            if checkcol:
                filtered_mslist.append(ms)
            else:
                print(f"Issue in : {ms}")

        mslist = filtered_mslist
        if len(mslist) == 0:
            print("No valid measurement set is found.")
            return 1, succeed, failed, total_images

        client_info = dask_client.scheduler_info()["workers"]
        njobs = len(client_info)
        worker_mem_list = []
        for addr, w in client_info.items():
            worker_mem_list.append(w["memory_limit"] / 1024**3)
        mem_limit = round(min(worker_mem_list), 3)
        n_threads = os.environ.get("OMP_NUM_THREADS")
        if n_threads is not None:
            n_threads = int(n_threads)
        else:
            n_threads = 1

        print("#################################")
        print(f"Total dask worker: {njobs}")
        print(f"CPU per worker: {n_threads}")
        print(f"Memory per worker: {mem_limit} GB")
        print("#################################")

        cutout_rsun = max(5, cutout_rsun) # Minimum 5 solar radii cutout, default is 10 solar radii
        
        tasks = []
        for i in range(len(mslist)):
            ms = mslist[i]
            num_pixel_in_psf = calc_npix_in_psf(weight, robust=robust)
            cellsize = calc_cellsize(ms, num_pixel_in_psf)
            instrument_fov = calc_field_of_view(ms, FWHM=False)
            msmd = msmetadata()
            msmd.open(ms)
            freqMHz = msmd.meanfreq(0, unit="MHz")
            msmd.close()
            cutout_rsun_arcsec = cutout_rsun*16*60
            fov = min(instrument_fov, 2*cutout_rsun_arcsec)
            imsize = int(fov / cellsize)
            imsize = get_fft_size(imsize)
            os.makedirs(workdir + "/logs", exist_ok=True)
            logfile = (
                workdir
                + "/logs/imaging_"
                + os.path.basename(ms).split(".ms")[0]
                + ".log"
            )
            tasks.append(
                delayed(perform_imaging)(
                    msname=ms,
                    workdir=workdir,
                    freqrange=freqrange,
                    timerange=timerange,
                    datacolumn=datacolumn,
                    imagedir=imagedir,
                    imsize=imsize,
                    cellsize=cellsize,
                    image_freqres=freqres,
                    image_timeres=timeres,
                    pol=pol,
                    weight=weight,
                    robust=robust,
                    minuv=minuv,
                    threshold=threshold,
                    use_multiscale=use_multiscale,
                    use_solar_mask=use_solar_mask,
                    savemodel=savemodel,
                    saveres=saveres,
                    cutout_rsun=cutout_rsun,
                    make_plots=make_plots,
                    ncpu=n_threads,
                    mem=mem_limit,
                    logfile=logfile,
                )
            )
        print(
            f"Starting imaging for ms : {ms}, Log file : {logfile}",
        )
        results = list(dask_client.gather(dask_client.compute(tasks)))
        all_image_list = []
        all_imaged_ms_list = []
        for i in range(len(results)):
            r = results[i][1]
            if len(r) == 0:
                print(
                    f"Imaging failed for ms : {mslist[i]}",
                )
            else:
                image_list = r["image"]
                if len(image_list) == 0:
                    print(
                        f"No image is made for ms : {mslist[i]}",
                    )
                else:
                    all_imaged_ms_list.append(mslist[i])
                    for image in image_list:
                        all_image_list.append(image)

        succeed = len(all_imaged_ms_list)
        failed = len(mslist) - succeed
        total_images = len(all_image_list)
        print(
            f"Numbers of input measurement sets : {len(mslist)}.",
        )
        print(
            f"Imaging successfully done for: {succeed} measurement sets.",
        )
        print(
            f"Imaging failed for: {failed} measurement sets.",
        )
        print(f"Total images made: {total_images}.")
        return 0, succeed, failed, total_images
    except Exception:
        traceback.print_exc()
        return 1, succeed, failed, total_images


def main(
    mslist,
    workdir,
    outdir,
    freqrange="",
    timerange="",
    datacolumn="CORRECTED_DATA",
    pol="I",
    freqres=-1,
    timeres=-1,
    weight="briggs",
    robust=0.0,
    minuv=0.0,
    threshold=1.0,
    cutout_rsun=10.0,
    use_multiscale=True,
    use_solar_mask=True,
    savemodel=True,
    saveres=True,
    make_plots=True,
    start_remote_log=False,
    cpu_frac=0.8,
    mem_frac=0.8,
    logfile=None,
    jobid=0,
    dask_client=None,
):
    """
    Perform distributed spectropolarimetric snapshot imaging on multiple measurement sets.

    Parameters
    ----------
    mslist : str
        Comma-separated list of measurement set paths to be imaged.
    workdir : str
        Directory for intermediate files, logs, and temporary outputs.
    outdir : str
        Directory where final images, models, and plots will be saved.
    freqrange : str, optional
        Frequency range to image (e.g., "500~1000MHz"). Default is "" (all frequencies).
    timerange : str, optional
        Time range to image (e.g., "09:30:00~09:40:00"). Default is "" (all times).
    datacolumn : str, optional
        Data column to image (e.g., "DATA", "CORRECTED_DATA"). Default is "CORRECTED_DATA".
    pol : str, optional
        Polarization product to image (e.g., "I", "XX", "RR", "QUV"). Default is "I".
    freqres : float, optional
        Frequency resolution in MHz for slicing the MS. Use -1 to disable. Default is -1.
    timeres : float, optional
        Time resolution in seconds for snapshot imaging. Use -1 to disable. Default is -1.
    weight : str, optional
        Weighting scheme for imaging ("natural", "uniform", "briggs"). Default is "briggs".
    robust : float, optional
        Robustness parameter for Briggs weighting. Default is 0.0.
    minuv : float, optional
        Minimum uv-distance (in wavelengths) to include in imaging. Default is 0.0.
    threshold : float, optional
        Cleaning threshold in sigma. Default is 1.0.
    cutout_rsun : float, optional
        Radius in solar radii to cut out around solar center. Default is 10.0.
    use_multiscale : bool, optional
        If True, enables multiscale CLEAN deconvolution. Default is True.
    use_solar_mask : bool, optional
        If True, applies a solar disk mask during CLEAN to reduce sidelobe artifacts. Default is True.
    savemodel : bool, optional
        If True, saves the CLEAN model images. Default is True.
    saveres : bool, optional
        If True, saves the residual images. Default is True.
    make_plots : bool, optional
        If True, generates diagnostic plots for each image. Default is True.
    start_remote_log : bool, optional
        Whether to enable remote logging using credentials in the workdir. Default is False.
    cpu_frac : float, optional
        Fraction of total CPUs to use per task. Default is 0.8.
    mem_frac : float, optional
        Fraction of total system memory to use per task. Default is 0.8.
    logfile : str, optional
        Log file
    jobid : int, optional
        Unique job identifier for logging and PID tracking. Default is 0.
    dask_client : dask.client, optional
        Dask client

    Returns
    -------
    int
        Success message
    int
        Succeeded ms number
    int
        Failed ms number
    int
        Total images
    """
    cpu_frac = min(0.8, abs(cpu_frac))
    mem_frac = min(0.8, abs(mem_frac))

    mslist = mslist.split(",")

    if workdir == "":
        workdir = os.path.dirname(os.path.abspath(mslist[0])) + "/workdir"
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)

    if outdir == "" or not os.path.exists(outdir):
        outdir = workdir
    outdir = outdir.rstrip("/")
    os.makedirs(outdir, exist_ok=True)

    if len(mslist) == 0:
        print("Please provide a valid measurement set list.")
        return 1, 0, 0, 0
    else:
        succeed = 0
        failed = len(mslist)
        total_images = 0

    ###########################
    # WSClean container
    ###########################
    container_name = "paircarswsclean"
    container_present = check_udocker_container(container_name)
    if not container_present:
        print(f"Initializing {container_name}...")
        container_name = initialize_wsclean_container(name=container_name, verbose=True)
        if container_name is None:
            print(
                f"Container {container_name} is not initiated. First initiate container and then run."
            )
            return 1, succeed, failed, total_images

    ############
    # Logger
    ############
    observer = None
    if (
        start_remote_log
        and os.path.exists(f"{workdir}/.jobname_password.npy")
        and logfile is not None
    ):
        time.sleep(5)
        jobname, password = np.load(
            f"{workdir}/.jobname_password.npy", allow_pickle=True
        )
        if os.path.exists(logfile):
            observer = init_logger(
                "all_imaging", logfile, jobname=jobname, password=password
            )
    if observer is None:
        print("Remote link or jobname is blank. Not transmiting to remote logger.")

    total_ncoarse = 0
    for msname in mslist:
        ncoarse = get_ncoarse(msname)
        total_ncoarse += ncoarse
    total_ncoarse = max(1, total_ncoarse)

    dask_cluster = None
    if dask_client is None:
        if mem_frac <= 0:
            mem_frac = 0.8
        if cpu_frac <= 0:
            cpu_frac = 0.8
        target_ms_sizes = [get_ms_size(msname) for msname in mslist]
        max_ms_size = max(target_ms_sizes)
        min_mem = round(10 * max_ms_size, 2)  # 10 times the size of the ms
        min_mem /= total_ncoarse

        result = get_local_dask_cluster(
            workdir,
            cpu_frac=cpu_frac,
            mem_frac=mem_frac,
            min_mem=min_mem,
            max_worker=len(mslist) + 1,
        )
        if result is None:
            print("Error occured in creating local cluster.")
            return 1, succeed, failed, total_images
        else:
            dask_client, dask_cluster, dask_dir, nworker = result
        scale_worker_and_wait(dask_cluster, dask_client, nworker)

    try:
        msg, succeed, failed, total_images = run_all_imaging(
            mslist,
            dask_client,
            workdir=workdir,
            outdir=outdir,
            freqrange=freqrange,
            timerange=timerange,
            datacolumn=datacolumn,
            freqres=freqres,
            timeres=timeres,
            weight=weight,
            robust=robust,
            minuv=minuv,
            threshold=threshold,
            use_multiscale=use_multiscale,
            use_solar_mask=use_solar_mask,
            pol=pol,
            make_plots=make_plots,
            cutout_rsun=cutout_rsun,
            savemodel=savemodel,
            saveres=saveres,
            cpu_frac=cpu_frac,
            mem_frac=mem_frac,
        )
    except Exception:
        traceback.print_exc()
        msg = 1
    finally:
        time.sleep(5)
        clean_shutdown(observer)
        for ms in mslist:
            drop_cache(ms)
        if dask_cluster is not None:
            dask_client.shutdown()
            dask_client.close()
            dask_cluster.close()
            drop_cache(workdir)
            drop_cache(outdir)
            os.system(f"rm -rf {dask_dir}")
    return msg, succeed, failed, total_images


def cli():
    parser = argparse.ArgumentParser(
        description="Perform spectropolarimetric snapshot imaging",
        formatter_class=SmartDefaultsHelpFormatter,
    )

    # Essential parameters
    basic_args = parser.add_argument_group(
        "###################\nEssential parameters\n###################"
    )
    basic_args.add_argument(
        "mslist",
        type=str,
        help="Comma-separated list of measurement sets (required)",
    )
    basic_args.add_argument(
        "--workdir",
        type=str,
        required=True,
        default="",
        help="Work directory for imaging",
    )
    basic_args.add_argument(
        "--outdir",
        type=str,
        default="",
        help="Output directory for imaging products",
    )

    # Advanced parameters
    adv_args = parser.add_argument_group(
        "###################\nAdvanced imaging parameters\n###################"
    )
    adv_args.add_argument(
        "--freqrange",
        type=str,
        default="",
        help="Frequency range to image",
    )
    adv_args.add_argument(
        "--timerange",
        type=str,
        default="",
        help="Time range to image",
    )
    adv_args.add_argument(
        "--datacolumn",
        type=str,
        default="CORRECTED_DATA",
        help="Data column to use for imaging",
    )
    adv_args.add_argument(
        "--pol",
        type=str,
        default="I",
        help="Stokes parameters to image",
    )
    adv_args.add_argument(
        "--freqres",
        type=float,
        default=-1,
        help="Frequency resolution per chunk in MHz (-1 for full)",
    )
    adv_args.add_argument(
        "--timeres",
        type=float,
        default=-1,
        help="Time resolution per chunk in seconds (-1 for full)",
    )
    adv_args.add_argument(
        "--weight",
        type=str,
        default="briggs",
        help="Imaging weighting scheme",
    )
    adv_args.add_argument(
        "--robust",
        type=float,
        default=0.0,
        help="Briggs robust parameter",
    )
    adv_args.add_argument(
        "--minuv_l",
        dest="minuv",
        type=float,
        default=0.0,
        help="Minimum UV distance in wavelengths",
    )
    adv_args.add_argument(
        "--threshold",
        type=float,
        default=1.0,
        help="CLEAN threshold in sigma",
    )
    adv_args.add_argument(
        "--cutout_rsun",
        type=float,
        default=10.0,
        help="Cutout radius for images (solar radii)",
    )
    adv_args.add_argument(
        "--no_multiscale",
        action="store_false",
        dest="use_multiscale",
        help="Do not use multiscale CLEAN",
    )
    adv_args.add_argument(
        "--no_solar_mask",
        action="store_false",
        dest="use_solar_mask",
        help="Do not use solar disk mask for CLEANing",
    )
    adv_args.add_argument(
        "--no_savemodel",
        action="store_false",
        dest="savemodel",
        help="Do no save model images",
    )
    adv_args.add_argument(
        "--no_saveres",
        action="store_false",
        dest="saveres",
        help="Do not save residual images",
    )
    adv_args.add_argument(
        "--no_make_plots",
        action="store_false",
        dest="make_plots",
        help="Do not generate helioprojective plots",
    )
    adv_args.add_argument(
        "--start_remote_log", action="store_true", help="Start remote logging"
    )

    # Resource management parameters
    hard_args = parser.add_argument_group(
        "###################\nHardware resource management parameters\n###################"
    )
    hard_args.add_argument(
        "--cpu_frac",
        type=float,
        default=0.8,
        help="Fraction of available CPU to use",
    )
    hard_args.add_argument(
        "--mem_frac",
        type=float,
        default=0.8,
        help="Fraction of available memory to use",
    )
    hard_args.add_argument(
        "--jobid",
        type=str,
        default="0",
        help="Job ID for process tracking and logging",
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        return 1

    args = parser.parse_args()

    msg, _, _, _ = main(
        mslist=args.mslist,
        workdir=args.workdir,
        outdir=args.outdir,
        freqrange=args.freqrange,
        timerange=args.timerange,
        datacolumn=args.datacolumn,
        pol=args.pol,
        freqres=args.freqres,
        timeres=args.timeres,
        weight=args.weight,
        robust=args.robust,
        minuv=args.minuv,
        threshold=args.threshold,
        cutout_rsun=args.cutout_rsun,
        use_multiscale=args.use_multiscale,
        use_solar_mask=args.use_solar_mask,
        savemodel=args.savemodel,
        saveres=args.saveres,
        make_plots=args.make_plots,
        start_remote_log=args.start_remote_log,
        cpu_frac=float(args.cpu_frac),
        mem_frac=float(args.mem_frac),
        jobid=args.jobid,
    )
    return msg
