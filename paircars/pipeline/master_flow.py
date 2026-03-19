import traceback
import time
import glob
import sys
import os
import socket
import requests
import getpass
import contextlib
import numpy as np
import argparse
from casatools import msmetadata
from astropy.io import fits
from datetime import datetime as dt
from multiprocessing import Event
from dask.distributed import get_client
from prefect import flow
from functools import partial
from prefect.context import get_run_context
from prefect_dask.task_runners import DaskTaskRunner
from prefect.settings import get_current_settings
from paircars.utils.basic_utils import (
    get_cachedir,
    internet_available,
)
from paircars.utils.calibration import (
    calc_bw_smearing_freqwidth,
    calc_time_smearing_timewidth,
    max_time_solar_smearing,
    interpolate_bpass,
    interpolate_quartical,
    get_caltable_metadata,
)
from paircars.utils.casatasks import reset_weights_and_flags
from paircars.utils.flagging import do_flag_backup, get_chans_flag
from paircars.utils.image_utils import filter_images
from paircars.utils.logger_utils import (
    SmartDefaultsHelpFormatter,
    clean_shutdown,
    generate_password,
    get_remote_logger_link,
    get_emails,
    init_logger,
)
from paircars.utils.ms_metadata import check_datacolumn_valid
from paircars.utils.mwa_ploting_utils import (
    plot_caltable_diagnostics,
    plot_quartical_tables,
)
from paircars.utils.mwa_utils import (
    get_ncoarse,
    get_MWA_OBSID,
    download_MWA_metafits,
    get_selfcal_ntimes,
    freq_to_MWA_coarse,
)
from paircars.utils.proc_manage_utils import (
    get_jobid,
    save_main_process_info,
    scale_worker_and_wait,
    get_local_dask_cluster,
    get_scheduler_name,
)
from paircars.utils.resource_utils import drop_cache
from paircars.data.sendmail import (
    send_paircars_notification as send_notification,
)
from paircars.clusterutils.slurm_cluster import (
    get_slurm_dask_cluster,
    is_slurm_job,
)
from paircars.utils.prefect_logger_utils import (
    start_flow_log_saver,
)
from paircars.pipeline import (
    move_solarcenter,
)
from paircars.pipeline.init_data import init_paircars_data
from paircars.pipeline.flows import basic_cal_subflow
from paircars.pipeline.tasks import (
    run_solar_phasecenter_jobs,
    run_ds_jobs,
    run_target_split_jobs,
    run_flag,
    run_import_model,
    run_basic_cal_jobs,
    run_apply_basiccal_sol,
    run_solar_siderealcor_jobs,
    run_selfcal_jobs,
    run_apply_selfcal_sol,
    run_imaging_jobs,
    run_apply_pbcor,
    run_make_overlay,
    run_make_msplot,
    send_task_notification,
)


@flow(
    name="P-AIRCARS Master control",
    version="3.0",
    description="Calibration and Imaging Pipeline for MWA Solar Observation",
    log_prints=True,
)
def master_control(
    target_datadir,
    workdir,
    outdir,
    dask_addr,
    # Metafits and calibrators
    target_metafits="",
    calibrator_datadir="",
    calibrator_metafits="",
    solar_data=True,
    # Pre-calibration
    do_forcereset_weightflag=False,
    do_cal_flag=True,
    do_import_model=True,
    # Basic calibration
    do_basic_cal=True,
    do_applycal=True,
    only_amplitude=False,
    redo_basic_cal=False,
    use_solarflagger=False,
    # Target data preparation
    freqrange="",
    timerange="",
    uvrange="",
    # Polarization self-calibration
    do_polcal=False,
    # Self-calibration
    do_selfcal=True,
    do_apply_selfcal=True,
    do_ap_selfcal=True,
    solar_selfcal=True,
    use_solar_mask=True,
    solint="60s",
    redo_selfcal=False,
    # Sidereal correction
    do_sidereal_cor=False,
    do_move_solarcenter=True,
    # Dynamic spectra
    make_ds=True,
    # Imaging
    do_imaging=True,
    do_pbcor=True,
    weight="briggs",
    robust=0.0,
    minuv=0,
    image_freqres=1.28,
    image_timeres=10.0,
    pol="IQUV",
    clean_threshold=1.0,
    use_multiscale=True,
    cutout_rsun=10.0,
    make_overlay=False,
    make_msplot=False,
    # Resource settings
    cpu_frac=0.8,
    mem_frac=0.8,
    max_worker=2,
    keep_backup=False,
    keep_calibrated_ms=True,
    # Remote logging
    masterlog=None,
    remote_logger=False,
    jobid=None,
    job_password=None,
    adaptive=False,
):
    """
    Master controller of the entire pipeline

    Parameters
    ----------
    target_datadir : str
        Target measurement set directory
    workdir : str
        Work directory path
    outdir : str
        Output directory

    target_metafits : str
        Target metafits file
    calibrator_datadir : str, optional
        Calibrator data directory
    calibrator_metafits : str, optional
        Calibrator metafits file
    solar_data : bool, optional
        Whether it is solar data or not

    do_forcereset_weightflag : bool, optional
        Reset weights and flags of the input ms
    do_cal_flag : bool, optional
        Perform flagging on calibrator
    do_import_model : bool, optional
        Import model visibilities of flux and polarization calibrators

    do_basic_cal : bool, optional
        Perform basic calibration
    do_applycal : bool, optional
        Apply basic calibration on target scans
    only_amplitude : bool, optional
        Apply only amplitude part of gain solution from calibrator
    redo_basic_cal : bool, optional
        Redo basic calibration
    use_solarflagger : bool, optional
        Use solar flagger or not

    freqrange : str, optional
        Frequency range to image in MHz (xx1~xx2,xx3~xx4,)
    timerange : str, optional
        Time range to image in YYYY/MM/DD/hh:mm:ss format (tt1~tt2,tt3~tt4,...)
    uvrange : str, optional
        UV-range for calibration

    do_polcal : bool, optional
        Perform full-polarization calibration and imaging

    do_selfcal : bool, optional
        Perform self-calibration
    do_apply_selfcal : bool, optional
        Apply self-calibration solutions
    do_ap_selfcal : bool, optional
        Perform amplitude-phase self-cal or not
    solint : str, optional
        Solution intervals in self-cal
    solar_selfcal : bool, optional
        Solar selfcal
    use_solar_mask : bool, optional
        Use solar mask or not
    redo_selfcal : bool, optional
        Redo self-calibration or not

    do_sidereal_cor : bool, optional
        Perform solar sidereal motion correction or not
    do_move_solarcenter: boo, optional
        Move phasecenter to solar center
    make_ds : bool, optional
        Make dynamic spectra

    do_imaging : bool, optional
        Perform final imaging
    do_pbcor : bool, optional
        Perform primary beam correction
    weight : str, optional
        Image weighting
    robust : float, optional
        Robust parameter for briggs weighting (-1 to 1)
    minuv : float, optional
        Minimum UV-lambda for final imaging
    image_freqres : float, optional
        Image frequency resolution in MHz (-1 means full bandwidth)
    image_timeres : float, optional
        Image temporal resolution in seconds (-1 means full scan duration)
    pol : str, optional
        Stokes parameters of final imaging
    clean_threshold : float, optional
        CLEAN threshold of final imaging in sigma
    use_multiscale : bool, optional
        Use multiscale scales or not
    cutout_rsun : float, optional
        Cutout image size from center in solar radii (default : 10.0 solar radii)
    make_overlay : bool, optional
        Make EUV MWA overlay for all images or not (default : per coarse channel images will be overlaid at 60s intervals)
    make_msplot : bool, optional
        Make diagnostic plots of measurement sets

    cpu_frac : float, optional
        CPU fraction to use
    mem_frac : float, optional
        Memory fraction to use
    max_worker: int, optional
        Maximum workers
    keep_backup : bool, optional
        Keep backup of of all intermediate data poducts, calibrated ms, self-cal rounds and including images, models and residual images
    keep_calibrated_ms : bool, optional
        Keep calibrated measurement sets or not

    masterlog : str, optional
        Master logfile
    remote_logger : bool, optional
        Enable remote logging of the pipeline status
    jobid : str, optional
        Job ID
    job_password : str, optional
        User specified job password for remote logger
    adaptive : bool, optional
        Whether do adaptive scaling or not

    Returns
    -------
    int
        Success message
    """
    print("P-AIRCARS workflow started...")
    emails = get_emails()

    #######################################################
    # Checking validity of target directories and metafits
    #######################################################
    if target_datadir.startswith("~"):
        print("Please provide full path of target directory.")
        return 1
    else:
        target_datadir = os.path.abspath(target_datadir)
        if os.path.exists(target_datadir) is False:
            print(
                f"Target data directory: {target_datadir} does not exist. Provide correct full path."
            )
            return 1

    #################################################
    # Checking workdir and outdir paths
    #################################################
    if workdir.startswith("~"):
        print("Please provide full path of work directory.")
        return 1
    else:
        workdir = os.path.abspath(workdir)
    if outdir.startswith("~"):
        print("Please provide full path of output directory.")
        return 1
    else:
        outdir = os.path.abspath(outdir)

    if jobid is None:
        jobid = get_jobid()

    print("#############################")
    print(f"P-AIRCARS Job ID: {jobid}")
    print("#############################")
    #########################################
    # Some validity checks for resources
    #########################################
    max_worker = max(2, max_worker)  # Minimum 2 workers are needed
    cpu_frac = min(0.8, abs(cpu_frac))
    mem_frac = min(0.8, abs(mem_frac))

    print("Sorting out target data....")
    #############################################
    # Listing target ms
    #############################################
    target_mslist = sorted(glob.glob(f"{target_datadir}/*.ms"))
    if len(target_mslist) == 0:
        print(
            f"No measurement set is present in target data directory: {target_datadir}"
        )
        if emails != "":
            email_msg = "No measurement set is present in the target data directory."
            send_task_notification(emails, email_msg, jobid, "N/A", "N/A")
        return 1

    #################################################
    # Verifying whether all target ms from same obsid
    #################################################
    target_ms_obsids = [get_MWA_OBSID(ms) for ms in target_mslist]
    all_same_obsids = all(x == target_ms_obsids[0] for x in target_ms_obsids)
    if not all_same_obsids:
        print(
            "All target measurement sets are not belong to same OBSID. Keep only measurement sets with same OBSID inside the target directory. P-AIRCARS has stopped."
        )
        return 1
    else:
        target_obsid = target_ms_obsids[0]

    ##############################################
    # Downloading target metafits if not exist
    ##############################################
    if target_metafits == "" or not os.path.exists(target_metafits):
        if os.path.exists(f"{target_datadir}/{target_obsid}.metafits"):
            target_metafits = f"{target_datadir}/{target_obsid}.metafits"
        else:
            try:
                target_metafits = download_MWA_metafits(
                    target_obsid, outdir=target_datadir
                )
            except Exception:
                traceback.print_exc()
                if emails != "":
                    email_msg = f"Target metafits for OBSID: {target_obsid} is not provided and also could not be downloaded. P-AIRCARS has stopped."
                    send_task_notification(emails, email_msg, jobid, "N/A", "N/A")
                print(
                    f"Target metafits for OBSID: {target_obsid} is not provided and also could not be downloaded. P-AIRCARS has stopped."
                )
                return 1

    ##################################################
    # Downloading target metafits if not match with ms
    ##################################################
    metafits_obsid = fits.getheader(target_metafits)["GPSTIME"]
    if metafits_obsid != target_obsid:
        print(
            "Mismatch between target ms OBSID: {target_obsid} and metafits OBSID: {metafits_obsid}. Downloading metafits for OBSID: {target_obsid}."
        )
        try:
            target_metafits = download_MWA_metafits(target_obsid, outdir=target_datadir)
        except Exception:
            traceback.print_exc()
            if emails != "":
                email_msg = f"Target metafits for OBSID: {target_obsid} could not be downloaded. P-AIRCARS has stopped."
                send_task_notification(emails, email_msg, jobid, "N/A", "N/A")
            print(
                f"Target metafits for OBSID: {target_obsid} could not be downloaded. P-AIRCARS has stopped."
            )
            return 1

    ################################################
    # Final target OBSID and frequency configuration
    ################################################
    target_header = fits.getheader(target_metafits)
    target_obsid = int(target_header["GPSTIME"])
    target_freq_config = target_header["CHANNELS"]
    target_coarse_chans = [int(c) for c in target_freq_config.split(",")]
    print(f"Target observation ID: {target_obsid}")
    print(f"Target coarse channels: {target_coarse_chans}")

    ################################################
    # Filtering calibrators
    ################################################
    print("Sorting out calibrator data....")
    cal_datadir_list = calibrator_datadir.split(",")
    cal_metafits_list = calibrator_metafits.split(",")
    final_cal_datadir_list = []
    final_cal_obsid_list = []
    final_cal_metafits_list = []
    final_cal_coarsechan_list = []

    for i in range(len(cal_datadir_list)):
        cal_datadir = cal_datadir_list[i]
        #############################################
        # Listing calibrator ms
        #############################################
        if cal_datadir != "" and os.path.exists(cal_datadir):
            cal_mslist = sorted(glob.glob(f"{cal_datadir}/*.ms"))
            if len(cal_mslist) == 0:
                print(
                    f"No measurement set is present in calibrator data directory: {cal_datadir}"
                )
                has_cal = False
            else:
                has_cal = True
        else:
            has_cal = False
        ######################################################
        # Verifying whether all calibraror ms from same obsid
        ######################################################
        if has_cal:
            cal_ms_obsids = [get_MWA_OBSID(ms) for ms in cal_mslist]
            all_same_obsids = all(x == cal_ms_obsids[0] for x in cal_ms_obsids)
            if not all_same_obsids:
                print(
                    "All calibrator measurement sets are not belong to same OBSID. Not using this calibrator."
                )
                has_cal = False
            else:
                cal_obsid = cal_ms_obsids[0]
        ##############################################
        # Searching for calibrator metafits
        ##############################################
        cal_metafits = None
        if has_cal and len(cal_metafits_list) > 0:
            for metafits in cal_metafits_list:
                if (metafits == "" or not os.path.exists(metafits)) and os.path.exists(
                    f"{cal_datadir}/{cal_obsid}.metafits"
                ):
                    metafits = f"{cal_datadir}/{cal_obsid}.metafits"
                metafits_obsid = fits.getheader(metafits)["GPSTIME"]
                if metafits_obsid == cal_obsid:
                    cal_metafits = metafits
                    break
        ######################################################
        # Downloading calibrator metafits if not match with ms
        ######################################################
        if has_cal and cal_metafits is None:
            try:
                cal_metafits = download_MWA_metafits(cal_obsid, outdir=cal_datadir)
            except Exception:
                traceback.print_exc()
                print(
                    f"Calibrator metafits for OBSID: {cal_obsid} could not be downloaded. Not using this calibrator."
                )
                has_cal = False
        ####################################################################
        # Including in final list if have overlapping frequency with target
        ####################################################################
        if has_cal and cal_metafits is not None:
            cal_header = fits.getheader(cal_metafits)
            cal_obsid = int(cal_header["GPSTIME"])
            if (
                abs(cal_obsid - target_obsid) < 12 * 3600
            ):  # Only if calibrator is 12 hours apart
                cal_freq_config = cal_header["CHANNELS"]
                cal_coarse_chans = [int(c) for c in cal_freq_config.split(",")]
                has_overlap = bool(set(cal_coarse_chans) & set(target_coarse_chans))
                if not has_overlap:
                    print(
                        f"Calibrator with OBSID: {cal_obsid} do not have frequency overlap with target."
                    )
                    print(f"Target coarse channels: {target_coarse_chans}")
                    print(f"Calibrator coarse channels: {cal_coarse_chans}")
                else:
                    final_cal_datadir_list.append(cal_datadir)
                    final_cal_obsid_list.append(cal_obsid)
                    final_cal_metafits_list.append(cal_metafits)
                    final_cal_coarsechan_list.append(cal_coarse_chans)

    ####################################################
    # Check whether there is calibrator available or not
    ####################################################
    calibrator_dic = {}
    if len(final_cal_datadir_list) > 0:
        ####################################################
        # Arranging in time
        ####################################################
        final_cal_datadir_list = np.array(final_cal_datadir_list)
        final_cal_obsid_list = np.array(final_cal_obsid_list)
        final_cal_metafits_list = np.array(final_cal_metafits_list)
        final_cal_coarsechan_list = np.array(final_cal_coarsechan_list)
        pos = np.argsort(abs(final_cal_obsid_list - target_obsid))
        final_cal_datadir_list = final_cal_datadir_list[pos].tolist()
        final_cal_obsid_list = final_cal_obsid_list[pos].tolist()
        final_cal_metafits_list = final_cal_metafits_list[pos].tolist()
        final_cal_coarsechan_list = final_cal_coarsechan_list[pos].tolist()
        all_coarse_chans = []
        for i in range(len(final_cal_datadir_list)):
            cal_obsid = final_cal_obsid_list[i]
            cal_datadir = final_cal_datadir_list[i]
            cal_metafits = final_cal_metafits_list[i]
            coarse_chans = final_cal_coarsechan_list[i]
            overlapping_coarse_chans = list(
                set(coarse_chans) & set(target_coarse_chans)
            )
            filtered_overlapping_coarse_chans = []
            for c in overlapping_coarse_chans:
                if c not in all_coarse_chans:
                    filtered_overlapping_coarse_chans.append(c)
                    all_coarse_chans.append(c)
            if len(filtered_overlapping_coarse_chans) > 0:
                calibrator_dic[cal_obsid] = [
                    cal_datadir,
                    cal_metafits,
                    filtered_overlapping_coarse_chans,
                ]

    if len(calibrator_dic) == 0:
        has_cal = False
        print("No calibrator data is found.")
    else:
        has_cal = True
        print(
            f"Total {len(calibrator_dic)} calibrator observations are sorted. Observation ID(s) are: {list(calibrator_dic.keys())}"
        )
        
    ######################################################
    # Making calibrator output directories
    ######################################################
    if has_cal:
        cal_outdir = f"{outdir}/calibrators"
        try:
            os.makedirs(cal_outdir, exist_ok=True)
        except Exception:
            print(
                f"Calibrator output directory: {cal_outdir} can not created. Please check the path carefully."
            )
            traceback.print_exc()
            has_cal = False
        basic_caldir = f"{cal_outdir}/caltables"
        os.makedirs(basic_caldir, exist_ok=True)

    #######################################
    # Preparing target working directories
    #######################################
    print("Preparing working directories....")
    if workdir == "":
        workdir = os.path.dirname(os.path.abspath(target_mslist[0])) + "/workdir"

    workdir = workdir.rstrip("/")
    if outdir == "":
        outdir = workdir

    workdir = f"{workdir}/{target_obsid}_{jobid}_target"
    try:
        os.makedirs(workdir, exist_ok=True)
    except Exception:
        print(
            f"Work directory: {workdir} can not be created. Please check the path carefully."
        )
        traceback.print_exc()
        return 1

    ####################################
    # Preparing target output directories
    ####################################
    outdir = outdir.rstrip("/")
    target_outdir = f"{outdir}/{target_obsid}_target"
    try:
        os.makedirs(target_outdir, exist_ok=True)
    except Exception:
        print(
            f"Output directory: {target_outdir} can not created. Please check the path carefully."
        )
        traceback.print_exc()
        return 1
    selfcaldir = f"{target_outdir}/caltables"
    os.makedirs(selfcaldir, exist_ok=True)

    #############################################
    # Change to workdir and determining scheduler
    #############################################
    os.chdir(workdir)
    scheduler_name = get_scheduler_name()

    #################################
    # Setup logger
    #################################
    logdir = f"{workdir}/logs"
    os.makedirs(logdir, exist_ok=True)

    if (
        scheduler_name != "local"
        or masterlog is None
        or os.path.exists(masterlog) is False
    ):
        master_logfile = f"{logdir}/main.log"
        ctx = get_run_context()
        flow_id = str(ctx.flow_run.id)
        flow_name = ctx.flow_run.name
        stop_event = Event()
        log_thread_flow = start_flow_log_saver(
            flow_id, flow_name, master_logfile, poll_interval=3, stop_event=stop_event
        )
        master_log_created = True
    else:
        master_log_created = False
        master_logfile = f"{logdir}/main.log"
        if os.path.exists(master_logfile):
            os.system(f"rm -rf {master_logfile}")
        os.symlink(masterlog, master_logfile)

    #####################################
    # Setup dask client
    #####################################
    dask_dir = None
    try:
        dask_client = get_client()
        dask_cluster = dask_client.cluster
    except Exception:
        dask_client, dask_cluster, dask_dir, nworker = get_local_dask_cluster(
            workdir,
            cpu_frac=cpu_frac,
            mem_frac=mem_frac,
            max_worker=max_worker,
        )
        if dask_client is None:
            print("Error occured in creating local cluster.")
            return 1

    #####################################
    # Initiating paircars data
    #####################################
    init_paircars_data()

    ############################################
    # Determine number of threads of main worker
    ############################################
    observer = None
    n_threads = os.environ.get("OMP_NUM_THREADS")
    if n_threads is None:
        n_threads = 1
    else:
        n_threads = max(1, int(n_threads))

    #########################################
    # Setup remote loggger and email notifier
    #########################################
    print("Setting up remote logger and email notifier...")
    try:
        #####################################
        # Reading remotelink and emails
        #####################################
        remote_link = ""
        internet_on = internet_available()
        if not internet_on:
            print("Internet connection is not available for remote logging.")
        else:
            if remote_logger:
                trial = 0
                while trial <= 5:
                    try:
                        remote_link = get_remote_logger_link()
                    except Exception:
                        traceback.print_exc()
                        pass
                    if remote_link != "":
                        break
                    else:
                        time.sleep(5)
                        trial += 1
                if remote_link == "":
                    print("Please provide a valid remote link.")
                    remote_logger = False

        if not remote_logger:
            timestamp = dt.utcnow().strftime("%Y-%m-%dT%H:%M:%S")
            internet_on = internet_available()
            if internet_on and emails != "":
                email_subject = (
                    f"P-AIRCARS Logger Details: {timestamp}, OBSID: {target_obsid}"
                )

                email_msg = f"P-AIRCARS Job ID: {jobid}"
                success_msg, error_msg = send_notification(
                    emails, email_subject, email_msg
                )
        else:
            ####################################
            # Job name and logging password
            ####################################
            hostname = socket.gethostname()
            timestamp = dt.utcnow().strftime("%Y-%m-%dT%H:%M:%S")
            username = getpass.getuser()
            jobname = f"{username}-{hostname}:{timestamp}:{target_obsid}"
            if job_password is None:
                password = generate_password()
            else:
                password = job_password
            np.save(
                f"{workdir}/.jobname_password.npy",
                np.array([jobname, password], dtype="object"),
            )
            ############
            # Logger
            ############
            if os.path.exists(f"{workdir}/.jobname_password.npy"):
                time.sleep(5)
                jobname, password = np.load(
                    f"{workdir}/.jobname_password.npy", allow_pickle=True
                )
                if master_logfile is not None and os.path.exists(master_logfile):
                    observer = init_logger(
                        "master_log", master_logfile, jobname=jobname, password=password
                    )
            if observer is None:
                print(
                    "Remote link or jobname is blank. Not transmiting to remote logger."
                )
            #####################
            # Notify over email
            #####################
            internet_on = internet_available()
            if emails != "" and internet_on:
                email_subject = (
                    f"P-AIRCARS Logger Details: {timestamp}, OBSID: {target_obsid}"
                )

                email_msg = (
                    f"P-AIRCARS Job ID: {jobid}\n"
                    f"Remote logger Job ID: {jobname}\n"
                    f"Remote access password: {password}"
                )
                success_msg, error_msg = send_notification(
                    emails, email_subject, email_msg
                )

        #####################################
        # Printing basic info of the pipeline
        #####################################
        print("###########################")
        print(f"Work directory: {workdir}")
        print(f"Final product directory: {outdir}")
        print("###########################")
        if remote_logger:
            print(
                "############################################################################"
            )
            print(remote_link)
            print(f"Remote Job ID: {jobname}")
            print(f"Remote access password: {password}")
            print(
                "#############################################################################"
            )

        if not has_cal:
            print(
                f"No suitable calibrators are available for target OBSID: {target_obsid}."
            )
            if emails != "":
                email_msg = f"No suitable calibrators are available for target OBSID: {target_obsid}."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )

        ###########################################
        # Setting up mutual conditions
        ###########################################
        # Move solar center, if any of these conditions are met
        if do_selfcal or do_applycal or do_apply_selfcal or do_imaging:
            if not do_move_solarcenter:
                do_move_solarcenter = True

        # Switch on cal flag and import model, if basic cal is needed
        if do_basic_cal:
            if not do_cal_flag:
                do_cal_flag = True
            if not do_import_model:
                do_import_model = True

        # Switch on applycal if selfcal is requested
        if do_selfcal:
            if not do_applycal:
                do_applycal = True

        # Switch on applycal and apply selfcal if imaging is requested
        if do_imaging:
            if not do_applycal:
                do_applycal = True
            if not do_apply_selfcal:
                do_apply_selfcal = True

        #####################################
        # Settings for solar data
        #####################################
        if solar_data:
            if not use_solar_mask:
                print("Use solar mask during CLEANing.")
                use_solar_mask = True
            if not solar_selfcal:
                solar_selfcal = True
            full_FoV = False
        else:
            if use_solar_mask:
                print("Stop using solar mask during CLEANing.")
                use_solar_mask = False
            if solar_selfcal:
                solar_selfcal = False
            full_FoV = True

        #####################################################################
        # Checking if ms is full pol for polarization calibration and imaging
        #####################################################################
        if do_polcal:
            print(
                "Checking measurement set suitability for polarization calibration...."
            )
            for msname in target_mslist:
                msmd = msmetadata()
                msmd.open(msname)
                npol = msmd.ncorrforpol()[0]
                msmd.close()
                if npol < 4:
                    print(
                        f"Measurement set: {msname} is not full-polar. Do not performing polarization analysis."
                    )
                    do_polcal = False
                    break

        #################################################
        # Determining maximum allowed frequency averaging
        #################################################
        print("Estimating optimal frequency averaging....")
        max_freqres_list = []
        freqres_list = []
        msmd = msmetadata()
        for msname in target_mslist:
            max_freqres = calc_bw_smearing_freqwidth(msname, full_FoV=full_FoV)
            max_freqres_list.append(max_freqres)
            msmd.open(msname)
            freqres = msmd.chanres(0, unit="MHz")[0]
            msmd.close()
            freqres_list.append(freqres)
        freqres = min(freqres_list)
        if freqres > 0.16:
            print(
                f"Frequency resolution: {round(freqres*1000,1)}kHz is more than 160kHz. Assuming channel flagging is already done before averaing."
            )
        max_freqres = min(max_freqres_list)
        if image_freqres > 0:
            image_freqres = max(image_freqres, freqres)
            freqavg = round(min(image_freqres, max_freqres), 2)
        else:
            freqavg = freqres
        freqavg = min(0.16, freqavg)
        image_freqres = round(image_freqres, 2)

        total_ncoarse = 0
        for msname in target_mslist:
            ncoarse = get_ncoarse(msname)
            total_ncoarse += ncoarse
        total_ncoarse = max(1, total_ncoarse)

        ################################################
        # Determining maximum allowed temporal averaging
        ################################################
        print("Estimating optimal temporal averaging....")
        max_timeres_list = []
        timeres_list = []
        for msname in target_mslist:
            if solar_data:  # For solar data, it is assumed Sun is tracked.
                max_timeres = calc_time_smearing_timewidth(msname)
            else:
                max_timeres = min(
                    calc_time_smearing_timewidth(msname),
                    max_time_solar_smearing(msname),
                )
            max_timeres_list.append(max_timeres)
            msmd.open(msname)
            times = msmd.timesforspws(0)
            timeres = np.nanmean(np.diff(times))
            msmd.close()
            timeres_list.append(timeres)
        timeres = min(timeres_list)
        quack_timestamps = int(4.0 / timeres)
        max_timeres = min(max_timeres_list)
        if image_timeres > (2 * 3660):  # If more than 2 hours
            print(
                "Image time integration is more than 2 hours, which may cause smearing due to solar differential rotation."
            )
        if image_timeres > 0:
            image_timeres = max(image_timeres, timeres)
            timeavg = round(min(image_timeres, max_timeres), 2)
        else:
            timeavg = timeres
        timeavg = min(2.0, timeavg)
        image_timeres = round(image_timeres, 2)
        print(f"Frequency resolution: {freqres}MHz, time resolution: {timeres}s.")
        print(f"Frequency averaging: {freqavg}MHz, time averaging: {timeavg}s.")
        print(
            f"Imaging frequency resolution: {image_freqres}MHz, time resolution: {image_timeres}s."
        )

        #############################
        # Reset any previous weights
        #############################
        print("Resetting previous flags and weights....")
        if len(target_mslist) > 0:
            for msname in target_mslist:
                reset_weights_and_flags(
                    msname, n_threads=n_threads, force_reset=do_forcereset_weightflag
                )
        if has_cal:
            cal_obsids = calibrator_dic.keys()
            for cal_obsid in cal_obsids:
                cal_datadir = calibrator_dic[cal_obsid][0]
                calibrator_mslist = glob.glob(f"{cal_datadir}/*.ms")
                if len(calibrator_mslist) > 0:
                    for msname in calibrator_mslist:
                        reset_weights_and_flags(
                            msname,
                            n_threads=n_threads,
                            force_reset=do_forcereset_weightflag,
                        )
        print("Reset is done.")
        #################################

        if (move_solarcenter or make_ds) and adaptive:
            scale_worker_and_wait(
                dask_cluster,
                dask_client,
                max(2, min(len(target_mslist) + 1, max_worker)),
            )
        ########################################
        # Moving phasecenter to the solar center
        ########################################
        if solar_data and do_move_solarcenter:
            if emails != "":
                email_msg = "Started moving phasecenter to solar center."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print("Starting task: Moving phasecenter to the Sun .....")
            print("###########################")
            future_movecenter = run_solar_phasecenter_jobs.with_options(
                task_run_name=f"moving_to_solar_center_{jobid}",
            ).submit(
                ",".join(target_mslist),
                workdir,
                jobid=jobid,
                cpu_frac=round(cpu_frac, 2),
                mem_frac=round(mem_frac, 2),
                remote_log=remote_logger,
            )
            try:
                msg, succeed, failed = future_movecenter.result()
                if emails != "":
                    email_msg = f"Moving phasecenter to solar center is done.\nSucceeded: {succeed}, failed: {failed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print("Finished task: Moving phasecenter to solar center is done.")
                print("###########################")
                filtered_ms = []
                for t_ms in target_mslist:
                    t_ms = t_ms.rstrip("/")
                    if os.path.exists(f"{t_ms}/.solarcenter_move_succeed"):
                        filtered_ms.append(t_ms)
                    else:
                        print(
                            f"Issue in moving phasecneter to solar center for measurement set: {t_ms}"
                        )
                if adaptive and len(filtered_ms) != len(target_mslist):
                    scale_worker_and_wait(
                        dask_cluster,
                        dask_client,
                        max(2, min(len(target_mslist) + 1, max_worker)),
                    )
                target_mslist = filtered_ms  # Filtered target mslist
            except Exception:
                print(
                    "Error in moving phasecenter to solar center. P-AIRCARS has stopped."
                )
                traceback.print_exc()
                if emails != "":
                    email_msg = "Error occured in moving phasecenter to solar center."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                return 1

        #######################################
        # Run dynamic spectra making
        #######################################
        if solar_data and make_ds:
            if emails != "":
                email_msg = "Started making solar dynamic spectra."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print("Starting task: Making dynamic spectra of solar target .....")
            print("###########################")
            future_maskms = run_ds_jobs.with_options(
                task_run_name=f"making_dynamic_spectra_{jobid}",
            ).submit(
                ",".join(target_mslist),
                target_metafits,
                workdir,
                target_outdir,
                jobid=jobid,
                cpu_frac=round(cpu_frac, 2),
                mem_frac=round(mem_frac, 2),
                remote_log=remote_logger,
            )
            try:
                msg, succeed, failed = future_maskms.result()
                if emails != "":
                    email_msg = f"Making solar dynamic spectra are done.\nSucceeded: {succeed}, failed: {failed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print("Finished task: Making solar dynamic spectra are done.")
                print("###########################")
            except Exception:
                print("!!! WARNING : Error in making dynamic spectra. !!!")
                traceback.print_exc()
                if emails != "":
                    email_msg = "Error occured in making dynamic spectra."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )

        ##########################################
        # Basic calibration flows
        ##########################################
        if has_cal:
            cal_obsids = list(calibrator_dic.keys())
            futures = []
            for cal_obsid in cal_obsids:
                cal_datadir, cal_metafits, coarse_chans = calibrator_dic[cal_obsid]
                future = basic_cal_subflow.with_options(
                    flow_run_name=f"basic_cal_{jobid}",
                    task_runner=DaskTaskRunner(address=dask_addr),
                )(
                    cal_obsid,
                    cal_datadir,
                    cal_metafits,
                    coarse_chans,
                    target_obsid,
                    workdir,
                    cal_outdir,
                    basic_caldir,
                    do_basic_cal,
                    redo_basic_cal,
                    do_cal_flag,
                    do_import_model,
                    do_polcal,
                    keep_backup,
                    quack_timestamps,
                    cpu_frac,
                    mem_frac,
                    jobid,
                    #timestamp,
                    #emails,
                    remote_logger,
                )
                futures.append(future)
            results = [f.result() for f in futures]

        ###################################################
        # Checking if selfcal tables already exist or not
        ###################################################
        if not redo_selfcal:
            if target_obsid is not None:
                print("Checking pre-existing self-calibration solutions...")
                selfcal_gaincal = sorted(
                    glob.glob(f"{selfcaldir}/selfcal_{target_obsid}*.gcal")
                )
                if len(selfcal_gaincal) > 0:
                    selfcal_bandpass = sorted(
                        glob.glob(f"{selfcaldir}/selfcal_{target_obsid}*.bcal")
                    )
                    if len(selfcal_bandpass) > 0:
                        selfcal_bandpass = interpolate_bpass(
                            selfcal_bandpass, overwrite=True
                        )
                    if do_polcal:
                        selfcal_leakages = sorted(
                            glob.glob(f"{selfcaldir}/selfcal_{target_obsid}*.dcal")
                        )
                        if len(selfcal_leakages) > 0:
                            selfcal_leakages = interpolate_quartical(
                                selfcal_leakages, overwrite=True
                            )
                            do_selfcal = False
                            print(
                                "Self-calibration solutions exist including polarisation calibration. Not performing self-calibration"
                            )
                            if emails != "":
                                email_msg = "Self-calibration solutions including polarisation for target are already present."
                                send_task_notification(
                                    emails, email_msg, jobid, target_obsid, timestamp
                                )
                    else:
                        print(
                            "Self-calibration solutions exist without polarisation calibration. Henc, performing self-calibration"
                        )

        ###################################################
        # Start spliting selfcal ms
        ###################################################
        if do_selfcal:
            ###############################################
            # Removing previous self-calibration artificats
            ###############################################
            print("Removing all previous self-calibration artificats...")
            os.system(
                f"rm -rf {workdir}/selfcal* {workdir}/.intselfcal* {workdir}/.polselfcal*"
            )
            prefix = "selfcal"
            try:
                time_interval = float(solint)
            except BaseException:
                if solint.endswith("s"):
                    time_interval = float(solint.split("s")[0])
                elif solint.endswith("min"):
                    time_interval = float(solint.split("min")[0]) * 60
                elif solint == "int":
                    time_interval = image_timeres
                else:
                    time_interval = 60.0

            if adaptive:
                scale_worker_and_wait(
                    dask_cluster,
                    dask_client,
                    max(2, min(total_ncoarse + 1, max_worker)),
                )

            ######################
            # Spliting
            ######################
            if emails != "":
                email_msg = "Started spliting of measurement sets for self-calibration."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print(f"Starting task: Spliting {prefix} .....")
            print("###########################")
            ntime = get_selfcal_ntimes(target_mslist[0])
            msmd.open(target_mslist[0])
            times = msmd.timesforspws(0)
            timeres = np.nanmean(np.diff(times))
            msmd.close()
            time_window = min(10, round(ntime * timeres, 1))  # Maximum 10s
            future_selfcal_split = run_target_split_jobs.with_options(
                task_run_name=f"spliting_{prefix}_{jobid}"
            ).submit(
                ",".join(target_mslist),
                target_metafits,
                workdir,
                datacolumn="data",
                timeres=timeavg,
                freqres=freqavg,
                prefix=prefix,
                force_split=True,
                only_disk=True,
                time_window=min(time_window, time_interval),
                time_interval=time_interval,
                quack_timestamps=quack_timestamps,
                jobid=jobid,
                cpu_frac=float(cpu_frac),
                mem_frac=float(mem_frac),
                remote_log=remote_logger,
            )
            print("Checking status of spliting of target for selfcal ...")
            try:
                msg, expected, succeed = future_selfcal_split.result()
                if emails != "":
                    email_msg = f"Spliting of measurement sets for self-calibration is done.\nExpected: {expected}, succeeded: {succeed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(
                    "Finished task: Spliting of measurement sets for self-calibration is done."
                )
                print("###########################")
            except Exception:
                print(
                    "!!!! WARNING: Error in running spliting target scans for selfcal. !!!!"
                )
                do_selfcal = False
                traceback.print_exc()
                if emails != "":
                    email_msg = "Error occured in spliting target measurement sets for self-calibration."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
            finally:
                if adaptive:
                    scale_worker_and_wait(dask_cluster, dask_client, 2)

        ######################################
        # Checking status of self-cal split
        ######################################
        if do_selfcal:
            print("Checking measurement sets before spawning self-calibrations....")
            ####################################
            # Filtering any corrupted ms
            #####################################
            selfcal_target_mslist = sorted(glob.glob(workdir + "/selfcal*_spw_*.ms"))
            if (selfcal_target_mslist) == 0:
                print(
                    "!!!! WARNING: Error in running spliting target scans for selfcal. !!!!"
                )
                do_selfcal = False
                if emails != "":
                    email_msg = "No splited measurement set is found for self-calibration. Not continuting for self-calibration."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )

            filtered_mslist = []  # Filtering in case any ms is corrupted
            for ms in selfcal_target_mslist:
                checkcol = check_datacolumn_valid(ms)
                if checkcol:
                    filtered_mslist.append(ms)
                else:
                    print(f"Issue in : {ms}")
                    os.system(f"rm -rf {ms}")
            selfcal_mslist = filtered_mslist
            if len(selfcal_mslist) == 0:
                print(
                    "No splited target scan ms are available in work directory for selfcal. Not continuing further for selfcal."
                )
                do_selfcal = False
                if emails != "":
                    email_msg = "No splited measurement set is found for self-calibration. Not continuting for self-calibration."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
            print(f"Selfcal mslist : {[os.path.basename(i) for i in selfcal_mslist]}")

        #########################################################
        # Flagging on targets for self-calibration
        #########################################################
        cal_applied = False
        if do_selfcal:
            if adaptive:
                scale_worker_and_wait(
                    dask_cluster,
                    dask_client,
                    max(2, min(len(selfcal_mslist) + 1, max_worker)),
                )

            ###################################
            # Apply basic calibration
            ###################################
            if has_cal:
                if emails != "":
                    email_msg = "Started applying basic calibration solution on self-calibration measurement sets."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(
                    "Starting task: Applying basic calibration on self-calibration measurement sets....."
                )
                print("###########################")
                future_apply_basical_selfcal = run_apply_basiccal_sol.with_options(
                    task_run_name=f"applying_basiccal_selfcal_{jobid}"
                ).submit(
                    ",".join(selfcal_mslist),
                    calibrator_metafits,
                    target_metafits,
                    workdir,
                    basic_caldir,
                    overwrite_datacolumn=False,
                    only_amplitude=only_amplitude,
                    applymode="calflag",
                    prefix="selfcal",
                    jobid=jobid,
                    cpu_frac=round(cpu_frac, 2),
                    mem_frac=round(mem_frac, 2),
                    remote_log=remote_logger,
                )
                try:
                    msg, succeed, failed = future_apply_basical_selfcal.result()
                    cal_applied = True
                    if emails != "":
                        email_msg = f"Applying basic calibration solution on self-calibration measurement sets are done.\nSucceeded: {succeed}, failed: {failed}."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    print("###########################")
                    print(
                        "Finished task: Applying basic calibration solution on self-calibration measurement sets are done."
                    )
                    print("###########################")
                except Exception:
                    print(
                        "!!!! WARNING: Error in applying basic calibration solutions on target. Continuing selfcal without basic calibration.!!!!"
                    )
                    traceback.print_exc()
                    cal_applied = False
                    do_applycal = False
                    if emails != "":
                        email_msg = "Error occured in applying basic calibration solutions on self-calibration measurement sets."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )

            if cal_applied:
                selfcal_applymode = "calonly"
                filtered_selfcalms_list = []
                for selfcalms in selfcal_mslist:
                    unflag_chans, flag_chans = get_chans_flag(
                        msname=selfcalms, n_threads=n_threads
                    )
                    if len(flag_chans) / (len(flag_chans) + len(unflag_chans)) <= 0.8:
                        filtered_selfcalms_list.append(selfcalms)
                    else:
                        print(
                            f"More than 80% channels are flagged for ms: {selfcalms}. Not using for self-calibration."
                        )
                selfcal_mslist = filtered_selfcalms_list
            else:
                selfcal_applymode = "calflag"

            ###############################################
            # Performing sidereal correction before selfcal
            ###############################################
            os.system(
                f"rm -rf {workdir}/*selfcal_int* {workdir}/*selfcal_pol* {workdir}/caltables/*selfcal*"
            )
            if adaptive:
                scale_worker_and_wait(
                    dask_cluster,
                    dask_client,
                    max(2, min(len(selfcal_mslist) + 1, max_worker)),
                )
            if do_sidereal_cor:
                if emails != "":
                    email_msg = "Started correcting for solar sidereal motion."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(
                    "Starting task: Sidereal motion correction for self-calibration measurement sets....."
                )
                print("###########################")
                future_sidereal_cor_selfcal = run_solar_siderealcor_jobs.with_options(
                    task_run_name=f"solar_sidereal_correction_{jobid}"
                ).submit(
                    ",".join(selfcal_mslist),
                    workdir,
                    prefix="selfcal",
                    jobid=jobid,
                    cpu_frac=round(cpu_frac, 2),
                    mem_frac=round(mem_frac, 2),
                    remote_log=remote_logger,
                )
                try:
                    msg, succeed, failed = future_sidereal_cor_selfcal.result()
                    if emails != "":
                        email_msg = f"Correction for solar sidereal motion is done.\nSucceeded: {succeed}, failed: {failed}."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    print("###########################")
                    print(
                        "Finished task: Correction for solar sidereal motion is done."
                    )
                    print("###########################")
                except Exception:
                    print("Sidereal correction is not successful.")
                    traceback.print_exc()
                    if emails != "":
                        email_msg = "Error occured in sidereal motion correction."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )

            ############################
            # Basic flagging for selfcal
            ############################
            if emails != "":
                email_msg = "Started flagging for self-calibration measurment sets."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print("Starting task: Flagging selfcal targets .....")
            print("###########################")
            future_flag = run_flag.with_options(
                task_run_name=f"flagging_selfcal_{jobid}"
            ).submit(
                ",".join(selfcal_mslist),
                target_metafits,
                workdir,
                target_outdir,
                flag_calibrators=False,
                flag_quack=False,
                datacolumn="corrected",
                run_solarflagger=use_solarflagger,
                jobid=jobid,
                cpu_frac=round(cpu_frac, 2),
                mem_frac=round(mem_frac, 2),
                remote_log=remote_logger,
            )
            try:
                msg, succeed, failed = future_flag.result()
                if emails != "":
                    email_msg = f"Flagging for self-calibration measurment sets are done.\nSucceeded: {succeed}, failed: {failed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                for s_ms in selfcal_mslist:
                    s_ms = s_ms.rstrip("/")
                    if os.path.exists(f"{s_ms}/.flag_failed"):
                        print(
                            f"Issue in flagging of measurement set: {s_ms}. Check calibration solutions carefully."
                        )
                print("###########################")
                print(
                    "Finished task: Flagging for self-calibration measurment sets are done."
                )
                print("###########################")
            except Exception:
                print(
                    "!!!! WARNING: Flagging error. Examine calibration solutions with caution. !!!!"
                )
                traceback.print_exc()
                if emails != "":
                    email_msg = (
                        "Error occured in flagging self-calibration measurement sets."
                    )
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )

            #############################
            # Self-calibration
            #############################
            if emails != "":
                email_msg = "Started self-calibration."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print("Starting task: Self-calibrations.....")
            print("###########################")
            if cal_applied:
                print("Calibrator solutions are applied.")
            else:
                print("Calibration solutions are not applied")
            future_selfcal = run_selfcal_jobs.with_options(
                task_run_name=f"selfcal_{jobid}"
            ).submit(
                ",".join(selfcal_mslist),
                workdir,
                selfcaldir,
                target_metafits,
                cal_applied,
                solint=solint,
                do_apcal=do_ap_selfcal,
                do_polcal=do_polcal,
                solar_selfcal=solar_selfcal,
                keep_backup=keep_backup,
                uvrange=uvrange,
                weight="briggs",
                robust=0.0,
                applymode=selfcal_applymode,
                use_solarflagger=use_solarflagger,
                jobid=jobid,
                cpu_frac=round(cpu_frac, 2),
                mem_frac=round(mem_frac, 2),
                remote_log=remote_logger,
            )
            try:
                (
                    msg,
                    int_succeed,
                    int_failed,
                    pol_succeed,
                    pol_failed,
                    int_DR,
                    pol_DR,
                ) = future_selfcal.result()
                if emails != "":
                    email_msg = f"Self-calibration is done.\nIntensity self-calibration, Succeeded: {int_succeed}, failed: {int_failed}, average DR: {int_DR}."
                    if do_polcal:
                        email_msg += f"\nPolarisation self-calibration, Succeeded: {pol_succeed}, failed: {pol_failed}, average DR; {pol_DR}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print("Finished task: Self-calibration is done.")
                print("###########################")
            except Exception:
                print(
                    "!!!! WARNING: Error in self-calibration on targets. Not applying self-calibration. !!!!"
                )
                do_apply_selfcal = False
                traceback.print_exc()
                if emails != "":
                    email_msg = "Error occured in self-calibration."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
            if adaptive:
                scale_worker_and_wait(dask_cluster, dask_client, 2)

        ########################################
        # Checking self-cal caltables
        ########################################
        print(
            f"Searching for self-calibration gaincal tables: {selfcaldir}/selfcal_{target_obsid}*.gcal"
        )
        selfcal_gaincal = sorted(
            glob.glob(f"{selfcaldir}/selfcal_{target_obsid}*.gcal")
        )
        if len(selfcal_gaincal) == 0:
            print(
                "Self-calibration is not performed and no self-calibration caltable is available."
            )
            do_apply_selfcal = False
            if emails != "":
                email_msg = "Self-calibration is not performed and no self-calibration caltable is available."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
        else:
            print("###################################################")
            print(
                f"Self-calibration gaincal tables in calibration directory: {selfcaldir}"
            )
            for gcal in selfcal_gaincal:
                print(f"{os.path.basename(gcal)}")
            print("####################################################")

            print(
                f"Searching for self-calibration bandpass tables: {selfcaldir}/selfcal_{target_obsid}*.bcal"
            )
            selfcal_bandpass = sorted(
                glob.glob(f"{selfcaldir}/selfcal_{target_obsid}*.bcal")
            )
            if len(selfcal_bandpass) > 0:
                print("###################################################")
                print(
                    f"Self-calibration bandpass tables in calibration directory: {selfcaldir}"
                )
                for bpass in selfcal_bandpass:
                    print(f"{os.path.basename(bpass)}")
                print("####################################################")
                selfcal_bandpass = interpolate_bpass(selfcal_bandpass, overwrite=True)
            if do_polcal:
                print(
                    f"Searching for self-calibration polarisation leakage tables: {selfcaldir}/selfcal_{target_obsid}*.dcal"
                )
                selfcal_leakages = sorted(
                    glob.glob(f"{selfcaldir}/selfcal_{target_obsid}*.dcal")
                )
                if len(selfcal_leakages) > 0:
                    print("###################################################")
                    print(
                        f"Self-calibration polarisation leakage tables in calibration directory: {selfcaldir}"
                    )
                    for dcal in selfcal_leakages:
                        print(f"{os.path.basename(dcal)}")
                    print("####################################################")
                    selfcal_leakages = interpolate_quartical(
                        selfcal_leakages, overwrite=True
                    )

            ###########################################
            # Plotting self-caltables
            ###########################################
            if do_selfcal and len(selfcal_gaincal) > 0:
                os.makedirs(f"{target_outdir}/diagnostic_plots", exist_ok=True)
                msg, gcal_plots = plot_caltable_diagnostics(
                    selfcal_gaincal,
                    f"{target_outdir}/diagnostic_plots/{target_obsid}_gcal",
                )
                if msg == 0:
                    print(
                        f"Diagnostic plots for self-calibration gaincal tables are saved in : {gcal_plots}."
                    )
                else:
                    print(
                        "Error in creating diagnostic plots for self-calibration gaincal tables."
                    )

            if do_selfcal and len(selfcal_bandpass) > 0:
                os.makedirs(f"{target_outdir}/diagnostic_plots", exist_ok=True)
                msg, bcal_plots = plot_caltable_diagnostics(
                    selfcal_bandpass,
                    f"{target_outdir}/diagnostic_plots/{target_obsid}_bcal",
                )
                if msg == 0:
                    print(
                        f"Diagnostic plots for self-calibration bandpass tables are saved in : {bcal_plots}."
                    )
                else:
                    print(
                        "Error in creating diagnostic plots for self-calibration bandpass tables."
                    )

            if do_selfcal and do_polcal:
                if len(selfcal_leakages) > 0:
                    os.makedirs(f"{target_outdir}/diagnostic_plots", exist_ok=True)
                    msg, dcal_plots = plot_quartical_tables(
                        selfcal_leakages,
                        f"{target_outdir}/diagnostic_plots/{target_obsid}_dcal",
                    )
                    if msg == 0:
                        print(
                            f"Diagnostic plots for self-calibration leakage tables are saved in : {dcal_plots}."
                        )
                    else:
                        print(
                            "Error in creating diagnostic plots for self-calibration leakage tables."
                        )

        #############################################
        # Spliting targets if not started already
        #############################################
        # If corrected data is requested or imaging is requested

        if do_applycal or do_apply_selfcal or do_imaging:
            if adaptive:
                scale_worker_and_wait(
                    dask_cluster,
                    dask_client,
                    max(2, min(total_ncoarse + 1, max_worker)),
                )
            prefix = "target"
            if emails != "":
                email_msg = "Started spliting target for final processing."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print(f"Starting task: Spliting {prefix} .....")
            print("###########################")
            future_split = run_target_split_jobs.with_options(
                task_run_name=f"spliting_{prefix}_{jobid}"
            ).submit(
                ",".join(target_mslist),
                target_metafits,
                workdir,
                datacolumn="data",
                force_split=True,
                freqres=freqavg,
                timeres=timeavg,
                quack_timestamps=quack_timestamps,
                prefix=prefix,
                jobid=jobid,
                cpu_frac=round(cpu_frac, 2),
                mem_frac=round(mem_frac, 2),
                remote_log=remote_logger,
            )
            ##########################################
            # Checking target spliting is done or not
            ##########################################
            print("Checking spliting of targets status...")
            try:
                msg, expected, succeed = future_split.result()
                if emails != "":
                    email_msg = f"Spliting target for final processing is done.\nExpected: {expected}, succeeded: {succeed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print("Finished task: Spliting target for final processing is done.")
                print("###########################")
            except Exception:
                print("!!!! WARNING: Error in spliting targets. !!!!")
                traceback.print_exc()
                if emails != "":
                    email_msg = "Error occured in spliting target for final processing."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                return 1
            finally:
                if adaptive:
                    scale_worker_and_wait(dask_cluster, dask_client, 2)

            ################################
            # Checking splited final ms list
            ################################
            split_target_mslist = sorted(glob.glob(workdir + "/target*_spw_*.ms"))
            if len(split_target_mslist) == 0:
                print("!!!! WARNING: No target ms are present. !!!!")
                if emails != "":
                    email_msg = (
                        "No target measurement set is present for final processing."
                    )
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                return 1
            ####################################
            # Filtering any corrupted ms
            #####################################
            print(
                "Checking final valid measurement sets before applying solutions and spawning imaging...."
            )
            filtered_mslist = []  # Filtering in case any ms is corrupted
            for ms in split_target_mslist:
                checkcol = check_datacolumn_valid(ms)
                if checkcol:
                    filtered_mslist.append(ms)
                else:
                    print(f"Issue in : {ms}")
                    os.system(f"rm -rf {ms}")
            split_target_mslist = filtered_mslist
            if len(split_target_mslist) == 0:
                print("No filtered target ms are available in work directory.")
                if emails != "":
                    email_msg = "No un-corrupted target measurement is present for final processing."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                return 1
            print(
                f"Target mslist : {[os.path.basename(i) for i in split_target_mslist]}"
            )

            #########################################################
            # Applying basic solutions on target scans
            #########################################################
            if adaptive:
                scale_worker_and_wait(
                    dask_cluster,
                    dask_client,
                    max(2, min(len(split_target_mslist) + 1, max_worker)),
                )

            ####################################
            # Applying basic calibration
            #####################################
            if (do_applycal or do_apply_selfcal) and has_cal:
                if emails != "":
                    email_msg = "Started applying basic calibration solutions on final target measurement sets."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(
                    "Starting task: Applying basic calibration on final target measurement sets....."
                )
                print("###########################")
                future_apply_basical = run_apply_basiccal_sol.with_options(
                    task_run_name=f"applying_basiccal_target_{jobid}"
                ).submit(
                    ",".join(split_target_mslist),
                    calibrator_metafits,
                    target_metafits,
                    workdir,
                    basic_caldir,
                    overwrite_datacolumn=True,
                    only_amplitude=only_amplitude,
                    applymode="calflag",
                    prefix="target",
                    jobid=jobid,
                    cpu_frac=round(cpu_frac, 2),
                    mem_frac=round(mem_frac, 2),
                    remote_log=remote_logger,
                )
                try:
                    msg, succeed, failed = future_apply_basical.result()
                    if emails != "":
                        email_msg = f"Applying basic calibration solutions on final target measurement sets are done.\nSucceeded: {succeed}, failed: {failed}."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    print("###########################")
                    print(
                        "Finished task: Applying basic calibration solutions on final target measurement sets are done."
                    )
                    print("###########################")
                except Exception:
                    print(
                        "!!!! WARNING: Error in applying basic calibration solutions on target scans. Not continuing further.!!!!"
                    )
                    traceback.print_exc()
                    if emails != "":
                        email_msg = "Error occured in applying basic calibration on final target measurement sets. P-AIRCARS has stopped."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    return 1

            ###################################
            # Correct sidereal motion
            ###################################
            if do_sidereal_cor:
                if emails != "":
                    email_msg = "Start correcting sidereal motion of the Sun on final target measurement sets."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(
                    "Starting task: Sidereal motion correction for final target measurement sets....."
                )
                print("###########################")
                future_sidereal_cor = run_solar_siderealcor_jobs.with_options(
                    task_run_name=f"solar_sidereal_correction_{jobid}"
                ).submit(
                    ",".join(split_target_mslist),
                    workdir,
                    prefix="target",
                    jobid=jobid,
                    cpu_frac=round(cpu_frac, 2),
                    mem_frac=round(mem_frac, 2),
                    remote_log=remote_logger,
                )
                try:
                    msg, succeed, failed = future_sidereal_cor.result()
                    if emails != "":
                        email_msg = f"Sidereal motion correction of the Sun on final target measurement sets are done.\nSucceeded: {succeed}, failed: {failed}."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    print("###########################")
                    print(
                        "Finished task: Sidereal motion correction of the Sun on final target measurement sets are done."
                    )
                    print("###########################")
                except Exception:
                    print("!!!! WARNING: Error in applying sidereal correction.!!!!")
                    traceback.print_exc()
                    if emails != "":
                        email_msg = "Error occured in sidereal motion correction on final target measurement sets."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )

            ########################################
            # Apply self-calibration
            ########################################
            if do_apply_selfcal:
                selfcal_applymode = "calonly"
                for msname in split_target_mslist:
                    if not os.path.exists(f"{msname}/.applied_sol"):
                        selfcal_applymode = "calflag"

                if emails != "":
                    email_msg = "Started applying self-calibration on final target measurement sets."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(
                    "Starting task: Applying self-calibration solutions on final target measurement sets....."
                )
                print("###########################")
                future_apply_selfcal = run_apply_selfcal_sol.with_options(
                    task_run_name=f"applying_selfcal_{jobid}"
                ).submit(
                    ",".join(split_target_mslist),
                    target_metafits,
                    workdir,
                    selfcaldir,
                    overwrite_datacolumn=False,
                    applymode=selfcal_applymode,
                    jobid=jobid,
                    cpu_frac=round(cpu_frac, 2),
                    mem_frac=round(mem_frac, 2),
                    remote_log=remote_logger,
                )
                try:
                    msg, gain_succeed, gain_failed, pol_succeed, pol_failed = (
                        future_apply_selfcal.result()
                    )
                    if emails != "":
                        email_msg = f"Applying self-calibration on final target measurement sets are done.\nGain solutions applied: Succeeded: {gain_succeed}, failed: {gain_failed}."
                        if do_polcal:
                            email_msg += f"\nPolarisation solution applied: Succeeded: {pol_succeed}, failed: {pol_failed}."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    print("###########################")
                    print(
                        "Finished task: Applying self-calibration on final target measurement sets are done."
                    )
                    print("###########################")
                except Exception:
                    print(
                        "!!!! WARNING: Error in applying self-calibration solutions on targets. !!!!"
                    )
                    traceback.print_exc()
                    if emails != "":
                        email_msg = "Error occured in applying self-calibration solutions on final target measurement sets."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )

            ############################
            # Basic flagging
            ############################
            if emails != "":
                email_msg = "Started flagging of final target measurement sets."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print("Starting task: Flagging final target measurement sets .....")
            print("###########################")
            if not use_solarflagger:
                dr_files = glob.glob(f"{selfcaldir}/selfcal_{target_obsid}*.DR")
                if len(dr_files) > 0:
                    int_DR_list = []
                    pol_DR_list = []
                    for dr_file in dr_files:
                        int_DR, pol_DR = np.load(dr_file, allow_pickle=True)
                        int_DR_list.append(int_DR)
                        pol_DR_list.append(pol_DR)
                    avg_int_DR = np.nanmedian(int_DR_list)
                    avg_pol_DR = np.nanmedian(pol_DR_list)
                    if avg_int_DR < 100 or avg_pol_DR < 100:
                        print(
                            f"Average intensity self-calibration dynamic range: {avg_int_DR} is smaller than 100."
                        )
                        print(
                            f"Average polarisation self-calibration dynamic range: {avg_pol_DR} is smaller than 100."
                        )
                        print("Using solar flagger.")
                        use_solarflagger = True

            future_flag = run_flag.with_options(
                task_run_name=f"flagging_target_{jobid}"
            ).submit(
                ",".join(split_target_mslist),
                target_metafits,
                workdir,
                target_outdir,
                flag_calibrators=False,
                flag_quack=False,
                datacolumn="corrected",
                run_solarflagger=use_solarflagger,
                jobid=jobid,
                cpu_frac=round(cpu_frac, 2),
                mem_frac=round(mem_frac, 2),
                remote_log=remote_logger,
            )
            try:
                msg, succeed, failed = future_flag.result()
                if emails != "":
                    email_msg = f"Flagging of final target measurement sets are done.\nSucceeded: {succeed}, failed: {failed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(
                    "Finished task: Flagging of final target measurement sets are done."
                )
                print("###########################")
            except Exception:
                print(
                    "!!!! WARNING: Flagging error. Examine calibration solutions with caution. !!!!"
                )
                traceback.print_exc()
                if emails != "":
                    email_msg = (
                        "Error occured in flagging of final target measurement sets."
                    )
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )

            ######################################
            # Imaging
            ######################################
            if do_imaging:
                if image_freqres > 0:
                    print(f"Image frequency resolution: {image_freqres} MHz.")
                else:
                    print("Image frequency resolution: entire corase channel.")
                if image_timeres > 0:
                    print(f"Image time resolution: {image_timeres} s.")
                else:
                    print("Imaging entire scan.")
                pol = pol.upper()
                if pol not in ["I", "IQUV"]:
                    pol = "IQUV"

                if (
                    not do_polcal
                ):  # Only if do_polcal is False, overwrite to make only Stokes I
                    pol = "I"

                if emails != "":
                    email_msg = "Started final imaging."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print("Starting task: Final imaging.....")
                print("###########################")
                future_imaging = run_imaging_jobs.with_options(
                    task_run_name=f"imaging_{jobid}"
                ).submit(
                    ",".join(split_target_mslist),
                    workdir,
                    target_outdir,
                    freqrange=freqrange,
                    timerange=timerange,
                    minuv=minuv,
                    weight=weight,
                    robust=float(robust),
                    pol=pol,
                    freqres=image_freqres,
                    timeres=image_timeres,
                    threshold=float(clean_threshold),
                    use_multiscale=use_multiscale,
                    use_solar_mask=use_solar_mask,
                    cutout_rsun=cutout_rsun,
                    savemodel=keep_backup,
                    saveres=keep_backup,
                    jobid=jobid,
                    cpu_frac=round(cpu_frac, 2),
                    mem_frac=round(mem_frac, 2),
                    remote_log=remote_logger,
                )
                try:
                    msg, succeed, failed, total_images = future_imaging.result()
                    if emails != "":
                        email_msg = f"Final imaging is done.\nSucceeded: {succeed}, failed: {failed}.\nTotal images made: {total_images}."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    print("###########################")
                    print("Finished task: Final imaging is done.")
                    print("###########################")
                except Exception:
                    print(
                        "!!!! WARNING: Final imaging on all measurement sets is not successful. Check the image directory. !!!!"
                    )
                    traceback.print_exc()
                    if emails != "":
                        email_msg = (
                            "Error occured in final imaging. P-AIRCARS has stopped."
                        )
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    return 1
            if adaptive:
                scale_worker_and_wait(dask_cluster, dask_client, 2)

        ########################################
        # Naming of image directory
        ########################################
        if weight == "briggs":
            weight_str = f"{weight}_{robust}"
        else:
            weight_str = weight
        if image_freqres == -1 and image_timeres == -1:
            imagedir = target_outdir + f"/imagedir_f_all_t_all_pol_{pol}_w_{weight_str}"
        elif image_freqres != -1 and image_timeres == -1:
            imagedir = (
                target_outdir
                + f"/imagedir_f_{image_freqres}_t_all_pol_{pol}_w_{weight_str}"
            )
        elif image_freqres == -1 and image_timeres != -1:
            imagedir = (
                target_outdir
                + f"/imagedir_f_all_t_{image_timeres}_pol_{pol}_w_{weight_str}"
            )
        else:
            imagedir = (
                target_outdir
                + f"/imagedir_f_{image_freqres}_t_{image_timeres}_pol_{pol}_w_{weight_str}"
            )

        os.system(f"rm -rf {imagedir}/images/*aia*.fits")
        os.system(f"rm -rf {imagedir}/images/*suvi*.fits")
        ###########################
        # Primary beam correction
        ###########################
        if do_pbcor:
            images = sorted(glob.glob(f"{imagedir}/images/*.fits"))
            if len(images) == 0:
                print(f"No image is present in image directory: {imagedir}/images")
                if emails != "":
                    email_msg = "No image is present in image directory for primary beam correction."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
            else:
                if adaptive:
                    scale_worker_and_wait(
                        dask_cluster,
                        dask_client,
                        max(2, min(len(images) + 1, max_worker)),
                    )
                if emails != "":
                    email_msg = "Started primary beam correction."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print("Starting task: Primary beam correction.....")
                print("###########################")
                future_pbcor = run_apply_pbcor.with_options(
                    task_run_name=f"applying_primary_beam_{jobid}"
                ).submit(
                    f"{imagedir}/images",
                    target_metafits,
                    workdir,
                    leakage_dir=selfcaldir,
                    jobid=jobid,
                    cpu_frac=round(cpu_frac, 2),
                    mem_frac=round(mem_frac, 2),
                    remote_log=remote_logger,
                )
                try:
                    msg, succeed, failed = future_pbcor.result()
                    if emails != "":
                        email_msg = f"Primary beam correction is done.\nSucceeded: {succeed}, failed: {failed}."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    print("###########################")
                    print("Finished task: Primary beam correction is done.")
                    print(f"Final image directory: {imagedir}/images")
                    print("###########################")
                except Exception:
                    print(
                        "!!!! WARNING: Primary beam corrections of the final images are not successful. !!!!"
                    )
                    traceback.print_exc()
                    if emails != "":
                        email_msg = "Error occured in primary beam correction. P-AIRCARS has stopped."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    return 1
                finally:
                    if adaptive:
                        scale_worker_and_wait(dask_cluster, dask_client, 2)

        ##############################################
        # Making diagnostic plots of measurement sets
        ##############################################
        if make_msplot:
            ###########################################
            # Ploting calibrator ms
            ###########################################
            split_cal_mslist = sorted(glob.glob(f"{workdir}/calibrator*_spw_*.ms"))
            if len(split_cal_mslist) == 0:
                print("No calibrator measurement set is present for ploting.")
            else:
                if adaptive:
                    scale_worker_and_wait(
                        dask_cluster,
                        dask_client,
                        max(2, min(len(split_cal_mslist) + 1, max_worker)),
                    )
                msplot_outdir = f"{cal_outdir}/ms_diagnostics_plots"
                os.makedirs(msplot_outdir, exist_ok=True)
                if len(split_cal_mslist) > 0:
                    if emails != "":
                        email_msg = "Started making diagnostic plots for calibrator measurement sets."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    print("###########################")
                    print(
                        "Starting task: Making diagnostic plots of calibrator measurement sets....."
                    )
                    print("###########################")
                    future_cal_plot = run_make_msplot.with_options(
                        task_run_name=f"making_msplot_cal_{jobid}"
                    ).submit(
                        ",".join(split_cal_mslist),
                        workdir,
                        msplot_outdir,
                        jobid=jobid,
                        cpu_frac=round(cpu_frac, 2),
                        mem_frac=round(mem_frac, 2),
                        remote_log=remote_logger,
                    )
                    try:
                        msg = future_cal_plot.result()
                        if emails != "":
                            email_msg = "Making diagnostic plots for calibrator measurement sets are done."
                            send_task_notification(
                                emails, email_msg, jobid, target_obsid, timestamp
                            )
                        print("###########################")
                        print(
                            "Finished task: Making diagnostic plots for calibrator measurment sets are done."
                        )
                        print("###########################")
                    except Exception:
                        print(
                            "!!!! WARNING: Diagnostic plot of calibrator measurment sets are not successful. !!!!"
                        )
                        traceback.print_exc()
                        if emails != "":
                            email_msg = "Error occured in making diagnostic plots of calibrator measurement sets."
                            send_task_notification(
                                emails, email_msg, jobid, target_obsid, timestamp
                            )

            ###########################################
            # Ploting target ms
            ###########################################
            split_target_mslist = sorted(glob.glob(f"{workdir}/target*_spw_*.ms"))
            if len(split_target_mslist) == 0:
                print("No target measurment set is present for ploting.")
            else:
                if adaptive:
                    scale_worker_and_wait(
                        dask_cluster,
                        dask_client,
                        max(2, min(len(split_target_mslist) + 1, max_worker)),
                    )
                msplot_outdir = f"{target_outdir}/ms_diagnostics_plots"
                os.makedirs(msplot_outdir, exist_ok=True)
                print("###########################")
                print(
                    "Starting task: Making diagnostic plots of target measurement sets....."
                )
                print("###########################")
                future_target_plot = run_make_msplot.with_options(
                    task_run_name=f"making_msplot_{jobid}"
                ).submit(
                    ",".join(split_target_mslist),
                    workdir,
                    msplot_outdir,
                    jobid=jobid,
                    cpu_frac=round(cpu_frac, 2),
                    mem_frac=round(mem_frac, 2),
                    remote_log=remote_logger,
                )
                try:
                    msg = future_target_plot.result()
                    if emails != "":
                        email_msg = "Making diagnostic plots for target measurement sets are done."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    print("###########################")
                    print(
                        "Finished task: Making diagnostic plots for target measurment sets are done."
                    )
                    print("###########################")
                except Exception:
                    print(
                        "!!!! WARNING: Diagnostic plot of target measurment sets are not successful. !!!!"
                    )
                    traceback.print_exc()
                    if emails != "":
                        email_msg = "Error occured in making diagnostic plots of target measurement sets."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
        if adaptive:
            scale_worker_and_wait(dask_cluster, dask_client, 2)
        #######################################
        # Make overlays
        #######################################
        images = sorted(glob.glob(f"{imagedir}/images/*.fits"))
        if len(images) == 0:
            print(
                f"No image is present in image directory: {imagedir}/images for making overlays"
            )
            if emails != "":
                email_msg = (
                    "No image is present in image directory for making overlays."
                )
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )

        #################################################################
        # Filtering only coarse channel images for default overlay mode
        #################################################################
        if make_overlay is False and len(images) > 0:
            images = filter_images(images, min_time_sep=60.0)
        internet_on = internet_available()
        if not internet_on:
            print("Internet connection is not available. Can not make overlays")
        elif len(images) > 0:
            if adaptive:
                scale_worker_and_wait(
                    dask_cluster, dask_client, max(2, min(len(images) + 1, max_worker))
                )
            #################################
            # Start overlays
            #################################
            if emails != "":
                email_msg = "Started making overlays."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print("Starting task: Making overlay on EUV images.....")
            print("###########################")
            future_overlay = run_make_overlay.with_options(
                task_run_name=f"making_overlay_{jobid}"
            ).submit(
                f"{imagedir}/images",
                f"{imagedir}/overlay_pngs",
                workdir=workdir,
                all_overlay=make_overlay,
                jobid=jobid,
                cpu_frac=round(cpu_frac, 2),
                remote_log=remote_logger,
            )
            try:
                msg, succeed, failed = future_overlay.result()
                if msg == 0:
                    if emails != "":
                        email_msg = f"Making overlays are done.\nSucceeded: {succeed}, failed: {failed}."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    print("###########################")
                    print("Finished task: Making overlays are done.")
                    print(f"Final image directory: {imagedir}/overlay_pngs")
                    print("###########################")
                else:
                    if emails != "":
                        email_msg = f"Making overlays are not successful EUV images could not be download.\nSucceeded: {succeed}, failed: {failed}."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    print("###########################")
                    print("Finished task: Making overlays are not successful.")
                    if len(glob.glob(f"{imagedir}/overlay_pngs/*.png")) == 0:
                        os.system(f"rm -rf {imagedir}/overlay_pngs")
                    else:
                        print(f"Final image directory: {imagedir}/overlay_pngs")
                    print("###########################")
            except Exception:
                print("!!!! WARNING: Overlay of the images are not successful. !!!!")
                traceback.print_exc()
                if emails != "":
                    email_msg = "Error occured in making overlays."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )

        ###########################################
        # Successful exit
        ###########################################
        print("P-AIRCARS calibration and imaging pipeline is successfully executed.")
        if emails != "":
            email_msg = "P-AIRCARS processing is done successfully."
            send_task_notification(emails, email_msg, jobid, target_obsid, timestamp)
        return 0
    except Exception as e:
        traceback.print_exc()
        if emails != "":
            email_msg = f"Error in running P-AIRCARS.\n{e}"
            send_task_notification(emails, email_msg, jobid, target_obsid, timestamp)
        print("###########################")
        print("Error occured in running P-AIRCARS.")
        print("###########################")
        return 1
    finally:
        time.sleep(5)
        datalist = sorted(glob.glob(f"{target_datadir}/*"))
        for data in datalist:
            drop_cache(data)
        callist = sorted(glob.glob(f"{calibrator_datadir}/*"))
        for cal in callist:
            drop_cache(cal)
        ######################################
        # Keeping flag backups
        ######################################
        # Flag backups of calibrator measurement sets
        ######################################
        final_cal_mslist = sorted(glob.glob(workdir + "/calibrator*_spw_*.ms"))
        if len(final_cal_mslist) > 0:
            os.makedirs(f"{cal_outdir}/ms_flags", exist_ok=True)
            print(
                f"Doing flag backup for calibrator measurement sets in: {cal_outdir}/ms_flags"
            )
            for cal_ms in final_cal_mslist:
                do_flag_backup(cal_ms, flagtype="finalflag")
                if os.path.exists(
                    f"{cal_outdir}/ms_flags/{os.path.basename(cal_ms)}.flagversions"
                ):
                    os.system(
                        f"rm -rf {cal_outdir}/ms_flags/{os.path.basename(cal_ms)}.flagversions"
                    )
                os.system(f"mv {cal_ms}.flagversions {cal_outdir}/ms_flags/")
                if keep_backup is False:
                    os.system(f"rm -rf {cal_ms}")
        ######################################
        # Flag backups of selfcal measurement sets
        ######################################
        final_selfcal_mslist = sorted(glob.glob(workdir + "/selfcal*_spw_*.ms"))
        if len(final_selfcal_mslist) > 0:
            os.makedirs(f"{target_outdir}/ms_flags", exist_ok=True)
            print(
                f"Doing flag backup of self-calibration measurement sets in: {target_outdir}/ms_flags"
            )
            for selfcal_ms in final_selfcal_mslist:
                do_flag_backup(selfcal_ms, flagtype="finalflag")
                if os.path.exists(
                    f"{target_outdir}/ms_flags/{os.path.basename(selfcal_ms)}.flagversions"
                ):
                    os.system(
                        f"rm -rf {target_outdir}/ms_flags/{os.path.basename(selfcal_ms)}.flagversions"
                    )
                os.system(f"mv {selfcal_ms}.flagversions {target_outdir}/ms_flags/")
                if keep_backup is False:
                    os.system(f"rm -rf {selfcal_ms}")
        ######################################
        # Flag backups of target measurement sets
        ######################################
        final_split_target_mslist = sorted(glob.glob(workdir + "/target*_spw_*.ms"))
        if len(final_split_target_mslist) > 0:
            os.makedirs(f"{target_outdir}/ms_flags", exist_ok=True)
            print(
                f"Doing flag backup target measurement sets in: {target_outdir}/ms_flags"
            )
            for target_ms in final_split_target_mslist:
                do_flag_backup(target_ms, flagtype="finalflag")
                if os.path.exists(
                    f"{target_outdir}/ms_flags/{os.path.basename(target_ms)}.flagversions"
                ):
                    os.system(
                        f"rm -rf {target_outdir}/ms_flags/{os.path.basename(target_ms)}.flagversions"
                    )
                os.system(f"mv {target_ms}.flagversions {target_outdir}/ms_flags/")
                if keep_calibrated_ms:
                    calibrated_msdir = f"{target_outdir}/calibrated_ms"
                    os.makedirs(calibrated_msdir, exist_ok=True)
                    os.system(f"mv {target_ms} {calibrated_msdir}")
                elif keep_backup is False:
                    os.system(f"rm -rf {target_ms}")
        time.sleep(5)
        drop_cache(workdir)
        drop_cache(outdir)
        if master_log_created:
            stop_event.set()
            log_thread_flow.join(timeout=5)
        if dask_dir is not None:
            os.system(f"rm -rf {dask_dir}")
        if observer is not None:
            clean_shutdown(observer)


def cli():
    parser = argparse.ArgumentParser(
        description="Run P-AIRCARS for calibration and imaging of solar observations.",
        formatter_class=SmartDefaultsHelpFormatter,
    )
    # === Essential parameters ===
    essential = parser.add_argument_group(
        "###################\nEssential parameters\n###################"
    )
    essential.add_argument(
        "target_datadir", type=str, help="Target measurement set directory"
    )
    essential.add_argument(
        "--workdir",
        type=str,
        dest="workdir",
        required=True,
        help="Working directory",
    )
    essential.add_argument(
        "--outdir",
        type=str,
        dest="outdir",
        required=True,
        help="Output products directory",
    )
    essential.add_argument(
        "--target_metafits",
        type=str,
        dest="target_metafits",
        help="Target metafits file",
    )
    essential.add_argument(
        "--cal_datadir",
        type=str,
        dest="cal_datadir",
        help="Calibrator measurement set directory",
    )
    essential.add_argument(
        "--cal_metafits",
        type=str,
        dest="cal_metafits",
        help="Calibrator metafits file",
    )

    # === Advanced calibration parameters ===
    advanced_cal = parser.add_argument_group(
        "###################\nAdvanced calibration parameters\n###################"
    )
    advanced_cal.add_argument(
        "--solint",
        type=str,
        default="60s",
        help="Solution interval for calibration (e.g. 'int', '10s', '5min', 'inf')",
    )
    advanced_cal.add_argument(
        "--cal_uvrange",
        type=str,
        default="",
        help="UV range to filter data for calibration (e.g. '>100klambda', '100~10000lambda')",
    )
    advanced_cal.add_argument(
        "--no_polcal",
        action="store_false",
        dest="do_polcal",
        help="Disable polarization calibration",
    )
    advanced_cal.add_argument(
        "--only_amplitude",
        action="store_true",
        help="Apply only amplitude part of gain solution from calibrator or not",
    )
    advanced_cal.add_argument(
        "--redo_basic_cal",
        action="store_true",
        help="Redo basic calibration or not",
    )
    advanced_cal.add_argument(
        "--redo_selfcal",
        action="store_true",
        help="Redo self-calibration",
    )
    advanced_cal.add_argument(
        "--use_solarflagger",
        action="store_true",
        help="Use solar flagger",
    )

    # === Advanced imaging parameters ===
    advanced_image = parser.add_argument_group(
        "###################\nAdvanced imaging parameters\n###################"
    )
    advanced_image.add_argument(
        "--freqrange",
        type=str,
        default="",
        help="Frequency range in MHz to select during imaging (comma-seperate, e.g. '100~110,130~140')",
    )
    advanced_image.add_argument(
        "--timerange",
        type=str,
        default="",
        help="Time range to select during imaging (comma-seperated, e.g. '2014/09/06/09:30:00~2014/09/06/09:45:00,2014/09/06/10:30:00~2014/09/06/10:45:00')",
    )
    advanced_image.add_argument(
        "--image_freqres",
        type=float,
        default=1.28,
        help="Output image frequency resolution in MHz (-1 = full)",
    )
    advanced_image.add_argument(
        "--image_timeres",
        type=float,
        default=10.0,
        help="Output image time resolution in seconds (-1 = full)",
    )
    advanced_image.add_argument(
        "--pol",
        type=str,
        default="IQUV",
        help="Stokes parameter(s) to image ('I' or 'IQUV')",
    )
    advanced_image.add_argument(
        "--minuv",
        type=float,
        default=0,
        help="Minimum baseline length (in wavelengths) to include in imaging",
    )
    advanced_image.add_argument(
        "--weight",
        type=str,
        default="briggs",
        help="Imaging weighting scheme (e.g. 'briggs', 'natural', 'uniform')",
    )
    advanced_image.add_argument(
        "--robust",
        type=float,
        default=0.0,
        help="Robust parameter for Briggs weighting (-2 to +2)",
    )
    advanced_image.add_argument(
        "--no_multiscale",
        action="store_false",
        dest="use_multiscale",
        help="Disable multiscale CLEAN for extended structures",
    )
    advanced_image.add_argument(
        "--clean_threshold",
        type=float,
        default=1.0,
        help="Clean threshold in sigma for final deconvolution",
    )
    advanced_image.add_argument(
        "--no_pbcor",
        action="store_false",
        dest="do_pbcor",
        help="Do not apply primary beam correction after imaging",
    )
    advanced_image.add_argument(
        "--cutout_rsun",
        type=float,
        default=10.0,
        help="Field of view cutout radius in solar radii",
    )
    advanced_image.add_argument(
        "--no_solar_mask",
        action="store_false",
        dest="use_solar_mask",
        help="Disable use solar disk mask during deconvolution",
    )
    advanced_image.add_argument(
        "--do_overlay",
        action="store_true",
        dest="make_overlay",
        help="Make overlay plot on EUV images for all images (default is to make overlays only one image per coarse channels at 10s intervals)",
    )
    advanced_image.add_argument(
        "--make_msplot",
        action="store_true",
        help="Make diagnostic plots of measurement sets",
    )

    # === Advanced options ===
    advanced = parser.add_argument_group(
        "###################\nAdvanced pipeline parameters\n###################"
    )
    advanced.add_argument(
        "--non_solar_data",
        action="store_false",
        dest="solar_data",
        help="Disable solar data mode",
    )
    advanced.add_argument(
        "--no_ds",
        action="store_false",
        dest="make_ds",
        help="Disable making solar dynamic spectra",
    )
    advanced.add_argument(
        "--do_forcereset_weightflag",
        action="store_true",
        help="Force reset of weights and flags (disabled by default)",
    )
    advanced.add_argument(
        "--no_cal_flag",
        action="store_false",
        dest="do_cal_flag",
        help="Disable initial flagging of calibrators",
    )
    advanced.add_argument(
        "--no_import_model",
        action="store_false",
        dest="do_import_model",
        help="Disable model import",
    )
    advanced.add_argument(
        "--no_basic_cal",
        action="store_false",
        dest="do_basic_cal",
        help="Disable basic gain calibration",
    )
    advanced.add_argument(
        "--do_sidereal_cor",
        action="store_true",
        dest="do_sidereal_cor",
        help="Sidereal motion correction for Sun (disabled by default)",
    )
    advanced.add_argument(
        "--no_solarcenter_move",
        action="store_false",
        dest="do_move_solarcenter",
        help="Disable moving phaseceneter to solar center",
    )
    advanced.add_argument(
        "--no_selfcal",
        action="store_false",
        dest="do_selfcal",
        help="Disable self-calibration",
    )
    advanced.add_argument(
        "--no_ap_selfcal",
        action="store_false",
        dest="do_ap_selfcal",
        help="Disable amplitude-phase self-calibration",
    )
    advanced.add_argument(
        "--no_solar_selfcal",
        action="store_false",
        dest="solar_selfcal",
        help="Disable solar-specific self-calibration parameters",
    )
    advanced.add_argument(
        "--no_applycal",
        action="store_false",
        dest="do_applycal",
        help="Disable application of basic calibration solutions",
    )
    advanced.add_argument(
        "--no_apply_selfcal",
        action="store_false",
        dest="do_apply_selfcal",
        help="Disable application of self-calibration solutions",
    )
    advanced.add_argument(
        "--no_imaging",
        action="store_false",
        dest="do_imaging",
        help="Disable final imaging",
    )

    # === Advanced local system/ per node hardware resource parameters ===
    advanced_resource = parser.add_argument_group(
        "###################\nAdvanced hardware resource parameters for local system or per node on HPC cluster\n###################"
    )
    advanced_resource.add_argument(
        "--cpu_frac",
        type=float,
        default=0.8,
        help="Fraction of CPU usuage per node",
    )
    advanced_resource.add_argument(
        "--mem_frac",
        type=float,
        default=0.8,
        help="Fraction of memory usuage per node",
    )
    advanced_resource.add_argument(
        "--max_worker",
        type=int,
        default=None,
        help="Maximum number of workers",
    )
    advanced_resource.add_argument(
        "--keep_backup",
        action="store_true",
        help="Keep backup of intermediate steps",
    )
    advanced_resource.add_argument(
        "--no_calibrated_ms",
        action="store_false",
        dest="keep_calibrated_ms",
        help="Keep calibrated measurement sets or not",
    )
    advanced_resource.add_argument(
        "--no_remote_logger",
        action="store_false",
        dest="remote_logger",
        help="Disable remote logger",
    )
    advanced_resource.add_argument(
        "--jobid",
        type=int,
        default=None,
        help="User provided P-AIRCARS job ID",
    )
    advanced_resource.add_argument(
        "--job_password",
        type=str,
        default=None,
        help="User specified job password",
    )
    advanced_resource.add_argument(
        "--masterlog",
        type=str,
        default=None,
        help="Master log file",
    )
    advanced_resource.add_argument(
        "--cluster",
        action="store_true",
        dest="cluster",
        help="Running in cluster environment",
    )
    advanced_resource.add_argument(
        "--adaptive",
        action="store_true",
        dest="adaptive",
        help="Whether do adaptive scaling or not in cluster environment (In local environment, it is always true)",
    )

    # === Advanced job scheduler parameters ===
    advanced_slurm = parser.add_argument_group(
        "###################\nAdvanced slurm cluster settings\n###################"
    )
    advanced_slurm.add_argument(
        "--partition",
        type=str,
        default=None,
        help="Partition name (Required)",
    )
    advanced_slurm.add_argument(
        "--account",
        type=str,
        default=None,
        help="Account name (If your cluster requires this, you should provide. Otherwise job can not be started)",
    )
    advanced_slurm.add_argument(
        "--walltime",
        type=str,
        default=None,
        help="Wall time, each slurm job can execute in maximum this time",
    )
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    scheduler_name = get_scheduler_name()
    ##################################################
    # Prefect settings check
    ##################################################
    prefect_settings = get_current_settings()
    prefect_env = prefect_settings.to_environment_variables()
    api_url = prefect_env.get("PREFECT_API_URL")
    f"{get_cachedir()}/prefect_{scheduler_name}"

    ######################################
    # Check connection to prefect server
    ######################################
    if api_url is not None:
        check_url = f"{api_url}/health"
        try:
            requests.get(check_url, timeout=60)
        except Exception:
            traceback.print_exc()
            print(f"Could not reach prefect server at: {api_url} from compute node.")
            return 1
    else:
        print("Prefecr server is not running.")
        if scheduler_name == "local":
            print("P-AIRCARS will use ephmeral mode.")
        else:
            print(
                f"First start prefect server to run P-AIRCARS in {scheduler_name} cluster."
            )
            return 1

    if args.jobid is None:
        jobid = get_jobid()
    else:
        jobid = args.jobid

    ###########################################################
    # Estimating jobs memory size (5 times each measurment set)
    ###########################################################
    if os.path.exists(args.target_datadir) is False:
        print(f"Target data directory: {args.target_datadir} does not exist.")
        return
    target_mslist = glob.glob(f"{args.target_datadir}/*.ms")
    if len(target_mslist) == 0:
        print(
            f"No measurement set is present in the target directory: {args.target_datadir}"
        )
        return

    total_ncoarse = 0
    for msname in target_mslist:
        ncoarse = get_ncoarse(msname)
        total_ncoarse += ncoarse
    total_ncoarse = max(1, total_ncoarse)

    if args.cal_datadir:
        if os.path.exists(args.cal_datadir):
            cal_mslist = glob.glob(f"{args.cal_datadir}/*.ms")
            if len(cal_mslist) == 0:
                print(
                    f"No calibrator measurement set is present in: {args.cal_datadir}"
                )
        else:
            print("Calibrator data direcotry does not exist.")

    ###############################################
    # Setup cluster environment
    ###############################################
    if scheduler_name == "local":
        pid = os.getpid()
    elif scheduler_name == "slurm":
        pid = os.environ.get("SLURM_JOB_ID")
    else:
        print("P-AIRCARS is only ready for local or slurm cluster.")
        return 1

    if args.cluster and scheduler_name == "local":
        print(
            "User wants to use cluster architechture, but no job scheduler is available. Stopping P-AIRCARS."
        )
        return

    if args.mem_frac <= 0:
        mem_frac = 0.8
    else:
        mem_frac = args.mem_frac
    if args.cpu_frac <= 0:
        cpu_frac = 0.8
    else:
        cpu_frac = args.cpu_frac

    if args.max_worker is None:
        max_worker = total_ncoarse + 1
    else:
        max_worker = max(int(args.max_worker), total_ncoarse + 1)
    max_worker = max(2, max_worker)  # Minimum 2 workers are needed

    slurm_job = is_slurm_job()
    if not args.cluster or scheduler_name == "local" or slurm_job is False:
        #######################################
        # Set up local cluster
        #######################################
        print("Setting up local cluster....")
        print(f"Maximum allowed worker: {max_worker}")
        dask_client, dask_cluster, dask_dir, nworker = get_local_dask_cluster(
            args.workdir,
            cpu_frac=cpu_frac,
            mem_frac=mem_frac,
            max_worker=max_worker,
        )
        if dask_client is None:
            print("Error occured in creating local cluster.")
            return 1

        scheduler_address = dask_client.scheduler.address
        save_main_process_info(
            pid,
            jobid,
            scheduler_address,
            args.target_datadir,
            os.path.abspath(args.workdir),
            os.path.abspath(args.outdir),
            args.cpu_frac,
            args.mem_frac,
        )
        adaptive = True  # For local cluster, always do adaptive scaling to avoid occupying resources
    else:
        if scheduler_name == "slurm":
            if args.partition is None:
                print("Please provide partition name to submit SLURM jobs.")
                return
            ########################################
            # Setting up slurm cluster
            ########################################
            print("Setting up slurm cluster....")
            print(f"Maximum allowed worker: {max_worker}")
            cluster_result = get_slurm_dask_cluster(
                args.workdir,
                jobid=jobid,
                cpu_frac=cpu_frac,
                mem_frac=mem_frac,
                max_worker=max_worker,
                min_mem=1,
                partition=args.partition,
                account=args.account,
                walltime=args.walltime,
            )
            if cluster_result is None:
                print("Error occured in creating slurm cluster.")
                return 1
            else:
                dask_client, dask_cluster, dask_dir, nworker = cluster_result
            scheduler_address = dask_client.scheduler.address
            save_main_process_info(
                pid,
                jobid,
                scheduler_address,
                args.target_datadir,
                os.path.abspath(args.workdir),
                os.path.abspath(args.outdir),
                args.cpu_frac,
                args.mem_frac,
            )
            adaptive = args.adaptive
            if not adaptive:
                nworker = min(total_ncoarse + 1, nworker)
                scale_worker_and_wait(dask_cluster, dask_client, max(2, nworker))
        else:
            print(
                f"P-AIRCARS is under development for job scheduler: {scheduler_name}. Stopping P-AIRCARS."
            )
            return

    ##########################################
    # Starting pipeline
    ##########################################
    try:
        dask_addr = dask_client.scheduler.address
        print("#########################################")
        print("Starting P-AIRCARS Pipeline....")
        print("#########################################")
        print(f"Total dask workers: {nworker}")
        msg = master_control.with_options(
            flow_run_name=f"paircars_{jobid}",
            task_runner=DaskTaskRunner(address=dask_addr),
        )(
            args.target_datadir,
            args.workdir,
            args.outdir,
            dask_addr,
            target_metafits=args.target_metafits,
            calibrator_datadir=args.cal_datadir,
            calibrator_metafits=args.cal_metafits,
            solar_data=args.solar_data,
            # Pre-calibration
            do_forcereset_weightflag=args.do_forcereset_weightflag,
            do_cal_flag=args.do_cal_flag,
            do_import_model=args.do_import_model,
            # Basic calibration
            do_basic_cal=args.do_basic_cal,
            do_applycal=args.do_applycal,
            only_amplitude=args.only_amplitude,
            redo_basic_cal=args.redo_basic_cal,
            use_solarflagger=args.use_solarflagger,
            # Target data preparation
            freqrange=args.freqrange,
            timerange=args.timerange,
            uvrange=args.cal_uvrange,
            # Polarization calibration
            do_polcal=args.do_polcal,
            # Self-calibration
            do_selfcal=args.do_selfcal,
            do_apply_selfcal=args.do_apply_selfcal,
            do_ap_selfcal=args.do_ap_selfcal,
            solar_selfcal=args.solar_selfcal,
            solint=args.solint,
            redo_selfcal=args.redo_selfcal,
            # Sidereal correction
            do_sidereal_cor=args.do_sidereal_cor,
            do_move_solarcenter=args.do_move_solarcenter,
            # Dynamic spectra
            make_ds=args.make_ds,
            # Imaging
            do_imaging=args.do_imaging,
            do_pbcor=args.do_pbcor,
            weight=args.weight,
            robust=args.robust,
            minuv=args.minuv,
            image_freqres=args.image_freqres,
            image_timeres=args.image_timeres,
            pol=args.pol,
            clean_threshold=args.clean_threshold,
            use_multiscale=args.use_multiscale,
            use_solar_mask=args.use_solar_mask,
            cutout_rsun=args.cutout_rsun,
            make_overlay=args.make_overlay,
            make_msplot=args.make_msplot,
            # Resource settings
            cpu_frac=args.cpu_frac,
            mem_frac=args.mem_frac,
            max_worker=nworker,
            keep_backup=args.keep_backup,
            keep_calibrated_ms=args.keep_calibrated_ms,
            # Remote logging
            masterlog=args.masterlog,
            remote_logger=args.remote_logger,
            jobid=jobid,
            job_password=args.job_password,
            adaptive=adaptive,
        )
        print("##########################################")
        if msg == 0:
            print("P-AIRCARS successfully executed.")
        else:
            print("Issued occured in P-AIRCARS execution.")
        print("##########################################")
    except Exception:
        traceback.print_exc()
    finally:
        time.sleep(5)
        print("Closing clusters...")
        with contextlib.suppress(Exception):
            dask_client.cancel(dask_client.futures)
        with contextlib.suppress(Exception):
            dask_client.close()
        with contextlib.suppress(Exception):
            dask_cluster.close()
        os.system(f"rm -rf {dask_dir}")
