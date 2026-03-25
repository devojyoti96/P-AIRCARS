import logging
import numpy as np
import argparse
import warnings
import time
import glob
import sys
import os
import subprocess
from astropy.io import fits
from astropy.wcs import FITSFixedWarning
from dask import delayed
from paircars.utils.basic_utils import (
    print_banner,
    capture_all_output,
)
from paircars.utils.image_utils import generate_tb_map
from paircars.utils.logger_utils import (
    SmartDefaultsHelpFormatter,
    clean_shutdown,
    init_logger,
    get_logger_safe,
)
from paircars.utils.mwa_utils import freq_to_MWA_coarse
from paircars.utils.mwa_ploting_utils import save_in_hpc, plot_in_hpc
from paircars.utils.proc_manage_utils import (
    scale_worker_and_wait,
    get_local_dask_cluster,
)
from paircars.utils.resource_utils import drop_cache

logging.getLogger("distributed").setLevel(logging.ERROR)
logging.getLogger("tornado.application").setLevel(logging.CRITICAL)
warnings.simplefilter("ignore", FITSFixedWarning)


def run_pbcor_wrapper(*args, **kwargs):
    with capture_all_output() as (out, err):
        result = run_pbcor(*args, **kwargs)
        return args[0], result, out.getvalue(), err.getvalue()


def get_fits_freq(image_file):
    hdr = fits.getheader(image_file)
    keys = hdr.keys()
    if "CTYPE3" in keys and hdr["CTYPE3"] == "FREQ":
        freq = hdr["CRVAL3"]
        return freq
    elif "CTYPE4" in keys and hdr["CTYPE4"] == "FREQ":
        freq = hdr["CRVAL4"]
        return freq
    else:
        print(f"No frequency axis in image: {image_file}.")
        return


def run_pbcor(
    imagename,
    metafits,
    pbdir,
    pbcor_dir,
    leakage_file="",
    restore=False,
    jobid=0,
    ncpu=1,
    verbose=False,
):
    """
    Run single image primary beam correction

    Parameters
    ----------
    imagename : str
        Imagename
    metafits : str
        Metafits file
    pbdir : str
        Primary beam directory
    pbcor_dir : str
        Primary beam corrected image directory
    leakage_file : str, optional
        Leakage information file
    restore : bool, optional
        Restore primary beam correction
    jobid : int, optional
        Job ID
    ncpu : int, optional
        Number of CPU threads to use
    verbose : bool, optional
        Verbose output

    Returns
    -------
    int
        Success message
    """
    ncpu = max(1, ncpu)
    freq = get_fits_freq(imagename)
    outfile = f"{pbcor_dir}/{os.path.basename(imagename).split('.fits')[0]}_pbcor.fits"
    pbfile = f"{pbdir}/freq_{freq}.npy"
    cmd = [
        "run-mwa-singlepbcor",
        "--num_threads",
        str(ncpu),
        "--interpolated",
    ]
    cmd.append("--pb_jones_file")
    cmd.append(pbfile)
    if os.path.exists(pbfile) is False:
        cmd.append("--save_pb")

    if restore:
        cmd.append("--restore")

    if leakage_file != "" and os.path.exists(leakage_file):
        cmd.append("--leakage_file")
        cmd.append(f"{leakage_file}")

    cmd.append(imagename)
    cmd.append(metafits)
    cmd.append(outfile)

    if verbose:
        print("Executing:", " ".join(cmd))

    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,  # Set to True if you want to raise on error
        )
        if verbose or result.returncode != 0:
            print(result.stdout)
        return result.returncode
    except Exception as e:
        print(f"Exception during primary beam correction: {e}")
        return 1


def get_leakage_file(image, leakage_dir):
    """
    Get leakage file for the image

    Parameters
    ----------
    image : str
        Imagename
    leakage_dir : str, optional
        Leakage file directory

    Returns
    -------
    str
        Leakage file name
    """
    header = fits.getheader(image)
    if header["CTYPE3"] == "FREQ":
        image_freq = round(float(header["CRVAL3"]) / 10**6, 3)
    elif header["CTYPE4"] == "FREQ":
        image_freq = round(float(header["CRVAL4"]) / 10**6, 3)
    else:
        image_freq = -1
    if image_freq > 0 and leakage_dir != 0 and os.path.exists(leakage_dir):
        image_coarse = freq_to_MWA_coarse(image_freq)
        leakage_file_list = glob.glob(
            f"{leakage_dir}/selfcal_*{image_coarse}_*.leakage"
        )
        if len(leakage_file_list) > 0:
            leakage_file = leakage_file_list[0]
        else:
            leakage_file = ""
    else:
        leakage_file = ""
    return leakage_file


def pbcor_all_images(
    imagedir,
    metafits,
    dask_client,
    leakage_dir="",
    make_TB=True,
    make_plots=True,
    restore=False,
    jobid=0,
    n_threads=1,
    mem_limit=1,
    njobs=1,
    logger=None,
):
    """
    Correct primary beam of MWA for images in a directory

    Parameters
    ----------
    imagedir : str
        Name of the image directory
    metafits : str
        Metafits file
    dask_client : dask.client
        Dask client
    leakage_dir : str, optional
        Leakage file directory
    make_TB : bool, optional
        Make brightness temperature map
    make_plots : bool, optional
        Make plots
    restore : bool, optional
        Restore primary beam correction
    jobid : int, optional
        Job ID
    n_threads : int, optional
        CPU threads to use
    mem_limit : float, optional
        Memory to use in GB
    njobs : int, optional
        Number of parallel jobs

    Returns
    -------
    int
        Success message
    int
        Succeeded image number
    int
        Failed image number
    """
    if logger is None:
        logger = get_logger_safe()

    imagedir = imagedir.rstrip("/")
    pbdir = f"{os.path.dirname(imagedir)}/pbdir"
    pbcor_dir = f"{os.path.dirname(imagedir)}/pbcor_images"
    os.makedirs(pbdir, exist_ok=True)
    os.makedirs(pbcor_dir, exist_ok=True)
    successful_pbcor = 0
    succeed = 0
    failed = 0
    try:
        images = glob.glob(f"{imagedir}/*.fits")
        if make_TB:
            tb_dir = f"{os.path.dirname(imagedir)}/tb_images"
            os.makedirs(tb_dir, exist_ok=True)
        if len(images) == 0:
            logger.critical(f"No image is present in image directory: {imagedir}")
            return 1, 0, 0
        else:
            succeed = 0
            failed = len(images)

        first_set = []
        remaining_set = []
        freqs = []
        for image in images:
            freq = get_fits_freq(image)
            if freq in freqs:
                remaining_set.append(image)
            else:
                freqs.append(freq)
                first_set.append(image)

        if len(first_set) > 0:
            tasks = []
            for image in first_set:
                leakage_file = get_leakage_file(image, leakage_dir=leakage_dir)
                task = delayed(run_pbcor_wrapper)(
                    image,
                    metafits,
                    pbdir,
                    pbcor_dir,
                    leakage_file=leakage_file,
                    restore=restore,
                    jobid=jobid,
                    ncpu=n_threads,
                )
                tasks.append(task)
            result_wrapper = []
            logger.info("Start correcting first set of images.")
            for i in range(0, len(tasks), njobs):
                batch = tasks[i : i + njobs]
                futures = dask_client.compute(batch)
                result_wrapper.extend(dask_client.gather(futures))
            result_wrapper = list(result_wrapper)
            results = []
            for r in result_wrapper:
                results.append(r[1])
                logger.debug("================")
                logger.debug(f"Worker log for: {os.path.basename(r[0])}")
                logger.debug("================")
                for line in r[2].splitlines():
                    logger.debug(line)
                for line in r[3].splitlines():
                    logger.debug(line)

            for r in results:
                if r == 0:
                    successful_pbcor += 1

        if len(remaining_set) > 0:
            tasks = []
            for image in remaining_set:
                leakage_file = get_leakage_file(image, leakage_dir=leakage_dir)
                task = delayed(run_pbcor_wrapper)(
                    image,
                    metafits,
                    pbdir,
                    pbcor_dir,
                    leakage_file=leakage_file,
                    restore=restore,
                    jobid=jobid,
                    ncpu=n_threads,
                )
                tasks.append(task)

            result_wrapper = []
            logger.info("Correcting remaining images of different timestamps.")
            for i in range(0, len(tasks), njobs):
                batch = tasks[i : i + njobs]
                futures = dask_client.compute(batch)
                result_wrapper.extend(dask_client.gather(futures))
            result_wrapper = list(result_wrapper)
            results = []
            for r in result_wrapper:
                results.append(r[1])
                logger.debug("================")
                logger.debug(f"Worker log for: {os.path.basename(r[0])}")
                logger.debug("================")
                for line in r[2].splitlines():
                    logger.debug(line)
                for line in r[3].splitlines():
                    logger.debug(line)

            for r in results:
                if r == 0:
                    successful_pbcor += 1

        ############################################
        # Saving fits in helioprojective coordinates
        ############################################
        if successful_pbcor > 0:
            hpcdir = f"{pbcor_dir}/hpcs"
            pbcor_images = glob.glob(f"{pbcor_dir}/*.fits")
            os.makedirs(hpcdir, exist_ok=True)
            logger.info(
                "Saving primary beam corrected images helioprojective coordinates."
            )
            for image in pbcor_images:
                save_in_hpc(image, outdir=hpcdir)
            if make_plots:
                logger.info(
                    "Making plots of primary beam corrected images in helioprojective coordinates."
                )
                pngdir = f"{pbcor_dir}/pngs"
                os.makedirs(pngdir, exist_ok=True)
                for image in pbcor_images:
                    try:
                        plot_in_hpc(
                            image,
                            draw_limb=True,
                            extensions=["png"],
                            outdirs=[pngdir],
                        )
                    except BaseException:
                        junkpng = f"{pngdir}/{os.path.basename(image).split('.fits')[0]}.png.junk"
                        os.system(f"touch {junkpng}")

        ####################################
        # Making brightness temperature maps
        ####################################
        if successful_pbcor > 0 and make_TB:
            logger.info("Making brightness temperature maps.")
            for pbcor_image in pbcor_images:
                tb_image = (
                    tb_dir
                    + "/"
                    + os.path.basename(pbcor_image).split(".fits")[0]
                    + "_TB.fits"
                )
                generate_tb_map(pbcor_image, outfile=tb_image)

            ############################################
            # Saving fits in helioprojective coordinates
            ###########################################
            hpcdir = f"{tb_dir}/hpcs"
            tb_images = glob.glob(f"{tb_dir}/*.fits")
            os.makedirs(hpcdir, exist_ok=True)
            logger.info(
                "Saving brightness temperature maps helioprojective coordinates."
            )
            for image in tb_images:
                save_in_hpc(image, outdir=hpcdir)

            if make_plots:
                logger.info("Making plots of brightness temperature maps.")
                pngdir = f"{tb_dir}/pngs"
                os.makedirs(pngdir, exist_ok=True)
                for image in tb_images:
                    plot_in_hpc(
                        image,
                        draw_limb=True,
                        extensions=["png"],
                        outdirs=[pngdir],
                    )

        #########################################
        # Final calculations
        #########################################
        logger.info(f"Total input images: {len(images)}")
        succeed = successful_pbcor
        failed = len(images) - succeed
        if successful_pbcor > 0:
            logger.info(f"Total primary beam corrected images: {len(pbcor_images)}")
            if make_TB:
                logger.info(f"Total brightness temperatures maps: {len(tb_images)}")
        else:
            logger.error("Total primary beam corrected images: 0")
        msg = 0
    except Exception:
        logger.exception("Exception occured in primary beam correction.", exc_info=True)
        msg = 1
    finally:
        os.system(f"rm -rf {pbdir}")
        return msg, succeed, failed


def main(
    imagedir,
    metafits,
    workdir="",
    leakage_dir="",
    make_TB=True,
    make_plots=True,
    restore=False,
    cpu_frac=0.8,
    mem_frac=0.8,
    logfile=None,
    jobid=0,
    verbose=False,
    start_remote_log=False,
    dask_client=None,
):
    """
    Primary beam correction of MWA for a sets of images in a directory

    Parameters
    ----------
    imagedir : str
        Image directory
    metafits : str
        Metafits file
    workdir : str, optional
        Work directory
    leakage_dir : str, optional
        Leakage file directory
    make_TB : bool, optional
        Make brightness temperature map or not
    make_plots : bool, optional
        Make png plots
    restore : bool, optional
        Restore primary beam correction
    cpu_frac : float,optional
        CPU fraction
    mem_frac : float
        Memory fraction
    logfile : str, optional
        Log file
    jobid : str, optional
        Job ID
    verbose : bool, optional
        Verbose logs
    start_remote_log : bool, optional
        Start remote logger
    dask_client : dask.client, optional
        Dask client

    Returns
    -------
    int
        Success message
    int
        Succeeded image number
    int
        Failed image number
    """
    logger = get_logger_safe()
    if verbose:
        logger.setLevel(logging.DEBUG)

    cpu_frac = min(0.8, abs(cpu_frac))
    mem_frac = min(0.8, abs(mem_frac))

    if workdir == "":
        workdir = imagedir + "/workdir"
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)
    logger.debug(f"Current working directory: {os.getcwd()}")

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
                "all_pbcor", logfile, jobname=jobname, password=password
            )

    dask_cluster = None
    if dask_client is None:
        dask_client, dask_cluster, dask_dir, nworker = get_local_dask_cluster(
            workdir,
            cpu_frac=cpu_frac,
            mem_frac=mem_frac,
        )
        if dask_client is None:
            logger.critical("Error occured in creating local cluster.")
            return 1
        scale_worker_and_wait(dask_cluster, dask_client, nworker)

    succeed = 0
    failed = 0
    try:
        for banner in print_banner(
            "Starting primary beam corrections.", no_print=True
        ).splitlines():
            logger.info(banner)
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

        logger.info("#################################")
        logger.info(f"Total dask worker: {njobs}")
        logger.info(f"CPU per worker: {n_threads}")
        logger.info(f"Memory per worker: {mem_limit} GB")
        logger.info("#################################")

        if os.path.exists(imagedir):
            msg, succeed, failed = pbcor_all_images(
                imagedir,
                metafits,
                dask_client,
                leakage_dir=leakage_dir,
                make_TB=make_TB,
                make_plots=make_plots,
                restore=restore,
                jobid=jobid,
                n_threads=n_threads,
                mem_limit=mem_limit,
                njobs=njobs,
                logger=logger,
            )
        else:
            logger.critical("Please provide correct image directory path.")
            msg = 1
    except Exception:
        logger.exception("Exception occured in primary beam correction.", exc_info=True)
        msg = 1
    finally:
        time.sleep(5)
        clean_shutdown(observer)
        if dask_cluster is not None:
            dask_client.shutdown()
            dask_client.close()
            dask_cluster.close()
            drop_cache(imagedir)
            drop_cache(workdir)
            os.system(f"rm -rf {dask_dir}")
    return msg, succeed, failed


def cli():
    parser = argparse.ArgumentParser(
        description="Correct all images for MWA full-pol averaged primary beam",
        formatter_class=SmartDefaultsHelpFormatter,
    )

    # Essential parameters
    basic_args = parser.add_argument_group(
        "###################\nEssential parameters\n###################"
    )
    basic_args.add_argument("imagedir", help="Path to image directory")
    basic_args.add_argument("--metafits", required=True, help="Metafits file")
    basic_args.add_argument("--workdir", default="", help="Path to work directory")

    # Advanced parameters
    adv_args = parser.add_argument_group(
        "###################\nAdvanced parameters\n###################"
    )
    adv_args.add_argument(
        "--leakage_dir",
        type=str,
        default="",
        help="Leakage file directory",
    )
    adv_args.add_argument(
        "--no_make_TB",
        action="store_false",
        dest="make_TB",
        help="Do not generate brightness temperature map",
    )
    adv_args.add_argument(
        "--no_make_plots",
        action="store_false",
        dest="make_plots",
        help="Do not make png plots",
    )
    adv_args.add_argument(
        "--restore",
        action="store_true",
        dest="restore",
        help="Restore primary beam correction",
    )
    adv_args.add_argument(
        "--verbose",
        action="store_true",
        help="Verbose logs",
    )
    adv_args.add_argument(
        "--jobid", type=int, default=0, help="Job ID for logging and PID tracking"
    )

    # Resource management parameters
    hard_args = parser.add_argument_group(
        "###################\nHardware resource management parameters\n###################"
    )
    hard_args.add_argument(
        "--cpu_frac", type=float, default=0.8, help="CPU usage fraction"
    )
    hard_args.add_argument(
        "--mem_frac", type=float, default=0.8, help="Memory usage fraction"
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        return 1

    args = parser.parse_args()

    msg, _, _ = main(
        args.imagedir,
        args.metafits,
        workdir=args.workdir,
        leakage_dir=args.leakage_dir,
        make_TB=args.make_TB,
        make_plots=args.make_plots,
        restore=args.restore,
        cpu_frac=args.cpu_frac,
        mem_frac=args.mem_frac,
        jobid=args.jobid,
        verbose=args.verbose,
    )
    return msg
