import logging
import psutil
import numpy as np
import argparse
import traceback
import time
import glob
import sys
import os
import dask
from astropy.io import fits
from paircars.utils.logger_utils import (
    SmartDefaultsHelpFormatter,
    clean_shutdown,
    init_logger,
)
from paircars.utils.basic_utils import timestamp_to_mjdsec
from paircars.utils.mwa_ploting_utils import make_mwa_overlay
from paircars.utils.resource_utils import drop_cache
from paircars.utils.proc_manage_utils import get_scheduler_name
from dask import delayed

logging.getLogger("distributed").setLevel(logging.ERROR)
logging.getLogger("tornado.application").setLevel(logging.CRITICAL)


def main(
    imagedir,
    outdir,
    workdir="",
    all_overlay=False,
    cpu_frac=0.8,
    logfile=None,
    jobid=0,
    start_remote_log=False,
    dask_client=None,
):
    """
    Run the EUV overlays

    Parameters
    ----------
    imagedir : str
        Image directory
    outdir : str
        Output directory
    workdir : str, optional
        Working directory
    all_overlay : bool, optional
        Whether make overlays of all images or not
    cpu_frac : float, optional
        Fraction of total CPU resources to use. Default is 0.8.
    logfile : str or None, optional
        Path to the log file for saving logs. If None, logging to file is skipped.
    jobid : int, optional
        Numeric job ID used for PID tracking. Default is 0.
    start_remote_log : bool, optional
        Whether to enable remote logging using credentials in the workdir. Default is False.
    dask_client : dask.client
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
    cpu_frac = min(0.8, abs(cpu_frac))

    if workdir == "":
        workdir = f"{imagedir}/workdir"
    os.makedirs(workdir, exist_ok=True)

    if outdir == "":
        outdir = workdir
    os.makedirs(outdir, exist_ok=True)

    ############
    # Logger
    ############
    observer = None
    if (
        start_remote_log
        and os.path.exists(f"{workdir}/jobname_password.npy")
        and logfile is not None
    ):
        time.sleep(1)
        jobname, password = np.load(
            f"{workdir}/jobname_password.npy", allow_pickle=True
        )
        if os.path.exists(logfile):
            observer = init_logger(
                "do_overlay", logfile, jobname=jobname, password=password
            )
    if observer == None:
        print("Remote link or jobname is blank. Not transmiting to remote logger.")

    imagelist = glob.glob(f"{imagedir}/*.fits")

    if len(imagelist) == 0:
        print("No image in the image directory.")
        return 1, 0, 0
    else:
        succeed = 0
        failed = len(imagelist)

    dask_cluster = None
    if dask_client is None:
        result = get_local_dask_cluster(
            workdir,
            cpu_frac=cpu_frac,
            mem_frac=cpu_frac,
        )
        if result is None:
            print("Error occured in creating local cluster.")
            return 1, succeed, failed
        else:
            dask_client, dask_cluster, dask_dir, nworker = result
        scale_worker_and_wait(dask_cluster, dask_client, nworker)

    try:
        ###############################################################################
        # Filtering only images with bandwidth of 1.28 MHz or more and at 10s intervals
        ###############################################################################
        if all_overlay is False:
            bws = []
            for image in imagelist:
                header = fits.getheader(image)
                keys = header.keys()
                if "CTYPE3" in keys and header["CTYPE3"] == "FREQ":
                    bw = round(float(header["CDELT3"]) / 10**6, 2)
                elif "CTYPE3" in keys and header["CTYPE4"] == "FREQ":
                    bw = round(float(header["CDELT4"]) / 10**6, 2)
                else:
                    bw = -1
                bws.append(bw)
            max_bw = max(bws)
            bws = np.array(bws)
            pos = np.where(bws == max_bw)
            imagelist = np.array(imagelist)
            filtered_imagelist = imagelist[pos]

            last_mjdsec = 0.0
            final_imagelist = []
            timelist = []
            for image in filtered_imagelist:
                header = fits.getheader(image)
                timeobs = header["DATE-OBS"].split(".")[0]
                mjdsec = timestamp_to_mjdsec(timeobs, date_format=1)
                if (mjdsec - last_mjdsec) >= 10.0:
                    final_imagelist.append(image)
                    timelist.append(mjdsec)
                    last_mjdsec = max(timelist)
            imagelist = final_imagelist

        if len(imagelist) > 0:
            print(f"Total images to overlay: {len(imagelist)}")
            """scheduler_name = get_scheduler_name()
            if scheduler_name == "local" or dask_client is None:
                ncpu = max(1, int(psutil.cpu_count() * cpu_frac))
                outimage_list = []
                for image in imagelist:
                    outimage = make_mwa_overlay(
                        image,
                        plot_file_prefix=os.path.basename(image).split(".fits")[0]
                        + "_euv_mwa_overlay",
                        extensions=["png"],
                        outdirs=[outdir],
                        keep_euv_fits=True,
                        ncpu=ncpu,
                        verbose=False,
                    )
            else:"""
            tasks = []
            ncpu = os.environ.get("OMP_NUM_THREADS")
            for image in imagelist:
                task = delayed(make_mwa_overlay)(
                    image,
                    plot_file_prefix=os.path.basename(image).split(".fits")[0]
                    + "_euv_mwa_overlay",
                    extensions=["png"],
                    outdirs=[outdir],
                    keep_euv_fits=True,
                    npcu=ncpu,
                    verbose=False,
                )
                tasks.append(task)
            futures = dask_client.compute(tasks)
            outimage_list = list(dask_client.gather(futures))
            outimage_list.append(outimage)
            if len(outimage_list) == 0:
                print("No overlay is made.")
                msg = 1
                succeed = 0
                failed = len(imagelist)
            else:
                print(f"Total images: {len(imagelist)}")
                print(f"Total overlays: {len(outimage_list)}")
                msg = 0
                succeed = len(outimage_list)
                failed = len(imagelist) - succeed
        else:
            msg = 1
    except Exception as e:
        traceback.print_exc()
        msg = 1
    finally:
        os.system(f"rm -rf {imagedir}/*aia*.fits")
        os.system(f"rm -rf {imagedir}/*suvi*.fits")
        time.sleep(1)
        drop_cache(imagedir)
        drop_cache(workdir)
        drop_cache(outdir)
        clean_shutdown(observer)
    return msg, succeed, failed


def cli():
    usage = "Overlay MWA images on EUV images"
    parser = argparse.ArgumentParser(
        description=usage, formatter_class=SmartDefaultsHelpFormatter
    )

    # Essential parameters
    basic_args = parser.add_argument_group(
        "###################\nEssential parameters\n###################"
    )
    basic_args.add_argument("imagedir", type=str, help="Image directory")
    basic_args.add_argument("outdir", type=str, help="Output directory")
    basic_args.add_argument(
        "--workdir", type=str, default="", help="Name of work directory"
    )

    # Advanced switches
    adv_args = parser.add_argument_group(
        "###################\nAdvanced parameters\n###################"
    )
    adv_args.add_argument(
        "--all_overlay", action="store_true", help="Make overlays of all images"
    )
    adv_args.add_argument(
        "--start_remote_log", action="store_true", help="Start remote logging"
    )

    # Resource management parameters
    hard_args = parser.add_argument_group(
        "###################\nHardware resource management parameters\n###################"
    )
    hard_args.add_argument(
        "--cpu_frac", type=float, default=0.8, help="CPU fraction to use"
    )
    hard_args.add_argument("--logfile", type=str, default=None, help="Log file")
    hard_args.add_argument("--jobid", type=int, default=0, help="Job ID")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        return 1

    args = parser.parse_args()

    msg, _, _ = main(
        args.imagedir,
        args.outdir,
        workdir=args.workdir,
        all_overlay=args.all_overlay,
        cpu_frac=args.cpu_frac,
        logfile=args.logfile,
        jobid=args.jobid,
        start_remote_log=args.start_remote_log,
    )
    return msg
