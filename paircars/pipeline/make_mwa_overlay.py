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
from paircars.utils.mwa_ploting_utils import make_mwa_overlay, get_all_euv_maps
from paircars.utils.resource_utils import drop_cache
from paircars.utils.image_utils import filter_images
from paircars.utils.proc_manage_utils import (
    scale_worker_and_wait,
    get_local_dask_cluster,
    get_scheduler_name,
)

logging.getLogger("distributed").setLevel(logging.ERROR)
logging.getLogger("tornado.application").setLevel(logging.CRITICAL)


def main(
    imagedir,
    outdir,
    workdir="",
    all_overlay=False,
    cpu_frac=0.8,
    mem_frac=0.8,
    logfile=None,
    jobid=0,
    start_remote_log=False,
    dask_client = None,
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
    mem_frac : float, optional
        Memory fraction to use.
    logfile : str or None, optional
        Path to the log file for saving logs. If None, logging to file is skipped.
    jobid : int, optional
        Numeric job ID used for PID tracking. Default is 0.
    start_remote_log : bool, optional
        Whether to enable remote logging using credentials in the workdir. Default is False.
    dask_client : dask. client
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
    os.chdir(workdir)

    if outdir == "":
        outdir = workdir
    os.makedirs(outdir, exist_ok=True)

    ############
    # Logger
    ############
    observer = None
    if (
        start_remote_log
        and os.path.exists(f"{workdir}/.jobname_password.npy")
        and logfile is not None
    ):
        time.sleep(1)
        jobname, password = np.load(
            f"{workdir}/.jobname_password.npy", allow_pickle=True
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
        
    ###############################
    # Creating dask client
    ###############################
    dask_cluster = None
    nworker=None
    if dask_client is None:
        if mem_frac <= 0:
            mem_frac = 0.8
        if cpu_frac <= 0:
            cpu_frac = 0.8
        image_sizes = [os.stat(image).st_size/1024**3 for image in imagelist]
        min_mem = max(image_sizes)*10
            
        result = get_local_dask_cluster(
            workdir,
            cpu_frac=cpu_frac,
            mem_frac=mem_frac,
            min_mem=min_mem,
            max_worker=len(imagelist) + 1,
        )
        if result is None:
            print("Error occured in creating local cluster.")
            return 1, succeed, failed
        else:
            dask_client, dask_cluster, dask_dir, nworker = result
        scale_worker_and_wait(dask_cluster, dask_client, nworker)
        nthreads = int(psutil.cpu_count()*cpu_frac)
        ncpu = max(1, int(nthreads/nworker))
    else:
        ncpu = os.environ["OMP_NUM_THREADS"]
        if ncpu is None:
            ncpu = 1
        else:
            ncpu = max(1, int(ncpu))
        client_info = dask_client.scheduler_info()["workers"]
        njobs = len(client_info)
        nthreads = ncpu * njobs
        
    try:
        ###############################################################################
        # Filtering only images with bandwidth of 1.28 MHz or more and at 60s intervals
        ###############################################################################
        if all_overlay is False:
            imagelist = filter_images(imagelist, min_time_sep=60.0)
        if len(imagelist) > 0:
            print(f"Total images to overlay: {len(imagelist)}")
            euv_maps = get_all_euv_maps(imagelist, workdir, wavelength=195, ncpu=nthreads)
            # Scatter once
            images_f = dask_client.scatter(imagelist, broadcast=False)
            maps_f = dask_client.scatter(euv_maps, broadcast=True)
            futures = []
            for i in range(len(imagelist)):
                futures.append(
                    dask_client.submit(
                        make_mwa_overlay,
                        images_f[i],
                        maps_f[i],
                        ncpu=ncpu,
                        plot_file_prefix=imagelist[i].split(".fits")[0],
                    )
                )
            results = dask_client.gather(futures)
            outimage_list = []
            for r in results:
                outimage_list.append(r[0])
                print(r[0])
                os.system(f"mv {r[0]} {outdir}")
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
        time.sleep(5)
        drop_cache(imagedir)
        clean_shutdown(observer)
        if dask_cluster is not None:
            dask_client.shutdown()
            dask_client.close()
            dask_cluster.close()
            drop_cache(workdir)
            os.system(f"rm -rf {dask_dir}")
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
    hard_args.add_argument(
        "--mem_frac", type=float, default=0.8, help="Memory fraction to use"
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
        mem_frac=args.mem_frac,
        logfile=args.logfile,
        jobid=args.jobid,
        start_remote_log=args.start_remote_log,
    )
    return msg
