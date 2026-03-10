import logging
import psutil
import numpy as np
import argparse
import traceback
import time
import glob
import sys
import os
import contextlib
from dask.distributed import as_completed
from paircars.utils.logger_utils import (
    SmartDefaultsHelpFormatter,
    clean_shutdown,
    init_logger,
)

from paircars.utils.mwa_ploting_utils import make_mwa_overlay, get_all_euv_maps
from paircars.utils.resource_utils import drop_cache
from paircars.utils.image_utils import filter_images
from paircars.utils.proc_manage_utils import (
    scale_worker_and_wait,
    get_local_dask_cluster,
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
    dask_client=None,
):
    cpu_frac = min(0.8, abs(cpu_frac))
    mem_frac = min(0.8, abs(mem_frac))

    if workdir == "":
        workdir = f"{imagedir}/workdir"
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)

    if outdir == "":
        outdir = workdir
    os.makedirs(outdir, exist_ok=True)

    ####################
    # Logger
    ####################
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

    if observer is None:
        print("Remote link or jobname is blank. Not transmiting to remote logger.")

    imagelist = glob.glob(f"{imagedir}/*.fits")

    if len(imagelist) == 0:
        print("No image in the image directory.")
        return 1, 0, 0

    succeed = 0
    failed = len(imagelist)
    
    ###############################################
    # Filter images
    ###############################################
    if not all_overlay:
        imagelist = filter_images(imagelist, min_time_sep=60.0)

    if len(imagelist) == 0:
        return 1, 0, 0

    print(f"Total images to overlay: {len(imagelist)}")
    
    try: 
        nthreads = int(os.environ.get("OMP_NUM_THREADS", 1))
        ###########################
        # Download EUV maps
        ###########################
        euv_fits_images = get_all_euv_maps(
            imagelist,
            workdir,
            wavelength=195,
            ncpu=nthreads,
        )
        if len(euv_fits_images)==0:
            print("No EUV images downloaded.")
            return 1, succeed, failed
    except Exception:
        print("Error occured in EUV fits downloading.")
        traceback.print_exc()
        return 1, succeed, failed
        

    ###############################
    # Dask cluster
    ###############################
    dask_cluster = None
    nworker = None
    if dask_client is None:
        image_sizes = [os.stat(image).st_size / 1024**3 for image in imagelist]
        max_image_size = max(image_sizes)
        min_mem = max(2, round(5 * max_image_size, 2))
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

        dask_client, dask_cluster, dask_dir, njobs = result
        scale_worker_and_wait(dask_cluster, dask_client, njobs)
        nthreads = int(psutil.cpu_count() * cpu_frac)
    else:
        ncpu = int(os.environ.get("OMP_NUM_THREADS", 1))
        client_info = dask_client.scheduler_info()["workers"]
        njobs = len(client_info)
        nthreads = ncpu * njobs

    try:
        ###########################
        # Overlay tasks
        ###########################
        print("Start making overlays....")
        results = []
        batch_size = max(1, njobs-1)
        for i in range(0, len(imagelist), batch_size):
            batch_imgs = imagelist[i:i+batch_size]
            batch_euv = euv_fits_images[i:i+batch_size]
            futures = []
            for img, euv_fits in zip(batch_imgs, batch_euv):
                futures.append(
                    dask_client.submit(
                        make_mwa_overlay,
                        img,
                        euv_fits,
                        workdir,
                        plot_file_prefix=os.path.basename(img).replace(".fits", ""),
                        outdirs=[outdir],
                        verbose=True,
                        pure=False,
                        retries=2,
                    )
                )

            ###########################
            # Collect batch results
            ###########################
            for f in as_completed(futures):
                try:
                    r = f.result(timeout=300)
                    results.append(r)
                except Exception as e:
                    print("Overlay failed:", e)

            # free worker memory
            dask_client.cancel(futures)
                
        ###########################
        # Move outputs
        ###########################
        outimage_list = []
        for r in results:
            if r is not None:
                outimage_list.append(r[0])
                if os.path.dirname(os.path.abspath(r[0]))!=os.path.abspath(outdir):
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
    except Exception:
        traceback.print_exc()
        msg = 1
    finally:
        os.system(f"rm -rf {imagedir}/*aia*.fits")
        os.system(f"rm -rf {imagedir}/*suvi*.fits")
        time.sleep(5)
        drop_cache(imagedir)
        clean_shutdown(observer)
        if dask_cluster is not None:
            with contextlib.suppress(Exception):
                dask_client.cancel(dask_client.futures)
            with contextlib.suppress(Exception):
                dask_client.close()
            with contextlib.suppress(Exception):
                dask_cluster.close()
            drop_cache(workdir)
            os.system(f"rm -rf {dask_dir}")
    return msg, succeed, failed


def cli():
    usage = "Overlay MWA images on EUV images"
    parser = argparse.ArgumentParser(
        description=usage,
        formatter_class=SmartDefaultsHelpFormatter,
    )
    basic_args = parser.add_argument_group("Essential parameters")
    basic_args.add_argument("imagedir", type=str)
    basic_args.add_argument("outdir", type=str)
    basic_args.add_argument("--workdir", type=str, default="")

    adv_args = parser.add_argument_group("Advanced parameters")
    adv_args.add_argument("--all_overlay", action="store_true")
    adv_args.add_argument("--start_remote_log", action="store_true")

    hard_args = parser.add_argument_group("Resource parameters")
    hard_args.add_argument("--cpu_frac", type=float, default=0.8)
    hard_args.add_argument("--mem_frac", type=float, default=0.8)
    hard_args.add_argument("--logfile", type=str, default=None)
    hard_args.add_argument("--jobid", type=int, default=0)

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
