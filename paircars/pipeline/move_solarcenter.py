import logging
import psutil
import numpy as np
import argparse
import traceback
import time
import glob
import sys
import os
from dask import delayed
from paircars.utils.logger_utils import (
    SmartDefaultsHelpFormatter,
    clean_shutdown,
    init_logger,
)
from paircars.utils.ms_metadata import get_ms_size
from paircars.utils.mwa_utils import get_ncoarse
from paircars.utils.proc_manage_utils import (
    scale_worker_and_wait,
    get_local_dask_cluster,
    get_scheduler_name,
)
from paircars.utils.resource_utils import drop_cache
from paircars.utils.sunpos_utils import move_to_sun

logging.getLogger("distributed").setLevel(logging.ERROR)
logging.getLogger("tornado.application").setLevel(logging.CRITICAL)


def main(
    mslist,
    workdir="",
    cpu_frac=0.8,
    mem_frac=0.8,
    logfile=None,
    jobid=0,
    start_remote_log=False,
    dask_client=None,
):
    """
    Run the flagging pipeline for a measurement set.

    Parameters
    ----------
    mslist : str
        List of measurement sets
    workdir : str, optional
        Working directory
    cpu_frac : float, optional
        Fraction of total CPU resources to use. Default is 0.8.
    mem_frac : float, optional
        Fraction of total memory to use. Default is 0.8.
    logfile : str or None, optional
        Path to the log file for saving logs. If None, logging to file is skipped.
    jobid : int, optional
        Numeric job ID used for PID tracking. Default is 0.
    start_remote_log : bool, optional
        Whether to enable remote logging using credentials in the workdir. Default is False.
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
    """
    cpu_frac = min(0.8, abs(cpu_frac))
    mem_frac = min(0.8, abs(mem_frac))

    mslist = mslist.split(",")

    if workdir == "":
        workdir = os.path.dirname(os.path.abspath(mslist[0])) + "/workdir"
    os.makedirs(workdir, exist_ok=True)

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
                "do_flagging", logfile, jobname=jobname, password=password
            )
    if observer == None:
        print("Remote link or jobname is blank. Not transmiting to remote logger.")

    if len(mslist) == 0:
        print("Please provide a valid measurement set list.")
        return 1, 0, 0
    else:
        succeed = 0
        failed = len(mslist)

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
            return 1, succeed, failed
        else:
            dask_client, dask_cluster, dask_dir, nworker = result
        scale_worker_and_wait(dask_cluster, dask_client, nworker)

    try:
        tasks = [delayed(move_to_sun)(msname) for msname in mslist]
        results = list(dask_client.gather(dask_client.compute(tasks)))
        failed = sum(results)
        succeed = len(mslist) - failed
        print(f"Total measurement sets: {len(mslist)}")
        print(f"Total success: {succeed}")
        print(f"Total failure: {failed}")
        if len(mslist) == failed:
            msg = 1
        else:
            msg = 0
    except Exception as e:
        traceback.print_exc()
        msg = 1
    finally:
        time.sleep(5)
        clean_shutdown(observer)
        if dask_cluster is not None:
            dask_client.shutdown()
            dask_client.close()
            dask_cluster.close()
            for ms in mslist:
                drop_cache(ms)
            drop_cache(workdir)
            os.system(f"rm -rf {dask_dir}")
    return msg, succeed, failed


def cli():
    usage = "Move phasecenter of the measurement set to the Sun"
    parser = argparse.ArgumentParser(
        description=usage, formatter_class=SmartDefaultsHelpFormatter
    )

    # Essential parameters
    basic_args = parser.add_argument_group(
        "###################\nEssential parameters\n###################"
    )
    basic_args.add_argument(
        "mslist", type=str, help="Measurement set list, comma seperated"
    )
    basic_args.add_argument(
        "--workdir", type=str, default="", help="Name of work directory"
    )

    # Advanced switches
    adv_args = parser.add_argument_group(
        "###################\nAdvanced parameters\n###################"
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
        args.mslist,
        workdir=args.workdir,
        cpu_frac=args.cpu_frac,
        mem_frac=args.mem_frac,
        logfile=args.logfile,
        jobid=args.jobid,
        start_remote_log=args.start_remote_log,
    )
    return msg
