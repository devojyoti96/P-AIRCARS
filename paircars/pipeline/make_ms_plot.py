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
from paircars.utils.logger_utils import (
    SmartDefaultsHelpFormatter,
    clean_shutdown,
    init_logger,
)
from paircars.utils.mwa_ploting_utils import plot_ms_diagnostics
from paircars.utils.resource_utils import drop_cache
from paircars.utils.proc_manage_utils import (
    scale_worker_and_wait,
    get_local_dask_cluster,
    get_scheduler_name,
)
from dask import delayed

logging.getLogger("distributed").setLevel(logging.ERROR)
logging.getLogger("tornado.application").setLevel(logging.CRITICAL)


def main(
    mslist,
    workdir,
    outdir,
    cpu_frac=0.8,
    mem_frac=0.8,
    logfile=None,
    jobid=0,
    start_remote_log=False,
    dask_client=None,
):
    """
    Run the measurement set plots

    Parameters
    ----------
    mslist : str
        Measurment set list (comma separated)
    workdir : str
        Working directory
    outdir : str
        Output directory
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
    """
    cpu_frac = min(0.8, cpu_frac)
    mem_frac = min(0.8, mem_frac)

    mslist = mslist.split(",")

    os.makedirs(workdir, exist_ok=True)
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
                "do_msplot", logfile, jobname=jobname, password=password
            )
    if observer == None:
        print("Remote link or jobname is blank. Not transmiting to remote logger.")

    if len(mslist) == 0:
        print("No measurement set is given.")
        return 1

    dask_cluster = None
    if dask_client is None:
        dask_client, dask_cluster, dask_dir = get_local_dask_cluster(
            workdir,
            mem_frac=mem_frac,
        )
        nworker = min(len(mslist), int(psutil.cpu_count() * cpu_frac) - 1)
        scale_worker_and_wait(dask_cluster, nworker + 1)

    try:
        scheduler_name = get_scheduler_name()
        if scheduler_name == "local":
            njobs = len(mslist)
            total_cpu = max(1, int(psutil.cpu_count() * cpu_frac))
            total_mem = (psutil.virtual_memory().available * mem_frac) / (
                1024**3
            )  # In GB
            n_threads = max(1, int(total_cpu / njobs))
            mem_limit = total_mem / njobs
            cpu_frac = -1
            mem_frac = -1
            print("#################################")
            print(f"Total dask worker: {njobs}")
            print(f"CPU per worker: {n_threads}")
            print(f"Memory per worker: {round(mem_limit,5)} GB")
            print("#################################")
        else:
            njobs = len(dask_client.scheduler_info()["workers"])
            n_threads = -1
            mem_limit = -1
            print("#################################")
            print(f"Total dask worker: {njobs}")
            print("#################################")
        tasks = [
            delayed(plot_ms_diagnostics)(
                msname,
                outdir=outdir,
                ncpu=n_threads,
                total_mem=mem_limit,
                cpu_frac=cpu_frac,
                mem_frac=mem_frac,
            )
            for msname in mslist
        ]
        results = list(dask_client.gather(dask_client.compute(tasks)))
        msg = 0
        final_plots = []
        for res in results:
            success_msg, plots = res
            msg += success_msg
            for p in plots:
                final_plots.append(p)
        print(f"Total measurment sets: {len(mslist)}.")
        print(f"Total successful measurement sets: {len(mslist)-msg}.")
        print(f"Total failed measurement sets: {msg}.")
        print(f"Total plots made: {len(final_plots)}.")
        if msg > 0:
            msg = 1
    except Exception as e:
        traceback.print_exc()
        msg = 1
    finally:
        time.sleep(1)
        for ms in mslist:
            drop_cache(ms)
        drop_cache(workdir)
        drop_cache(outdir)
        clean_shutdown(observer)
    return msg


def cli():
    usage = "Make diagnostic plots of measurement sets"
    parser = argparse.ArgumentParser(
        description=usage, formatter_class=SmartDefaultsHelpFormatter
    )

    # Essential parameters
    basic_args = parser.add_argument_group(
        "###################\nEssential parameters\n###################"
    )
    basic_args.add_argument(
        "mslist", type=str, help="Measurement set list (comma separated)"
    )
    basic_args.add_argument("workdir", type=str, help="Name of work directory")
    basic_args.add_argument("outdir", type=str, help="Output directory")

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

    msg = main(
        args.mslist,
        args.workdir,
        args.outdir,
        cpu_frac=args.cpu_frac,
        mem_frac=args.mem_frac,
        logfile=args.logfile,
        jobid=args.jobid,
        start_remote_log=args.start_remote_log,
    )
    return msg


if __name__ == "__main__":
    result = cli()
    print(
        "\n###################\nPloting measurement set diagnostics are done.\n###################\n"
    )
    os._exit(result)
