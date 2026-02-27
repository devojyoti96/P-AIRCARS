import logging
import psutil
import dask
import numpy as np
import argparse
import traceback
import time
import sys
import os
from casatools import msmetadata
from dask import delayed
from paircars.utils.casatasks import single_mstransform
from paircars.utils.logger_utils import (
    SmartDefaultsHelpFormatter,
    clean_shutdown,
    init_logger,
)
from paircars.utils.ms_metadata import get_timeranges
from paircars.utils.mwa_utils import get_MWA_coarse_bands
from paircars.utils.proc_manage_utils import (
    scale_worker_and_wait,
    get_local_dask_cluster,
    get_scheduler_name,
)
from paircars.utils.resource_utils import drop_cache


logging.getLogger("distributed").setLevel(logging.ERROR)
logging.getLogger("tornado.application").setLevel(logging.CRITICAL)


def chanlist_to_str(lst):
    lst = sorted(lst)
    ranges = []
    start = lst[0]
    for i in range(1, len(lst)):
        if lst[i] != lst[i - 1] + 1:
            if lst[i - 1] > start:
                ranges.append(f"{start}~{lst[i - 1]}")
            elif lst[i - 1] == start:
                ranges.append(f"{start}")
            start = lst[i]
    if lst[-1] > start:
        ranges.append(f"{start}~{lst[-1]}")
    elif lst[-1] == start:
        ranges.append(f"{start}")
    return ";".join(ranges)


def split_target_scans(
    msname,
    workdir,
    timeres,
    freqres,
    datacolumn,
    scan=1,
    prefix="targets",
    time_interval=-1,
    time_window=-1,
    quack_timestamps=-1,
    force_split=False,
    n_threads=-1,
):
    """
    Split target scans

    Parameters
    ----------
    msname : str
        Measurement set
    workdir : str
        Work directory
    timeres : float
        Time resolution in seconds
    freqres : float
        Frequency resolution in MHz
    datacolumn : str
        Data column to split
    scan : int
        Scan to split
    prefix : str, optional
        Splited ms prefix
    time_interval : float
        Time interval in seconds
    time_window : float
        Time window in seconds
    quack_timestamps : int, optional
        Number of timestamps ignored at the start and end of each scan
    force_split : bool, optional
        Force split
    n_threads : int, optional
        Number of threads to use

    Returns
    -------
    list
        Splited ms list
    """
    n_threads = max(1, n_threads)
    print(f"Spliting measurement set: {msname}")
    try:
        os.chdir(workdir)
        #######################################
        # Extracting time frequency information
        #######################################
        msmd = msmetadata()
        msmd.open(msname)
        chanres = msmd.chanres(0, unit="MHz")[0]
        freqs = msmd.chanfreqs(0, unit="MHz")
        bw = max(freqs) - min(freqs)
        nchan = msmd.nchan(0)
        msmd.close()
        if freqres > 0:  # Image resolution is in MHz
            chanwidth = int(freqres / chanres)
            if chanwidth < 1:
                chanwidth = 1
        else:
            chanwidth = 1
        if timeres > 0:  # Image resolution is in seconds
            timebin = str(timeres) + "s"
        else:
            timebin = ""

        #############################
        # Making spectral chunks
        #############################
        coarse_channel_bands = get_MWA_coarse_bands(msname)
        chanlist = []
        for chan in coarse_channel_bands:
            start_chan = chan[0]
            end_chan = chan[-1]
            if end_chan > start_chan:
                chanlist.append(f"{start_chan}~{end_chan-1}")
            elif start_chan == end_chan:
                chanlist.append(f"{start_chan}")

        splited_ms_list = []
        timerange_list = get_timeranges(
            msname,
            time_interval,
            time_window,
            quack_timestamps=quack_timestamps,
        )
        timerange = ",".join(timerange_list)
        for chanrange in chanlist:
            outputvis = f"{workdir}/{prefix}_{os.path.basename(msname).split('.ms')[0]}_spw_{chanrange}.ms"
            if os.path.exists(f"{outputvis}/.splited") and force_split is False:
                print(f"{outputvis} is already splited successfully.")
                splited_ms_list.append(outputvis)
            else:
                if os.path.exists(outputvis):
                    os.system(f"rm -rf {outputvis}")
                if os.path.exists(f"{outputvis}.flagversions"):
                    os.system(f"rm -rf {outputvis}.flagversions")
                splited_ms = single_mstransform(
                    msname=msname,
                    outputms=outputvis,
                    width=chanwidth,
                    timebin=timebin,
                    datacolumn=datacolumn,
                    spw="0:" + chanrange,
                    corr="",
                    timerange=timerange,
                    n_threads=n_threads,
                )
                if os.path.exists(splited_ms):
                    splited_ms_list.append(splited_ms)
        if len(splited_ms_list) == 0:
            print(f"Spliting of measurement set: {msname} is unsuccessful.")
            return 1, []
        else:
            for splited_ms in splited_ms_list:
                drop_cache(splited_ms)
            print(f"Spliting of measurement set: {msname} is done successfully.")
            return 0, splited_ms_list
    except Exception as e:
        traceback.print_exc()
        print(f"Spliting of measurement set: {msname} is unsuccessful.")
        return 1, []
    finally:
        time.sleep(1)
        drop_cache(msname)


def main(
    mslist,
    workdir="",
    datacolumn="data",
    scan=1,
    time_window=-1,
    time_interval=-1,
    quack_timestamps=-1,
    freqres=-1,
    timeres=-1,
    prefix="targets",
    force_split=False,
    cpu_frac=0.8,
    mem_frac=0.8,
    logfile=None,
    jobid=0,
    start_remote_log=False,
    dask_client=None,
):
    """
    Split target scans from a measurement set into smaller chunks for parallel processing.

    Parameters
    ----------
    mslist : str
        Measurement sets (comma separated).
    workdir : str, optional
        Working directory for intermediate and output products. If empty, defaults to `<msname>/workdir`.
    datacolumn : str, optional
        Column of the MS to use for splitting (e.g., "DATA", "CORRECTED"). Default is "data".
    scan : int, optional
        Scan numbers to split.
    time_window : float, optional
        Time window in seconds for a single time chunk. Set -1 to disable. Default is -1.
    time_interval : float, optional
        Time interval in seconds between two time chunks. Set -1 to disable. Default is -1.
    quack_timestamps : int, optional
       Number of timestamps to flag at the beginning and end of each scan ("quack"). -1 to disable. Default is -1.
    freqres : float, optional
        Frequency resolution in MHz for spectral averaging. Set -1 to disable. Default is -1.
    timeres : float, optional
        Time resolution in seconds for time averaging. Set -1 to disable. Default is -1.
    prefix : str, optional
        Prefix for the output split MS files. Default is "targets".
    force_split : bool, optional
        Force to split
    cpu_frac : float, optional
        Fraction of available CPUs to allocate per task. Default is 0.8.
    mem_frac : float, optional
        Fraction of available memory to allocate per task. Default is 0.8.
    logfile : str or None, optional
        Path to log file. If None, logging to file is disabled. Default is None.
    jobid : int, optional
        Job identifier for tracking and PID storage. Default is 0.
    start_remote_log : bool, optional
        If True, enables remote logging using credentials stored in workdir. Default is False.
    dask_client : dask.client, optional
        Dask client

    Returns
    -------
    int
        Success message
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
        time.sleep(5)
        jobname, password = np.load(
            f"{workdir}/jobname_password.npy", allow_pickle=True
        )
        if os.path.exists(logfile):
            observer = init_logger(
                "do_target_split", logfile, jobname=jobname, password=password
            )
    if observer == None:
        print("Remote link or jobname is blank. Not transmiting to remote logger.")

    dask_cluster = None
    if dask_client is None:
        if mem_frac <= 0:
            mem_frac = 0.8
        result = get_local_dask_cluster(
            workdir,
            mem_frac=mem_frac,
        )
        if result is None:
            print("Error occured in creating local cluster.")
            return 1
        else:
            dask_client, dask_cluster, dask_dir = result
        nworker = max(2, int(psutil.cpu_count() * cpu_frac))
        scale_worker_and_wait(dask_cluster, dask_client, nworker)

    try:
        if len(mslist) > 0:
            print("###################################")
            print(f"Start spliting measurement sets in coarse frequency bands.")
            print("###################################")
            ##################################
            # Parallel spliting
            ##################################
            scheduler_name = get_scheduler_name()
            if scheduler_name == "local":
                total_cpu = max(1, int(psutil.cpu_count() * cpu_frac))
                total_mem = (psutil.virtual_memory().available * mem_frac) / (
                    1024**3
                )  # In GB
                njobs = max(1, total_cpu)
                n_threads = max(1, int(total_cpu / njobs))
            else:
                client_info = dask_client.scheduler_info()["workers"]
                njobs = len(client_info)
                n_threads = os.environ.get("OMP_NUM_THREADS")
                if n_threads is not None:
                    n_threads = int(n_threads)
                else:
                    n_threads = 1

            print(f"Total dask worker: {njobs}")
            print(f"CPU per worker: {n_threads}")
            tasks = [
                delayed(split_target_scans)(
                    msname,
                    workdir,
                    float(timeres),
                    float(freqres),
                    datacolumn,
                    time_window=float(time_window),
                    time_interval=float(time_interval),
                    quack_timestamps=int(quack_timestamps),
                    force_split=force_split,
                    scan=scan,
                    prefix=prefix,
                    n_threads=n_threads,
                )
                for msname in mslist
            ]

            future = dask_client.compute(tasks)
            results = dask_client.gather(future)
            succeed = 0
            failed = 0
            for r in results:
                msg, splited_ms = r
                if msg == 0:
                    succeed += 1
                else:
                    failed += 1
            print("########################################")
            print(f"Total measurement sets: {len(mslist)}")
            print(f"Total successful spliting: {succeed}")
            print(f"Total failed spliting: {failed}")
            print("#########################################")
            if failed == len(mslist):
                msg = 1
            else:
                msg = 0
        else:
            print("Please provide correct measurement set list.")
            msg = 1
    except Exception as e:
        traceback.print_exc()
        msg = 1
    finally:
        time.sleep(5)
        for msname in mslist:
            drop_cache(msname)
        drop_cache(workdir)
        clean_shutdown(observer)
        if dask_cluster is not None:
            dask_client.shutdown()
            dask_client.close()
            dask_cluster.close()
            os.system(f"rm -rf {dask_dir}")
        if msg == 0:
            print("All measurement sets are splited successfully.")
        else:
            print("Error occured in spliting measurement sets.")
        return msg


def cli():
    parser = argparse.ArgumentParser(
        description="Split measurement set into coarse channels",
        formatter_class=SmartDefaultsHelpFormatter,
    )

    # Essential parameters
    basic_args = parser.add_argument_group(
        "###################\nEssential parameters\n###################"
    )
    basic_args.add_argument(
        "mslist",
        type=str,
        help="Name of measurement sets (required positional argument)",
    )
    basic_args.add_argument(
        "--workdir",
        type=str,
        default="",
        help="Name of work directory",
    )

    # Advanced parameters
    adv_args = parser.add_argument_group(
        "###################\nAdvanced parameters\n###################"
    )
    adv_args.add_argument(
        "--datacolumn",
        type=str,
        default="data",
        help="Data column to split",
    )
    adv_args.add_argument(
        "--scan",
        type=int,
        default=1,
        help="Target scan to split",
    )
    adv_args.add_argument(
        "--time_window",
        type=float,
        default=-1,
        help="Time window in seconds of a single time chunk",
    )
    adv_args.add_argument(
        "--time_interval",
        type=float,
        default=-1,
        help="Time interval in seconds between two time chunks",
    )
    adv_args.add_argument(
        "--quack_timestamps",
        type=int,
        default=-1,
        help="Time stamps to ignore at the start and end of the each scan",
    )
    adv_args.add_argument(
        "--freqres",
        type=float,
        default=-1,
        help="Frequency to average in MHz",
        metavar="Float",
    )
    adv_args.add_argument(
        "--timeres",
        type=float,
        default=-1,
        help="Time bin to average in seconds",
        metavar="Float",
    )
    adv_args.add_argument(
        "--prefix",
        type=str,
        default="targets",
        help="Splited ms prefix name",
    )
    adv_args.add_argument("--force_split", action="store_true", help="Force to split")
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
        help="CPU fraction to use",
        metavar="Float",
    )
    hard_args.add_argument(
        "--mem_frac",
        type=float,
        default=0.8,
        help="Memory fraction to use",
        metavar="Float",
    )
    hard_args.add_argument("--logfile", type=str, default=None, help="Log file")
    hard_args.add_argument("--jobid", type=int, default=0, help="Job ID")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        return 1

    args = parser.parse_args()

    msg = main(
        mslist=args.mslist,
        workdir=args.workdir,
        datacolumn=args.datacolumn,
        scan=args.scan,
        time_window=args.time_window,
        time_interval=args.time_interval,
        quack_timestamps=args.quack_timestamps,
        force_split=args.force_split,
        freqres=args.freqres,
        timeres=args.timeres,
        prefix=args.prefix,
        cpu_frac=args.cpu_frac,
        mem_frac=args.mem_frac,
        logfile=args.logfile,
        jobid=args.jobid,
        start_remote_log=args.start_remote_log,
    )
    return msg
