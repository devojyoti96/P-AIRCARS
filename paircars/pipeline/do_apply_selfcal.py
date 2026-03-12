import logging
import numpy as np
import argparse
import traceback
import time
import glob
import sys
import os
from casatools import msmetadata
from dask import delayed
from astropy.io import fits
from paircars.utils.logger_utils import (
    SmartDefaultsHelpFormatter,
    clean_shutdown,
    init_logger,
)
from paircars.utils.ms_metadata import check_datacolumn_valid, get_ms_size
from paircars.utils.mwa_utils import freq_to_MWA_coarse, get_ncoarse
from paircars.utils.proc_manage_utils import (
    scale_worker_and_wait,
    get_local_dask_cluster,
)
from paircars.utils.resource_utils import drop_cache
from paircars.pipeline.do_apply_basiccal import applysol

logging.getLogger("distributed").setLevel(logging.ERROR)
logging.getLogger("tornado.application").setLevel(logging.CRITICAL)


def run_all_applysol(
    mslist,
    metafits,
    dask_client,
    workdir,
    caldir,
    overwrite_datacolumn=False,
    applymode="calonly",
    force_apply=False,
    cpu_frac=0.8,
    mem_frac=0.8,
):
    """
    Apply self-calibrator solutions on all target scans

    Parameters
    ----------
    mslist : list
        Measurement set list
    metafits: str
        Metafits file
    dask_client : dask.client
        Dask client
    workdir : str
        Working directory
    caldir : str
        Calibration directory
    overwrite_datacolumn : bool, optional
        Overwrite data column or not
    applymode : str, optional
        Apply mode
    force_apply : bool, optional
        Force to apply solutions even already applied
    cpu_frac : float, optional
        CPU fraction to use
    mem_frac : float, optional
        Memory fraction to use

    Returns
    --------
    int
        Succeeded gain solution ms number
    int
        Failed gain solution ms number
    int
        Succeeded polarisation solution ms number
    int
        Failed polarisation solution ms number
    """
    cpu_frac = min(0.8, abs(cpu_frac))
    mem_frac = min(0.8, abs(mem_frac))

    if len(mslist) == 0:
        print("Please provide a valid measurement set list.")
        return 0, 0, 0, 0
    else:
        gain_succeed = 0
        gain_failed = len(mslist)
        pol_succeed = 0
        pol_failed = len(mslist)

    try:
        os.chdir(workdir)
        mslist = np.unique(mslist).tolist()
        header = fits.getheader(metafits)
        obsid = header["GPSTIME"]
        selfcal_tables = sorted(glob.glob(f"{caldir}/selfcal_{obsid}_coarsechan*.gcal"))
        sorted(glob.glob(f"{caldir}/selfcal_{obsid}_coarsechan*.bcal"))
        sorted(glob.glob(f"{caldir}/selfcal_{obsid}_coarsechan*.dcal"))
        if len(selfcal_tables) == 0:
            print(f"No self-cal caltable is present in {caldir}.")
            return gain_succeed, gain_failed, pol_succeed, pol_failed
        selfcal_tables_start_chans = np.array(
            [
                int(
                    os.path.basename(i)
                    .split(".gcal")[0]
                    .split("coarsechan_")[-1]
                    .split("_")[0]
                )
                for i in selfcal_tables
            ]
        )
        selfcal_tables_end_chans = np.array(
            [
                int(
                    os.path.basename(i)
                    .split(".gcal")[0]
                    .split("coarsechan_")[-1]
                    .split("_")[-1]
                )
                for i in selfcal_tables
            ]
        )
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
            print("No valid measurement set.")
            return gain_succeed, gain_failed, pol_succeed, pol_failed

        ####################################
        # Applycal jobs
        ####################################
        print(f"Total ms list: {len(mslist)}")
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

        tasks = []
        msmd = msmetadata()
        for ms in mslist:
            msmd.open(ms)
            freqs = msmd.chanfreqs(0, unit="MHz")
            start_freq = np.nanmin(freqs)
            end_freq = np.nanmax(freqs)
            start_coarse_chan = freq_to_MWA_coarse(start_freq)
            end_coarse_chan = freq_to_MWA_coarse(end_freq)
            msmd.close()
            gaintable = []
            quartical_table = []
            for i in range(len(selfcal_tables_start_chans)):
                s = selfcal_tables_start_chans[i]
                e = selfcal_tables_end_chans[i]
                if start_coarse_chan >= s and end_coarse_chan <= e:
                    gaintable = selfcal_tables[i]
                    gaintable_prefix = gaintable.split(".gcal")[0]
                    gaintable = [f"{gaintable_prefix}.gcal", f"{gaintable_prefix}.bcal"]
                    quartical_table = [f"{gaintable_prefix}.dcal"]

            if len(gaintable) == 0:
                print(
                    f"Measurement set coarse channel : {start_coarse_chan} to {end_coarse_chan}. Corresponding self-calibration table is not present."
                )
            else:
                if len(quartical_table) == 0:
                    print(
                        f"Measurement set coarse channel : {start_coarse_chan} to {end_coarse_chan}. Corresponding polarisation self-calibration table is not present."
                    )
                    os.system(f"touch {ms}/.nopolselfcal")
                tasks.append(
                    delayed(applysol)(
                        ms,
                        workdir,
                        gaintable=gaintable,
                        quartical_table=quartical_table,
                        overwrite_datacolumn=overwrite_datacolumn,
                        applymode=applymode,
                        interp=["linear,linear"],
                        n_threads=n_threads,
                        mem_limit=mem_limit,
                        force_apply=force_apply,
                        soltype="selfcal",
                    )
                )
        if len(tasks) > 0:
            print("Applying solutions...")
            results = list(dask_client.gather(dask_client.compute(tasks)))
            gain_msg = []
            pol_msg = []
            for r in results:
                gain_msg.append(r[0])
                pol_msg.append(r[1])
            gain_failed = sum(gain_msg)
            pol_failed = sum(pol_msg)
            gain_succeed = len(mslist) - gain_failed
            pol_succeed = len(mslist) - pol_failed
            print("##################")
            print(f"Total measurement sets: {len(mslist)}")
            print(f"Gain solution applied, Succeeded: {gain_succeed}")
            print(f"Gain solution applied, Failed: {gain_failed}")
            print(f"Polarisation solution applied, Succeeded: {pol_succeed}")
            print(f"Polarisation solution applied, Failed: {pol_failed}")
            if gain_failed == 0 and pol_failed == 0:
                print(
                    "Applying gain and polarisation self-calibration solutions for targets are done successfully."
                )
            elif pol_failed == 0:
                print(
                    "Applying gain self-calibration solutions for targets are done successfully, but failed for polarisation solutions."
                )
            else:
                print(
                    "Applying self-calibration solutions for targets are not done successfully."
                )
            print("##################")
        else:
            print("##################")
            print(
                "Applying self-calibration solutions for targets are not done successfully. No suitable calibration solutions are found."
            )
            print("##################")
    except Exception as e:
        traceback.print_exc()
        print("##################")
        print(
            "Applying self-calibration solutions for targets are not done successfully."
        )
        print("##################")
    finally:
        return gain_succeed, gain_failed, pol_succeed, pol_failed


def main(
    mslist,
    metafits,
    workdir,
    caldir,
    applymode="calonly",
    overwrite_datacolumn=False,
    force_apply=False,
    start_remote_log=False,
    cpu_frac=0.8,
    mem_frac=0.8,
    logfile=None,
    jobid=0,
    dask_client=None,
):
    """
    Apply calibration solutions to a list of measurement sets.

    Parameters
    ----------
    mslist : str
        Comma-separated list of measurement set paths to apply calibration to.
    metafits : str
        Metafits file
    workdir : str
        Directory for logs, intermediate files, and PID tracking.
    caldir : str
        Directory containing calibration tables (e.g., gain, bandpass, polarization).
    applymode : str, optional
        Mode for applying calibration (e.g., "calonly", "calflag", "flagonly"). Default is "calonly".
    overwrite_datacolumn : bool, optional
        If True, overwrites the existing corrected data column. Default is False.
    force_apply : bool, optional
        If True, applies calibration even if it appears to have been applied already. Default is False.
    start_remote_log : bool, optional
        Whether to enable remote logging using credentials found in `workdir`. Default is False.
    cpu_frac : float, optional
        Fraction of available CPU resources to allocate per task. Default is 0.8.
    mem_frac : float, optional
        Fraction of system memory to allocate per task. Default is 0.8.
    logfile : str or None, optional
        Path to the logfile for saving logs. If None, logging to file is disabled. Default is None.
    jobid : int, optional
        Job ID for PID tracking and logging. Default is 0.
    dask_client : dask.client, optional
        Dask client address

    Returns
    --------
    int
        Succeeded gain solution ms number
    int
        Failed gain solution ms number
    int
        Succeeded polarisation solution ms number
    int
        Failed polarisation solution ms number
    """
    cpu_frac = min(0.8, abs(cpu_frac))
    mem_frac = min(0.8, abs(mem_frac))

    # Get first MS from mslist for fallback directory creation
    mslist = mslist.split(",")

    if workdir == "":
        workdir = os.path.dirname(os.path.abspath(mslist[0])) + "/workdir"
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)

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
                "apply_selfcal", logfile, jobname=jobname, password=password
            )
    if observer is None:
        print("Remote link or jobname is blank. Not transmiting to remote logger.")

    if len(mslist) == 0:
        print("Please provide a valid measurement set list.")
        return 0, 0, 0, 0
    else:
        gain_succeed = 0
        gain_failed = len(mslist)
        pol_succeed = 0
        pol_failed = len(mslist)

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
            return 0, 0, 0, 0
        else:
            dask_client, dask_cluster, dask_dir, nworker = result
        scale_worker_and_wait(dask_cluster, dask_client, nworker)

    try:
        print("###################################")
        print("Starting applying solutions...")
        print("###################################")

        if caldir == "" or not os.path.exists(caldir):
            print("Provide existing caltable directory.")
        else:
            gain_succeed, gain_failed, pol_succeed, pol_failed = run_all_applysol(
                mslist,
                metafits,
                dask_client,
                workdir,
                caldir,
                overwrite_datacolumn=overwrite_datacolumn,
                applymode=applymode,
                force_apply=force_apply,
                cpu_frac=cpu_frac,
                mem_frac=mem_frac,
            )
    except Exception:
        traceback.print_exc()
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
            os.system(f"rm -rf {dask_dir}")
    return gain_succeed, gain_failed, pol_succeed, pol_failed


def cli():
    parser = argparse.ArgumentParser(
        description="Apply self-calibration solutions to target scans",
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
        "metafits",
        type=str,
        help="Metafits file (required)",
    )
    basic_args.add_argument(
        "--workdir",
        type=str,
        default="",
        required=True,
        help="Working directory for intermediate files",
    )
    basic_args.add_argument(
        "--caldir",
        type=str,
        default="",
        required=True,
        help="Directory containing self-calibration tables",
    )

    # Advanced parameters
    adv_args = parser.add_argument_group(
        "###################\nAdvanced parameters\n###################"
    )
    adv_args.add_argument(
        "--applymode",
        type=str,
        default="calonly",
        help="Applycal mode (e.g. 'calonly', 'calflag')",
    )
    adv_args.add_argument(
        "--overwrite_datacolumn",
        action="store_true",
        help="Overwrite corrected data column in MS",
    )
    adv_args.add_argument(
        "--force_apply",
        action="store_true",
        help="Force apply calibration even if already applied",
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
    hard_args.add_argument(
        "--logfile", type=str, default=None, help="Optional path to log file"
    )
    hard_args.add_argument(
        "--jobid", type=str, default="0", help="Job ID for logging and PID tracking"
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        return 1

    args = parser.parse_args()

    gain_succeed, gain_failed, pol_succeed, pol_failed = main(
        args.mslist,
        args.metafits,
        args.workdir,
        args.caldir,
        applymode=args.applymode,
        overwrite_datacolumn=args.overwrite_datacolumn,
        force_apply=args.force_apply,
        start_remote_log=args.start_remote_log,
        cpu_frac=float(args.cpu_frac),
        mem_frac=float(args.mem_frac),
        logfile=args.logfile,
        jobid=args.jobid,
    )
    if gain_failed == 0 and pol_failed == 0:
        msg = 0
    else:
        msg = 1
    return msg
