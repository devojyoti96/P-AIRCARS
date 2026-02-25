import os
import glob
import numpy as np
import time
import sys
import dask
import psutil
import traceback
import logging
import argparse
import subprocess
from dask import delayed
from casatasks import setjy
from casatools import table as casatable, msmetadata
from paircars.utils.basic_utils import suppress_output, get_datadir
from paircars.utils.logger_utils import (
    SmartDefaultsHelpFormatter,
    clean_shutdown,
    init_logger,
)
from paircars.utils.ms_metadata import get_ms_size
from paircars.utils.proc_manage_utils import (
    scale_worker_and_wait,
    get_local_dask_cluster,
    get_scheduler_name,
)
from paircars.utils.resource_utils import drop_cache
from paircars.utils.udocker_utils import (
    run_hyperdrive,
    initialize_hyperdrive_container,
    check_udocker_container,
)


logging.getLogger("distributed").setLevel(logging.ERROR)
logging.getLogger("tornado.application").setLevel(logging.CRITICAL)
datadir = get_datadir()


def import_hyperdrive_model(
    msname, metafits, beamfile="", sourcelist="", ncpu=1, verbose=False
):
    """
    Simulate visibilities and import in the measurement set

    Parameters
    ----------
    msname : str
        Name of the measurement set
    metafits : str
        Name of the metafits file
    beamfile : str, optional
        Beam file name
    sourcelist : str, optional
        Source file name
    ncpu : int, optional
        Number of cpu threads to use
    verbose : bool, optional
        Verbose output or not
    """
    ncpu = max(1, ncpu)

    msname = msname.rstrip("/")
    msname = os.path.abspath(msname)
    print(
        "#######################\nImporting model for ms:"
        + msname
        + "\n###################\n"
    )
    if beamfile == "" or os.path.exists(beamfile) is False:
        with suppress_output():
            msmd = msmetadata()
            msmd.open(msname)
            freqres = msmd.chanres(0, unit="kHz")[0]
            msmd.close()
        beam_files = glob.glob(f"{datadir}/mwa_full_embedded_element_pattern*.h5")
        beam_files_freqs = []
        for beamfile in beam_files:
            if os.path.basename(beamfile) == "mwa_full_embedded_element_pattern.h5":
                beam_file_freq = 1280.0
            else:
                beam_file_freq = float(
                    os.path.basename(beamfile)
                    .split(".h5")[0]
                    .split("mwa_full_embedded_element_pattern_")[-1]
                )
            beam_files_freqs.append(beam_file_freq)
        beam_files_freqs = np.array(beam_files_freqs)
        pos = np.argmin(np.abs(beam_files_freqs - freqres))
        beamfile = beam_files[pos]
    if sourcelist == "" or os.path.exists(sourcelist) is not True:
        sourcelist = f"{datadir}/GGSM.txt"
    model_msname = msname.split(".ms")[0] + "_model.ms"
    try:
        starttime = time.time()
        with suppress_output():
            msmd = msmetadata()
            msmd.open(msname)
            nchan = msmd.nchan(0)
            mid_freq = msmd.meanfreq(0, unit="MHz")
            freqres = msmd.chanres(0, unit="kHz")[0]
            npol = msmd.ncorrforpol()[0]
            nant = msmd.nantennas()
            times = msmd.timesforfield(0)
            ntime = len(times)
            timeres = msmd.exposuretime(scan=1)["value"]
            nrow = msmd.nrows()
            msmd.close()

        hyperdrive_cmd_args = [
            f"hyperdrive",
            "vis-simulate",
            "-m",
            metafits,
            "--beam-file",
            beamfile,
            "--middle-freq",
            str(mid_freq),
            "--freq-res",
            str(freqres),
            "--time-res",
            str(timeres),
            "--source-dist-cutoff",
            "180",
            "-s",
            sourcelist,
            "-n",
            "2000",
            "--output-model-files",
            f"{model_msname}",
            "--output-model-freq-average",
            f"{freqres}kHz",
            "--num-fine-channels",
            str(nchan),
            "--num-timesteps",
            str(ntime),
            "--output-model-time-average",
            f"{timeres}s",
        ]
        hyperdrive_cmd = " ".join(hyperdrive_cmd_args)
        result = run_hyperdrive(hyperdrive_cmd, ncpu=ncpu, verbose=verbose)
        if result != 0:
            print("Error occured in hyperdrive.")
            return 1

        ########################
        # Importing model
        ########################
        with suppress_output():
            data_table = casatable()
            data_table.open(msname, nomodify=False)
            column_names = data_table.colnames()
            if "MODEL_DATA" not in column_names:
                data_table.close()
                setjy(
                    vis=msname,
                    standard="manual",
                    fluxdensity=[1, 0, 0, 0],
                    usescratch=True,
                )
                data_table.open(msname, nomodify=False)
            model_table = casatable()
            model_table.open(model_msname, nomodify=False)
            baselines = [
                *zip(data_table.getcol("ANTENNA1"), data_table.getcol("ANTENNA2"))
            ]
            m_array = model_table.getcol("DATA")
            pos = np.array([i[0] != i[1] for i in baselines])
            model_array = np.empty((npol, nchan, len(baselines)), dtype="complex")
            model_array[..., pos] = m_array
            model_array[..., ~pos] = 0.0
            data_table.putcol("MODEL_DATA", model_array)
            data_table.close()
            model_table.close()
        del m_array, model_array
        print(f"Model import done in: {round(time.time()-starttime,2)}s")
        return 0
    except Exception as e:
        print(f"Model simulation and import failed for: {msname}.")
        traceback.print_exc()
        return 1
    finally:
        os.system(f"rm -rf {model_msname}")


def run_all_modeling(
    mslist, dask_client, metafits, beamfile, sourcelist, ncpu, verbose
):
    """
    Run all modeling

    Parameters
    ----------
    mslist : list
        Measurement set list
    dask_client : dask. client
        Dask client
    metafits : str
        Metafits file
    beamfile : str
        MWA primary beam file
    sourcelist : str
        Source list file
    ncpu : int
        Number of CPU threads
    verbose : bool
        Verbose output

    Returns
    -------
    int
        Total failure
    """
    msg = 0
    try:
        if len(mslist) > 0:
            tasks = []
            for msname in mslist:
                tasks.append(
                    delayed(import_hyperdrive_model)(
                        msname,
                        metafits,
                        beamfile=beamfile,
                        sourcelist=sourcelist,
                        ncpu=ncpu,
                        verbose=verbose,
                    )
                )
            print("Start import modeling...")
            futures = dask_client.compute(tasks)
            results = dask_client.gather(futures)
            for i in range(len(results)):
                if results[i] != 0:
                    print(f"Error in model import for ms: {mslist[i]}.")
                    msg += 1
        else:
            print("Please provide a valid measurement set list.")
            msg = -1
    except Exception as e:
        traceback.print_exc()
        msg = -1
    return msg


def main(
    mslist,
    metafits,
    workdir,
    beamfile="",
    sourcelist="",
    verbose=False,
    cpu_frac=0.8,
    mem_frac=0.8,
    logfile=None,
    jobid="0",
    start_remote_log=False,
    dask_client=None,
):
    """
    Make dynamic spectra

    Parameters
    ----------
    mslist : str
        Measurement set list (comma separated)
    metafits : str
        Metafits file
    workdir : str
        Work directory
    beamfile : str, optional
        MWA beam file
    sourcelist : str, optional
        MWA global sky model (fits or ascii in wsclean format)
    verbose : bool, optional
        Verbose output or not
    cpu_frac : float, optional
        CPU fraction
    mem_frac : float, optional
        Memory fraction
    logfile : str, optional
        Log file
    jobid : str, optional
        Job ID
    start_remote_log : bool, optional
        Start remote log
    dask_client: dask.client, optional
        Dask client

    Returns
    -------
    int
        Success messsage
    """
    cpu_frac = min(0.8, cpu_frac)
    mem_frac = min(0.8, mem_frac)

    mslist = mslist.split(",")

    if workdir == "":
        workdir = os.path.dirname(os.path.abspath(mslist[0])) + "/workdir"
    os.makedirs(workdir, exist_ok=True)

    ###########################
    # Hyperdrive container
    ###########################
    container_name = "paircarshyperdrive"
    container_present = check_udocker_container(container_name)
    if not container_present:
        container_name = initialize_hyperdrive_container(name=container_name)
        if container_name is None:
            print(
                f"Container {container_name} is not initiated. First initiate container and then run."
            )
            return 1

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
                "ds_plot", logfile, jobname=jobname, password=password
            )
    if observer == None:
        print("Remote link or jobname is blank. Not transmiting to remote logger.")

    if dask_client is None:
        scheduler_name = "local"
    else:
        scheduler_name = get_scheduler_name()

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
        nworker = min(len(mslist), int(psutil.cpu_count() * cpu_frac) - 1)
        scale_worker_and_wait(dask_cluster, dask_client, nworker + 1)

    #################################################
    # Number of jobs in local and cluster environment
    ##################################################
    if scheduler_name == "local":
        ms_sizes = [get_ms_size(ms) for ms in mslist]
        per_job_mem = 2 * max(ms_sizes)
        mem_limit = (psutil.virtual_memory().available * mem_frac) / (1024**3)
        max_njobs = int(mem_limit / per_job_mem)
        njobs = max(1, min(max_njobs, len(mslist)))
        ncpu = max(1, int(psutil.cpu_count() * cpu_frac / njobs))
    else:
        client_info = dask_client.scheduler_info()["workers"]
        njobs = len(client_info)
        worker_mem_list = []
        for addr, w in client_info.items():
            worker_mem_list.append(w["memory_limit"] / 1024**3)
        mem_limit = round(min(worker_mem_list), 3)
        ncpu = os.environ.get("OMP_NUM_THREADS")
        if ncpu is not None:
            ncpu = int(ncpu)
        else:
            ncpu = 1

    print("#################################")
    print(f"Total dask worker: {njobs}")
    print(f"CPU per worker: {ncpu}")
    print(f"Memory per worker: {mem_limit} GB")
    print("#################################")

    try:
        msg = run_all_modeling(
            mslist, dask_client, metafits, beamfile, sourcelist, ncpu, verbose
        )
        if msg < 0:
            msg = 1
        elif msg > 0:
            print(f"Total model import failure: {msg}")
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
    return msg


################################
# CLI interface
################################
def cli():
    parser = argparse.ArgumentParser(description="Simulate and import MWA visibilities")

    # Essential parameters
    basic_args = parser.add_argument_group(
        "###################\nEssential parameters\n###################"
    )
    basic_args.add_argument(
        "mslist",
        type=str,
        help="Name of the measurement sets (comma seperated)",
    )
    basic_args.add_argument(
        "metafits",
        type=str,
        help="Name of the metafits file",
    )
    basic_args.add_argument(
        "--workdir",
        type=str,
        required=True,
        help="Work directory",
    )

    # Advanced parameters
    adv_args = parser.add_argument_group(
        "###################\nAdvanced parameters\n###################"
    )
    adv_args.add_argument(
        "--beamfile",
        type=str,
        default="",
        help="Name of the MWA PB file",
    )
    adv_args.add_argument(
        "--sourcelist",
        type=str,
        default="",
        help="Source model file",
    )
    adv_args.add_argument("--verbose", action="store_true", help="Verbose output")
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
        help="CPU fraction",
    )
    hard_args.add_argument(
        "--mem_frac",
        type=float,
        default=0.8,
        help="Memory fraction",
    )
    hard_args.add_argument("--logfile", type=str, default=None, help="Log file")
    hard_args.add_argument("--jobid", type=int, default=0, help="Job ID")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        return 1

    args = parser.parse_args()

    msg = main(
        args.mslist,
        args.metafits,
        args.workdir,
        beamfile=args.beamfile,
        sourcelist=args.sourcelist,
        verbose=args.verbose,
        cpu_frac=float(args.cpu_frac),
        mem_frac=float(args.mem_frac),
        logfile=args.logfile,
        jobid=args.jobid,
        start_remote_log=args.start_remote_log,
    )
    return msg


if __name__ == "__main__":
    result = cli()
    print(
        "\n###################\Visibility simulation is finished.\n###################\n"
    )
    os._exit(result)
