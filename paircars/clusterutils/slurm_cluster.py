import types
import dask
import gc
import time
import os
import subprocess
import sys
import traceback
import logging
import shlex
import re
from dask.distributed import Client
from dask_jobqueue import SLURMCluster
from paircars.utils.basic_utils import *
from paircars.utils.proc_manage_utils import *


def get_slurm_node_resources(partition=None, cpu_frac=0.8, mem_frac=0.8):
    """
    Get node resources for SLURM cluster

    Parameters
    ----------
    partition : str, optional
        Partition name
    cpu_frac : float, optional
        CPU fraction to use
    mem_frac : float, optional
        Memory fraction to use

    Returns
    -------
    int
        Number of CPU threads
    float
        Memory in GB
    """
    if partition is not None:
        cmd = ["sinfo", "-h", "-p", partition, "-o", "%c %m"]
    else:
        cmd = ["sinfo", "-h", "-o", "%c %m"]
    out = subprocess.check_output(cmd).decode().strip().split("\n")
    cores = []
    mems = []
    for line in out:
        c, m = line.split()
        cores.append(int(c.rstrip("+")))
        mems.append(int(m.rstrip("+")))
    total_cpu = min(cores)
    total_mem = min(mems) / (1024)  # In GB
    cpu_frac = min(0.8, cpu_frac)
    mem_frac = min(0.8, mem_frac)
    ncpu = max(1, int(total_cpu * cpu_frac))
    mem = round(total_mem * mem_frac, 1)
    return ncpu, mem


def get_slurm_dask_cluster(
    dask_dir,
    jobid=None,
    cpu_frac=0.8,
    mem_frac=0.8,
    max_mem=16,
    partition=None,
    account=None,
    walltime="24:00:00",
    python_path=None,
    spill_frac=0.7,
    verbose=True,
):
    """
    Launch a SLURMCluster using a YAML configuration and return a connected Dask client.

    Parameters
    ----------
    dask_dir : str
        Dask working directory (for temporary files)
    jobid : int
        JobID of P-AIRCARS to avoid mixup of cluster configurration with other P-AIRCARS jobs.
    cpu_frac : float, optional
        CPU fraction to use
    mem_frac : float, optional
        Memory fraction to use
    max_mem : float, optional
        Maximum job memory in GB
    partition : str, optional
        SLURM partition name
        Note: If your cluster requires this, you should provide. Otherwise, error will occur.
    account : str, optional
        SLURM account name
        Note: If your cluster requires this, you should provide. Otherwise, error will occur.
    walltime : str, optional
        Job walltime, maximum time the SLURM job can run (HH:MM:SS)
    spill_frac : float
        Fraction of memory to spill to disk
    verbose : bool
        Print Dask dashboard URL and diagnostics

    Returns
    -------
    client : dask.distributed.Client
        Connected Dask client
    cluster : dask_jobqueue.SLURMCluster
        SLURM Dask cluster
    str
        Dask directory used
    """
    logging.getLogger("distributed").setLevel(logging.ERROR)
    scheduler_name = get_scheduler_name()
    if scheduler_name != "slurm":
        print("SLURM is not avilable as job scheduler in your cluster.")
        return

    cpu_frac = min(0.8, cpu_frac)
    mem_frac = min(0.8, mem_frac)

    if jobid is None:
        jobid = get_jobid()

    os.makedirs(dask_dir, exist_ok=True)
    log_dir = f"{dask_dir}/slurm_log_{jobid}"
    os.makedirs(log_dir, exist_ok=True)

    dask_dir = os.path.join(dask_dir.rstrip("/"), f"dask_{int(time.time())}")
    dask_dir_tmp = os.path.join(dask_dir, "tmp")
    os.makedirs(dask_dir_tmp, exist_ok=True)

    try:
        dask.config.set(
            {
                "temporary-directory": dask_dir_tmp,
                "distributed.worker.memory.target": spill_frac,
                "distributed.worker.memory.spill": spill_frac + 0.1,
                "distributed.worker.memory.pause": spill_frac + 0.2,
                "distributed.worker.memory.terminate": spill_frac + 0.25,
            }
        )
        ncpu, mem = get_slurm_node_resources(
            partition=partition, cpu_frac=cpu_frac, mem_frac=mem_frac
        )
        if python_path is None:
            python_path = sys.executable
        interface = detect_best_interface()

        job_extra = [
            f"--nodes=1",
            f"--ntasks=1",
            f"--cpus-per-task={ncpu}",
            f"--mem={min(max_mem,mem)}G",
            f"--exclusive",
            f"--output={log_dir}/paircars_{jobid}-%j.out",
            f"--error={log_dir}/paircars_{jobid}-%j.err",
        ]

        cluster = SLURMCluster(
            queue=partition,
            account=account,
            cores=ncpu,
            n_workers=1,
            walltime=walltime,
            memory=f"{min(max_mem,mem)}G",
            processes=1,
            interface=interface,
            python=python_path,
            local_directory=dask_dir_tmp,
            death_timeout=60,
            log_directory=log_dir,
            name=f"paircars_{jobid}",
            shared_temp_directory=dask_dir_tmp,
            env_extra=[
                "OMP_NUM_THREADS=1",
                "MKL_NUM_THREADS=1",
                "OPENBLAS_NUM_THREADS=1",
                "NUMEXPR_NUM_THREADS=1",
                "MALLOC_TRIM_THRESHOLD_=0",
                f"TMPDIR={dask_dir_tmp}",
                f"TMP={dask_dir_tmp}",
                f"TEMP={dask_dir_tmp}",
                f"DASK_TEMPORARY_DIRECTORY={dask_dir_tmp}",
                "PYTHONWARNINGS=ignore::UserWarning:contextlib",
            ],
        )

        cluster.scale(1)
        client = Client(cluster, heartbeat_interval="5s")
        client.run_on_scheduler(gc.collect)
        if verbose:
            print("####################################################")
            print(f"Dask dashboard available at: {client.dashboard_link}")
            print(f"CPU per worker: {ncpu}")
            print(f"Memory per worker: {min(max_mem,usable_mem)}GB")
            print("####################################################")

        return client, cluster, dask_dir
    except Exception as e:
        print("Error occured in creating SLURM cluster.")
        traceback.print_exc()
        os.system(f"rm -rf {output_path} {log_dir} {dask_dir}")


def slurm_time_to_seconds(timestr):
    """
    Convert SLURM time format (D-HH:MM:SS or HH:MM:SS) to seconds.

    Parameters
    ----------
    timestr : str
        Time string in SLURM format

    Returns
    -------
    float
        Time in seconds
    """
    if timestr.lower() in ["infinite", "unlimited"]:
        return float("inf")
    if "-" in timestr:
        days, hms = timestr.split("-")
        h, m, s = map(int, hms.split(":"))
        return int(days) * 86400 + h * 3600 + m * 60 + s
    else:
        h, m, s = map(int, timestr.split(":"))
        return h * 3600 + m * 60 + s


def get_max_walltime(partition):
    """
    Get maximum wall time for the partition

    Parameters
    ----------
    partition : str
        Partition name

    Returns
    -------
    str
        Maximum wall time
    """
    result = subprocess.run(
        ["scontrol", "show", "partition"],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError("Failed to query SLURM partitions.")
    output = result.stdout
    partitions = {}
    blocks = output.split("\n\n")
    for block in blocks:
        name_match = re.search(r"PartitionName=(\S+)", block)
        time_match = re.search(r"MaxTime=(\S+)", block)
        if name_match and time_match:
            part_name = name_match.group(1)
            max_time = time_match.group(1)
            partitions[part_name] = max_time
    if partition not in partitions:
        raise ValueError(f"Partition {partition} not found.")
    max_time = partitions[partition]
    return max_time, slurm_time_to_seconds(max_time)


def submit_master_flow(args, jobid):
    """
    Submit P-AIRCARS master flow to a slurm job

    Parameters
    ----------
    args : dict
        Arparser dictionary
    jobid : int
        P-AIRCARS jobid

    Returns
    -------
    int
        Success message
    """
    scheduler_name = get_scheduler_name()
    if scheduler_name is not "slurm":
        print("SLURM job scheduler is not available.")
        return 1
    cli_cmd = " ".join(shlex.quote(arg) for arg in sys.argv[1:])
    if args.partition and args.partition is not None:
        max_time, max_time_seconds = get_max_walltime(args.partition)
    else:
        print("Please provide partition name to run SLURM jobs.")
        return 1

    try:
        #################################
        # Determining wall time
        #################################
        if args.walltime is None:
            walltime = max_time
        else:
            wall_time_second = slurm_time_to_seconds(args.walltime)
            if wall_time_seconod > max_time_second:
                print(
                    f"Walltime : {args.walltime} is larger than maximum allowed time: {max_time}."
                )
                walltime = max_time
            else:
                walltime = args.walltime
        #############################
        # Determining cpu and memory
        #############################
        ncpu, mem = get_slurm_node_resources(
            partition=args.partition, cpu_frac=args.cpu_frac, mem_frac=args.mem_frac
        )
        script = f"""#!/bin/bash
        #SBATCH --job-name=paircars_{jobid}
        #SBATCH --time={walltime}
        #SBATCH --output={args.workdir}/paircars_{jobid}_%j.out
        #SBATCH --output={args.workdir}/paircars_{jobid}_%j.err
        #SBATCH --partition={args.partition}
        #SBATCH --partition={args.partition}
        #SBATCH --nodes=1
        #SBATCH --ntasks=1
        #SBATCH --cpus-per-task={min(8,ncpu)}
        #SBATCH --mem={min(16,mem)}G
        """
        if args.account:
            script += f"#SBATCH --account={args.account}\n"
        os.makedirs(args.workdir, exist_ok=True)
        script_path = os.path.join(args.workdir, f"paircars_slurm_{jobid}.sh")
        with open(script_path, "w") as f:
            f.write(script)
        subprocess.run(["sbatch", script_path], check=True)
        return 0
    except Exception as e:
        traceback.print_exc()
        return 1


# Exposing only functions
__all__ = [
    name
    for name, obj in globals().items()
    if isinstance(obj, types.FunctionType) and obj.__module__ == __name__
]
