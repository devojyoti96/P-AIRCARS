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
from paircars.utils.basic_utils import get_cachedir
from paircars.utils.proc_manage_utils import (
    get_scheduler_name,
    detect_best_interface,
    get_jobid,
)


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
    max_worker=1,
    partition=None,
    account=None,
    walltime="24:00:00",
    python_path=None,
    spill_frac=0.7,
    env=[],
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
    max_worker : float, optional
        Maximum number of worker
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
    env : list, optional
        List of environment variables
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
    log_dir = f"{dask_dir}/slurm_logs"
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
        ncpu = max(1, int(ncpu/max_worker))
        
        if python_path is None:
            python_path = sys.executable
        interface = detect_best_interface()
        mem_limit = round(min(max_mem, mem), 2)
        job_extra = [
            f"--nodes=1",
            f"--ntasks=1",
            f"--cpus-per-task={ncpu}",
            f"--mem={mem_limit}G",
            f"--output={log_dir}/paircars_{jobid}-%j.out",
            f"--error={log_dir}/paircars_{jobid}-%j.err",
        ]

        cluster = SLURMCluster(
            queue=partition,
            account=account,
            cores=1,
            n_workers=1,
            walltime=walltime,
            memory=f"{min(max_mem,mem)}G",
            processes=1,
            interface=interface,
            python=python_path,
            local_directory=dask_dir_tmp,
            death_timeout=300,
            log_directory=log_dir,
            name=f"paircars_{jobid}",
            shared_temp_directory=dask_dir_tmp,
            env_extra=[
                "PYTHONUNBUFFERED=1",
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
            ]
            + env,
        )

        cluster.scale(1)
        client = Client(cluster, heartbeat_interval="5s")
        client.run_on_scheduler(gc.collect)
        if verbose:
            print("####################################################")
            print(f"Dask dashboard available at: {client.dashboard_link}")
            print(f"CPU per worker: {ncpu}")
            print(f"Memory per worker: {mem_limit}GB")
            print("####################################################")
        return client, cluster, dask_dir
    except Exception as e:
        print("Error occured in creating SLURM cluster.")
        traceback.print_exc()
        os.system(f"rm -rf {log_dir} {dask_dir}")
        return


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


def submit_slurm_master_flow(args, jobid):
    """
    Submit P-AIRCARS master flow to a slurm cluster

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
    cli_cmd = (
        "run-mwa-masterflow "
        + " ".join(shlex.quote(arg) for arg in sys.argv[1:])
        + f" --jobid {jobid}"
    )
    if hasattr(args, "partition") and args.partition is not None:
        max_time, max_time_seconds = get_max_walltime(args.partition)
    else:
        print("Please provide partition name to run SLURM jobs.")
        return 1
    if hasattr(args, "workdir") and args.workdir is not None:
        os.makedirs(args.workdir, exist_ok=True)
    else:
        print("Please provide a work directory.")
        return 1

    cachedir = f"{get_cachedir()}/prefect_{scheduler_name}"

    prefect_env_list = [
        f"PREFECT_HOME={cachedir}/prefect_home",
        "PREFECT_API_MODE=server",
        f"PREFECT_API_DATABASE_CONNECTION_URL=sqlite+aiosqlite:///{cachedir}/prefect_home/prefect.db",
        "PREFECT_SERVER_ALLOW_EPHEMERAL_MODE=false",
        "PREFECT_API_URL=http://127.0.0.1:4260/api",
        f"PREFECT_PROFILE=paircarspipe_{scheduler_name}",
        f"PREFECT_PROFILES_PATH={cachedir}/prefect_home/profiles.toml",
        f"PREFECT_LOCAL_STORAGE_PATH={cachedir}/prefect_home/storage",
        f"PREFECT_LOGGING_SETTINGS_PATH={cachedir}/prefect_home/logging.yml",
        f"PREFECT_MEMO_STORE_PATH={cachedir}/prefect_home/memo_store.toml",
    ]

    try:
        #################################
        # Determining wall time
        #################################
        if hasattr(args, "walltime"):
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
        else:
            walltime = max_time
        #############################
        # Determining cpu and memory
        #############################
        if hasattr(args, "cpu_frac") is False:
            cpu_frac = 0.8
        else:
            cpu_frac = args.cpu_frac
        if hasattr(args, "mem_frac") is False:
            mem_frac = 0.8
        else:
            mem_frac = args.mem_frac
        ncpu, mem = get_slurm_node_resources(
            partition=args.partition, cpu_frac=cpu_frac, mem_frac=mem_frac
        )
        script_args = [
            "#!/bin/bash",
            f"#SBATCH --job-name=paircars_{jobid}",
            f"#SBATCH --time={walltime}",
            f"#SBATCH --output={args.workdir}/paircars_{jobid}.log",
            f"#SBATCH --error={args.workdir}/paircars_{jobid}.log",
            f"#SBATCH --partition={args.partition}",
            "#SBATCH --nodes=1",
            "#SBATCH --ntasks=1",
            f"#SBATCH --cpus-per-task={min(8,ncpu)}",
            f"#SBATCH --mem={min(16,mem)}G",
        ]
        if hasattr(args, "account") and args.account is not None:
            script_args.append(f"#SBATCH --account={args.account}\n")
        script_args.append("init-paircars-prefect start\n")
        if len(prefect_env_list) > 0:
            for i in prefect_env_list:
                script_args.append(f"export {i}")
        script_args.append("export PYTHONUNBUFFERED=1\n")
        script_args.append(cli_cmd)
        script_path = os.path.join(args.workdir, f"paircars_slurm_{jobid}.sh")
        with open(script_path, "w") as f:
            for script_arg in script_args:
                f.write(f"{script_arg}\n")
        print("######################################################")
        print(f"P-AIRCARS Job ID: {jobid}")
        print(f"Batch script: {script_path} is ready for submission.")
        print("######################################################")
        result = subprocess.run(["sbatch", script_path], stderr=subprocess.DEVNULL)
        exit_code = result.returncode
        return 0 if exit_code == 0 else 1
    except Exception as e:
        traceback.print_exc()
        return 1
