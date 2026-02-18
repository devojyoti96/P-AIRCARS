import types
import resource
import psutil
import dask
import numpy as np
import warnings
import gc
import logging
import time
import glob
import os
import subprocess
import sys
import tempfile
import shutil
import yaml
import socket
import traceback
from pathlib import Path
from dask import delayed, compute, config
from dask.distributed import Client, LocalCluster
from dask_jobqueue import SLURMCluster
from datetime import datetime as dt, timedelta
from .basic_utils import *


#################################
# Process management
#################################
def get_nprocess_paircarspipe(jobid):
    """
    Get numbers of processes currently running

    Parameters
    ----------
    workdir : str
        Work directory name
    jobid : int
        Job ID

    Returns
    -------
    int
        Number of running processes
    """
    cachedir = get_cachedir()
    pid_file = f"{cachedir}/pids/pids_{jobid}.txt"
    pids = np.loadtxt(pid_file, unpack=True)
    n_process = 0
    for pid in pids:
        if psutil.pid_exists(int(pid)):
            n_process += 1
    return n_process


def get_jobid():
    """
    Get Job ID with millisecond-level uniqueness.

    Returns
    -------
    int
        Job ID in the format YYYYMMDDHHMMSSmmm (milliseconds)
    """
    cachedir = get_cachedir()
    jobid_file = os.path.join(cachedir, "jobids.txt")
    if os.path.exists(jobid_file):
        prev_jobids = np.loadtxt(jobid_file, unpack=True, dtype="int64")
        if prev_jobids.size == 0:
            prev_jobids = []
        elif prev_jobids.size == 1:
            prev_jobids = [str(prev_jobids)]
        else:
            prev_jobids = [str(jid) for jid in prev_jobids]
    else:
        prev_jobids = []

    if len(prev_jobids) > 0:
        FORMAT = "%Y%m%d%H%M%S%f"
        CUTOFF = dt.utcnow() - timedelta(days=15)
        filtered_prev_jobids = []
        for job_id in prev_jobids:
            job_time = dt.strptime(job_id.ljust(20, "0"), FORMAT)  # pad if truncated
            if job_time >= CUTOFF or job_id == 0:  # Job ID 0 is always kept
                filtered_prev_jobids.append(job_id)
        prev_jobids = filtered_prev_jobids

    now = dt.utcnow()
    cur_jobid = (
        now.strftime("%Y%m%d%H%M%S") + f"{int(now.microsecond/1000):03d}"
    )  # ms = first 3 digits of microseconds
    prev_jobids.append(cur_jobid)

    job_ids_int = np.array(prev_jobids, dtype=np.int64)
    np.savetxt(jobid_file, job_ids_int, fmt="%d")

    return int(cur_jobid)


def save_main_process_info(pid, jobid, msdir, workdir, outdir, cpu_frac, mem_frac):
    """
    Save main processes info

    Parameters
    ----------
    pid : int
        Main job process id
    jobid : int
        Job ID
    msdir : str
        Measurement set directory
    workdir : str
        Work directory
    outdir : str
        Output directory
    cpu_frac : float
        CPU fraction of the job
    mem_frac : float
        Memory fraction of the job

    Returns
    -------
    str
        Job info file name
    """
    cachedir = get_cachedir()
    prev_main_pids = glob.glob(f"{cachedir}/main_pids_*.txt")
    prev_jobids = [
        str(os.path.basename(i).rstrip(".txt").split("main_pids_")[-1])
        for i in prev_main_pids
    ]
    if len(prev_jobids) > 0:
        FORMAT = "%Y%m%d%H%M%S%f"
        CUTOFF = dt.utcnow() - timedelta(days=15)
        filtered_prev_jobids = []
        for i in range(len(prev_jobids)):
            job_id = prev_jobids[i]
            job_time = dt.strptime(job_id.ljust(20, "0"), FORMAT)  # pad if truncated
            if job_time < CUTOFF or job_id == 0:  # Job ID 0 is always kept
                filtered_prev_jobids.append(job_id)
            else:
                os.system(f"rm -rf {prev_main_pids[i]}")
                if os.path.exists(f"{cachedir}/pids/pids_{job_id}.txt"):
                    os.system(f"rm -rf {cachedir}/pids/pids_{job_id}.txt")
    main_job_file = f"{cachedir}/main_pids_{jobid}.txt"
    main_str = f"{jobid} {pid} {msdir} {workdir} {outdir} {cpu_frac} {mem_frac}"
    with open(main_job_file, "w") as f:
        f.write(main_str)
    return main_job_file


def save_pid(pid, pid_file):
    """
    Save PID

    Parameters
    ----------
    pid : int
        Process ID
    pid_file : str
        File to save
    """
    try:
        pid_file = os.path.abspath(pid_file)
        os.makedirs(os.path.dirname(pid_file), exist_ok=True)
        pids = []
        if os.path.exists(pid_file):
            with open(pid_file, "r") as f:
                pids = [int(line.strip()) for line in f if line.strip()]
        pids.append(int(pid))
        with open(pid_file, "w") as f:
            for p in pids:
                f.write(f"{p}\n")
    except:
        pass


def generate_activate_env(outfile="activate_env.sh"):
    """
    Generate a shell script that activates the current Python environment.

    This works for both Conda and virtualenv environments and is safe for use in
    non-interactive shells (e.g., Slurm batch jobs) by explicitly sourcing `conda.sh`.

    If conda is not found in $PATH, it will try loading either `anaconda` or `anaconda3` module.

    Parameters
    ----------
    outfile : str
        Path to the shell script to write (default: ./activate_env.sh).

    Returns
    -------
    str
        Output file name
    """
    outfile = Path(outfile).expanduser().resolve()
    putfile = os.path.abspath(outfile)
    lines = ["#!/bin/bash", ""]

    def module_exists(name):
        """Check if a module exists using 'module avail'."""
        try:
            subprocess.run(
                ["module", "avail", name],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=False,
            )
            return True
        except Exception:
            return False

    # Conda-based environment
    if "CONDA_DEFAULT_ENV" in os.environ:
        conda_env = os.environ["CONDA_DEFAULT_ENV"]
        lines.append("# === Activate Conda Environment Safely ===")
        lines.append("if ! command -v conda >/dev/null 2>&1; then")
        if module_exists("anaconda"):
            lines.append("    module load anaconda")
        elif module_exists("anaconda3"):
            lines.append("    module load anaconda3")
        else:
            lines.append("    echo 'No Conda module found (anaconda or anaconda3)'")
            lines.append("    exit 1")
        lines.append("fi")
        lines.append("source $(conda info --base)/etc/profile.d/conda.sh")
        lines.append(f"conda activate {conda_env}")
    # Virtualenv-based environment
    elif "VIRTUAL_ENV" in os.environ:
        venv_path = os.environ["VIRTUAL_ENV"]
        lines.append("# === Activate Virtualenv ===")
        lines.append(f"source {venv_path}/bin/activate")
    else:
        python_path = sys.executable
        lines.append(
            "# === No Conda/Virtualenv Detected — Using current Python directly ==="
        )
        lines.append(f"echo 'No Conda or virtualenv detected; using: {python_path}'")
        lines.append(f"export PATH={os.path.dirname(python_path)}:$PATH")
    # Write file
    with open(outfile, "w") as f:
        f.write("\n".join(lines) + "\n")
    os.chmod(outfile, 0o755)
    print(f"Created activation script at: {outfile}")
    return outfile


def get_total_worker(cluster):
    """
    Get total workers in the cluster

    Parameters
    ----------
    cluster : dask.cluster
        Dask cluster

    Returns
    -------
    int
        Number of workers
    """
    return len(cluster.workers)


def scale_worker_and_wait(dask_cluster, nworker, timeout=60, poll_interval=1):
    """
    Scale worker and wait until it is done

    Parameters
    ----------
    dask_cluster : dask.cluster
        Dask cluster
    nworker : int
        Number of worker
    timeout : float, optional
        Timeout, show a warning and move
    poll_interval : float, optional
        Check interval in seconds
    """
    print(f"Start scaling to {nworker} workers")
    dask_cluster.scale(nworker)
    timeout = 60
    c = 0
    while c < timeout:
        if get_total_worker(dask_cluster) == nworker:
            print(f"Successfully scaled to {nworker} workers")
            return 0
        else:
            time.sleep(poll_interval)
            c += poll_interval
    print(f"Dask cluster did not scale to {nworker} within {timeout} seconds.")
    return 1


def wait_for_dask_workers(client, min_worker=1, timeout=60):
    """
    Wait until the Dask cluster has a minimum number of total and/or new workers.

    Parameters
    ----------
    client : dask.distributed.Client
        Dask client
    min_worker : int, optional
        Minimum new connected workers (default: 1)
    timeout : float, optional
        Maximum time to wait in seconds (default: 60)

    Raises
    ------
    TimeoutError
        If the required number of workers do not connect in time.
    """
    client.wait_for_workers(n_workers=min_worker, timeout=timeout)


def get_local_dask_cluster(
    dask_dir,
    mem_frac=0.8,
    ncpu=-1,
    mem=-1,
    spill_frac=0.7,
    verbose=True,
):
    """
    Create a local Dask cluster

    Parameters
    ----------
    dask_dir : str
        Dask temporary directory
    mem_frac : float, optional
        Fraction of total memory to use
    spill_frac : float, optional
        Spill to disk at this fraction
    verbose : bool, optional
        Verbose (details of cluster)

    Returns
    -------
    client : dask.distributed.Client
        Dask client
    cluster : dask.distributed.LocalCluster
        Dask cluster
    str
        Dask directory
    """
    logging.getLogger("distributed").setLevel(logging.ERROR)
    print("Creating local cluster on the current node.")
    # Set up Dask working directories
    dask_dir = os.path.join(dask_dir.rstrip("/"), f"dask_{int(time.time())}")
    dask_dir_tmp = os.path.join(dask_dir, "tmp")
    os.makedirs(dask_dir_tmp, exist_ok=True)
    try:
        total_mem = psutil.virtual_memory().total / 1024**3  # In GB
        mem_frac = min(mem_frac, 0.8)
        usable_mem = total_mem * mem_frac
        # Raise file descriptor limit
        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        if soft < int(hard * 0.8):
            resource.setrlimit(resource.RLIMIT_NOFILE, (int(hard * 0.8), hard))
        dask.config.set(
            {
                "temporary-directory": dask_dir,
                "distributed.worker.memory.target": spill_frac,
                "distributed.worker.memory.spill": spill_frac + 0.1,
                "distributed.worker.memory.pause": spill_frac + 0.2,
                "distributed.worker.memory.terminate": spill_frac + 0.25,
            }
        )
        cluster = LocalCluster(
            n_workers=1,
            threads_per_worker=1,
            memory_limit=f"{usable_mem}GB",
            local_directory=dask_dir,
            dashboard_address=":0",
            processes=True,
            env={
                "TMPDIR": dask_dir_tmp,
                "TMP": dask_dir_tmp,
                "TEMP": dask_dir_tmp,
                "DASK_TEMPORARY_DIRECTORY": dask_dir_tmp,
                "MALLOC_TRIM_THRESHOLD_": "0",
                "PYTHONWARNINGS": "ignore::UserWarning:contextlib",
            },
        )
        client = Client(cluster, heartbeat_interval="5s")
        client.run_on_scheduler(gc.collect)
        if verbose:
            print("####################################################")
            print(f"Dask dashboard available at: {client.dashboard_link}")
            print("####################################################")
        return client, cluster, dask_dir
    except Exception as e:
        print ("Error occured in creating local cluster.")
        traceback.print_exc()
        os.system(f"rm -rf {dask_dir_tmp}")


def detect_best_interface():
    """
    Automatically detect best network interface for Dask.

    Returns
    -------
    str
        Best interface
    """
    interfaces = psutil.net_if_addrs().keys()
    for iface in interfaces:
        if iface.startswith("ib"):
            return iface
    for iface in interfaces:
        if iface.startswith(("eth", "en")):
            return iface
    for iface in interfaces:
        if (
            iface != "lo"
            and not iface.startswith("wl")
            and not iface.startswith("docker")
            and not iface.startswith("veth")
            and not iface.startswith("br")
        ):
            return iface
    return None


def get_scheduler_name():
    """
    Get job scheduler available

    Returns
    -------
    str
        Scheduler name (local, pbs, slurm)
    """
    if shutil.which("sbatch"):
        return "slurm"
    elif shutil.which("bsub"):
        return "lsf"
    elif shutil.which("qhost"):
        return "sge"
    elif shutil.which("qsub"):
        return "pbs"
    elif shutil.which("condor_submit"):
        return "htcondor"
    elif shutil.which("msub"):
        return "mab"
    elif shutil.which("oarsub"):
        return "oar"
    else:
        return "local"


def get_total_nodes(partition=None):
    """
    Get total nodes

    Parameters
    ----------
    partitiion : str, optional
        Partition or queue (depending on type of scheduler)

    Returns
    -------
    int
        Total node number
    """
    if partition is None:
        print("No partition is given. Providing nodes of entire cluster.")
    scheduler_name = get_scheduler_name()
    if scheduler_name == "slurm":
        if partition is None:
            cmd = f"sinfo -h -o '%D'"
        else:
            cmd = f"sinfo -p {partition} -h -o '%D'"
        output = subprocess.check_output(cmd, shell=True).decode().strip().split()
        return sum(int(x) for x in output)
    elif scheduler_name == "pbs":
        if partition is None:
            cmd = "pbsnodes -a | grep 'Mom =' | wc -l"
        else:
            cmd = f"pbsnodes -a | grep 'queue = {partition}' | wc -l"
        output = subprocess.check_output(cmd, shell=True).decode().strip()
        return int(output)
    elif scheduler_name == "lsf":
        cmd = "bhosts -noheader | wc -l"
        output = subprocess.check_output(cmd, shell=True).decode().strip()
        return int(output)
    elif scheduler_name == "sge":
        cmd = "qhost | grep lx | wc -l"
        output = subprocess.check_output(cmd, shell=True).decode().strip()
        return int(output)
    elif scheduler_name == "htcondor":
        cmd = "condor_status -noheader | wc -l"
        output = subprocess.check_output(cmd, shell=True).decode().strip()
        return int(output)
    elif scheduler_name == "oar":
        if partition:
            cmd = f"oarnodes -l | grep {partition} | wc -l"
        else:
            cmd = "oarnodes -s | grep Alive | wc -l"

        output = subprocess.check_output(cmd, shell=True).decode().strip()
        return int(output)
    elif scheduler_name == "moab":
        cmd = "mdiag -n"
        output = subprocess.check_output(cmd, shell=True).decode()
        for line in output.splitlines():
            if "Total Nodes" in line:
                return int(line.split(":")[1].strip())
        return
    elif scheduler_name == "local":
        return 1
    else:
        return None


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
    mem = round(total_mem * mem_frac,1)
    return ncpu, mem


def create_slurm_config(
    output_path,
    dask_dir,
    log_dir,
    cpu_frac=0.8,
    mem_frac=0.8,
    partition=None,
    account=None,
    project=None,
    walltime="24:00:00",
    job_name="paircars",
    python_path=None,
    exclusive=True,
):
    """
    Create a SLURMCluster Dask YAML configuration file.

    Parameters
    ----------
    output_path : str
        Path where YAML file will be written
    dask_dir : str
        Slurm/Dask local directory to spill
    log_dir : str
        SLURM log directory (Not task log directory)
    cpu_frac : float, optional
        CPU fraction to use
    mem_frac : float, optional
        Memory fraction to use
    partition : str, optional
        SLURM partition name
        Note: If your cluster requires this, you should provide. Otherwise, error will occur.
    account : str, optional
        SLURM account name
        Note: If your cluster requires this, you should provide. Otherwise, error will occur.
    project : str, optional
        SLURM project code (if required)
    walltime : str, optional
        Job walltime, maximum time the SLURM job can run (HH:MM:SS)
    job_name : str, optional
        SLURM job name
    log_dir : str, optional
        Directory for SLURM logs
    local_directory : str, optional
        Worker scratch directory
    interface : str, optional
        Network interface (e.g., 'ib0', 'eth0')
    python_path : str, optional
        Explicit python executable (default: current python)
    exclusive : bool
        Whether to request exclusive node allocation

    Returns
    -------
    str
        Path to generated YAML file
    """
    cpu_frac = min(0.8, cpu_frac)
    mem_frac = min(0.8, mem_frac)
    ncpu, mem = get_slurm_node_resources(
        partition=partition, cpu_frac=cpu_frac, mem_frac=mem_frac
    )
    if python_path is None:
        python_path = sys.executable
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)
    job_extra = [
        f"--nodes=1",
        f"--ntasks=1",
        f"--cpus-per-task={ncpu}",
        f"--mem={mem}G",
    ]
    if exclusive:
        job_extra.append("--exclusive")
    if log_dir:
        job_extra.extend(
            [
                f"--output={log_dir}/dask-%j.out",
                f"--error={log_dir}/dask-%j.err",
            ]
        )
    config = {
        "jobqueue": {
            "slurm": {
                "queue": partition,
                "cores": ncpu,
                "processes": 1,
                "memory": mem,
                "walltime": walltime,
                "name": job_name,
                "python": python_path,
                "job-extra": job_extra,
                "death-timeout": 60,
                "shared-temp-directory": dask_dir,
                "env-extra": [
                    "OMP_NUM_THREADS=1",
                    "MKL_NUM_THREADS=1",
                    "OPENBLAS_NUM_THREADS=1",
                    "NUMEXPR_NUM_THREADS=1",
                    "MALLOC_TRIM_THRESHOLD_=0",
                ],
            }
        }
    }
    interface = detect_best_interface()
    if account:
        config["jobqueue"]["slurm"]["account"] = account
    if project:
        config["jobqueue"]["slurm"]["project"] = project
    if interface:
        config["jobqueue"]["slurm"]["interface"] = interface
    if dask_dir:
        config["jobqueue"]["slurm"]["local-directory"] = dask_dir
    with open(output_path, "w") as f:
        yaml.dump(config, f, sort_keys=False)
    return output_path


def get_slurm_dask_cluster(
    dask_dir,
    jobid=None,
    cpu_frac=0.8,
    mem_frac=0.8,
    partition=None,
    account=None,
    project=None,
    walltime="24:00:00",
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
    partition : str, optional
        SLURM partition name
        Note: If your cluster requires this, you should provide. Otherwise, error will occur.
    account : str, optional
        SLURM account name
        Note: If your cluster requires this, you should provide. Otherwise, error will occur.
    project : str, optional
        SLURM project code (if required)
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
    scheduler_name = get_scheduler_name()
    if scheduler_name!="slurm":
        print ("SLURM is not avilable as job scheduler in your cluster.")
        return 
    logging.getLogger("distributed").setLevel(logging.ERROR)

    cpu_frac = min(0.8, cpu_frac)
    mem_frac = min(0.8, mem_frac)
    if jobid is None:
        jobid = get_jobid()

    output_path = f"{dask_dir}/slurm_config_{jobid}.yaml"
    log_dir = f"{dask_dir}/slurm_log_{jobid}"
    os.makedirs(log_dir, exist_ok=True)
    
    dask_dir = os.path.join(dask_dir.rstrip("/"), f"dask_{int(time.time())}")
    dask_dir_tmp = os.path.join(dask_dir, "tmp")
    os.makedirs(dask_dir_tmp, exist_ok=True)

    try:
        slurm_config_yaml = create_slurm_config(
            output_path,
            dask_dir,
            log_dir,
            cpu_frac=cpu_frac,
            mem_frac=mem_frac,
            partition=partition,
            account=account,
            project=project,
            walltime=walltime,
            job_name=f"paircars_{jobid}",
            exclusive=True,
        )

        dask.config.set(
            {
                "temporary-directory": dask_dir,
                "distributed.worker.memory.target": spill_frac,
                "distributed.worker.memory.spill": spill_frac + 0.1,
                "distributed.worker.memory.pause": spill_frac + 0.2,
                "distributed.worker.memory.terminate": spill_frac + 0.25,
            }
        )

        with open(slurm_config_yaml, "r") as f:
            cluster_config = yaml.safe_load(f)
        dask.config.set(cluster_config)

        cluster = SLURMCluster(
            local_directory=dask_dir_tmp,
            env_extra=[
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
            print("####################################################")

        return client, cluster, dask_dir
    except Exception as e:
        print ("Error occured in creating SLURM cluster.")
        traceback.print_exc()
        os.system(f"rm -rf {output_path} {log_dir} {dask_dir_tmp}")


# Exposing only functions
__all__ = [
    name
    for name, obj in globals().items()
    if isinstance(obj, types.FunctionType) and obj.__module__ == __name__
]
