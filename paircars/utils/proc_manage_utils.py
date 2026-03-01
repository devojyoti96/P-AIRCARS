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
import socket
import shlex
import traceback
from pathlib import Path
from dotenv import load_dotenv
from dask.distributed import Client, LocalCluster
from datetime import datetime as dt, timedelta
from .basic_utils import get_cachedir


#################################
# Process management
#################################
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


def save_main_process_info(
    pid, jobid, scheduler_address, msdir, workdir, outdir, cpu_frac, mem_frac
):
    """
    Save main processes info

    Parameters
    ----------
    pid : int
        Main job process id for local cluster or scheduler job id
    jobid : int
        Job ID
    scheduler_address : str
        Dask scheduler address
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
            if job_time > CUTOFF or job_id == 0:  # Job ID 0 is always kept
                filtered_prev_jobids.append(job_id)
            else:
                os.system(f"rm -rf {prev_main_pids[i]}")
    main_job_file = f"{cachedir}/main_pids_{jobid}.txt"
    main_str = f"{jobid} {pid} {scheduler_address} {msdir} {workdir} {outdir} {cpu_frac} {mem_frac}"
    with open(main_job_file, "w") as f:
        f.write(main_str)
    return main_job_file


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


def get_total_worker(client):
    """
    Get total workers in the cluster

    Parameters
    ----------
    client : dask.client
        Dask client for the cluster

    Returns
    -------
    int
        Number of workers
    """
    return len(client.scheduler_info()["workers"])


def scale_worker_and_wait(
    dask_cluster, dask_client, nworker, timeout=60, poll_interval=1
):
    """
    Scale worker and wait until it is done

    Parameters
    ----------
    dask_cluster : dask.cluster
        Dask cluster
    dask_client : dask.client
        Dask client for the same cluster
    nworker : int
        Number of worker
    timeout : float, optional
        Timeout, show a warning and move
    poll_interval : float, optional
        Check interval in seconds
    """
    print(f"Start scaling to {nworker} workers")
    dask_cluster.scale(nworker)
    time.sleep(1)
    timeout = 60
    c = 0
    while c < timeout:
        running_worker = get_total_worker(dask_client)
        if running_worker == nworker:
            print(f"Successfully scaled to {running_worker} workers")
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
    cpu_frac=0.8,
    mem_frac=0.8,
    min_mem=1,
    max_worker=-1,
    spill_frac=0.7,
    verbose=True,
):
    """
    Create a local Dask cluster

    Parameters
    ----------
    dask_dir : str
        Dask temporary directory
    cpu_frac : float, optional
        CPU fraction to use
    mem_frac : float, optional
        Fraction of total memory to use
    min_mem : float, optional
        Minimum required per job memory in GB
    max_worker : int, optional
        Maximum worker
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
    int
        Number of workers
    """
    cpu_frac = min(abs(cpu_frac), 0.8)
    mem_frac = min(abs(mem_frac), 0.8)
    logging.getLogger("distributed").setLevel(logging.ERROR)
    print("Creating local cluster on the current node.")
    # Set up Dask working directories
    dask_dir = os.path.join(dask_dir.rstrip("/"), f"dask_{int(time.time())}")
    dask_dir_tmp = os.path.join(dask_dir, "tmp")
    os.makedirs(dask_dir_tmp, exist_ok=True)
    try:
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
        usable_cpu = max(1, int(psutil.cpu_count() * cpu_frac))
        total_mem = psutil.virtual_memory().total / 1024**3  # In GB
        usable_mem = round(total_mem * mem_frac, 2)
        n_worker_mem = int(usable_mem / min_mem)
        if n_worker_mem < 2:
            print("Minimum available memory is not sufficient for at-least 2 workers.")
            return
        n_worker_cpu = usable_cpu
        n_worker = min(n_worker_cpu, n_worker_mem)
        if max_worker > 0:
            n_worker = min(n_worker, max_worker)
            n_worker = max(2, n_worker)

        mem_limit = round(usable_mem / n_worker, 2)
        ncpu = max(1, int(usable_cpu / n_worker))

        cluster = LocalCluster(
            n_workers=1,
            threads_per_worker=1,
            memory_limit=f"{mem_limit}GiB",
            local_directory=dask_dir,
            dashboard_address=":0",
            processes=True,
            env={
                "PYTHONUNBUFFERED": "1",
                "OMP_NUM_THREADS": f"{ncpu}",
                "MKL_NUM_THREADS": f"{ncpu}",
                "OPENBLAS_NUM_THREADS": f"{ncpu}",
                "NUMEXPR_NUM_THREADS": f"{ncpu}",
                "RAYON_NUM_THREADS": f"{ncpu}",
                "MALLOC_TRIM_THRESHOLD_": "0",
                "TMPDIR": f"{dask_dir_tmp}",
                "TMP": f"{dask_dir_tmp}",
                "TEMP": f"{dask_dir_tmp}",
                "DASK_TEMPORARY_DIRECTORY": f"{dask_dir_tmp}",
                "PYTHONWARNINGS": "ignore::UserWarning:contextlib",
            },
        )
        client = Client(cluster, heartbeat_interval="5s")
        client.run_on_scheduler(gc.collect)
        if verbose:
            print("####################################################")
            print(f"Dask dashboard available at: {client.dashboard_link}")
            print(f"Total usable cpu: {usable_cpu}")
            print(f"Total usable memory: {usable_mem} GB")
            print(f"CPU per worker: {ncpu}")
            print(f"Memory per worker: {mem_limit}GB")
            print(f"Maximum number of workers: {n_worker}")
            print("####################################################")
        return client, cluster, dask_dir, n_worker
    except Exception as e:
        print("Error occured in creating local cluster.")
        traceback.print_exc()
        os.system(f"rm -rf {dask_dir}")
        os.system(f"rm -rf {dask_dir_tmp}")
        return


def submit_local_master_flow(args, jobid):
    """
    Submit P-AIRCARS master flow to a local cluster

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
    if scheduler_name != "local":
        print(
            f"Job scheduler is not local. Available job scheduler is : {scheduler_name}"
        )
        return 1
    cli_cmd = "run-mwa-masterflow " + " ".join(shlex.quote(arg) for arg in sys.argv[1:])
    if hasattr(args, "workdir") and args.workdir is not None:
        os.makedirs(args.workdir, exist_ok=True)
    else:
        print("Please provide a work directory.")
        return 1

    cachedir = f"{get_cachedir()}/prefect_{scheduler_name}"
    config_file = f"{cachedir}/prefect.config.npy"
    prefect_env_list = []
    if os.path.exists(config_file):
        config = np.load(config_file, allow_pickle=True).all()
        load_dotenv(dotenv_path=config["ENV_FILE"], override=True)
        envlist = os.environ
        envlist["PREFECT_API_URL"] = config["NODE_URL"]
        for env in envlist:
            if "PREFECT" in env:
                prefect_env_list.append(f"export {env}={envlist.get(env)}")
    log_file = f"{args.workdir}/main_paircars_{jobid}.log"
    cli_cmd.append("--masterlog")
    cli_cmd.append(f"{log_file}")
    try:
        script_args = ["#!/bin/bash\n"]
        if len(prefect_env_list) > 0:
            for i in prefect_env_list:
                script_args.append(f"{i}")
        script_args.append("export PYTHONUNBUFFERED=1\n")
        script_args.append(cli_cmd)
        script_path = os.path.join(args.workdir, f"paircars_local_{jobid}.sh")
        with open(script_path, "w") as f:
            for script_arg in script_args:
                f.write(f"{script_arg}\n")        
        with open(log_file, "w") as log:
            result = subprocess.run(
                ["bash", script_path],
                stdout=log,
                stderr=log,
                text=True
            )
        exit_code = result.returncode
        return 0 if exit_code == 0 else 1
    except Exception as e:
        traceback.print_exc()
        return 1


##############################################
# Scheduler and hardware architecture related
##############################################
def detect_best_interface(scheduler_ip=None):
    """
    Automatically detect best IPv4 network interface for Dask.

    Parameters
    ----------
    scheduler_ip : str, optional
        If provided, ensures selected interface can route to this IP.

    Returns
    -------
    str or None
        Best interface name or None if not found.
    """
    stats = psutil.net_if_stats()
    addrs = psutil.net_if_addrs()
    candidates = []
    for iface, iface_addrs in addrs.items():
        if iface == "lo":
            continue
        if iface.startswith(("docker", "veth", "br", "wl")):
            continue
        # Interface must be UP
        if iface not in stats or not stats[iface].isup:
            continue
        # Must have IPv4
        ipv4 = None
        for addr in iface_addrs:
            if addr.family == socket.AF_INET:
                if not addr.address.startswith("127."):
                    ipv4 = addr.address
                    break
        if ipv4:
            candidates.append((iface, ipv4))
    if not candidates:
        return None
    # If scheduler_ip provided, check routing
    if scheduler_ip:
        for iface, ip in candidates:
            try:
                s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                s.bind((ip, 0))
                s.settimeout(1)
                s.connect((scheduler_ip, 80))
                s.close()
                return iface
            except Exception:
                continue
    # Prefer InfiniBand if available
    for iface, _ in candidates:
        if iface.startswith("ib"):
            return iface
    # Otherwise return first valid interface
    return candidates[0][0]


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
