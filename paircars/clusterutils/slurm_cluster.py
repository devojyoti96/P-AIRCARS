import types
import dask
import gc
import time
import os
import subprocess
import sys
import traceback
import logging
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
            f"--mem={mem}G",
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
            memory=f"{mem}G",
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
            print("####################################################")

        return client, cluster, dask_dir
    except Exception as e:
        print("Error occured in creating SLURM cluster.")
        traceback.print_exc()
        os.system(f"rm -rf {output_path} {log_dir} {dask_dir}")
        
# Exposing only functions
__all__ = [
    name
    for name, obj in globals().items()
    if isinstance(obj, types.FunctionType) and obj.__module__ == __name__
]

