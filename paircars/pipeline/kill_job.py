import psutil
import numpy as np
import argparse
import time
import sys
import os
import signal
import traceback
import subprocess
from distributed import Client
from paircars.utils import get_cachedir, drop_cache, get_scheduler_name


def terminate_process_and_children(pid, grace_period=3.0):
    """
    Try to gracefully terminate a process tree, then force kill if needed.
    """
    try:
        parent = psutil.Process(pid)
        children = parent.children(recursive=True)
        for child in children:
            try:
                child.terminate()
            except Exception:
                pass
        parent.terminate()
        gone, alive = psutil.wait_procs([parent] + children, timeout=grace_period)
        for p in alive:
            try:
                p.kill()
            except Exception:
                pass
    except (psutil.NoSuchProcess, ProcessLookupError):
        pass


def kill_localscheduler(jobid):
    """
    Kill local scheduler
    
    Parameters
    ----------
    jobid : int
        P-AIRCARS job ID
    """
    try:
        cachedir = get_cachedir()
        jobfile_name = f"{cachedir}/main_pids_{jobid}.txt"
        if os.path.exists(jobfile_name) is False:
            print(f"No P-AIRCARS job information available for job ID; {jobfile_name}")
            return 
        try:
            results = np.loadtxt(jobfile_name, dtype="str", unpack=True)
            main_pid = int(results[1])
            scheduler_address = str(results[2])
            msdir = str(results[3])
            workdir = str(results[4])
            outdir = str(results[5])
        except Exception as e:
            print(f"Could not read job file.")
            traceback.print_exc()
            return
                    
        print(f"Attempting to terminate main PID: {main_pid}")
        terminate_process_and_children(main_pid)
        os.system(f"rm -rf {workdir}/tmp_paircars_*")

        try:
            print("Closing dask cluster....")
            client = Client(address=scheduler_address,timeout=30)
            client.shutdown()
            client.close()
        except:
            print (f"Dask cluster at: {scheduler_address} is already closed.")
            
        print("Dropping caches...")
        drop_cache(msdir)
        drop_cache(workdir)
        drop_cache(outdir)
        drop_cache(cachedir)
        print("Cleanup complete.")     
        return 
    except Exception as e:
        print(f"Error in killing P-AIRCARS job: {jobid}")  
        traceback.print_exc()
        return  

def kill_slurmscheduler(jobid):
    """
    Kill local scheduler
    
    Parameters
    ----------
    jobid : int
        P-AIRCARS job ID
    """
    try:
        cachedir = get_cachedir()
        jobfile_name = f"{cachedir}/main_pids_{jobid}.txt"
        if os.path.exists(jobfile_name) is False:
            print(f"No P-AIRCARS job information available for job ID; {jobfile_name}")
            return 
        try:
            results = np.loadtxt(jobfile_name, dtype="str", unpack=True)
            main_jobid = int(results[1])
            scheduler_address = str(results[2])
            msdir = str(results[3])
            workdir = str(results[4])
            outdir = str(results[5])
        except Exception as e:
            print(f"Could not read job file.")
            traceback.print_exc()
            return

        print(f"Attempting to terminate main slurm jobid: {main_jobid}")
        subprocess.run(["scancel", main_jobid])
        os.system(f"rm -rf {workdir}/tmp_paircars_*")
        
        try:
            print("Closing dask cluster....")
            client = Client(address=scheduler_address,timeput=30)
            client.shutdown()
            client.close()
        except:
            print (f"Dask cluster at: {scheduler_address} is already closed.")

        print("Dropping caches...")
        drop_cache(msdir)
        drop_cache(workdir)
        drop_cache(outdir)
        drop_cache(cachedir)
        print("Cleanup complete.")     
        return 
    except Exception as e:
        print(f"Error in killing P-AIRCARS job: {jobid}")  
        traceback.print_exc()
        return 


def kill_paircarsjob():
    """
    Gracefully terminate all processes related to a P-AIRCARS job.
    """
    parser = argparse.ArgumentParser(description="Kill P-AIRCARS Job")
    parser.add_argument(
        "--jobid", type=str, required=True, help="P-AIRCARS Job ID to kill"
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    
    scheduler_name = get_scheduler_name()
    
    if scheduler_name=="local":
        kill_localscheduler(args.jobid)
    elif scheduler_name=="slurm":
        kill_slurmscheduler(args.jobid)    
    
   
if __name__ == "__main__":
    kill_paircarsjob()
