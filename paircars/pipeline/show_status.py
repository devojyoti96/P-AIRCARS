import psutil
import argparse
import traceback
import glob
import sys
import os
import subprocess
from paircars.utils.basic_utils import get_cachedir
from paircars.utils.resource_utils import drop_cache
from paircars.utils.logger_utils import SmartDefaultsHelpFormatter
from paircars.utils.proc_manage_utils import get_scheduler_name

def is_slurm_job_running(job_id):
    result = subprocess.run(
        ["squeue", "-j", str(job_id), "-h"],
        capture_output=True,
        text=True
    )
    return result.stdout.strip() != ""


def show_local_job_status(clean_old_jobs=False):
    """
    Show P-AIRCARS local cluster jobs status

    Parameters
    ----------
    clean_old_jobs : bool, optional
        Clean old informations for stopped jobs
        
    Returns
    -------
    int
        Number of jobs running
    """
    cachedir = get_cachedir()
    msg=0
    try:
        main_pid_files = glob.glob(f"{cachedir}/main_pids_*.txt")
        if len(main_pid_files) == 0:
            print("No P-AIRCARS jobs is running.")
        else:
            print("####################")
            print("P-AIRCARS Job status")
            print("####################")
            for pid_file in main_pid_files:
                with open(pid_file, "r") as f:
                    line = f.read().split(" ")
                jobid = line[0]
                pid = line[1]
                workdir = line[4]
                outdir = line[5]
                if psutil.pid_exists(int(pid)):
                    running = "Running/Waiting"
                    msg+=1
                else:
                    running = "Done/Stopped"
                print(
                    f"Job ID: {jobid}, Work direcory: {workdir}, Output directory: {outdir}, Status: {running}"
                )
                print(
                    "#########################################################################################"
                )
                if clean_old_jobs and running == "Done/Stopped":
                    os.system(f"rm -rf {pid_file}")
    except Exception as e:
        traceback.print_exc()
    finally:
        return msg


def show_slurm_job_status(clean_old_jobs=False):
    """
    Show P-AIRCARS slurm cluster jobs status

    Parameters
    ----------
    clean_old_jobs : bool, optional
        Clean old informations for stopped jobs
        
    Returns
    -------
    int
        Number of jobs running
    """
    cachedir = get_cachedir()
    msg=0
    try:
        main_pid_files = glob.glob(f"{cachedir}/main_pids_*.txt")
        if len(main_pid_files) == 0:
            print("No P-AIRCARS jobs is running.")
        else:
            print("####################")
            print("P-AIRCARS Job status")
            print("####################")
            for pid_file in main_pid_files:
                with open(pid_file, "r") as f:
                    line = f.read().split(" ")
                jobid = line[0]
                pid = line[1]
                workdir = line[4]
                outdir = line[5]
                if is_slurm_job_running(int(pid)):
                    running = "Running/Waiting"
                    msg+=1
                else:
                    running = "Done/Stopped"
                print(
                    f"Job ID: {jobid}, Work direcory: {workdir}, Output directory: {outdir}, Status: {running}"
                )
                print(
                    "#########################################################################################"
                )
                if clean_old_jobs and running == "Done/Stopped":
                    os.system(f"rm -rf {pid_file}")
    except Exception as e:
        traceback.print_exc()
    finally:
        return msg
        
        
def cli():
    parser = argparse.ArgumentParser(
        description="Show P-AIRCARS jobs status.",
        formatter_class=SmartDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--show",
        action="store_true",
        dest="show",
        help="Show job status",
    )
    parser.add_argument(
        "--clean_old_jobs",
        action="store_true",
        help="Clean old jobs",
    )
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    scheduler_name = get_scheduler_name()
    try:
        args = parser.parse_args()
        if args.show:
            if scheduler_name=="local":
                show_local_job_status(clean_old_jobs=args.clean_old_jobs)
            elif scheduler_name=="slurm":
                show_slurm_job_status(clean_old_jobs=args.clean_old_jobs)    
            else:
                print (f"P-AIRCARS is not ready for job scheduler: {scheduler_name}")
    except Exception as e:
        traceback.print_exc()


if __name__ == "__main__":
    cli()
