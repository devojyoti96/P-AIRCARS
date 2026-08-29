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
import numpy as np
from dotenv import load_dotenv
from dask.distributed import Client
from dask_jobqueue import PBSCluster
from pyfiglet import Figlet
from collections import deque
from paircars.utils.basic_utils import get_cachedir
from paircars.utils.proc_manage_utils import (
    get_scheduler_name,
    detect_best_interface,
    get_jobid,
    get_total_nodes,
)

def is_pbs_job():
    """
    Check whether the current process is running as a PBS job.
    """
    return any(
        var in os.environ
        for var in [
            "PBS_JOBID",
            "PBS_JOBNAME",
            "PBS_NODEFILE",
        ]
    )
    
def get_available_nodes(queue=None):
    """
    Get available nodes of a PBS queue.

    Parameters
    ----------
    queue : str, optional
        PBS queue name.

    Returns
    -------
    list
        Available node names.
    """
    cmd = ["pbsnodes", "-a", "-S"]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True,
    )
    available = []
    for line in result.stdout.splitlines():
        fields = line.split()
        # Skip header/empty lines
        if not fields or fields[0].lower() in ["vnode", "node"]:
            continue
        name = fields[0]
        state = fields[1] if len(fields) > 2 else ""
        q = fields[5]
        if state in ["free", "job-busy"]:
            if queue is None or q==queue:
                available.append(name)
    return available
    
    
