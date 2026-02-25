import argparse
import sys
from paircars.utils.proc_manage_utils import get_scheduler_name
from paircars.utils.logger_utils import SmartDefaultsHelpFormatter
from paircars.utils.killjob_utils import (
    terminate_process_and_children,
    kill_localscheduler,
    kill_slurmscheduler,
)


def kill_paircarsjob():
    """
    Gracefully terminate all processes related to a P-AIRCARS job.
    """
    parser = argparse.ArgumentParser(
        description="Kill P-AIRCARS Job", formatter_class=SmartDefaultsHelpFormatter
    )
    parser.add_argument(
        "--jobid", type=str, required=True, help="P-AIRCARS Job ID to kill"
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    scheduler_name = get_scheduler_name()

    if scheduler_name == "local":
        kill_localscheduler(args.jobid)
    elif scheduler_name == "slurm":
        kill_slurmscheduler(args.jobid)


if __name__ == "__main__":
    kill_paircarsjob()
