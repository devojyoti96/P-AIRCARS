import argparse
import sys
import os
import numpy as np
from paircars.utils.basic_utils import get_cachedir
from paircars.utils.logger_utils import SmartDefaultsHelpFormatter
from paircars.utils.prefect_setup_utils import (
    start_server,
    stop_prefect_server,
    prefect_server_status,
    prefect_config,
    save_prefect_env_to_file,
    show_prefect_config,
    get_free_port,
)


def cli():
    parser = argparse.ArgumentParser(
        description="Manage a local Prefect server. Only for single-node work station.",
        formatter_class=SmartDefaultsHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", help="Sub-command help")
    # Start
    start_parser = subparsers.add_parser("start", help="Start the Prefect server")
    start_parser.add_argument(
        "--show-config",
        action="store_true",
        help="Display Prefect config after startup",
    )
    # Stop
    subparsers.add_parser("stop", help="Stop the Prefect server")

    # Status
    subparsers.add_parser("status", help="Check if the Prefect server is running")

    # Env
    subparsers.add_parser(
        "save_env", help="Save the Prefect environment to a .env file"
    )

    # Config
    subparsers.add_parser("config", help="Print the current Prefect config")

    parser.add_argument("--jobid", default="local", help="Job ID")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()

    if args.command == "start":
        port = get_free_port()
        msg, config_file, profile_path, env_file, dashboard, pid_file = start_server(
            port, show_config=args.show_config, jobid=args.jobid
        )
        if msg == 0:
            print(f"Prefect server started at port: {port}")
            print(f"Profile file: {profile_path}")
            print(f"Environment file: {env_file}")
            print(f"Dashboard file: {dashboard}")
            print(f"Server process ID file: {pid_file}")
        else:
            print(f"Error in starting prefect server at port: {port}")
    elif args.command == "save_env":
        profile_path, env_file, dashboard = save_prefect_env_to_file(jobid=args.jobid)
        print(f"Profile file: {profile_path}")
        print(f"Environment file: {env_file}")
        print(f"Dashboard file: {dashboard}")
    elif args.command == "config":
        show_prefect_config(jobid=args.jobid)
    else:
        cachedir = f"{get_cachedir()}/prefect_{args.jobid}"
        os.makedirs(cachedir, exist_ok=True)
        config_file = f"{cachedir}/prefect.config.npy"
        if os.path.exists(config_file) is False:
            print(f"Configuration file for job ID: {args.jobid} does not exist.")
            return 1
        config = np.load(config_file, allow_pickle=True).all()
        if args.command == "stop":
            msg = stop_prefect_server(jobid=args.jobid)
            if msg == 0:
                print(f"Prefect server stopped at: {config['SERVER_DASHBOARD']}.")
            else:
                print(
                    f"Error in stopping prefect server at: {config['SERVER_DASHBOARD']}."
                )
        elif args.command == "status":
            if prefect_server_status(jobid=args.jobid):
                print(f"Prefect server is running at: {config['SERVER_DASHBOARD']}")
            else:
                print(
                    f"Prefect server is not running at: {config['SERVER_DASHBOARD']}."
                )
        else:
            parser.print_help()


if __name__ == "__main__":
    cli()
