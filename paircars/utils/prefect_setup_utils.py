import os
import subprocess
import time
import socket
import signal
import argparse
import toml
import traceback
import numpy as np
from pathlib import Path
from dotenv import load_dotenv
from .basic_utils import get_cachedir


def get_free_port(start_port=4200, end_port=4300):
    """
    Get free port

    Parameters
    ----------
    start_port : int, optional
        Start port range
    end_port : int, optional
        End port range

    Returns
    -------
    int
        Free port
    """
    for port in range(4200, 4301):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind(("0.0.0.0", port))
                print(f"Free port: {port}")
                return port
            except OSError:
                continue


# === CONFIG ===
def prefect_config(port, jobid="local"):
    """
    Configure prefect

    Parameters
    ----------
    port : int
        Free port
    jobid : int, optional
        P-AIRCARS Job ID

    Returns
    -------
    str
        Configuration file name
    dict
        Configuration dictionary
    """
    cachedir = f"{get_cachedir()}/prefect_{jobid}"
    os.makedirs(cachedir, exist_ok=True)
    config_file = f"{cachedir}/prefect.config.npy"
    PREFECT_HOME = f"{cachedir}/prefect_home"
    os.makedirs(PREFECT_HOME, exist_ok=True)
    DB_URL = f"sqlite+aiosqlite:///{PREFECT_HOME}/prefect.db"
    LOG_FILE = os.path.join(PREFECT_HOME, "server.log")
    profile_path = os.path.join(PREFECT_HOME, "profiles.toml")
    memo_path = os.path.join(PREFECT_HOME, "memo_store.toml")
    storage = os.path.join(PREFECT_HOME, "storage")
    os.makedirs(storage, exist_ok=True)
    ENV_FILE = os.path.join(cachedir, "paircars_prefect.env")
    SERVER_HOST = "127.0.0.1"
    SERVER_PORT = f"{port}"
    #hostname = socket.gethostname()
    SERVER_URL = f"http://{SERVER_HOST}:{SERVER_PORT}/api"
    SERVER_DASHBOARD = f"http://{SERVER_HOST}:{SERVER_PORT}/dashboard"
    profile_name = f"paircarspipe_{jobid}"
    pid_file = os.path.join(PREFECT_HOME, "server.pid")
    logging_path = os.path.join(PREFECT_HOME, "logging.yml")
    config = {
        "CACHEDIR": cachedir,
        "PREFECT_HOME": PREFECT_HOME,
        "DB_URL": DB_URL,
        "LOG_FILE": LOG_FILE,
        "PROFILE_PATH": profile_path,
        "MEMO_PATH": memo_path,
        "STORAGE": storage,
        "ENV_FILE": ENV_FILE,
        "SERVER_HOST": SERVER_HOST,
        "SERVER_PORT": SERVER_PORT,
        "SERVER_URL": SERVER_URL,
        "SERVER_DASHBOARD": SERVER_DASHBOARD,
        "PROFILE_NAME": profile_name,
        "PID_FILE": pid_file,
        "LOGGING_PATH": logging_path,
    }
    np.save(config_file, config)
    return config_file, config


####################################
# Start and save
####################################
def write_prefect_profile(jobid="local"):
    """
    Save prefect profile

    Parameters
    ----------
    jobid : int, optional
        P-AIRCARS Job ID

    Returns
    -------
    str
        Profile file
    """
    cachedir = f"{get_cachedir()}/prefect_{jobid}"
    os.makedirs(cachedir, exist_ok=True)
    config_file = f"{cachedir}/prefect.config.npy"
    if os.path.exists(config_file) is False:
        print(f"Configuration file for job ID: {jobid} does not exist.")
        return False
    config = np.load(config_file, allow_pickle=True).all()
    # Load existing TOML config or start new
    profile_path = config["PROFILE_PATH"]
    if os.path.exists(profile_path):
        data = toml.load(profile_path)
    else:
        data = {}
    # Set active profile
    profile_name = config["PROFILE_NAME"]
    data["active"] = profile_name
    # Set config under [profiles.<profile_name>]
    if "profiles" not in data:
        data["profiles"] = {}
    data["profiles"][profile_name] = {
        "PREFECT_API_URL": config["SERVER_URL"],
        "PREFECT_HOME": config["PREFECT_HOME"],
        "PREFECT_API_DATABASE_CONNECTION_URL": config["DB_URL"],
    }
    with open(profile_path, "w") as f:
        toml.dump(data, f)
    print(f"Prefect profile '{profile_name}' written to {profile_path}")
    return profile_path


def save_prefect_env_to_file(jobid="local"):
    """
    Save current Prefect server env config to a .env file for reuse.

    Parameters
    ----------
    jobid : int, optional
        P-AIRCARS Job ID

    Returns
    -------
    str
        Profile file
    str
        Environment file
    str
        Dashboard file
    """
    cachedir = f"{get_cachedir()}/prefect_{jobid}"
    os.makedirs(cachedir, exist_ok=True)
    config_file = f"{cachedir}/prefect.config.npy"
    if os.path.exists(config_file) is False:
        print(f"Configuration file for job ID: {jobid} does not exist.")
        return False
    config = np.load(config_file, allow_pickle=True).all()
    cachedir = config["CACHEDIR"]
    env_file = config["ENV_FILE"]
    dashboard = f"{cachedir}/prefect.dashboard"
    with open(env_file, "w") as f:
        f.write(f"PREFECT_HOME={config['PREFECT_HOME']}\n")
        f.write("PREFECT_API_MODE=server\n")
        f.write(f"PREFECT_API_DATABASE_CONNECTION_URL={config['DB_URL']}\n")
        f.write("PREFECT_SERVER_ALLOW_EPHEMERAL_MODE=false\n")
        f.write(f"PREFECT_API_URL={config['SERVER_URL']}\n")
        f.write(f"PREFECT_PROFILE={config['PROFILE_NAME']}\n")
        f.write(f"PREFECT_PROFILES_PATH={config['PROFILE_PATH']}\n")
        f.write(f"PREFECT_LOCAL_STORAGE_PATH={config['STORAGE']}\n")
        f.write(f"PREFECT_LOGGING_SETTINGS_PATH={config['LOGGING_PATH']}\n")
        f.write(f"PREFECT_MEMO_STORE_PATH={config['MEMO_PATH']}\n")
    print(f"Saved Prefect server environment to {env_file}")
    if os.path.exists(dashboard) is not True:
        with open(dashboard, "w") as f:
            f.write(f"{config['SERVER_DASHBOARD']}")
    profile_path = write_prefect_profile(jobid=jobid)
    return profile_path, env_file, dashboard


def start_server(port, show_config=False, jobid="local"):
    """
    Start prefect server if it is not running

    Parameters
    ----------
    port : int
        Free port number
    show_config : bool, optional
        Show configuration of prefect server
    jobid : int, optional
        P-AIRCARS job ID

    Returns
    -------
    0, config_file, profile_path, env_file, dashboard, pid_file
    int
        Success message
    str
        Configuration file
    str
        Profile file
    str
        Environment file
    str
        Dashboard file
    str
        Server process ID file
    """
    config_file, config = prefect_config(port, jobid=jobid)
    cachedir = config["CACHEDIR"]
    pid_file = config["PID_FILE"]
    print("Starting Prefect server...")
    if prefect_server_status(jobid=jobid):
        stop_prefect_server(jobid=jobid)
    env = get_prefect_env(jobid=jobid)
    os.makedirs(config["PREFECT_HOME"], exist_ok=True)
    profile_path, env_file, dashboard = save_prefect_env_to_file(jobid=jobid)
    with open(config["LOG_FILE"], "w") as f:
        server_proc = subprocess.Popen(
            [
                "prefect",
                "server",
                "start",
                "--host",
                config["SERVER_HOST"],
                "--port",
                config["SERVER_PORT"],
            ],
            stdout=f,
            stderr=subprocess.STDOUT,
            env=env,
        )
    server_started = False
    for _ in range(1800):  # wait up to 1800s for the server to respond
        if prefect_server_status(jobid=jobid):
            if show_config:
                show_prefect_config(jobid=jobid)
            server_started = True
            break
        else:
            time.sleep(5)
    if server_started:
        with open(pid_file, "w") as pf:
            pf.write(str(server_proc.pid))
        print(f"Prefect server is now running at {config['SERVER_DASHBOARD']}")
        if os.path.exists(dashboard) is not True:
            with open(dashboard, "w") as f:
                f.write(f"{config['SERVER_DASHBOARD']}")
        return 0, config_file, profile_path, env_file, dashboard, pid_file
    else:
        print(
            f"Server did not respond within 30 minutes. Check logs at {config['LOG_FILE']} for more details"
        )
        return 0, config_file, profile_path, env_file, dashboard, pid_file


#########################################
# Stop prefect server
##########################################
def kill_port(port):
    """
    Kill a running port

    Parameters
    ----------
    port : int
        Port number
    """
    print(f"Closing previous prefect server at port : {port}.")
    result = subprocess.run(
        ["lsof", "-t", f"-i:{port}"],
        capture_output=True,
        text=True,
    )
    for pid in result.stdout.split():
        os.kill(int(pid), signal.SIGKILL)


def stop_prefect_server(jobid="local"):
    """
    Stop prefect server running in the current installation
    Note: it will only stop prefect server which is running from the current installation
    For this pipeline, a free port between 4200 to 4300 is chosen.

    Parameters
    ----------
    jobid : int, optional
        P-AIRCARS job ID

    Returns
    -------
    int
        Success message
    """
    cachedir = f"{get_cachedir()}/prefect_{jobid}"
    os.makedirs(cachedir, exist_ok=True)
    config_file = f"{cachedir}/prefect.config.npy"
    if os.path.exists(config_file) is False:
        print(f"Configuration file for job ID: {jobid} does not exist.")
        os.system(f"rm -rf {cachedir}")
        return 1
    config = np.load(config_file, allow_pickle=True).all()
    pid_file = config["PID_FILE"]
    cachedir = config["CACHEDIR"]
    try:
        if not os.path.exists(pid_file):
            try:
                kill_port(config["SERVER_PORT"])
                msg = 0
            except:
                msg = 1
        else:
            with open(pid_file, "r") as f:
                pid = int(f.read().strip())
            print(f"Stopping Prefect server with PID {pid} ...")
            os.kill(pid, signal.SIGTERM)
            print(f"Server stopped and {cachedir} removed.")
            msg = 0
    except ProcessLookupError:
        print(f"No such process with PID {pid}. Removing stale {cachedir} directory.")
        msg = 0
    except Exception as e:
        print(f"Error stopping server")
        traceback.print_exc()
        msg = 1
    finally:
        os.system(f"rm -rf {cachedir}")
        return msg


############################################
# Prefect server status
############################################
def prefect_server_status(jobid="local"):
    """
    Get prefect server status

    Parameters
    ----------
    jobid : int, optional
        P-AIRCARS job ID
    """
    cachedir = f"{get_cachedir()}/prefect_{jobid}"
    os.makedirs(cachedir, exist_ok=True)
    config_file = f"{cachedir}/prefect.config.npy"
    if os.path.exists(config_file) is False:
        print(f"Configuration file for job ID: {jobid} does not exist.")
        return False
    config = np.load(config_file, allow_pickle=True).all()
    try:
        with socket.create_connection(
            (config["SERVER_HOST"], config["SERVER_PORT"]), timeout=2
        ):
            return True
    except OSError:
        return False


def get_prefect_env(jobid="local"):
    """
    Get environment variables of prefect

    Parameters
    ----------
    jobid : int, optional
        P-AIRCARS job ID

    Returns
    -------
    dict
        Environment dictionary
    """
    cachedir = f"{get_cachedir()}/prefect_{jobid}"
    os.makedirs(cachedir, exist_ok=True)
    config_file = f"{cachedir}/prefect.config.npy"
    if os.path.exists(config_file) is False:
        print(f"Configuration file for job ID: {jobid} does not exist.")
        return
    config = np.load(config_file, allow_pickle=True).all()
    env = os.environ.copy()
    env["PREFECT_HOME"] = config["PREFECT_HOME"]
    env["PREFECT_API_MODE"] = "server"
    env["PREFECT_API_DATABASE_CONNECTION_URL"] = config["DB_URL"]
    env["PREFECT_SERVER_ALLOW_EPHEMERAL_MODE"] = "false"
    env["PREFECT_API_URL"] = config["SERVER_URL"]
    env["PREFECT_PROFILE"] = config["PROFILE_NAME"]
    env["PREFECT_PROFILES_PATH"] = config["PROFILE_PATH"]
    env["PREFECT_LOCAL_STORAGE_PATH"] = config["STORAGE"]
    env["PREFECT_LOGGING_SETTINGS_PATH"] = config["LOGGING_PATH"]
    env["PREFECT_MEMO_STORE_PATH"] = config["MEMO_PATH"]
    return env


def show_prefect_config(jobid="local"):
    """
    Print the effective Prefect config in this environment.

    Parameters
    ----------
    jobid : int, optional
        P-AIRCARS job ID
    """
    cachedir = f"{get_cachedir()}/prefect_{jobid}"
    os.makedirs(cachedir, exist_ok=True)
    config_file = f"{cachedir}/prefect.config.npy"
    if os.path.exists(config_file) is False:
        print(f"Configuration file for job ID: {jobid} does not exist.")
        return
    config = np.load(config_file, allow_pickle=True).all()
    load_dotenv(dotenv_path=config["ENV_FILE"], override=True)
    env = os.environ.copy()
    print("Prefect config in current environment ...")
    subprocess.run(["prefect", "config", "view"], env=env)
