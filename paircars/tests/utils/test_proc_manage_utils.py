import pytest
import os
import sys
import tempfile
import numpy as np
import time
from pathlib import Path
from dask import delayed, compute
from dask.distributed import Client, LocalCluster
from datetime import datetime as dt
from unittest.mock import patch, MagicMock, mock_open, call
from itertools import chain, repeat
from paircars.utils.proc_manage_utils import *


@patch("paircars.utils.proc_manage_utils.np.savetxt")
@patch("paircars.utils.proc_manage_utils.np.loadtxt")
@patch("paircars.utils.proc_manage_utils.os.path.exists")
@patch("paircars.utils.proc_manage_utils.get_cachedir")
@patch("paircars.utils.proc_manage_utils.dt")
def test_get_jobid(mock_dt, mock_getdir, mock_exists, mock_loadtxt, mock_savetxt):
    fake_time = dt(2025, 7, 1, 15, 30, 45, 123456)
    mock_dt.utcnow.return_value = fake_time
    mock_getdir.return_value = "/mock/.paircars"
    mock_exists.return_value = False
    mock_loadtxt.return_value = []
    jobid = get_jobid()
    expected = int("20250701153045123")
    assert jobid == expected
    mock_savetxt.assert_called_once()


@patch("paircars.utils.proc_manage_utils.dt")
@patch("builtins.open", new_callable=mock_open)
@patch("paircars.utils.proc_manage_utils.os.system")
@patch("paircars.utils.proc_manage_utils.os.path.exists")
@patch("paircars.utils.proc_manage_utils.glob.glob")
@patch(
    "paircars.utils.proc_manage_utils.get_cachedir",
    return_value="/mock/.paircars",
)
def test_save_main_process_info(
    mock_get_cachedir,
    mock_glob,
    mock_exists,
    mock_system,
    mock_openfile,
    mock_dt,
):
    mock_glob.return_value = ["/mock/.paircars/main_pids_20250625000000000000.txt"]
    mock_exists.return_value = True
    fake_now = dt(2025, 7, 1, 0, 0, 0)
    mock_dt.utcnow.return_value = fake_now
    mock_dt.strptime.side_effect = lambda s, fmt: dt.strptime(s, fmt)
    result = save_main_process_info(
        1234,
        "20250701010101010101",
        "scheduler",
        "/mock/workdir",
        "/mock/workdir",
        "/mock/outdir",
        0.5,
        0.6,
    )
    expected_file = "/mock/.paircars/main_pids_20250701010101010101.txt"
    assert result == expected_file
    mock_openfile().write.assert_called_once_with(
        "20250701010101010101 1234 scheduler /mock/workdir /mock/workdir /mock/outdir 0.5 0.6"
    )
    mock_glob.return_value = ["/mock/.paircars/main_pids_20250625000000000000.txt"]


@pytest.mark.parametrize(
    "scheduler_info_return, expected",
    [
        ({"workers": {}}, 0),
        ({"workers": {"tcp://1": {}, "tcp://2": {}}}, 2),
        ({"workers": {"w1": {}, "w2": {}, "w3": {}}}, 3),
    ],
)
def test_get_total_worker(scheduler_info_return, expected):
    mock_client = MagicMock()
    mock_client.scheduler_info.return_value = scheduler_info_return
    result = get_total_worker(mock_client)
    assert result == expected
    mock_client.scheduler_info.assert_called_once()


@pytest.mark.parametrize(
    "worker_sequence, expected_return",
    [
        ([0, 1, 2, 3], 0),
    ],
)
@patch("paircars.utils.proc_manage_utils.time.sleep", return_value=None)
@patch("paircars.utils.proc_manage_utils.get_total_worker")
def test_scale_worker_and_wait(
    mock_get_total_worker,
    mock_sleep,
    worker_sequence,
    expected_return,
):
    target_workers = 3
    # Mock cluster
    mock_cluster = MagicMock()
    mock_client = MagicMock()
    # Configure get_total_worker to simulate scaling progression
    mock_get_total_worker.side_effect = worker_sequence
    result = scale_worker_and_wait(
        mock_cluster,
        mock_client,
        target_workers,
        timeout=10,
        poll_interval=1,
    )
    # Check return code
    assert result == expected_return
    # Ensure scale() was called correctly
    mock_cluster.scale.assert_called_once_with(target_workers)
    # Ensure polling happened
    assert mock_get_total_worker.called


@pytest.mark.parametrize(
    "total_cpu,total_mem,min_mem,max_worker,raise_cpu_error,expect_none",
    [
        # Normal case
        (16, 64, 4, -1, False, False),
        # Low memory → early return (n_worker_mem < 2)
        (4, 4, 4, -1, False, True),
        # max_worker restriction
        (16, 64, 4, 2, False, False),
        # Exception branch
        (None, None, 4, -1, True, True),
    ],
)
@patch("paircars.utils.proc_manage_utils.os.system")
@patch("paircars.utils.proc_manage_utils.traceback.print_exc")
@patch("paircars.utils.proc_manage_utils.Client")
@patch("paircars.utils.proc_manage_utils.LocalCluster")
@patch("paircars.utils.proc_manage_utils.psutil.virtual_memory")
@patch("paircars.utils.proc_manage_utils.psutil.cpu_count")
@patch("paircars.utils.proc_manage_utils.dask.config.set")
@patch("paircars.utils.proc_manage_utils.resource.setrlimit")
@patch("paircars.utils.proc_manage_utils.resource.getrlimit")
@patch("paircars.utils.proc_manage_utils.os.makedirs")
@patch("paircars.utils.proc_manage_utils.time.time", return_value=123456)
def test_get_local_dask_cluster(
    mock_time,
    mock_makedirs,
    mock_getrlimit,
    mock_setrlimit,
    mock_dask_config,
    mock_cpu_count,
    mock_virtual_memory,
    mock_localcluster,
    mock_client,
    mock_traceback,
    mock_system,
    total_cpu,
    total_mem,
    min_mem,
    max_worker,
    raise_cpu_error,
    expect_none,
):
    if raise_cpu_error:
        mock_cpu_count.side_effect = Exception("CPU fail")
    else:
        mock_cpu_count.return_value = total_cpu
        mock_virtual_memory.return_value = MagicMock(total=total_mem * 1024**3)
    mock_getrlimit.return_value = (1024, 4096)
    fake_cluster = MagicMock()
    fake_client = MagicMock()
    fake_client.dashboard_link = "http://localhost:8787"
    mock_localcluster.return_value = fake_cluster
    mock_client.return_value = fake_client
    result = get_local_dask_cluster(
        dask_dir="/tmp/test",
        min_mem=min_mem,
        max_worker=max_worker,
        verbose=False,
    )
    if expect_none:
        assert result is None
        if raise_cpu_error:
            mock_traceback.assert_called_once()
            mock_system.assert_called_once()
        else:
            mock_localcluster.assert_not_called()
        return
    client, cluster, dask_dir, n_worker = result
    assert client == fake_client
    assert cluster == fake_cluster
    assert "dask_123456" in dask_dir
    assert n_worker >= 2
    mock_localcluster.assert_called_once()
    mock_client.assert_called_once_with(fake_cluster, heartbeat_interval="5s")


@pytest.mark.parametrize(
    "scheduler_name,workdir_exists,log2term,config_exists,subprocess_fail,outer_fail,expected_return",
    [
        ("slurm", True, False, False, False, False, 1),
        ("local", False, False, False, False, False, 1),
        ("local", True, False, False, False, False, 0),
        ("local", True, True, True, False, False, 0),
        ("local", True, False, False, True, False, 1),
        ("local", True, False, False, False, True, 1),
    ],
)
@patch("paircars.utils.proc_manage_utils.traceback.print_exc")
@patch("paircars.utils.proc_manage_utils.subprocess.Popen")
@patch("paircars.utils.proc_manage_utils.load_dotenv")
@patch("paircars.utils.proc_manage_utils.np.load")
@patch("paircars.utils.proc_manage_utils.os.path.exists")
@patch("paircars.utils.proc_manage_utils.get_cachedir")
@patch("paircars.utils.proc_manage_utils.get_scheduler_name")
@patch("paircars.utils.proc_manage_utils.os.makedirs")
@patch("paircars.utils.proc_manage_utils.open", new_callable=mock_open)
@patch("paircars.utils.proc_manage_utils.Figlet")
def test_submit_local_master_flow(
    mock_figlet,
    mock_open_file,
    mock_makedirs,
    mock_get_scheduler,
    mock_get_cachedir,
    mock_exists,
    mock_npload,
    mock_load_dotenv,
    mock_popen,
    mock_traceback,
    scheduler_name,
    workdir_exists,
    log2term,
    config_exists,
    subprocess_fail,
    outer_fail,
    expected_return,
):
    mock_get_scheduler.return_value = scheduler_name
    args = MagicMock()
    args.workdir = "/tmp/testdir" if workdir_exists else None
    args.log2term = log2term
    jobid = 123

    # -------------------------
    # Config mocking
    # -------------------------
    mock_get_cachedir.return_value = "/tmp/cache"
    mock_exists.return_value = config_exists

    if config_exists:
        mock_npload.return_value.all.return_value = {
            "ENV_FILE": "/tmp/.env",
            "NODE_URL": "http://localhost:4200",
        }

    # -------------------------
    # Mock Figlet
    # -------------------------
    mock_figlet.return_value.renderText.return_value = "ASCII"

    # -------------------------
    # Mock subprocess
    # -------------------------
    if subprocess_fail:
        mock_popen.side_effect = Exception("Subprocess failed")
    else:
        fake_process = MagicMock()
        fake_process.stdout = ["Task run started\n", "Flow run done\n"]
        mock_popen.return_value = fake_process

    # -------------------------
    # Force outer failure
    # -------------------------
    if outer_fail:
        mock_open_file.side_effect = Exception("File write fail")

    # -------------------------
    # Run
    # -------------------------
    result = submit_local_master_flow(args, jobid)

    # -------------------------
    # Assertions
    # -------------------------
    assert result == expected_return

    if scheduler_name != "local":
        return

    if not workdir_exists:
        return

    if outer_fail:
        mock_traceback.assert_called_once()
        return

    if subprocess_fail:
        assert result == 1
    elif expected_return == 0:
        mock_popen.assert_called()
