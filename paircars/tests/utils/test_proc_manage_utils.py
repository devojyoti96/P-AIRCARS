import pytest
from datetime import datetime as dt
from unittest.mock import patch, MagicMock, mock_open
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
    )
    # Check return code
    assert result == expected_return
    # Ensure scale() was called correctly
    mock_cluster.scale.assert_called_once_with(target_workers)


