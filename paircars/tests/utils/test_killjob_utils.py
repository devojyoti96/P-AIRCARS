import pytest
from unittest.mock import patch, MagicMock, Mock
import psutil
from paircars.utils.killjob_utils import *


@pytest.mark.parametrize(
    "port, port_status, lsof_output, expected_kill_calls",
    [
        # Case 1: Port closed → kill PIDs
        (8000, False, "1234\n5678\n", [1234, 5678]),
        # Case 2: Port open → do nothing
        (8000, True, "1234\n5678\n", []),
        # Case 3: No processes found
        (8000, False, "", []),
    ],
)
@patch("paircars.utils.killjob_utils.os.kill")
@patch("paircars.utils.killjob_utils.subprocess.run")
@patch("paircars.utils.killjob_utils.check_port_status")
def test_kill_port(
    mock_check_port,
    mock_subprocess,
    mock_kill,
    port,
    port_status,
    lsof_output,
    expected_kill_calls,
):
    mock_check_port.return_value = port_status

    mock_result = MagicMock()
    mock_result.stdout = lsof_output
    mock_subprocess.return_value = mock_result

    kill_port(port)

    if expected_kill_calls:
        actual_calls = [call.args[0] for call in mock_kill.call_args_list]
        assert actual_calls == expected_kill_calls
    else:
        mock_kill.assert_not_called()


@pytest.mark.parametrize("has_process", [True, False])
@patch("paircars.utils.killjob_utils.psutil.wait_procs")
@patch("paircars.utils.killjob_utils.psutil.Process")
def test_terminate_process_and_children(mock_process_cls, mock_wait_procs, has_process):
    mock_parent = MagicMock()
    mock_children = [MagicMock(), MagicMock()]
    mock_process_cls.return_value = mock_parent
    mock_parent.children.return_value = mock_children
    mock_wait_procs.return_value = ([], mock_children)

    if not has_process:
        mock_process_cls.side_effect = psutil.NoSuchProcess(9999)

    terminate_process_and_children(9999)

    if has_process:
        assert mock_parent.terminate.call_count == 1
        assert all(child.terminate.call_count == 1 for child in mock_children)
        assert all(child.kill.call_count == 1 for child in mock_children)
    else:
        mock_process_cls.assert_called_once()


@pytest.mark.parametrize(
    "file_exists, load_error, dask_error, outer_error",
    [
        (False, False, False, False),  # job file missing
        (True, True, False, False),  # loadtxt fails
        (True, False, False, False),  # normal execution
        (True, False, True, False),  # dask shutdown fails
        (True, False, False, True),  # outer exception
    ],
)
def test_kill_localscheduler(
    mocker,
    file_exists,
    load_error,
    dask_error,
    outer_error,
):
    from paircars.utils import killjob_utils

    jobid = 123

    mocker.patch("paircars.utils.killjob_utils.print")
    mocker.patch("paircars.utils.killjob_utils.traceback.print_exc")

    mocker.patch(
        "paircars.utils.killjob_utils.get_cachedir",
        side_effect=Exception("outer") if outer_error else lambda: "/mock/cache",
    )

    mocker.patch(
        "paircars.utils.killjob_utils.os.path.exists",
        return_value=file_exists,
    )

    if file_exists and not load_error and not outer_error:
        mocker.patch(
            "paircars.utils.killjob_utils.np.loadtxt",
            return_value=[
                "0",
                "9999",
                "tcp://scheduler:8786",
                "/mock/ms",
                "/mock/work",
                "/mock/out",
            ],
        )
    elif file_exists and load_error:
        mocker.patch(
            "paircars.utils.killjob_utils.np.loadtxt",
            side_effect=Exception("read error"),
        )

    mocker.patch(
        "paircars.utils.killjob_utils.terminate_process_and_children",
        Mock(),
    )

    if not outer_error:
        mock_client = MagicMock()
        if dask_error:
            mock_client.shutdown.side_effect = Exception("dask closed")
        mocker.patch(
            "paircars.utils.killjob_utils.Client",
            return_value=mock_client,
        )

    mocker.patch(
        "paircars.utils.killjob_utils.drop_cache",
        Mock(),
    )

    killjob_utils.kill_localscheduler(jobid)

    if not file_exists or outer_error:
        killjob_utils.terminate_process_and_children.assert_not_called()
    elif load_error:
        killjob_utils.terminate_process_and_children.assert_not_called()
    else:
        killjob_utils.terminate_process_and_children.assert_called_once_with(9999)


@pytest.mark.parametrize(
    "file_exists, load_error, dask_error, outer_error",
    [
        (False, False, False, False),
        (True, True, False, False),
        (True, False, False, False),
        (True, False, True, False),
        (True, False, False, True),
    ],
)
def test_kill_slurmscheduler(
    mocker,
    file_exists,
    load_error,
    dask_error,
    outer_error,
):
    from paircars.utils import killjob_utils

    jobid = 456

    mocker.patch("paircars.utils.killjob_utils.print")
    mocker.patch("paircars.utils.killjob_utils.traceback.print_exc")

    mocker.patch(
        "paircars.utils.killjob_utils.get_cachedir",
        side_effect=Exception("outer") if outer_error else lambda: "/mock/cache",
    )

    mocker.patch(
        "paircars.utils.killjob_utils.os.path.exists",
        return_value=file_exists,
    )

    if file_exists and not load_error and not outer_error:
        mocker.patch(
            "paircars.utils.killjob_utils.np.loadtxt",
            return_value=[
                "0",
                "7777",
                "tcp://scheduler:8786",
                "/mock/ms",
                "/mock/work",
                "/mock/out",
            ],
        )
    elif file_exists and load_error:
        mocker.patch(
            "paircars.utils.killjob_utils.np.loadtxt",
            side_effect=Exception("read error"),
        )

    mock_run = mocker.patch(
        "paircars.utils.killjob_utils.subprocess.run",
        Mock(),
    )

    if not outer_error:
        mock_client = MagicMock()
        if dask_error:
            mock_client.shutdown.side_effect = Exception("closed")
        mocker.patch(
            "paircars.utils.killjob_utils.Client",
            return_value=mock_client,
        )

    mock_drop = mocker.patch(
        "paircars.utils.killjob_utils.drop_cache",
        Mock(),
    )

    killjob_utils.kill_slurmscheduler(jobid)

    if not file_exists or outer_error or load_error:
        mock_run.assert_not_called()
    else:
        mock_run.assert_called_once_with(["scancel", "7777"])
        assert mock_drop.call_count == 4
