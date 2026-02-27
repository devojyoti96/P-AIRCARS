import pytest
from unittest.mock import patch, MagicMock, Mock, call
import psutil
import numpy as np

from paircars.utils.killjob_utils import (
    terminate_process_and_children,
)


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
        mock_run.assert_called_once_with(["scancel", 7777])
        assert mock_drop.call_count == 4
