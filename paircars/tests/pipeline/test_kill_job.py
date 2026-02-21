import pytest
from unittest.mock import patch, MagicMock, Mock, call
import psutil
import numpy as np

from paircars.pipeline.kill_job import (
    terminate_process_and_children,
    kill_paircarsjob,
)


@pytest.mark.parametrize("has_process", [True, False])
@patch("paircars.pipeline.kill_job.psutil.wait_procs")
@patch("paircars.pipeline.kill_job.psutil.Process")
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
        (True, True, False, False),    # loadtxt fails
        (True, False, False, False),   # normal execution
        (True, False, True, False),    # dask shutdown fails
        (True, False, False, True),    # outer exception
    ],
)
def test_kill_localscheduler(
    mocker,
    file_exists,
    load_error,
    dask_error,
    outer_error,
):
    from paircars.pipeline import kill_job

    jobid = 123

    mocker.patch("paircars.pipeline.kill_job.print")
    mocker.patch("paircars.pipeline.kill_job.traceback.print_exc")

    mocker.patch(
        "paircars.pipeline.kill_job.get_cachedir",
        side_effect=Exception("outer") if outer_error else lambda: "/mock/cache",
    )

    mocker.patch(
        "paircars.pipeline.kill_job.os.path.exists",
        return_value=file_exists,
    )

    if file_exists and not load_error and not outer_error:
        mocker.patch(
            "paircars.pipeline.kill_job.np.loadtxt",
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
            "paircars.pipeline.kill_job.np.loadtxt",
            side_effect=Exception("read error"),
        )

    mocker.patch(
        "paircars.pipeline.kill_job.terminate_process_and_children",
        Mock(),
    )

    if not outer_error:
        mock_client = MagicMock()
        if dask_error:
            mock_client.shutdown.side_effect = Exception("dask closed")
        mocker.patch(
            "paircars.pipeline.kill_job.Client",
            return_value=mock_client,
        )

    mocker.patch(
        "paircars.pipeline.kill_job.drop_cache",
        Mock(),
    )

    kill_job.kill_localscheduler(jobid)

    if not file_exists or outer_error:
        kill_job.terminate_process_and_children.assert_not_called()
    elif load_error:
        kill_job.terminate_process_and_children.assert_not_called()
    else:
        kill_job.terminate_process_and_children.assert_called_once_with(9999)
        
        
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
    from paircars.pipeline import kill_job

    jobid = 456

    mocker.patch("paircars.pipeline.kill_job.print")
    mocker.patch("paircars.pipeline.kill_job.traceback.print_exc")

    mocker.patch(
        "paircars.pipeline.kill_job.get_cachedir",
        side_effect=Exception("outer") if outer_error else lambda: "/mock/cache",
    )

    mocker.patch(
        "paircars.pipeline.kill_job.os.path.exists",
        return_value=file_exists,
    )

    if file_exists and not load_error and not outer_error:
        mocker.patch(
            "paircars.pipeline.kill_job.np.loadtxt",
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
            "paircars.pipeline.kill_job.np.loadtxt",
            side_effect=Exception("read error"),
        )

    mock_run = mocker.patch(
        "paircars.pipeline.kill_job.subprocess.run",
        Mock(),
    )

    if not outer_error:
        mock_client = MagicMock()
        if dask_error:
            mock_client.shutdown.side_effect = Exception("closed")
        mocker.patch(
            "paircars.pipeline.kill_job.Client",
            return_value=mock_client,
        )

    mock_drop = mocker.patch(
        "paircars.pipeline.kill_job.drop_cache",
        Mock(),
    )

    kill_job.kill_slurmscheduler(jobid)

    if not file_exists or outer_error or load_error:
        mock_run.assert_not_called()
    else:
        mock_run.assert_called_once_with(["scancel", 7777])
        assert mock_drop.call_count == 4
        
        
@pytest.mark.parametrize(
    "argv, scheduler_name, expect_exit, expect_local, expect_slurm",
    [
        (["prog"], None, True, False, False),
        (["prog", "--jobid", "123"], "local", False, True, False),
        (["prog", "--jobid", "123"], "slurm", False, False, True),
    ],
)
def test_kill_paircarsjob(
    mocker,
    monkeypatch,
    argv,
    scheduler_name,
    expect_exit,
    expect_local,
    expect_slurm,
):
    from paircars.pipeline import kill_job

    monkeypatch.setattr(kill_job.sys, "argv", argv)

    mocker.patch("paircars.pipeline.kill_job.print")

    if expect_exit:
        mock_exit = mocker.patch("paircars.pipeline.kill_job.sys.exit", side_effect=SystemExit)
        with pytest.raises(SystemExit):
            kill_job.kill_paircarsjob()
        mock_exit.assert_called_once_with(1)
        return

    mocker.patch(
        "paircars.pipeline.kill_job.get_scheduler_name",
        return_value=scheduler_name,
    )

    mock_local = mocker.patch(
        "paircars.pipeline.kill_job.kill_localscheduler",
        Mock(),
    )

    mock_slurm = mocker.patch(
        "paircars.pipeline.kill_job.kill_slurmscheduler",
        Mock(),
    )

    kill_job.kill_paircarsjob()

    if expect_local:
        mock_local.assert_called_once_with("123")
        mock_slurm.assert_not_called()
    elif expect_slurm:
        mock_slurm.assert_called_once_with("123")
        mock_local.assert_not_called()
        
