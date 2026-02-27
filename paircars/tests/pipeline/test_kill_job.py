import pytest
from unittest.mock import patch, MagicMock, Mock, call
import psutil
import numpy as np

from paircars.pipeline.kill_job import (
    kill_paircarsjob,
)


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
        mock_exit = mocker.patch(
            "paircars.pipeline.kill_job.sys.exit", side_effect=SystemExit
        )
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
