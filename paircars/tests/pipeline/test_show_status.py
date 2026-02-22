import pytest
from unittest.mock import patch, MagicMock, mock_open
from paircars.pipeline.show_status import *


@pytest.mark.parametrize(
    "pid_alive, clean_old_jobs, expect_rm",
    [
        (True, False, False),
        (False, False, False),
        (False, True, True),
    ],
)
@patch("paircars.pipeline.show_status.drop_cache")
@patch("paircars.pipeline.show_status.os.path.exists", return_value=True)
@patch("paircars.pipeline.show_status.os.system")
@patch(
    "paircars.pipeline.show_status.open",
    new_callable=mock_open,
    read_data="1234 5678 dummy workdir outdir",
)
@patch("paircars.pipeline.show_status.glob.glob")
@patch("paircars.pipeline.show_status.get_cachedir", return_value="/mock/cache")
def test_show_job_status(
    mock_cachedir,
    mock_glob,
    mock_open_func,
    mock_system,
    mock_exists,
    mock_drop_cache,
    pid_alive,
    clean_old_jobs,
    expect_rm,
):
    # Mock job file present
    mock_glob.return_value = ["/mock/cache/main_pids_1234.txt"]
    show_job_status(clean_old_jobs=clean_old_jobs)


@pytest.mark.parametrize(
    "argv_args, expect_show_called, expect_exit_called",
    [
        (["show_paircars_status"], False, True),  # No args, should exit
        (["show_paircars_status", "--show"], True, False),  # Show only
        (
            ["show_paircars_status", "--show", "--clean_old_jobs"],
            True,
            False,
        ),  # Show + clean
    ],
)
@patch("paircars.pipeline.show_status.show_job_status")
@patch("paircars.pipeline.show_status.sys.exit")
@patch("paircars.pipeline.show_status.sys.argv", new_callable=list)
def test_cli_show_job_status(
    mock_argv,
    mock_exit,
    mock_show_status,
    argv_args,
    expect_show_called,
    expect_exit_called,
):
    # Patch argv directly
    sys.argv[:] = argv_args

    cli()

    if expect_show_called:
        mock_show_status.assert_called_once()
    else:
        mock_show_status.assert_not_called()

    if expect_exit_called:
        mock_exit.assert_called_once_with(1)
    else:
        mock_exit.assert_not_called()
