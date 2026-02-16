import pytest
import psutil
import traceback
import tempfile
import os
from unittest.mock import patch, MagicMock
from paircars.utils.udocker_utils import *


@pytest.mark.parametrize(
    "datadir, exists_datadir, exists_tarball, should_set",
    [
        (None, False, False, False),          # datadir None
        ("/mock/data", False, False, False),  # datadir missing
        ("/mock/data", True, False, False),   # tarball missing
        ("/mock/data", True, True, True),     # success case
    ],
)
@patch("paircars.utils.udocker_utils.os.makedirs")
@patch("paircars.utils.udocker_utils.os.path.exists")
@patch("paircars.utils.udocker_utils.get_datadir")
def test_set_udocker_env(
    mock_get_datadir,
    mock_exists,
    mock_makedirs,
    datadir,
    exists_datadir,
    exists_tarball,
    should_set,
):
    with patch.dict(os.environ, {}, clear=True):

        mock_get_datadir.return_value = datadir

        def exists_side_effect(path):
            if datadir is None:
                return False
            if path == datadir:
                return exists_datadir
            if path == f"{datadir}/udocker-englib-1.2.11.tar.gz":
                return exists_tarball
            return False

        mock_exists.side_effect = exists_side_effect

        result = set_udocker_env()

        if should_set:
            mock_makedirs.assert_called_once_with(
                f"{datadir}/udocker",
                exist_ok=True,
            )
            assert result == datadir
            assert os.environ["UDOCKER_DIR"] == f"{datadir}/udocker"
            assert (
                os.environ["UDOCKER_TARBALL"]
                == f"{datadir}/udocker-englib-1.2.11.tar.gz"
            )
        else:
            mock_makedirs.assert_not_called()
            assert result is None
            assert "UDOCKER_DIR" not in os.environ
            assert "UDOCKER_TARBALL" not in os.environ


@patch("paircars.utils.udocker_utils.set_udocker_env")
def test_init_udocker(mock_env):
    init_udocker()


@pytest.mark.parametrize(
    "system_return, expected",
    [
        (0, True),  # Container present
        (1, False),  # Container absent
    ],
)
@patch("paircars.utils.udocker_utils.os.system")
@patch("paircars.utils.udocker_utils.set_udocker_env")
def test_check_udocker_container(mock_env, system_mock, system_return, expected):
    # First call: udocker inspect, Second call: cleanup
    system_mock.side_effect = [system_return, None]
    result = check_udocker_container("test_container")
    assert result is expected
    assert system_mock.call_count == 2


@pytest.mark.parametrize(
    "check_container, container_present, expected_return",
    [
        (True, False, 1),  # container check fails, fallback fails
        (False, True, 0),  # skip check, run successfully
    ],
)
@patch("paircars.utils.udocker_utils.traceback.print_exc")
@patch("paircars.utils.udocker_utils.psutil.Process")
@patch("paircars.utils.udocker_utils.os.system")
@patch("paircars.utils.udocker_utils.initialize_wsclean_container")
@patch("paircars.utils.udocker_utils.check_udocker_container")
@patch("paircars.utils.udocker_utils.tempfile.mkdtemp", return_value="/mock/temp")
@patch("paircars.utils.udocker_utils.os.getcwd", return_value="/mock")
@patch(
    "paircars.utils.udocker_utils.os.path.abspath", side_effect=lambda x: f"/abs/{x}"
)
@patch("paircars.utils.udocker_utils.os.path.dirname", side_effect=lambda x: "/abs")
@patch("paircars.utils.udocker_utils.set_udocker_env")
def test_run_wsclean_param_cases(
    mock_env,
    mock_dirname,
    mock_abspath,
    mock_getcwd,
    mock_mkdtemp,
    mock_check,
    mock_init,
    mock_system,
    mock_process,
    mock_traceback,
    check_container,
    container_present,
    expected_return,
):
    mock_check.return_value = container_present
    mock_init.return_value = None if not container_present else "paircarswsclean"
    mock_system.return_value = 0
    mock_process.return_value.memory_info.return_value.rss = 2.5 * 1024**3  # 2.5 GB
    result = run_wsclean(
        "wsclean -name mock test.ms",
        container_name="paircarswsclean",
        check_container=check_container,
        verbose=False,
    )
    assert result == expected_return


@pytest.mark.parametrize(
    "container_present, expected",
    [
        (False, 0),  # Container not found, init fails
        (True, 0),  # Normal run success
    ],
)
@patch("paircars.utils.udocker_utils.traceback.print_exc")
@patch("paircars.utils.udocker_utils.psutil.Process")
@patch("paircars.utils.udocker_utils.os.system")
@patch("paircars.utils.udocker_utils.initialize_wsclean_container")
@patch("paircars.utils.udocker_utils.check_udocker_container")
@patch("paircars.utils.udocker_utils.tempfile.mkdtemp", return_value="/mock/temp")
@patch("paircars.utils.udocker_utils.os.getcwd", return_value="/mock")
@patch(
    "paircars.utils.udocker_utils.os.path.abspath", side_effect=lambda x: f"/abs/{x}"
)
@patch("paircars.utils.udocker_utils.os.path.dirname", side_effect=lambda x: "/abs")
@patch("paircars.utils.udocker_utils.set_udocker_env")
def test_run_solar_sidereal_cor(
    mock_env,
    mock_dirname,
    mock_abspath,
    mock_getcwd,
    mock_mkdtemp,
    mock_check,
    mock_init,
    mock_system,
    mock_process,
    mock_traceback,
    container_present,
    expected,
):
    mock_check.return_value = container_present
    mock_init.return_value = None if not container_present else "paircarswsclean"
    mock_system.return_value = 0
    mock_process.return_value.memory_info.return_value.rss = 2.5 * 1024**3
    result = run_solar_sidereal_cor(
        msname="test.ms",
        only_uvw=False,
        container_name="paircarswsclean",
        verbose=False,
    )
    assert result == expected


@pytest.mark.parametrize(
    "container_present, expected",
    [
        (False, 0),  # container missing, init fails
        (True, 0),  # normal run, successful
    ],
)
@patch("paircars.utils.udocker_utils.traceback.print_exc")
@patch("paircars.utils.udocker_utils.psutil.Process")
@patch("paircars.utils.udocker_utils.os.system")
@patch("paircars.utils.udocker_utils.initialize_wsclean_container")
@patch("paircars.utils.udocker_utils.check_udocker_container")
@patch("paircars.utils.udocker_utils.tempfile.mkdtemp", return_value="/mock/temp")
@patch("paircars.utils.udocker_utils.os.getcwd", return_value="/mock")
@patch(
    "paircars.utils.udocker_utils.os.path.abspath", side_effect=lambda x: f"/abs/{x}"
)
@patch("paircars.utils.udocker_utils.os.path.dirname", side_effect=lambda x: "/abs")
@patch("paircars.utils.udocker_utils.set_udocker_env")
def test_run_chgcenter_param_cases(
    mock_env,
    mock_dirname,
    mock_abspath,
    mock_getcwd,
    mock_mkdtemp,
    mock_check,
    mock_init,
    mock_system,
    mock_process,
    mock_traceback,
    container_present,
    expected,
):
    mock_check.return_value = container_present
    mock_init.return_value = None if not container_present else "paircarswsclean"
    mock_system.return_value = 0
    mock_process.return_value.memory_info.return_value.rss = 2.5 * 1024**3
    result = run_chgcenter(
        msname="test.ms",
        ra="00:00:00.0",
        dec="-30:00:00.0",
        only_uvw=False,
        container_name="paircarswsclean",
        verbose=False,
    )
    assert result == expected
