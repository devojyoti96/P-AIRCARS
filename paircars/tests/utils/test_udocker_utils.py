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
        (None, False, False, False),  # datadir None
        ("/mock/data", False, False, False),  # datadir missing
        ("/mock/data", True, False, False),  # tarball missing
        ("/mock/data", True, True, True),  # success case
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
    "returncode, side_effect, expected",
    [
        (0, None, True),                 # container exists
        (1, None, False),                # container does not exist
        (None, Exception("error"), False),  # subprocess raises exception
    ],
)
def test_check_udocker_container(mocker, returncode, side_effect, expected):
    mocker.patch("paircars.utils.udocker_utils.set_udocker_env")
    mock_run = mocker.patch("paircars.utils.udocker_utils.subprocess.run")
    if side_effect:
        mock_run.side_effect = side_effect
    else:
        mock_result = MagicMock()
        mock_result.returncode = returncode
        mock_run.return_value = mock_result
    result = check_udocker_container("test_container")
    assert result is expected


@pytest.mark.parametrize(
    "image_exists, update, verbose, pull_rc, expected",
    [
        (1, False, False, 0, "test_container"),
        (0, False, False, 0, "test_container"),
        (0, True, False, 0, "test_container"),
        (1, False, False, 1, None),
    ],
)
def test_initialize_container(
    mocker, image_exists, update, verbose, pull_rc, expected
):
    mocker.patch("paircars.utils.udocker_utils.set_udocker_env")
    mocker.patch(
        "paircars.utils.udocker_utils.os.system",
        return_value=image_exists,
    )
    mock_run = mocker.patch(
        "paircars.utils.udocker_utils.subprocess.run"
    )
    mock_result = MagicMock()
    mock_result.returncode = pull_rc
    mock_run.return_value = mock_result
    mocker.patch("paircars.utils.udocker_utils.print")
    result = initialize_container(
        image_name="test_image",
        name="test_container",
        update=update,
        verbose=verbose,
    )
    assert result == expected
    
    
@pytest.mark.parametrize(
    "update, verbose, returned_value",
    [
        (False, False, "paircarswsclean"),
        (True, False, "paircarswsclean"),
        (False, True, "paircarswsclean"),
        (True, True, None),  # simulate failure
    ],
)
def test_initialize_wsclean_container(
    mocker, update, verbose, returned_value
):
    mock_init = mocker.patch(
        "paircars.utils.udocker_utils.initialize_container",
        return_value=returned_value,
    )
    mocker.patch("paircars.utils.udocker_utils.print")
    result = initialize_wsclean_container(
        name="paircarswsclean",
        update=update,
        verbose=verbose,
    )
    mock_init.assert_called_once_with(
        "devojyoti96/wsclean-solar:latest",
        "paircarswsclean",
        update=update,
        verbose=verbose,
    )
    assert result == returned_value
    
    
@pytest.mark.parametrize(
    "update, verbose, returned_value",
    [
        (False, False, "paircarsquartical"),
        (True, False, "paircarsquartical"),
        (False, True, "paircarsquartical"),
        (True, True, None),  # simulate failure
    ],
)
def test_initialize_quartical_container(
    mocker, update, verbose, returned_value
):
    mock_init = mocker.patch(
        "paircars.utils.udocker_utils.initialize_container",
        return_value=returned_value,
    )
    mocker.patch("paircars.utils.udocker_utils.print")
    result = initialize_quartical_container(
        name="paircarsquartical",
        update=update,
        verbose=verbose,
    )
    mock_init.assert_called_once_with(
        "devojyoti96/quartical:0.2.6",
        "paircarsquartical",
        update=update,
        verbose=verbose,
    )
    assert result == returned_value
    
    
@pytest.mark.parametrize(
    "update, verbose, returned_value",
    [
        (False, False, "paircarsshadems"),
        (True, False, "paircarsshadems"),
        (False, True, "paircarsshadems"),
        (True, True, None),  # simulate failure
    ],
)
def test_initialize_shadems_container(
    mocker, update, verbose, returned_value
):
    mock_init = mocker.patch(
        "paircars.utils.udocker_utils.initialize_container",
        return_value=returned_value,
    )
    mocker.patch("paircars.utils.udocker_utils.print")
    result = initialize_shadems_container(
        name="paircarsshadems",
        update=update,
        verbose=verbose,
    )
    mock_init.assert_called_once_with(
        "devojyoti96/shadems:v0.5.4",
        "paircarsshadems",
        update=update,
        verbose=verbose,
    )
    assert result == returned_value
    
    
@pytest.mark.parametrize(
    "container_present, init_return, run_rc, raise_exc, expected",
    [
        (True, None, 0, False, 0),
        (True, None, 1, False, 1),
        (False, "paircarswsclean", 0, False, 0),
        (False, None, 0, False, 1),
        (True, None, 0, True, 1),
    ],
)
def test_run_wsclean(
    mocker,
    container_present,
    init_return,
    run_rc,
    raise_exc,
    expected,
):
    mocker.patch("paircars.utils.udocker_utils.set_udocker_env")
    mocker.patch("paircars.utils.udocker_utils.print")
    mocker.patch("paircars.utils.udocker_utils.traceback.print_exc")
    mocker.patch(
        "paircars.utils.udocker_utils.check_udocker_container",
        return_value=container_present,
    )
    mocker.patch(
        "paircars.utils.udocker_utils.initialize_wsclean_container",
        return_value=init_return,
    )
    mocker.patch(
        "paircars.utils.udocker_utils.tempfile._get_candidate_names",
        return_value=iter(["abc123"]),
    )
    mock_run = mocker.patch(
        "paircars.utils.udocker_utils.subprocess.run"
    )
    if raise_exc:
        mock_run.side_effect = Exception("run failed")
    else:
        mock_result = MagicMock()
        mock_result.returncode = run_rc
        mock_run.return_value = mock_result
    cmd = "wsclean -size 512 512 test.ms"
    result = run_wsclean(
        wsclean_cmd=cmd,
        container_name="paircarswsclean",
        check_container=True,
        verbose=False,
    )
    assert result == expected
    

@pytest.mark.parametrize(
    "container_present, init_return, only_uvw, run_rc, raise_exc, expected",
    [
        (True, None, False, 0, False, 0),
        (True, None, True, 0, False, 0),
        (True, None, False, 1, False, 1),
        (False, "paircarswsclean", False, 0, False, 0),
        (False, None, False, 0, False, 1),
        (True, None, False, 0, True, 1),
    ],
)
def test_run_solar_sidereal_cor(
    mocker,
    container_present,
    init_return,
    only_uvw,
    run_rc,
    raise_exc,
    expected,
):
    mocker.patch("paircars.utils.udocker_utils.set_udocker_env")
    mocker.patch("paircars.utils.udocker_utils.print")
    mocker.patch("paircars.utils.udocker_utils.traceback.print_exc")
    mocker.patch(
        "paircars.utils.udocker_utils.check_udocker_container",
        return_value=container_present,
    )
    mocker.patch(
        "paircars.utils.udocker_utils.initialize_wsclean_container",
        return_value=init_return,
    )
    mocker.patch(
        "paircars.utils.udocker_utils.tempfile._get_candidate_names",
        return_value=iter(["abc123"]),
    )
    mock_run = mocker.patch(
        "paircars.utils.udocker_utils.subprocess.run"
    )
    if raise_exc:
        mock_run.side_effect = Exception("failure")
    else:
        mock_result = MagicMock()
        mock_result.returncode = run_rc
        mock_run.return_value = mock_result
    result = run_solar_sidereal_cor(
        msname="test.ms",
        only_uvw=only_uvw,
        container_name="paircarswsclean",
        check_container=True,
        verbose=False,
    )
    assert result == expected
    
    
@pytest.mark.parametrize(
    "container_present, init_return, only_uvw, run_rc, raise_exc, expected",
    [
        (True, None, False, 0, False, 0),
        (True, None, True, 0, False, 0),
        (True, None, False, 1, False, 1),
        (False, "paircarswsclean", False, 0, False, 0),
        (False, None, False, 0, False, 1),
        (True, None, False, 0, True, 1),
    ],
)
def test_run_chgcenter(
    mocker,
    container_present,
    init_return,
    only_uvw,
    run_rc,
    raise_exc,
    expected,
):
    mocker.patch("paircars.utils.udocker_utils.set_udocker_env")
    mocker.patch("paircars.utils.udocker_utils.print")
    mocker.patch("paircars.utils.udocker_utils.traceback.print_exc")

    mocker.patch(
        "paircars.utils.udocker_utils.check_udocker_container",
        return_value=container_present,
    )

    mocker.patch(
        "paircars.utils.udocker_utils.initialize_wsclean_container",
        return_value=init_return,
    )

    mocker.patch(
        "paircars.utils.udocker_utils.tempfile._get_candidate_names",
        return_value=iter(["abc123"]),
    )

    mock_run = mocker.patch(
        "paircars.utils.udocker_utils.subprocess.run"
    )

    if raise_exc:
        mock_run.side_effect = Exception("failure")
    else:
        mock_result = MagicMock()
        mock_result.returncode = run_rc
        mock_run.return_value = mock_result

    result = run_chgcenter(
        msname="test.ms",
        ra="00h00m00.0s",
        dec="00d00m00.0s",
        only_uvw=only_uvw,
        container_name="paircarswsclean",
        check_container=True,
        verbose=False,
    )

    assert result == expected
    
    
@pytest.mark.parametrize(
    "container_present, init_return, cmd, run_rc, raise_exc, expected",
    [
        (True, None, "shadems test.ms", 0, False, 0),
        (True, None, "shadems test.ms", 1, False, 1),
        (True, None, "shadems --help", 0, False, 0),
        (False, "paircarsshadems", "shadems test.ms", 0, False, 0),
        (False, None, "shadems test.ms", 0, False, 1),
        (True, None, "shadems test.ms", 0, True, 1),
    ],
)
def test_run_shadems(
    mocker,
    container_present,
    init_return,
    cmd,
    run_rc,
    raise_exc,
    expected,
):
    mocker.patch("paircars.utils.udocker_utils.set_udocker_env")
    mocker.patch("paircars.utils.udocker_utils.print")
    mocker.patch("paircars.utils.udocker_utils.traceback.print_exc")
    mocker.patch(
        "paircars.utils.udocker_utils.check_udocker_container",
        return_value=container_present,
    )
    mocker.patch(
        "paircars.utils.udocker_utils.initialize_shadems_container",
        return_value=init_return,
    )
    mocker.patch(
        "paircars.utils.udocker_utils.tempfile._get_candidate_names",
        return_value=iter(["abc123"]),
    )
    mocker.patch(
        "paircars.utils.udocker_utils.os.getcwd",
        return_value="/tmp",
    )

    mock_run = mocker.patch(
        "paircars.utils.udocker_utils.subprocess.run"
    )

    if raise_exc:
        mock_run.side_effect = Exception("failure")
    else:
        mock_result = MagicMock()
        mock_result.returncode = run_rc
        mock_run.return_value = mock_result

    result = run_shadems(
        cmd=cmd,
        container_name="paircarsshadems",
        check_container=True,
        verbose=False,
    )

    assert result == expected
    
    
@pytest.mark.parametrize(
    "container_present, init_return, cmd, run_rc, raise_exc, expected",
    [
        (True, None, "goquartical", 0, False, 1),
        (
            True,
            None,
            "quartical input_ms.path=test.ms output.gain_directory=cal "
            "output.log_directory=log load_from=/other/path/gain/table",
            0,
            False,
            0,
        ),
        (
            True,
            None,
            "quartical input_ms.path=test.ms",
            1,
            False,
            1,
        ),
        (False, "paircarsquartical", "goquartical", 0, False, 1),
        (False, None, "goquartical", 0, False, 1),
        (True, None, "goquartical", 0, True, 1),
        (True, None, "", 0, False, 1),
    ],
)
def test_run_quartical(
    mocker,
    container_present,
    init_return,
    cmd,
    run_rc,
    raise_exc,
    expected,
):
    mocker.patch("paircars.utils.udocker_utils.set_udocker_env")
    mocker.patch("paircars.utils.udocker_utils.print")
    mocker.patch("paircars.utils.udocker_utils.traceback.print_exc")
    mocker.patch(
        "paircars.utils.udocker_utils.check_udocker_container",
        return_value=container_present,
    )
    mocker.patch(
        "paircars.utils.udocker_utils.initialize_quartical_container",
        return_value=init_return,
    )
    mocker.patch(
        "paircars.utils.udocker_utils.tempfile._get_candidate_names",
        return_value=iter(["abc123"]),
    )
    mocker.patch(
        "paircars.utils.udocker_utils.os.getcwd",
        return_value="/tmp",
    )
    mocker.patch(
        "paircars.utils.udocker_utils.os.system",
        return_value=0,
    )

    mock_run = mocker.patch(
        "paircars.utils.udocker_utils.subprocess.run"
    )

    if raise_exc:
        mock_run.side_effect = Exception("failure")
    else:
        mock_result = MagicMock()
        mock_result.returncode = run_rc
        mock_run.return_value = mock_result

    result = run_quartical(
        cmd=cmd,
        container_name="paircarsquartical",
        check_container=True,
        verbose=False,
    )

    assert result == expected
    
