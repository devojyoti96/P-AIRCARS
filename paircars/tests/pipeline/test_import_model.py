import pytest
import numpy as np
from unittest.mock import patch, MagicMock
from paircars.pipeline.import_model import *


@pytest.mark.parametrize("raise_error", [False, True])
def test_import_hyperdrive(tmp_path, monkeypatch, raise_error):
    msname = "test.ms"
    metafits = "test.metafits"
    monkeypatch.setattr(
        "paircars.pipeline.import_model.datadir",
        str(tmp_path),
    )
    monkeypatch.setattr(
        "glob.glob",
        lambda x: [str(tmp_path / "mwa_full_embedded_element_pattern_150.h5")]
    )
    monkeypatch.setattr("os.path.exists", lambda x: True)
    if raise_error:
        def fake_run(*args, **kwargs):
            raise RuntimeError("hyperdrive failed")
    else:
        def fake_run(*args, **kwargs):
            return None
    monkeypatch.setattr("subprocess.run", fake_run)
    fake_msmd = MagicMock()
    fake_msmd.nchan.return_value = 4
    fake_msmd.meanfreq.return_value = 150.0
    fake_msmd.chanres.return_value = [40.0]
    fake_msmd.ncorrforpol.return_value = [2]
    fake_msmd.nantennas.return_value = 4
    fake_msmd.timesforfield.return_value = [1, 2, 3]
    fake_msmd.exposuretime.return_value = {"value": 2.0}
    fake_msmd.nrows.return_value = 10

    monkeypatch.setattr(
        "paircars.pipeline.import_model.msmetadata",
        lambda: fake_msmd
    )
    fake_table = MagicMock()
    fake_table.colnames.return_value = ["DATA", "MODEL_DATA"]
    fake_table.getcol.side_effect = [
        np.array([0, 1]),  # ANTENNA1
        np.array([1, 2]),  # ANTENNA2
        np.ones((2, 4, 2), dtype=complex),  # model DATA
    ]
    monkeypatch.setattr(
        "paircars.pipeline.import_model.casatable",
        lambda: fake_table
    )
    monkeypatch.setattr(
        "paircars.pipeline.import_model.setjy",
        lambda **kwargs: None
    )
    monkeypatch.setattr("os.system", lambda x: 0)
    result = import_hyperdrive_model(msname, metafits)
    if raise_error:
        assert result == 1
    else:
        assert result == 0


@patch("paircars.pipeline.import_model.clean_shutdown")
@patch("paircars.pipeline.import_model.drop_cache")
@patch("paircars.pipeline.import_model.time.sleep", return_value=None)
@patch("paircars.pipeline.import_model.os.system")
@patch("paircars.pipeline.import_model.os.makedirs")
@patch("paircars.pipeline.import_model.os.path.exists")
@patch("paircars.pipeline.import_model.get_cachedir", return_value="/tmp")
@patch("paircars.pipeline.import_model.save_pid")
@patch("paircars.pipeline.import_model.get_local_dask_cluster")
@patch("paircars.pipeline.import_model.scale_worker_and_wait")
@patch("paircars.pipeline.import_model.psutil.cpu_count", return_value=8)
@patch("paircars.pipeline.import_model.psutil.virtual_memory")
@patch("paircars.pipeline.import_model.np.load", return_value=("job", "pass"))
@patch("paircars.pipeline.import_model.init_logger")
@patch("paircars.pipeline.import_model.get_ms_size", return_value=1)
@patch("paircars.pipeline.import_model.import_hyperdrive_model")
@patch("paircars.pipeline.import_model.delayed")
def test_main_import_model(
    mock_delayed,
    mock_import_model,
    mock_get_ms_size,
    mock_init_logger,
    mock_np_load,
    mock_virtual_mem,
    mock_cpu_count,
    mock_scale,
    mock_get_cluster,
    mock_save_pid,
    mock_get_cachedir,
    mock_path_exists,
    mock_makedirs,
    mock_os_system,
    mock_sleep,
    mock_drop_cache,
    mock_shutdown,
):
    def fake_exists(path):
        if "jobname_password.npy" in path:
            return True
        if "log.txt" in path:
            return True
        return True

    mock_path_exists.side_effect = fake_exists
    mock_virtual_mem.return_value.available = 16 * 1024**3  # 16GB
    mock_delayed.side_effect = lambda fn: fn
    mock_client = MagicMock()
    mock_cluster = MagicMock()
    mock_client.compute.side_effect = lambda tasks: tasks
    mock_client.gather.side_effect = lambda tasks: tasks
    mock_get_cluster.return_value = (mock_client, mock_cluster, "/tmp/dask")
    mock_import_model.return_value = 0
    msg = main(
        mslist="mock1.ms,mock2.ms",
        metafits="meta.fits",
        workdir="/tmp/work",
    )
    assert msg == 0
    assert mock_import_model.call_count == 2
    mock_import_model.side_effect = [0, 1]

    msg = main(
        mslist="mock1.ms,mock2.ms",
        metafits="meta.fits",
        workdir="/tmp/work",
    )
    assert msg == 1

    msg = main(
        mslist="",
        metafits="meta.fits",
        workdir="/tmp/work",
    )
    assert msg == 1

    mock_import_model.side_effect = Exception("Simulated failure")
    msg = main(
        mslist="mock.ms",
        metafits="meta.fits",
        workdir="/tmp/work",
    )

    assert msg == 1
    mock_import_model.side_effect = None
    mock_import_model.return_value = 0

    msg = main(
        mslist="mock.ms",
        metafits="meta.fits",
        workdir="/tmp/work",
        start_remote_log=True,
        logfile="log.txt",
    )

    assert msg == 0
    mock_init_logger.assert_called()
    mock_import_model.return_value = 0

    external_client = MagicMock()

    msg = main(
        mslist="mock.ms",
        metafits="meta.fits",
        workdir="/tmp/work",
        dask_client=external_client,
    )

    assert msg == 0


@pytest.mark.parametrize(
    "argv, should_exit",
    [
        (["prog.py"], True),  # Missing required args
        (["prog.py", "mock.ms", "--workdir", "/mock/work"], False),  # Valid
    ],
)
@patch("paircars.pipeline.import_model.main", return_value=0)
@patch("paircars.pipeline.import_model.sys.exit")
@patch("paircars.pipeline.import_model.argparse.ArgumentParser.print_help")
def test_cli(
    mock_print_help,
    mock_exit,
    mock_main,
    argv,
    should_exit,
):
    with patch("sys.argv", argv):
        from paircars.pipeline import import_model

        result = import_model.cli()
        assert result == should_exit
