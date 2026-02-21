import pytest
from unittest.mock import patch, MagicMock
from paircars.pipeline.do_target_split import *


def test_chanlist_to_str():
    result = chanlist_to_str([0, 1, 2, 10, 45])
    assert result == "0~2;10;45"


@patch("paircars.pipeline.do_target_split.drop_cache")
@patch("paircars.pipeline.do_target_split.time.sleep", return_value=None)
@patch("paircars.pipeline.do_target_split.os.chdir")
@patch("paircars.pipeline.do_target_split.os.path.exists", return_value=False)
@patch("paircars.pipeline.do_target_split.get_timeranges", return_value=["0~60"])
@patch("paircars.pipeline.do_target_split.get_MWA_coarse_bands")
@patch("paircars.pipeline.do_target_split.single_mstransform")
@patch("paircars.pipeline.do_target_split.msmetadata")
@patch("paircars.pipeline.do_target_split.psutil.virtual_memory")
@patch("paircars.pipeline.do_target_split.psutil.cpu_count", return_value=8)
@patch("paircars.pipeline.do_target_split.delayed")
def test_split_target_scans(
    mock_delayed,
    mock_cpu,
    mock_virtual_mem,
    mock_msmetadata,
    mock_single_mstransform,
    mock_coarse_bands,
    mock_timeranges,
    mock_exists,
    mock_chdir,
    mock_sleep,
    mock_drop_cache,
):
    mock_delayed.side_effect = lambda fn: fn
    mock_virtual_mem.return_value.available = 16 * 1024**3  # 16GB
    mock_msmd = MagicMock()
    mock_msmd.open.return_value = None
    mock_msmd.close.return_value = None
    mock_msmd.chanres.return_value = [0.04]
    mock_msmd.chanfreqs.return_value = [100, 101, 102]
    mock_msmd.nchan.return_value = 3
    mock_msmetadata.return_value = mock_msmd
    mock_client = MagicMock()
    mock_client.compute.side_effect = lambda tasks: tasks
    mock_client.gather.side_effect = lambda tasks: tasks
    mock_coarse_bands.return_value = [(0, 5), (5, 10)]
    mock_single_mstransform.side_effect = lambda **kwargs: kwargs["outputms"]
    status, result = split_target_scans(
        msname="mock.ms",
        dask_client=mock_client,
        workdir="/tmp",
        timeres=10,
        freqres=1.0,
        datacolumn="DATA",
    )
    assert status == 0
    assert len(result) == 2
    mock_exists.return_value = True
    mock_coarse_bands.return_value = [(0, 5)]
    status, result = split_target_scans(
        msname="mock.ms",
        dask_client=mock_client,
        workdir="/tmp",
        timeres=10,
        freqres=1.0,
        datacolumn="DATA",
    )
    assert status == 0
    assert len(result) == 1
    mock_exists.return_value = False
    mock_coarse_bands.return_value = []
    status, result = split_target_scans(
        msname="mock.ms",
        dask_client=mock_client,
        workdir="/tmp",
        timeres=10,
        freqres=1.0,
        datacolumn="DATA",
    )
    assert status == 0
    mock_msmd.chanres.return_value = [5.0]  # larger than freqres
    mock_coarse_bands.return_value = [(0, 5)]
    status, result = split_target_scans(
        msname="mock.ms",
        dask_client=mock_client,
        workdir="/tmp",
        timeres=10,
        freqres=1.0,
        datacolumn="DATA",
    )
    assert status == 0
    mock_single_mstransform.side_effect = Exception("Simulated failure")
    status, result = split_target_scans(
        msname="mock.ms",
        dask_client=mock_client,
        workdir="/tmp",
        timeres=10,
        freqres=1.0,
        datacolumn="DATA",
    )
    assert result == []


@patch("paircars.pipeline.do_target_split.clean_shutdown")
@patch("paircars.pipeline.do_target_split.drop_cache")
@patch("paircars.pipeline.do_target_split.time.sleep", return_value=None)
@patch("paircars.pipeline.do_target_split.os.system")
@patch("paircars.pipeline.do_target_split.os.makedirs")
@patch("paircars.pipeline.do_target_split.os.path.exists", return_value=True)
@patch("paircars.pipeline.do_target_split.get_local_dask_cluster")
@patch("paircars.pipeline.do_target_split.scale_worker_and_wait")
@patch("paircars.pipeline.do_target_split.psutil.cpu_count", return_value=8)
@patch("paircars.pipeline.do_target_split.init_logger")
@patch("paircars.pipeline.do_target_split.split_target_scans")
def test_main_split_target_scans(
    mock_split_target_scans,
    mock_init_logger,
    mock_cpu_count,
    mock_scale,
    mock_get_cluster,
    mock_path_exists,
    mock_makedirs,
    mock_os_system,
    mock_sleep,
    mock_drop_cache,
    mock_shutdown,
):
    mock_client = MagicMock()
    mock_cluster = MagicMock()
    mock_get_cluster.return_value = (mock_client, mock_cluster, "/tmp/dask")
    mock_split_target_scans.return_value = (0, ["chunk1.ms", "chunk2.ms"])
    msg = main(
        mslist="mock1.ms,mock2.ms",
        workdir="/tmp/work",
        datacolumn="DATA",
    )
    assert msg == 0
    assert mock_split_target_scans.call_count == 2
    mock_split_target_scans.return_value = (1, [])
    msg = main(
        mslist="mock.ms",
        workdir="/tmp/work",
        datacolumn="DATA",
    )
    assert msg == 1
    msg = main(
        mslist="",
        workdir="/tmp/work",
        datacolumn="DATA",
    )
    assert msg == 1
    mock_split_target_scans.side_effect = Exception("Simulated failure")
    msg = main(
        mslist="mock.ms",
        workdir="/tmp/work",
        datacolumn="DATA",
    )
    assert msg == 1


@pytest.mark.parametrize(
    "argv, should_exit",
    [
        (["prog.py"], True),  # Missing args
        (
            [
                "prog.py",
                "mock.ms",
                "--workdir",
                "/mock/work",
                "--scans",
                "1,2",
                "--prefix",
                "targets",
            ],
            False,
        ),  # Normal CLI call
    ],
)
@patch("paircars.pipeline.do_target_split.main", return_value=0)
@patch("paircars.pipeline.do_target_split.sys.exit")
@patch("paircars.pipeline.do_target_split.argparse.ArgumentParser.print_help")
def test_cli_split_target_scans(
    mock_print_help,
    mock_exit,
    mock_main,
    argv,
    should_exit,
):
    with patch("sys.argv", argv):
        from paircars.pipeline import do_target_split

        result = do_target_split.cli()
        assert result == should_exit
