import pytest
from unittest.mock import patch, MagicMock
from casatools import table
from paircars.pipeline.flagging import *


def test_single_ms_flag(dummy_msname):
    result = single_ms_flag(
        msname=dummy_msname,
        badspw="0:0;1",
        bad_ants_str="1,2",
        datacolumn="data",
        use_tfcrop=True,
        use_rflag=True,
        flagdimension="freqtime",
        flag_autocorr=True,
        n_threads=-1,
        mem_limit=-1,
    )
    assert result == 0
    tb = table()
    tb.open(dummy_msname, nomodify=False)
    flag = tb.getcol("FLAG")
    flag *= False
    tb.putcol("FLAG", flag)
    tb.flush()
    tb.close()
    os.system(f"rm -rf {dummy_msname}.flagversions")
    assert os.path.exists(f"{dummy_msname}.flagversions") == False


@patch("paircars.pipeline.flagging.flagsummary")
@patch("paircars.pipeline.flagging.single_ms_flag")
@patch("paircars.pipeline.flagging.do_flag_backup")
@patch("paircars.pipeline.flagging.get_mwa_bad_ants", return_value="1,2")
@patch("paircars.pipeline.flagging.get_bad_chans", return_value="0:1~10")
@patch("paircars.pipeline.flagging.suppress_output")
@patch("paircars.pipeline.flagging.os.chdir")
@patch("paircars.pipeline.flagging.delayed")
@patch("casatasks.flagdata",return_value=0)
def test_do_flagging(
    mock_flagdata,
    mock_delayed,
    mock_chdir,
    mock_suppress,
    mock_get_bad_chans,
    mock_get_bad_ants,
    mock_flag_backup,
    mock_single_ms_flag,
    mock_flagsummary,
    dummy_metafits,
):
    mock_delayed.side_effect = lambda fn: fn
    mock_client = MagicMock()
    mock_client.compute.side_effect = lambda tasks: tasks
    mock_client.gather.side_effect = lambda tasks: tasks
    mock_single_ms_flag.return_value = 0
    msg = do_flagging(
        mslist=["mock1.ms", "mock2.ms"],
        metafits=dummy_metafits,
        dask_client=mock_client,
        workdir="/tmp",
        outdir="/tmp/out",
    )
    msg = do_flagging(
        mslist=["mock.ms"],
        metafits=dummy_metafits,
        dask_client=mock_client,
        workdir="/tmp",
        outdir="/tmp/out",
        flag_bad_ants=False,
        flag_bad_spw=False,
    )
    mock_single_ms_flag.side_effect = Exception("Simulated failure")
    msg = do_flagging(
        mslist=["mock.ms"],
        metafits=dummy_metafits,
        dask_client=mock_client,
        workdir="/tmp",
        outdir="/tmp/out",
    )


@patch("paircars.pipeline.flagging.clean_shutdown")
@patch("paircars.pipeline.flagging.drop_cache")
@patch("paircars.pipeline.flagging.time.sleep", return_value=None)
@patch("paircars.pipeline.flagging.os.system")
@patch("paircars.pipeline.flagging.os.makedirs")
@patch("paircars.pipeline.flagging.os.path.exists")
@patch("paircars.pipeline.flagging.get_local_dask_cluster")
@patch("paircars.pipeline.flagging.scale_worker_and_wait")
@patch("paircars.pipeline.flagging.np.load", return_value=("job", "pass"))
@patch("paircars.pipeline.flagging.init_logger")
@patch("paircars.pipeline.flagging.do_flagging")
@patch("paircars.pipeline.flagging.os.chdir",return_value=True)
@patch("paircars.pipeline.flagging.get_ncoarse",return_value=1)
@patch("casatasks.flagdata",return_value=1)
def test_main_flagging(
    mock_flagdata,
    mock_ncoarse,
    mock_chdir,
    mock_do_flagging,
    mock_init_logger,
    mock_np_load,
    mock_scale,
    mock_get_cluster,
    mock_path_exists,
    mock_makedirs,
    mock_os_system,
    mock_sleep,
    mock_drop_cache,
    mock_shutdown,
    dummy_metafits,
):
    def fake_exists(path):
        if "jobname_password.npy" in path:
            return True
        if "log.txt" in path:
            return True
        return True

    mock_path_exists.side_effect = fake_exists
    mock_client = MagicMock()
    mock_cluster = MagicMock()
    mock_get_cluster.return_value = (mock_client, mock_cluster, "/tmp/dask", 1)
    mock_do_flagging.return_value = (0, 2, 0)
    msg, succeed, failed = main(
        mslist="mock1.ms,mock2.ms",
        metafits=dummy_metafits,
        workdir="/tmp/work",
        outdir="/tmp/out",
    )
    assert msg == 0
    mock_do_flagging.assert_called_once()
    mock_do_flagging.return_value = (1, 0, 1)
    msg, succeed, failed = main(
        mslist="mock.ms",
        metafits=dummy_metafits,
        workdir="/tmp/work",
        outdir="/tmp/out",
    )
    assert msg == 1
    msg, succeed, failed = main(
        mslist="",
        metafits=dummy_metafits,
        workdir="/tmp/work",
        outdir="/tmp/out",
    )
    assert msg == 1
    mock_do_flagging.side_effect = Exception("Simulated failure")
    msg, succeed, failed = main(
        mslist="mock.ms",
        metafits=dummy_metafits,
        workdir="/tmp/work",
        outdir="/tmp/out",
    )
    assert msg == 1
    mock_do_flagging.side_effect = None
    mock_do_flagging.return_value = (0, 1, 0)
    msg, succeed, failed = main(
        mslist="mock.ms",
        metafits=dummy_metafits,
        workdir="/tmp/work",
        outdir="/tmp/out",
        start_remote_log=True,
        logfile="log.txt",
    )
    assert msg == 0
    mock_init_logger.assert_called()


@pytest.mark.parametrize(
    "argv_args, expect_main_called, expected_exit",
    [
        (["prog"], False, 1),  # No args: expect sys.exit(1)
        (
            ["prog", "mock.ms", "--workdir", "mockdir"],
            True,
            0,  # Valid: expect main() call and return value
        ),
    ],
)
@patch("paircars.pipeline.flagging.main", return_value= (0, 1, 0))
@patch("paircars.pipeline.flagging.sys.exit")
@patch("paircars.pipeline.flagging.argparse.ArgumentParser.print_help")
def test_cli_flagging(
    mock_print_help,
    mock_exit,
    mock_main,
    argv_args,
    expect_main_called,
    expected_exit,
):
    with patch("sys.argv", argv_args):
        from paircars.pipeline import flagging

        result = flagging.cli()
        assert result == expected_exit
