import pytest
from unittest.mock import patch, MagicMock
from paircars.pipeline.do_apply_selfcal import *


@pytest.mark.parametrize(
    "has_tables, valid_ms, coarse_match, dask_sum, raise_exc",
    [
        (False, True, True, 0, False),  # no selfcal tables
        (True, False, True, 0, False),  # no valid MS
        (True, True, False, 0, False),  # no matching coarse
        (True, True, True, 0, False),  # success
        (True, True, True, 1, False),  # failure
        (True, True, True, 0, True),  # exception branch
    ],
)
def test_run_all_applysol_selfcal(
    has_tables,
    valid_ms,
    coarse_match,
    dask_sum,
    raise_exc,
):
    mslist = ["a.ms"]
    fake_client = MagicMock()
    with (
        patch("paircars.pipeline.do_apply_selfcal.psutil.cpu_count", return_value=8),
        patch("paircars.pipeline.do_apply_selfcal.psutil.virtual_memory") as m_mem,
        patch("paircars.pipeline.do_apply_selfcal.os.chdir"),
        patch("paircars.pipeline.do_apply_selfcal.np.unique", return_value=mslist),
        patch("paircars.pipeline.do_apply_selfcal.fits.getheader") as m_header,
        patch("paircars.pipeline.do_apply_selfcal.glob.glob") as m_glob,
        patch("paircars.pipeline.do_apply_selfcal.check_datacolumn_valid") as m_check,
        patch("paircars.pipeline.do_apply_selfcal.msmetadata") as m_msmd,
        patch("paircars.pipeline.do_apply_selfcal.freq_to_MWA_coarse", return_value=10),
        patch("paircars.pipeline.do_apply_selfcal.delayed", side_effect=lambda f: f),
        patch("paircars.pipeline.do_apply_selfcal.applysol") as m_apply,
        patch("paircars.pipeline.do_apply_selfcal.os.system") as m_system,
        patch("paircars.pipeline.do_apply_selfcal.traceback.print_exc"),
    ):
        mem_mock = MagicMock()
        mem_mock.available = 16 * 1024**3
        m_mem.return_value = mem_mock
        m_header.return_value = {"GPSTIME": "123"}
        if not has_tables:
            m_glob.return_value = []
        else:
            m_glob.side_effect = [
                ["/cal/selfcal_123_coarsechan_0_20.gcal"],  # gcal
                ["/cal/selfcal_123_coarsechan_0_20.bcal"],  # bcal
                ["/cal/selfcal_123_coarsechan_0_20.dcal"],  # dcal
            ]
        m_check.return_value = valid_ms
        msmd_inst = MagicMock()
        m_msmd.return_value = msmd_inst
        msmd_inst.chanfreqs.return_value = [150]
        msmd_inst.open.return_value = None
        msmd_inst.close.return_value = None
        if not coarse_match:
            with patch(
                "paircars.pipeline.do_apply_selfcal.freq_to_MWA_coarse",
                return_value=1000,
            ):
                pass
        if raise_exc:
            fake_client.compute.side_effect = Exception("boom")
        else:
            fake_client.compute.side_effect = lambda x: x
            fake_client.gather.side_effect = lambda x: [0] if dask_sum == 0 else [1]
        result = run_all_applysol(
            mslist=mslist,
            metafits="meta.fits",
            dask_client=fake_client,
            workdir="/tmp",
            caldir="/cal",
        )
        if raise_exc:
            assert result == 1
            m_system.assert_any_call("rm -rf casa*log")
            return

        if not has_tables:
            assert result == 1
            return

        if not valid_ms:
            assert result == 1
            return

        if not coarse_match:
            assert result == 1
            return


@pytest.mark.parametrize(
    "start_remote_log, provide_dask, caldir_exists, raise_exc",
    [
        (False, False, True, False),  # normal local dask
        (True, False, True, False),  # remote logging
        (False, True, True, False),  # external dask client
        (False, False, False, False),  # invalid caldir
        (False, False, True, True),  # exception branch
    ],
)
def test_main_apply_selfcal(
    start_remote_log,
    provide_dask,
    caldir_exists,
    raise_exc,
):
    ms_input = "a.ms,b.ms"
    fake_client = MagicMock()
    fake_cluster = MagicMock()
    with (
        patch("paircars.pipeline.do_apply_selfcal.get_cachedir", return_value="/tmp"),
        patch("paircars.pipeline.do_apply_selfcal.save_pid"),
        patch("paircars.pipeline.do_apply_selfcal.get_local_dask_cluster") as m_cluster,
        patch("paircars.pipeline.do_apply_selfcal.scale_worker_and_wait"),
        patch("paircars.pipeline.do_apply_selfcal.run_all_applysol") as m_run,
        patch("paircars.pipeline.do_apply_selfcal.drop_cache") as m_drop,
        patch("paircars.pipeline.do_apply_selfcal.clean_shutdown") as m_clean,
        patch("paircars.pipeline.do_apply_selfcal.init_logger") as m_logger,
        patch("paircars.pipeline.do_apply_selfcal.os.makedirs"),
        patch("paircars.pipeline.do_apply_selfcal.os.path.exists") as m_exists,
        patch("paircars.pipeline.do_apply_selfcal.os.system"),
        patch("paircars.pipeline.do_apply_selfcal.psutil.cpu_count", return_value=16),
        patch("paircars.pipeline.do_apply_selfcal.time.sleep"),
        patch("paircars.pipeline.do_apply_selfcal.traceback.print_exc"),
    ):
        m_cluster.return_value = (fake_client, fake_cluster, "/tmp/daskdir")

        def exists_side_effect(path):
            if "jobname_password.npy" in path:
                return start_remote_log
            if path == "/cal":
                return caldir_exists
            if path == "/tmp/log.txt":
                return True
            return True

        m_exists.side_effect = exists_side_effect
        if start_remote_log:
            with patch(
                "paircars.pipeline.do_apply_selfcal.np.load",
                return_value=("job", "pass"),
            ):
                m_logger.return_value = MagicMock()
                m_run.return_value = 0

                result = main(
                    mslist=ms_input,
                    metafits="meta.fits",
                    workdir="",
                    caldir="/cal",
                    start_remote_log=True,
                    logfile="/tmp/log.txt",
                    dask_client=None if not provide_dask else fake_client,
                )
        else:
            if raise_exc:
                m_run.side_effect = Exception("boom")
            else:
                m_run.return_value = 0

            result = main(
                mslist=ms_input,
                metafits="meta.fits",
                workdir="",
                caldir="/cal",
                start_remote_log=False,
                logfile=None,
                dask_client=None if not provide_dask else fake_client,
            )
        if not caldir_exists:
            assert result == 1
            return
        if raise_exc:
            assert result == 1
        else:
            assert result == 0
        if caldir_exists:
            assert m_run.called
        assert m_drop.called
        assert m_clean.called
        if not provide_dask:
            fake_client.close.assert_called()
            fake_cluster.close.assert_called()


@pytest.mark.parametrize(
    "argv, should_exit",
    [
        (["prog.py"], True),
        (
            [
                "prog.py",
                "ms1.ms,ms2.ms",
                "--workdir",
                "/mock/work",
                "--caldir",
                "/mock/caldir",
                "--overwrite_datacolumn",
                "--force_apply",
                "--cpu_frac",
                "0.6",
                "--mem_frac",
                "0.7",
                "--jobid",
                "321",
            ],
            False,
        ),
    ],
)
@patch("paircars.pipeline.do_apply_selfcal.main", return_value=0)
@patch("paircars.pipeline.do_apply_selfcal.sys.exit")
@patch("paircars.pipeline.do_apply_selfcal.argparse.ArgumentParser.print_help")
def test_cli_apply_selfcal(mock_print_help, mock_exit, mock_main, argv, should_exit):
    with patch("sys.argv", argv):
        from paircars.pipeline import do_apply_selfcal

        result = do_apply_selfcal.cli()
        assert result == should_exit
