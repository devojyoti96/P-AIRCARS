import pytest
from unittest.mock import patch, MagicMock, call
from paircars.pipeline.basic_cal import *


@patch("paircars.pipeline.basic_cal.limit_threads")
@patch("paircars.pipeline.basic_cal.suppress_output")
@patch("casatasks.flagdata")
@patch("casatasks.bandpass")
def test_run_bandpass(
    mock_bandpass,
    mock_flagdata,
    mock_suppress_output,
    mock_limit_threads,
):
    """Test full bandpass calibration with mocks"""
    msname = "/mock/path/test.ms"
    workdir = "/mock/workdir"
    expected_caltable = "/mock/workdir/test.bcal"

    result = run_bandpass(
        msname,
        workdir,
        uvrange=">100lambda",
        refant="1",
        solint="int",
        combine="scan",
        gaintable=[],
        gainfield=[],
        interp=[],
        n_threads=2,
    )

    expected_calls = [call(n_threads=2)]
    mock_limit_threads.assert_has_calls(expected_calls)
    mock_bandpass.assert_called_once_with(
        vis=msname,
        caltable=expected_caltable,
        uvrange=">100lambda",
        refant="1",
        solint="int",
        combine="scan",
        gaintable=[],
        gainfield=[],
        interp=[],
    )
    mock_flagdata.assert_called_once_with(
        vis=expected_caltable,
        mode="rflag",
        datacolumn="CPARAM",
        flagbackup=False,
    )


@patch("paircars.pipeline.basic_cal.limit_threads")
@patch("paircars.pipeline.basic_cal.suppress_output")
@patch("paircars.pipeline.basic_cal.crossphasecal")
def test_run_crossphasecal(
    mock_crossphasecal,
    mock_suppress_output,
    mock_limit_threads,
):
    msname = "/mock/path/test.ms"
    workdir = "/mock/workdir"
    expected_caltable = "/mock/workdir/test.kcrosscal"
    gaintable = ["/mock/path/test.bcal"]

    result = run_crossphasecal(
        msname,
        workdir,
        uvrange=">100lambda",
        gaintable=gaintable,
        n_threads=2,
    )

    expected_calls = [call(n_threads=2)]
    mock_limit_threads.assert_has_calls(expected_calls)
    mock_crossphasecal.assert_called_once_with(
        msname,
        expected_caltable,
        uvrange=">100lambda",
        gaintable=gaintable[0],
        n_threads=2,
    )


@patch("casatasks.applycal")
@patch("paircars.pipeline.basic_cal.suppress_output")
@patch("paircars.pipeline.basic_cal.limit_threads")
def test_run_applycal(
    mock_limit_threads,
    mock_suppress_output,
    mock_applycal,
):
    msname = "mock.ms"
    field = "1"
    scan = "2"
    gaintable = ["mock.kcal", "mock.bcal"]
    gainfield = ["1", "1"]
    interp = ["", ""]
    calwt = [False, False]
    result = run_applycal(
        msname=msname,
        applymode="calonly",
        flagbackup=True,
        gaintable=gaintable,
        gainfield=gainfield,
        interp=interp,
        calwt=calwt,
        n_threads=2,
    )
    assert result is None
    mock_applycal.assert_called_once()


@pytest.mark.parametrize(
    "msg_return, expect_print",
    [
        (0, False),  # normal case
        (2, True),  # issue case
    ],
)
def test_run_postcal_flag(
    msg_return,
    expect_print,
    mocker,
    capsys,
):
    msname = "test.ms"
    m_limit = mocker.patch("paircars.pipeline.basic_cal.limit_threads")
    m_flag = mocker.patch(
        "paircars.pipeline.basic_cal.single_ms_flag",
        return_value=msg_return,
    )
    result = run_postcal_flag(
        msname=msname,
        datacolumn="residual",
        threshold=5.0,
        n_threads=4,
        mem_limit=8,
    )
    m_limit.assert_called_once_with(n_threads=4)
    m_flag.assert_called_once_with(
        msname=msname,
        badspw="",
        bad_ants_str="",
        datacolumn="residual",
        use_tfcrop=True,
        use_rflag=True,
        flagdimension="freqtime",
        flag_autocorr=False,
        threshold=5.0,
        n_threads=4,
        mem_limit=8,
    )
    captured = capsys.readouterr()
    if expect_print:
        assert "Issue in post-calibration flagging" in captured.out
    else:
        assert captured.out == ""
    assert result is None


@pytest.mark.parametrize(
    "npol, do_polcal, applysol, do_postcal_flag, bandpass_ok, raise_exc",
    [
        (4, True, True, True, True, False),  # Full success
        (2, True, True, True, True, False),  # Not full-pol
        (4, True, True, True, False, False),  # Bandpass fails
        (4, True, False, True, True, False),  # applysol=False
        (4, True, True, False, True, False),  # do_postcal_flag=False
        (4, True, True, True, True, True),  # Exception branch
    ],
)
def test_single_ms_cal_and_flag(
    npol,
    do_polcal,
    applysol,
    do_postcal_flag,
    bandpass_ok,
    raise_exc,
    tmp_path,
):
    msname = "/mock/test.ms"
    workdir = str(tmp_path)

    with (
        patch("paircars.pipeline.basic_cal.msmetadata") as m_msmd,
        patch("paircars.pipeline.basic_cal.os.path.exists") as m_exists,
        patch("paircars.pipeline.basic_cal.os.system"),
        patch("paircars.pipeline.basic_cal.run_bandpass") as m_bandpass,
        patch("paircars.pipeline.basic_cal.run_crossphasecal") as m_cross,
        patch("paircars.pipeline.basic_cal.run_applycal") as m_apply,
        patch("paircars.pipeline.basic_cal.run_postcal_flag") as m_postflag,
        patch("paircars.pipeline.basic_cal.do_flag_backup") as m_backup,
        patch("paircars.pipeline.basic_cal.drop_cache") as m_drop,
        patch("paircars.pipeline.basic_cal.time.sleep"),
    ):

        m_msmd_inst = MagicMock()
        m_msmd.return_value = m_msmd_inst
        m_msmd_inst.ncorrforpol.return_value = [npol]
        if raise_exc:
            m_bandpass.side_effect = Exception("boom")
        else:
            if bandpass_ok:
                m_bandpass.return_value = f"{workdir}/test_caltable.bcal"
            else:
                m_bandpass.return_value = None

        m_cross.return_value = f"{workdir}/test_caltable.kcrosscal"

        def exists_side_effect(path):
            if raise_exc:
                return False
            if not bandpass_ok:
                return False
            return True

        m_exists.side_effect = exists_side_effect

        result = single_ms_cal_and_flag(
            msname=msname,
            workdir=workdir,
            cal_round=1,
            refant="1",
            uvrange=">100lambda",
            do_polcal=do_polcal,
            applysol=applysol,
            do_postcal_flag=do_postcal_flag,
            flag_threshold=5.0,
            n_threads=2,
            mem_limit=1024,
        )

        m_drop.assert_called_once_with(msname)

        if raise_exc:
            assert result == []
            return

        if not bandpass_ok:
            assert result == []
            m_bandpass.assert_called_once()
            return
        m_bandpass.assert_called_once()
        if do_polcal and npol == 4:
            m_cross.assert_called_once()
        else:
            m_cross.assert_not_called()
        if applysol:
            m_apply.assert_called_once()
        else:
            m_apply.assert_not_called()
        if applysol and do_postcal_flag:
            m_backup.assert_called_once()
            m_postflag.assert_called_once()
        else:
            m_postflag.assert_not_called()
        assert isinstance(result, list)
        assert len(result) == 2


@pytest.mark.parametrize(
    "cpu_frac, mem_frac, mslist, results",
    [
        # Normal case
        (
            0.5,
            0.5,
            ["a.ms", "b.ms"],
            [["a.bcal", "a.kcal"], ["b.bcal", "b.kcal"]],
        ),
        # cpu_frac > 0.8 clamp
        (
            0.95,
            0.5,
            ["a.ms"],
            [["a.bcal", None]],
        ),
        # mem_frac > 0.8 clamp
        (
            0.5,
            0.95,
            ["a.ms"],
            [["a.bcal", "a.kcal"]],
        ),
        # Empty result branch
        (
            0.5,
            0.5,
            ["a.ms"],
            [[]],
        ),
        # All None branch
        (
            0.5,
            0.5,
            ["a.ms"],
            [[None, None]],
        ),
    ],
)
def test_single_round_cal_and_flag(cpu_frac, mem_frac, mslist, results):
    mslist = ["a.ms", "b.ms"]
    fake_client = MagicMock()

    with (
        patch("paircars.pipeline.basic_cal.psutil.cpu_count", return_value=16),
        patch("paircars.pipeline.basic_cal.psutil.virtual_memory") as m_mem,
        patch("paircars.pipeline.basic_cal.delayed", side_effect=lambda f: f),
        patch("paircars.pipeline.basic_cal.single_ms_cal_and_flag") as m_single,
    ):
        mem_mock = MagicMock()
        mem_mock.available = 32 * 1024**3  # 32 GB
        m_mem.return_value = mem_mock
        m_single.side_effect = [
            ["a.bcal", None],
            ["b.bcal", "b.kcal"],
        ]
        fake_client.compute.side_effect = lambda x: x
        fake_client.gather.side_effect = lambda x: x

        output = single_round_cal_and_flag(
            mslist,
            fake_client,
            "workdir",
            1,
            refant=1,
            uvrange="",
            do_polcal=True,
            applysol=True,
            do_postcal_flag=True,
            flag_threshold=5.0,
            cpu_frac=cpu_frac,
            mem_frac=mem_frac,
        )

        assert output == {
            "a.ms": ["a.bcal"],
            "b.ms": ["b.bcal", "b.kcal"],
        }


@pytest.mark.parametrize(
    "npol, keep_backup, raise_exc",
    [
        (4, True, False),  # full-pol, backup
        (2, False, False),  # non full-pol
        (4, False, True),  # exception branch
    ],
)
def test_run_basic_cal_rounds(npol, keep_backup, raise_exc):

    mslist = ["a.ms", "b.ms"]
    fake_client = MagicMock()

    with (
        patch("paircars.pipeline.basic_cal.msmetadata") as m_msmd,
        patch("paircars.pipeline.basic_cal.get_unflagged_antennas") as m_unflag,
        patch("paircars.pipeline.basic_cal.get_gleam_uvrange") as m_gleam,
        patch("paircars.pipeline.basic_cal.get_uvrange_exclude") as m_uv_ex,
        patch("paircars.pipeline.basic_cal.single_round_cal_and_flag") as m_single,
        patch("paircars.pipeline.basic_cal.flagsummary") as m_flagsummary,
        patch("paircars.pipeline.basic_cal.delayed", side_effect=lambda f: f),
        patch("paircars.pipeline.basic_cal.os.chdir"),
        patch("paircars.pipeline.basic_cal.os.makedirs"),
        patch("paircars.pipeline.basic_cal.os.path.exists", return_value=True),
        patch("paircars.pipeline.basic_cal.os.system"),
        patch("paircars.pipeline.basic_cal.traceback.print_exc"),
    ):
        with patch.dict("sys.modules", {"casatasks": MagicMock(flagdata=MagicMock())}):
            m_msmd_inst = MagicMock()
            m_msmd.return_value = m_msmd_inst
            m_msmd_inst.ncorrforpol.return_value = [npol]
            m_msmd_inst.antennaids.return_value = [0]

            m_unflag.return_value = (["ant1"], [0.1])
            m_gleam.return_value = "100~200lambda"
            m_uv_ex.return_value = ["<50lambda"]
            if raise_exc:
                m_single.side_effect = Exception("boom")
            else:
                m_single.return_value = {"a.ms": ["a.bcal"], "b.ms": ["b.bcal"]}
            fake_client.compute.side_effect = lambda x: x
            fake_client.gather.side_effect = lambda x: x
            status, caltables = run_basic_cal_rounds(
                mslist=mslist,
                dask_client=fake_client,
                workdir="/tmp",
                outdir="/tmp",
                refant="",
                uvrange="",
                keep_backup=keep_backup,
                perform_polcal=True,
            )
            if raise_exc:
                assert status == 1
                assert caltables == []
                return

            assert status == 0
            assert isinstance(caltables, list)
            if npol == 4:
                assert m_single.call_count == 3
            else:
                assert m_single.call_count == 2
            if keep_backup:
                assert m_single.called


@pytest.mark.parametrize(
    "start_remote_log, provide_dask, empty_ms, raise_exc",
    [
        (False, False, False, False),  # normal local cluster
        (True, False, False, False),  # remote log branch
        (False, True, False, False),  # external dask client
        (False, False, True, False),  # empty mslist branch
        (False, False, False, True),  # exception branch
    ],
)
def test_main(
    start_remote_log,
    provide_dask,
    empty_ms,
    raise_exc,
):

    ms_input = "" if empty_ms else "a.ms,b.ms"

    fake_client = MagicMock()
    fake_cluster = MagicMock()

    with (
        patch("paircars.pipeline.basic_cal.get_MWA_OBSID", return_value="123"),
        patch("paircars.pipeline.basic_cal.get_local_dask_cluster") as m_cluster,
        patch("paircars.pipeline.basic_cal.scale_worker_and_wait"),
        patch("paircars.pipeline.basic_cal.run_basic_cal_rounds") as m_run,
        patch("paircars.pipeline.basic_cal.drop_cache") as m_drop,
        patch("paircars.pipeline.basic_cal.clean_shutdown") as m_clean,
        patch("paircars.pipeline.basic_cal.init_logger") as m_logger,
        patch("paircars.pipeline.basic_cal.os.makedirs"),
        patch("paircars.pipeline.basic_cal.os.path.exists", return_value=True),
        patch("paircars.pipeline.basic_cal.os.system"),
        patch("paircars.pipeline.basic_cal.psutil.cpu_count", return_value=16),
        patch("paircars.pipeline.basic_cal.time.sleep"),
        patch("paircars.pipeline.basic_cal.traceback.print_exc"),
    ):

        if not provide_dask:
            m_cluster.return_value = (fake_client, fake_cluster, "/tmp/daskdir")

        if start_remote_log:
            with patch(
                "paircars.pipeline.basic_cal.np.load", return_value=("job", "pass")
            ):

                m_logger.return_value = MagicMock()

                m_run.return_value = (0, ["a.bcal"])

                msg = main(
                    mslist=ms_input,
                    workdir="",
                    outdir="",
                    start_remote_log=True,
                    logfile="/tmp/log.txt",
                    dask_client=None if not provide_dask else fake_client,
                )
        else:
            if raise_exc:
                m_run.side_effect = Exception("boom")
            else:
                m_run.return_value = (0, ["a.bcal"])

            msg = main(
                mslist=ms_input,
                workdir="",
                outdir="",
                start_remote_log=False,
                logfile=None,
                dask_client=None if not provide_dask else fake_client,
            )

        if raise_exc:
            assert msg == 1
        else:
            assert msg == 0
        if not empty_ms:
            assert m_run.called
        assert m_drop.called
        assert m_clean.called
        if not provide_dask:
            fake_client.close.assert_called()
            fake_cluster.close.assert_called()


@pytest.mark.parametrize(
    "argv_args, expect_main_called, expected_exit",
    [
        (["prog", "--mslist", "a.ms"], True, 0),
        (["prog"], False, 1),
    ],
)
@patch("paircars.pipeline.basic_cal.main", return_value=0)
@patch("paircars.pipeline.basic_cal.sys.exit")
@patch("paircars.pipeline.basic_cal.argparse.ArgumentParser.print_help")
def test_cli(
    mock_print_help,
    mock_exit,
    mock_main,
    argv_args,
    expect_main_called,
    expected_exit,
):
    with patch("sys.argv", argv_args):
        from paircars.pipeline import basic_cal

        result = basic_cal.cli()

        if expect_main_called:
            mock_main.assert_called()
        else:
            mock_print_help.assert_called()

        assert result == expected_exit
