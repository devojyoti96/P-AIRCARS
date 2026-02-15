import pytest
from unittest.mock import patch, MagicMock
from paircars.pipeline.master_flow import *
from paircars.pipeline.init_data import init_paircars_data


@pytest.mark.parametrize("mock_msg,raises", [(0, False), (1, True)])
@patch("paircars.pipeline.master_flow.get_run_context")
@patch("paircars.pipeline.master_flow.start_log_task_saver")
@patch("paircars.pipeline.master_flow.os.makedirs")
@patch("paircars.pipeline.master_flow.os.path.exists")
@patch("paircars.pipeline.master_flow.os.remove")
@patch("paircars.pipeline.move_solarcenter.main")
@patch("paircars.pipeline.master_flow.get_dask_client")
def test_run_solar_phasecenter_jobs(
    mock_dask_client,
    mock_main,
    mock_remove,
    mock_exists,
    mock_makedirs,
    mock_log_task_saver,
    mock_get_ctx,
    mock_msg,
    raises,
):
    mock_main.return_value = mock_msg
    mock_exists.return_value = True
    mock_ctx = MagicMock()
    mock_ctx.task_run.id = "abc123"
    mock_ctx.task_run.name = "mock_task"
    mock_get_ctx.return_value = mock_ctx
    # Mock log thread
    mock_thread = MagicMock()
    mock_log_task_saver.return_value = mock_thread
    # Mock Dask client context manager
    mock_client = MagicMock()
    mock_dask_client.return_value.__enter__.return_value = mock_client
    # Parameters
    kwargs = dict(
        mslist="mock1.ms,mock2.ms",
        workdir="/mock/workdir",
        prefix="target",
        jobid=42,
        cpu_frac=0.5,
        mem_frac=0.5,
        remote_log=True,
    )
    if raises:
        with pytest.raises(
            RuntimeError, match="Moving phasecenter to solar center is failed."
        ):
            run_solar_phasecenter_jobs.fn(**kwargs)
    else:
        result = run_solar_phasecenter_jobs.fn(**kwargs)
        assert result == 0
    mock_makedirs.assert_called_once_with("/mock/workdir/logs", exist_ok=True)
    # Logfile removed if exists
    mock_remove.assert_called_once_with("/mock/workdir/logs/cor_phasecenter_target.log")
    # Log saver started
    mock_log_task_saver.assert_called_once()
    # Thread join always happens
    mock_thread.join.assert_called_once_with(timeout=5)
    # move_solarcenter.main called
    mock_main.assert_called_once()


import pytest
from unittest.mock import patch, MagicMock


@pytest.mark.parametrize("mock_msg,raises", [(0, False), (1, True)])
@patch("paircars.pipeline.master_flow.get_run_context")
@patch("paircars.pipeline.master_flow.start_log_task_saver")
@patch("paircars.pipeline.master_flow.os.makedirs")
@patch("paircars.pipeline.master_flow.os.path.exists")
@patch("paircars.pipeline.master_flow.os.remove")
@patch("paircars.pipeline.mwa_make_ds.main")
@patch("paircars.pipeline.master_flow.get_dask_client")
def test_run_ds_jobs(
    mock_dask_client,
    mock_main,
    mock_remove,
    mock_exists,
    mock_makedirs,
    mock_log_task_saver,
    mock_get_ctx,
    mock_msg,
    raises,
):
    mock_main.return_value = mock_msg
    mock_exists.return_value = True

    mock_ctx = MagicMock()
    mock_ctx.task_run.id = "abc123"
    mock_ctx.task_run.name = "mock_task"
    mock_get_ctx.return_value = mock_ctx

    mock_thread = MagicMock()
    mock_log_task_saver.return_value = mock_thread

    mock_client = MagicMock()
    mock_dask_client.return_value.__enter__.return_value = mock_client

    kwargs = dict(
        mslist="mock1.ms,mock2.ms",
        metafits="mock.metafits",
        workdir="/mock/workdir",
        outdir="/mock/outdir",
        plot_quantity="TB",
        jobid=42,
        cpu_frac=0.5,
        mem_frac=0.5,
        remote_log=True,
    )

    if raises:
        with pytest.raises(RuntimeError, match="Dynamic spectrum making is failed."):
            run_ds_jobs.fn(**kwargs)
    else:
        result = run_ds_jobs.fn(**kwargs)
        assert result == 0

    mock_makedirs.assert_called_once_with("/mock/workdir/logs", exist_ok=True)
    mock_remove.assert_called_once_with("/mock/workdir/logs/ds_target.log")
    mock_log_task_saver.assert_called_once()
    mock_thread.join.assert_called_once_with(timeout=5)

    mock_main.assert_called_once_with(
        "mock1.ms,mock2.ms",
        "mock.metafits",
        "/mock/workdir",
        "/mock/outdir",
        plot_quantity="TB",
        cpu_frac=0.5,
        mem_frac=0.5,
        logfile="/mock/workdir/logs/ds_target.log",
        jobid=42,
        start_remote_log=True,
        dask_client=mock_client,
    )


@pytest.mark.parametrize("mock_msg,raises", [(0, False), (1, True)])
@patch("paircars.pipeline.master_flow.get_run_context")
@patch("paircars.pipeline.master_flow.start_log_task_saver")
@patch("paircars.pipeline.master_flow.os.makedirs")
@patch("paircars.pipeline.master_flow.os.path.exists")
@patch("paircars.pipeline.master_flow.os.remove")
@patch("paircars.pipeline.do_target_split.main")
@patch("paircars.pipeline.master_flow.get_dask_client")
def test_run_target_split_jobs(
    mock_dask_client,
    mock_main,
    mock_remove,
    mock_exists,
    mock_makedirs,
    mock_log_task_saver,
    mock_get_ctx,
    mock_msg,
    raises,
):
    mock_main.return_value = mock_msg
    mock_exists.return_value = True

    mock_ctx = MagicMock()
    mock_ctx.task_run.id = "abc123"
    mock_ctx.task_run.name = "mock_task"
    mock_get_ctx.return_value = mock_ctx

    mock_thread = MagicMock()
    mock_log_task_saver.return_value = mock_thread

    mock_client = MagicMock()
    mock_dask_client.return_value.__enter__.return_value = mock_client

    kwargs = dict(
        mslist="mock1.ms,mock2.ms",
        workdir="/mock/workdir",
        datacolumn="data",
        timeres=2.0,
        freqres=0.1,
        prefix="target",
        time_window=10,
        time_interval=5,
        quack_timestamps=2,
        jobid=42,
        cpu_frac=0.5,
        mem_frac=0.5,
        remote_log=True,
        force_split=True,
    )

    if raises:
        with pytest.raises(
            RuntimeError,
            match="Spliting measurement set into coarse channels is failed.",
        ):
            run_target_split_jobs.fn(**kwargs)
    else:
        result = run_target_split_jobs.fn(**kwargs)
        assert result == 0

    mock_makedirs.assert_called_once_with("/mock/workdir/logs", exist_ok=True)
    mock_remove.assert_called_once_with("/mock/workdir/logs/split_target.log")
    mock_log_task_saver.assert_called_once()
    mock_thread.join.assert_called_once_with(timeout=5)

    mock_main.assert_called_once_with(
        "mock1.ms,mock2.ms",
        workdir="/mock/workdir",
        datacolumn="data",
        time_window=10,
        time_interval=5,
        freqres=0.1,
        timeres=2.0,
        quack_timestamps=2,
        prefix="target",
        cpu_frac=0.5,
        mem_frac=0.5,
        logfile="/mock/workdir/logs/split_target.log",
        jobid=42,
        start_remote_log=True,
        force_split=True,
        dask_client=mock_client,
    )


@pytest.mark.parametrize(
    "flag_calibrators,mock_msg,raises",
    [
        (True, 0, False),
        (False, 0, False),
        (True, 1, True),
    ],
)
@patch("paircars.pipeline.master_flow.get_run_context")
@patch("paircars.pipeline.master_flow.start_log_task_saver")
@patch("paircars.pipeline.master_flow.os.makedirs")
@patch("paircars.pipeline.master_flow.os.path.exists")
@patch("paircars.pipeline.master_flow.os.remove")
@patch("paircars.pipeline.flagging.main")
@patch("paircars.pipeline.master_flow.get_dask_client")
def test_run_flag(
    mock_dask_client,
    mock_main,
    mock_remove,
    mock_exists,
    mock_makedirs,
    mock_log_task_saver,
    mock_get_ctx,
    flag_calibrators,
    mock_msg,
    raises,
):
    mock_main.return_value = mock_msg
    mock_exists.return_value = True

    mock_ctx = MagicMock()
    mock_ctx.task_run.id = "abc123"
    mock_ctx.task_run.name = "mock_task"
    mock_get_ctx.return_value = mock_ctx

    mock_thread = MagicMock()
    mock_log_task_saver.return_value = mock_thread

    mock_client = MagicMock()
    mock_dask_client.return_value.__enter__.return_value = mock_client

    kwargs = dict(
        mslist="mock1.ms,mock2.ms",
        metafits="mock.metafits",
        workdir="/mock/workdir",
        outdir="/mock/outdir",
        datacolumn="DATA",
        flag_calibrators=flag_calibrators,
        flag_quack=True,
        restore_flag=True,
        jobid=42,
        cpu_frac=0.5,
        mem_frac=0.5,
        remote_log=True,
    )

    if raises:
        with pytest.raises(RuntimeError, match="Flagging is failed."):
            run_flag.fn(**kwargs)
    else:
        result = run_flag.fn(**kwargs)
        assert result == 0

    if flag_calibrators:
        expected_log = "/mock/workdir/logs/flagging_cal_calibrator.log"
        expected_flagdimension = "freqtime"
        expected_use_tfcrop = True
    else:
        expected_log = "/mock/workdir/logs/flagging_target_target.log"
        expected_flagdimension = "freq"
        expected_use_tfcrop = False

    mock_makedirs.assert_called_once_with("/mock/workdir/logs", exist_ok=True)
    mock_remove.assert_called_once_with(expected_log)
    mock_log_task_saver.assert_called_once()
    mock_thread.join.assert_called_once_with(timeout=5)

    mock_main.assert_called_once_with(
        "mock1.ms,mock2.ms",
        "mock.metafits",
        workdir="/mock/workdir",
        outdir="/mock/outdir",
        datacolumn="DATA",
        flag_bad_ants=True,
        flag_bad_spw=True,
        use_tfcrop=expected_use_tfcrop,
        flag_autocorr=True,
        flag_quack=True,
        flagdimension=expected_flagdimension,
        restore_flag=True,
        flagbackup=False,
        cpu_frac=0.5,
        mem_frac=0.5,
        logfile=expected_log,
        jobid=42,
        start_remote_log=True,
        dask_client=mock_client,
    )


@pytest.mark.parametrize("mock_msg,raises", [(0, False), (1, True)])
@patch("paircars.pipeline.master_flow.get_run_context")
@patch("paircars.pipeline.master_flow.start_log_task_saver")
@patch("paircars.pipeline.master_flow.os.makedirs")
@patch("paircars.pipeline.master_flow.os.path.exists")
@patch("paircars.pipeline.master_flow.os.remove")
@patch("paircars.pipeline.import_model.main")
@patch("paircars.pipeline.master_flow.get_dask_client")
def test_run_import_model(
    mock_dask_client,
    mock_main,
    mock_remove,
    mock_exists,
    mock_makedirs,
    mock_log_task_saver,
    mock_get_ctx,
    mock_msg,
    raises,
):
    mock_main.return_value = mock_msg
    mock_exists.return_value = True

    mock_ctx = MagicMock()
    mock_ctx.task_run.id = "abc123"
    mock_ctx.task_run.name = "mock_task"
    mock_get_ctx.return_value = mock_ctx

    mock_thread = MagicMock()
    mock_log_task_saver.return_value = mock_thread

    mock_client = MagicMock()
    mock_dask_client.return_value.__enter__.return_value = mock_client

    kwargs = dict(
        mslist="mock1.ms,mock2.ms",
        metafits="mock.metafits",
        workdir="/mock/workdir",
        jobid=42,
        cpu_frac=0.5,
        mem_frac=0.5,
        remote_log=True,
    )

    if raises:
        with pytest.raises(
            RuntimeError,
            match="Importing calibrator model is failed.",
        ):
            run_import_model.fn(**kwargs)
    else:
        result = run_import_model.fn(**kwargs)
        assert result == 0

    mock_makedirs.assert_called_once_with("/mock/workdir/logs", exist_ok=True)
    mock_remove.assert_called_once_with("/mock/workdir/logs/modeling.log")
    mock_log_task_saver.assert_called_once()
    mock_thread.join.assert_called_once_with(timeout=5)

    mock_main.assert_called_once_with(
        "mock1.ms,mock2.ms",
        "mock.metafits",
        "/mock/workdir",
        cpu_frac=0.5,
        mem_frac=0.5,
        logfile="/mock/workdir/logs/modeling.log",
        jobid=42,
        start_remote_log=True,
        dask_client=mock_client,
    )


@pytest.mark.parametrize(
    "mock_msg,raises",
    [
        (0, False),
        (1, True),
    ],
)
@patch("paircars.pipeline.master_flow.get_run_context")
@patch("paircars.pipeline.master_flow.start_log_task_saver")
@patch("paircars.pipeline.master_flow.os.makedirs")
@patch("paircars.pipeline.master_flow.os.path.exists")
@patch("paircars.pipeline.master_flow.os.remove")
@patch("paircars.pipeline.basic_cal.main")
@patch("paircars.pipeline.master_flow.get_dask_client")
def test_run_basic_cal_jobs(
    mock_dask_client,
    mock_main,
    mock_remove,
    mock_exists,
    mock_makedirs,
    mock_log_task_saver,
    mock_get_ctx,
    mock_msg,
    raises,
):
    mock_main.return_value = mock_msg
    mock_exists.return_value = True

    mock_ctx = MagicMock()
    mock_ctx.task_run.id = "abc123"
    mock_ctx.task_run.name = "mock_task"
    mock_get_ctx.return_value = mock_ctx

    mock_thread = MagicMock()
    mock_log_task_saver.return_value = mock_thread

    mock_client = MagicMock()
    mock_dask_client.return_value.__enter__.return_value = mock_client

    kwargs = dict(
        mslist="mock1.ms,mock2.ms",
        workdir="/mock/workdir",
        outdir="/mock/outdir",
        perform_polcal=True,
        jobid=42,
        cpu_frac=0.5,
        mem_frac=0.5,
        keep_backup=True,
        remote_log=True,
    )

    if raises:
        with pytest.raises(RuntimeError, match="Basic calibration is failed."):
            run_basic_cal_jobs.fn(**kwargs)
    else:
        result = run_basic_cal_jobs.fn(**kwargs)
        assert result == 0

    mock_makedirs.assert_called_once_with("/mock/workdir/logs", exist_ok=True)
    mock_remove.assert_called_once_with("/mock/workdir/logs/basic_cal.log")
    mock_log_task_saver.assert_called_once()
    mock_thread.join.assert_called_once_with(timeout=5)

    mock_main.assert_called_once_with(
        "mock1.ms,mock2.ms",
        "/mock/workdir",
        "/mock/outdir",
        perform_polcal=True,
        keep_backup=True,
        start_remote_log=True,
        cpu_frac=0.5,
        mem_frac=0.5,
        logfile="/mock/workdir/logs/basic_cal.log",
        jobid=42,
        dask_client=mock_client,
    )


import pytest
from unittest.mock import patch, MagicMock


@pytest.mark.parametrize("mock_msg,raises", [(0, False), (1, True)])
@patch("paircars.pipeline.master_flow.get_run_context")
@patch("paircars.pipeline.master_flow.start_log_task_saver")
@patch("paircars.pipeline.master_flow.os.makedirs")
@patch("paircars.pipeline.master_flow.os.path.exists")
@patch("paircars.pipeline.master_flow.os.remove")
@patch("paircars.pipeline.do_apply_basiccal.main")
@patch("paircars.pipeline.master_flow.get_dask_client")
def test_run_apply_basiccal_sol(
    mock_dask_client,
    mock_main,
    mock_remove,
    mock_exists,
    mock_makedirs,
    mock_log_task_saver,
    mock_get_ctx,
    mock_msg,
    raises,
):
    mock_main.return_value = mock_msg
    mock_exists.return_value = True

    mock_ctx = MagicMock()
    mock_ctx.task_run.id = "abc123"
    mock_ctx.task_run.name = "mock_task"
    mock_get_ctx.return_value = mock_ctx

    mock_thread = MagicMock()
    mock_log_task_saver.return_value = mock_thread

    mock_client = MagicMock()
    mock_dask_client.return_value.__enter__.return_value = mock_client

    kwargs = dict(
        mslist="mock1.ms,mock2.ms",
        calibrator_metafits="cal.meta",
        target_metafits="tar.meta",
        workdir="/mock/workdir",
        caldir="/mock/caldir",
        overwrite_datacolumn=False,
        only_amplitude=True,
        applymode="calflag",
        prefix="target",
        jobid=42,
        cpu_frac=0.5,
        mem_frac=0.5,
        remote_log=True,
    )

    if raises:
        with pytest.raises(
            RuntimeError,
            match="Applying basic calibration solutions is failed.",
        ):
            run_apply_basiccal_sol.fn(**kwargs)
    else:
        result = run_apply_basiccal_sol.fn(**kwargs)
        assert result == 0

    mock_makedirs.assert_called_once_with("/mock/workdir/logs", exist_ok=True)
    mock_remove.assert_called_once_with("/mock/workdir/logs/apply_basiccal_target.log")
    mock_log_task_saver.assert_called_once()
    mock_thread.join.assert_called_once_with(timeout=5)

    mock_main.assert_called_once_with(
        "mock1.ms,mock2.ms",
        "cal.meta",
        "tar.meta",
        "/mock/workdir",
        "/mock/caldir",
        applymode="calflag",
        overwrite_datacolumn=False,
        only_amplitude=True,
        start_remote_log=True,
        cpu_frac=0.5,
        mem_frac=0.5,
        logfile="/mock/workdir/logs/apply_basiccal_target.log",
        jobid=42,
        dask_client=mock_client,
    )


@pytest.mark.parametrize("mock_msg,raises", [(0, False), (1, True)])
@patch("paircars.pipeline.master_flow.get_run_context")
@patch("paircars.pipeline.master_flow.start_log_task_saver")
@patch("paircars.pipeline.master_flow.os.makedirs")
@patch("paircars.pipeline.master_flow.os.path.exists")
@patch("paircars.pipeline.master_flow.os.remove")
@patch("paircars.pipeline.do_sidereal_cor.main")
@patch("paircars.pipeline.master_flow.get_dask_client")
def test_run_solar_siderealcor_jobs(
    mock_dask_client,
    mock_main,
    mock_remove,
    mock_exists,
    mock_makedirs,
    mock_log_task_saver,
    mock_get_ctx,
    mock_msg,
    raises,
):
    mock_main.return_value = mock_msg
    mock_exists.return_value = True

    mock_ctx = MagicMock()
    mock_ctx.task_run.id = "abc123"
    mock_ctx.task_run.name = "mock_task"
    mock_get_ctx.return_value = mock_ctx

    mock_thread = MagicMock()
    mock_log_task_saver.return_value = mock_thread

    mock_client = MagicMock()
    mock_dask_client.return_value.__enter__.return_value = mock_client

    kwargs = dict(
        mslist="mock1.ms,mock2.ms",
        workdir="/mock/workdir",
        prefix="target",
        jobid=42,
        cpu_frac=0.5,
        mem_frac=0.5,
        remote_log=True,
    )

    if raises:
        with pytest.raises(
            RuntimeError,
            match="Solar sidereal motion correction is failed.",
        ):
            run_solar_siderealcor_jobs.fn(**kwargs)
    else:
        result = run_solar_siderealcor_jobs.fn(**kwargs)
        assert result == 0

    mock_makedirs.assert_called_once_with("/mock/workdir/logs", exist_ok=True)
    mock_remove.assert_called_once_with("/mock/workdir/logs/cor_sidereal_target.log")
    mock_log_task_saver.assert_called_once()
    mock_thread.join.assert_called_once_with(timeout=5)

    mock_main.assert_called_once_with(
        "mock1.ms,mock2.ms",
        workdir="/mock/workdir",
        cpu_frac=0.5,
        mem_frac=0.5,
        logfile="/mock/workdir/logs/cor_sidereal_target.log",
        jobid=42,
        start_remote_log=True,
        dask_client=mock_client,
    )


@pytest.mark.parametrize("mock_msg,raises", [(0, False), (1, True)])
@patch("paircars.pipeline.master_flow.get_run_context")
@patch("paircars.pipeline.master_flow.start_log_task_saver")
@patch("paircars.pipeline.master_flow.os.makedirs")
@patch("paircars.pipeline.master_flow.os.path.exists")
@patch("paircars.pipeline.master_flow.os.remove")
@patch("paircars.pipeline.do_selfcal.main")
@patch("paircars.pipeline.master_flow.get_dask_client")
def test_run_selfcal_jobs(
    mock_dask_client,
    mock_main,
    mock_remove,
    mock_exists,
    mock_makedirs,
    mock_log_task_saver,
    mock_get_ctx,
    mock_msg,
    raises,
):
    mock_main.return_value = mock_msg
    mock_exists.return_value = True

    mock_ctx = MagicMock()
    mock_ctx.task_run.id = "abc123"
    mock_ctx.task_run.name = "mock_task"
    mock_get_ctx.return_value = mock_ctx

    mock_thread = MagicMock()
    mock_log_task_saver.return_value = mock_thread

    mock_client = MagicMock()
    mock_dask_client.return_value.__enter__.return_value = mock_client

    kwargs = dict(
        mslist="mock1.ms,mock2.ms",
        workdir="/mock/workdir",
        caldir="/mock/caldir",
        metafits="mock.metafits",
        cal_applied=True,
        start_thresh=5.0,
        stop_thresh=3.0,
        max_iter=100,
        max_DR=100000,
        min_iter=5,
        conv_frac=0.3,
        solint="30s",
        do_apcal=True,
        do_polcal=True,
        solar_selfcal=True,
        keep_backup=False,
        uvrange="",
        minuv=0,
        weight="briggs",
        robust=0.0,
        applymode="calonly",
        min_tol_factor=1.0,
        jobid=42,
        cpu_frac=0.5,
        mem_frac=0.5,
        remote_log=True,
    )

    if raises:
        with pytest.raises(RuntimeError, match="Self-calibration is failed."):
            run_selfcal_jobs.fn(**kwargs)
    else:
        result = run_selfcal_jobs.fn(**kwargs)
        assert result == 0

    mock_makedirs.assert_called_once_with("/mock/workdir/logs", exist_ok=True)
    mock_remove.assert_called_once_with("/mock/workdir/logs/selfcal_target.log")
    mock_log_task_saver.assert_called_once()
    mock_thread.join.assert_called_once_with(timeout=5)

    mock_main.assert_called_once_with(
        "mock1.ms,mock2.ms",
        "mock.metafits",
        "/mock/workdir",
        "/mock/caldir",
        cal_applied=True,
        start_thresh=5.0,
        stop_thresh=3.0,
        max_iter=100.0,
        max_DR=100000.0,
        min_iter=5.0,
        conv_frac=0.3,
        solint="30s",
        uvrange="",
        minuv=0.0,
        weight="briggs",
        robust=0.0,
        applymode="calonly",
        min_tol_factor=1.0,
        do_apcal=True,
        do_polcal=True,
        solar_selfcal=True,
        keep_backup=False,
        cpu_frac=0.5,
        mem_frac=0.5,
        logfile="/mock/workdir/logs/selfcal_target.log",
        jobid=42,
        start_remote_log=True,
        dask_client=mock_client,
    )


@pytest.mark.parametrize("mock_msg,raises", [(0, False), (1, True)])
@patch("paircars.pipeline.master_flow.get_run_context")
@patch("paircars.pipeline.master_flow.start_log_task_saver")
@patch("paircars.pipeline.master_flow.os.makedirs")
@patch("paircars.pipeline.master_flow.os.path.exists")
@patch("paircars.pipeline.master_flow.os.remove")
@patch("paircars.pipeline.do_apply_selfcal.main")
@patch("paircars.pipeline.master_flow.get_dask_client")
def test_run_apply_selfcal_sol(
    mock_dask_client,
    mock_main,
    mock_remove,
    mock_exists,
    mock_makedirs,
    mock_log_task_saver,
    mock_get_ctx,
    mock_msg,
    raises,
):
    mock_main.return_value = mock_msg
    mock_exists.return_value = True

    mock_ctx = MagicMock()
    mock_ctx.task_run.id = "abc123"
    mock_ctx.task_run.name = "mock_task"
    mock_get_ctx.return_value = mock_ctx

    mock_thread = MagicMock()
    mock_log_task_saver.return_value = mock_thread

    mock_client = MagicMock()
    mock_dask_client.return_value.__enter__.return_value = mock_client

    kwargs = dict(
        mslist="mock1.ms,mock2.ms",
        metafits="mock.metafits",
        workdir="/mock/workdir",
        caldir="/mock/caldir",
        overwrite_datacolumn=False,
        applymode="calflag",
        jobid=42,
        cpu_frac=0.5,
        mem_frac=0.5,
        remote_log=True,
    )

    if raises:
        with pytest.raises(
            RuntimeError,
            match="Applying self-calibration solutions is failed.",
        ):
            run_apply_selfcal_sol.fn(**kwargs)
    else:
        result = run_apply_selfcal_sol.fn(**kwargs)
        assert result == 0

    mock_makedirs.assert_called_once_with("/mock/workdir/logs", exist_ok=True)
    mock_remove.assert_called_once_with("/mock/workdir/logs/apply_selfcal.log")
    mock_log_task_saver.assert_called_once()
    mock_thread.join.assert_called_once_with(timeout=5)

    mock_main.assert_called_once_with(
        "mock1.ms,mock2.ms",
        "mock.metafits",
        "/mock/workdir",
        "/mock/caldir",
        applymode="calflag",
        overwrite_datacolumn=False,
        start_remote_log=True,
        cpu_frac=0.5,
        mem_frac=0.5,
        logfile="/mock/workdir/logs/apply_selfcal.log",
        jobid=42,
        dask_client=mock_client,
    )


@pytest.mark.parametrize("mock_msg,raises", [(0, False), (1, True)])
@patch("paircars.pipeline.master_flow.get_run_context")
@patch("paircars.pipeline.master_flow.start_log_task_saver")
@patch("paircars.pipeline.master_flow.os.makedirs")
@patch("paircars.pipeline.master_flow.os.path.exists")
@patch("paircars.pipeline.master_flow.os.remove")
@patch("paircars.pipeline.do_imaging.main")
@patch("paircars.pipeline.master_flow.get_dask_client")
def test_run_imaging_jobs(
    mock_dask_client,
    mock_main,
    mock_remove,
    mock_exists,
    mock_makedirs,
    mock_log_task_saver,
    mock_get_ctx,
    mock_msg,
    raises,
):
    mock_main.return_value = mock_msg
    mock_exists.return_value = True

    mock_ctx = MagicMock()
    mock_ctx.task_run.id = "abc123"
    mock_ctx.task_run.name = "mock_task"
    mock_get_ctx.return_value = mock_ctx

    mock_thread = MagicMock()
    mock_log_task_saver.return_value = mock_thread

    mock_client = MagicMock()
    mock_dask_client.return_value.__enter__.return_value = mock_client

    kwargs = dict(
        mslist="mock1.ms,mock2.ms",
        workdir="/mock/workdir",
        outdir="/mock/outdir",
        freqrange="100~200",
        timerange="2023/01/01/00:00:00~2023/01/01/00:10:00",
        minuv=10,
        weight="briggs",
        robust=0.5,
        pol="IQUV",
        freqres=1.28,
        timeres=10.0,
        threshold=1.0,
        use_multiscale=True,
        use_solar_mask=True,
        cutout_rsun=4.0,
        savemodel=True,
        saveres=True,
        jobid=42,
        cpu_frac=0.5,
        mem_frac=0.5,
        remote_log=True,
    )

    if raises:
        with pytest.raises(RuntimeError, match="Imaging is failed."):
            run_imaging_jobs.fn(**kwargs)
    else:
        result = run_imaging_jobs.fn(**kwargs)
        assert result == 0

    mock_makedirs.assert_called_once_with("/mock/workdir/logs", exist_ok=True)
    mock_remove.assert_called_once_with("/mock/workdir/logs/imaging_target.log")
    mock_log_task_saver.assert_called_once()
    mock_thread.join.assert_called_once_with(timeout=5)

    mock_main.assert_called_once_with(
        "mock1.ms,mock2.ms",
        "/mock/workdir",
        "/mock/outdir",
        freqrange="100~200",
        timerange="2023/01/01/00:00:00~2023/01/01/00:10:00",
        pol="IQUV",
        freqres=1.28,
        timeres=10.0,
        weight="briggs",
        robust=0.5,
        minuv=10.0,
        threshold=1.0,
        cutout_rsun=4.0,
        use_multiscale=True,
        use_solar_mask=True,
        savemodel=True,
        saveres=True,
        start_remote_log=True,
        cpu_frac=0.5,
        mem_frac=0.5,
        jobid=42,
        logfile="/mock/workdir/logs/imaging_target.log",
        dask_client=mock_client,
    )


import pytest
from unittest.mock import patch, MagicMock


@pytest.mark.parametrize("mock_msg,raises", [(0, False), (1, True)])
@patch("paircars.pipeline.master_flow.get_run_context")
@patch("paircars.pipeline.master_flow.start_log_task_saver")
@patch("paircars.pipeline.master_flow.os.makedirs")
@patch("paircars.pipeline.master_flow.os.path.exists")
@patch("paircars.pipeline.master_flow.os.remove")
@patch("paircars.pipeline.mwa_pbcor.main")
@patch("paircars.pipeline.master_flow.get_dask_client")
def test_run_apply_pbcor(
    mock_dask_client,
    mock_main,
    mock_remove,
    mock_exists,
    mock_makedirs,
    mock_log_task_saver,
    mock_get_ctx,
    mock_msg,
    raises,
):
    mock_main.return_value = mock_msg
    mock_exists.return_value = True

    mock_ctx = MagicMock()
    mock_ctx.task_run.id = "abc123"
    mock_ctx.task_run.name = "mock_task"
    mock_get_ctx.return_value = mock_ctx

    mock_thread = MagicMock()
    mock_log_task_saver.return_value = mock_thread

    mock_client = MagicMock()
    mock_dask_client.return_value.__enter__.return_value = mock_client

    kwargs = dict(
        imagedir="/mock/imagedir",
        metafits="mock.metafits",
        workdir="/mock/workdir",
        jobid=42,
        cpu_frac=0.5,
        mem_frac=0.5,
        remote_log=True,
    )

    if raises:
        with pytest.raises(
            RuntimeError,
            match="Primary beam correction is failed.",
        ):
            run_apply_pbcor.fn(**kwargs)
    else:
        result = run_apply_pbcor.fn(**kwargs)
        assert result == 0

    mock_makedirs.assert_called_once_with("/mock/workdir/logs", exist_ok=True)
    mock_remove.assert_called_once_with("/mock/workdir/logs/apply_pbcor.log")
    mock_log_task_saver.assert_called_once()
    mock_thread.join.assert_called_once_with(timeout=5)

    mock_main.assert_called_once_with(
        "/mock/imagedir",
        "mock.metafits",
        workdir="/mock/workdir",
        cpu_frac=0.5,
        mem_frac=0.5,
        logfile="/mock/workdir/logs/apply_pbcor.log",
        jobid=42,
        start_remote_log=True,
        dask_client=mock_client,
    )


@patch("paircars.pipeline.master_flow.get_run_context")
@patch("paircars.pipeline.master_flow.start_log_task_saver")
@patch("paircars.pipeline.master_flow.os.makedirs")
@patch("paircars.pipeline.master_flow.os.path.exists")
@patch("paircars.pipeline.master_flow.os.remove")
@patch("paircars.pipeline.make_mwa_overlay.main")
def test_run_make_overlay(
    mock_main,
    mock_remove,
    mock_exists,
    mock_makedirs,
    mock_log_task_saver,
    mock_get_ctx,
):
    mock_exists.return_value = True

    mock_ctx = MagicMock()
    mock_ctx.task_run.id = "abc123"
    mock_ctx.task_run.name = "mock_task"
    mock_get_ctx.return_value = mock_ctx

    mock_thread = MagicMock()
    mock_log_task_saver.return_value = mock_thread

    kwargs = dict(
        imagedir="/mock/imagedir",
        outdir="/mock/outdir",
        workdir="/mock/workdir",
        jobid=42,
        cpu_frac=0.5,
        mem_frac=0.5,
        remote_log=True,
    )

    result = run_make_overlay.fn(**kwargs)
    assert result == 0

    mock_makedirs.assert_any_call("/mock/workdir/logs", exist_ok=True)
    mock_makedirs.assert_any_call("/mock/outdir", exist_ok=True)
    mock_remove.assert_called_once_with("/mock/workdir/logs/overlay.log")
    mock_log_task_saver.assert_called_once()
    mock_thread.join.assert_called_once_with(timeout=5)

    mock_main.assert_called_once_with(
        "/mock/imagedir",
        "/mock/outdir",
        workdir="/mock/workdir",
        cpu_frac=0.5,
        logfile="/mock/workdir/logs/overlay.log",
        jobid=42,
        start_remote_log=True,
    )
