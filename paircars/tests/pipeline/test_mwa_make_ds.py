import pytest, shutil
from unittest.mock import patch, MagicMock, call, ANY
from paircars.pipeline.mwa_make_ds import *


@pytest.mark.parametrize(
    "cpu_frac, mem_frac, raise_exception",
    [
        (0.5, 0.5, False),  # Normal case
        (0.9, 0.9, False),  # Clamp cpu/mem > 0.8
        (0.5, 0.5, True),  # Exception path
    ],
)
@patch("paircars.pipeline.mwa_make_ds.drop_cache")
@patch("paircars.pipeline.mwa_make_ds.time.sleep")
@patch("paircars.pipeline.mwa_make_ds.traceback.print_exc")
@patch("paircars.pipeline.mwa_make_ds.os.system")
@patch("paircars.pipeline.mwa_make_ds.glob.glob")
@patch("paircars.pipeline.mwa_make_ds.make_ds_plot")
@patch("paircars.pipeline.mwa_make_ds.get_MWA_OBSID")
@patch("paircars.pipeline.mwa_make_ds.get_ms_size")
@patch("paircars.pipeline.mwa_make_ds.psutil.virtual_memory")
@patch("paircars.pipeline.mwa_make_ds.psutil.cpu_count")
@patch("paircars.pipeline.mwa_make_ds.delayed")
def test_make_solar_DS(
    m_delayed,
    m_cpu_count,
    m_virtual_mem,
    m_get_ms_size,
    m_get_obsid,
    m_make_ds_plot,
    m_glob,
    m_system,
    m_print_exc,
    m_sleep,
    m_drop_cache,
    cpu_frac,
    mem_frac,
    raise_exception,
):
    mslist = ["a.ms", "b.ms"]
    m_cpu_count.return_value = 8

    class FakeMem:
        available = 8 * 1024**3  # 8 GB

    m_virtual_mem.return_value = FakeMem()

    m_get_ms_size.return_value = 1  # 1 GB per MS
    m_get_obsid.return_value = "123456"
    fake_client = MagicMock()

    if raise_exception:
        fake_client.compute.side_effect = Exception("Dask failure")
    else:
        fake_client.compute.return_value = ["future"]
        fake_client.gather.return_value = [("file1.ds", None)]

    m_delayed.side_effect = lambda f: lambda *args, **kwargs: ("file1.ds", None)

    m_make_ds_plot.return_value = "final_plot.png"
    m_glob.return_value = ["dummy.nc"]
    result = make_solar_DS(
        mslist=mslist,
        dask_client=fake_client,
        metafits="meta.fits",
        workdir="workdir",
        outdir="outdir",
        cpu_frac=cpu_frac,
        mem_frac=mem_frac,
    )

    if raise_exception:
        assert result is None
        m_print_exc.assert_called_once()
    else:
        assert result == "final_plot.png"

        # Plot called
        m_make_ds_plot.assert_called_once()

        # GOES cleanup
        m_system.assert_called_with("rm -rf dummy.nc")

    # drop_cache always called for each MS
    assert m_drop_cache.call_count == len(mslist) + 1
    shutil.rmtree("outdir", ignore_errors=True)

    # sleep always called in finally
    m_sleep.assert_called_once_with(5)


@pytest.mark.parametrize(
    "mslist, start_remote_log, inject_client, ds_returns, raise_exception",
    [
        # Normal success
        ("a.ms,b.ms", False, False, "plot.png", False),
        # DS returns None → failure
        ("a.ms,b.ms", False, False, None, False),
        # Empty mslist
        ("", False, False, "plot.png", False),
        # Injected dask client (no cluster creation)
        ("a.ms", False, True, "plot.png", False),
        # Remote logger branch
        ("a.ms", True, False, "plot.png", False),
        # Exception branch
        ("a.ms", False, False, "plot.png", True),
    ],
)
@patch("paircars.pipeline.mwa_make_ds.os.system")
@patch("paircars.pipeline.mwa_make_ds.clean_shutdown")
@patch("paircars.pipeline.mwa_make_ds.drop_cache")
@patch("paircars.pipeline.mwa_make_ds.time.sleep")
@patch("paircars.pipeline.mwa_make_ds.traceback.print_exc")
@patch("paircars.pipeline.mwa_make_ds.make_solar_DS")
@patch("paircars.pipeline.mwa_make_ds.scale_worker_and_wait")
@patch("paircars.pipeline.mwa_make_ds.get_local_dask_cluster")
@patch("paircars.pipeline.mwa_make_ds.init_logger")
@patch("paircars.pipeline.mwa_make_ds.np.load")
@patch("paircars.pipeline.mwa_make_ds.os.path.exists")
@patch("paircars.pipeline.mwa_make_ds.os.makedirs")
@patch("paircars.pipeline.mwa_make_ds.psutil.cpu_count")
def test_main_ds(
    m_cpu_count,
    m_makedirs,
    m_exists,
    m_np_load,
    m_init_logger,
    m_get_cluster,
    m_scale,
    m_make_solar_DS,
    m_print_exc,
    m_sleep,
    m_drop_cache,
    m_clean_shutdown,
    m_system,
    mslist,
    start_remote_log,
    inject_client,
    ds_returns,
    raise_exception,
):
    m_cpu_count.return_value = 8
    m_exists.return_value = True
    m_np_load.return_value = ("jobname", "password")

    fake_client = MagicMock()
    fake_cluster = MagicMock()
    m_get_cluster.return_value = (fake_client, fake_cluster, "/daskdir")

    if raise_exception:
        m_make_solar_DS.side_effect = Exception("failure")
    else:
        m_make_solar_DS.return_value = ds_returns

    dask_client = fake_client if inject_client else None

    msg = main(
        mslist=mslist,
        metafits="meta.fits",
        workdir="",
        outdir="",
        logfile="test.log",
        jobid="42",
        start_remote_log=start_remote_log,
        dask_client=dask_client,
    )

    # Workdir created
    m_makedirs.assert_called()

    # Remote logger
    if start_remote_log:
        m_init_logger.assert_called()
    else:
        m_init_logger.assert_not_called()

    # Cluster creation
    if inject_client:
        m_get_cluster.assert_not_called()
    else:
        m_get_cluster.assert_called()

    # Return logic
    if raise_exception:
        assert msg == 1
        m_print_exc.assert_called_once()
    elif ds_returns is None:
        assert msg == 1
    else:
        assert msg == 0

    # clean shutdown always called
    m_clean_shutdown.assert_called_once()

    # Cluster close + cleanup only if auto-created
    if not inject_client:
        fake_client.close.assert_called_once()
        fake_cluster.close.assert_called_once()
        m_system.assert_called_with("rm -rf /daskdir")

    # sleep in finally
    m_sleep.assert_called()
