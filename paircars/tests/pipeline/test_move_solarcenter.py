import pytest
from unittest.mock import patch, MagicMock, call, ANY
from paircars.pipeline.move_solarcenter import *


@pytest.mark.parametrize(
    "start_remote_log, logfile, dask_client_injected, raise_exception",
    [
        # Normal run, no remote log, create local cluster
        (False, None, False, False),
        # Remote log enabled and logfile exists
        (True, "test.log", False, False),
        # Injected dask client (skip cluster creation)
        (False, None, True, False),
        # Exception during compute
        (False, None, False, True),
    ],
)
@patch("paircars.pipeline.move_solarcenter.clean_shutdown")
@patch("paircars.pipeline.move_solarcenter.drop_cache")
@patch("paircars.pipeline.move_solarcenter.scale_worker_and_wait")
@patch("paircars.pipeline.move_solarcenter.get_local_dask_cluster")
@patch("paircars.pipeline.move_solarcenter.init_logger")
@patch("paircars.pipeline.move_solarcenter.np.load")
@patch("paircars.pipeline.move_solarcenter.os.path.exists")
@patch("paircars.pipeline.move_solarcenter.os.makedirs")
@patch("paircars.pipeline.move_solarcenter.get_MWA_OBSID")
@patch("paircars.pipeline.move_solarcenter.psutil.cpu_count")
@patch("paircars.pipeline.move_solarcenter.delayed")
def test_main(
    m_delayed,
    m_cpu_count,
    m_get_obsid,
    m_makedirs,
    m_exists,
    m_np_load,
    m_init_logger,
    m_get_cluster,
    m_scale,
    m_drop_cache,
    m_clean_shutdown,
    start_remote_log,
    logfile,
    dask_client_injected,
    raise_exception,
):
    # -------------------------
    # Basic mocks
    # -------------------------
    m_get_obsid.return_value = "12345"
    m_cpu_count.return_value = 8
    m_exists.return_value = True
    m_np_load.return_value = ("jobname", "password")
    fake_client = MagicMock()
    fake_cluster = MagicMock()
    if raise_exception:
        fake_client.compute.side_effect = Exception("Dask failure")
    else:
        fake_client.compute.return_value = ["future"]
        fake_client.gather.return_value = ["done"]
    m_get_cluster.return_value = (fake_client, fake_cluster, "/daskdir")
    # delayed(move_to_sun)(ms) → just return ms
    m_delayed.side_effect = lambda f: lambda x: x
    # Inject client optionally
    dask_client = fake_client if dask_client_injected else None
    # -------------------------
    # Run
    # -------------------------
    msg = main(
        mslist="a.ms,b.ms",
        workdir="",
        cpu_frac=0.5,
        mem_frac=0.5,
        logfile=logfile,
        jobid=42,
        start_remote_log=start_remote_log,
        dask_client=dask_client,
    )
    # Workdir created
    m_makedirs.assert_called()
    # Cluster created only if not injected
    if not dask_client_injected:
        m_get_cluster.assert_called_once()
        m_scale.assert_called()
    else:
        m_get_cluster.assert_not_called()
    # Remote logger only if enabled and logfile provided
    if start_remote_log and logfile:
        m_init_logger.assert_called()
    else:
        m_init_logger.assert_not_called()
    # Cleanup always happens
    assert m_drop_cache.call_count >= 2
    m_clean_shutdown.assert_called_once()
    # Return code
    if raise_exception:
        assert msg == 1
    else:
        assert msg == 0
