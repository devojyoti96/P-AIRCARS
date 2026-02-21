import pytest
from unittest.mock import patch, MagicMock, ANY
from paircars.pipeline.mwa_pbcor import *


@pytest.mark.parametrize(
    "header_dict, expected_output, expect_warning",
    [
        ({"CTYPE3": "FREQ", "CRVAL3": 1.4e9}, 1.4e9, False),
        ({"CTYPE4": "FREQ", "CRVAL4": 1.5e9}, 1.5e9, False),
        ({"CTYPE3": "STOKES"}, None, True),
    ],
)
@patch("paircars.utils.image_utils.fits.getheader")
def test_get_fits_freq(
    mock_getheader, header_dict, expected_output, expect_warning, capsys
):
    mock_hdr = MagicMock()
    mock_hdr.keys.return_value = header_dict.keys()
    mock_hdr.__getitem__.side_effect = header_dict.__getitem__
    mock_getheader.return_value = mock_hdr
    result = get_fits_freq("mock.fits")

    assert result == expected_output

    if expect_warning:
        captured = capsys.readouterr()
        assert "No frequency axis" in captured.out


@pytest.mark.parametrize(
    "pb_exists, restore, verbose, returncode, raise_exception",
    [
        # Normal success
        (True, False, False, 0, False),
        # pbfile missing → --save_pb
        (False, False, False, 0, False),
        # restore branch
        (True, True, False, 0, False),
        # verbose branch
        (True, False, True, 0, False),
        # non-zero return code
        (True, False, False, 1, False),
        # subprocess exception
        (True, False, False, 0, True),
    ],
)
@patch("paircars.pipeline.mwa_pbcor.subprocess.run")
@patch("paircars.pipeline.mwa_pbcor.os.path.exists")
@patch("paircars.pipeline.mwa_pbcor.get_fits_freq")
def test_run_pbcor(
    m_get_freq,
    m_exists,
    m_subprocess,
    pb_exists,
    restore,
    verbose,
    returncode,
    raise_exception,
):
    m_get_freq.return_value = 150.0
    m_exists.return_value = pb_exists

    fake_result = MagicMock()
    fake_result.returncode = returncode
    fake_result.stdout = "mock output"

    if raise_exception:
        m_subprocess.side_effect = Exception("boom")
    else:
        m_subprocess.return_value = fake_result
    result = run_pbcor(
        imagename="image.fits",
        metafits="meta.fits",
        pbdir="pbdir",
        pbcor_dir="outdir",
        restore=restore,
        ncpu=4,
        verbose=verbose,
    )
    if raise_exception:
        assert result == 1
        return
    expected_pbfile = "pbdir/freq_150.0.npy"
    expected_outfile = "outdir/image_pbcor.fits"

    called_cmd = m_subprocess.call_args[0][0]

    assert "run-mwa-singlepbcor" in called_cmd
    assert "--num_threads" in called_cmd
    assert "4" in called_cmd
    assert "--interpolated" in called_cmd
    assert "--pb_jones_file" in called_cmd
    assert expected_pbfile in called_cmd
    assert "image.fits" in called_cmd
    assert "meta.fits" in called_cmd
    assert expected_outfile in called_cmd

    if not pb_exists:
        assert "--save_pb" in called_cmd
    else:
        assert "--save_pb" not in called_cmd

    if restore:
        assert "--restore" in called_cmd
    else:
        assert "--restore" not in called_cmd
    assert result == returncode


@pytest.mark.parametrize(
    "images, run_results, make_TB, make_plots, plot_raises, expect_return",
    [
        # No images
        ([], [], True, True, False, 1),
        # Normal success, no TB
        (["a.fits", "b.fits"], [0, 0], False, True, False, 0),
        # Success with TB + plots
        (["a.fits", "b.fits"], [0, 0], True, True, False, 0),
        # Plot raises exception → junk branch
        (["a.fits"], [0], True, True, True, 0),
        # run_pbcor returns failures
        (["a.fits"], [1], True, True, False, 0),
        # Exception branch
        (["a.fits"], None, True, True, False, 1),
    ],
)
@patch("paircars.pipeline.mwa_pbcor.os.system")
@patch("paircars.pipeline.mwa_pbcor.traceback.print_exc")
@patch("paircars.pipeline.mwa_pbcor.plot_in_hpc")
@patch("paircars.pipeline.mwa_pbcor.save_in_hpc")
@patch("paircars.pipeline.mwa_pbcor.generate_tb_map")
@patch("paircars.pipeline.mwa_pbcor.run_pbcor")
@patch("paircars.pipeline.mwa_pbcor.get_fits_freq")
@patch("paircars.pipeline.mwa_pbcor.os.path.getsize")
@patch("paircars.pipeline.mwa_pbcor.glob.glob")
@patch("paircars.pipeline.mwa_pbcor.os.makedirs")
@patch("paircars.pipeline.mwa_pbcor.psutil.virtual_memory")
@patch("paircars.pipeline.mwa_pbcor.psutil.cpu_count")
@patch("paircars.pipeline.mwa_pbcor.delayed")
def test_pbcor_all_images(
    m_delayed,
    m_cpu,
    m_vm,
    m_makedirs,
    m_glob,
    m_getsize,
    m_getfreq,
    m_run_pbcor,
    m_generate_tb,
    m_save_hpc,
    m_plot_hpc,
    m_print_exc,
    m_system,
    images,
    run_results,
    make_TB,
    make_plots,
    plot_raises,
    expect_return,
):
    m_cpu.return_value = 8

    class FakeMem:
        available = 16 * 1024**3  # 16GB

    m_vm.return_value = FakeMem()
    m_getsize.return_value = 1 * 1024**3  # 1GB per image
    m_glob.return_value = images
    m_getfreq.side_effect = lambda x: 150 if "a" in x else 151
    fake_client = MagicMock()

    if run_results is None:
        m_run_pbcor.side_effect = Exception("boom")
    else:
        # delayed returns function passthrough
        m_delayed.side_effect = lambda f: lambda *a, **k: run_results.pop(0)
        fake_client.compute.side_effect = lambda batch: batch
        fake_client.gather.side_effect = lambda batch: batch
    # Plot exception branch
    if plot_raises:
        m_plot_hpc.side_effect = Exception("plot fail")
    result = pbcor_all_images(
        imagedir="images",
        metafits="meta.fits",
        dask_client=fake_client,
        make_TB=make_TB,
        make_plots=make_plots,
        restore=False,
    )


@pytest.mark.parametrize(
    "imagedir_exists, start_remote_log, inject_client, pb_return, raise_exception",
    [
        # Normal success
        (True, False, False, 0, False),
        # pbcor returns failure
        (True, False, False, 1, False),
        # Directory does not exist
        (False, False, False, 0, False),
        # Injected dask client
        (True, False, True, 0, False),
        # Remote logger branch
        (True, True, False, 0, False),
        # Exception branch
        (True, False, False, 0, True),
    ],
)
@patch("paircars.pipeline.mwa_pbcor.os.system")
@patch("paircars.pipeline.mwa_pbcor.clean_shutdown")
@patch("paircars.pipeline.mwa_pbcor.drop_cache")
@patch("paircars.pipeline.mwa_pbcor.time.sleep")
@patch("paircars.pipeline.mwa_pbcor.traceback.print_exc")
@patch("paircars.pipeline.mwa_pbcor.pbcor_all_images")
@patch("paircars.pipeline.mwa_pbcor.scale_worker_and_wait")
@patch("paircars.pipeline.mwa_pbcor.get_local_dask_cluster")
@patch("paircars.pipeline.mwa_pbcor.init_logger")
@patch("paircars.pipeline.mwa_pbcor.np.load")
@patch("paircars.pipeline.mwa_pbcor.os.path.exists")
@patch("paircars.pipeline.mwa_pbcor.os.makedirs")
@patch("paircars.pipeline.mwa_pbcor.psutil.cpu_count")
def test_main_pbcor(
    m_cpu_count,
    m_makedirs,
    m_exists,
    m_np_load,
    m_init_logger,
    m_get_cluster,
    m_scale,
    m_pbcor_all_images,
    m_print_exc,
    m_sleep,
    m_drop_cache,
    m_clean_shutdown,
    m_system,
    imagedir_exists,
    start_remote_log,
    inject_client,
    pb_return,
    raise_exception,
):
    m_cpu_count.return_value = 8
    m_exists.return_value = imagedir_exists
    m_np_load.return_value = ("job", "pass")

    fake_client = MagicMock()
    fake_cluster = MagicMock()
    m_get_cluster.return_value = (fake_client, fake_cluster, "/daskdir")

    if raise_exception:
        m_pbcor_all_images.side_effect = Exception("boom")
    else:
        m_pbcor_all_images.return_value = pb_return

    dask_client = fake_client if inject_client else None

    result = main(
        imagedir="images",
        metafits="meta.fits",
        workdir="",
        make_TB=True,
        make_plots=True,
        restore=False,
        cpu_frac=0.5,
        mem_frac=0.5,
        logfile="log.txt",
        jobid=42,
        start_remote_log=start_remote_log,
        dask_client=dask_client,
    )

    # Workdir created
    m_makedirs.assert_called()

    # Remote logger
    if start_remote_log and imagedir_exists:
        m_init_logger.assert_called()
    else:
        m_init_logger.assert_not_called()

    # Cluster creation
    if inject_client:
        m_get_cluster.assert_not_called()
    else:
        m_get_cluster.assert_called()
        fake_client.close.assert_called_once()
        fake_cluster.close.assert_called_once()
        m_system.assert_called_with("rm -rf /daskdir")

    # Return logic
    if raise_exception:
        assert result == 1
        m_print_exc.assert_called_once()
    elif not imagedir_exists:
        assert result == 1
    else:
        assert result == pb_return

    # Cleanup always runs
    assert m_drop_cache.call_count == 2  # imagedir + workdir
    m_clean_shutdown.assert_called_once()
    m_sleep.assert_called()


@pytest.mark.parametrize(
    "argv, expect_exit",
    [
        (["prog.py"], True),  # No arguments: help and exit
        (["prog.py", "mockdir", "--no_make_TB"], False),  # Normal run
    ],
)
@patch("paircars.pipeline.mwa_pbcor.main", return_value=0)
@patch("paircars.pipeline.mwa_pbcor.sys.exit")
@patch("paircars.pipeline.mwa_pbcor.argparse.ArgumentParser.print_help")
def test_cli_function(
    mock_print_help,
    mock_exit,
    mock_main,
    argv,
    expect_exit,
):
    with patch("sys.argv", argv):
        from paircars.pipeline import mwa_pbcor

        result = mwa_pbcor.cli()
        assert result == expect_exit
