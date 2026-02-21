import pytest
from unittest.mock import patch, MagicMock, call
from paircars.pipeline.make_mwa_overlay import *


@patch("paircars.pipeline.make_mwa_overlay.clean_shutdown")
@patch("paircars.pipeline.make_mwa_overlay.drop_cache")
@patch("paircars.pipeline.make_mwa_overlay.time.sleep", return_value=None)
@patch("paircars.pipeline.make_mwa_overlay.os.system")
@patch("paircars.pipeline.make_mwa_overlay.os.makedirs")
@patch("paircars.pipeline.make_mwa_overlay.os.path.exists")
@patch("paircars.pipeline.make_mwa_overlay.psutil.cpu_count", return_value=8)
@patch("paircars.pipeline.make_mwa_overlay.np.load", return_value=("job", "pass"))
@patch("paircars.pipeline.make_mwa_overlay.init_logger")
@patch("paircars.pipeline.make_mwa_overlay.glob.glob")
@patch("paircars.pipeline.make_mwa_overlay.make_mwa_overlay")
def test_main_make_mwa_overlay(
    mock_make_overlay,
    mock_glob,
    mock_init_logger,
    mock_np_load,
    mock_cpu_count,
    mock_path_exists,
    mock_makedirs,
    mock_os_system,
    mock_sleep,
    mock_drop_cache,
    mock_shutdown,
):

    # -----------------------------
    # Fake os.path.exists behavior
    # -----------------------------
    def fake_exists(path):
        if "jobname_password.npy" in path:
            return True
        if "log.txt" in path:
            return True
        return True

    mock_path_exists.side_effect = fake_exists

    # =========================================================
    # CASE 1: Successful overlays
    # =========================================================
    mock_glob.return_value = ["img1.fits", "img2.fits"]
    mock_make_overlay.return_value = "overlay.png"

    msg = main(
        imagedir="/tmp/images",
        outdir="/tmp/out",
    )

    assert msg == 0
    assert mock_make_overlay.call_count == 2

    # =========================================================
    # CASE 2: No images in directory
    # =========================================================
    mock_glob.return_value = []

    msg = main(
        imagedir="/tmp/images",
        outdir="/tmp/out",
    )

    assert msg == 1

    # =========================================================
    # CASE 3: Overlay returns None (simulate failure)
    # =========================================================
    mock_glob.return_value = ["img1.fits"]
    mock_make_overlay.return_value = None

    msg = main(
        imagedir="/tmp/images",
        outdir="/tmp/out",
    )

    assert msg == 0  # still considered success because list not empty

    # =========================================================
    # CASE 4: Exception branch
    # =========================================================
    mock_make_overlay.side_effect = Exception("Simulated failure")

    msg = main(
        imagedir="/tmp/images",
        outdir="/tmp/out",
    )

    assert msg == 1

    # =========================================================
    # CASE 5: Remote logging branch
    # =========================================================
    mock_make_overlay.side_effect = None
    mock_make_overlay.return_value = "overlay.png"
    mock_glob.return_value = ["img1.fits"]

    msg = main(
        imagedir="/tmp/images",
        outdir="/tmp/out",
        start_remote_log=True,
        logfile="log.txt",
    )

    assert msg == 0
    mock_init_logger.assert_called()
