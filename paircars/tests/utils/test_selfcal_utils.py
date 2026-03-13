import pytest
import numpy as np
import os
from unittest.mock import patch, MagicMock
from paircars.utils.selfcal_utils import *


def test_determine_disk_visibility(dummy_msname):
    chans, timestamps, detected = determine_disk_visibility(dummy_msname)


def test_flag_non_disk(dummy_msname):
    flag_non_disk(dummy_msname)


def test_get_quiet_sun_flux():
    flux = get_quiet_sun_flux(100)
    assert flux > 10**4


def test_make_qs_model(dummy_msname):
    qs_model = make_qs_model(dummy_msname)
    assert os.path.exists(qs_model) is True
    assert qs_model == "quiet_sun.cl"
    os.system(f"rm -rf {qs_model}")


@pytest.mark.parametrize(
    "exists_side_effect, make_model_side_effect, expected_msg, expect_applycal",
    [
        ([False, True], None, 1, True),
        ([False, False], None, 1, False),
        ([False], RuntimeError("Boom"), 1, False),
    ],
)
@patch("paircars.utils.selfcal_utils.suppress_output")
@patch("paircars.utils.selfcal_utils.make_qs_model")
@patch("paircars.utils.selfcal_utils.os.system")
@patch("paircars.utils.selfcal_utils.os.path.exists")
@patch("casatasks.applycal")
@patch("casatasks.gaincal")
@patch("casatasks.ft")
@patch("casatasks.delmod")
@patch("casatasks.flagmanager")
def test_quiet_sun_selfcal(
    m_flagmanager,
    m_delmod,
    m_ft,
    m_gaincal,
    m_applycal,
    m_exists,
    m_system,
    m_make_qs_model,
    m_suppress,
    exists_side_effect,
    make_model_side_effect,
    expected_msg,
    expect_applycal,
    dummy_msname,
):
    """
    Unified test covering:
    - success
    - no gain solutions
    - exception
    """

    logger = MagicMock()

    # Mock suppress_output context manager
    m_suppress.return_value.__enter__.return_value = None
    m_suppress.return_value.__exit__.return_value = None

    # Mock make_qs_model behaviour
    if make_model_side_effect:
        m_make_qs_model.side_effect = make_model_side_effect
    else:
        m_make_qs_model.return_value = "fake_qs.cl"

    # Mock os.path.exists behaviour
    m_exists.side_effect = exists_side_effect

    msg, caltable = quiet_sun_selfcal(
        msname=dummy_msname,
        logger=logger,
        selfcaldir="/tmp",
        refant="1",
        solint="60s",
    )

    # -----------------------------------
    # Assertions
    # -----------------------------------
    assert msg == expected_msg


def test_check_valid_image(dummy_image):
    valid_image = check_valid_image(dummy_image)
    assert valid_image is True


def test_calc_leakage(dummy_image):
    q, u, v, _, _, _ = calc_leakage(dummy_image)
    assert q <= 1
    assert u <= 1
    assert v <= 1


def test_correct_image_leakage(dummy_image):
    cor_image, cor_model = correct_image_leakage(dummy_image)
    assert os.path.exists(cor_image) is True
    assert cor_model is None
    os.system(f"rm -rf {cor_image}")


@pytest.mark.parametrize(
    "pbcor, leakagecor, pbuncor",
    [
        (False, False, False),  # no corrections
        (True, False, False),  # pb only
        (True, True, False),  # pb + leakage
        (True, True, True),  # full pipeline
    ],
)
@patch("paircars.utils.selfcal_utils.subprocess.run")
@patch("paircars.utils.selfcal_utils.os.path.exists")
@patch("paircars.utils.selfcal_utils.correct_image_leakage")
@patch("paircars.utils.selfcal_utils.calc_leakage")
@patch("paircars.utils.selfcal_utils.fits.getheader")
def test_correct_pbcor_leakage(
    m_getheader,
    m_calc_leakage,
    m_correct_leakage,
    m_exists,
    m_run,
    pbcor,
    leakagecor,
    pbuncor,
):
    """
    Test correct_pbcor_leakage with all combinations of flags.
    """
    m_getheader.return_value = {"CRVAL3": 150e6}
    m_exists.return_value = True
    m_run.return_value = MagicMock(returncode=0)
    m_calc_leakage.return_value = (0.1, 0.2, 0.3, 0.01, 0.02, 0.03)  # q,u,v  # errors
    m_correct_leakage.return_value = (
        "image_leakcor.fits",
        "model_leakcor.fits",
    )
    out_img, out_model, leak_info = correct_pbcor_leakage(
        imagename="test.fits",
        modelname="model.fits",
        metafits="meta.fits",
        pbcor=pbcor,
        leakagecor=leakagecor,
        pbuncor=pbuncor,
        ncpu=2,
    )
    assert isinstance(out_img, str)
    assert isinstance(out_model, str)
    assert isinstance(leak_info, list)
    if leakagecor:
        m_calc_leakage.assert_called_once()
        m_correct_leakage.assert_called_once()
        assert len(leak_info) == 6
    else:
        m_calc_leakage.assert_not_called()
        m_correct_leakage.assert_not_called()
        assert leak_info == []
    if pbcor or pbuncor:
        assert m_run.called


@pytest.mark.parametrize("valid_image", [True, False])
@patch("paircars.utils.selfcal_utils.fits.writeto")
@patch("paircars.utils.selfcal_utils.fits.getdata")
@patch("paircars.utils.selfcal_utils.fits.getheader")
@patch("paircars.utils.selfcal_utils.correct_pbcor_leakage")
@patch("paircars.utils.selfcal_utils.check_valid_image")
def test_single_image_update_leakage(
    m_check_valid,
    m_correct_pb,
    m_getheader,
    m_getdata,
    m_writeto,
    valid_image,
):
    """
    Test single_image_update_leakage for valid and invalid image cases.
    """
    m_check_valid.return_value = valid_image
    # fake cube data (4 Stokes)
    cube_data = np.ones((4, 1, 10, 10))
    wsclean_data = np.zeros((1, 1, 10, 10))

    m_correct_pb.return_value = (
        "corrected_image.fits",
        "corrected_model.fits",
        [0.1, 0.2, 0.3, 0.01, 0.02, 0.03],
    )

    m_getheader.return_value = {"TEST": "HDR"}

    # FITS behavior
    def getdata_side_effect(name):
        if name in ["corrected_image.fits", "corrected_model.fits"]:
            return cube_data
        return wsclean_data.copy()

    m_getdata.side_effect = getdata_side_effect

    wsclean_images = ["imgI.fits", "imgQ.fits", "imgU.fits", "imgV.fits"]
    wsclean_models = ["modI.fits", "modQ.fits", "modU.fits", "modV.fits"]

    result = single_image_update_leakage(
        wsclean_images=wsclean_images,
        wsclean_models=wsclean_models,
        image_cube="cube.fits",
        model_cube="model_cube.fits",
        metafits="meta.fits",
    )

    if valid_image:
        m_correct_pb.assert_called_once()
        assert m_writeto.call_count == 8  # 4 images + 4 models
        assert isinstance(result, list)
        assert len(result) == 6
    else:
        m_correct_pb.assert_not_called()
        m_writeto.assert_not_called()
        assert result is None


@pytest.mark.parametrize(
    "valid_image, use_logger, leakage_return",
    [
        (True, True, [1, 2, 3]),  # normal case
        (True, False, [1, 2, 3]),  # print branch
        (True, True, None),  # no leakage appended
        (False, True, [1, 2, 3]),  # invalid image skipped
    ],
)
@patch("paircars.utils.selfcal_utils.os.system")
@patch("paircars.utils.selfcal_utils.single_image_update_leakage")
@patch("paircars.utils.selfcal_utils.check_valid_image")
def test_correct_spectrosnap_pbleak(
    m_check_valid,
    m_single_update,
    m_system,
    valid_image,
    use_logger,
    leakage_return,
):
    """
    Test correct_spectrosnap_pbleak covering all branches.
    """
    m_check_valid.return_value = valid_image
    m_single_update.return_value = leakage_return
    logger = MagicMock() if use_logger else None
    # include MFS and non-MFS image
    image_dic = {
        "test_image.fits": ["imgI", "imgQ"],
    }

    model_dic = {
        "test_image.fits": ["modI", "modQ"],
    }
    result = correct_spectrosnap_pbleak(
        image_dic=image_dic,
        model_dic=model_dic,
        metafits="meta.fits",
        logger=logger,
        pbcor=True,
        leakagecor=True,
        pbuncor=True,
        ncpu=2,
    )
    # MFS image should always be skipped → only one call possible
    if valid_image:
        m_single_update.assert_called_once()

        if leakage_return is not None:
            assert result == [leakage_return]
        else:
            assert result == []
    else:
        m_single_update.assert_not_called()
        assert result == []

    # logger vs print branch
    if valid_image and use_logger:
        logger.info.assert_called_once()

    # cleanup command always executed
    m_system.assert_called_once()
    assert "rm -rf" in m_system.call_args[0][0]
