import pytest
from unittest.mock import patch, MagicMock
from paircars.pipeline.single_image_mwapbcor import *


@pytest.mark.parametrize(
    "restore, stokes_axis, load_pb, shape_match, save_pb, expect_return",
    [
        # Normal correction, full pol axis=3
        (False, 3, False, True, False, True),
        # Restore branch
        (True, 3, False, True, False, True),
        # Load PB file (shape matches)
        (False, 3, True, True, False, True),
        # Load PB file (shape mismatch → re-estimate)
        (False, 3, True, False, False, True),
        # Save PB branch
        (False, 3, False, True, True, True),
        # Stokes axis=4
        (False, 4, False, True, False, True),
        # Stokes I only
        (False, 1, False, True, False, True),
        # Exception branch
        (False, 3, False, True, False, False),
    ],
)
@patch("paircars.pipeline.single_image_mwapbcor.traceback.print_exc")
@patch("paircars.pipeline.single_image_mwapbcor.calc_leakage")
@patch("paircars.pipeline.single_image_mwapbcor.fits.open")
@patch("paircars.pipeline.single_image_mwapbcor.fits.writeto")
@patch("paircars.pipeline.single_image_mwapbcor.inv", side_effect=lambda x: x)
@patch("paircars.pipeline.single_image_mwapbcor.B2IQUV")
@patch("paircars.pipeline.single_image_mwapbcor.get_inst_pols")
@patch("paircars.pipeline.single_image_mwapbcor.get_IQUV")
@patch("paircars.pipeline.single_image_mwapbcor.get_jones_array")
@patch("paircars.pipeline.single_image_mwapbcor.get_azza_from_fits")
@patch("paircars.pipeline.single_image_mwapbcor.np.save")
@patch("paircars.pipeline.single_image_mwapbcor.np.load")
@patch("paircars.pipeline.single_image_mwapbcor.glob.glob")
@patch("paircars.pipeline.single_image_mwapbcor.fits.getdata")
@patch("paircars.pipeline.single_image_mwapbcor.fits.getheader")
@patch("paircars.pipeline.single_image_mwapbcor.os.system")
@patch("paircars.pipeline.single_image_mwapbcor.os.path.exists")
def test_get_pbcor_image(
    m_exists,
    m_system,
    m_getheader,
    m_getdata,
    m_glob,
    m_np_load,
    m_np_save,
    m_get_azza,
    m_get_jones,
    m_get_IQUV,
    m_get_inst_pols,
    m_B2IQUV,
    m_inv,
    m_writeto,
    m_fits_open,
    m_calc_leakage,
    m_print_exc,
    restore,
    stokes_axis,
    load_pb,
    shape_match,
    save_pb,
    expect_return,
):
    # ----------------------------
    # FITS header setup
    # ----------------------------
    header = {
        "CTYPE3": "FREQ" if stokes_axis != 3 else "STOKES",
        "CTYPE4": "STOKES" if stokes_axis == 4 else "FREQ",
        "CRVAL3": 150e6,
        "CRVAL4": 150e6,
        "CDELT3": 1000.0,
        "CDELT4": 1000.0,
        "GRIDNUM": 1,
    }
    m_getheader.return_value = header
    m_getdata.return_value = np.zeros((1, 4, 2, 2))

    # Beam files
    m_glob.return_value = ["mwa_full_embedded_element_pattern_150.0.h5"]

    # Sweet spots
    m_np_load.return_value = {1: [None, None, None, 5]}

    # Az/ZA
    fake_azza = {
        "za_rad": np.zeros((2, 2)),
        "astro_az_rad": np.zeros((2, 2)),
    }
    m_get_azza.return_value = fake_azza

    # Jones
    fake_jones = np.ones((4, 2, 2))
    m_get_jones.return_value = fake_jones

    # PB load branch
    if load_pb:
        pb_array = np.array(
            [False, fake_jones if shape_match else np.ones((1, 2, 2))],
            dtype="object",
        )
        m_np_load.return_value = pb_array
        m_exists.side_effect = lambda x: True
    else:
        m_exists.side_effect = lambda x: False

    # IQUV + inst pol
    fake_stokes = {"I": np.ones((2, 2))}
    m_get_IQUV.return_value = fake_stokes
    m_get_inst_pols.return_value = np.ones((2, 2, 2, 2))

    m_B2IQUV.return_value = {
        "I": np.ones((2, 2)),
        "Q": np.ones((2, 2)),
        "U": np.ones((2, 2)),
        "V": np.ones((2, 2)),
    }

    # Leakage
    m_calc_leakage.return_value = (0.1, 0.1, 0.1, 0, 0, 0)

    # FITS open mock
    fake_hdul = MagicMock()
    m_fits_open.return_value.__enter__.return_value = fake_hdul

    if not expect_return:
        m_get_jones.side_effect = Exception("boom")

    result = get_pbcor_image(
        imagename="test.fits",
        outfile="out.fits",
        metafits="meta.fits",
        pb_jones_file="pb.npy" if load_pb else "",
        save_pb=save_pb,
        restore=restore,
        nthreads=2,
    )
