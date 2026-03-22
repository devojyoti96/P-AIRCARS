import pytest
import numpy as np
import os
from astropy.io import fits
from unittest.mock import patch, MagicMock
from paircars.utils.image_utils import *


@patch("paircars.utils.image_utils.fits.writeto")
@patch("paircars.utils.image_utils.fits.getheader")
@patch("paircars.utils.image_utils.fits.getdata")
@patch("paircars.utils.image_utils.os.path.exists", return_value=True)
@patch("paircars.utils.image_utils.os.system")
@patch("paircars.utils.image_utils.run_wsclean", return_value=0)
def test_create_circular_mask_success(
    mock_run_wsclean,
    mock_system,
    mock_exists,
    mock_getdata,
    mock_getheader,
    mock_writeto,
    tmp_path,
):
    # Prepare mock FITS data
    imsize = 100
    dummy_data = MagicMock()
    dummy_data.__getitem__.return_value = np.zeros((imsize, imsize), dtype=float)
    mock_getdata.return_value = np.zeros((1, 1, imsize, imsize), dtype=float)
    mock_getheader.return_value = {"SIMPLE": True}

    # Create fake MS path
    ms_path = tmp_path / "dummy.ms"
    ms_path.mkdir()

    # Run the function
    mask_file = create_circular_mask(
        str(ms_path), cellsize=5.0, imsize=imsize, mask_radius=10
    )

    # Assertions
    assert mask_file.endswith("-mask.fits")
    mock_run_wsclean.assert_called_once()
    mock_writeto.assert_called_once()
    assert mock_exists.called


@pytest.mark.parametrize(
    "radius",
    [
        (10),
        (90),
        (5),
        (60),
    ],
)
def test_create_circular_mask_array(radius):
    mock_data = np.empty((200, 200))
    masked_array = create_circular_mask_array(mock_data, radius)
    max_pix = radius**2
    assert np.nansum(masked_array) > max_pix


def test_calc_solar_image_stat(dummy_image):
    maxval, minval, rms, total_val, mean_val, median_val, rms_dyn, minmax_dyn = (
        calc_solar_image_stat(dummy_image, disc_size=18)
    )


def test_calc_dyn_range(dummy_image):
    flux, dr, rms = calc_dyn_range(dummy_image, dummy_image, dummy_image)


def test_generate_tb_map(dummy_image):
    outfile = dummy_image.split(".fits")[0] + "_TB.fits"
    assert generate_tb_map(dummy_image) == outfile
    header = fits.getheader(outfile)
    assert header["BUNIT"] == "K"
    os.system(f"rm -rf {outfile}")


@pytest.mark.parametrize(
    "cutout_size",
    [
        (0.1),
        (1.0),
        (0.2),
    ],
)
def test_cutout_image(dummy_image, cutout_size):
    output_file = dummy_image.split(".fits")[0] + "-cutout.fits"
    assert cutout_image(dummy_image, output_file, x_deg=cutout_size) == output_file
    header = fits.getheader(output_file)
    cdelt = abs(header["CDELT1"])
    npix = header["NAXIS1"]
    imsize = round(cdelt * npix, 1)
    assert imsize == cutout_size
    os.system(f"rm -rf {output_file}")


def test_make_timeavg_image(dummy_image):
    outfile_name = dummy_image.split(".fits")[0] + "-tavg.fits"
    assert (
        make_timeavg_image(
            [dummy_image, dummy_image, dummy_image],
            outfile_name,
            keep_wsclean_images=True,
        )
        == outfile_name
    )
    os.system(f"rm -rf {outfile_name}")


def test_make_freqavg_image(dummy_image):
    outfile_name = dummy_image.split(".fits")[0] + "-favg.fits"
    assert (
        make_freqavg_image(
            [dummy_image, dummy_image, dummy_image],
            outfile_name,
            keep_wsclean_images=True,
        )
        == outfile_name
    )
    os.system(f"rm -rf {outfile_name}")


def test_make_stokes_wsclean_imagecube(dummy_image):
    outfile_name = dummy_image.split(".fits")[0] + "-StokesI.fits"
    result = make_stokes_wsclean_imagecube([dummy_image], outfile_name)
    assert result == outfile_name
    os.system(f"rm -rf {outfile_name}")


@pytest.mark.parametrize(
    "ctype_key",
    ["CTYPE3", "CTYPE4"],
)
@patch("paircars.utils.image_utils.timestamp_to_mjdsec")
@patch("paircars.utils.image_utils.fits.getheader")
def test_filter_images(mock_getheader, mock_mjdsec, ctype_key):
    headers = {
        "img1.fits": {
            ctype_key: "FREQ",
            "CDELT3" if ctype_key == "CTYPE3" else "CDELT4": 2e6,
            "CRVAL3" if ctype_key == "CTYPE3" else "CRVAL4": 100e6,
            "DATE-OBS": "2023-01-01T00:00:00.000",
        },
        "img2.fits": {
            ctype_key: "FREQ",
            "CDELT3" if ctype_key == "CTYPE3" else "CDELT4": 2e6,
            "CRVAL3" if ctype_key == "CTYPE3" else "CRVAL4": 100e6,
            "DATE-OBS": "2023-01-01T00:01:30.000",
        },
        "img3.fits": {
            ctype_key: "FREQ",
            "CDELT3" if ctype_key == "CTYPE3" else "CDELT4": 1e6,  # lower BW
            "CRVAL3" if ctype_key == "CTYPE3" else "CRVAL4": 100e6,
            "DATE-OBS": "2023-01-01T00:02:00.000",
        },
    }

    def header_side_effect(image):
        return headers[image]

    mock_getheader.side_effect = header_side_effect
    times = {
        "2023-01-01T00:00:00": 0,
        "2023-01-01T00:01:30": 90,
        "2023-01-01T00:02:00": 120,
    }
    mock_mjdsec.side_effect = lambda t, date_format=1: times[t]
    imagelist = ["img1.fits", "img2.fits", "img3.fits"]
    result = filter_images(imagelist, min_time_sep=60)
    # img3 should be removed due to smaller BW
    assert "img3.fits" not in result
    # img1 and img2 should remain (90 sec separation)
    assert "img1.fits" in result
    assert "img2.fits" in result
