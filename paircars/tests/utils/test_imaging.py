import pytest
from paircars.utils.imaging import *


@pytest.mark.parametrize(
    "n,expected",
    [
        (1, True),
        (2, True),
        (3, True),
        (5, True),
        (7, True),
        (30, True),  # 2*3*5
        (49, True),  # 7^2
        (11, False),
        (13, False),
        (121, False),  # 11^2
    ],
)
def test_is_fft_good(n, expected):
    assert is_fft_good(n) == expected


@pytest.mark.parametrize(
    "n,expected",
    [
        (1, 2),  # edge case
        (2, 2),
        (3, 4),  # next power-of-two (<128 rule)
        (10, 16),
        (31, 32),
        (64, 64),
        (127, 128),
        (128, 128),
        (129, 140),  # next FFT-good number using 2,3,5,7 factors
    ],
)
def test_get_fft_size(n, expected):
    result = get_fft_size(n)
    assert result >= n
    assert result % 2 == 0


def test_calc_sun_dia():
    assert calc_sun_dia(1000.0) == 34.2


def test_calc_maxuv(dummy_msname):
    maxuv, maxuv_l = calc_maxuv(dummy_msname)
    assert maxuv == 2610.63
    assert maxuv_l == 1154.31


def test_calc_minuv(dummy_msname):
    minuv, minuv_l = calc_minuv(dummy_msname)
    assert minuv == 5.36
    assert minuv_l == 2.37


def test_calc_field_of_view(dummy_msname):
    assert calc_field_of_view(dummy_msname, FWHM=True) == 142281.87
    assert calc_field_of_view(dummy_msname, FWHM=False) == 237913.94


def test_get_optimal_image_interval(dummy_msname):
    ntime, nchan = get_optimal_image_interval(
        dummy_msname,
        temporal_tol_factor=0.1,
        spectral_tol_factor=0.1,
    )
    assert ntime == 1
    assert nchan == 1
    ntime, nchan = get_optimal_image_interval(
        dummy_msname,
        temporal_tol_factor=1.0,
        spectral_tol_factor=0.001,
    )
    assert ntime == 1
    assert nchan == 1


def test_calc_psf(dummy_msname):
    assert calc_psf(dummy_msname) == 214.43


def test_calc_npix_in_psf():
    assert calc_npix_in_psf("natural") == 5.0
    assert calc_npix_in_psf("uniform") == 3.0
    assert calc_npix_in_psf("briggs", robust=0.0) == 4.0


def test_calc_cellsize(dummy_msname):
    assert calc_cellsize(dummy_msname, 3) == 71.5
    assert calc_cellsize(dummy_msname, 5) == 42.9


def test_calc_multiscale_scales(dummy_msname):
    scales = calc_multiscale_scales(dummy_msname, 3, max_scale=16)
    assert scales == [0, 3, 6, 12, 13]
    scales = calc_multiscale_scales(dummy_msname, 3, max_scale=8)
    assert scales == [0, 3, 6]


def test_get_multiscale_bias():
    assert get_multiscale_bias(100) == 0.6
    assert get_multiscale_bias(200) == 0.9
    assert get_multiscale_bias(150) == 0.775
