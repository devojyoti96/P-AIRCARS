import pytest
import os
import traceback
from paircars.utils.imaging import *


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
