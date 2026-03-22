import traceback
import os
import numpy as np
import pytest
from unittest.mock import patch, MagicMock
from casatasks import casalog
from casatools import table
from paircars.utils.calibration import *

try:
    casalogfile = casalog.logfile()
    os.system("rm -rf " + casalogfile)
except BaseException:
    traceback.print_exc()
    pass


def test_fill_nan_gains():
    x = np.array([1, 2, 3, 4, 5])
    data = np.array([1.0, np.nan, 3.0, np.nan, 5.0])
    result = fill_nan_gains(x, data)
    expected = np.array([1, 2, 3, 4, 5], dtype=float)
    np.testing.assert_allclose(result, expected)


def test_fluxcal_caltable(dummy_caltable):
    scaled_caltable = dummy_caltable.split(".bcal")[0] + "_scaled.bcal"
    os.system(f"cp -r {dummy_caltable} {scaled_caltable}")
    scaled_caltable = fluxcal_caltable(scaled_caltable)
    tb = table()
    tb.open(scaled_caltable)
    scaled_gain = tb.getcol("CPARAM")
    tb.close()
    tb.open(dummy_caltable)
    gain = tb.getcol("CPARAM")
    tb.close()
    os.system(f"rm -rf {scaled_caltable}")
    scale_factor = np.abs(np.nanmedian(scaled_gain / gain))
    assert scale_factor < 1.0


def test_merge_caltables(dummy_caltables, tmp_path):
    merged = tmp_path / "merged.K"
    result = merge_caltables(
        dummy_caltables.copy(), str(merged), append=False, keepcopy=True
    )
    assert os.path.exists(result)
    tb = table()
    tb.open(result)
    merged_rows = tb.nrows()
    tb.close()
    tb.open(dummy_caltables[0])
    single_rows = tb.nrows()
    tb.close()
    assert merged_rows == 2 * single_rows


@pytest.mark.parametrize("overwrite", [True, False])
@patch("paircars.utils.calibration.os.system")
@patch("paircars.utils.calibration.os.path.exists")
@patch("paircars.utils.calibration.fill_nan_gains")
@patch("paircars.utils.calibration.table")
def test_interpolate_bpass(
    mock_table,
    mock_fill_nan,
    mock_exists,
    mock_system,
    overwrite,
):
    tb = MagicMock()
    mock_table.return_value = tb
    freqs = np.array([[100.0, 110.0, 120.0]])
    gains = np.ones((2, 3, 2), dtype=complex)
    flags = np.zeros((2, 3, 2), dtype=bool)

    def getcol_side_effect(name):
        if name == "CHAN_FREQ":
            return freqs
        if name == "CPARAM":
            return gains.copy()
        if name == "FLAG":
            return flags.copy()

    tb.getcol.side_effect = getcol_side_effect
    mock_fill_nan.side_effect = lambda x, y: y
    mock_exists.return_value = False
    caltables = ["cal1", "cal2"]
    result = interpolate_bpass(caltables, overwrite=overwrite)
    assert len(result) == 2
    if overwrite:
        assert result == ["cal1", "cal2"]
    else:
        assert result == ["cal1.interp", "cal2.interp"]
    assert tb.open.called
    assert tb.getcol.called
    assert tb.putcol.called


@pytest.mark.parametrize("overwrite", [True, False])
@patch("paircars.utils.calibration.os.system")
@patch("paircars.utils.calibration.dask.compute")
@patch("paircars.utils.calibration.xds_to_zarr")
@patch("paircars.utils.calibration.xds_from_zarr")
@patch("paircars.utils.calibration.get_quartical_soltype")
@patch("paircars.utils.calibration.fill_nan_gains")
def test_interpolate_quartical(
    mock_fill,
    mock_soltype,
    mock_xds_from,
    mock_xds_to,
    mock_compute,
    mock_system,
    overwrite,
):
    mock_soltype.return_value = ["G"]
    ntime = 2
    nchan = 3
    nant = 2
    ndir = 1
    npol = 4
    freqs = np.array([100.0, 110.0, 120.0])
    gain_data = np.ones((ntime, nchan, nant, ndir, npol), dtype=complex)
    gain_flags = np.zeros((ntime, nchan, nant, ndir), dtype=bool)
    ds = MagicMock()
    ds.gain_freq.to_numpy.return_value = freqs
    ds.gains.to_numpy.return_value = gain_data.copy()
    ds.gain_flags.to_numpy.return_value = gain_flags.copy()
    ds.gain_flags.values = gain_flags.copy()
    ds.update = MagicMock()
    mock_xds_from.return_value = [ds]
    mock_fill.side_effect = lambda x, y: y
    caltables = ["cal1", "cal2"]
    result = interpolate_quartical(caltables, overwrite=overwrite)
    # ---- Assertions ----
    assert len(result) == 2
    if overwrite:
        assert result == ["cal1", "cal2"]
    else:
        assert result == ["cal1.interp", "cal2.interp"]
    assert mock_xds_from.called
    assert mock_xds_to.called
    assert mock_compute.called


@patch("paircars.utils.calibration.table")
def test_get_cal_flag_info(mock_table):
    tb = MagicMock()
    mock_table.return_value = tb
    npol = 2
    nchan = 3
    ntime = 2
    nant = 2
    # CASA stores flags as (npol, nchan, nrows)
    nrows = ntime * nant
    flags = np.zeros((npol, nchan, nrows), dtype=bool)
    # Fully flag channel 1
    flags[:, 1, :] = True
    # Fully flag antenna 0
    flags[:, :, 0::2] = True
    # Fully flag time index 1
    flags[:, :, 2:] = True
    times = np.array([1, 1, 2, 2])

    def getcol_side_effect(col):
        if col == "FLAG":
            return flags
        elif col == "TIME":
            return times

    tb.getcol.side_effect = getcol_side_effect
    result = get_cal_flag_info("fake_caltable")
    flag_chans, flag_ants, flag_times, flag_frac, chan_frac, ant_frac, time_frac = (
        result
    )
    # ---- Assertions ----
    assert isinstance(flag_chans, list)
    assert isinstance(flag_ants, list)
    assert isinstance(flag_times, list)
    assert 1 in flag_chans
    assert flag_frac >= 0
    assert 0 <= chan_frac <= 1
    assert 0 <= ant_frac <= 1
    assert 0 <= time_frac <= 1

    assert tb.open.called
    assert tb.close.called


def test_get_psf_size(dummy_msname):
    assert get_psf_size(dummy_msname) == 214.43


def test_calc_bw_smearing_freqwidth(dummy_msname):
    assert calc_bw_smearing_freqwidth(dummy_msname, full_FoV=False, FWHM=False) == 3.04
    assert calc_bw_smearing_freqwidth(dummy_msname, full_FoV=False, FWHM=True) == 3.04
    assert calc_bw_smearing_freqwidth(dummy_msname, full_FoV=True, FWHM=True) == 0.16
    assert calc_bw_smearing_freqwidth(dummy_msname, full_FoV=True, FWHM=False) == 0.16


def test_calc_time_smearing_timewidth(dummy_msname):
    assert (
        calc_time_smearing_timewidth(dummy_msname, full_FoV=False, FWHM=False) == 84.0
    )
    assert calc_time_smearing_timewidth(dummy_msname, full_FoV=False, FWHM=True) == 84.0
    assert calc_time_smearing_timewidth(dummy_msname, full_FoV=True, FWHM=True) == 4.0
    assert calc_time_smearing_timewidth(dummy_msname, full_FoV=True, FWHM=False) == 2.0


def test_max_time_solar_smearing(dummy_msname):
    assert max_time_solar_smearing(dummy_msname) == 2573.16


def test_get_caltable_metadata(dummy_caltable):
    metadata = get_caltable_metadata(dummy_caltable)
    assert metadata["JonesType"] == "B Jones"


def test_get_nearest_bandpass_table(dummy_caltable):
    nearest_caltable = get_nearest_bandpass_table([dummy_caltable], 132)
    assert nearest_caltable == dummy_caltable


def test_get_nearest_gaincal_table(dummy_caltable):
    nearest_caltable = get_nearest_gaincal_table(
        [dummy_caltable], "2015/01/01/00:00:00"
    )
    assert nearest_caltable == dummy_caltable


def test_get_gleam_uvrange(dummy_msname):
    uvrange = get_gleam_uvrange(dummy_msname)
    assert uvrange == "49.5~1104.6lambda"


def test_uvrange_casa_to_quartical(dummy_msname):
    minuv, maxuv = uvrange_casa_to_quartical(dummy_msname)
    assert minuv == 0.0
    assert maxuv == 0.0
    minuv, maxuv = uvrange_casa_to_quartical(dummy_msname, uvrange="1~100lambda")
    assert minuv == 2.3
    assert maxuv == 226.3


def test_solint_in_float():
    assert solint_in_float("20s") == 20.0
    assert solint_in_float("1min") == 60.0


def test_quartical_matrix_normalize(dummy_quartical_table):
    norm_table = quartical_matrix_normalize(dummy_quartical_table, overwrite=False)
    assert norm_table == f"{dummy_quartical_table}.poldist"
    os.system(f"rm -rf {dummy_quartical_table}.poldist")


def test_get_quartical_table_metadata(dummy_quartical_table):
    metadata = get_quartical_table_metadata(dummy_quartical_table)
    assert metadata["JonesType"] == "complex"
