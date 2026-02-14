import pytest
import traceback
import os
import numpy as np
from casatasks import casalog
from casatools import table
from paircars.utils.calibration import *
from unittest.mock import MagicMock, patch

try:
    casalogfile = casalog.logfile()
    os.system("rm -rf " + casalogfile)
except BaseException:
    traceback.print_exc()
    pass


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


def test_merge_caltables_from_fixture(dummy_caltables, tmp_path):
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
