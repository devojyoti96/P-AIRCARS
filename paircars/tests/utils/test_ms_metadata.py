import pytest
import os
from unittest.mock import patch
from paircars.utils.ms_metadata import *


@pytest.mark.parametrize(
    "uvrange_input, expected_output, expect_error",
    [
        (">200lambda", ["<200lambda"], False),
        ("<100lambda", [">100lambda"], False),
        ("10~1000lambda", ["<10lambda", ">1000lambda"], False),
        ("   >300lambda  ", ["<300lambda"], False),
        ("50~500lambda", ["<50lambda", ">500lambda"], False),
        ("500", None, True),
        ("20-100lambda", None, True),
        ("100~lambda", None, True),
        ("lambda~200", None, True),
        ("", None, True),
        ("300~10lambda", None, True),
    ],
)
def test_get_uvrange_exclude(uvrange_input, expected_output, expect_error):
    if expect_error:
        with pytest.raises(ValueError):
            get_uvrange_exclude(uvrange_input)
    else:
        assert get_uvrange_exclude(uvrange_input) == expected_output


def test_get_phasecenter(dummy_msname):
    ra, dec = get_phasecenter(dummy_msname, fieldID=0)
    assert ra == 5.58562
    assert dec == 2.4158


def test_get_observatory_name(dummy_msname):
    assert get_observatory_name(dummy_msname) == "MWA"


def test_get_observatory_coord(dummy_msname):
    lat, lon, height = get_observatory_coord(dummy_msname)
    assert lat == -26.703
    assert lon == 116.671
    assert height == 377.83


def test_get_timeranges(dummy_msname):
    t = get_timeranges(dummy_msname, 5, 60)
    assert t == ["2015/03/27/06:57:37.00~2015/03/27/06:57:39.00"]


def test_calc_fractional_bandwidth(dummy_msname):
    assert calc_fractional_bandwidth(dummy_msname) == 0.0


def test_baseline_names(dummy_msname):
    bs_names = baseline_names(dummy_msname)
    for bs in bs_names:
        assert "&&" in bs


def test_get_ms_size(dummy_msname):
    autocor_size = get_ms_size(dummy_msname, only_autocorr=True)
    noautocor_size = get_ms_size(dummy_msname, only_autocorr=False)
    assert noautocor_size > autocor_size


def test_get_column_size(dummy_msname):
    autocor_size = get_column_size(dummy_msname, only_autocorr=True)
    noautocor_size = get_column_size(dummy_msname, only_autocorr=False)
    assert noautocor_size > autocor_size


def test_get_ms_scan_size(dummy_msname):
    size = get_ms_scan_size(dummy_msname, 1)
    assert isinstance(size, float)
    assert size > 0


def test_get_chunk_size(dummy_msname):
    assert get_chunk_size(dummy_msname, memory_limit=1) == 1


def test_check_datacolumn_valid(dummy_msname):
    assert check_datacolumn_valid(dummy_msname, datacolumn="DATA") == True


def test_get_bad_ants(dummy_msname):
    ant_list, ant_str = get_bad_ants(dummy_msname)
    assert ant_list == []
    assert ant_str == ""


def test_get_common_spw():
    assert get_common_spw("0:0~100", "0:50~70") == "0:50~70"


def test_scans_in_timerange(dummy_msname):
    assert scans_in_timerange(
        dummy_msname, timerange="2015/03/27/06:57:37.00~2015/03/27/06:57:39.00"
    ) == {1: "2015/03/27/06:57:37.00~2015/03/27/06:57:39.00"}


def test_get_refant(dummy_msname):
    assert get_refant(dummy_msname) == "0"


def test_get_ms_scans(dummy_msname):
    assert get_ms_scans(dummy_msname) == [1]


def test_get_pol_names(dummy_msname):
    pollist = get_pol_names(dummy_msname)
    assert pollist == ["XX", "XY", "YX", "YY"]
