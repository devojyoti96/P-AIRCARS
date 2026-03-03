import pytest
import traceback
import os
import numpy as np
from casatasks import casalog
from casatools import table
from paircars.utils.ds_utils import *
from unittest.mock import MagicMock, patch

try:
    casalogfile = casalog.logfile()
    os.system("rm -rf " + casalogfile)
except BaseException:
    traceback.print_exc()
    pass


def test_calc_T_rec():
    result = calc_T_rec(50.1336)
    assert result == 2363.92


def test_calc_T_pickup():
    result = calc_T_pickup(75.0)
    assert result == 20.142


def test_cal_sun_solid_angle():
    result = cal_sun_solid_angle(100.0)
    assert result == 0.00044164209135847456


def test_cal_norm_crosscorr(dummy_msname):
    result = cal_norm_crosscorr(dummy_msname, 0, 1)
    assert np.nanmedian(result[0]) <= 1


def test_get_short_baselines(dummy_msname):
    result = get_short_baselines(dummy_msname)
    assert result[0] == [0, 1]


@patch("paircars.utils.ds_utils.np.save")
@patch("paircars.utils.ds_utils.radec_sun", return_value=("SUN", 0, 0, 0.0, 0.0))
@patch("paircars.utils.ds_utils.get_bad_chans", return_value="")
@patch(
    "paircars.utils.ds_utils.get_pb_radec", return_value=[None, None, None, 1.0, 1.0]
)
@patch(
    "paircars.utils.ds_utils.make_primarybeammap",
    return_value=(1, 1, 10, 1, 1, 1, 10, 1, 1, 1),
)
@patch(
    "paircars.utils.ds_utils.cal_norm_crosscorr", return_value=np.ones((4, 10)) * 0.1
)
@patch("paircars.utils.ds_utils.cal_sun_solid_angle", return_value=1.0)
@patch("paircars.utils.ds_utils.calc_T_pickup", return_value=5.0)
@patch("paircars.utils.ds_utils.calc_T_rec", return_value=50.0)
@patch("paircars.utils.ds_utils.get_short_baselines", return_value=[[0, 1]])
@patch("paircars.utils.ds_utils.msmetadata")
def test_calc_dynamic_spectrum(
    m_msmd,
    m_baselines,
    m_trec,
    m_tpick,
    m_sunomega,
    m_rn,
    m_pbmap,
    m_pb,
    m_badchans,
    m_radec,
    m_save,
):
    msmd = MagicMock()
    msmd.chanfreqs.return_value = np.linspace(100, 110, 5)
    msmd.ncorrforpol.return_value = [4]
    msmd.meanfreq.return_value = 105
    msmd.chanres.return_value = [320, 320, 320]
    msmd.timesforspws.return_value = np.linspace(0, 10, 10)
    msmd.open.return_value = None
    msmd.close.return_value = None
    m_msmd.return_value = msmd

    ds, rn = calc_dynamic_spectrum(
        "test.ms",
        "test.metafits",
        "/tmp",
        n_threads=1,
    )

    assert ds.endswith("_ds.npy")
    assert rn.endswith("_rn.npy")
    assert m_save.call_count == 2
