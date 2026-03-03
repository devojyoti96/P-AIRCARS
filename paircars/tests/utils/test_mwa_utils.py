import pytest
import psutil
import numpy as np
import os
import traceback
from casatasks import casalog
from casatools import ms as casamstool, table
from unittest.mock import patch, MagicMock, mock_open, call
from paircars.utils.mwa_utils import *

try:
    casalogfile = casalog.logfile()
    os.system("rm -rf " + casalogfile)
except BaseException:
    traceback.print_exc()
    pass


def test_get_MWA_OBSID(dummy_msname):
    obsid = get_MWA_OBSID(dummy_msname)
    assert obsid == 1111474560


def test_get_ncoarse(dummy_msname):
    ncoarse = get_ncoarse(dummy_msname)
    assert ncoarse == 1


def test_freq_to_MWA_coarse():
    coarse = freq_to_MWA_coarse(128)
    assert coarse == 100


def test_get_MWA_coarse_chan(dummy_msname):
    coarse = get_MWA_coarse_chan(dummy_msname)
    assert coarse == 104


def test_get_MWA_coarse_bands(dummy_msname):
    coarse_bands = get_MWA_coarse_bands(dummy_msname)
    assert coarse_bands == [(0, 0, [0])]


def test_get_bad_chans(dummy_msname):
    assert get_bad_chans(dummy_msname) == ""


def test_get_good_chans(dummy_msname):
    assert get_good_chans(dummy_msname) == "0:0"


def test_get_mwa_bad_ants(dummy_metafits):
    bad_ants = get_mwa_bad_ants(dummy_metafits)
    assert bad_ants == "Tile103,Tile108"


def test_download_MWA_metafits():
    metafits = download_MWA_metafits(1111474560)
    assert os.path.exists(metafits) is True
    assert metafits == "./1111474560.metafits"
    os.system(f"rm -rf {metafits}")
