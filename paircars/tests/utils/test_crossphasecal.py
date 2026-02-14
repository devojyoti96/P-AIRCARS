import pytest
import traceback
import os
import numpy as np
from casatasks import casalog
from casatools import table
from paircars.utils.crossphasecal import *
from unittest.mock import MagicMock, patch

try:
    casalogfile = casalog.logfile()
    os.system("rm -rf " + casalogfile)
except BaseException:
    traceback.print_exc()
    pass


def test_create_blank_table(dummy_msname):
    caltable_name = "test.kcross"
    result = create_blank_table(dummy_msname, caltable_name)
    assert result == caltable_name
    os.system(f"rm -rf {result}")


def test_create_crossphase_table(dummy_msname):
    caltable_name = "test.kcross"
    freqs = np.array([132, 133, 134]) * 10**6
    crossphase = np.ones(3)
    flags = crossphase.astype("bool")
    result = create_crossphase_table(
        dummy_msname, caltable_name, freqs, crossphase, flags
    )
    assert result == caltable_name
    tb = table()
    tb.open(result)
    gain = tb.getcol("CPARAM")
    tb.close()
    phase = round(np.nanmedian(np.angle(gain[0, ...], deg=True)), 1)
    assert phase == 1.0
    os.system(f"rm -rf {result}")


def test_fitted_crossphase():
    freqs = np.arange(10) * 10**6
    crossphase = np.ones(10)
    fit_crossphase = fitted_crossphase(freqs, crossphase)
    assert len(fit_crossphase) == 10
