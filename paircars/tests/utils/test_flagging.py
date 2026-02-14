import pytest
import traceback
import os
from casatasks import casalog
from casatools import table
from unittest.mock import patch, MagicMock, call
from paircars.utils.flagging import *

try:
    casalogfile = casalog.logfile()
    os.system("rm -rf " + casalogfile)
except BaseException:
    traceback.print_exc()
    pass


def test_flagsummary(dummy_msname):
    summary_file = "sum.txt"
    result = flagsummary(dummy_msname, summary_file)
    assert result == summary_file
    os.system(f"rm -rf {result}")


def test_do_flag_backup(dummy_msname):
    from casatasks import flagmanager

    do_flag_backup(dummy_msname, flagtype="test_flagdata")
    flags = flagmanager(vis=dummy_msname, mode="list")
    flagged = False
    for f in flags:
        if f != "MS":
            ver_name = flags[f]["name"]
            if "test_flagdata" in ver_name:
                flagmanager(vis=dummy_msname, mode="delete", versionname=ver_name)
                if flagged is not True:
                    flagged = True
    assert flagged == True


def test_get_unflagged_antennas(dummy_msname):
    tb = table()
    tb.open(dummy_msname, nomodify=False)
    flag = tb.getcol("FLAG")
    flag *= False
    tb.putcol("FLAG", flag)
    tb.flush()
    tb.close()
    antlist, fraclist = get_unflagged_antennas(dummy_msname)
    antlist = antlist.tolist()
    fraclist = fraclist.tolist()
    assert fraclist == [0] * len(antlist)


@pytest.mark.parametrize(
    "mode,flagbackup",
    [
        ("rflag", True),
        ("tfcrop", True),
        ("rflag", False),
        ("tfcrop", False),
    ],
)
@patch("paircars.utils.flagging.traceback.print_exc")
@patch("casatasks.flagmanager")
@patch("casatasks.flagdata")
@patch("paircars.utils.flagging.do_flag_backup")
@patch("paircars.utils.flagging.calc_maxuv")
@patch("paircars.utils.flagging.suppress_output")
def test_uvbin_flag(
    mock_suppress,
    mock_calc_maxuv,
    mock_do_flag_backup,
    mock_flagdata,
    mock_flagmanager,
    mock_print_exc,
    mode,
    flagbackup,
):
    # setup
    mock_calc_maxuv.return_value = (1000.0, 100)
    mock_suppress.return_value.__enter__.return_value = None
    mock_suppress.return_value.__exit__.return_value = None

    ret = uvbin_flag(
        "test.ms",
        uvbin_size=10,
        mode=mode,
        threshold=5.0,
        flagbackup=flagbackup,
    )

    assert ret == 0
    mock_calc_maxuv.assert_called_once_with("test.ms")

    if flagbackup:
        mock_do_flag_backup.assert_called_once_with(
            "test.ms", flagtype="uvbin_flagdata"
        )
    else:
        mock_do_flag_backup.assert_not_called()

    # flagdata should be called multiple times
    assert mock_flagdata.call_count > 0

    # no restore/delete on success
    mock_flagmanager.assert_not_called()
    mock_print_exc.assert_not_called()


def test_get_chans_flag(dummy_msname):
    unflag_chans, flag_chans = get_chans_flag(dummy_msname)
    assert unflag_chans == [0]
    assert flag_chans == []


def test_calc_flag_fraction(dummy_msname):
    frac = calc_flag_fraction(dummy_msname)
    assert frac == 0.0


@pytest.mark.parametrize(
    "uvrange", ["100~200", ">100", "<200", "200~300lambda", ">200lambda", "<300lambda"]
)
def test_flag_outside_uvrange(dummy_msname, uvrange):
    assert flag_outside_uvrange(dummy_msname, uvrange) == 0
    tb = table()
    tb.open(dummy_msname, nomodify=False)
    flag = tb.getcol("FLAG")
    flag *= False
    tb.putcol("FLAG", flag)
    tb.flush()
    tb.close()


def test_flag_quartical_table(dummy_quartical_table):
    result = flag_quartical_table(dummy_quartical_table, threshold=1000.0)
    assert result == dummy_quartical_table
