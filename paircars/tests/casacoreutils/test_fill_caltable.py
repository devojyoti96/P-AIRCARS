from unittest.mock import MagicMock, patch

from paircars.casacoreutils.fill_caltable import make_caltable_columns


@patch("paircars.casacoreutils.fill_caltable.makearrcoldesc")
@patch("paircars.casacoreutils.fill_caltable.makecoldesc")
@patch("paircars.casacoreutils.fill_caltable.maketabdesc")
@patch("paircars.casacoreutils.fill_caltable.table")
def test_make_caltable_columns(
    m_table,
    m_maketabdesc,
    m_makecoldesc,
    m_makearrcoldesc,
):
    """
    Test make_caltable_columns with mocked casacore tables.
    """

    # ----------------------------
    # Mock ANTENNA table
    # ----------------------------
    antenna_table = MagicMock()
    antenna_table.__len__.return_value = 3  # 3 antennas
    antenna_table.__enter__.return_value = antenna_table

    # ----------------------------
    # Mock SPECTRAL_WINDOW table
    # ----------------------------
    spw_table = MagicMock()
    spw_table.__len__.return_value = 2  # 2 SPWs
    spw_table.__enter__.return_value = spw_table

    # ----------------------------
    # Mock output caltable
    # ----------------------------
    cal_table = MagicMock()
    cal_table.__len__.return_value = 0

    # table() called 3 times:
    # 1. ANTENNA
    # 2. SPECTRAL_WINDOW
    # 3. output caltable
    m_table.side_effect = [antenna_table, spw_table, cal_table]

    # dummy descriptors
    m_makecoldesc.return_value = MagicMock()
    m_makearrcoldesc.return_value = MagicMock()
    m_maketabdesc.return_value = MagicMock()

    # ----------------------------
    # Run function
    # ----------------------------
    result = make_caltable_columns(
        msname="test.ms",
        caltable="test.cal",
        nchan=4,
    )

    # ----------------------------
    # Assertions
    # ----------------------------

    # return value
    assert result == "test.cal"

    # rows added = nants × nspw = 3 × 2 = 6
    assert cal_table.addrows.call_count == 6

    # verify column writes happened
    assert cal_table.putcell.called

    # verify table closed
    cal_table.close.assert_called()
