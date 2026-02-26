import pytest
import os
from unittest.mock import patch, MagicMock, mock_open, call
from paircars.utils.basic_utils import *


def test_suppress_output_fd():
    with suppress_output():
        os.write(1, b"This should not appear\n")
        os.write(2, b"This error should not appear\n")


def test_get_cachedir(mocker):
    dummy_home = "/dummy/home"
    dummy_user = "dummyuser"
    expected_cachedir = f"{dummy_home}/.paircarspipe"
    mocker.patch.dict("os.environ", {"HOME": dummy_home})
    mocker.patch("os.getlogin", return_value=dummy_user)
    makedirs_mock = mocker.patch("os.makedirs")
    cachedir = get_cachedir()
    assert cachedir == expected_cachedir
    makedirs_mock.assert_any_call(expected_cachedir, exist_ok=True)


@pytest.mark.parametrize(
    "input_datadir, cachedir, expected_datadir",
    [
        ("", "/mock/cache", "/mock/cache/paircarspipe_data"),  # default
        ("/custom/data", "/mock/cache", "/custom/data/paircarspipe_data"),  # user-provided
    ],
)
@patch("paircars.utils.basic_utils.open", new_callable=mock_open)
@patch("paircars.utils.basic_utils.os.makedirs")
@patch("paircars.utils.basic_utils.get_cachedir")
def test_create_datadir(
    mock_get_cachedir,
    mock_makedirs,
    mock_open_file,
    input_datadir,
    cachedir,
    expected_datadir,
):
    mock_get_cachedir.return_value = cachedir

    create_datadir(datadir=input_datadir)

    # Check directory creation
    mock_makedirs.assert_called_once_with(expected_datadir, exist_ok=True)

    # Check file write
    mock_open_file.assert_called_once_with(
        f"{cachedir}/paircarspipe_data_dir.txt", "w"
    )
    mock_open_file().write.assert_called_once_with(expected_datadir + "\n")


@pytest.mark.parametrize(
    "file_exists, file_contents, expected_result",
    [
        (True, "/custom/data\n", "/custom/data"),
        (False, "", None),
    ],
)
@patch("paircars.utils.basic_utils.get_cachedir", return_value="/mock/cache")
@patch("paircars.utils.basic_utils.os.makedirs")
@patch("paircars.utils.basic_utils.os.path.exists")
def test_get_datadir(
    mock_exists,
    mock_makedirs,
    mock_get_cachedir,
    file_exists,
    file_contents,
    expected_result,
):
    open_path = "paircars.utils.basic_utils.open"
    with patch(open_path, mock_open(read_data=file_contents)) as mock_open_func:
        mock_exists.return_value = file_exists

        result = get_datadir()

        mock_get_cachedir.assert_called_once()
        mock_exists.assert_called_once_with("/mock/cache/paircarspipe_data_dir.txt")

        if file_exists:
            mock_open_func.assert_called_once_with(
                "/mock/cache/paircarspipe_data_dir.txt", "r"
            )
            mock_open_func().read.assert_called_once()
            mock_makedirs.assert_called_once_with(expected_result, exist_ok=True)
            assert result == expected_result
        else:
            mock_open_func.assert_not_called()
            mock_makedirs.assert_not_called()
            assert result is None

@pytest.mark.parametrize(
    "time_side_effect, conn_side_effect, expected_result, expected_calls, expected_sleep_calls",
    [
        # Immediate success
        (
            [0, 1],  # time progression
            [MagicMock()],  # connection succeeds immediately
            True,
            1,  # create_connection calls
            0,  # sleep calls
        ),
        # Success after one retry
        (
            [0, 1, 3],
            [OSError("Not ready"), MagicMock()],
            True,
            2,
            1,
        ),
        # Timeout
        (
            [0, 30, 61],
            OSError("Still closed"),
            False,
            1,  # at least one attempt
            1,
        ),
    ],
)
@patch("paircars.utils.basic_utils.time.sleep")
@patch("paircars.utils.basic_utils.time.time")
@patch("paircars.utils.basic_utils.socket.create_connection")
def test_wait_for_port(
    mock_create_conn,
    mock_time,
    mock_sleep,
    time_side_effect,
    conn_side_effect,
    expected_result,
    expected_calls,
    expected_sleep_calls,
):
    mock_time.side_effect = time_side_effect
    mock_create_conn.side_effect = conn_side_effect

    result = wait_for_port("localhost", 1234, timeout=60)

    assert result is expected_result
    assert mock_create_conn.call_count == expected_calls
    assert mock_sleep.call_count == expected_sleep_calls
    

@pytest.mark.parametrize(
    "bind_side_effects, expected_port",
    [
        # First port is free
        (
            [None],  # bind succeeds immediately
            4200,
        ),
        # First busy, second free
        (
            [OSError("Busy"), None],
            4201,
        ),
        # All busy
        (
            [OSError("Busy"), OSError("Busy")],
            None,
        ),
    ],
)
@patch("paircars.utils.basic_utils.socket.socket")
def test_get_free_port(
    mock_socket_class,
    bind_side_effects,
    expected_port,
):
    # Create mock socket instance
    mock_socket = MagicMock()
    mock_socket.__enter__.return_value = mock_socket
    mock_socket.__exit__.return_value = None

    # Set bind behavior
    mock_socket.bind.side_effect = bind_side_effects

    mock_socket_class.return_value = mock_socket

    result = get_free_port(start_port=4200, end_port=4202)

    assert result == expected_port
    
    
@pytest.mark.parametrize(
    "bind_side_effect, expected_result",
    [
        # Port free
        (None, True),
        # Port occupied
        (OSError("Address already in use"), False),
    ],
)
@patch("paircars.utils.basic_utils.socket.socket")
def test_check_port_status(
    mock_socket_class,
    bind_side_effect,
    expected_result,
):
    # Mock socket instance
    mock_socket = MagicMock()
    mock_socket.__enter__.return_value = mock_socket
    mock_socket.__exit__.return_value = None

    mock_socket.bind.side_effect = bind_side_effect

    mock_socket_class.return_value = mock_socket

    result = check_port_status(4200)

    assert result is expected_result
    mock_socket.bind.assert_called_once_with(("127.0.0.1", 4200))
    
        
def test_ra_dec_to_deg():
    radeg, decdeg = ra_dec_to_deg("00h00m00s", "00d00m00s")
    assert radeg == 0
    assert decdeg == 0
    radeg, decdeg = ra_dec_to_deg("01h30m03s", "+12d02m30s")
    assert radeg == 22.5125
    assert decdeg == 12.0417
    radeg, decdeg = ra_dec_to_deg("01h30m03s", "-12d02m30s")
    assert radeg == 22.5125
    assert decdeg == -12.0417


def test_ra_dec_to_hms_dms():
    ra, dec = ra_dec_to_hms_dms(0.0, 0.0)
    assert ra == "0h00m00s"
    assert dec == "0d00m00s"
    ra, dec = ra_dec_to_hms_dms(22.5125, -12.0417)
    assert ra == "1h30m03s"
    assert dec == "-12d02m30.12s"
    a, dec = ra_dec_to_hms_dms(22.5125, 12.0417)
    assert ra == "1h30m03s"
    assert dec == "12d02m30.12s"


@pytest.mark.parametrize(
    "lst,target_chunk_size,result",
    [
        ([1, 2, 3, 4, 5], 2, [[1, 2, 3], [4, 5]]),
        ([1, 2, 3, 4, 5], 3, [[1, 2, 3], [4, 5]]),
        ([1, 2, 3], 1, [[1], [2], [3]]),
        ([1], 1, [[1]]),
        ([], 1, [[]]),
    ],
)
def test_split_into_chunks(lst, target_chunk_size, result):
    assert split_into_chunks(lst, target_chunk_size) == result


@pytest.mark.parametrize(
    "timestamps,expected",
    [
        (
            ["2014-05-01T00:00:00", "2014-05-01T01:00:00", "2014-05-01T02:00:00"],
            "2014-05-01T01:00:00",
        ),
        (
            ["2024-02-29T00:00:00", "2024-02-28T23:59:59", "2024-03-01T00:00:01"],
            "2024-02-29T08:00:00",
        ),
        ([], ""),
    ],
)
def test_average_timestamp(timestamps, expected):
    assert average_timestamp(timestamps) == expected


@pytest.mark.parametrize("n,base,result", [(10, 2, 12), (16, 3, 18), (19, 5, 20)])
def test_ceil_to_multiple(n, base, result):
    assert ceil_to_multiple(n, base) == result


@pytest.mark.parametrize(
    "ra1, dec1, ra2, dec2, expected",
    [
        (10.0, -30.0, 10.0, -30.0, 0.0),
        (0.0, 0.0, 90.0, 0.0, 90.0),
        (0.0, 0.0, 180.0, 0.0, 180.0),
        (10.0, 20.0, 30.0, 40.0, 26.33),
        (120.0, -45.0, 130.0, -44.0, 7.2),
    ],
)
def test_angular_separation_equatorial(ra1, dec1, ra2, dec2, expected):
    result = angular_separation_equatorial(ra1, dec1, ra2, dec2)
    assert result == expected


@pytest.mark.parametrize(
    "timestamp, date_format",
    [
        ("2024/06/30/12:34:56", 0),
        ("2024-06-30T12:34:56", 1),
        ("2024-06-30 12:34:56", 2),
        ("2024_06_30_12_34_56", 3),
    ],
)
def test_timestamp_to_mjdsec(timestamp, date_format):
    expected = 5226467696.0
    result = timestamp_to_mjdsec(timestamp, date_format)
    assert result == expected


def test_mjdsec_to_timestamp():
    assert mjdsec_to_timestamp(5226467696.0, 0) == "2024-06-30T12:34:56.00"
    assert mjdsec_to_timestamp(5226467696.0, 1) == "2024/06/30/12:34:56.00"
    assert mjdsec_to_timestamp(5226467696.0, 2) == "2024-06-30 12:34:56.00"
