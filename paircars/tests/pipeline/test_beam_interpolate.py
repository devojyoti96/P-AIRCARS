import pytest
from unittest.mock import patch, MagicMock
from paircars.pipeline.beam_interpolate import *


@pytest.mark.parametrize(
    "file_exists, raise_exc",
    [
        (True, False),
        (False, False),
        (False, True),
    ],
)
def test_do_beam_interpolate(file_exists, raise_exc):
    fake_h5_read = MagicMock()
    fake_h5_write = MagicMock()
    fake_keys = ["X1_1000000"]
    fake_h5_read.keys.return_value = fake_keys

    def fake_getitem(key):
        if key == "modes":
            return np.ones((1, 1))
        else:
            return np.ones((2, 3))

    fake_h5_read.__getitem__.side_effect = fake_getitem
    with (
        patch(
            "paircars.pipeline.beam_interpolate.os.path.exists",
            return_value=file_exists,
        ),
        patch("paircars.pipeline.beam_interpolate.os.system") as m_system,
        patch("paircars.pipeline.beam_interpolate.h5py.File") as m_h5,
        patch("paircars.pipeline.beam_interpolate.np.memmap") as m_memmap,
        patch("paircars.pipeline.beam_interpolate.CubicSpline") as m_spline,
        patch("paircars.pipeline.beam_interpolate.time.time", return_value=0),
        patch("paircars.pipeline.beam_interpolate.traceback.print_exc"),
    ):
        if raise_exc:
            m_h5.side_effect = Exception("boom")
        else:
            m_h5.side_effect = [fake_h5_read, fake_h5_write]
            m_memmap.return_value = np.zeros((2, 16, 1, 2, 3))
            spline_instance = MagicMock()
            spline_instance.return_value = np.ones(1)
            m_spline.return_value = spline_instance
        result = do_beam_interpolate("test.h5", new_freq_res=160)


@pytest.mark.parametrize(
    "argv_args, expect_main_called, expected_exit",
    [
        (["prog", "--mslist", "a.ms"], True, 0),
        (["prog"], False, 1),
    ],
)
@patch("paircars.pipeline.beam_interpolate.do_beam_interpolate", return_value=0)
@patch("paircars.pipeline.beam_interpolate.sys.exit")
@patch("paircars.pipeline.beam_interpolate.argparse.ArgumentParser.print_help")
def test_cli(
    mock_print_help,
    mock_exit,
    mock_main,
    argv_args,
    expect_main_called,
    expected_exit,
):
    with patch("sys.argv", argv_args):
        from paircars.pipeline import beam_interpolate

        result = beam_interpolate.cli()

        if expect_main_called:
            mock_main.assert_called()
        else:
            mock_print_help.assert_called()

        assert result == expected_exit
