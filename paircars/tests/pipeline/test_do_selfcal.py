import pytest
from unittest.mock import patch, MagicMock
from itertools import cycle
from paircars.pipeline.do_selfcal import *


@patch("paircars.pipeline.do_selfcal.drop_cache")
@patch("paircars.pipeline.do_selfcal.clean_shutdown")
@patch("paircars.pipeline.do_selfcal.time.sleep", return_value=None)
@patch("paircars.pipeline.do_selfcal.os.chdir")
@patch("paircars.pipeline.do_selfcal.os.makedirs")
@patch("paircars.pipeline.do_selfcal.os.path.exists", return_value=False)
@patch("paircars.pipeline.do_selfcal.os.system")
@patch(
    "paircars.pipeline.do_selfcal.create_logger",
    return_value=(MagicMock(), "log.log"),
)
@patch("paircars.pipeline.do_selfcal.init_logger")
@patch(
    "paircars.pipeline.do_selfcal.get_unflagged_antennas",
    return_value=(["ant1", "ant2"], [0.1, 0.1]),
)
@patch("paircars.pipeline.do_selfcal.calc_cellsize", return_value=5.0)
@patch("paircars.pipeline.do_selfcal.calc_field_of_view", return_value=1200)
@patch("paircars.pipeline.do_selfcal.flag_non_disk", return_value=0)
@patch("paircars.pipeline.do_selfcal.uvbin_flag", return_value=0)
@patch("paircars.pipeline.do_selfcal.check_datacolumn_valid", return_value=True)
@patch("paircars.pipeline.do_selfcal.msmetadata")
@patch("casatasks.flagmanager", return_value={0: {"name": "applycal"}})
@patch("casatasks.flagdata")
@patch("casatasks.split")
@patch("paircars.pipeline.do_selfcal.limit_threads")
@patch("paircars.pipeline.do_selfcal.quiet_sun_selfcal", return_value=(0, "g0.cal"))
@patch("paircars.pipeline.do_selfcal.selfcal_round")
def test_do_selfcal_function(
    mock_selfcal_round,
    mock_quiet_sun,
    mock_limit_threads,
    mock_split,
    mock_flagdata,
    mock_flagmanager,
    mock_msmetadata,
    mock_check_data,
    mock_uvbin,
    mock_flag_non_disk,
    mock_fov,
    mock_cellsize,
    mock_unflagged,
    mock_init_logger,
    mock_create_logger,
    mock_os_system,
    mock_path_exists,
    mock_makedirs,
    mock_chdir,
    mock_sleep,
    mock_shutdown,
    mock_drop_cache,
):
    mock_msmd = MagicMock()
    mock_msmd.open.return_value = None
    mock_msmd.close.return_value = None
    mock_msmd.scannumbers.return_value = [0]
    mock_msmd.fieldsforscan.return_value = [0]
    mock_msmd.meanfreq.return_value = 100.0
    mock_msmetadata.return_value = mock_msmd
    mock_selfcal_round.side_effect = cycle(
        [
            (0, ["g0.cal"], 100.0, 0.01, "", "", "", None),
            (0, ["g1.cal"], 110.0, 0.009, "", "", "", None),
            (0, ["g2.cal"], 111.0, 0.008, "", "", "", None),
        ]
    )
    status, msname, caltable, _ = do_selfcal(
        msname="mock.ms",
        workdir="/tmp",
        selfcaldir="/tmp",
        metafits="meta.fits",
        min_iter=1,
        max_iter=3,
    )
    assert status == 1
    assert isinstance(caltable, list)
    mock_selfcal_round.side_effect = [
        (1, [], 0, 0, "", "", "", None),
        (1, [], 0, 0, "", "", "", None),
    ]
    status, msname, caltable, _ = do_selfcal(
        msname="mock.ms",
        workdir="/tmp",
        selfcaldir="/tmp",
        metafits="meta.fits",
    )
    assert status == 1
    assert caltable == []
    mock_selfcal_round.side_effect = [
        (2, [], 0, 0, "", "", "", None),
    ]
    status, msname, caltable, _ = do_selfcal(
        msname="mock.ms",
        workdir="/tmp",
        selfcaldir="/tmp",
        metafits="meta.fits",
    )
    assert status == 1
    assert caltable == []
    mock_selfcal_round.side_effect = cycle(
        [
            (0, ["g0.cal"], 100.0, 0.01, "", "", "", None),
            (0, ["g1.cal"], 150.0, 0.009, "", "", "", None),
            (0, ["g2.cal"], 80.0, 0.009, "", "", "", None),
        ]
    )
    status, msname, caltable, _ = do_selfcal(
        msname="mock.ms",
        workdir="/tmp",
        selfcaldir="/tmp",
        metafits="meta.fits",
        min_iter=1,
        do_apcal=True,
    )
    assert status == 1
    assert isinstance(caltable, list)
    mock_selfcal_round.side_effect = cycle(
        [(0, [f"g{i}.cal"], 100.0, 0.01, "", "", "", None) for i in range(20)]
    )
    status, msname, caltable, _ = do_selfcal(
        msname="mock.ms",
        workdir="/tmp",
        selfcaldir="/tmp",
        metafits="meta.fits",
        max_iter=5,
        min_iter=1,
    )
    assert status == 1
    assert isinstance(caltable, list)
    mock_selfcal_round.side_effect = Exception("Simulated failure")
    status, msname, caltable, _ = do_selfcal(
        msname="mock.ms",
        workdir="/tmp",
        selfcaldir="/tmp",
        metafits="meta.fits",
    )
    assert status == 1
    assert caltable == []


@patch("paircars.pipeline.do_selfcal.clean_shutdown")
@patch("paircars.pipeline.do_selfcal.time.sleep", return_value=None)
@patch("paircars.pipeline.do_selfcal.os.chdir")
@patch("paircars.pipeline.do_selfcal.os.makedirs")
@patch("paircars.pipeline.do_selfcal.os.path.exists", return_value=False)
@patch("paircars.pipeline.do_selfcal.os.system")
@patch(
    "paircars.pipeline.do_selfcal.create_logger",
    return_value=(MagicMock(), "log.log"),
)
@patch("paircars.pipeline.do_selfcal.init_logger")
@patch(
    "paircars.pipeline.do_selfcal.get_unflagged_antennas",
    return_value=(["ant1", "ant2"], [0.1, 0.1]),
)
@patch("paircars.pipeline.do_selfcal.calc_cellsize", return_value=5.0)
@patch("paircars.pipeline.do_selfcal.calc_field_of_view", return_value=1200)
@patch("paircars.pipeline.do_selfcal.uvbin_flag", return_value=0)
@patch("paircars.pipeline.do_selfcal.check_datacolumn_valid", return_value=True)
@patch("paircars.pipeline.do_selfcal.msmetadata")
@patch("casatasks.flagmanager",return_value=True)
@patch("casatasks.flagdata")
@patch("casatasks.split")
@patch("paircars.pipeline.do_selfcal.limit_threads")
@patch("paircars.pipeline.do_selfcal.selfcal_round")
@patch("paircars.pipeline.do_selfcal.do_flag_backup",return_value=True)
@patch("paircars.pipeline.do_selfcal.flag_non_disk",return_value=0)
@patch("paircars.pipeline.do_selfcal.get_chans_flag",return_value= ([0],[]))
@patch("paircars.pipeline.do_selfcal.weighted_mean",return_value= (1,0.1))
def test_do_polselfcal(
    mock_mean,
    mock_get_chans_flag,
    mock_flag_non_disk,
    mock_flagbackup,
    mock_selfcal_round,
    mock_limit_threads,
    mock_split,
    mock_flagdata,
    mock_flagmanager,
    mock_msmetadata,
    mock_check_data,
    mock_uvbin,
    mock_fov,
    mock_cellsize,
    mock_unflagged,
    mock_init_logger,
    mock_create_logger,
    mock_os_system,
    mock_path_exists,
    mock_makedirs,
    mock_chdir,
    mock_sleep,
    mock_shutdown,
):
    mock_msmd = MagicMock()
    mock_msmd.open.return_value = None
    mock_msmd.close.return_value = None
    mock_msmd.scannumbers.return_value = [0]
    mock_msmd.fieldsforscan.return_value = [0]
    mock_msmd.meanfreq.return_value = 100.0
    mock_msmd.timesforspws.return_value = [0, 10]
    mock_msmd.exposuretime.return_value = {"value": 10}
    mock_msmetadata.return_value = mock_msmd
    mock_selfcal_round.side_effect = cycle(
        [
            (0, ["g0.cal"], 100.0, 0.01, "", "", "", [[0.1, 0.1, 0.1, 0, 0, 0]]),
            (0, ["g1.cal"], 110.0, 0.009, "", "", "", [[0.05, 0.05, 0.05, 0, 0, 0]]),
            (0, ["g2.cal"], 111.0, 0.008, "", "", "", [[0.01, 0.01, 0.01, 0, 0, 0]]),
        ]
    )
    status, msname, caltable, leakagetable = do_polselfcal(
        msname="mock.ms",
        workdir="/tmp",
        selfcaldir="/tmp",
        metafits="meta.fits",
        min_iter=1,
        max_iter=3,
    )
    assert status == 0
    assert isinstance(caltable, list)
    mock_selfcal_round.side_effect = [(1, [], 0, 0, "", "", "", [[0, 0, 0, 0, 0, 0]])]
    status, msname, caltable, leakagetable = do_polselfcal(
        msname="mock.ms",
        workdir="/tmp",
        selfcaldir="/tmp",
        metafits="meta.fits",
    )
    assert status == 1
    assert caltable == []
    mock_selfcal_round.side_effect = [(2, [], 0, 0, "", "", "", [[0, 0, 0, 0, 0, 0]])]
    status, msname, caltable, leakagetable = do_polselfcal(
        msname="mock.ms",
        workdir="/tmp",
        selfcaldir="/tmp",
        metafits="meta.fits",
    )
    assert status == 1
    assert caltable == []
    mock_selfcal_round.side_effect = cycle(
        [
            (0, ["g0.cal"], 100.0, 0.01, "", "", "", [[0.01, 0.01, 0.01, 0, 0, 0]]),
            (0, ["g1.cal"], 200.0, 0.009, "", "", "", [[0.01, 0.01, 0.01, 0, 0, 0]]),
            (0, ["g2.cal"], 80.0, 0.009, "", "", "", [[0.01, 0.01, 0.01, 0, 0, 0]]),
        ]
    )
    status, msname, caltable, leakagetable = do_polselfcal(
        msname="mock.ms",
        workdir="/tmp",
        selfcaldir="/tmp",
        metafits="meta.fits",
        min_iter=1,
    )
    assert status == 0
    assert isinstance(caltable, list)
    mock_selfcal_round.side_effect = cycle(
        [
            (0, ["g0.cal"], 1000000.0, 0.01, "", "", "", [[0.01, 0.01, 0.01, 0, 0, 0]]),
        ]
    )
    status, msname, caltable, leakagetable = do_polselfcal(
        msname="mock.ms",
        workdir="/tmp",
        selfcaldir="/tmp",
        metafits="meta.fits",
        min_iter=1,
    )
    assert status == 0
    assert isinstance(caltable, list)
    mock_selfcal_round.side_effect = cycle(
        [
            (0, ["gX.cal"], 100.0, 0.01, "", "", "", [[0.02, 0.02, 0.02, 0, 0, 0]])
            for _ in range(20)
        ]
    )
    status, msname, caltable, leakagetable = do_polselfcal(
        msname="mock.ms",
        workdir="/tmp",
        selfcaldir="/tmp",
        metafits="meta.fits",
        max_iter=3,
        min_iter=1,
    )
    assert status == 0
    assert isinstance(caltable, list)
    mock_selfcal_round.side_effect = Exception("Simulated failure")
    status, msname, caltable, leakagetable = do_polselfcal(
        msname="mock.ms",
        workdir="/tmp",
        selfcaldir="/tmp",
        metafits="meta.fits",
    )
    assert status == 1
    assert caltable == []


@patch("paircars.pipeline.do_selfcal.do_polselfcal")
@patch("paircars.pipeline.do_selfcal.do_selfcal")
@patch("paircars.pipeline.do_selfcal.msmetadata")
@patch(
    "paircars.pipeline.do_selfcal.get_unflagged_antennas",
    return_value=(["ant1", "ant2"], [0.1, 0.1]),
)
def test_do_full_selfcal(
    mock_unflagged,
    mock_msmetadata,
    mock_do_selfcal,
    mock_do_polselfcal,
):
    mock_msmd = MagicMock()
    mock_msmd.open.return_value = None
    mock_msmd.close.return_value = None
    mock_msmd.antennaids.return_value = [0]
    mock_msmetadata.return_value = mock_msmd
    mock_do_selfcal.return_value = (1, "int.ms", [], False)
    result = do_full_selfcal(
        msname="mock.ms",
        workdir="/tmp",
        selfcaldir="/tmp/selfcal",
        metafits="meta.fits",
    )
    assert result == (1, 1, [], [], "")
    mock_do_selfcal.return_value = (0, "int.ms", ["int.cal"], True)
    result = do_full_selfcal(
        msname="mock.ms",
        workdir="/tmp",
        selfcaldir="/tmp/selfcal",
        metafits="meta.fits",
        do_polcal=False,
    )
    assert result == (0, 1, ["int.cal"], [], "")
    mock_do_selfcal.return_value = (0, "int.ms", ["int.cal"], True)
    mock_do_polselfcal.return_value = (0, "pol.ms", ["pol.cal"], "leakage.npy")
    result = do_full_selfcal(
        msname="mock.ms",
        workdir="/tmp",
        selfcaldir="/tmp/selfcal",
        metafits="meta.fits",
    )
    assert result == (0, 0, ["int.cal"], ["pol.cal"], "leakage.npy")
    mock_do_selfcal.return_value = (0, "int.ms", ["int.cal"], True)
    mock_do_polselfcal.return_value = (2, "pol.ms", [], "")
    result = do_full_selfcal(
        msname="mock.ms",
        workdir="/tmp",
        selfcaldir="/tmp/selfcal",
        metafits="meta.fits",
    )
    assert result == (0, 2, ["int.cal"], [], "")


@patch("paircars.pipeline.do_selfcal.clean_shutdown")
@patch("paircars.pipeline.do_selfcal.drop_cache")
@patch("paircars.pipeline.do_selfcal.time.sleep", return_value=None)
@patch("paircars.pipeline.do_selfcal.os.system")
@patch("paircars.pipeline.do_selfcal.os.makedirs")
@patch("paircars.pipeline.do_selfcal.os.path.exists", return_value=True)
@patch("paircars.pipeline.do_selfcal.get_local_dask_cluster")
@patch("paircars.pipeline.do_selfcal.scale_worker_and_wait")
@patch("paircars.pipeline.do_selfcal.check_udocker_container", return_value=True)
@patch("paircars.pipeline.do_selfcal.initialize_wsclean_container")
@patch("paircars.pipeline.do_selfcal.fits.getheader", return_value={"GPSTIME": 12345})
@patch("paircars.pipeline.do_selfcal.check_datacolumn_valid", return_value=True)
@patch("paircars.pipeline.do_selfcal.msmetadata")
@patch("paircars.pipeline.do_selfcal.get_caltable_metadata")
@patch("paircars.pipeline.do_selfcal.get_quartical_table_metadata")
@patch("paircars.pipeline.do_selfcal.freq_to_MWA_coarse", return_value=10)
@patch("paircars.pipeline.do_selfcal.do_full_selfcal")
@patch("paircars.pipeline.do_selfcal.os.chdir",return_value=True)
@patch("paircars.pipeline.do_selfcal.get_ncoarse",return_value=1)
def test_main_selfcal(
    mock_ncoarse,
    mock_chdir,
    mock_do_full_selfcal,
    mock_freq_to_coarse,
    mock_get_quart_meta,
    mock_get_cal_meta,
    mock_msmetadata,
    mock_check_data,
    mock_getheader,
    mock_init_container,
    mock_check_container,
    mock_scale,
    mock_get_cluster,
    mock_path_exists,
    mock_makedirs,
    mock_os_system,
    mock_sleep,
    mock_drop_cache,
    mock_shutdown,
):
    mock_client = MagicMock()
    mock_cluster = MagicMock()
    mock_client.compute.side_effect = lambda tasks: tasks
    mock_client.gather.side_effect = lambda tasks: tasks
    mock_get_cluster.return_value = (mock_client, mock_cluster, "/tmp/dask", 1)
    mock_msmd = MagicMock()
    mock_msmd.open.return_value = None
    mock_msmd.close.return_value = None
    mock_msmd.timesforspws.return_value = [0, 10]
    mock_msmetadata.return_value = mock_msmd
    mock_get_cal_meta.return_value = {
        "Channel 0 frequency (MHz)": 100,
        "Bandwidth (MHz)": 10,
    }
    mock_get_quart_meta.return_value = {
        "Channel 0 frequency (MHz)": 100,
        "Bandwidth (MHz)": 10,
    }
    mock_do_full_selfcal.return_value = (
        0,
        0,
        ["g.cal", "b.cal"],
        ["d.cal"],
    )
    msg = main(
        mslist="mock.ms",
        metafits="meta.fits",
        workdir="/tmp/work",
        caldir="/tmp/cal",
    )
    mock_do_full_selfcal.return_value = (
        1,
        1,
        [],
        [],
    )
    msg = main(
        mslist="mock.ms",
        metafits="meta.fits",
        workdir="/tmp/work",
        caldir="/tmp/cal",
    )
    mock_check_data.return_value = False
    msg = main(
        mslist="mock.ms",
        metafits="meta.fits",
        workdir="/tmp/work",
        caldir="/tmp/cal",
    )
    mock_check_container.return_value = False
    mock_init_container.return_value = None
    msg = main(
        mslist="mock.ms",
        metafits="meta.fits",
        workdir="/tmp/work",
        caldir="/tmp/cal",
    )
    mock_do_full_selfcal.side_effect = Exception("Simulated failure")
    msg = main(
        mslist="mock.ms",
        metafits="meta.fits",
        workdir="/tmp/work",
        caldir="/tmp/cal",
    )


@pytest.mark.parametrize(
    "argv, should_exit",
    [
        (["prog.py"], True),  # No args → help
        (
            [
                "prog.py",
                "ms1.ms,ms2.ms",
                "--workdir",
                "/mock/work",
                "--caldir",
                "/mock/caltables",
                "--start_thresh",
                "5",
                "--stop_thresh",
                "3",
                "--no_apcal",
                "--keep_backup",
            ],
            False,
        ),
    ],
)
@patch("paircars.pipeline.do_selfcal.main", return_value= (0, 0, 0, 0, 0))
@patch("paircars.pipeline.do_selfcal.sys.exit")
@patch("paircars.pipeline.do_selfcal.argparse.ArgumentParser.print_help")
def test_cli_selfcal(mock_print_help, mock_exit, mock_main, argv, should_exit):
    with patch("sys.argv", argv):
        from paircars.pipeline import do_selfcal

        result = do_selfcal.cli()
        assert result == should_exit
