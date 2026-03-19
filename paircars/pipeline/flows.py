import traceback
import glob
import os
from prefect import flow
from paircars.utils.calibration import (
    interpolate_bpass,
    get_caltable_metadata,
)
from paircars.utils.mwa_ploting_utils import (
    plot_caltable_diagnostics,
)
from paircars.utils.mwa_utils import (
    freq_to_MWA_coarse,
)
from paircars.pipeline.tasks import (
    run_solar_phasecenter_jobs,
    run_ds_jobs,
    run_target_split_jobs,
    run_flag,
    run_import_model,
    run_basic_cal_jobs,
    run_apply_basiccal_sol,
    run_solar_siderealcor_jobs,
    run_selfcal_jobs,
    run_apply_selfcal_sol,
    run_imaging_jobs,
    run_apply_pbcor,
    run_make_overlay,
    run_make_msplot,
    #send_task_notification,
)


@flow(
    name="Basic calibration",
    description="Perform basic calibration using calibrator observations",
    log_prints=True,
)
def basic_cal_subflow(
    # Core observational inputs
    cal_obsid,
    cal_datadir,
    cal_metafits,
    coarse_chans,
    target_obsid,
    # I/O and workspace
    workdir,
    cal_outdir,
    basic_caldir,
    # Calibration controls
    do_basic_cal,
    redo_basic_cal,
    do_cal_flag,
    do_import_model,
    do_polcal,
    keep_backup,
    # Data conditioning
    quack_timestamps,
    # Resource management
    cpu_frac,
    mem_frac,
    # Logging / metadata
    jobid,
    #timestamp,
    #emails,
    remote_logger,
):
    """
    Basic calibration sub flow
    """
    ##########################################
    # Checking presence of basic caltables
    ##########################################
    print(f"Performing calibration for calibrator with OBSID: {cal_obsid}")
    if not redo_basic_cal:
        print(
            f"Searching for existing bandpass tables: {basic_caldir}/calibrator_{cal_obsid}*.bcal"
        )
        bandpass_tables = sorted(
            glob.glob(f"{basic_caldir}/calibrator_{cal_obsid}*.bcal")
        )
        print(
            f"Searching for existing crossphase tables: {basic_caldir}/calibrator_{cal_obsid}*.kcrossscal"
        )
        crossphase_tables = sorted(
            glob.glob(f"{basic_caldir}/calibrator_{cal_obsid}*.kcrosscal")
        )
        if len(bandpass_tables) < len(coarse_chans):
            bpass_coarse_chans = []
            for bpass in bandpass_tables:
                cal_metadata = get_caltable_metadata(bpass)
                freqMHz = cal_metadata["Channel 0 frequency (MHz)"]
                bpass_coarse_chans.append(freq_to_MWA_coarse(freqMHz))
            coarse_chans = [x for x in coarse_chans if x not in set(bpass_coarse_chans)]
        else:
            print("###################################################")
            print(
                f"Bandpass tables are already present in calibration directory: {basic_caldir}"
            )
            for bpass in bandpass_tables:
                print(f"{os.path.basename(bpass)}")
            if len(crossphase_tables) < len(coarse_chans):
                kcross_coarse_chans = []
                for kcross in crossphase_tables:
                    cal_metadata = get_caltable_metadata(kcross)
                    freqMHz = cal_metadata["Channel 0 frequency (MHz)"]
                    kcross_coarse_chans.append(freq_to_MWA_coarse(freqMHz))
                coarse_chans = [
                    x for x in coarse_chans if x not in set(kcross_coarse_chans)
                ]
            else:
                print("####################################################")
                print(
                    f"Crosshand phase tables are already present in calibration directory: {basic_caldir}"
                )
                for kcross in crossphase_tables:
                    print(f"{os.path.basename(kcross)}")
                print("####################################################")
                #if emails != "":
                #    email_msg = f"All gain solutions from calibrator for calibrator OBSID: {cal_obsid} are already present."
                #    send_task_notification(
                #        emails, email_msg, jobid, target_obsid, timestamp
                 #   )
                return 0, bandpass_tables, crossphase_tables

    ############################
    # Calibrator ms list
    ############################
    cal_mslist = glob.glob(f"{cal_datadir}/*.ms")
    if len(cal_mslist) == 0 or len(coarse_chans) == 0:
        print(
            f"No calibrator measurement set with coarse channels: {coarse_chans} is present in: {cal_datadir}"
        )
        return 1, [], []

    print (cal_mslist)
   
