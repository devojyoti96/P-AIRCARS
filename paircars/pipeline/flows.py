import traceback
import glob
import os
import numpy as np
from prefect import flow
from astropy.io import fits
from casatools import msmetadata
from paircars.utils.basic_utils import print_banner
from paircars.utils.calibration import (
    interpolate_bpass,
    interpolate_quartical,
    get_caltable_metadata,
    scale_bandpass,
)
from paircars.utils.flagging import get_chans_flag
from paircars.utils.mwa_ploting_utils import (
    plot_caltable_diagnostics,
    plot_quartical_tables,
)
from paircars.utils.mwa_utils import (
    freq_to_MWA_coarse,
    get_selfcal_ntimes,
)
from paircars.utils.ms_metadata import check_datacolumn_valid
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
    send_task_notification,
)
from prefect.context import get_run_context
from multiprocessing import Event
from paircars.utils.prefect_logger_utils import (
    start_flow_log_saver,
)

#########################
# Pre-processing subflow
#########################
@flow(
    name="Pre-processing target",
    description="Perform pre-processing on target measurement sets",
    log_prints=True,
)
def pre_process_subflow(
    # Core observational inputs
    target_mslist,
    target_metafits,
    target_obsid,
    solar_data,
    # I/O and workspace
    workdir,
    target_outdir,
    # Processing controls
    do_move_solarcenter,
    make_ds,
    # Resource management
    cpu_frac,
    mem_frac,
    # Logging / metadata
    jobid,
    timestamp,
    emails,
    remote_logger,
):
    """
    Pre-processing of target measurement set subflow
    
    Returns
    -------
    int
        Flow success message
    list
        Filtered target measurement set list
    """
    logdir = f"{workdir}/logs"
    os.makedirs(logdir, exist_ok=True)
    pre_process_logfile = f"{logdir}/preprocess_subflow_{target_obsid}.log"
    ctx = get_run_context()
    flow_id = str(ctx.flow_run.id)
    flow_name = ctx.flow_run.name
    stop_event = Event()
    log_thread_flow = start_flow_log_saver(
        flow_id, flow_name, pre_process_logfile, poll_interval=3, stop_event=stop_event
    )
    try:
        ########################################
        # Moving phasecenter to the solar center
        ########################################
        if solar_data and do_move_solarcenter:
            if emails != "":
                email_msg = f"[{target_obsid}] Started moving phasecenter to solar center."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                )
            print_banner("Starting task: Moving phasecenter to the Sun.")
            future_movecenter = run_solar_phasecenter_jobs.with_options(
                task_run_name=f"move_solarcenter_{target_obsid}",
            ).submit(
                ",".join(target_mslist),
                workdir,
                jobid=jobid,
                cpu_frac=round(cpu_frac, 2),
                mem_frac=round(mem_frac, 2),
                remote_log=remote_logger,
            )
            try:
                msg, succeed, failed = future_movecenter.result()
                if emails != "":
                    email_msg = f"[{target_obsid}] Moving phasecenter to solar center is done.\nSucceeded: {succeed}, failed: {failed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
                print_banner("Finished task: Moving phasecenter to solar center is done.")
                filtered_ms = []
                for t_ms in target_mslist:
                    t_ms = t_ms.rstrip("/")
                    if os.path.exists(f"{t_ms}/.solarcenter_move_succeed"):
                        filtered_ms.append(t_ms)
                    else:
                        print(
                            f"Issue in moving phasecneter to solar center: {t_ms}"
                        )
                target_mslist = filtered_ms  # Filtered target mslist
            except Exception:
                print(
                    "Error in moving phasecenter to solar center. P-AIRCARS has stopped."
                )
                traceback.print_exc()
                if emails != "":
                    email_msg = f"[{target_obsid}] Error occured in moving phasecenter to solar center."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
                return 1, []
        #######################################
        # Run dynamic spectra making
        #######################################
        if solar_data and make_ds:
            if emails != "":
                email_msg = f"[{target_obsid}] Started making solar dynamic spectra."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                )
            print_banner("Starting task: Making dynamic spectra of solar target.")
            future_maskms = run_ds_jobs.with_options(
                task_run_name=f"make_ds_{target_obsid}",
            ).submit(
                ",".join(target_mslist),
                target_metafits,
                workdir,
                target_outdir,
                jobid=jobid,
                cpu_frac=round(cpu_frac, 2),
                mem_frac=round(mem_frac, 2),
                remote_log=remote_logger,
            )
            try:
                msg, succeed, failed = future_maskms.result()
                if emails != "":
                    email_msg = f"[{target_obsid}] Making solar dynamic spectra are done.\nSucceeded: {succeed}, failed: {failed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
                print_banner("Finished task: Making solar dynamic spectra are done.")
            except Exception:
                print("!!! WARNING : Error in making dynamic spectra. !!!")
                traceback.print_exc()
                if emails != "":
                    email_msg = f"[{target_obsid}] Error occured in making dynamic spectra."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
        
        return 0, target_mslist 
    except Exception:
        traceback.print_exc()
        return 1, []
    finally:
        stop_event.set()
        log_thread_flow.join(timeout=5)


############################
# Basic calibration subflow
############################
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
    target_metafits,
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
    timestamp,
    emails,
    remote_logger,
):
    """
    Basic calibration sub flow
    """
    logdir = f"{workdir}/logs"
    os.makedirs(logdir, exist_ok=True)
    basic_cal_logfile = f"{logdir}/basiccal_subflow_{cal_obsid}.log"
    ctx = get_run_context()
    flow_id = str(ctx.flow_run.id)
    flow_name = ctx.flow_run.name
    stop_event = Event()
    log_thread_flow = start_flow_log_saver(
        flow_id, flow_name, basic_cal_logfile, poll_interval=3, stop_event=stop_event
    )
    try:
        ##########################################
        # Checking presence of basic caltables
        ##########################################
        if not redo_basic_cal:
            print(
                f"Searching for existing bandpass tables:\n{basic_caldir}/calibrator_{cal_obsid}*.bcal"
            )
            bandpass_tables = sorted(
                glob.glob(f"{basic_caldir}/calibrator_{cal_obsid}*.bcal")
            )
            print(
                f"Searching for existing crossphase tables:\n{basic_caldir}/calibrator_{cal_obsid}*.kcrossscal"
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
                print_banner(
                    f"Bandpass tables are already present.\nCalibration directory: {basic_caldir}"
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
                    print_banner(
                        f"Crosshand phase tables are already present.\nCalibration directory: {basic_caldir}"
                    )
                    for kcross in crossphase_tables:
                        print(f"{os.path.basename(kcross)}")
                    if emails != "":
                        email_msg = f"[{cal_obsid}] All gain solutions from calibrator are already present."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                        )
                    return 0, bandpass_tables, crossphase_tables

        ############################
        # Calibrator ms list
        ############################
        cal_mslist = glob.glob(f"{cal_datadir}/*.ms")
        if len(cal_mslist) == 0 or len(coarse_chans) == 0:
            print_banner(
                f"No calibrator measurement set present.\nCoarse channels: {coarse_chans}.\nCalibrator directory: {cal_datadir}"
            )
            if emails != "":
                email_msg = f"[{cal_obsid}] No calibrator measurement set with coarse channels: {coarse_chans} is present in: {cal_datadir}."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                )
            return 1, [], []

        ##############################
        # Run spliting jobs
        ##############################
        # If basic calibration is requested and calibrator ms and metafits are present
        if do_cal_flag or do_import_model or do_basic_cal:
            prefix = "calibrator"
            if emails != "":
                email_msg = f"[{cal_obsid}] Started spliting of calibrator measurement sets."
                send_task_notification(emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}")
            print_banner(
                "Starting task: Spliting of calibrator measurement sets."
            )
            future_cal_split = run_target_split_jobs.with_options(
                task_run_name=f"split_{cal_obsid}"
            ).submit(
                ",".join(cal_mslist),
                cal_metafits,
                workdir,
                datacolumn="data",
                split_coarse_chans=coarse_chans,
                timeres=10.0,
                freqres=0.16,
                prefix=prefix,
                force_split=False,
                time_window=-1,
                time_interval=-1,
                quack_timestamps=quack_timestamps,
                jobid=jobid,
                cpu_frac=float(cpu_frac),
                mem_frac=float(mem_frac),
                remote_log=remote_logger,
            )
            try:
                msg, expected, succeed = future_cal_split.result()
                if emails != "":
                    email_msg = f"[{cal_obsid}] Spliting of calibrator measurement sets are done.\nExpected: {expected}, succeeded: {succeed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
                print_banner(
                    "Finished task: Spliting of calibrator measurement sets are done."
                )
            except Exception:
                print_banner(
                    "!!!! WARNING: Error in spliting calibrator measurement sets. !!!!"
                )
                traceback.print_exc()
                if emails != "":
                    email_msg = f"[{cal_obsid}] Spliting calibrator measurement set is failed."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
                return 1, [], []

        if do_cal_flag or do_import_model or do_basic_cal:
            split_cal_mslist = sorted(
                glob.glob(f"{workdir}/calibrator_{cal_obsid}*_ch_*.ms")
            )
            if len(split_cal_mslist) == 0:
                print_banner(
                    "No splited measurement set is present for basic calibration."
                )
                if emails != "":
                    email_msg = f"[{cal_obsid}] No splited measurement set is present for basic calibration."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
                return 1, [], []

        ##################################
        # Run flagging jobs on calibrators
        ##################################
        # Only if basic calibration is requested
        if do_cal_flag:
            if emails != "":
                email_msg = f"[{cal_obsid}] Started flagging of calibrators."
                send_task_notification(emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}")
            print_banner("Starting task: Flagging calibrators.")
            future_flag = run_flag.with_options(
                task_run_name=f"flag_{cal_obsid}"
            ).submit(
                ",".join(split_cal_mslist),
                cal_metafits,
                workdir,
                cal_outdir,
                flag_calibrators=True,
                jobid=jobid,
                flag_quack=False,
                cpu_frac=round(cpu_frac, 2),
                mem_frac=round(mem_frac, 2),
                remote_log=remote_logger,
            )
            try:
                msg, succeed, failed = future_flag.result()
                if emails != "":
                    email_msg = f"[{cal_obsid}] Flagging of calibrator is done.\nSucceeded: {succeed}, failed: {failed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
                filtered_ms = []
                for c_ms in split_cal_mslist:
                    c_ms = c_ms.rstrip("/")
                    if os.path.exists(f"{c_ms}/.flag_succeed"):
                        filtered_ms.append(c_ms)
                    else:
                        print(f"Issue in flagging of measurement set: {c_ms}")
                split_cal_mslist = filtered_ms  # Filtered target mslist
                print_banner(
                    "Finished task: Flagging of calibrator is done."
                )
            except Exception:
                print_banner(
                    "!!!! WARNING: Flagging error for calibrator. !!!!"
                )
                traceback.print_exc()
                if emails != "":
                    email_msg = f"[{cal_obsid}] Error in flagging calibrators."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )

        #################################
        # Import model
        #################################
        if do_import_model:
            if emails != "":
                email_msg = (
                    f"[{cal_obsid}] Started importing sky model."
                )
                send_task_notification(emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}")
            print_banner("Starting task: Importing model visibilities.")
            future_import_model = run_import_model.with_options(
                task_run_name=f"model_{cal_obsid}"
            ).submit(
                ",".join(split_cal_mslist),
                cal_metafits,
                workdir,
                jobid=jobid,
                cpu_frac=round(cpu_frac, 2),
                mem_frac=round(mem_frac, 2),
                remote_log=remote_logger,
            )
            try:
                msg, succeed, failed = future_import_model.result()
                if emails != "":
                    email_msg = f"[{cal_obsid}] Model import for calibrator is done.\nSucceeded: {succeed}, failed: {failed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
                print_banner(
                    "Finished task: Model import for calibrator is done."
                )
                filtered_ms = []
                for c_ms in split_cal_mslist:
                    c_ms = c_ms.rstrip("/")
                    if os.path.exists(f"{c_ms}/.modeling_succeed"):
                        filtered_ms.append(c_ms)
                    else:
                        print(
                            f"Issue in importing calibrator sky model: {c_ms}"
                        )
                split_cal_mslist = filtered_ms  # Filtered target mslist
            except Exception:
                print_banner(
                    "!!!! WARNING: Error in importing calibrator models.\nNot continuing calibration. !!!!"
                )
                traceback.print_exc()
                if emails != "":
                    email_msg = f"[{cal_obsid}] Error occured in importing model for calibrators.\nNot using calibrator solutions."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
                return 1, [], []

        ###############################
        # Run basic calibration
        ###############################
        if do_basic_cal:
            if emails != "":
                email_msg = f"[{cal_obsid}] Started basic calibration."
                send_task_notification(emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}")
            print_banner("Starting task: Performing basic calibration.")
            future_basical = run_basic_cal_jobs.with_options(
                task_run_name=f"calibration_{cal_obsid}"
            ).submit(
                ",".join(split_cal_mslist),
                cal_metafits,
                workdir,
                cal_outdir,
                perform_polcal=do_polcal,
                jobid=jobid,
                cpu_frac=round(cpu_frac, 2),
                mem_frac=round(mem_frac, 2),
                keep_backup=keep_backup,
                remote_log=remote_logger,
            )
            try:
                msg, succeed, failed = future_basical.result()
                if emails != "":
                    email_msg = f"[{cal_obsid}] Basic calibration is done.\nSucceeded: {succeed}, failed: {failed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
                print_banner("Finished task: Basic calibration is done.")
            except Exception:
                print_banner(
                    "!!!! WARNING: Error in basic calibration. !!!!"
                )
                traceback.print_exc()
                if emails != "":
                    email_msg = f"[{cal_obsid}] Error occured in basic calibration."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
                return 1, [], []

        ##################################################################
        # Checking and interpolating bandpass tables
        ##################################################################
        print(f"Searching for bandpass tables:\n{basic_caldir}/calibrator_{cal_obsid}*.bcal")
        bandpass_tables = sorted(glob.glob(f"{basic_caldir}/calibrator_{cal_obsid}*.bcal"))
        if len(bandpass_tables) == 0:
            print(
                f"No bandpass table is present.\nCalibration directory : {basic_caldir}."
            )
            if emails != "":
                email_msg = f"[{cal_obsid}] No bandpass calibration table is found."
                send_task_notification(emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}")
            return 1, [], []
        bandpass_tables = interpolate_bpass(bandpass_tables, overwrite=True)
        
        ################################
        # Scale bandpass for attenuators
        ################################
        calibrator_header = fits.getheader(cal_metafits)
        cal_attn = calibrator_header["ATTEN_DB"]
        target_header = fits.getheader(target_metafits)
        target_attn = target_header["ATTEN_DB"]
        for bpass_table in bandpass_tables:
            print(f"Scaling for attenuation: {bpass_table}")
            scale_bandpass(bpass_table, cal_attn, target_attn)
            
        print_banner(
            f"Bandpass tables in calibration directory:\n{basic_caldir}"
        )
        for bpass in bandpass_tables:
            print(f"{os.path.basename(bpass)}")
            
        ######################################
        # Checking crossphase tables
        ######################################
        print(
            f"Searching for crossphase tables:\n{basic_caldir}/calibrator_{cal_obsid}*.kcrossscal"
        )
        crossphase_tables = sorted(
            glob.glob(f"{basic_caldir}/calibrator_{cal_obsid}*.kcrosscal")
        )
        if len(crossphase_tables) > 0:
            crossphase_tables = interpolate_bpass(crossphase_tables, overwrite=True)
            print_banner(
                f"Crosshand phase tables in calibration directory: \n{basic_caldir}"
            )
            for kcross in crossphase_tables:
                print(f"{os.path.basename(kcross)}")

        ###############################################
        # Making diagnostic plots
        ###############################################
        if len(bandpass_tables) > 0 and do_basic_cal:
            os.makedirs(f"{cal_outdir}/diagnostic_plots", exist_ok=True)
            msg, bpass_plots = plot_caltable_diagnostics(
                bandpass_tables,
                f"{cal_outdir}/diagnostic_plots/{cal_obsid}_bcal",
            )
            if msg == 0:
                print_banner(
                    f"Diagnostic plots for bandpass tables are saved in:\n{bpass_plots}."
                )
            else:
                print(
                    "Error in creating diagnostic plots for bandpass tables."
                )
        if len(crossphase_tables) > 0 and do_basic_cal:
            os.makedirs(f"{cal_outdir}/diagnostic_plots", exist_ok=True)
            msg, kcross_plots = plot_caltable_diagnostics(
                crossphase_tables,
                f"{cal_outdir}/diagnostic_plots/{cal_obsid}_kcrosscal",
                quantities=["phase"],
                plot_all_ants=False,
            )
            if msg == 0:
                print_banner(
                    f"Diagnostic plots for crosshand phase tables are saved in:\n{kcross_plots}."
                )
            else:
                print(
                    "Error in creating diagnostic plots for crosshand phase tables."
                )
        return 0, bandpass_tables, crossphase_tables
    except Exception:
        traceback.print_exc()
        return 1, [], []
    finally:
        stop_event.set()
        log_thread_flow.join(timeout=5)
        

########################################################
# Self-calibration subflows
########################################################        
@flow(
    name="Self-calibration",
    description="Perform self-calibration on target measurement sets",
    log_prints=True,
)
def selfcal_subflow(
    # Core observational inputs
    target_mslist,
    target_metafits,
    target_obsid,
    # I/O and workspace
    workdir,
    basic_caldir,
    selfcaldir,
    target_outdir,
    # Processing controls
    redo_selfcal,
    do_selfcal,
    has_cal,
    solar_selfcal,
    do_sidereal_cor,
    use_solarflagger,
    keep_backup,
    # Selfcal parameters
    solint,
    timeavg, 
    freqavg,
    image_timeres,
    image_freqres,
    quack_timestamps,
    only_amplitude,
    do_ap_selfcal,
    do_polcal,
    uvrange,
    # Resource management
    cpu_frac,
    mem_frac,
    # Logging / metadata
    jobid,
    timestamp,
    emails,
    remote_logger,
):
    """
    Self-calibration subflow
    
    Returns
    -------
    int
        Flow success message
    list
        Self-calibration gaincal tables
    list
        Self-calibration bandpass tables
    list
        Self-calibration polcal leakage tables
    """
    logdir = f"{workdir}/logs"
    os.makedirs(logdir, exist_ok=True)
    pre_process_logfile = f"{logdir}/selfcal_subflow_{target_obsid}.log"
    ctx = get_run_context()
    flow_id = str(ctx.flow_run.id)
    flow_name = ctx.flow_run.name
    stop_event = Event()
    log_thread_flow = start_flow_log_saver(
        flow_id, flow_name, pre_process_logfile, poll_interval=3, stop_event=stop_event
    )
    try:
        ###################################################
        # Checking if selfcal tables already exist or not
        ###################################################
        if not redo_selfcal:
            print("Checking pre-existing self-calibration solutions.")
            selfcal_gaincal = sorted(
                glob.glob(f"{selfcaldir}/selfcal_{target_obsid}*.gcal")
            )
            if len(selfcal_gaincal) > 0:
                selfcal_bandpass = sorted(
                    glob.glob(f"{selfcaldir}/selfcal_{target_obsid}*.bcal")
                )
                if len(selfcal_bandpass) > 0:
                    selfcal_bandpass = interpolate_bpass(
                        selfcal_bandpass, overwrite=True
                    )
                if do_polcal:
                    selfcal_leakage = sorted(
                        glob.glob(f"{selfcaldir}/selfcal_{target_obsid}*.dcal")
                    )
                    if len(selfcal_leakage) > 0:
                        selfcal_leakage = interpolate_quartical(
                            selfcal_leakage, overwrite=True
                        )
                        print_banner(
                            "Self-calibration solutions exist including polarisation calibration.\nNot performing self-calibration"
                        )
                        if emails != "":
                            email_msg = f"[{target_obsid}] Self-calibration solutions including polarisation for target are already present."
                            send_task_notification(
                                emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                            )
                        return 0, selfcal_gaincal, selfcal_bandpass, selfcal_leakage
                    else:
                        print_banner(
                            "Self-calibration solutions exist without polarisation calibration.\nHence, performing self-calibration"
                        )
                else:
                    print_banner(
                        "Self-calibration solutions exist without polarisation calibration.\nPolarisation calibration is not requested."
                    )
                    return 0, selfcal_gaincal, selfcal_bandpass, []

        ###################################################
        # Start spliting selfcal ms
        ###################################################
        if not do_selfcal:
            print_banner("Self-calibration is not requested and previous self-calibration tables are also not present.")
            if emails != "":
                email_msg = f"[{target_obsid}] Self-calibration is not requested and previous self-calibration tables are also not present."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                )
            return 1, [], [], []
        else:
            ###############################################
            # Removing previous self-calibration artificats
            ###############################################
            print("Removing all previous self-calibration artificats.")
            os.system(
                f"rm -rf {workdir}/selfcal* {workdir}/.intselfcal* {workdir}/.polselfcal*"
            )
            prefix = "selfcal"
            try:
                time_interval = float(solint)
            except BaseException:
                if solint.endswith("s"):
                    time_interval = float(solint.split("s")[0])
                elif solint.endswith("min"):
                    time_interval = float(solint.split("min")[0]) * 60
                elif solint == "int":
                    time_interval = image_timeres
                else:
                    time_interval = 60.0

            ######################
            # Spliting
            ######################
            if emails != "":
                email_msg = f"[{target_obsid}] Started spliting of measurement sets for self-calibration."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                )
            print_banner(f"Starting task: Spliting {prefix}.")
            ntime = get_selfcal_ntimes(target_mslist[0])
            msmd = msmetadata()
            msmd.open(target_mslist[0])
            times = msmd.timesforspws(0)
            timeres = np.nanmean(np.diff(times))
            msmd.close()
            time_window = min(10, round(ntime * timeres, 1))  # Maximum 10s
            future_selfcal_split = run_target_split_jobs.with_options(
                task_run_name=f"split_{target_obsid}"
            ).submit(
                ",".join(target_mslist),
                target_metafits,
                workdir,
                datacolumn="data",
                timeres=timeavg,
                freqres=freqavg,
                prefix=prefix,
                force_split=True,
                only_disk=True,
                time_window=min(time_window, time_interval),
                time_interval=time_interval,
                quack_timestamps=quack_timestamps,
                jobid=jobid,
                cpu_frac=float(cpu_frac),
                mem_frac=float(mem_frac),
                remote_log=remote_logger,
            )
            try:
                msg, expected, succeed = future_selfcal_split.result()
                if emails != "":
                    email_msg = f"[{target_obsid}] Spliting of measurement sets for self-calibration is done.\nExpected: {expected}, succeeded: {succeed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
                print_banner(
                    "Finished task: Spliting of measurement sets for self-calibration is done."
                )
            except Exception:
                print_banner(
                    "!!!! WARNING: Error in running spliting target scans for selfcal. !!!!"
                )
                traceback.print_exc()
                if emails != "":
                    email_msg = f"[{target_obsid}] Error occured in spliting target measurement sets for self-calibration."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
                return 1, [], [], []

            ######################################
            # Checking status of self-cal split
            ######################################
            print("Checking measurement sets before spawning self-calibrations.")
            ####################################
            # Filtering any corrupted ms
            #####################################
            selfcal_target_mslist = sorted(glob.glob(workdir + "/selfcal*_ch_*.ms"))
            if (selfcal_target_mslist) == 0:
                print_banner(
                    "!!!! WARNING: Error in running spliting target scans for selfcal. !!!!"
                )
                if emails != "":
                    email_msg = f"[{target_obsid}] No splited measurement set is found for self-calibration.\nNot continuting for self-calibration."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
                return 1, [], [], []

            filtered_mslist = []  # Filtering in case any ms is corrupted
            for ms in selfcal_target_mslist:
                checkcol = check_datacolumn_valid(ms)
                if checkcol:
                    filtered_mslist.append(ms)
                else:
                    print(f"Issue in : {ms}")
                    os.system(f"rm -rf {ms}")
            selfcal_mslist = filtered_mslist
            if len(selfcal_mslist) == 0:
                print_banner(
                    "No splited target scan ms are available in work directory for selfcal.\nNot continuing further for selfcal."
                )
                if emails != "":
                    email_msg = f"[{target_obsid}] No splited measurement set is found for self-calibration. Not continuting for self-calibration."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
                return 1, [], [], []
                
            print_banner("Selfcal measurement set list:")
            for ms in [os.path.basename(i) for i in selfcal_mslist]:
                print(ms)

            #########################################################
            # Flagging on targets for self-calibration
            #########################################################
            cal_applied = False
            ###################################
            # Apply basic calibration
            ###################################
            if has_cal:
                if emails != "":
                    email_msg = f"[{target_obsid}] Started applying basic calibration solution on self-calibration measurement sets."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
                print_banner(
                    "Starting task: Applying basic calibration on self-calibration measurement sets."
                )
                future_apply_basical_selfcal = run_apply_basiccal_sol.with_options(
                    task_run_name=f"apply_basic_cal_{target_obsid}"
                ).submit(
                    ",".join(selfcal_mslist),
                    target_metafits,
                    workdir,
                    basic_caldir,
                    overwrite_datacolumn=False,
                    only_amplitude=only_amplitude,
                    applymode="calflag",
                    prefix="selfcal",
                    jobid=jobid,
                    cpu_frac=round(cpu_frac, 2),
                    mem_frac=round(mem_frac, 2),
                    remote_log=remote_logger,
                )
                try:
                    msg, succeed, failed = future_apply_basical_selfcal.result()
                    cal_applied = True
                    if emails != "":
                        email_msg = f"[{target_obsid}] Applying basic calibration solution on self-calibration measurement sets are done.\nSucceeded: {succeed}, failed: {failed}."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                        )
                    print_banner(
                        "Finished task: Applying basic calibration solution on self-calibration measurement sets are done."
                    )
                except Exception:
                    print_banner(
                        "!!!! WARNING: Error in applying basic calibration solutions on target.\nContinuing selfcal without basic calibration.!!!!"
                    )
                    traceback.print_exc()
                    if emails != "":
                        email_msg = f"[{target_obsid}] Error occured in applying basic calibration solutions on self-calibration measurement sets."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                        )
                        
            ########################################
            # Filtering out for self-calibration
            ########################################
            if cal_applied:
                selfcal_applymode = "calonly"
                filtered_selfcalms_list = []
                for selfcalms in selfcal_mslist:
                    unflag_chans, flag_chans = get_chans_flag(
                        msname=selfcalms, n_threads=1
                    )
                    if len(flag_chans) / (len(flag_chans) + len(unflag_chans)) <= 0.8:
                        filtered_selfcalms_list.append(selfcalms)
                    else:
                        print(
                            f"More than 80% channels are flagged for ms: {selfcalms}. Not using for self-calibration."
                        )
                if len(filtered_selfcalms_list)==0:
                    print_banner("No measurement set is present with unflagged data for self-calibration after applying basic-calibration.")
                    if emails != "":
                        email_msg = f"[{target_obsid}] No measurement set is present with unflagged data for self-calibration after applying basic-calibration."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                        )
                    return 1, [], [], []
                else:
                    selfcal_mslist = filtered_selfcalms_list
            else:
                selfcal_applymode = "calflag"
            
            ###############################################
            # Performing sidereal correction before selfcal
            ###############################################
            os.system(
                f"rm -rf {workdir}/*selfcal_int* {workdir}/*selfcal_pol* {workdir}/caltables/*selfcal*"
            )
            if do_sidereal_cor:
                if emails != "":
                    email_msg = f"[{target_obsid}] Started correcting for solar sidereal motion."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
                print_banner(
                    "Starting task: Sidereal motion correction for self-calibration measurement sets."
                )
                future_sidereal_cor_selfcal = run_solar_siderealcor_jobs.with_options(
                    task_run_name=f"sidereal_cor_{target_obsid}"
                ).submit(
                    ",".join(selfcal_mslist),
                    workdir,
                    prefix="selfcal",
                    jobid=jobid,
                    cpu_frac=round(cpu_frac, 2),
                    mem_frac=round(mem_frac, 2),
                    remote_log=remote_logger,
                )
                try:
                    msg, succeed, failed = future_sidereal_cor_selfcal.result()
                    if emails != "":
                        email_msg = f"[{target_obsid}] Correction for solar sidereal motion is done.\nSucceeded: {succeed}, failed: {failed}."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                        )
                    print_banner(
                        "Finished task: Correction for solar sidereal motion is done."
                    )
                except Exception:
                    print_banner("!!! WARNING : Sidereal correction is not successful. !!!")
                    traceback.print_exc()
                    if emails != "":
                        email_msg = f"[{target_obsid}] Error occured in sidereal motion correction."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                        )

            ############################
            # Basic flagging for selfcal
            ############################
            if emails != "":
                email_msg = f"[{target_obsid}] Started flagging for self-calibration measurment sets."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                )
            print_banner("Starting task: Flagging selfcal targets.")
            future_flag = run_flag.with_options(
                task_run_name=f"flag_{target_obsid}"
            ).submit(
                ",".join(selfcal_mslist),
                target_metafits,
                workdir,
                target_outdir,
                flag_calibrators=False,
                flag_quack=False,
                datacolumn="corrected",
                run_solarflagger=use_solarflagger,
                jobid=jobid,
                cpu_frac=round(cpu_frac, 2),
                mem_frac=round(mem_frac, 2),
                remote_log=remote_logger,
            )
            try:
                msg, succeed, failed = future_flag.result()
                if emails != "":
                    email_msg = f"[{target_obsid}] Flagging for self-calibration measurment sets are done.\nSucceeded: {succeed}, failed: {failed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
                for s_ms in selfcal_mslist:
                    s_ms = s_ms.rstrip("/")
                    if os.path.exists(f"{s_ms}/.flag_failed"):
                        print(
                            f"Issue in flagging: {s_ms}. Check calibration solutions carefully."
                        )
                print_banner(
                    "Finished task: Flagging for self-calibration measurment sets are done."
                )
            except Exception:
                print_banner(
                    "!!!! WARNING: Flagging error. Examine calibration solutions with caution. !!!!"
                )
                traceback.print_exc()
                if emails != "":
                    email_msg = (
                        f"[{target_obsid}] Error occured in flagging self-calibration measurement sets."
                    )
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )

            #############################
            # Self-calibration
            #############################
            if emails != "":
                email_msg = f"[{target_obsid}] Started self-calibration."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                )
            print_banner("Starting task: Self-calibrations.")
            if cal_applied:
                print("Calibrator solutions are applied.")
            else:
                print("Calibration solutions are not applied")
            future_selfcal = run_selfcal_jobs.with_options(
                task_run_name=f"selfcal_{jobid}"
            ).submit(
                ",".join(selfcal_mslist),
                workdir,
                selfcaldir,
                target_metafits,
                cal_applied,
                solint=solint,
                do_apcal=do_ap_selfcal,
                do_polcal=do_polcal,
                solar_selfcal=solar_selfcal,
                keep_backup=keep_backup,
                uvrange=uvrange,
                weight="briggs",
                robust=0.0,
                applymode=selfcal_applymode,
                use_solarflagger=use_solarflagger,
                jobid=jobid,
                cpu_frac=round(cpu_frac, 2),
                mem_frac=round(mem_frac, 2),
                remote_log=remote_logger,
            )
            try:
                (
                    msg,
                    int_succeed,
                    int_failed,
                    pol_succeed,
                    pol_failed,
                    int_DR,
                    pol_DR,
                ) = future_selfcal.result()
                if emails != "":
                    email_msg = f"[{target_obsid}] Self-calibration is done.\nIntensity self-calibration, Succeeded: {int_succeed}, failed: {int_failed}, average DR: {int_DR}."
                    if do_polcal:
                        email_msg += f"\nPolarisation self-calibration, Succeeded: {pol_succeed}, failed: {pol_failed}, average DR; {pol_DR}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
                print_banner("Finished task: Self-calibration is done.")
            except Exception:
                print_banner(
                    "!!!! WARNING: Error in self-calibration on targets. Not applying self-calibration. !!!!"
                )
                traceback.print_exc()
                if emails != "":
                    email_msg = f"[{target_obsid}] Error occured in self-calibration."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
                return 1, [], [], []

            ########################################
            # Checking self-cal caltables
            ########################################
            print(
                f"Searching for self-calibration gaincal tables:\n{selfcaldir}/selfcal_{target_obsid}*.gcal"
            )
            selfcal_gaincal = sorted(
                glob.glob(f"{selfcaldir}/selfcal_{target_obsid}*.gcal")
            )
            if len(selfcal_gaincal) == 0:
                print_banner(
                    "Self-calibration is not performed and no self-calibration caltable is available."
                )
                if emails != "":
                    email_msg = f"[{target_obsid}] Self-calibration is not performed and no self-calibration caltable is available."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp, flow_name=f"Subflow {flow_name}",
                    )
            else:
                print_banner(
                    f"Self-calibration gaincal tables in calibration directory:\n{selfcaldir}"
                )
                for gcal in selfcal_gaincal:
                    print(f"{os.path.basename(gcal)}")
                print(
                    f"Searching for self-calibration bandpass tables:\n{selfcaldir}/selfcal_{target_obsid}*.bcal"
                )
                selfcal_bandpass = sorted(
                    glob.glob(f"{selfcaldir}/selfcal_{target_obsid}*.bcal")
                )
                if len(selfcal_bandpass) > 0:
                    print_banner(
                        f"Self-calibration bandpass tables in calibration directory:\n{selfcaldir}"
                    )
                    for bpass in selfcal_bandpass:
                        print(f"{os.path.basename(bpass)}")
                    selfcal_bandpass = interpolate_bpass(selfcal_bandpass, overwrite=True)
                if do_polcal:
                    print(
                        f"Searching for self-calibration polarisation leakage tables:\n{selfcaldir}/selfcal_{target_obsid}*.dcal"
                    )
                    selfcal_leakage = sorted(
                        glob.glob(f"{selfcaldir}/selfcal_{target_obsid}*.dcal")
                    )
                    if len(selfcal_leakage) > 0:
                        print_banner(
                            f"Self-calibration polarisation leakage tables in calibration directory:\n{selfcaldir}"
                        )
                        for dcal in selfcal_leakage:
                            print(f"{os.path.basename(dcal)}")
                        selfcal_leakage = interpolate_quartical(
                            selfcal_leakage, overwrite=True
                        )

                ###########################################
                # Plotting self-caltables
                ###########################################
                if len(selfcal_gaincal) > 0:
                    os.makedirs(f"{target_outdir}/diagnostic_plots", exist_ok=True)
                    msg, gcal_plots = plot_caltable_diagnostics(
                        selfcal_gaincal,
                        f"{target_outdir}/diagnostic_plots/{target_obsid}_gcal",
                    )
                    if msg == 0:
                        print(
                            f"Diagnostic plots for self-calibration gaincal tables are saved in:\n{gcal_plots}."
                        )
                    else:
                        print(
                            "Error in creating diagnostic plots for self-calibration gaincal tables."
                        )

                if len(selfcal_bandpass) > 0:
                    os.makedirs(f"{target_outdir}/diagnostic_plots", exist_ok=True)
                    msg, bcal_plots = plot_caltable_diagnostics(
                        selfcal_bandpass,
                        f"{target_outdir}/diagnostic_plots/{target_obsid}_bcal",
                    )
                    if msg == 0:
                        print(
                            f"Diagnostic plots for self-calibration bandpass tables are saved in:\n{bcal_plots}."
                        )
                    else:
                        print(
                            "Error in creating diagnostic plots for self-calibration bandpass tables."
                        )

                if do_polcal:
                    if len(selfcal_leakage) > 0:
                        os.makedirs(f"{target_outdir}/diagnostic_plots", exist_ok=True)
                        msg, dcal_plots = plot_quartical_tables(
                            selfcal_leakage,
                            f"{target_outdir}/diagnostic_plots/{target_obsid}_dcal",
                        )
                        if msg == 0:
                            print(
                                f"Diagnostic plots for self-calibration leakage tables are saved in: \n{dcal_plots}."
                            )
                        else:
                            print(
                                "Error in creating diagnostic plots for self-calibration leakage tables."
                            )
                return 0, selfcal_gaincal, selfcal_bandpass, selfcal_leakage
    except Exception:
        traceback.print_exc()
        return 1, [], [], []
    finally:
        stop_event.set()
        log_thread_flow.join(timeout=5)
            
                 

