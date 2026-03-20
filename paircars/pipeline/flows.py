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
    send_task_notification,
)
from prefect.context import get_run_context
from multiprocessing import Event
from paircars.utils.prefect_logger_utils import (
    start_flow_log_saver,
)

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
    pre_process_logfile = f"{logdir}/pre_process_{target_obsid}.log"
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
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print("Starting task: Moving phasecenter to the Sun.")
            print("###########################")
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
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print("Finished task: Moving phasecenter to solar center is done.")
                print("###########################")
                filtered_ms = []
                for t_ms in target_mslist:
                    t_ms = t_ms.rstrip("/")
                    if os.path.exists(f"{t_ms}/.solarcenter_move_succeed"):
                        filtered_ms.append(t_ms)
                    else:
                        print(
                            f"Issue in moving phasecneter to solar center for measurement set: {t_ms}"
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
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                return 1, []
        #######################################
        # Run dynamic spectra making
        #######################################
        if solar_data and make_ds:
            if emails != "":
                email_msg = f"[{target_obsid}] Started making solar dynamic spectra."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print("Starting task: Making dynamic spectra of solar target.")
            print("###########################")
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
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print("Finished task: Making solar dynamic spectra are done.")
                print("###########################")
            except Exception:
                print("!!! WARNING : Error in making dynamic spectra. !!!")
                traceback.print_exc()
                if emails != "":
                    email_msg = f"[{target_obsid}] Error occured in making dynamic spectra."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
        
        return 0, target_mslist 
    except Exception:
        traceback.print_exc()
        return 1, []
    finally:
        stop_event.set()
        log_thread_flow.join(timeout=5)


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
    timestamp,
    emails,
    remote_logger,
):
    """
    Basic calibration sub flow
    """
    logdir = f"{workdir}/logs"
    os.makedirs(logdir, exist_ok=True)
    basic_cal_logfile = f"{logdir}/basic_cal_{cal_obsid}.log"
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
                    if emails != "":
                        email_msg = f"[{cal_obsid}] All gain solutions from calibrator are already present."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    return 0, bandpass_tables, crossphase_tables

        ############################
        # Calibrator ms list
        ############################
        cal_mslist = glob.glob(f"{cal_datadir}/*.ms")
        if len(cal_mslist) == 0 or len(coarse_chans) == 0:
            print(
                f"No calibrator measurement set with coarse channels: {coarse_chans} is present in: {cal_datadir}"
            )
            if emails != "":
                email_msg = f"[{cal_obsid}] No calibrator measurement set with coarse channels: {coarse_chans} is present in: {cal_datadir}."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
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
                send_task_notification(emails, email_msg, jobid, target_obsid, timestamp)
            print("###########################")
            print(
                "Starting task: Spliting of calibrator measurement sets."
            )
            print("###########################")
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
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(
                    "Finished task: Spliting of calibrator measurement sets are done."
                )
                print("###########################")
            except Exception:
                print(
                    "!!!! WARNING: Error in spliting calibrator measurement sets. !!!!"
                )
                traceback.print_exc()
                if emails != "":
                    email_msg = f"[{cal_obsid}] Spliting calibrator measurement set is failed."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                return 1, [], []

        if do_cal_flag or do_import_model or do_basic_cal:
            split_cal_mslist = sorted(
                glob.glob(f"{workdir}/calibrator_{cal_obsid}*spw_*.ms")
            )
            if len(split_cal_mslist) == 0:
                print(
                    "No splited measurement set is present for basic calibration."
                )
                if emails != "":
                    email_msg = f"[{cal_obsid}] No splited measurement set is present for basic calibration."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                return 1, [], []

        ##################################
        # Run flagging jobs on calibrators
        ##################################
        # Only if basic calibration is requested
        if do_cal_flag:
            if emails != "":
                email_msg = f"[{cal_obsid}] Started flagging of calibrators."
                send_task_notification(emails, email_msg, jobid, target_obsid, timestamp)
            print("###########################")
            print("Starting task: Flagging calibrators.")
            print("###########################")
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
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                filtered_ms = []
                for c_ms in split_cal_mslist:
                    c_ms = c_ms.rstrip("/")
                    if os.path.exists(f"{c_ms}/.flag_succeed"):
                        filtered_ms.append(c_ms)
                    else:
                        print(f"Issue in flagging of measurement set: {c_ms}")
                split_cal_mslist = filtered_ms  # Filtered target mslist
                print("###########################")
                print(
                    "Finished task: Flagging of calibrator is done."
                )
                print("###########################")
            except Exception:
                print(
                    "!!!! WARNING: Flagging error for calibrator. !!!!"
                )
                traceback.print_exc()
                if emails != "":
                    email_msg = f"[{cal_obsid}] Error in flagging calibrators."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )

        #################################
        # Import model
        #################################
        if do_import_model:
            if emails != "":
                email_msg = (
                    f"[{cal_obsid}] Started importing sky model."
                )
                send_task_notification(emails, email_msg, jobid, target_obsid, timestamp)
            print("###########################")
            print("Starting task: Importing model visibilities.")
            print("###########################")
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
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(
                    "Finished task: Model import for calibrator is done."
                )
                print("###########################")
                filtered_ms = []
                for c_ms in split_cal_mslist:
                    c_ms = c_ms.rstrip("/")
                    if os.path.exists(f"{c_ms}/.modeling_succeed"):
                        filtered_ms.append(c_ms)
                    else:
                        print(
                            f"Issue in importing calibrator sky model of measurement set: {c_ms}"
                        )
                split_cal_mslist = filtered_ms  # Filtered target mslist
            except Exception:
                print(
                    "!!!! WARNING: Error in importing calibrator models. Not continuing calibration. !!!!"
                )
                traceback.print_exc()
                if emails != "":
                    email_msg = f"[{cal_obsid}] Error occured in importing model for calibrators. Not using calibrator solutions."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                return 1, [], []

        ###############################
        # Run basic calibration
        ###############################
        if do_basic_cal:
            if emails != "":
                email_msg = f"[{cal_obsid}] Started basic calibration."
                send_task_notification(emails, email_msg, jobid, target_obsid, timestamp)
            print("###########################")
            print("Starting task: Performing basic calibration.")
            print("###########################")
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
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print("Finished task: Basic calibration is done.")
                print("###########################")
            except Exception:
                print(
                    "!!!! WARNING: Error in basic calibration. !!!!"
                )
                traceback.print_exc()
                if emails != "":
                    email_msg = f"[{cal_obsid}] Error occured in basic calibration."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                return 1, [], []

        ##################################################################
        # Checking presence of necessary caltables if not checked already
        #################################################################
        print(f"Searching for bandpass tables: {basic_caldir}/calibrator_{cal_obsid}*.bcal")
        bandpass_tables = sorted(glob.glob(f"{basic_caldir}/calibrator_{cal_obsid}*.bcal"))
        if len(bandpass_tables) > 0:
            bandpass_tables = interpolate_bpass(bandpass_tables, overwrite=True)
        print(
            f"Searching for crossphase tables: {basic_caldir}/calibrator_{cal_obsid}*.kcrossscal"
        )
        crossphase_tables = sorted(
            glob.glob(f"{basic_caldir}/calibrator_{cal_obsid}*.kcrosscal")
        )
        if len(crossphase_tables) > 0:
            crossphase_tables = interpolate_bpass(crossphase_tables, overwrite=True)
        if len(bandpass_tables) == 0:
            print(
                f"No bandpass table is present in calibration directory : {basic_caldir}."
            )
            if emails != "":
                email_msg = f"[{cal_obsid}] No bandpass calibration table is found."
                send_task_notification(emails, email_msg, jobid, target_obsid, timestamp)
            return 1, [], []
        else:
            print("###################################################")
            print(
                f"Bandpass tables in calibration directory: {basic_caldir}"
            )
            for bpass in bandpass_tables:
                print(f"{os.path.basename(bpass)}")
            print("####################################################")
            if len(crossphase_tables) > 0:
                print(
                    f"Crosshand phase tables in calibration directory: {basic_caldir}"
                )
                for kcross in crossphase_tables:
                    print(f"{os.path.basename(kcross)}")
                print("####################################################")

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
                print(
                    f"Diagnostic plots for bandpass tables are saved in : {bpass_plots}."
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
                print(
                    f"Diagnostic plots for crosshand phase tables are saved in : {kcross_plots}."
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
