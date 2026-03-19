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
    send_task_notification,
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

    ##############################
    # Run spliting jobs
    ##############################
    # If basic calibration is requested and calibrator ms and metafits are present
    if do_cal_flag or do_import_model or do_basic_cal:
        prefix = "calibrator"
        #if emails != "":
        #    email_msg = f"Started spliting of calibrator measurement sets for OBSID {cal_obsid}."
        #    send_task_notification(emails, email_msg, jobid, target_obsid, timestamp)
        print("###########################")
        print(
            f"Starting task: Spliting of calibrator measurement sets for OBSID {cal_obsid}......"
        )
        print("###########################")
        future_cal_split = run_target_split_jobs.with_options(
            task_run_name=f"spliting_{prefix}_{cal_obsid}_{jobid}"
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
            #if emails != "":
            #    email_msg = f"Spliting of calibrator measurement sets for OBSID: {cal_obsid} are done.\nExpected: {expected}, succeeded: {succeed}."
            #    send_task_notification(
            #        emails, email_msg, jobid, target_obsid, timestamp
            #    )
            print("###########################")
            print(
                f"Finished task: Spliting of calibrator measurement sets for OBSID {cal_obsid} are done."
            )
            print("###########################")
        except Exception:
            print(
                f"!!!! WARNING: Error in spliting calibrator measurement sets for OBSID {cal_obsid}. !!!!"
            )
            traceback.print_exc()
            #if emails != "":
            #    email_msg = f"Spliting calibrator measurement set for OBSID {cal_obsid} is failed."
            #   send_task_notification(
             #       emails, email_msg, jobid, target_obsid, timestamp
             #   )
            return 1, [], []

    if do_cal_flag or do_import_model or do_basic_cal:
        split_cal_mslist = sorted(
            glob.glob(f"{workdir}/calibrator_{cal_obsid}*spw_*.ms")
        )
        if len(split_cal_mslist) == 0:
            print(
                f"No splited measurement set is present for basic calibration for OBSID {cal_obsid}."
            )
            return 1, [], []

    ##################################
    # Run flagging jobs on calibrators
    ##################################
    # Only if basic calibration is requested
    if do_cal_flag:
        #if emails != "":
        #    email_msg = f"Started flagging of calibrators for OBSID {cal_obsid}."
        #    send_task_notification(emails, email_msg, jobid, target_obsid, timestamp)
        print("###########################")
        print(f"Starting task: Flagging calibrators for OBSID {cal_obsid}....")
        print("###########################")
        future_flag = run_flag.with_options(
            task_run_name=f"flagging_cal_{cal_obsid}_{jobid}"
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
            #if emails != "":
             #   email_msg = f"Flagging of calibrator for OBSID {cal_obsid} is done.\nSucceeded: {succeed}, failed: {failed}."
             #   send_task_notification(
             #       emails, email_msg, jobid, target_obsid, timestamp
              #  )
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
                f"Finished task: Flagging of calibrator is done for OBSID {cal_obsid}."
            )
            print("###########################")
        except Exception:
            print(
                f"!!!! WARNING: Flagging error for calibrator with OBSID {cal_obsid}. !!!!"
            )
            traceback.print_exc()
            '''if emails != "":
                email_msg = f"Error in flagging calibrators for OBSID {cal_obsid}."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )'''
            return 1, [], []

    #################################
    # Import model
    #################################
    if do_import_model:
        '''if emails != "":
            email_msg = (
                f"Started importing sky model for calibrator for OBSID {cal_obsid}."
            )
            send_task_notification(emails, email_msg, jobid, target_obsid, timestamp)'''
        print("###########################")
        print(f"Starting task: Importing model visibilities for OBSID {cal_obsid}....")
        print("###########################")
        future_import_model = run_import_model.with_options(
            task_run_name=f"importing_model_visibilities_{cal_obsid}_{jobid}"
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
            '''if emails != "":
                email_msg = f"Model import for calibrator for OBSID {cal_obsid} is done.\nSucceeded: {succeed}, failed: {failed}."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )'''
            print("###########################")
            print(
                f"Finished task: Model import for calibrator for OBSID {cal_obsid} is done."
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
                f"!!!! WARNING: Error in importing calibrator models for OBSID {cal_obsid}. Not continuing calibration. !!!!"
            )
            traceback.print_exc()
            '''if emails != "":
                email_msg = f"Error occured in importing model for calibrators for OBSID {cal_obsid}. Not using calibrator solutions."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )'''
            return 1, [], []

    ###############################
    # Run basic calibration
    ###############################
    if do_basic_cal:
        '''if emails != "":
            email_msg = f"Started basic calibration for OBSID {cal_obsid}."
            send_task_notification(emails, email_msg, jobid, target_obsid, timestamp)'''
        print("###########################")
        print(f"Starting task: Performing basic calibration for OBSID {cal_obsid}.....")
        print("###########################")
        future_basical = run_basic_cal_jobs.with_options(
            task_run_name=f"basic_calibration_{cal_obsid}_{jobid}"
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
            '''if emails != "":
                email_msg = f"Basic calibration is done for OBSID {cal_obsid}.\nSucceeded: {succeed}, failed: {failed}."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )'''
            print("###########################")
            print(f"Finished task: Basic calibration is done for OBSID {cal_obsid}.")
            print("###########################")
        except Exception:
            print(
                f"!!!! WARNING: Error in basic calibration for OBSID {cal_obsid}. !!!!"
            )
            traceback.print_exc()
            '''if emails != "":
                email_msg = f"Error occured in basic calibration for OBSID {cal_obsid}."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )'''
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
            f"No bandpass table is present for OBSID {cal_obsid} in calibration directory : {basic_caldir}."
        )
        '''if emails != "":
            email_msg = f"No bandpass calibration table is found for OBSID {cal_obsid}."
            send_task_notification(emails, email_msg, jobid, target_obsid, timestamp)'''
        return 1, [], []
    else:
        print("###################################################")
        print(
            f"Bandpass tables for OBSID {cal_obsid} in calibration directory: {basic_caldir}"
        )
        for bpass in bandpass_tables:
            print(f"{os.path.basename(bpass)}")
        print("####################################################")
        if len(crossphase_tables) > 0:
            print(
                f"Crosshand phase tables for OBSID {cal_obsid} in calibration directory: {basic_caldir}"
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
                f"Diagnostic plots for bandpass tables for OBSID {cal_obsid} are saved in : {bpass_plots}."
            )
        else:
            print(
                f"Error in creating diagnostic plots for bandpass tables for OBSID {cal_obsid}."
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
                f"Diagnostic plots for crosshand phase tables for OBSID {cal_obsid} are saved in : {kcross_plots}."
            )
        else:
            print(
                f"Error in creating diagnostic plots for crosshand phase tables for OBSID {cal_obsid}."
            )
    return 0, bandpass_tables, crossphase_tables
