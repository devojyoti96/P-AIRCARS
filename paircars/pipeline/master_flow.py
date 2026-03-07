import os
import psutil
import numpy as np
import argparse
import traceback
import copy
import time
import glob
import sys
import os
import socket
import requests
import getpass
from collections import Counter
from casatools import msmetadata
from astropy.io import fits
from datetime import datetime as dt
from multiprocessing import Process, Event
from dask.distributed import get_client
from dotenv import load_dotenv
from pyfiglet import Figlet
from prefect import flow, task
from prefect.context import get_run_context
from prefect_dask.task_runners import DaskTaskRunner
from prefect_dask import get_dask_client
from prefect.settings import get_current_settings
from paircars.utils.basic_utils import (
    get_cachedir,
    timestamp_to_mjdsec,
)
from paircars.utils.calibration import (
    calc_bw_smearing_freqwidth,
    calc_time_smearing_timewidth,
    max_time_solar_smearing,
    interpolate_bpass,
    interpolate_quartical,
)
from paircars.utils.casatasks import reset_weights_and_flags
from paircars.utils.flagging import do_flag_backup, get_chans_flag
from paircars.utils.image_utils import filter_images
from paircars.utils.logger_utils import (
    SmartDefaultsHelpFormatter,
    clean_shutdown,
    generate_password,
    get_remote_logger_link,
    get_emails,
    init_logger,
)
from paircars.utils.ms_metadata import get_ms_size, check_datacolumn_valid
from paircars.utils.mwa_ploting_utils import (
    plot_caltable_diagnostics,
    plot_quartical_tables,
)
from paircars.utils.mwa_utils import (
    get_ncoarse,
    get_MWA_coarse_chan,
    get_MWA_OBSID,
    download_MWA_metafits,
)
from paircars.utils.proc_manage_utils import (
    get_jobid,
    save_main_process_info,
    get_total_worker,
    scale_worker_and_wait,
    get_total_nodes,
    get_total_nodes,
    get_local_dask_cluster,
    get_scheduler_name,
)
from paircars.utils.resource_utils import drop_cache
from paircars.data.sendmail import (
    send_paircars_notification as send_notification,
)
from paircars.clusterutils.slurm_cluster import (
    get_slurm_dask_cluster,
    get_slurm_node_resources,
    is_slurm_job,
)
from paircars.utils.prefect_logger_utils import (
    start_log_task_saver,
    start_flow_log_saver,
)
from paircars.pipeline import (
    mwa_make_ds,
    do_target_split,
    flagging,
    import_model,
    basic_cal,
    do_apply_basiccal,
    do_sidereal_cor,
    do_selfcal,
    do_apply_selfcal,
    do_imaging,
    mwa_pbcor,
    make_mwa_overlay,
    move_solarcenter,
    make_ms_plot,
)
from paircars.pipeline.init_data import init_paircars_data


@task(
    name="moving_to_solar_center",
    retries=2,
    retry_delay_seconds=60,
    log_prints=True,
)
def run_solar_phasecenter_jobs(
    mslist,
    workdir,
    prefix="target",
    jobid=0,
    cpu_frac=0.8,
    mem_frac=0.8,
    remote_log=False,
):
    """
    Move phase center to the Sun

    Parameters
    ----------
    mslist: str
        List of the measurement sets (comma separated)
    workdir : str
        Work directory
    prefix : str, optional
        Measurement set prefix
    cpu_frac : float, optional
        CPU fraction to use
    mem_frac : float, optional
        Memory fraction to use
    remote_log: bool, optional
        Start remote logger

    Returns
    -------
    int
        Success message
    int
        Succeeded ms number
    int
        Failed ms number
    """
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)
    phasecor_basename = f"cor_phasecenter_{prefix}"
    logdir = f"{workdir}/logs"
    os.makedirs(logdir, exist_ok=True)
    logfile = f"{logdir}/{phasecor_basename}.log"
    if os.path.exists(logfile):
        os.remove(logfile)
    ctx = get_run_context()
    task_id = str(ctx.task_run.id)
    task_name = ctx.task_run.name
    stop_event = Event()
    log_thread_sidereal = start_log_task_saver(
        task_id, task_name, logfile, poll_interval=3, stop_event=stop_event
    )
    try:
        #######################
        print("###########################")
        print("Moving phasecenter to the Sun .....")
        print("###########################")
        #######################
        # Moving phasecenter motion correction
        #######################
        with get_dask_client() as dask_client:
            msg, succeed, failed = move_solarcenter.main(
                mslist,
                workdir=workdir,
                cpu_frac=float(cpu_frac),
                mem_frac=float(mem_frac),
                logfile=logfile,
                jobid=jobid,
                start_remote_log=remote_log,
                dask_client=dask_client,
            )
    finally:
        stop_event.set()
        log_thread_sidereal.join(timeout=5)
    if msg != 0:
        raise RuntimeError("Moving phasecenter to solar center is failed.")
    else:
        return msg, succeed, failed


@task(
    name="making_dynamic_spectra",
    retries=2,
    retry_delay_seconds=60,
    log_prints=True,
)
def run_ds_jobs(
    mslist,
    metafits,
    workdir,
    outdir,
    plot_quantity="TB",
    jobid=0,
    cpu_frac=0.8,
    mem_frac=0.8,
    remote_log=False,
):
    """
    Make dynamic spectra of the solar target

    Parameters
    ----------
    mslist : str
        Measurement sets (comma separated)
    metafits : str
        Metafits file
    workdir : str
        Name of the work directory
    outdir : str
        Name of the output directory
    plot_quantity : str, optional
        Plot quantity (TB or flux)
    cpu_frac : float, optional
        CPU fraction to use
    mem_frac : float, optional
        Memory fraction to use
    remote_log: bool, optional
        Start remote logger

    Returns
    -------
    int
        Success message
    int
        Succeeded ms number
    int
        Failed ms number
    """
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)
    ds_basename = "ds_target"
    logdir = f"{workdir}/logs"
    os.makedirs(logdir, exist_ok=True)
    logfile = f"{logdir}/{ds_basename}.log"
    if os.path.exists(logfile):
        os.remove(logfile)
    ctx = get_run_context()
    task_id = str(ctx.task_run.id)
    task_name = ctx.task_run.name
    stop_event = Event()
    log_thread_ds = start_log_task_saver(
        task_id, task_name, logfile, poll_interval=3, stop_event=stop_event
    )
    try:
        ##################
        print("###########################")
        print("Making dynamic spectra of solar target .....")
        print("###########################")
        ##########################
        # Making dynamic spectrum
        ##########################
        with get_dask_client() as dask_client:
            msg, succeed, failed = mwa_make_ds.main(
                mslist,
                metafits,
                workdir,
                outdir,
                plot_quantity=plot_quantity,
                cpu_frac=float(cpu_frac),
                mem_frac=float(mem_frac),
                logfile=logfile,
                jobid=jobid,
                start_remote_log=remote_log,
                dask_client=dask_client,
            )
    finally:
        stop_event.set()
        log_thread_ds.join(timeout=5)
    if msg != 0:
        raise RuntimeError("Dynamic spectrum making is failed.")
    else:
        return msg, succeed, failed


@task(
    name="spliting_ms",
    retries=2,
    retry_delay_seconds=60,
    log_prints=True,
)
def run_target_split_jobs(
    mslist,
    metafits,
    workdir,
    datacolumn="data",
    timeres=-1,
    freqres=-1,
    prefix="target",
    time_window=-1,
    time_interval=-1,
    quack_timestamps=-1,
    force_split=False,
    only_disk=False,
    jobid=0,
    cpu_frac=0.8,
    mem_frac=0.8,
    remote_log=False,
):
    """
    Split measurement set

    Parameters
    ----------
    mslist: str
        Name of the measurement sets (comma separated)
    metafits : str
        Metafits file
    workdir : str
        Working directory
    datacolumn : str, optional
        Data column
    timeres : float, optional
        Time bin to average in seconds
    freqres : float, optional
        Frequency averaging in MHz
    prefix : str, optional
        Prefix of splited targets
    time_window : float, optional
        Time window in seconds
    time_interval : float, optional
        Time interval in seconds
    quack_timestamps: int, optional
        Number of timestamps to flag at the beginning and end of each scan ("quack").
    force_split : bool, optional
        Force to split
    only_disk : bool, optional
        Split only disk times
    cpu_frac : float, optional
        CPU fraction to use
    mem_frac : float, optional
        Memory fraction to use
    remote_log: bool, optional
        Start remote logger

    Returns
    -------
    int
        Success message for spliting measurement set
    int
        Expected splited ms number
    int
        Succeeded splited ms number
    """
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)
    split_basename = f"split_{prefix}"
    logdir = f"{workdir}/logs"
    os.makedirs(logdir, exist_ok=True)
    logfile = f"{logdir}/{split_basename}.log"
    if os.path.exists(logfile):
        os.remove(logfile)
    ctx = get_run_context()
    task_id = str(ctx.task_run.id)
    task_name = ctx.task_run.name
    stop_event = Event()
    log_thread_split = start_log_task_saver(
        task_id, task_name, logfile, poll_interval=3, stop_event=stop_event
    )
    try:
        ############
        print("###########################")
        print(f"Spliting {prefix} .....")
        print("###########################")
        ##################
        # Spliting ms
        ##################
        with get_dask_client() as dask_client:
            msg, expected, succeed = do_target_split.main(
                mslist,
                metafits,
                workdir=workdir,
                datacolumn=datacolumn,
                time_window=time_window,
                time_interval=time_interval,
                freqres=freqres,
                timeres=timeres,
                quack_timestamps=quack_timestamps,
                force_split=force_split,
                only_disk=only_disk,
                prefix=prefix,
                cpu_frac=float(cpu_frac),
                mem_frac=float(mem_frac),
                logfile=logfile,
                jobid=jobid,
                start_remote_log=remote_log,
                dask_client=dask_client,
            )
    finally:
        stop_event.set()
        log_thread_split.join(timeout=5)
    if msg != 0:
        raise RuntimeError("Spliting measurement set into coarse channels is failed.")
    else:
        return msg, expected, succeed


@task(
    name="flagging",
    retries=2,
    retry_delay_seconds=60,
    log_prints=True,
)
def run_flag(
    mslist,
    metafits,
    workdir,
    outdir,
    datacolumn="DATA",
    flag_calibrators=True,
    flag_quack=True,
    restore_flag=True,
    jobid=0,
    cpu_frac=0.8,
    mem_frac=0.8,
    remote_log=False,
):
    """
    Run flagging jobs

    Parameters
    ----------
    mslist: str
        Name of the measurement sets (comma separted)
    metafits : str
        Metafits file
    workdir : str
        Working directory
    outdir : str
        Output directory
    datacolumn : str, optional
        Data column
    flag_calibrators : bool, optional
        Flag calibrator fields
    flag_quack : bool, optional
        Flag quack timestamps
    restore_flag : bool, optional
        Restore flags or not
    jobid : int, optional
        Job ID
    cpu_frac : float, optional
        CPU fraction to use
    mem_frac : float, optional
        Memory fraction to use
    remote_log: bool, optional
        Start remote logger

    Returns
    -------
    int
        Success message
    int
        Succeeded ms number
    int
        Failed ms number
    """
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)
    if flag_calibrators:
        flagdimension = "freqtime"
        flagfield_type = "cal"
        use_tfcrop = True
        flag_basename = f"flagging_{flagfield_type}_calibrator"
    else:
        flagdimension = "freq"
        flagfield_type = "target"
        use_tfcrop = False
        flag_basename = f"flagging_{flagfield_type}_target"
    logdir = f"{workdir}/logs"
    os.makedirs(logdir, exist_ok=True)
    logfile = f"{logdir}/{flag_basename}.log"
    if os.path.exists(logfile):
        os.remove(logfile)
    ctx = get_run_context()
    task_id = str(ctx.task_run.id)
    task_name = ctx.task_run.name
    stop_event = Event()
    log_thread_flag = start_log_task_saver(
        task_id, task_name, logfile, poll_interval=3, stop_event=stop_event
    )
    try:
        ##############
        print("###########################")
        print("Flagging ....")
        print("###########################")
        ########################
        # Calibrator ms flagging
        ########################
        with get_dask_client() as dask_client:
            msg, succeed, failed = flagging.main(
                mslist,
                metafits,
                workdir=workdir,
                outdir=outdir,
                datacolumn=datacolumn,
                flag_bad_ants=True,
                flag_bad_spw=False,
                use_tfcrop=use_tfcrop,
                flag_autocorr=True,
                flag_quack=flag_quack,
                flagdimension=flagdimension,
                restore_flag=restore_flag,
                flagbackup=False,
                cpu_frac=float(cpu_frac),
                mem_frac=float(mem_frac),
                logfile=logfile,
                jobid=jobid,
                start_remote_log=remote_log,
                dask_client=dask_client,
            )
    finally:
        stop_event.set()
        log_thread_flag.join(timeout=5)
    if msg != 0:
        raise RuntimeError("Flagging is failed.")
    else:
        return msg, succeed, failed


@task(
    name="importing_model_visibilities",
    retries=2,
    retry_delay_seconds=60,
    log_prints=True,
)
def run_import_model(
    mslist,
    metafits,
    workdir,
    jobid=0,
    cpu_frac=0.8,
    mem_frac=0.8,
    remote_log=False,
):
    """
    Importing calibrator models

    Parameters
    ----------
    mslist : str
        Name of the measurement sets (comma separated)
    metafits : str
        Metafits file
    workdir : str
        Working directory
    cpu_frac : float, optional
        CPU fraction to use
    mem_frac : float, optional
        Memory fraction to use
    remote_log: bool, optional
        Start remote logger

    Returns
    -------
    int
        Success message
    int
        Succeeded ms number
    int
        Failed ms number
    """
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)
    model_basename = "modeling"
    logdir = f"{workdir}/logs"
    os.makedirs(logdir, exist_ok=True)
    logfile = f"{logdir}/{model_basename}.log"
    if os.path.exists(logfile):
        os.remove(logfile)
    ctx = get_run_context()
    task_id = str(ctx.task_run.id)
    task_name = ctx.task_run.name
    stop_event = Event()
    log_thread_model = start_log_task_saver(
        task_id, task_name, logfile, poll_interval=3, stop_event=stop_event
    )
    try:
        ##############
        print("###########################")
        print("Importing model visibilities ....")
        print("###########################")
        ###################################
        # Calibrator ms visibility import
        ###################################
        with get_dask_client() as dask_client:
            msg, succeed, failed = import_model.main(
                mslist,
                metafits,
                workdir,
                cpu_frac=float(cpu_frac),
                mem_frac=float(mem_frac),
                logfile=logfile,
                jobid=jobid,
                start_remote_log=remote_log,
                dask_client=dask_client,
            )
    finally:
        stop_event.set()
        log_thread_model.join(timeout=5)
    if msg != 0:
        raise RuntimeError("Importing calibrator model is failed.")
    else:
        return msg, succeed, failed


@task(
    name="basic_calibration",
    retries=2,
    retry_delay_seconds=60,
    log_prints=True,
)
def run_basic_cal_jobs(
    mslist,
    metafits,
    workdir,
    outdir,
    perform_polcal=False,
    jobid=0,
    cpu_frac=0.8,
    mem_frac=0.8,
    keep_backup=False,
    remote_log=False,
):
    """
    Perform basic calibration

    Parameters
    ----------
    mslist: str
        Name of the measurement sets (comma seperated)
    metafits: str
        Metafits file
    workdir : str
        Working directory
    outdir : str
        Output directory
    perform_polcal : bool, optional
        Perform full polarization calibration
    cpu_frac : float, optional
        CPU fraction to use
    mem_frac : float, optional
        Memory fraction to use
    keep_backup : bool, optional
        Keep backups
    remote_log: bool, optional
        Start remote logger

    Returns
    -------
    int
        Success message for basic calibration
    int
        Succeeded ms number
    int
        Failed ms number
    """
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)
    cal_basename = "basic_cal"
    logdir = f"{workdir}/logs"
    os.makedirs(logdir, exist_ok=True)
    logfile = f"{logdir}/{cal_basename}.log"
    if os.path.exists(logfile):
        os.remove(logfile)
    ctx = get_run_context()
    task_id = str(ctx.task_run.id)
    task_name = ctx.task_run.name
    stop_event = Event()
    log_thread_cal = start_log_task_saver(
        task_id, task_name, logfile, poll_interval=3, stop_event=stop_event
    )
    try:
        ##############
        print("###########################")
        print("Performing basic calibration .....")
        print("###########################")
        ########################
        # Basic calibration
        ########################
        with get_dask_client() as dask_client:
            msg, succeed, failed = basic_cal.main(
                mslist,
                metafits,
                workdir,
                outdir,
                perform_polcal=perform_polcal,
                keep_backup=keep_backup,
                start_remote_log=remote_log,
                cpu_frac=float(cpu_frac),
                mem_frac=float(mem_frac),
                logfile=logfile,
                jobid=jobid,
                dask_client=dask_client,
            )
    finally:
        stop_event.set()
        log_thread_cal.join(timeout=5)
    if msg != 0:
        raise RuntimeError("Basic calibration is failed.")
    else:
        return msg, succeed, failed


@task(
    name="applying_basic_calibration",
    retries=2,
    retry_delay_seconds=60,
    log_prints=True,
)
def run_apply_basiccal_sol(
    mslist,
    calibrator_metafits,
    target_metafits,
    workdir,
    caldir,
    overwrite_datacolumn=True,
    only_amplitude=False,
    applymode="calflag",
    prefix="target",
    jobid=0,
    cpu_frac=0.8,
    mem_frac=0.8,
    remote_log=False,
):
    """
    Apply basic calibration solutions on splited target scans

    Parameters
    ----------
    mslist: str
        Target measurement set list (comma separated)
    calibrator_metafits : str
        Calibrator metafits
    target_metafits : str
        Target metafits
    workdir : str
        Working directory
    caldir : str
        Caltable directory
    overwrite_datacolumn : bool
        Overwrite data column or not
    only_amplitude : bool, optional
        Apply only amplitude part of the solution
    applymode : str, optional
        Applycal mode
    prefix : str, optional
        Applying on target of selfcal ms
    cpu_frac : float, optional
        CPU fraction to use
    mem_frac : float, optional
        Memory fraction to use
    remote_log: bool, optional
        Start remote logger

    Returns
    -------
    int
        Success message for applying calibration solutions and spliting target scans
    int
        Succeeded ms number
    int
        Failed ms number
    """
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)
    applycal_basename = f"apply_basiccal_{prefix}"
    logdir = f"{workdir}/logs"
    os.makedirs(logdir, exist_ok=True)
    logfile = f"{logdir}/{applycal_basename}.log"
    if os.path.exists(logfile):
        os.remove(logfile)
    ctx = get_run_context()
    task_id = str(ctx.task_run.id)
    task_name = ctx.task_run.name
    stop_event = Event()
    log_thread_apply = start_log_task_saver(
        task_id, task_name, logfile, poll_interval=3, stop_event=stop_event
    )
    try:
        ######################
        print("###########################")
        print("Applying basic calibration solutions on solar target .....")
        print("###########################")
        ######################
        # Applying basic calibration
        ######################
        with get_dask_client() as dask_client:
            msg, succeed, failed = do_apply_basiccal.main(
                mslist,
                calibrator_metafits,
                target_metafits,
                workdir,
                caldir,
                applymode=applymode,
                overwrite_datacolumn=overwrite_datacolumn,
                only_amplitude=only_amplitude,
                start_remote_log=remote_log,
                cpu_frac=float(cpu_frac),
                mem_frac=float(mem_frac),
                logfile=logfile,
                jobid=jobid,
                dask_client=dask_client,
            )
    finally:
        stop_event.set()
        log_thread_apply.join(timeout=5)
    if msg != 0:
        raise RuntimeError("Applying basic calibration solutions is failed.")
    else:
        return msg, succeed, failed


@task(
    name="solar_sidereal_correction",
    retries=2,
    retry_delay_seconds=60,
    log_prints=True,
)
def run_solar_siderealcor_jobs(
    mslist,
    workdir,
    prefix="target",
    jobid=0,
    cpu_frac=0.8,
    mem_frac=0.8,
    remote_log=False,
):
    """
    Apply sidereal motion correction of the Sun

    Parameters
    ----------
    mslist: str
        List of the measurement sets (comma separated)
    workdir : str
        Work directory
    prefix : str, optional
        Measurement set prefix
    cpu_frac : float, optional
        CPU fraction to use
    mem_frac : float, optional
        Memory fraction to use
    remote_log: bool, optional
        Start remote logger

    Returns
    -------
    int
        Success message
    int
        Succeeded ms number
    int
        Failed ms number
    """
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)
    sidereal_basename = f"cor_sidereal_{prefix}"
    logdir = f"{workdir}/logs"
    os.makedirs(logdir, exist_ok=True)
    logfile = f"{logdir}/{sidereal_basename}.log"
    if os.path.exists(logfile):
        os.remove(logfile)
    ctx = get_run_context()
    task_id = str(ctx.task_run.id)
    task_name = ctx.task_run.name
    stop_event = Event()
    log_thread_sidereal = start_log_task_saver(
        task_id, task_name, logfile, poll_interval=3, stop_event=stop_event
    )
    try:
        #######################
        print("###########################")
        print("Correcting sidereal motion .....")
        print("###########################")
        #######################
        # Sidereal motion correction
        #######################
        with get_dask_client() as dask_client:
            msg, succeed, failed = do_sidereal_cor.main(
                mslist,
                workdir=workdir,
                cpu_frac=float(cpu_frac),
                mem_frac=float(mem_frac),
                logfile=logfile,
                jobid=jobid,
                start_remote_log=remote_log,
                dask_client=dask_client,
            )
    finally:
        stop_event.set()
        log_thread_sidereal.join(timeout=5)
    if msg != 0:
        raise RuntimeError("Solar sidereal motion correction is failed.")
    else:
        return msg, succeed, failed


@task(
    name="selfcal",
    log_prints=True,
)
def run_selfcal_jobs(
    mslist,
    workdir,
    caldir,
    metafits,
    cal_applied,
    start_thresh=5.0,
    stop_thresh=3.0,
    max_iter=100,
    max_DR=100000,
    min_iter=5,
    conv_frac=0.3,
    solint="60s",
    do_apcal=True,
    do_polcal=True,
    solar_selfcal=True,
    keep_backup=False,
    uvrange="",
    minuv=0,
    weight="briggs",
    robust=0.0,
    applymode="calonly",
    min_tol_factor=1.0,
    jobid=0,
    cpu_frac=0.8,
    mem_frac=0.8,
    remote_log=False,
):
    """
    Self-calibration on target scans

    Parameters
    ----------
    mslist: str
        Target measurement set list (comma separated)
    workdir : str
        Working directory
    caldir : str
        Caltable directory
    metafits : str
        Metafits file
    cal_applied : bool
        Whether calibration solutions are applied or not
    cpu_frac : float, optional
        CPU fraction to use
    mem_frac : float, optional
        Memory fraction to use
    start_threshold : int, optional
        Start CLEAN threhold
    end_threshold : int, optional
        End CLEAN threshold
    max_iter : int, optional
        Maximum numbers of selfcal iterations
    max_DR : float, optional
        Maximum dynamic range
    min_iter : int, optional
        Minimum numbers of seflcal iterations at different stages
    conv_frac : float, optional
        Dynamic range fractional change to consider as converged
    uvrange : str, optional
        UV-range for calibration
    minuv : float, optionial
        Minimum UV-lambda to use in imaging
    weight : str, optional
        Image weighitng scheme
    robust : float, optional
        Robustness parameter for briggs weighting
    solint : str, optional
        Solutions interval
    do_apcal : bool, optional
        Perform ap-selfcal or not
    do_polcal : bool, optional
        Perform polarisation selfcal or not
    min_tol_factor : float, optional
        Minimum tolerance in temporal variation in imaging
    applymode : str, optional
        Solution apply mode
    solar_selfcal : bool, optional
        Whether is is solar selfcal or not
    remote_log: bool, optional
        Start remote logger

    Returns
    -------
    int
        Success message for self-calibration
    int
        Intensity self-calibration succeeded ms number
    int
        Intensity self-calibration failed ms number
    int
        Polarisation self-calibration succeed ms number
    int
        Polarisation self-calibration failed ms number
    """
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)
    selfcal_basename = "selfcal_target"
    logdir = f"{workdir}/logs"
    os.makedirs(logdir, exist_ok=True)
    logfile = f"{logdir}/{selfcal_basename}.log"
    if os.path.exists(logfile):
        os.remove(logfile)
    ctx = get_run_context()
    task_id = str(ctx.task_run.id)
    task_name = ctx.task_run.name
    stop_event = Event()
    log_thread_selfcal = start_log_task_saver(
        task_id, task_name, logfile, poll_interval=3, stop_event=stop_event
    )
    try:
        ########################
        print("###########################")
        print("Performing self-calibration of solar targets .....")
        print("###########################")
        ########################
        # Selfcal jobs
        ########################
        with get_dask_client() as dask_client:
            msg, int_succeed, int_failed, pol_succeed, pol_failed = do_selfcal.main(
                mslist,
                metafits,
                workdir,
                caldir,
                cal_applied=cal_applied,
                start_thresh=float(start_thresh),
                stop_thresh=float(stop_thresh),
                max_iter=float(max_iter),
                max_DR=float(max_DR),
                min_iter=float(min_iter),
                conv_frac=float(conv_frac),
                solint=solint,
                uvrange=uvrange,
                minuv=float(minuv),
                weight=weight,
                robust=float(robust),
                applymode=applymode,
                min_tol_factor=float(min_tol_factor),
                do_apcal=do_apcal,
                do_polcal=do_polcal,
                solar_selfcal=solar_selfcal,
                keep_backup=keep_backup,
                cpu_frac=float(cpu_frac),
                mem_frac=float(mem_frac),
                logfile=logfile,
                jobid=jobid,
                start_remote_log=remote_log,
                dask_client=dask_client,
            )
    finally:
        stop_event.set()
        log_thread_selfcal.join(timeout=5)
    if msg != 0:
        raise RuntimeError("Self-calibration is failed.")
    else:
        return msg, int_succeed, int_failed, pol_succeed, pol_failed


@task(
    name="applying_self-calibration",
    retries=2,
    retry_delay_seconds=60,
    log_prints=True,
)
def run_apply_selfcal_sol(
    mslist,
    metafits,
    workdir,
    caldir,
    overwrite_datacolumn=True,
    applymode="calflag",
    jobid=0,
    cpu_frac=0.8,
    mem_frac=0.8,
    remote_log=False,
):
    """
    Apply self-calibration solutions on splited target scans

    Parameters
    ----------
    mslist: str
        Target measurement set list (comma separated)
    metafits : str
        Metafits file
    workdir : str
        Working directory
    caldir : str
        Caltable directory
    applymode : str, optional
        Applycal mode
    cpu_frac : float, optional
        CPU fraction to use
    mem_frac : float, optional
        Memory fraction to use
    overwrite_datacolumn : bool
        Overwrite data column or not
    remote_log: bool, optional
        Start remote logger

    Returns
    -------
    int
        Success message for applying calibration solutions and spliting target scans
    int
        Succeeded gain solution ms number
    int
        Failed gain solution ms number
    int
        Succeeded polarisation solution ms number
    int
        Failed polarisation solution ms number
    """
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)
    applycal_basename = "apply_selfcal"
    logdir = f"{workdir}/logs"
    os.makedirs(logdir, exist_ok=True)
    logfile = f"{logdir}/{applycal_basename}.log"
    if os.path.exists(logfile):
        os.remove(logfile)
    ctx = get_run_context()
    task_id = str(ctx.task_run.id)
    task_name = ctx.task_run.name
    stop_event = Event()
    log_thread_applyselfcal = start_log_task_saver(
        task_id, task_name, logfile, poll_interval=3, stop_event=stop_event
    )
    try:
        ##################
        print("###########################")
        print("Applying self-calibration solutions on targets .....")
        print("###########################")
        ########################
        # Applying self-calibration
        ########################
        with get_dask_client() as dask_client:
            gain_succeed, gain_failed, pol_succeed, pol_failed = do_apply_selfcal.main(
                mslist,
                metafits,
                workdir,
                caldir,
                applymode=applymode,
                overwrite_datacolumn=overwrite_datacolumn,
                start_remote_log=remote_log,
                cpu_frac=float(cpu_frac),
                mem_frac=float(mem_frac),
                logfile=logfile,
                jobid=jobid,
                dask_client=dask_client,
            )
        if gain_failed == 0:
            msg = 0
        else:
            msg = 1
    finally:
        stop_event.set()
        log_thread_applyselfcal.join(timeout=5)
    if msg != 0:
        raise RuntimeError("Applying self-calibration solutions is failed.")
    else:
        return msg, gain_succeed, gain_failed, pol_succeed, pol_failed


@task(
    name="imaging",
    log_prints=True,
)
def run_imaging_jobs(
    mslist,
    workdir,
    outdir,
    freqrange="",
    timerange="",
    minuv=-1,
    weight="briggs",
    robust=0.0,
    pol="IQUV",
    freqres=1.28,
    timeres=10.0,
    threshold=1.0,
    use_multiscale=True,
    use_solar_mask=True,
    cutout_rsun=10.0,
    savemodel=False,
    saveres=False,
    jobid=0,
    cpu_frac=0.8,
    mem_frac=0.8,
    remote_log=False,
):
    """
    Imaging on target scans

    Parameters
    ----------
    mslist: str
        Target measurement set list (comma separated)
    workdir : str
        Working directory
    outdir : str
        Output image directory
    freqrange : str, optional
        Frequency range to image in MHz
    timerange : str, optional
        Time range to image (YYYY/MM/DD/hh:mm:ss~YYYY/MM/DD/hh:mm:ss)
    cpu_frac : float, optional
        CPU fraction to use
    mem_frac : float, optional
        Memory fraction to use
    minuv : float, optionial
        Minimum UV-lambda to use in imaging
    weight : str, optional
        Imaging weighting
    robust : float, optional
        Briggs weighting robust parameter (-1 to 1)
    pol : str, optional
        Stokes parameters to image
    freqres : float, optional
        Frequency resolution of spectral chunk in MHz (default : -1, no spectral chunking)
    timeres : float, optional
        Time resolution of temporal chunks in MHz (default : -1, no temporal chunking)
    threshold : float, optional
        CLEAN threshold in sigma
    use_multiscale : bool, optional
        Use multiscale or not
    use_solar_mask : bool, optional
        Use solar mask or not
    cutout_rsun : float, optional
        Cutout image size from center in solar radii (default : 10.0 solar radii)
    savemodel : bool, optional
        Save model images or not
    saveres : bool, optional
        Save residual images or not
    remote_log: bool, optional
        Start remote logger

    Returns
    -------
    int
        Success message for imaging
    int
        Succeeded ms number
    int
        Failed ms number
    int
        Total images
    """
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)
    imaging_basename = "imaging_target"
    logdir = f"{workdir}/logs"
    os.makedirs(logdir, exist_ok=True)
    logfile = f"{logdir}/{imaging_basename}.log"
    if os.path.exists(logfile):
        os.remove(logfile)
    ctx = get_run_context()
    task_id = str(ctx.task_run.id)
    task_name = ctx.task_run.name
    stop_event = Event()
    log_thread_imaging = start_log_task_saver(
        task_id, task_name, logfile, poll_interval=3, stop_event=stop_event
    )
    try:
        ######################
        print("###########################")
        print("Performing imaging of target scans .....")
        print("###########################")
        #######################
        # Performing imaging
        #######################
        with get_dask_client() as dask_client:
            msg, succeed, failed, total_images = do_imaging.main(
                mslist,
                workdir,
                outdir,
                freqrange=freqrange,
                timerange=timerange,
                pol=pol,
                freqres=float(freqres),
                timeres=float(timeres),
                weight=weight,
                robust=float(robust),
                minuv=float(minuv),
                threshold=float(threshold),
                cutout_rsun=float(cutout_rsun),
                use_multiscale=use_multiscale,
                use_solar_mask=use_solar_mask,
                savemodel=savemodel,
                saveres=saveres,
                start_remote_log=remote_log,
                cpu_frac=float(cpu_frac),
                mem_frac=float(mem_frac),
                jobid=jobid,
                logfile=logfile,
                dask_client=dask_client,
            )
    finally:
        stop_event.set()
        log_thread_imaging.join(timeout=5)
    if msg != 0:
        raise RuntimeError("Imaging is failed.")
    else:
        return msg, succeed, failed, total_images


@task(
    name="applying_primary_beam",
    log_prints=True,
)
def run_apply_pbcor(
    imagedir,
    metafits,
    workdir,
    leakage_dir="",
    jobid=0,
    cpu_frac=0.8,
    mem_frac=0.8,
    remote_log=False,
):
    """
    Apply primary beam corrections on all images

    Parameters
    ----------
    imagedir: str
        Image directory name
    metafits : str
        Metafits file
    workdir : str
        Work directory
    leakage_dir : str, optional
        Leakage dile directory
    cpu_frac : float, optional
        CPU fraction to use
    mem_frac : float, optional
        Memory fraction to use
    remote_log: bool, optional
        Start remote logger

    Returns
    -------
    int
        Success message for applying primary beam correction on all images
    int
        Succeeded image number
    int
        Failed image number
    """
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)
    applypbcor_basename = "apply_pbcor"
    logdir = f"{workdir}/logs"
    os.makedirs(logdir, exist_ok=True)
    logfile = f"{logdir}/{applypbcor_basename}.log"
    if os.path.exists(logfile):
        os.remove(logfile)
    ctx = get_run_context()
    task_id = str(ctx.task_run.id)
    task_name = ctx.task_run.name
    stop_event = Event()
    log_thread_pbcor = start_log_task_saver(
        task_id, task_name, logfile, poll_interval=3, stop_event=stop_event
    )
    try:
        ###################
        print("###########################")
        print("Applying primary beam corrections on all images .....")
        print("###########################")
        #####################
        # Applying primary beam correction
        #####################
        with get_dask_client() as dask_client:
            msg, succeed, failed = mwa_pbcor.main(
                imagedir,
                metafits,
                leakage_dir=leakage_dir,
                workdir=workdir,
                cpu_frac=float(cpu_frac),
                mem_frac=float(mem_frac),
                logfile=logfile,
                jobid=jobid,
                start_remote_log=remote_log,
                dask_client=dask_client,
            )
    finally:
        stop_event.set()
        log_thread_pbcor.join(timeout=5)
    if msg != 0:
        raise RuntimeError("Primary beam correction is failed.")
    else:
        return msg, succeed, failed


@task(
    name="making_overlay",
    log_prints=True,
)
def run_make_overlay(
    imagedir,
    outdir,
    workdir="",
    all_overlay=False,
    jobid=0,
    cpu_frac=0.8,
    remote_log=False,
):
    """
    Making overlays of all images on EUV images

    Parameters
    ----------
    imagedir : str
        Image directory name
    outdir : str
        Output directory
    workdir : str, optional
        Work directory
    all_overlay : bool, optional
        Whether to make overlays for all images or not
    cpu_frac : float, optional
        CPU fraction to use
    remote_log: bool, optional
        Start remote logger

    Returns
    -------
    int
        Success message for EUV overlays
    int
        Succeeded image number
    int
        Failed image number
    """
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)
    overlay_basename = "do_overlay"
    logdir = f"{workdir}/logs"
    os.makedirs(logdir, exist_ok=True)
    logfile = f"{logdir}/{overlay_basename}.log"
    if os.path.exists(logfile):
        os.remove(logfile)
    ctx = get_run_context()
    task_id = str(ctx.task_run.id)
    task_name = ctx.task_run.name
    stop_event = Event()
    log_thread_overlay = start_log_task_saver(
        task_id, task_name, logfile, poll_interval=3, stop_event=stop_event
    )
    os.makedirs(outdir, exist_ok=True)
    try:
        ###################
        print("###########################")
        print("Making overlays of images .....")
        print("###########################")
        #####################
        # Making overlays
        #####################
        with get_dask_client() as dask_client:
            msg, succeed, failed = make_mwa_overlay.main(
                imagedir,
                outdir,
                workdir=workdir,
                all_overlay=all_overlay,
                cpu_frac=float(cpu_frac),
                logfile=logfile,
                jobid=jobid,
                start_remote_log=remote_log,
                dask_client=dask_client,
            )
    finally:
        stop_event.set()
        log_thread_overlay.join(timeout=5)
    if msg != 0:
        raise RuntimeError("EUV overlay is failed.")
    else:
        return msg, succeed, failed


@task(
    name="making_msplot",
    retries=2,
    retry_delay_seconds=60,
    log_prints=True,
)
def run_make_msplot(
    mslist,
    workdir,
    outdir,
    jobid=0,
    cpu_frac=0.8,
    mem_frac=0.8,
    remote_log=False,
):
    """
    Making diagnostic plots of measurement sets

    Parameters
    ----------
    mslist : str
        Measurement set list (comma separated)
    workdir : str, optional
        Work directory
    outdir : str
        Output directory
    cpu_frac : float, optional
        CPU fraction to use
    mem_frac : float, optional
        Memory fraction to use
    remote_log: bool, optional
        Start remote logger

    Returns
    -------
    int
        Success message for measurement set ploting
    """
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)
    msplot_basename = "do_msplot"
    logdir = f"{workdir}/logs"
    os.makedirs(logdir, exist_ok=True)
    logfile = f"{logdir}/{msplot_basename}.log"
    if os.path.exists(logfile):
        os.remove(logfile)
    ctx = get_run_context()
    task_id = str(ctx.task_run.id)
    task_name = ctx.task_run.name
    stop_event = Event()
    log_thread_overlay = start_log_task_saver(
        task_id, task_name, logfile, poll_interval=3, stop_event=stop_event
    )
    os.makedirs(outdir, exist_ok=True)
    try:
        ###################
        print("###########################")
        print("Making diagnostic plots of all measurement sets .....")
        print("###########################")
        #####################
        # Making plots
        #####################
        with get_dask_client() as dask_client:
            msg = make_ms_plot.main(
                mslist,
                workdir,
                outdir,
                cpu_frac=float(cpu_frac),
                mem_frac=float(mem_frac),
                logfile=logfile,
                jobid=jobid,
                start_remote_log=remote_log,
                dask_client=dask_client,
            )
    finally:
        stop_event.set()
        log_thread_overlay.join(timeout=5)
    if msg != 0:
        raise RuntimeError("Measurement set diagnostic ploting is failed.")
    else:
        return msg


def send_task_notification(emails, msg, jobid, obsid, logger_timestamp):
    """
    Send notification after each task is finished

    Parameters
    ----------
    emails : str
        E-mail ids
    msg : str
        Notification message
    jobid : int
        JobID
    obsid : int
        Observation ID
    logger_timestamp : str
        Logger timestamp
    """
    email_subject = f"P-AIRCARS Logger Details: {logger_timestamp}, OBSID: {obsid}"
    email_msg = f"{msg}"
    success_msg, error_msg = send_notification(emails, email_subject, email_msg)


@flow(
    name="P-AIRCARS Master control",
    version="3.0",
    description="Calibration and Imaging Pipeline for MWA Solar Observation",
    log_prints=True,
)
def master_control(
    target_datadir,
    target_metafits,
    workdir,
    outdir,
    calibrator_datadir="",
    calibrator_metafits="",
    solar_data=True,
    # Pre-calibration
    do_forcereset_weightflag=False,
    do_cal_flag=True,
    do_import_model=True,
    # Basic calibration
    do_basic_cal=True,
    do_applycal=True,
    only_amplitude=False,
    # Target data preparation
    freqrange="",
    timerange="",
    uvrange="",
    # Polarization self-calibration
    do_polcal=False,
    # Self-calibration
    do_selfcal=True,
    do_apply_selfcal=True,
    do_ap_selfcal=True,
    solar_selfcal=True,
    use_solar_mask=True,
    solint="60s",
    # Sidereal correction
    do_sidereal_cor=False,
    do_move_solarcenter=True,
    # Dynamic spectra
    make_ds=True,
    # Imaging
    do_imaging=True,
    do_pbcor=True,
    weight="briggs",
    robust=0.0,
    minuv=0,
    image_freqres=1.28,
    image_timeres=10.0,
    pol="IQUV",
    clean_threshold=1.0,
    use_multiscale=True,
    cutout_rsun=10.0,
    make_overlay=False,
    make_msplot=False,
    # Resource settings
    cpu_frac=0.8,
    mem_frac=0.8,
    max_worker=2,
    keep_backup=False,
    keep_calibrated_ms=True,
    # Remote logging
    masterlog=None,
    remote_logger=False,
    jobid=None,
    job_password=None,
    adaptive=False,
):
    """
    Master controller of the entire pipeline

    Parameters
    ----------
    target_datadir : str
        Target measurement set directory
    target_metafits : str
        Target metafits file
    workdir : str
        Work directory path
    outdir : str
        Output directory
    calibrator_datadir : str, optional
        Calibrator data directory
    calibrator_metafits : str, optional
        Calibrator metafits file
    solar_data : bool, optional
        Whether it is solar data or not

    do_forcereset_weightflag : bool, optional
        Reset weights and flags of the input ms
    do_cal_flag : bool, optional
        Perform flagging on calibrator
    do_import_model : bool, optional
        Import model visibilities of flux and polarization calibrators

    do_basic_cal : bool, optional
        Perform basic calibration
    do_applycal : bool, optional
        Apply basic calibration on target scans
    only_amplitude : bool, optional
        Apply only amplitude part of gain solution from calibrator

    freqrange : str, optional
        Frequency range to image in MHz (xx1~xx2,xx3~xx4,)
    timerange : str, optional
        Time range to image in YYYY/MM/DD/hh:mm:ss format (tt1~tt2,tt3~tt4,...)
    uvrange : str, optional
        UV-range for calibration

    do_polcal : bool, optional
        Perform full-polarization calibration and imaging

    do_selfcal : bool, optional
        Perform self-calibration
    do_apply_selfcal : bool, optional
        Apply self-calibration solutions
    do_ap_selfcal : bool, optional
        Perform amplitude-phase self-cal or not
    solint : str, optional
        Solution intervals in self-cal
    solar_selfcal : bool, optional
        Solar selfcal
    use_solar_mask : bool, optional
        Use solar mask or not

    do_sidereal_cor : bool, optional
        Perform solar sidereal motion correction or not
    do_move_solarcenter: boo, optional
        Move phasecenter to solar center
    make_ds : bool, optional
        Make dynamic spectra

    do_imaging : bool, optional
        Perform final imaging
    do_pbcor : bool, optional
        Perform primary beam correction
    weight : str, optional
        Image weighting
    robust : float, optional
        Robust parameter for briggs weighting (-1 to 1)
    minuv : float, optional
        Minimum UV-lambda for final imaging
    image_freqres : float, optional
        Image frequency resolution in MHz (-1 means full bandwidth)
    image_timeres : float, optional
        Image temporal resolution in seconds (-1 means full scan duration)
    pol : str, optional
        Stokes parameters of final imaging
    clean_threshold : float, optional
        CLEAN threshold of final imaging in sigma
    use_multiscale : bool, optional
        Use multiscale scales or not
    cutout_rsun : float, optional
        Cutout image size from center in solar radii (default : 10.0 solar radii)
    make_overlay : bool, optional
        Make EUV MWA overlay for all images or not (default : per coarse channel images will be overlaid at 60s intervals)
    make_msplot : bool, optional
        Make diagnostic plots of measurement sets

    cpu_frac : float, optional
        CPU fraction to use
    mem_frac : float, optional
        Memory fraction to use
    max_worker: int, optional
        Maximum workers
    keep_backup : bool, optional
        Keep backup of of all intermediate data poducts, calibrated ms, self-cal rounds and including images, models and residual images
    keep_calibrated_ms : bool, optional
        Keep calibrated measurement sets or not

    masterlog : str, optional
        Master logfile
    remote_logger : bool, optional
        Enable remote logging of the pipeline status
    jobid : str, optional
        Job ID
    job_password : str, optional
        User specified job password for remote logger
    adaptive : bool, optional
        Whether do adaptive scaling or not

    Returns
    -------
    int
        Success message
    """
    print("P-AIRCARS workfkow started...")
    if target_datadir.startswith("~"):
        print("Please provide full path of target directory.")
        return 1
    else:
        target_datadir = os.path.abspath(target_datadir)
        if os.path.exists(target_datadir) is False:
            print(
                f"Target data directory: {target_datadir} does not exist. Provide correct full path."
            )
            return 1
    if target_metafits.startswith("~"):
        print("Please provide full path of target metafits.")
        return 1
    else:
        target_metafits = os.path.abspath(target_metafits)
        if os.path.exists(target_metafits) is False:
            print(f"Target metafits: {target_metafits} does not exist.")
            return 1
    if calibrator_datadir.startswith("~"):
        print("Please provide full path of calibrator data directory.")
        return 1
    else:
        if calibrator_datadir != "":
            calibrator_datadir = os.path.abspath(calibrator_datadir)
    if calibrator_metafits.startswith("~"):
        print("Please provide full path of calibrator metafits.")
        return 1
    else:
        if calibrator_metafits != "":
            calibrator_metafits = os.path.abspath(calibrator_metafits)
    if workdir.startswith("~"):
        print("Please provide full path of work directory.")
        return 1
    else:
        workdir = os.path.abspath(workdir)
    if outdir.startswith("~"):
        print("Please provide full path of output directory.")
        return 1
    else:
        outdir = os.path.abspath(outdir)

    if jobid is None:
        jobid = get_jobid()
    #############################################
    # Listing target ms
    #############################################
    target_mslist = sorted(glob.glob(f"{target_datadir}/*.ms"))
    if len(target_mslist) == 0:
        print(
            f"No measurement set is present in target data directory: {target_datadir}"
        )
        if emails != "":
            email_msg = "No measurement set is present in the target data directory."
            send_task_notification(emails, email_msg, jobid, target_obsid, timestamp)
    test_msname = target_mslist[0]
    if os.path.exists(target_metafits) is False:
        target_obsid = get_MWA_OBSID(test_msname)
        try:
            target_metafits = download_MWA_metafits(
                target_obsid, outdir=os.path.dirname(test_msname)
            )
        except Exception:
            tracebcak.print_exc()
            target_metafits = None
        if target_metafits is None or os.path.exists(target_metafits) is False:
            print(
                f"Target metafits {target_metafits} does not exist. P-AIRCARS has stopped."
            )
            if emails != "":
                email_msg = "Target metafits file does not exist."
                send_task_notification(emails, email_msg, jobid, "N/A", timestamp)
            return 1
    target_header = fits.getheader(target_metafits)
    target_obsid = target_header["GPSTIME"]
    target_freq_config = target_header["CHANNELS"]
    target_coarse_chans = [get_MWA_coarse_chan(ms) for ms in target_mslist]

    #########################################
    # Verifying target obsid
    #########################################
    target_ms_obsids = [get_MWA_OBSID(ms) for ms in target_mslist]
    all_same_obsids = all(x == target_ms_obsids[0] for x in target_ms_obsids)
    if all_same_obsids is False:
        print(
            "All target measurement sets are not belong to same OBSID. Keep only measurement sets with same OBSID inside the target directory. P-AIRCARS has stopped."
        )
        return 1
    else:
        target_ms_obsid = target_ms_obsids[0]
        if target_ms_obsid != target_obsid:
            print(
                f"Target measurement set OBSID: {target_ms_obsid} is different from metafits provided OBSID: {target_obsid}. Provide correct metafits file. P-AIRCARS has stopped."
            )
            return 1

    ###################################
    # Preparing working directories
    ###################################
    print("Preparing working directories....")
    if workdir == "":
        workdir = os.path.dirname(os.path.abspath(target_mslist[0])) + "/workdir"
    workdir = workdir.rstrip("/")
    if outdir == "":
        outdir = workdir
    workdir = f"{workdir}/{target_obsid}"
    try:
        os.makedirs(workdir, exist_ok=True)
    except Exception:
        print(
            f"Work directory: {workdir} can not be created. Please check the path carefully."
        )
        traceback.print_exc()
        return 1

    os.chdir(workdir)
    scheduler_name = get_scheduler_name()
    #################################
    # Setup logger
    #################################
    logdir = f"{workdir}/logs"
    os.makedirs(logdir, exist_ok=True)
    if (
        scheduler_name != "local"
        or masterlog is None
        or os.path.exists(masterlog) is False
    ):
        master_logfile = f"{logdir}/main.log"
        ctx = get_run_context()
        flow_id = str(ctx.flow_run.id)
        flow_name = ctx.flow_run.name
        stop_event = Event()
        log_thread_flow = start_flow_log_saver(
            flow_id, flow_name, master_logfile, poll_interval=3, stop_event=stop_event
        )
        master_log_created = True
    else:
        master_log_created = False
        master_logfile = f"{logdir}/main.log"
        if os.path.exists(master_logfile):
            os.system(f"rm -rf {master_logfile}")
        os.symlink(masterlog, master_logfile)

    dask_dir = None
    try:
        dask_client = get_client()
        dask_cluster = dask_client.cluster
    except:
        if mem_frac <= 0:
            mem_frac = 0.8
        result = get_local_dask_cluster(
            workdir,
            mem_frac=mem_frac,
        )
        if result is None:
            print("Error occured in creating local cluster.")
            return 1
        else:
            dask_client, dask_cluster, dask_dir = result

    #####################################
    # Initiating paircars data
    #####################################
    init_paircars_data()
    observer = None

    ###################################################
    # Measurement set check and other working directory
    ###################################################
    outdir = outdir.rstrip("/")
    outdir = f"{outdir}/{target_obsid}"
    caldir = f"{outdir}/caltables"
    caldir = caldir.rstrip("/")
    try:
        os.makedirs(outdir, exist_ok=True)
    except Exception:
        print(
            f"Output directory: {outdir} can not created. Please check the path carefully."
        )
        traceback.print_exc()
        return 1
    os.makedirs(caldir, exist_ok=True)

    max_worker = max(2, max_worker)  # Minimum 2 workers are needed
    cpu_frac = min(0.8, abs(cpu_frac))
    mem_frac = min(0.8, abs(mem_frac))

    n_threads = os.environ.get("OMP_NUM_THREADS")
    if n_threads is None:
        n_threads = 1
    else:
        n_threads = max(1, int(n_threads))

    observer = None
    try:
        #####################################
        # Reading remotelink and emails
        #####################################
        remote_link = ""
        if remote_logger:
            trial = 0
            while trial <= 5:
                try:
                    remote_link = get_remote_logger_link()
                except Exception:
                    traceback.print_exc()
                    pass
                if remote_link != "":
                    break
                else:
                    time.sleep(5)
                    trial += 1
            if remote_link == "":
                print("Please provide a valid remote link.")
                remote_logger = False
    except Exception:
        traceback.print_exc()
        remote_logger = False

    ###############################
    # Setting up email notification
    ###############################
    try:
        emails = get_emails()
        if not remote_logger:
            timestamp = dt.utcnow().strftime("%Y-%m-%dT%H:%M:%S")
            if emails != "":
                email_subject = (
                    f"P-AIRCARS Logger Details: {timestamp}, OBSID: {target_obsid}"
                )

                email_msg = f"P-AIRCARS Job ID: {jobid}"
                success_msg, error_msg = send_notification(
                    emails, email_subject, email_msg
                )
        else:
            ####################################
            # Job name and logging password
            ####################################
            hostname = socket.gethostname()
            timestamp = dt.utcnow().strftime("%Y-%m-%dT%H:%M:%S")
            username = getpass.getuser()
            jobname = f"{username}-{hostname}:{timestamp}:{target_obsid}"
            timestamp1 = dt.utcnow().strftime("%Y%m%dT%H%M%S")
            remote_job_id = f"{hostname}_{timestamp1}_{target_obsid}"
            if job_password is None:
                password = generate_password()
            else:
                password = job_password
            np.save(
                f"{workdir}/.jobname_password.npy",
                np.array([jobname, password], dtype="object"),
            )
            ############
            # Logger
            ############
            if os.path.exists(f"{workdir}/.jobname_password.npy"):
                time.sleep(5)
                jobname, password = np.load(
                    f"{workdir}/.jobname_password.npy", allow_pickle=True
                )
                if master_logfile is not None and os.path.exists(master_logfile):
                    observer = init_logger(
                        "master_log", master_logfile, jobname=jobname, password=password
                    )
            if observer == None:
                print(
                    "Remote link or jobname is blank. Not transmiting to remote logger."
                )
    except Exception:
        traceback.print_exc()
        emails = ""

    try:
        #####################
        # Notify over email
        #####################
        if emails != "":
            email_subject = (
                f"P-AIRCARS Logger Details: {timestamp}, OBSID: {target_obsid}"
            )

            email_msg = (
                f"P-AIRCARS Job ID: {jobid}\n"
                f"Remote logger Job ID: {jobname}\n"
                f"Remote access password: {password}"
            )
            success_msg, error_msg = send_notification(emails, email_subject, email_msg)
        #####################################
        # Printing basic info of the pipeline
        #####################################
        print("###########################")
        print(f"P-AIRCARS Job ID: {jobid}")
        print(f"Work directory: {workdir}")
        print(f"Final product directory: {outdir}")
        print("###########################")
        if remote_logger:
            print(
                "############################################################################"
            )
            print(remote_link)
            print(f"Job ID: {jobname}")
            print(f"Remote access password: {password}")
            print(
                "#############################################################################"
            )

        ############################################
        # Determining where to use calibrator or not
        #############################################
        calibrator_mslist = sorted(glob.glob(f"{calibrator_datadir}/*.ms"))
        calibrator_obsid = None
        if len(calibrator_mslist) == 0:
            print(
                f"No calibrator observation is provided. Continuing based on self-calibration."
            )
            has_cal = False
        ######################################################
        # Downloading calibrator metafits if it does not exist
        ######################################################
        if os.path.exists(calibrator_metafits) is False:
            test_cal_ms = calibrator_mslist[0]
            cal_obsid = get_MWA_OBSID(test_cal_ms)
            try:
                calibrator_metafits = download_MWA_metafits(
                    cal_obsid, outdir=os.path.dirname(test_cal_ms)
                )
            except Exception:
                tracebcak.print_exc()
                calibrator_metafits = None
        if calibrator_metafits is not None and os.path.exists(calibrator_metafits):
            calibrator_header = fits.getheader(calibrator_metafits)
            calibrator_obsid = calibrator_header["GPSTIME"]
            calibrator_freq_config = calibrator_header["CHANNELS"]
            if np.abs(calibrator_obsid - target_obsid) > 12 * 3600:
                print("Calibrator observations were taken 12 hours apart.")
                has_cal = True
            elif target_freq_config != calibrator_freq_config:
                print(f"Target coarse channels: {target_freq_config}.")
                print(f"Calibrator coarse channels: {calibrator_freq_config}.")
                print("Calibrator and target frequency configuration is different.")
                has_cal = False
            else:
                has_cal = True
        else:
            print(
                f"Calibrator ms is available, however, calibrator metafits is not available."
            )
            has_cal = False

        #########################################
        # Verifying calibrator obsid
        #########################################
        if has_cal:
            cal_ms_obsids = [get_MWA_OBSID(ms) for ms in calibrator_mslist]
            all_same_obsids = all(x == cal_ms_obsids[0] for x in cal_ms_obsids)
            if all_same_obsids is False:
                print(
                    "All calibrator measurement sets are not belong to same OBSID. Keep only measurement sets with same OBSID inside the calibrator directory."
                )
                print("P-AIRCARS will not use calibrators.")
                has_cal = False
            else:
                cal_ms_obsid = cal_ms_obsids[0]
                if cal_ms_obsid != calibrator_obsid:
                    print(
                        f"Calibrator measurement set OBSID: {cal_ms_obsid} is different from metafits provided OBSID: {calibrator_obsid}. Provide correct metafits file."
                    )
                    print("P-AIRCARS will not use calibrators.")
                    has_cal = False

        ######################################################
        # Filtering only matching coarse channel calibrator ms
        ######################################################
        if has_cal:
            print("Filtering calibrator measurement sets...")
            filtered_calms = []
            for ms in calibrator_mslist:
                coarse_chan = get_MWA_coarse_chan(ms)
                if coarse_chan in target_coarse_chans:
                    filtered_calms.append(ms)
                    print(
                        f"Coarse channel: {coarse_chan} of calibrator measurement set: {ms} is used."
                    )
            calibrator_mslist = filtered_calms

        #####################################
        # Settings for solar data
        #####################################
        if solar_data:
            if not use_solar_mask:
                print("Use solar mask during CLEANing.")
                use_solar_mask = True
            if not solar_selfcal:
                solar_selfcal = True
            full_FoV = False
        else:
            if use_solar_mask:
                print("Stop using solar mask during CLEANing.")
                use_solar_mask = False
            if solar_selfcal:
                solar_selfcal = False
            full_FoV = True

        #####################################################################
        # Checking if ms is full pol for polarization calibration and imaging
        #####################################################################
        if do_polcal:
            print(
                "Checking measurement set suitability for polarization calibration...."
            )
            for msname in target_mslist:
                msmd = msmetadata()
                msmd.open(msname)
                npol = msmd.ncorrforpol()[0]
                msmd.close()
                if npol < 4:
                    print(
                        f"Measurement set: {ms} is not full-polar. Do not performing polarization analysis."
                    )
                    do_polcal = False
                    break

        #################################################
        # Determining maximum allowed frequency averaging
        #################################################
        print("Estimating optimal frequency averaging....")
        max_freqres_list = []
        freqres_list = []
        msmd = msmetadata()
        for msname in target_mslist:
            max_freqres = calc_bw_smearing_freqwidth(msname, full_FoV=full_FoV)
            max_freqres_list.append(max_freqres)
            msmd.open(msname)
            freqres = msmd.chanres(0, unit="MHz")[0]
            msmd.close()
            freqres_list.append(freqres)
        freqres = min(freqres_list)
        if freqres > 0.16:
            print(
                f"Frequency resolution: {round(freqres*1000,1)}kHz is more than 160kHz. Assuming channel flagging is already done before averaing."
            )
        max_freqres = min(max_freqres_list)
        if image_freqres > 0:
            image_freqres = max(image_freqres, freqres)
            freqavg = round(min(image_freqres, max_freqres), 2)
        else:
            freqavg = freqres
        freqavg = min(0.16, freqavg)
        image_freqres = round(image_freqres, 2)

        total_ncoarse = 0
        for msname in target_mslist:
            ncoarse = get_ncoarse(msname)
            total_ncoarse += ncoarse
        total_ncoarse = max(1, total_ncoarse)

        ################################################
        # Determining maximum allowed temporal averaging
        ################################################
        print("Estimating optimal temporal averaging....")
        max_timeres_list = []
        timeres_list = []
        for msname in target_mslist:
            if solar_data:  # For solar data, it is assumed Sun is tracked.
                max_timeres = calc_time_smearing_timewidth(msname)
            else:
                max_timeres = min(
                    calc_time_smearing_timewidth(msname),
                    max_time_solar_smearing(msname),
                )
            max_timeres_list.append(max_timeres)
            msmd.open(msname)
            times = msmd.timesforspws(0)
            timeres = np.nanmean(np.diff(times))
            msmd.close()
            timeres_list.append(timeres)
        timeres = min(timeres_list)
        quack_timestamps = int(4.0 / timeres)
        max_timeres = min(max_timeres_list)
        if image_timeres > (2 * 3660):  # If more than 2 hours
            print(
                "Image time integration is more than 2 hours, which may cause smearing due to solar differential rotation."
            )
        if image_timeres > 0:
            image_timeres = max(image_timeres, timeres)
            timeavg = round(min(image_timeres, max_timeres), 2)
        else:
            timeavg = timeres
        timeavg = min(2.0, timeavg)
        image_timeres = round(image_timeres, 2)
        print(f"Frequency resolution: {freqres}MHz, time resolution: {timeres}s.")
        print(f"Frequency averaging: {freqavg}MHz, time averaging: {timeavg}s.")
        print(
            f"Imaging frequency resolution: {image_freqres}MHz, time resolution: {image_timeres}s."
        )

        #############################
        # Reset any previous weights
        ############################
        print("Resetting previous flags and weights....")
        for msname in target_mslist:
            reset_weights_and_flags(
                msname, n_threads=n_threads, force_reset=do_forcereset_weightflag
            )
        for msname in calibrator_mslist:
            reset_weights_and_flags(
                msname, n_threads=n_threads, force_reset=do_forcereset_weightflag
            )
        print("Reset is done.")

        if (move_solarcenter or make_ds) and adaptive:
            scale_worker_and_wait(
                dask_cluster, dask_client, min(len(target_mslist) + 1, max_worker)
            )

        ########################################
        # Moving phasecenter to the solar center
        ########################################
        if solar_data and do_move_solarcenter:
            if emails != "":
                email_msg = "Started moving phasecenter to solar center."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print("Starting task: Moving phasecenter to the Sun .....")
            print("###########################")
            future_movecenter = run_solar_phasecenter_jobs.with_options(
                task_run_name=f"moving_to_solar_center_{jobid}",
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
                    email_msg = f"Moving phasecenter to solar center is done.\nSucceeded: {succeed}, failed: {failed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(f"Finished task: Moving phasecenter to solar center is done.")
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
                if adaptive and len(filtered_ms) != len(target_mslist):
                    scale_worker_and_wait(
                        dask_cluster,
                        dask_client,
                        min(len(target_mslist) + 1, max_worker),
                    )
                target_mslist = filtered_ms  # Filtered target mslist
            except Exception:
                print(
                    "Error in moving phasecenter to solar center. P-AIRCARS has stopped."
                )
                traceback.print_exc()
                if emails != "":
                    email_msg = "Error occured in moving phasecenter to solar center."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                return 1

        #######################################
        # Run dynamic spectra making
        #######################################
        if solar_data and make_ds:
            if emails != "":
                email_msg = "Started making solar dynamic spectra."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print("Starting task: Making dynamic spectra of solar target .....")
            print("###########################")
            future_maskms = run_ds_jobs.with_options(
                task_run_name=f"making_dynamic_spectra_{jobid}",
            ).submit(
                ",".join(target_mslist),
                target_metafits,
                workdir,
                outdir,
                jobid=jobid,
                cpu_frac=round(cpu_frac, 2),
                mem_frac=round(mem_frac, 2),
                remote_log=remote_logger,
            )
            try:
                msg, succeed, failed = future_maskms.result()
                if emails != "":
                    email_msg = f"Making solar dynamic spectra are done.\nSucceeded: {succeed}, failed: {failed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(f"Finished task: Making solar dynamic spectra are done.")
                print("###########################")
            except Exception:
                print("!!! WARNING : Error in making dynamic spectra. !!!")
                traceback.print_exc()
                if emails != "":
                    email_msg = "Error occured in making dynamic spectra."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )

        ##############################
        # Run spliting jobs
        ##############################
        # If basic calibration is requested and calibrator ms and metafits are present
        if (do_basic_cal or do_cal_flag or do_import_model) and has_cal:
            if adaptive:
                scale_worker_and_wait(
                    dask_cluster, dask_client, min(total_ncoarse + 1, max_worker)
                )
            prefix = "calibrator"
            if emails != "":
                email_msg = "Started spliting of calibrator measurement sets."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print(f"Starting task: Spliting {prefix} .....")
            print("###########################")
            future_cal_split = run_target_split_jobs.with_options(
                task_run_name=f"spliting_{prefix}_{jobid}"
            ).submit(
                ",".join(calibrator_mslist),
                calibrator_metafits,
                workdir,
                datacolumn="data",
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
                    email_msg = f"Spliting of calibrator measurement sets are done.\nExpected: {expected}, succeeded: {succeed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(
                    f"Finished task: Spliting of calibrator measurement sets are done."
                )
                print("###########################")
            except Exception:
                print(
                    "!!!! WARNING: Error in spliting calibrator measurement sets. !!!!"
                )
                traceback.print_exc()
                if emails != "":
                    email_msg = "Spliting calibrator measurement set is failed."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                has_cal = False

        if (do_cal_flag or do_import_model or do_basic_cal) and has_cal:
            split_cal_mslist = sorted(glob.glob(f"{workdir}/calibrator*_spw_*.ms"))
            if len(split_cal_mslist) == 0:
                print("No splited measurement set is present for basic calibration.")
                has_cal = False
                if adaptive:
                    scale_worker_and_wait(
                        dask_cluster,
                        dask_client,
                        1,
                    )
            else:
                if adaptive:
                    scale_worker_and_wait(
                        dask_cluster,
                        dask_client,
                        min(len(split_cal_mslist) + 1, max_worker),
                    )

        ##################################
        # Run flagging jobs on calibrators
        ##################################
        # Only if basic calibration is requested
        if do_cal_flag and has_cal:
            if emails != "":
                email_msg = "Started flagging of calibrators."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print("Starting task: Flagging calibrators ....")
            print("###########################")
            future_flag = run_flag.with_options(
                task_run_name=f"flagging_cal_{jobid}"
            ).submit(
                ",".join(split_cal_mslist),
                calibrator_metafits,
                workdir,
                outdir,
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
                    email_msg = f"Flagging of calibrator is done.\nSucceeded: {succeed}, failed: {failed}."
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
                if adaptive and len(filtered_ms) != len(split_cal_mslist):
                    scale_worker_and_wait(
                        dask_cluster,
                        dask_client,
                        min(len(split_cal_mslist) + 1, max_worker),
                    )
                split_cal_mslist = filtered_ms  # Filtered target mslist
                print("###########################")
                print(f"Finished task: Flagging of calibrator is done.")
                print("###########################")
            except Exception:
                print("!!!! WARNING: Flagging error. P-AIRCARS has stopped. !!!!")
                traceback.print_exc()
                if emails != "":
                    email_msg = "Error in flagging calibrators."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                return 1

        #################################
        # Import model
        #################################
        if do_import_model and has_cal:
            if emails != "":
                email_msg = "Started importing sky model for calibrator."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print("Starting task: Importing model visibilities ....")
            print("###########################")
            future_import_model = run_import_model.with_options(
                task_run_name=f"importing_model_visibilities_{jobid}"
            ).submit(
                ",".join(split_cal_mslist),
                calibrator_metafits,
                workdir,
                jobid=jobid,
                cpu_frac=round(cpu_frac, 2),
                mem_frac=round(mem_frac, 2),
                remote_log=remote_logger,
            )
            try:
                msg, succeed, failed = future_import_model.result()
                if emails != "":
                    email_msg = f"Model import for calibrator is done.\nSucceeded: {succeed}, failed: {failed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(f"Finished task: Model import for calibrator is done.")
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
                if adaptive and len(filtered_ms) != len(split_cal_mslist):
                    scale_worker_and_wait(
                        dask_cluster,
                        dask_client,
                        min(len(split_cal_mslist) + 1, max_worker),
                    )
                split_cal_mslist = filtered_ms  # Filtered target mslist
            except Exception:
                print(
                    "!!!! WARNING: Error in importing calibrator models. Not continuing calibration. !!!!"
                )
                traceback.print_exc()
                if emails != "":
                    email_msg = "Error occured in importing model for calibrators. Not using calibrator solutions."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                has_cal = False
                if do_selfcal is False:
                    print(
                        "Self-calibration is also switched off. P-AIRCARS has stopped."
                    )
                    return 1

        ###############################
        # Run basic calibration
        ###############################
        if do_basic_cal and has_cal:
            if emails != "":
                email_msg = "Started basic calibration."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print("Starting task: Performing basic calibration .....")
            print("###########################")
            future_basical = run_basic_cal_jobs.with_options(
                task_run_name=f"basic_calibration_{jobid}"
            ).submit(
                ",".join(split_cal_mslist),
                calibrator_metafits,
                workdir,
                outdir,
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
                    email_msg = f"Basic calibration is done.\nSucceeded: {succeed}, failed: {failed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(f"Finished task: Basic calibration is done.")
                print("###########################")
            except Exception:
                print(
                    "!!!! WARNING: Error in basic calibration. Starting without basic calibration. !!!!"
                )
                traceback.print_exc()
                if emails != "":
                    email_msg = "Error occured in basic calibration."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                has_cal = False

        if (do_cal_flag or do_import_model or do_basic_cal) and adaptive:
            scale_worker_and_wait(dask_cluster, dask_client, 1)

        ##########################################
        # Checking presence of necessary caltables
        ##########################################
        if calibrator_obsid is not None:
            print(
                f"Searching for bandpass tables: {caldir}/calibrator_{calibrator_obsid}*.bcal"
            )
            bandpass_tables = sorted(
                glob.glob(f"{caldir}/calibrator_{calibrator_obsid}*.bcal")
            )
            bandpass_tables = interpolate_bpass(bandpass_tables, overwrite=True)
            print(
                f"Searching for crossphase tables: {caldir}/calibrator_{calibrator_obsid}*.kcrossscal"
            )
            crossphase_tables = sorted(
                glob.glob(f"{caldir}/calibrator_{calibrator_obsid}*.kcrosscal")
            )
            crossphase_tables = interpolate_bpass(crossphase_tables, overwrite=True)
            if len(bandpass_tables) == 0:
                print(
                    f"No bandpass table is present in calibration directory : {caldir}."
                )
                has_cal = False
                if emails != "":
                    email_msg = "No bandpass calibration table is found."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
            else:
                has_cal = True
                print("###################################################")
                print(f"Bandpass tables in calibration directory: {caldir}")
                for bpass in bandpass_tables:
                    print(f"{os.path.basename(bpass)}")
                print("####################################################")
                print(f"Crosshand phase tables in calibration directory: {caldir}")
                for kcross in crossphase_tables:
                    print(f"{os.path.basename(kcross)}")
                print("####################################################")
                caltables = bandpass_tables + crossphase_tables
        else:
            has_cal = False

        ###############################################
        # Making diagnostic plots
        ###############################################
        if has_cal and len(bandpass_tables) > 0 and do_basic_cal:
            os.makedirs(f"{outdir}/diagnostic_plots", exist_ok=True)
            msg, bpass_plots = plot_caltable_diagnostics(
                bandpass_tables, f"{outdir}/diagnostic_plots/{calibrator_obsid}_bcal"
            )
            if msg == 0:
                print(
                    f"Diagnostic plots for bandpass tables are saved in : {bpass_plots}."
                )
            else:
                print(f"Error in creating diagnostic plots for bandpass tables.")
        if has_cal and len(crossphase_tables) > 0 and do_basic_cal:
            os.makedirs(f"{outdir}/diagnostic_plots", exist_ok=True)
            msg, kcross_plots = plot_caltable_diagnostics(
                crossphase_tables,
                f"{outdir}/diagnostic_plots/{calibrator_obsid}_kcrosscal",
                plot_all_ants=False,
            )
            if msg == 0:
                print(
                    f"Diagnostic plots for crosshand phase tables are saved in : {kcross_plots}."
                )
            else:
                print(f"Error in creating diagnostic plots for crosshand phase tables.")

        ###################################################
        # Start spliting selfcal ms
        ###################################################
        if do_selfcal:
            ###############################################
            # Removing previous self-calibration artificats
            ###############################################
            print("Removing all previous self-calibration artificats...")
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

            if adaptive:
                scale_worker_and_wait(
                    dask_cluster, dask_client, min(total_ncoarse + 1, max_worker)
                )

            ######################
            # Spliting
            ######################
            if emails != "":
                email_msg = "Started spliting of measurement sets for self-calibration."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print(f"Starting task: Spliting {prefix} .....")
            print("###########################")
            future_selfcal_split = run_target_split_jobs.with_options(
                task_run_name=f"spliting_{prefix}_{jobid}"
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
                time_window=min(1.0, time_interval),
                time_interval=time_interval,
                quack_timestamps=quack_timestamps,
                jobid=jobid,
                cpu_frac=float(cpu_frac),
                mem_frac=float(mem_frac),
                remote_log=remote_logger,
            )
            print("Checking status of spliting of target for selfcal ...")
            try:
                msg, expected, succeed = future_selfcal_split.result()
                if emails != "":
                    email_msg = f"Spliting of measurement sets for self-calibration is done.\nExpected: {expected}, succeeded: {succeed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(
                    f"Finished task: Spliting of measurement sets for self-calibration is done."
                )
                print("###########################")
            except Exception:
                print(
                    "!!!! WARNING: Error in running spliting target scans for selfcal. !!!!"
                )
                do_selfcal = False
                traceback.print_exc()
                if emails != "":
                    email_msg = "Error occured in spliting target measurement sets for self-calibration."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
            finally:
                if adaptive:
                    scale_worker_and_wait(dask_cluster, dask_client, 1)

        ######################################
        # Checking status of self-cal split
        ######################################
        if do_selfcal:
            print("Checking measurement sets before spawning self-calibrations....")
            ####################################
            # Filtering any corrupted ms
            #####################################
            selfcal_target_mslist = sorted(glob.glob(workdir + "/selfcal*_spw_*.ms"))
            if (selfcal_target_mslist) == 0:
                print(
                    "!!!! WARNING: Error in running spliting target scans for selfcal. !!!!"
                )
                do_selfcal = False
                if emails != "":
                    email_msg = "No splited measurement set is found for self-calibration. Not continuting for self-calibration."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )

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
                print(
                    "No splited target scan ms are available in work directory for selfcal. Not continuing further for selfcal."
                )
                do_selfcal = False
                if emails != "":
                    email_msg = "No splited measurement set is found for self-calibration. Not continuting for self-calibration."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
            print(f"Selfcal mslist : {[os.path.basename(i) for i in selfcal_mslist]}")

        #########################################################
        # Flagging on targets for self-calibration
        #########################################################
        cal_applied = False
        if do_selfcal:
            if adaptive:
                scale_worker_and_wait(
                    dask_cluster, dask_client, min(len(selfcal_mslist) + 1, max_worker)
                )
            ############################
            # Basic flagging for selfcal
            ############################
            if emails != "":
                email_msg = "Started flagging for self-calibration measurment sets."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print("Starting task: Flagging selfcal targets .....")
            print("###########################")
            future_flag = run_flag.with_options(
                task_run_name=f"flagging_selfcal_{jobid}"
            ).submit(
                ",".join(selfcal_mslist),
                target_metafits,
                workdir,
                outdir,
                flag_calibrators=False,
                flag_quack=False,
                jobid=jobid,
                cpu_frac=round(cpu_frac, 2),
                mem_frac=round(mem_frac, 2),
                remote_log=remote_logger,
            )
            try:
                msg, succeed, failed = future_flag.result()
                if emails != "":
                    email_msg = f"Flagging for self-calibration measurment sets are done.\nSucceeded: {succeed}, failed: {failed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                for s_ms in selfcal_mslist:
                    s_ms = s_ms.rstrip("/")
                    if os.path.exists(f"{s_ms}/.flag_failed"):
                        print(
                            f"Issue in flagging of measurement set: {s_ms}. Check calibration solutions carefully."
                        )
                print("###########################")
                print(
                    f"Finished task: Flagging for self-calibration measurment sets are done."
                )
                print("###########################")
            except Exception:
                print(
                    "!!!! WARNING: Flagging error. Examine calibration solutions with caution. !!!!"
                )
                traceback.print_exc()
                if emails != "":
                    email_msg = (
                        "Error occured in flagging self-calibration measurement sets."
                    )
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )

            ###################################
            # Apply basic calibration
            ###################################
            if has_cal:
                if emails != "":
                    email_msg = "Started applying basic calibration solution on self-calibration measurement sets."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(
                    "Starting task: Applying basic calibration on self-calibration measurement sets....."
                )
                print("###########################")
                future_apply_basical_selfcal = run_apply_basiccal_sol.with_options(
                    task_run_name=f"applying_basiccal_selfcal_{jobid}"
                ).submit(
                    ",".join(selfcal_mslist),
                    calibrator_metafits,
                    target_metafits,
                    workdir,
                    caldir,
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
                        email_msg = f"Applying basic calibration solution on self-calibration measurement sets are done.\nSucceeded: {succeed}, failed: {failed}."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    print("###########################")
                    print(
                        f"Finished task: Applying basic calibration solution on self-calibration measurement sets are done."
                    )
                    print("###########################")
                except Exception:
                    print(
                        "!!!! WARNING: Error in applying basic calibration solutions on target. Continuing selfcal without basic calibration.!!!!"
                    )
                    traceback.print_exc()
                    cal_applied = False
                    do_applycal = False
                    if emails != "":
                        email_msg = "Error occured in applying basic calibration solutions on self-calibration measurement sets."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )

            if cal_applied:
                selfcal_applymode = "calonly"
                filtered_selfcalms_list = []
                for selfcalms in selfcal_mslist:
                    unflag_chans, flag_chans = get_chans_flag(
                        msname=selfcalms, n_threads=n_threads
                    )
                    if len(flag_chans) / (len(flag_chans) + len(unflag_chans)) <= 0.8:
                        filtered_selfcalms_list.append(selfcalms)
                    else:
                        print(
                            f"More than 80% channels are flagged for ms: {selfcalms}. Not using for self-calibration."
                        )
                selfcal_mslist = filtered_selfcalms_list
            else:
                selfcal_applymode = "calflag"

            ###############################################
            # Performing sidereal correction before selfcal
            ###############################################
            os.system(
                f"rm -rf {workdir}/*selfcal_int* {workdir}/*selfcal_pol* {workdir}/caltables/*selfcal*"
            )
            if adaptive:
                scale_worker_and_wait(
                    dask_cluster,
                    dask_client,
                    min(len(selfcal_mslist) + 1, max_worker),
                )
            if do_sidereal_cor:
                if emails != "":
                    email_msg = "Started correcting for solar sidereal motion."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(
                    "Starting task: Sidereal motion correction for self-calibration measurement sets....."
                )
                print("###########################")
                future_sidereal_cor_selfcal = run_solar_siderealcor_jobs.with_options(
                    task_run_name=f"solar_sidereal_correction_{jobid}"
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
                        email_msg = f"Correction for solar sidereal motion is done.\nSucceeded: {succeed}, failed: {failed}."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    print("###########################")
                    print(
                        f"Finished task: Correction for solar sidereal motion is done."
                    )
                    print("###########################")
                except Exception:
                    print("Sidereal correction is not successful.")
                    traceback.print_exc()
                    if emails != "":
                        email_msg = "Error occured in sidereal motion correction."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )

            #############################
            # Self-calibration
            #############################
            if emails != "":
                email_msg = "Started self-calibration."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print("Starting task: Self-calibrations.....")
            print("###########################")
            if cal_applied:
                print("Calibrator solutions are applied.")
            else:
                print("Calibration soplutions are not applied")
            future_selfcal = run_selfcal_jobs.with_options(
                task_run_name=f"selfcal_{jobid}"
            ).submit(
                ",".join(selfcal_mslist),
                workdir,
                caldir,
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
                jobid=jobid,
                cpu_frac=round(cpu_frac, 2),
                mem_frac=round(mem_frac, 2),
                remote_log=remote_logger,
            )
            try:
                msg, int_succeed, int_failed, pol_succeed, pol_failed = (
                    future_selfcal.result()
                )
                if emails != "":
                    email_msg = f"Self-calibration is done.\nIntensity self-calibration, Succeeded: {int_succeed}, failed: {int_failed}."
                    if do_polcal:
                        email_msg += f"\nPolarisation self-calibration, Succeeded: {pol_succeed}, failed: {pol_failed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(f"Finished task: Self-calibration is done.")
                print("###########################")
            except Exception:
                print(
                    "!!!! WARNING: Error in self-calibration on targets. Not applying self-calibration. !!!!"
                )
                do_apply_selfcal = False
                traceback.print_exc()
                if emails != "":
                    email_msg = "Error occured in self-calibration."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
            if adaptive:
                scale_worker_and_wait(dask_cluster, dask_client, 1)

        ########################################
        # Checking self-cal caltables
        ########################################
        selfcal_gaincal = sorted(glob.glob(f"{caldir}/selfcal_{target_obsid}*.gcal"))
        selfcal_bandpass = sorted(glob.glob(f"{caldir}/selfcal_{target_obsid}*.bcal"))
        selfcal_bandpass = interpolate_bpass(selfcal_bandpass, overwrite=True)
        if do_polcal:
            selfcal_leakages = sorted(
                glob.glob(f"{caldir}/selfcal_{target_obsid}*.dcal")
            )
            selfcal_leakages = interpolate_quartical(selfcal_leakages, overwrite=True)

        if len(selfcal_gaincal) == 0:
            print(
                "Self-calibration is not performed and no self-calibration caltable is available."
            )
            do_apply_selfcal = False
            if emails != "":
                email_msg = "Self-calibration is not performed and no self-calibration caltable is available."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )

        ###########################################
        # Plotting self-caltables
        ###########################################
        if do_selfcal and len(selfcal_gaincal) > 0:
            os.makedirs(f"{outdir}/diagnostic_plots", exist_ok=True)
            msg, gcal_plots = plot_caltable_diagnostics(
                selfcal_gaincal, f"{outdir}/diagnostic_plots/{target_obsid}_gcal"
            )
            if msg == 0:
                print(
                    f"Diagnostic plots for self-calibration gaincal tables are saved in : {gcal_plots}."
                )
            else:
                print(
                    f"Error in creating diagnostic plots for self-calibration gaincal tables."
                )

        if do_selfcal and len(selfcal_bandpass) > 0:
            os.makedirs(f"{outdir}/diagnostic_plots", exist_ok=True)
            msg, bcal_plots = plot_caltable_diagnostics(
                selfcal_bandpass, f"{outdir}/diagnostic_plots/{target_obsid}_bcal"
            )
            if msg == 0:
                print(
                    f"Diagnostic plots for self-calibration bandpass tables are saved in : {bcal_plots}."
                )
            else:
                print(
                    f"Error in creating diagnostic plots for self-calibration bandpass tables."
                )

        if do_selfcal and do_polcal and len(selfcal_leakages) > 0:
            os.makedirs(f"{outdir}/diagnostic_plots", exist_ok=True)
            msg, dcal_plots = plot_quartical_tables(
                selfcal_leakages, f"{outdir}/diagnostic_plots/{target_obsid}_dcal"
            )
            if msg == 0:
                print(
                    f"Diagnostic plots for self-calibration leakage tables are saved in : {dcal_plots}."
                )
            else:
                print(
                    f"Error in creating diagnostic plots for self-calibration leakage tables."
                )

        #############################################
        # Spliting targets if not started already
        #############################################
        # If corrected data is requested or imaging is requested
        if do_applycal or do_apply_selfcal or do_imaging:
            if adaptive:
                scale_worker_and_wait(
                    dask_cluster, dask_client, min(total_ncoarse + 1, max_worker)
                )
            prefix = "target"
            if emails != "":
                email_msg = "Started spliting target for final processing."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print(f"Starting task: Spliting {prefix} .....")
            print("###########################")
            future_split = run_target_split_jobs.with_options(
                task_run_name=f"spliting_{prefix}_{jobid}"
            ).submit(
                ",".join(target_mslist),
                target_metafits,
                workdir,
                datacolumn="data",
                force_split=True,
                freqres=freqavg,
                timeres=timeavg,
                quack_timestamps=quack_timestamps,
                prefix=prefix,
                jobid=jobid,
                cpu_frac=round(cpu_frac, 2),
                mem_frac=round(mem_frac, 2),
                remote_log=remote_logger,
            )
            ##########################################
            # Checking target spliting is done or not
            ##########################################
            print("Checking spliting of targets status...")
            try:
                msg, expected, succeed = future_split.result()
                if emails != "":
                    email_msg = f"Spliting target for final processing is done.\nExpected: {expected}, succeeded: {succeed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(f"Finished task: Spliting target for final processing is done.")
                print("###########################")
            except Exception:
                print("!!!! WARNING: Error in spliting targets. !!!!")
                traceback.print_exc()
                if emails != "":
                    email_msg = "Error occured in spliting target for final processing."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                return 1
            finally:
                if adaptive:
                    scale_worker_and_wait(dask_cluster, dask_client, 1)

            ################################
            # Checking splited final ms list
            ################################
            split_target_mslist = sorted(glob.glob(workdir + "/target*_spw_*.ms"))
            if len(split_target_mslist) == 0:
                print("!!!! WARNING: No target ms are present. !!!!")
                if emails != "":
                    email_msg = (
                        "No target measurement set is present for final processing."
                    )
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                return 1
            ####################################
            # Filtering any corrupted ms
            #####################################
            print(
                "Checking final valid measurement sets before applying solutions and spawning imaging...."
            )
            filtered_mslist = []  # Filtering in case any ms is corrupted
            for ms in split_target_mslist:
                checkcol = check_datacolumn_valid(ms)
                if checkcol:
                    filtered_mslist.append(ms)
                else:
                    print(f"Issue in : {ms}")
                    os.system(f"rm -rf {ms}")
            split_target_mslist = filtered_mslist
            if len(split_target_mslist) == 0:
                print("No filtered target ms are available in work directory.")
                if emails != "":
                    email_msg = "No un-corrupted target measurement is present for final processing."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                return 1
            print(
                f"Target mslist : {[os.path.basename(i) for i in split_target_mslist]}"
            )

            #########################################################
            # Applying basic solutions on target scans
            #########################################################
            if adaptive:
                scale_worker_and_wait(
                    dask_cluster,
                    dask_client,
                    min(len(split_target_mslist) + 1, max_worker),
                )
            ############################
            # Basic flagging
            ############################
            if emails != "":
                email_msg = "Started flagging of final target measurement sets."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print("Starting task: Flagging final target measurement sets .....")
            print("###########################")
            future_flag = run_flag.with_options(
                task_run_name=f"flagging_target_{jobid}"
            ).submit(
                ",".join(split_target_mslist),
                target_metafits,
                workdir,
                outdir,
                flag_calibrators=False,
                flag_quack=False,
                jobid=jobid,
                cpu_frac=round(cpu_frac, 2),
                mem_frac=round(mem_frac, 2),
                remote_log=remote_logger,
            )
            try:
                msg, succeed, failed = future_flag.result()
                if emails != "":
                    email_msg = f"Flagging of final target measurement sets are done.\nSucceeded: {succeed}, failed: {failed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(
                    f"Finished task: Flagging of final target measurement sets are done."
                )
                print("###########################")
            except Exception:
                print(
                    "!!!! WARNING: Flagging error. Examine calibration solutions with caution. !!!!"
                )
                traceback.print_exc()
                if emails != "":
                    email_msg = (
                        "Error occured in flagging of final target measurement sets."
                    )
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )

            ####################################
            # Applying basic calibration
            #####################################
            if (do_applycal or do_apply_selfcal) and has_cal:
                if emails != "":
                    email_msg = "Started applying basic calibration solutions on final target measurement sets."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(
                    "Starting task: Applying basic calibration on final target measurement sets....."
                )
                print("###########################")
                future_apply_basical = run_apply_basiccal_sol.with_options(
                    task_run_name=f"applying_basiccal_target_{jobid}"
                ).submit(
                    ",".join(split_target_mslist),
                    calibrator_metafits,
                    target_metafits,
                    workdir,
                    caldir,
                    overwrite_datacolumn=True,
                    only_amplitude=only_amplitude,
                    applymode="calflag",
                    prefix="target",
                    jobid=jobid,
                    cpu_frac=round(cpu_frac, 2),
                    mem_frac=round(mem_frac, 2),
                    remote_log=remote_logger,
                )
                try:
                    msg, succeed, failed = future_apply_basical.result()
                    if emails != "":
                        email_msg = f"Applying basic calibration solutions on final target measurement sets are done.\nSucceeded: {succeed}, failed: {failed}."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    print("###########################")
                    print(
                        f"Finished task: Applying basic calibration solutions on final target measurement sets are done."
                    )
                    print("###########################")
                except Exception:
                    print(
                        "!!!! WARNING: Error in applying basic calibration solutions on target scans. Not continuing further.!!!!"
                    )
                    traceback.print_exc()
                    if emails != "":
                        email_msg = "Error occured in applying basic calibration on final target measurement sets. P-AIRCARS has stopped."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    return 1

            ###################################
            # Correct sidereal motion
            ###################################
            if do_sidereal_cor:
                if emails != "":
                    email_msg = "Start correcting sidereal motion of the Sun on final target measurement sets."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(
                    "Starting task: Sidereal motion correction for final target measurement sets....."
                )
                print("###########################")
                future_sidereal_cor = run_solar_siderealcor_jobs.with_options(
                    task_run_name=f"solar_sidereal_correction_{jobid}"
                ).submit(
                    ",".join(split_target_mslist),
                    workdir,
                    prefix="target",
                    jobid=jobid,
                    cpu_frac=round(cpu_frac, 2),
                    mem_frac=round(mem_frac, 2),
                    remote_log=remote_logger,
                )
                try:
                    msg, succeed, failed = future_sidereal_cor.result()
                    if emails != "":
                        email_msg = f"Sidereal motion correction of the Sun on final target measurement sets are done.\nSucceeded: {succeed}, failed: {failed}."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    print("###########################")
                    print(
                        f"Finished task: Sidereal motion correction of the Sun on final target measurement sets are done."
                    )
                    print("###########################")
                except Exception:
                    print("!!!! WARNING: Error in applying sidereal correction.!!!!")
                    traceback.print_exc()
                    if emails != "":
                        email_msg = "Error occured in sidereal motion correction on final target measurement sets."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )

            ########################################
            # Apply self-calibration
            ########################################
            if do_apply_selfcal:
                selfcal_applymode = "calonly"
                for msname in split_target_mslist:
                    if not os.path.exists(f"{msname}/.applied_sol"):
                        selfcal_applymode = "calflag"

                if emails != "":
                    email_msg = "Started applying self-calibration on final target measurement sets."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(
                    "Starting task: Applying self-calibration solutions on final target measurement sets....."
                )
                print("###########################")
                future_apply_selfcal = run_apply_selfcal_sol.with_options(
                    task_run_name=f"applying_selfcal_{jobid}"
                ).submit(
                    ",".join(split_target_mslist),
                    target_metafits,
                    workdir,
                    caldir,
                    overwrite_datacolumn=False,
                    applymode=selfcal_applymode,
                    jobid=jobid,
                    cpu_frac=round(cpu_frac, 2),
                    mem_frac=round(mem_frac, 2),
                    remote_log=remote_logger,
                )
                try:
                    msg, gain_succeed, gain_failed, pol_succeed, pol_failed = (
                        future_apply_selfcal.result()
                    )
                    if emails != "":
                        email_msg = f"Applying self-calibration on final target measurement sets are done.\nGain solutions applied: Succeeded: {gain_succeed}, failed: {gain_failed}."
                        if do_polcal:
                            email_msg += f"\nPolarisation solution applied: Succeeded: {pol_succeed}, failed: {pol_failed}."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    print("###########################")
                    print(
                        f"Finished task: Applying self-calibration on final target measurement sets are done."
                    )
                    print("###########################")
                except Exception:
                    print(
                        "!!!! WARNING: Error in applying self-calibration solutions on targets. !!!!"
                    )
                    traceback.print_exc()
                    if emails != "":
                        email_msg = "Error occured in applying self-calibration solutions on final target measurement sets."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )

            ######################################
            # Imaging
            ######################################
            if do_imaging:
                if image_freqres > 0:
                    print(f"Image frequency resolution: {image_freqres} MHz.")
                else:
                    print(f"Image frequency resolution: entire corase channel.")
                if image_timeres > 0:
                    print(f"Image time resolution: {image_timeres} s.")
                else:
                    print("Imaging entire scan.")
                pol = pol.upper()
                if pol not in ["I", "IQUV"]:
                    pol = "IQUV"

                if (
                    do_polcal == False
                ):  # Only if do_polcal is False, overwrite to make only Stokes I
                    pol = "I"

                if emails != "":
                    email_msg = "Started final imaging."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print("Starting task: Final imaging.....")
                print("###########################")
                future_imaging = run_imaging_jobs.with_options(
                    task_run_name=f"imaging_{jobid}"
                ).submit(
                    ",".join(split_target_mslist),
                    workdir,
                    outdir,
                    freqrange=freqrange,
                    timerange=timerange,
                    minuv=minuv,
                    weight=weight,
                    robust=float(robust),
                    pol=pol,
                    freqres=image_freqres,
                    timeres=image_timeres,
                    threshold=float(clean_threshold),
                    use_multiscale=use_multiscale,
                    use_solar_mask=use_solar_mask,
                    cutout_rsun=cutout_rsun,
                    savemodel=keep_backup,
                    saveres=keep_backup,
                    jobid=jobid,
                    cpu_frac=round(cpu_frac, 2),
                    mem_frac=round(mem_frac, 2),
                    remote_log=remote_logger,
                )
                try:
                    msg, succeed, failed, total_images = future_imaging.result()
                    if emails != "":
                        email_msg = f"Final imaging is done.\nSucceeded: {succeed}, failed: {failed}.\nTotal images made: {total_images}."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    print("###########################")
                    print(f"Finished task: Final imaging is done.")
                    print("###########################")
                except Exception:
                    print(
                        "!!!! WARNING: Final imaging on all measurement sets is not successful. Check the image directory. !!!!"
                    )
                    traceback.print_exc()
                    if emails != "":
                        email_msg = (
                            "Error occured in final imaging. P-AIRCARS has stopped."
                        )
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    return 1
            if adaptive:
                scale_worker_and_wait(dask_cluster, dask_client, 1)

        ########################################
        # Naming of image directory
        ########################################
        if weight == "briggs":
            weight_str = f"{weight}_{robust}"
        else:
            weight_str = weight
        if image_freqres == -1 and image_timeres == -1:
            imagedir = outdir + f"/imagedir_f_all_t_all_pol_{pol}_w_{weight_str}"
        elif image_freqres != -1 and image_timeres == -1:
            imagedir = (
                outdir + f"/imagedir_f_{image_freqres}_t_all_pol_{pol}_w_{weight_str}"
            )
        elif image_freqres == -1 and image_timeres != -1:
            imagedir = (
                outdir + f"/imagedir_f_all_t_{image_timeres}_pol_{pol}_w_{weight_str}"
            )
        else:
            imagedir = (
                outdir
                + f"/imagedir_f_{image_freqres}_t_{image_timeres}_pol_{pol}_w_{weight_str}"
            )

        os.system(f"rm -rf {imagedir}/images/*aia*.fits")
        os.system(f"rm -rf {imagedir}/images/*suvi*.fits")
        ###########################
        # Primary beam correction
        ###########################
        if do_pbcor:
            images = sorted(glob.glob(f"{imagedir}/images/*.fits"))
            if len(images) == 0:
                print(f"No image is present in image directory: {imagedir}/images")
                if emails != "":
                    email_msg = "No image is present in image directory for primary beam correction."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
            else:
                if adaptive:
                    scale_worker_and_wait(
                        dask_cluster, dask_client, min(len(images) + 1, max_worker)
                    )
                if emails != "":
                    email_msg = "Started primary beam correction."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print("Starting task: Primary beam correction.....")
                print("###########################")
                future_pbcor = run_apply_pbcor.with_options(
                    task_run_name=f"applying_primary_beam_{jobid}"
                ).submit(
                    f"{imagedir}/images",
                    target_metafits,
                    workdir,
                    leakage_dir=caldir,
                    jobid=jobid,
                    cpu_frac=round(cpu_frac, 2),
                    mem_frac=round(mem_frac, 2),
                    remote_log=remote_logger,
                )
                try:
                    msg, succeed, failed = future_pbcor.result()
                    if emails != "":
                        email_msg = f"Primary beam correction is done.\nSucceeded: {succeed}, failed: {failed}."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    print("###########################")
                    print(f"Finished task: Primary beam correction is done.")
                    print(f"Final image directory: {imagedir}/images")
                    print("###########################")
                except Exception:
                    print(
                        "!!!! WARNING: Primary beam corrections of the final images are not successful. !!!!"
                    )
                    traceback.print_exc()
                    if emails != "":
                        email_msg = "Error occured in primary beam correction. P-AIRCARS has stopped."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    return 1
                finally:
                    if adaptive:
                        scale_worker_and_wait(dask_cluster, dask_client, 1)

        ##############################################
        # Making diagnostic plots of measurement sets
        ##############################################
        if make_msplot:
            ###########################################
            # Ploting calibrator ms
            ###########################################
            split_cal_mslist = sorted(glob.glob(f"{workdir}/calibrator*_spw_*.ms"))
            if len(split_cal_mslist) == 0:
                print("No calibrator measurement set is present for ploting.")
            else:
                if adaptive:
                    scale_worker_and_wait(
                        dask_cluster,
                        dask_client,
                        min(len(split_cal_mslist) + 1, max_worker),
                    )
                msplot_outdir = f"{outdir}/ms_diagnostics_plots"
                os.makedirs(msplot_outdir, exist_ok=True)
                if len(split_cal_mslist) > 0:
                    if emails != "":
                        email_msg = "Started making diagnostic plots for calibrator measurement sets."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    print("###########################")
                    print(
                        "Starting task: Making diagnostic plots of calibrator measurement sets....."
                    )
                    print("###########################")
                    future_cal_plot = run_make_msplot.with_options(
                        task_run_name=f"making_msplot_cal_{jobid}"
                    ).submit(
                        ",".join(split_cal_mslist),
                        workdir,
                        msplot_outdir,
                        jobid=jobid,
                        cpu_frac=round(cpu_frac, 2),
                        mem_frac=round(mem_frac, 2),
                        remote_log=remote_logger,
                    )
                    try:
                        msg = future_cal_plot.result()
                        if emails != "":
                            email_msg = "Making diagnostic plots for calibrator measurement sets are done."
                            send_task_notification(
                                emails, email_msg, jobid, target_obsid, timestamp
                            )
                        print("###########################")
                        print(
                            f"Finished task: Making diagnostic plots for calibrator measurment sets are done."
                        )
                        print("###########################")
                    except Exception:
                        print(
                            "!!!! WARNING: Diagnostic plot of calibrator measurment sets are not successful. !!!!"
                        )
                        traceback.print_exc()
                        if emails != "":
                            email_msg = "Error occured in making diagnostic plots of calibrator measurement sets."
                            send_task_notification(
                                emails, email_msg, jobid, target_obsid, timestamp
                            )

            ###########################################
            # Ploting target ms
            ###########################################
            split_target_mslist = sorted(glob.glob(f"{workdir}/target*_spw_*.ms"))
            if len(split_target_mslist) == 0:
                print("No target measurment set is present for ploting.")
            else:
                if adaptive:
                    scale_worker_and_wait(
                        dask_cluster,
                        dask_client,
                        min(len(split_target_mslist) + 1, max_worker),
                    )
                msplot_outdir = f"{outdir}/ms_diagnostics_plots"
                os.makedirs(msplot_outdir, exist_ok=True)
                print("###########################")
                print(
                    "Starting task: Making diagnostic plots of target measurement sets....."
                )
                print("###########################")
                future_target_plot = run_make_msplot.with_options(
                    task_run_name=f"making_msplot_{jobid}"
                ).submit(
                    ",".join(split_target_mslist),
                    workdir,
                    msplot_outdir,
                    jobid=jobid,
                    cpu_frac=round(cpu_frac, 2),
                    mem_frac=round(mem_frac, 2),
                    remote_log=remote_logger,
                )
                try:
                    msg = future_target_plot.result()
                    if emails != "":
                        email_msg = "Making diagnostic plots for target measurement sets are done."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )
                    print("###########################")
                    print(
                        f"Finished task: Making diagnostic plots for target measurment sets are done."
                    )
                    print("###########################")
                except Exception:
                    print(
                        "!!!! WARNING: Diagnostic plot of target measurment sets are not successful. !!!!"
                    )
                    traceback.print_exc()
                    if emails != "":
                        email_msg = "Error occured in making diagnostic plots of target measurement sets."
                        send_task_notification(
                            emails, email_msg, jobid, target_obsid, timestamp
                        )

        scale_worker_and_wait(dask_cluster, dask_client, 1)
        #######################################
        # Make overlays
        #######################################
        images = sorted(glob.glob(f"{imagedir}/images/*.fits"))
        if len(images) == 0:
            print(
                f"No image is present in image directory: {imagedir}/images for making overlays"
            )
            if emails != "":
                email_msg = (
                    "No image is present in image directory for making overlays."
                )
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )

        #################################################################
        # Filtering only coarse channel images for default overlay mode
        #################################################################
        if make_overlay is False:
            images = filter_images(images, min_time_sep=60.0)
        if len(images) > 0:
            if adaptive:
                scale_worker_and_wait(
                    dask_cluster, dask_client, min(len(images) + 1, max_worker)
                )
            #################################
            # Start overlays
            #################################
            if emails != "":
                email_msg = "Started making overlays."
                send_task_notification(
                    emails, email_msg, jobid, target_obsid, timestamp
                )
            print("###########################")
            print("Starting task: Making overlay on EUV images.....")
            print("###########################")
            future_overlay = run_make_overlay.with_options(
                task_run_name=f"making_overlay_{jobid}"
            ).submit(
                f"{imagedir}/images",
                f"{imagedir}/overlay_pngs",
                workdir=workdir,
                all_overlay=make_overlay,
                jobid=jobid,
                cpu_frac=round(cpu_frac, 2),
                remote_log=remote_logger,
            )
            try:
                msg, succeed, failed = future_overlay.result()
                if emails != "":
                    email_msg = f"Making overlays are done.\nSucceeded: {succeed}, failed: {failed}."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )
                print("###########################")
                print(f"Finished task: Making overlays are done.")
                print(f"Final image directory: {imagedir}/overlay_pngs")
                print("###########################")
            except Exception:
                print("!!!! WARNING: Overlay of the images are not successful. !!!!")
                traceback.print_exc()
                if emails != "":
                    email_msg = "Error occured in making overlays."
                    send_task_notification(
                        emails, email_msg, jobid, target_obsid, timestamp
                    )

        ###########################################
        # Successful exit
        ###########################################
        print(f"P-AIRCARS calibration and imaging pipeline is successfully executed.")
        if emails != "":
            email_msg = "P-AIRCARS processing is done successfully."
            send_task_notification(emails, email_msg, jobid, target_obsid, timestamp)
        return 0
    except Exception as e:
        traceback.print_exc()
        if emails != "":
            email_msg = f"Error in running P-AIRCARS.\n{e}"
            send_task_notification(emails, email_msg, jobid, target_obsid, timestamp)
        print("###########################")
        print(f"Error occured in running P-AIRCARS.")
        print("###########################")
        return 1
    finally:
        time.sleep(5)
        datalist = sorted(glob.glob(f"{target_datadir}/*"))
        for data in datalist:
            drop_cache(data)
        callist = sorted(glob.glob(f"{calibrator_datadir}/*"))
        for cal in callist:
            drop_cache(cal)
        ######################################
        # Keeping flag backups
        ######################################
        # Flag backups of calibrator measurement sets
        ######################################
        final_cal_mslist = sorted(glob.glob(workdir + "/calibrator*_spw_*.ms"))
        if len(final_cal_mslist) > 0:
            os.makedirs(f"{outdir}/ms_flags", exist_ok=True)
            print(f"Doing flag backup in: {outdir}/ms_flags")
            for cal_ms in final_cal_mslist:
                do_flag_backup(cal_ms, flagtype="finalflag")
                if os.path.exists(
                    f"{outdir}/ms_flags/{os.path.basename(cal_ms)}.flagversions"
                ):
                    os.system(
                        f"rm -rf {outdir}/ms_flags/{os.path.basename(cal_ms)}.flagversions"
                    )
                os.system(f"mv {cal_ms}.flagversions {outdir}/ms_flags/")
                if keep_backup is False:
                    os.system(f"rm -rf {cal_ms}")
        ######################################
        # Flag backups of selfcal measurement sets
        ######################################
        final_selfcal_mslist = sorted(glob.glob(workdir + "/selfcal*_spw_*.ms"))
        if len(final_selfcal_mslist) > 0:
            os.makedirs(f"{outdir}/ms_flags", exist_ok=True)
            print(f"Doing flag backup in: {outdir}/ms_flags")
            for selfcal_ms in final_selfcal_mslist:
                do_flag_backup(selfcal_ms, flagtype="finalflag")
                if os.path.exists(
                    f"{outdir}/ms_flags/{os.path.basename(selfcal_ms)}.flagversions"
                ):
                    os.system(
                        f"rm -rf {outdir}/ms_flags/{os.path.basename(selfcal_ms)}.flagversions"
                    )
                os.system(f"mv {selfcal_ms}.flagversions {outdir}/ms_flags/")
                if keep_backup is False:
                    os.system(f"rm -rf {selfcal_ms}")
        ######################################
        # Flag backups of target measurement sets
        ######################################
        final_split_target_mslist = sorted(glob.glob(workdir + "/target*_spw_*.ms"))
        if len(final_split_target_mslist) > 0:
            os.makedirs(f"{outdir}/ms_flags", exist_ok=True)
            print(f"Doing flag backup in: {outdir}/ms_flags")
            for target_ms in final_split_target_mslist:
                do_flag_backup(target_ms, flagtype="finalflag")
                if os.path.exists(
                    f"{outdir}/ms_flags/{os.path.basename(target_ms)}.flagversions"
                ):
                    os.system(
                        f"rm -rf {outdir}/ms_flags/{os.path.basename(target_ms)}.flagversions"
                    )
                os.system(f"mv {target_ms}.flagversions {outdir}/ms_flags/")
                if keep_calibrated_ms:
                    calibrated_msdir = f"{outdir}/calibrated_ms"
                    os.makedirs(calibrated_msdir, exist_ok=True)
                    os.system(f"mv {target_ms} {calibrated_msdir}")
                elif keep_backup is False:
                    os.system(f"rm -rf {target_ms}")
        time.sleep(5)
        drop_cache(workdir)
        drop_cache(outdir)
        if master_log_created:
            stop_event.set()
            log_thread_flow.join(timeout=5)
        if dask_dir is not None:
            os.system(f"rm -rf {dask_dir}")
        if observer is not None:
            clean_shutdown(observer)


def cli():
    parser = argparse.ArgumentParser(
        description="Run P-AIRCARS for calibration and imaging of solar observations.",
        formatter_class=SmartDefaultsHelpFormatter,
    )
    # === Essential parameters ===
    essential = parser.add_argument_group(
        "###################\nEssential parameters\n###################"
    )
    essential.add_argument(
        "target_datadir", type=str, help="Target measurement set directory"
    )
    essential.add_argument("target_metafits", type=str, help="Target metafits file")
    essential.add_argument(
        "--workdir",
        type=str,
        dest="workdir",
        required=True,
        help="Working directory",
    )
    essential.add_argument(
        "--outdir",
        type=str,
        dest="outdir",
        required=True,
        help="Output products directory",
    )
    essential.add_argument(
        "--cal_datadir",
        type=str,
        dest="cal_datadir",
        help="Calibrator measurement set directory",
    )
    essential.add_argument(
        "--cal_metafits",
        type=str,
        dest="cal_metafits",
        help="Calibrator metafits file",
    )

    # === Advanced calibration parameters ===
    advanced_cal = parser.add_argument_group(
        "###################\nAdvanced calibration parameters\n###################"
    )
    advanced_cal.add_argument(
        "--solint",
        type=str,
        default="60s",
        help="Solution interval for calibration (e.g. 'int', '10s', '5min', 'inf')",
    )
    advanced_cal.add_argument(
        "--cal_uvrange",
        type=str,
        default="",
        help="UV range to filter data for calibration (e.g. '>100klambda', '100~10000lambda')",
    )
    advanced_cal.add_argument(
        "--no_polcal",
        action="store_false",
        dest="do_polcal",
        help="Disable polarization calibration",
    )
    advanced_cal.add_argument(
        "--only_amplitude",
        action="store_true",
        help="Apply only amplitude part of gain solution from calibrator or not",
    )

    # === Advanced imaging parameters ===
    advanced_image = parser.add_argument_group(
        "###################\nAdvanced imaging parameters\n###################"
    )
    advanced_image.add_argument(
        "--freqrange",
        type=str,
        default="",
        help="Frequency range in MHz to select during imaging (comma-seperate, e.g. '100~110,130~140')",
    )
    advanced_image.add_argument(
        "--timerange",
        type=str,
        default="",
        help="Time range to select during imaging (comma-seperated, e.g. '2014/09/06/09:30:00~2014/09/06/09:45:00,2014/09/06/10:30:00~2014/09/06/10:45:00')",
    )
    advanced_image.add_argument(
        "--image_freqres",
        type=float,
        default=1.28,
        help="Output image frequency resolution in MHz (-1 = full)",
    )
    advanced_image.add_argument(
        "--image_timeres",
        type=float,
        default=10.0,
        help="Output image time resolution in seconds (-1 = full)",
    )
    advanced_image.add_argument(
        "--pol",
        type=str,
        default="IQUV",
        help="Stokes parameter(s) to image ('I' or 'IQUV')",
    )
    advanced_image.add_argument(
        "--minuv",
        type=float,
        default=0,
        help="Minimum baseline length (in wavelengths) to include in imaging",
    )
    advanced_image.add_argument(
        "--weight",
        type=str,
        default="briggs",
        help="Imaging weighting scheme (e.g. 'briggs', 'natural', 'uniform')",
    )
    advanced_image.add_argument(
        "--robust",
        type=float,
        default=0.0,
        help="Robust parameter for Briggs weighting (-2 to +2)",
    )
    advanced_image.add_argument(
        "--no_multiscale",
        action="store_false",
        dest="use_multiscale",
        help="Disable multiscale CLEAN for extended structures",
    )
    advanced_image.add_argument(
        "--clean_threshold",
        type=float,
        default=1.0,
        help="Clean threshold in sigma for final deconvolution",
    )
    advanced_image.add_argument(
        "--no_pbcor",
        action="store_false",
        dest="do_pbcor",
        help="Do not apply primary beam correction after imaging",
    )
    advanced_image.add_argument(
        "--cutout_rsun",
        type=float,
        default=10.0,
        help="Field of view cutout radius in solar radii",
    )
    advanced_image.add_argument(
        "--no_solar_mask",
        action="store_false",
        dest="use_solar_mask",
        help="Disable use solar disk mask during deconvolution",
    )
    advanced_image.add_argument(
        "--do_overlay",
        action="store_true",
        dest="make_overlay",
        help="Make overlay plot on EUV images for all images (default is to make overlays only one image per coarse channels at 10s intervals)",
    )
    advanced_image.add_argument(
        "--make_msplot",
        action="store_true",
        help="Make diagnostic plots of measurement sets",
    )

    # === Advanced options ===
    advanced = parser.add_argument_group(
        "###################\nAdvanced pipeline parameters\n###################"
    )
    advanced.add_argument(
        "--non_solar_data",
        action="store_false",
        dest="solar_data",
        help="Disable solar data mode",
    )
    advanced.add_argument(
        "--no_ds",
        action="store_false",
        dest="make_ds",
        help="Disable making solar dynamic spectra",
    )
    advanced.add_argument(
        "--do_forcereset_weightflag",
        action="store_true",
        help="Force reset of weights and flags (disabled by default)",
    )
    advanced.add_argument(
        "--no_cal_flag",
        action="store_false",
        dest="do_cal_flag",
        help="Disable initial flagging of calibrators",
    )
    advanced.add_argument(
        "--no_import_model",
        action="store_false",
        dest="do_import_model",
        help="Disable model import",
    )
    advanced.add_argument(
        "--no_basic_cal",
        action="store_false",
        dest="do_basic_cal",
        help="Disable basic gain calibration",
    )
    advanced.add_argument(
        "--do_sidereal_cor",
        action="store_true",
        dest="do_sidereal_cor",
        help="Sidereal motion correction for Sun (disabled by default)",
    )
    advanced.add_argument(
        "--no_solarcenter_move",
        action="store_false",
        dest="do_move_solarcenter",
        help="Disable moving phaseceneter to solar center",
    )
    advanced.add_argument(
        "--no_selfcal",
        action="store_false",
        dest="do_selfcal",
        help="Disable self-calibration",
    )
    advanced.add_argument(
        "--no_ap_selfcal",
        action="store_false",
        dest="do_ap_selfcal",
        help="Disable amplitude-phase self-calibration",
    )
    advanced.add_argument(
        "--no_solar_selfcal",
        action="store_false",
        dest="solar_selfcal",
        help="Disable solar-specific self-calibration parameters",
    )
    advanced.add_argument(
        "--no_applycal",
        action="store_false",
        dest="do_applycal",
        help="Disable application of basic calibration solutions",
    )
    advanced.add_argument(
        "--no_apply_selfcal",
        action="store_false",
        dest="do_apply_selfcal",
        help="Disable application of self-calibration solutions",
    )
    advanced.add_argument(
        "--no_imaging",
        action="store_false",
        dest="do_imaging",
        help="Disable final imaging",
    )

    # === Advanced local system/ per node hardware resource parameters ===
    advanced_resource = parser.add_argument_group(
        "###################\nAdvanced hardware resource parameters for local system or per node on HPC cluster\n###################"
    )
    advanced_resource.add_argument(
        "--cpu_frac",
        type=float,
        default=0.8,
        help="Fraction of CPU usuage per node",
    )
    advanced_resource.add_argument(
        "--mem_frac",
        type=float,
        default=0.8,
        help="Fraction of memory usuage per node",
    )
    advanced_resource.add_argument(
        "--max_worker",
        type=int,
        default=None,
        help="Maximum number of workers",
    )
    advanced_resource.add_argument(
        "--keep_backup",
        action="store_true",
        help="Keep backup of intermediate steps",
    )
    advanced_resource.add_argument(
        "--no_calibrated_ms",
        action="store_false",
        dest="keep_calibrated_ms",
        help="Keep calibrated measurement sets or not",
    )
    advanced_resource.add_argument(
        "--no_remote_logger",
        action="store_false",
        dest="remote_logger",
        help="Disable remote logger",
    )
    advanced_resource.add_argument(
        "--jobid",
        type=int,
        default=None,
        help="User provided P-AIRCARS job ID",
    )
    advanced_resource.add_argument(
        "--job_password",
        type=str,
        default=None,
        help="User specified job password",
    )
    advanced_resource.add_argument(
        "--masterlog",
        type=str,
        default=None,
        help="Master log file",
    )
    advanced_resource.add_argument(
        "--cluster",
        action="store_true",
        dest="cluster",
        help="Running in cluster environment",
    )
    advanced_resource.add_argument(
        "--adaptive",
        action="store_true",
        dest="adaptive",
        help="Whether do adaptive scaling or not in cluster environment (In local environment, it is always true)",
    )

    # === Advanced job scheduler parameters ===
    advanced_slurm = parser.add_argument_group(
        "###################\nAdvanced slurm cluster settings\n###################"
    )
    advanced_slurm.add_argument(
        "--partition",
        type=str,
        default=None,
        help="Partition name (Required)",
    )
    advanced_slurm.add_argument(
        "--account",
        type=str,
        default=None,
        help="Account name (If your cluster requires this, you should provide. Otherwise job can not be started)",
    )
    advanced_slurm.add_argument(
        "--walltime",
        type=str,
        default=None,
        help="Wall time, each slurm job can execute in maximum this time",
    )
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    scheduler_name = get_scheduler_name()
    ##################################################
    # Prefect settings check
    ##################################################
    prefect_settings = get_current_settings()
    prefect_env = prefect_settings.to_environment_variables()
    api_url = prefect_env.get("PREFECT_API_URL")
    cachedir = f"{get_cachedir()}/prefect_{scheduler_name}"
    config_file = f"{cachedir}/prefect.config.npy"

    ######################################
    # Check connection to prefect server
    ######################################
    if api_url is not None:
        check_url = f"{api_url}/health"
        try:
            r = requests.get(check_url, timeout=60)
        except Exception:
            traceback.print_exc()
            print(f"Could not reach prefect server at: {api_url} from compute node.")
            return 1
    else:
        print("Prefecr server is not running.")
        if scheduler_name == "local":
            print("P-AIRCARS will use ephmeral mode.")
        else:
            print(
                f"First start prefect server to run P-AIRCARS in {scheduler_name} cluster."
            )
            return 1

    f = Figlet(font="big")
    print(f.renderText("P-AIRCARS"))

    if args.jobid is None:
        jobid = get_jobid()
    else:
        jobid = args.jobid

    os.system(f"rm -rf {args.workdir}/dask_*")

    ###########################################################
    # Estimating jobs memory size (5 times each measurment set)
    ###########################################################
    if os.path.exists(args.target_datadir) is False:
        print(f"Target data directory: {args.target_datadir} does not exist.")
        return
    target_mslist = glob.glob(f"{args.target_datadir}/*.ms")
    if len(target_mslist) == 0:
        print(
            f"No measurement set is present in the target directory: {args.target_datadir}"
        )
        return

    total_ncoarse = 0
    for msname in target_mslist:
        ncoarse = get_ncoarse(msname)
        total_ncoarse += ncoarse
    total_ncoarse = max(1, total_ncoarse)

    target_ms_sizes = [get_ms_size(target_msname) for target_msname in target_mslist]
    total_ms_size_target = sum(target_ms_sizes)
    min_mem_target = round(10 * total_ms_size_target / total_ncoarse, 2)
    min_mem_cal = 0.0

    if args.cal_datadir:
        if os.path.exists(args.cal_datadir):
            cal_mslist = glob.glob(f"{args.cal_datadir}/*.ms")
            if len(cal_mslist) == 0:
                print(
                    f"No calibrator measurement set is present in: {args.cal_datadir}"
                )
            else:
                cal_ms_sizes = [get_ms_size(cal_msname) for cal_msname in cal_mslist]
                total_ms_size_cal = sum(cal_ms_sizes)
                min_mem_cal = round(10 * total_ms_size_cal / total_ncoarse, 2)
        else:
            print(f"Calibrator data direcotry does not exist.")

    min_mem = max(min_mem_target, min_mem_cal)

    ###############################################
    # Setup cluster environment
    ###############################################
    if scheduler_name == "local":
        pid = os.getpid()
    elif scheduler_name == "slurm":
        pid = os.environ.get("SLURM_JOB_ID")
    else:
        print("P-AIRCARS is only ready for local or slurm cluster.")
        return 1

    if args.cluster is True and scheduler_name == "local":
        print(
            "User wants to use cluster architechture, but no job scheduler is available. Stopping P-AIRCARS."
        )
        return

    if args.mem_frac <= 0:
        mem_frac = 0.8
    else:
        mem_frac = args.mem_frac
    if args.cpu_frac <= 0:
        cpu_frac = 0.8
    else:
        cpu_frac = args.cpu_frac

    if args.max_worker is None:
        max_worker = total_ncoarse + 1
    else:
        max_worker = max(int(args.max_worker), total_ncoarse + 1)
    max_worker = max(2, max_worker)  # Minimum 2 workers are needed

    slurm_job = is_slurm_job()
    if args.cluster is not True or scheduler_name == "local" or slurm_job is False:
        #######################################
        # Set up local cluster
        #######################################
        print("Setting up local cluster....")
        print(f"Maximum allowed worker: {max_worker}")
        print(f"Minimum per worker memory: {min_mem}GB")
        result = get_local_dask_cluster(
            args.workdir,
            cpu_frac=cpu_frac,
            mem_frac=mem_frac,
            min_mem=min_mem,
            max_worker=max_worker,
        )
        if result is None:
            print("Error occured in creating local cluster.")
            return 1
        else:
            dask_client, dask_cluster, dask_dir, nworker = result

        scheduler_address = dask_client.scheduler.address
        main_job_file = save_main_process_info(
            pid,
            jobid,
            scheduler_address,
            args.target_datadir,
            os.path.abspath(args.workdir),
            os.path.abspath(args.outdir),
            args.cpu_frac,
            args.mem_frac,
        )
        adaptive = True  # For local cluster, always do adaptive scaling to avoid occupying resources
    else:
        if scheduler_name == "slurm":
            if args.partition is None:
                print("Please provide partition name to submit SLURM jobs.")
                return
            ########################################
            # Setting up slurm cluster
            ########################################
            print("Setting up slurm cluster....")
            print(f"Maximum allowed worker: {max_worker}")
            print(f"Minimum per worker memory: {min_mem}GB")
            cluster_result = get_slurm_dask_cluster(
                args.workdir,
                jobid=jobid,
                cpu_frac=cpu_frac,
                mem_frac=mem_frac,
                max_worker=max_worker,
                min_mem=min_mem,
                partition=args.partition,
                account=args.account,
                walltime=args.walltime,
            )
            if cluster_result is None:
                print("Error occured in creating slurm cluster.")
                return 1
            else:
                dask_client, dask_cluster, dask_dir, nworker = cluster_result
            scheduler_address = dask_client.scheduler.address
            main_job_file = save_main_process_info(
                pid,
                jobid,
                scheduler_address,
                args.target_datadir,
                os.path.abspath(args.workdir),
                os.path.abspath(args.outdir),
                args.cpu_frac,
                args.mem_frac,
            )
            adaptive = args.adaptive
            if not adaptive:
                nworker = min(total_ncoarse + 1, nworker)
                scale_worker_and_wait(dask_cluster, dask_client, nworker)
        else:
            print(
                f"P-AIRCARS is under development for job scheduler: {scheduler_name}. Stopping P-AIRCARS."
            )
            return

    ##########################################
    # Starting pipeline
    ##########################################
    try:
        dask_addr = dask_client.scheduler.address
        print("#########################################")
        print("Starting P-AIRCARS Pipeline....")
        print("#########################################")
        print(f"Total dask workers: {nworker}")
        msg = master_control.with_options(
            flow_run_name=f"paircars_{jobid}",
            task_runner=DaskTaskRunner(address=dask_addr),
        )(
            args.target_datadir,
            args.target_metafits,
            args.workdir,
            args.outdir,
            calibrator_datadir=args.cal_datadir,
            calibrator_metafits=args.cal_metafits,
            solar_data=args.solar_data,
            # Pre-calibration
            do_forcereset_weightflag=args.do_forcereset_weightflag,
            do_cal_flag=args.do_cal_flag,
            do_import_model=args.do_import_model,
            # Basic calibration
            do_basic_cal=args.do_basic_cal,
            do_applycal=args.do_applycal,
            # Target data preparation
            freqrange=args.freqrange,
            timerange=args.timerange,
            uvrange=args.cal_uvrange,
            # Polarization calibration
            do_polcal=args.do_polcal,
            # Self-calibration
            do_selfcal=args.do_selfcal,
            do_apply_selfcal=args.do_apply_selfcal,
            only_amplitude=args.only_amplitude,
            do_ap_selfcal=args.do_ap_selfcal,
            solar_selfcal=args.solar_selfcal,
            solint=args.solint,
            # Sidereal correction
            do_sidereal_cor=args.do_sidereal_cor,
            do_move_solarcenter=args.do_move_solarcenter,
            # Dynamic spectra
            make_ds=args.make_ds,
            # Imaging
            do_imaging=args.do_imaging,
            do_pbcor=args.do_pbcor,
            weight=args.weight,
            robust=args.robust,
            minuv=args.minuv,
            image_freqres=args.image_freqres,
            image_timeres=args.image_timeres,
            pol=args.pol,
            clean_threshold=args.clean_threshold,
            use_multiscale=args.use_multiscale,
            use_solar_mask=args.use_solar_mask,
            cutout_rsun=args.cutout_rsun,
            make_overlay=args.make_overlay,
            make_msplot=args.make_msplot,
            # Resource settings
            cpu_frac=args.cpu_frac,
            mem_frac=args.mem_frac,
            max_worker=nworker,
            keep_backup=args.keep_backup,
            keep_calibrated_ms=args.keep_calibrated_ms,
            # Remote logging
            masterlog=args.masterlog,
            remote_logger=args.remote_logger,
            jobid=jobid,
            job_password=args.job_password,
            adaptive=adaptive,
        )
        print("##########################################")
        if msg == 0:
            print("P-AIRCARS successfully executed.")
        else:
            print("Issued occured in P-AIRCARS execution.")
        print("##########################################")
    except Exception:
        traceback.print_exc()
    finally:
        time.sleep(5)
        print("Clearning caches...")
        drop_cache(args.target_datadir)
        drop_cache(args.cal_datadir)
        drop_cache(args.workdir)
        drop_cache(args.outdir)
        dask_client.shutdown()
        dask_client.close()
        dask_cluster.close()
        os.system(f"rm -rf {dask_dir}")
