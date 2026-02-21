import psutil
import traceback
import tempfile
import time
import glob
import os
import subprocess
from .basic_utils import get_datadir

####################
# uDOCKER related
####################


def set_udocker_env():
    datadir = get_datadir()
    if (
        datadir is None
        or os.path.exists(datadir) is False
        or os.path.exists(f"{datadir}/udocker-englib-1.2.11.tar.gz") is False
    ):
        print("P-AIRCARS data directory and docker environment is not setup yet")
        return
    udocker_dir = f"{datadir}/udocker"
    os.makedirs(udocker_dir, exist_ok=True)
    os.environ["UDOCKER_DIR"] = udocker_dir
    os.environ["UDOCKER_TARBALL"] = f"{datadir}/udocker-englib-1.2.11.tar.gz"
    return datadir


def init_udocker():
    set_udocker_env()
    os.system("udocker install")


def check_udocker_container(name):
    """
    Check whether a docker container is present or not

    Parameters
    ----------
    name : str
        Container name

    Returns
    -------
    bool
        Whether present or not
    """
    set_udocker_env()
    try:
        result = subprocess.run(
            ["udocker", "--insecure", "--quiet", "inspect", name],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return result.returncode == 0
    except Exception:
        return False

def initialize_container(image_name, name, update=False, verbose=False):
    """
    Initialize container
    
    Parameters
    ----------
    image_name: str
        Docker image name
    name : str
        Container name
    update : bool, optional
        Update or not
    verbose : bool, optional
        Verbose output
        
    Returns
    -------
    bool
        Whether initialized successfully or not
    """
    set_udocker_env()
    check_cmd = f"udocker images | grep -q {image_name}"
    image_exists = os.system(check_cmd)
    if image_exists != 0:
        if verbose:
            result = subprocess.run(
                ["udocker", "pull", f"{image_name}"],
            )
        else:
            result = subprocess.run(
                ["udocker", "pull", f"{image_name}"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        a = result.returncode
    else:
        if update:
            if verbose:
                subprocess.run(
                    ["udocker", "rm", f"{name}"],
                )
                subprocess.run(
                    ["udocker", "rmi", f"{image_name}"],
                )
            else:
                subprocess.run(
                    ["udocker", "rm", f"{name}"],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
                subprocess.run(
                    ["udocker", "rmi", f"{image_name}"],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
            print("Re-downloading docker image.")
            if verbose:
                result = subprocess.run(
                    ["udocker", "pull", f"{image_name}"],
                )
            else:
                result = subprocess.run(
                    ["udocker", "pull", f"{image_name}"],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
            a = result.returncode
            if a == 0:
                print("Re-downloaded docker image.")
            else:
                print("Re-downloading container image is failed.")
                return
        else:
            print(f"Image {image_name} already present.")
            a = 0
    if a == 0:
        if verbose:
            result = subprocess.run(
                ["udocker", "pull", f"{image_name}"],
            )
        else:
            result = subprocess.run(
                ["udocker", "create", f"--name={name}", f"{image_name}"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        a = result.returncode
        print(f"Container started with name : {name}")
        return name
    else:
        print(f"Container could not be created with name : {name}")
        return

def initialize_wsclean_container(name="paircarswsclean", update=False):
    """
    Initialize WSClean container

    Parameters
    ----------
    name : str, optional
        Name of the container
    update : bool, optional
        Update container

    Returns
    -------
    bool
        Whether initialized successfully or not
    """
    print("Initializing wsclean container.")
    image_name = "devojyoti96/wsclean-solar:latest"
    msg = initialize_container(image_name, name, update=update, verbose=verbose)
    return msg


def initialize_quartical_container(name="paircarsquartical", update=False, verbose=True):
    """
    Initialize quartical container

    Parameters
    ----------
    name : str, optional
        Name of the container
    update : bool, optional
        Update container
    verbose : bool, optional
        Verbose output 

    Returns
    -------
    bool
        Whether initialized successfully or not
    """
    print("Initializing quartical container.")
    image_name = "devojyoti96/quartical:0.2.6"
    msg = initialize_container(image_name, name, update=update, verbose=verbose)
    return msg


def initialize_shadems_container(name="paircarsshadems", update=False, verbose=False):
    """
    Initialize shadems container

    Parameters
    ----------
    name : str, optional
        Name of the container
    update : bool, optional
        Update container
    verbose : bool, optional
        Verbose output

    Returns
    -------
    bool
        Whether initialized successfully or not
    """
    print("Initializing shadems container.")
    image_name = "devojyoti96/shadems:v0.5.4"
    msg = initialize_container(image_name, name, update=update, verbose=verbose)
    return msg


def run_wsclean(
    wsclean_cmd,
    container_name="paircarswsclean",
    check_container=False,
    verbose=False,
):
    """
    Run WSClean inside a udocker container (no root permission required).

    Parameters
    ----------
    wsclean_cmd : str
        Full WSClean command as a string.
    container_name : str, optional
        Container name
    check_container : bool, optional
        Check container presence or not
    verbose : bool, optional
        Verbose output or not

    Returns
    -------
    int
        Success message
    """
    set_udocker_env()
    def show_file(path):
        try:
            print(open(path).read())
        except Exception as e:
            print(f"{e}")

    if check_container:
        container_present = check_udocker_container(container_name)
        if not container_present:
            container_name = initialize_wsclean_container(name=container_name)
            if container_name is None:
                print(
                    f"Container {container_name} is not initiated. First initiate container and then run."
                )
                return 1
    msname = wsclean_cmd.split(" ")[-1]
    msname = os.path.abspath(msname)
    mspath = os.path.dirname(msname)
    temp_name = "wsclean_udocker_" + next(tempfile._get_candidate_names())
    temp_docker_path = os.path.join(mspath, temp_name)
    wsclean_cmd_args = wsclean_cmd.split(" ")[:-1]
    if "-fits-mask" in wsclean_cmd_args:
        index = wsclean_cmd_args.index("-fits-mask")
        name = wsclean_cmd_args[index + 1]
        namedir = os.path.dirname(os.path.abspath(name))
        basename = os.path.basename(os.path.abspath(name))
        wsclean_cmd_args.remove(name)
        wsclean_cmd_args.insert(index + 1, temp_docker_path + "/" + basename)
    if "-name" not in wsclean_cmd_args:
        wsclean_cmd_args.append(
            "-name " + temp_docker_path + "/" + os.path.basename(msname).split(".ms")[0]
        )
    else:
        index = wsclean_cmd_args.index("-name")
        name = wsclean_cmd_args[index + 1]
        namedir = os.path.dirname(os.path.abspath(name))
        basename = os.path.basename(os.path.abspath(name))
        wsclean_cmd_args.remove(name)
        wsclean_cmd_args.insert(index + 1, temp_docker_path + "/" + basename)
    if "-temp-dir" not in wsclean_cmd_args:
        wsclean_cmd_args.append("-temp-dir " + temp_docker_path)
    else:
        index = wsclean_cmd_args.index("-temp-dir")
        name = os.path.abspath(wsclean_cmd_args[index + 1])
        wsclean_cmd_args.remove(name)
        wsclean_cmd_args.insert(index + 1, temp_docker_path)
    wsclean_cmd = (
        " ".join(wsclean_cmd_args)
        + " "
        + temp_docker_path
        + "/"
        + os.path.basename(msname)
    )
    wsclean_cmd_args = wsclean_cmd.split(" ")
    try:
        full_command = [
            "udocker",
            "run",
            "--nobanner",
            f"--volume={mspath}:{temp_docker_path}",
            "--workdir",
            f"{temp_docker_path}",
            f"{container_name}",
        ] + wsclean_cmd_args
        if verbose:
            print(f"{wsclean_cmd}\n")
            result = subprocess.run(
                full_command,
            )
        else:
            result = subprocess.run(
                full_command,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        exit_code = result.returncode
        if exit_code != 0:
            print("##########################")
            print(os.path.basename(msname))
            print("##########################")
        return 0 if exit_code == 0 else 1
    except Exception as e:
        traceback.print_exc()
        return 1


def run_solar_sidereal_cor(
    msname="",
    only_uvw=False,
    container_name="paircarswsclean",
    check_container=False,
    verbose=False,
):
    """
    Run chgcenter inside a udocker container to correct solar sidereal motion (no root permission required).

    Parameters
    ----------
    msname : str
        Name of the measurement set
    only_uvw : bool, optional
        Update only UVW values
        Note: This is required when visibilities are properly phase rotated in correlator to track the Sun,
        but while creating the MS, UVW values are estimated using the first phasecenter of the Sun.
    check_container : bool, optional
        Check container
    container_name : str, optional
        Container name
    verbose : bool, optional
        Verbose output or not

    Returns
    -------
    int
        Success message
    """
    set_udocker_env()
    if check_container:
        container_present = check_udocker_container(container_name)
        if not container_present:
            container_name = initialize_wsclean_container(name=container_name)
            if container_name is None:
                print(
                    f"Container {container_name} is not initiated. First initiate container and then run."
                )
                return 1
    msname = os.path.abspath(msname)
    mspath = os.path.dirname(msname)
    temp_name = "chgcenter_udocker_" + next(tempfile._get_candidate_names())
    temp_docker_path = os.path.join(mspath, temp_name)
    if only_uvw:
        cmd = (
            "chgcentre -only-uvw -solarcenter "
            + temp_docker_path
            + "/"
            + os.path.basename(msname)
        )
    else:
        cmd = (
            "chgcentre -solarcenter "
            + temp_docker_path
            + "/"
            + os.path.basename(msname)
        )
    cmd_args = cmd.split(" ")
    try:
        full_command = [
            "udocker",
            "--quiet",
            "run",
            "--nobanner",
            f"--volume={mspath}:{temp_docker_path}",
            "--workdir",
            f"{temp_docker_path}",
            "paircarswsclean",
        ] + cmd_args
        if verbose:
            print(f"{cmd}\n")
            result = subprocess.run(
                full_command,
            )
        else:
            result = subprocess.run(
                full_command,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        exit_code = result.returncode
        return 0 if exit_code == 0 else 1
    except Exception as e:
        traceback.print_exc()
        return 1


def run_chgcenter(
    msname,
    ra,
    dec,
    only_uvw=False,
    container_name="paircarswsclean",
    check_container=False,
    verbose=False,
):
    """
    Run chgcenter inside a udocker container (no root permission required).

    Parameters
    ----------
    msname : str
        Name of the measurement set
    ra : str
        RA can either be 00h00m00.0s or 00:00:00.0
    dec : str
        Dec can either be 00d00m00.0s or 00.00.00.0
    only_uvw : bool, optional
        Update only UVW values
        Note: This is required when visibilities are properly phase rotated in correlator,
        but while creating the MS, UVW values are estimated using a wrong phase center.
    check_container : bool, optional
        Check container
    container_name : str, optional
        Container name
    verbose : bool, optional
        Verbose output

    Returns
    -------
    int
        Success message
    """
    set_udocker_env()
    if check_container:
        container_present = check_udocker_container(container_name)
        if not container_present:
            container_name = initialize_wsclean_container(name=container_name)
            if container_name is None:
                print(
                    f"Container {container_name} is not initiated. First initiate container and then run."
                )
                return 1
    msname = os.path.abspath(msname)
    mspath = os.path.dirname(msname)
    temp_name = "chgcenter_udocker_" + next(tempfile._get_candidate_names())
    temp_docker_path = os.path.join(mspath, temp_name)
    if only_uvw:
        cmd = (
            "chgcentre -only-uvw "
            + temp_docker_path
            + "/"
            + os.path.basename(msname)
            + " "
            + ra
            + " "
            + dec
        )
    else:
        cmd = (
            "chgcentre "
            + temp_docker_path
            + "/"
            + os.path.basename(msname)
            + " "
            + ra
            + " "
            + dec
        )
    cmd_args = cmd.split(" ")
    try:
        full_command = [
            "udocker",
            "--quiet",
            "run",
            "--nobanner",
            f"--volume={mspath}:{temp_docker_path}",
            "--workdir",
            f"{temp_docker_path}",
            f"{container_name}",
        ] + cmd_args
        if verbose:
            print(f"{cmd}\n")
            result = subprocess.run(
                full_command,
            )
        else:
            result = subprocess.run(
                full_command,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        exit_code = result.returncode
        return 0 if exit_code == 0 else 1
    except Exception as e:
        traceback.print_exc()
        return 1


def run_shadems(
    cmd,
    container_name="paircarsshadems",
    check_container=False,
    verbose=False,
):
    """
    Run shadems inside a udocker container (no root permission required).

    Parameters
    ----------
    cmd : str
        Shadems command
    container_name : str, optional
        Container name
    check_container : bool, optional
        Check container
    verbose : bool, optional
        Verbose output

    Returns
    -------
    int
        Success message
    """
    set_udocker_env()
    if check_container:
        container_present = check_udocker_container(container_name)
        if not container_present:
            container_name = initialize_shadems_container(name=container_name)
            if container_name is None:
                print(
                    f"Container {container_name} is not initiated. First initiate container and then run."
                )
                return 1
    splited_cmd = cmd.split(" ")
    if splited_cmd[-1] in ["-h", "--help"]:
        verbose = True
        datapath = os.getcwd()
    else:
        msname = splited_cmd[-1]
        datapath = os.path.dirname(os.path.abspath(msname))
    temp_name = "shadems_udocker_" + next(tempfile._get_candidate_names())
    temp_docker_path = os.path.join(datapath, temp_name)
    if splited_cmd[-1] not in ["-h", "--help"]:
        cmd = f"{' '.join(splited_cmd[:-1])} {temp_docker_path}/{os.path.basename(msname)}"
    cmd_args = cmd.split(" ")
    try:
        full_command = [
            "udocker",
            "--quiet",
            "run",
            "--nobanner",
            f"--volume={datapath}:{temp_docker_path}",
            "--workdir",
            f"{temp_docker_path}",
            f"{container_name}",
        ] + cmd_args
        if verbose:
            print(f"{cmd}\n")
            result = subprocess.run(
                full_command,
            )
        else:
            result = subprocess.run(
                full_command,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        exit_code = result.returncode
        return 0 if exit_code == 0 else 1
    except Exception as e:
        traceback.print_exc()
        return 1


def run_quartical(
    cmd,
    container_name="paircarsquartical",
    check_container=False,
    verbose=False,
):
    """
    Run quartical inside a udocker container (no root permission required).

    Parameters
    ----------
    cmd : str
        Quartical command
    container_name : str, optional
        Container name
    check_container : bool, optional
        Check container
    verbose : bool, optional
        Verbose output

    Returns
    -------
    int
        Success message
    """
    set_udocker_env()
    if check_container:
        container_present = check_udocker_container(container_name)
        if not container_present:
            container_name = initialize_quartical_container(name=container_name)
            if container_name is None:
                print(
                    f"Container {container_name} is not initiated. First initiate container and then run."
                )
                return 1
    splited_cmd = cmd.split(" ")
    if len(splited_cmd) == 1 and "goquartical" in cmd:
        verbose = True
        datapath = os.getcwd()
        temp_name = "quartical_udocker_" + next(tempfile._get_candidate_names())
        temp_docker_path = os.path.join(datapath, temp_name)
    elif len(splited_cmd) > 1:
        for i in range(len(splited_cmd)):
            cmd_arg = splited_cmd[i]
            if "input_ms.path" in cmd_arg:
                msname = cmd_arg.split("input_ms.path=")[-1]
                datapath = os.path.dirname(os.path.abspath(msname))
                temp_name = "quartical_udocker_" + next(tempfile._get_candidate_names())
                temp_docker_path = os.path.join(datapath, temp_name)
                temp_msname = f"{temp_docker_path}/{os.path.basename(msname)}"
                cmd_arg = f"input_ms.path={temp_msname}"
                splited_cmd[i] = cmd_arg
            if "output.gain_directory" in cmd_arg:
                caltable = cmd_arg.split("output.gain_directory=")[-1]
                temp_caltable = f"{temp_docker_path}/{os.path.basename(caltable)}"
                cmd_arg = f"output.gain_directory={temp_caltable}"
                splited_cmd[i] = cmd_arg
            if "output.log_directory" in cmd_arg:
                log = cmd_arg.split("output.log_directory=")[-1]
                temp_log = f"{temp_docker_path}/{os.path.basename(log)}"
                cmd_arg = f"output.log_directory={temp_log}"
                splited_cmd[i] = cmd_arg
            if "load_from" in cmd_arg:
                gaintable = cmd_arg.split("load_from=")[-1]
                gain_path = os.path.dirname(os.path.dirname(gaintable))
                if gain_path != datapath:
                    temp_gain_path = (
                        f"{datapath}/{os.path.basename(os.path.dirname(gaintable))}"
                    )
                    os.system(f"rm -rf {temp_gain_path}")
                    os.system(f"cp -r {os.path.dirname(gaintable)} {temp_gain_path}")
                gaintable = gaintable.split("/")[-2:]
                gaintable = "/".join(gaintable)
                temp_gaintable = f"{temp_docker_path}/{gaintable}"
                cmd_arg = f"{cmd_arg.split('=')[0]}={temp_gaintable}"
                splited_cmd[i] = cmd_arg
        cmd = " ".join(splited_cmd)
    else:
        print("Please provide valid command.")
        return 1
    cmd_args = cmd.split(" ")
    try:
        full_command = [
            "udocker",
            "--quiet",
            "run",
            "--nobanner",
            f"--volume={datapath}:{temp_docker_path}",
            "--workdir",
            f"{temp_docker_path}",
            f"{container_name}",
        ] + cmd_args
        if verbose:
            print(f"{cmd}\n")
            result = subprocess.run(
                full_command,
            )
        else:
            result = subprocess.run(
                full_command,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        exit_code = result.returncode
        if "load_from" in cmd_arg and gain_path != datapath:
            os.system(f"rm -rf {temp_gain_path}")
        return 0 if exit_code == 0 else 1
    except Exception as e:
        traceback.print_exc()
        return 1

