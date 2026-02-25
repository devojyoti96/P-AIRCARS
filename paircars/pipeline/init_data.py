import logging
import psutil
import argparse
import requests
import sys
import os
import getpass
from datetime import datetime as dt
from parfive import Downloader
from paircars.utils.basic_utils import (
    create_datadir,
    get_datadir,
    get_cachedir,
    check_port_status,
    get_free_port,
)
from paircars.utils.logger_utils import SmartDefaultsHelpFormatter, clean_shutdown
from paircars.utils.prefect_setup_utils import start_server
from paircars.utils.resource_utils import has_space
from paircars.utils.proc_manage_utils import get_scheduler_name
from paircars.utils.udocker_utils import (
    init_udocker,
    initialize_wsclean_container,
    initialize_quartical_container,
    initialize_shadems_container,
    initialize_hyperdrive_container,
)
from paircars.pipeline.beam_interpolate import *

logging.getLogger("distributed").setLevel(logging.ERROR)
logging.getLogger("tornado.application").setLevel(logging.CRITICAL)

all_filenames = [
    "udocker-englib-1.2.11.tar.gz",
    "de440s.bsp",
    "GGSM.txt",
    "haslam_map.fits",
    "MWA_sweet_spots.npy",
    "Ref_mean_bandpass_final.npy",
    "postgres_credentials.npy",
    "mwa_full_embedded_element_pattern.h5",
]


def get_zenodo_file_urls(record_id):
    url = f"https://zenodo.org/api/records/{record_id}"
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()
    return [(f["links"]["self"], f["key"]) for f in data.get("files", [])]


def download_with_parfive(record_id, update=False, output_dir="zenodo_download"):
    print("####################################")
    print("Downloading P-AIRCARS data files ...")
    print("####################################")
    urls = get_zenodo_file_urls(record_id)
    urls.append(
        (
            "http://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5",
            "mwa_full_embedded_element_pattern.h5",
        )
    )
    os.makedirs(output_dir, exist_ok=True)
    total_cpu = psutil.cpu_count()
    dl = Downloader(max_conn=min(total_cpu, len(all_filenames) + 1))
    for file_url, filename in urls:
        if filename in all_filenames:
            if os.path.exists(f"{output_dir}/{filename}") == False or update:
                if os.path.exists(f"{output_dir}/{filename}"):
                    os.system(f"rm -rf {output_dir}/{filename}")
                dl.enqueue_file(file_url, path=output_dir, filename=filename)
    results = dl.download()
    for f in results:
        os.chmod(f, 0o755)


def init_paircars_data(update=False, remote_link=None, emails=None):
    """
    Initiate P-AIRCARS data

    Parameters
    ----------
    update : bool, optional
        Update data, if already exists
    remote_link : str, optional
        Remote logger link to save in database
    emails : str, optional
        Email addresses to send remote logger JobID and password
    """
    datadir = get_datadir()
    os.makedirs(datadir, exist_ok=True)
    cachedir = get_cachedir()
    username = getpass.getuser()
    linkfile = f"{cachedir}/remotelink_{username}.txt"
    emailfile = f"{cachedir}/emails_{username}.txt"
    if not os.path.exists(linkfile):
        with open(linkfile, "w") as f:
            f.write("")

    if remote_link is not None:
        with open(linkfile, "w") as f:
            f.write(str(remote_link))

    if emails is not None:
        with open(emailfile, "w") as f:
            f.write(str(emails))

    unavailable_files = [
        f for f in all_filenames if not os.path.exists(f"{datadir}/{f}")
    ]

    if unavailable_files or update:
        record_id = "18640418"
        download_with_parfive(record_id, update=update, output_dir=datadir)
        timestr = dt.utcnow().strftime("%Y-%m-%d %H:%M:%S")
        print(f"P-AIRCARS data are updated in: {datadir} at time: {timestr}")
    freqres_list = [40, 80, 160, 320, 640]
    mwapb_file = f"{datadir}/mwa_full_embedded_element_pattern.h5"
    for freqres in freqres_list:
        outfile = mwapb_file.split(".h5")[0] + f"_{freqres}.h5"
        if os.path.exists(outfile) is False or update:
            print(f"Making interpolated beam at frequency resolution: {freqres} kHz")
            do_beam_interpolate(mwapb_file, new_freq_res=int(freqres))


def main(
    init=False,
    datadir="",
    update=False,
    link=None,
    emails=None,
):
    """
    Initiate P-AIRCARS setup

    Parameters
    ----------
    init : bool, optional
        Initiate setup
    datadir : str, optional
        User provided custom data directory
    update : bool, optional
        Update existing data (if corrupted by somehow)
    link : str, optional
        Remote link
    emails : str, optional
        E-mails for notifications
    """
    required_gb = 20

    port = 4260
    postgres_port = 5260

    if check_port_status(port) is False:
        if scheduler_name != "local":
            port = get_free_port(start_port=4260, end_port=5250)

    if check_port_status(postgres_port) is False:
        if scheduler_name != "local":
            postgres_portport = get_free_port(start_port=5260, end_port=6250)

    if init:
        ######################################
        # Downloading data files
        ######################################
        create_datadir(datadir=datadir)
        datadir = get_datadir()
        print(f"P-AIRCARS data directory: {datadir}")
        if has_space(datadir, required_gb) is False:
            print(
                f"Minimum {required_gb}GB disk space is required in data directory: {datadir}. Please check disk space."
            )
            return 1
        init_paircars_data(update=update, remote_link=link, emails=emails)
        print(f"P-AIRCARS data are initiated.")

        #########################################
        # Docker containers initiation
        #########################################
        init_udocker()
        print("uDOCKER is inititalized")
        wsclean_container_name = initialize_wsclean_container(
            update=update, verbose=True
        )
        if (
            wsclean_container_name is not None
            and wsclean_container_name == "paircarswsclean"
        ):
            print("WSClean container is initialized")
        else:
            return 1
        quartical_container_name = initialize_quartical_container(
            update=update, verbose=True
        )
        if (
            quartical_container_name is not None
            and quartical_container_name == "paircarsquartical"
        ):
            print("Quartical container is initialized")
        else:
            return 1
        shadems_container_name = initialize_shadems_container(
            update=update, verbose=True
        )
        if (
            shadems_container_name is not None
            and shadems_container_name == "paircarsshadems"
        ):
            print("Shadems container is initialized")
        else:
            return 1
        hyperdrive_container_name = initialize_hyperdrive_container(
            update=update, verbose=True
        )
        if (
            hyperdrive_container_name is not None
            and hyperdrive_container_name == "paircarshyperdrive"
        ):
            print("Hyperdrive container is initialized")
        else:
            return 1
        postgres_container_name = initialize_postgres_container(
            update=update, verbose=True
        )
        if (
            postgres_container_name is not None
            and postgres_container_name == "paircarspostgres"
        ):
            print("PostgreSQL container is initialized")
        else:
            return 1

        #########################################
        # prefect server setup
        #########################################
        print ("Prefect setup....")
        scheduler_name = get_scheduler_name()
        msg, config_file, profile_path, env_file, dashboard, pid_file = start_server(
            port, postgres_port, scheduler_name=scheduler_name
        )
        config = np.load(config_file, allow_pickle=True).all()
        if msg != 0:
            if scheduler_name != "local":
                print(
                    f"P-AIRCARS will not work in prefect ephemeral mode in cluster environment with job scheduler: {scheduler_name}"
                )
                return 1
            else:
                print(
                    f"Error in starting prefect server at port. P-AIRCARS will use ephemeral mode in local cluster."
                )
        return 0
    else:
        return 1


def cli():
    usage = "Initiate P-AIRCARS data"
    parser = argparse.ArgumentParser(
        description=usage, formatter_class=SmartDefaultsHelpFormatter
    )
    parser.add_argument("--init", action="store_true", help="Initiate data")
    parser.add_argument(
        "--datadir", type=str, default="", help="User provided data directory"
    )
    parser.add_argument("--update", action="store_true", help="Update existing data")
    parser.add_argument(
        "--remotelink", dest="link", default=None, help="Set remote log link"
    )
    parser.add_argument(
        "--emails",
        dest="emails",
        default=None,
        help="Email addresses (comma seperated) to send Job ID and password for remote logger",
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        return 1

    args = parser.parse_args()

    msg = main(
        init=args.init,
        datadir=args.datadir,
        update=args.update,
        link=args.link,
        emails=args.emails,
    )
    if msg != 0:
        print("Error in initial setup.")
    return msg


if __name__ == "__main__":
    msg = cli()
    os._exit(msg)
