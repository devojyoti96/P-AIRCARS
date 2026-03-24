import subprocess
import pickle
import os
import tempfile
from .udocker_utils import (
    initialize_hyperbeam_container,
    check_udocker_container,
    init_udocker,
)


class FEEBeam:
    def __init__(self, pbfile, check_container=False):
        """
        MWA FEE beam class

        Parameters
        ----------
        pbfile : str
            MWA PB HDF5 file
        check_container : bool, optional
            Check hyperbeam container
        """
        self.pbfile = pbfile
        self.check_container = check_container
        init_udocker()
        self.env = os.environ.copy()
        self.container_name = "paircarshyperbeam"
        if self.check_container:
            container_present = check_udocker_container(self.container_name)
            if not container_present:
                print(f"Initializing {self.container_name}...")
                self.container_name = initialize_hyperbeam_container(
                    name=self.container_name, verbose=True
                )
                if self.container_name is None:
                    print(
                        f"Container {self.container_name} is not initiated. First initiate container and then run."
                    )
                    return

    def calc_jones_array(
        self,
        az_rad,
        za_rad,
        freq,
        delay,
        amps,
        norm,
        lat,
        iau_order,
    ):
        """
        Calc primary beam jones array

        Parameters
        ----------
        az_rad : numpy.array
            Azimuth array in radians
        za_rad : numpy.array
            Zenith angle array in radians
        freq : float
            Frequency in Hz
        delay : list
            List of beamformer delays
        amps : list
            List of dipole amplitudes
        norm : bool
            Zenith normalization or not
        lat : float
            Latitude of the array
        iau_order : bool
            MWA beam in IAU order or not

        Returns
        -------
        numpy.array
            Beam Jones (shape: coodinates, 4 components)
        """
        env_vars = [
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
            "RAYON_NUM_THREADS",
        ]
        pbdir = os.path.dirname(os.path.abspath(self.pbfile))
        temp_name = "mwapb_udocker_" + next(tempfile._get_candidate_names())
        temp_pbdir_path = os.path.join(pbdir, temp_name)
        data = dict(
            az_rad=az_rad,
            za_rad=za_rad,
            freq=freq,
            delay=delay,
            amps=amps,
            norm=norm,
            lat=lat,
            iau_order=iau_order,
            pbfile=f"{temp_pbdir_path}/{os.path.basename(self.pbfile)}",
        )
        full_command = ["udocker", "--quiet", "run", "--nobanner"]
        env_keys = list(self.env.keys())
        for var in env_vars:
            if var in env_keys:
                full_command.append(f"--env={var}={self.env[var]}")
        full_command = full_command + [
            f"--volume={pbdir}:{temp_pbdir_path}",
            "--workdir",
            f"{temp_pbdir_path}",
            f"{self.container_name}",
            "python",
            "/app/hyperbeam_array.py",
        ]
        proc = subprocess.run(
            full_command,
            env=self.env,
            input=pickle.dumps(data),
            capture_output=True,
        )
        return pickle.loads(proc.stdout)

    def calc_jones(
        self,
        az_rad,
        za_rad,
        freq,
        delay,
        amps,
        norm,
        lat,
        iau_order,
    ):
        """
        Calc primary beam jones array

        Parameters
        ----------
        az_rad : float
            Azimuth in radian
        za_rad : float
            Zenith angle in radian
        freq : float
            Frequency in Hz
        delay : list
            List of beamformer delays
        amps : list
            List of dipole amplitudes
        norm : bool
            Zenith normalization or not
        lat : float
            Latitude of the array
        iau_order : bool
            MWA beam in IAU order or not

        Returns
        -------
        numpy.array
            Beam Jones (shape: coodinates, 4 components)
        """
        pbdir = os.path.dirname(os.path.abspath(self.pbfile))
        temp_name = "mwapb_udocker_" + next(tempfile._get_candidate_names())
        temp_pbdir_path = os.path.join(pbdir, temp_name)
        data = dict(
            az_rad=az_rad,
            za_rad=za_rad,
            freq=freq,
            delay=delay,
            amps=amps,
            norm=norm,
            lat=lat,
            iau_order=iau_order,
            pbfile=f"{temp_pbdir_path}/{os.path.basename(self.pbfile)}",
        )
        proc = subprocess.run(
            [
                "udocker",
                "run",
                "--nobanner",
                f"--volume={pbdir}:{temp_pbdir_path}",
                "--workdir",
                f"{temp_pbdir_path}",
                f"{self.container_name}",
                "python",
                "/app/hyperbeam_single.py",
            ],
            env=self.env,
            input=pickle.dumps(data),
            capture_output=True,
        )
        return pickle.loads(proc.stdout)
