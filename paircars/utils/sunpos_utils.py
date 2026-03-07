import astropy.units as u
import glob
import os
import numpy as np
from astropy.time import Time
from astropy.coordinates import (
    EarthLocation,
    AltAz,
    get_sun,
    solar_system_ephemeris,
)
from astropy.io import fits
from astropy.wcs import WCS
from casatools import msmetadata
from .basic_utils import get_datadir, mjdsec_to_timestamp
from .udocker_utils import run_solar_sidereal_cor, run_chgcenter
from .image_utils import create_circular_mask_array

#####################################
# Sun position related
#####################################
datadir = get_datadir()
try:
    solar_system_ephemeris.set(f"{datadir}/de440s")
except:
    solar_system_ephemeris.set("builtin")


def get_solar_elevation(lat, lon, elev, date_time):
    """
    Get solar elevation

    Parameters
    ----------
    lat : float
        Latitude in degrees
    lon : float
        Longitude in degrees
    elev : float
        Elevation in degrees
    date_time : str
        Date time in YYYY-MM-DDThh:mm:ss (ISOT) format, default : present time


    Returns
    -------
    float
        Solar elevation in degree
    """
    latitude = lat * u.deg  # In degree
    longitude = lon * u.deg  # In degree
    elevation = elev * u.m  # In meter
    if date_time == "":
        astro_time = Time.now()
    else:
        astro_time = Time(date_time)
    location = EarthLocation(lat=latitude, lon=longitude, height=elevation)
    sun_coords = get_sun(astro_time)  # In GCRS (geocentric frame)
    altaz_frame = AltAz(obstime=astro_time, location=location)
    sun_altaz = sun_coords.transform_to(altaz_frame)
    solar_elevation = sun_altaz.alt.deg
    return round(solar_elevation, 3)


def radec_sun(msname):
    """
    RA DEC of the Sun at the midpoint of scan (offline version)

    Parameters
    ----------
    msname : str
        Name of the measurement set

    Returns
    -------
    str
        RA DEC of the Sun in J2000
    str
        RA string
    str
        DEC string
    float
        RA in degree
    float
        DEC in degree
    """
    msmd = msmetadata()
    msmd.open(msname)
    times = msmd.timesforspws(0)
    msmd.close()
    msmd.done()
    mid_time = times[int(len(times) / 2)]
    mid_timestamp = mjdsec_to_timestamp(mid_time)
    astro_time = Time(mid_timestamp, scale="utc")
    sun_coord = get_sun(astro_time)  # In GCRS (geocentric frame)
    sun_ra = (
        f"{int(sun_coord.ra.hms.h)}h"
        f"{int(sun_coord.ra.hms.m)}m"
        f"{round(sun_coord.ra.hms.s, 2)}s"
    )
    sun_dec = (
        f"{int(sun_coord.dec.dms.d)}d"
        f"{abs(int(sun_coord.dec.dms.m))}m"
        f"{abs(round(sun_coord.dec.dms.s, 2))}s"
    )
    sun_radec_string = f"J2000 {sun_ra} {sun_dec}"
    return (
        sun_radec_string,
        sun_ra,
        sun_dec,
        sun_coord.ra.deg,
        sun_coord.dec.deg,
    )


def radec_sun_at_time(timestamp):
    """
    RA DEC of the Sun a given time

    Parameters
    ----------
    timestamp : str
        Time in format dd-mm-yyyyThh:mm:ss

    Returns
    -------
    str
        RA DEC of the Sun in J2000
    str
        RA string
    str
        DEC string
    float
        RA in degree
    float
        DEC in degree
    """
    astro_time = Time(timestamp, scale="utc")
    sun_coord = get_sun(astro_time)  # In GCRS (geocentric frame)
    sun_ra = (
        f"{int(sun_coord.ra.hms.h)}h"
        f"{int(sun_coord.ra.hms.m)}m"
        f"{round(sun_coord.ra.hms.s, 2)}s"
    )
    sun_dec = (
        f"{int(sun_coord.dec.dms.d)}d"
        f"{abs(int(sun_coord.dec.dms.m))}m"
        f"{abs(round(sun_coord.dec.dms.s, 2))}s"
    )
    sun_radec_string = f"J2000 {sun_ra} {sun_dec}"
    return (
        sun_radec_string,
        sun_ra,
        sun_dec,
        sun_coord.ra.deg,
        sun_coord.dec.deg,
    )


def move_to_sun(msname, only_uvw=False):
    """
    Move the phasecenter of the measurement set at the center of the Sun (Assuming ms has one scan)

    Parameters
    ----------
    msname : str
        Name of the measurement set
    only_uvw : bool, optional
        Note: This is required when visibilities are properly phase rotated in correlator to track the Sun,
        but while creating the MS, UVW values are estimated using a wrong phase center at the start of solar center at the start.

    Returns
    -------
    int
        Success message
    """
    msname = msname.rstrip("/")
    os.system(f"rm -rf {msname}/.solarcenter_move_*")
    print(f"Moving phasecenter to solar center for measurement set: {msname}")
    sun_radec_string, sunra, sundec, sunra_deg, sundec_deg = radec_sun(msname)
    msg = run_chgcenter(
        msname, sunra, sundec, only_uvw=only_uvw, container_name="paircarswsclean"
    )
    if msg != 0:
        print("Phasecenter could not be shifted.")
        os.system(f"touch {msname}/.solarcenter_move_failed")
    else:
        os.system(f"touch {msname}/.solarcenter_move_succeed")
    return msg


def cal_solar_phaseshift(imagename, sigma=10):
    """
    Calculate the difference between solar center and phase center of the image

    Parameters
    ----------
    imagename : str
            Name of the image
    sigma : float
            If Gaussian fitting is not used, threshold for estimating center of mass as solar center (default =10)

    Returns
    -------
    float
            RA of the solar center in degree
    float
            DEC of the solarcenter in degree
    bool
            Whther phase shift required or not. Not required if less than image pixel size
    """
    from scipy.ndimage import center_of_mass

    data = fits.getdata(imagename)
    header = fits.getheader(imagename)
    obstime = header["DATE-OBS"]
    (
        _,
        _,
        _,
        sun_radeg,
        sun_decdeg,
    ) = radec_sun_at_time(obstime)
    cellsize = float(header["CDELT1"]) * 3600.0  # In arcsec
    pix_radius = int((4 * 16 * 60) / cellsize)  # 4 solar radii
    circular_mask = create_circular_mask_array(data[0, 0, ...], pix_radius)
    I_rms = data[0, 0, ...].copy()
    I_rms[circular_mask] = np.nan
    rms = np.nanstd(I_rms)
    I = data[0, 0, ...].copy()
    I[I >= (sigma * rms)] = 1
    I[I < (sigma * rms)] = 0
    cx, cy = center_of_mass(I)
    w = WCS(imagename).celestial
    result = w.array_index_to_world(int(cy), int(cx))
    x_cen = result[0].ra.deg
    y_cen = result[0].dec.deg
    ra = float(x_cen)
    dec = float(y_cen)
    if np.sqrt((ra - sun_radeg) ** 2 + (dec - sun_decdeg) ** 2) < cellsize / 3600.0:
        msg = False
    else:
        msg = True
    return ra, dec, msg


def shift_solarcenter(imagename, sigma=10, overwrite=True):
    """
    Function to shift solar center to image phase center

    Parameters
    ----------
    imagename : str
        Name of the image
    sigma : float, optional
        Sigma threshold for masking solar disk
    overwrite : bool, optional
        Overwrite existing image or not

    Returns
    -------
    int
        Success code 0: Successfully shifted, 1: Shifting is not required, 2: Error in shifting
    """
    ra, dec, shiftsun = cal_solar_phaseshift(imagename, sigma=sigma)
    try:
        if shiftsun:
            w = WCS(imagename).celestial
            pix = w.all_world2pix(np.array([[ra, dec]]), 0)
            ra_pix = int(pix[0][0])
            dec_pix = int(pix[0][1])
            data = fits.getdata(imagename)
            header = fits.getheader(imagename)
            header["CRPIX1"] = float(ra_pix)
            header["CRPIX2"] = float(dec_pix)
            if overwrite:
                fits.writeto(imagename, data=data, header=header, overwrite=True)
            else:
                fits.writeto(imagename.split(".fits")[0]+"_centered.fits", data=data, header=header, overwrite=True)
            msg = 0
        else:
            msg = 1
    except Exception:
        msg = 2
        traceback.print_exc()
    finally:
        return msg


def correct_solar_sidereal_motion(msname="", verbose=False):
    """
    Correct sodereal motion of the Sun

    Parameters
    ----------
    msname : str
        Name of the measurement set

    Returns
    -------
    int
        Success message
    """
    print(f"Correcting sidereal motion for ms: {msname}\n")
    if os.path.exists(msname + "/.sidereal_cor") == False:
        msg = run_solar_sidereal_cor(
            msname=msname, container_name="paircarswsclean", verbose=verbose
        )
        if msg != 0:
            print("Sidereal motion correction is not successful.")
        else:
            os.system("touch " + msname + "/.sidereal_cor")
        return msg
    else:
        print(f"Sidereal motion correction is already done for ms: {msname}")
        return 0
