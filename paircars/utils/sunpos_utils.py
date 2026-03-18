import astropy.units as u
import os
import traceback
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
from .imaging import calc_sun_dia

#####################################
# Sun position related
#####################################
datadir = get_datadir()
try:
    solar_system_ephemeris.set(f"{datadir}/de440s")
except Exception:
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
        Whether phase shift required or not. Not required if less than image pixel size
    """
    def gaussian_2d(xy, amplitude, x0, y0, sigma_x, sigma_y, offset):
        x, y = xy
        g = offset + amplitude * np.exp(
            -(
                ((x - x0) ** 2) / (2 * sigma_x**2)
                + ((y - y0) ** 2) / (2 * sigma_y**2)
            )
        )
        return g.ravel()
    data = fits.getdata(imagename)
    header = fits.getheader(imagename)
    obstime = header["DATE-OBS"]
    if header["CTYPE3"]=="FREQ":
        freqMHz = float(header["CRVAL3"])/10**6 # In MHz
        sun_dia = calc_sun_dia(freqMHz) # In arcmin
    elif header["CTYPE4"]=="FREQ":
        freqMHz = float(header["CRVAL4"])/10**6 # In MHz
        sun_dia = calc_sun_dia(freqMHz) # In arcmin
    else:
        sun_dia = 32 # In arcmin
    (
        _,
        _,
        _,
        sun_radeg,
        sun_decdeg,
    ) = radec_sun_at_time(obstime)
    cellsize = float(abs(header["CDELT1"])) * 3600.0  # In arcsec
    imsize = int(header["NAXIS1"]) # Image size
    pix_radius = min(imsize, int((4 * 16 * 60) / cellsize))  # 4 solar radii
    if data.ndim == 4:
        data2d = data[0, 0, ...]
    elif data.ndim == 3:
        data2d = data[0,...]
    else:
        data2d = data        
    circular_mask = create_circular_mask_array(data, pix_radius)
    try:
        from scipy.optimize import curve_fit
        from scipy.ndimage import gaussian_filter
        data2d = gaussian_filter(data2d, sigma=3)
        max_pos = np.where(data2d==np.nanmax(data2d))
        y0, x0 = max_pos[0][0], max_pos[1][0]  
        y_min = max(0, y0 - pix_radius)
        y_max = min(data2d.shape[0], y0 + pix_radius)
        x_min = max(0, x0 - pix_radius)
        x_max = min(data2d.shape[1], x0 + pix_radius)
        y_grid, x_grid = np.mgrid[y_min:y_max, x_min:x_max]
        subdata = data2d[y_min:y_max, x_min:x_max]
        base_mean = np.nanmean(data2d[~circular_mask])
        sigma = int((sun_dia/2)*60.0/cellsize) 
        p0 = [np.nanmax(subdata), x0, y0, sigma, sigma, base_mean]
        popt, pcov = curve_fit(gaussian_2d,(x_grid, y_grid),data_region.ravel(),p0=p0,maxfev=5000)
        apparent_pix_x = int(popt[1])
        apparent_pix_y = int(popt[2])
    except Exception:
        traceback.print_exc()
        print("Using imsmooth")
        from casatasks import imsmooth, exportfits
        imsmooth(imagename=imagename,outfile=f"{imagename}.smoothed",targetres=True,beam={"major":f"{sun_dia}arcmin","minor":f"{sun_dia}arcmin","pa":"0deg"},overwrite=True)
        exportfits(imagename=f"{imagename}.smoothed",fitsimage=f"{imagename}.smoothed.fits",overwrite=True)
        os.system(f"rm -rf {imagename}.smoothed")
        data_smoothed = fits.getdata(f"{imagename}.smoothed.fits")
        os.system(f"rm -rf {imagename}.smoothed.fits")
        if data_smoothed.ndim == 4:
            data2d_smoothed = data_smoothed[0, 0, ...]
        elif data.ndim == 3:
            data2d_smoothed = data_smoothed[0,...]
        else:
            data2d_smoothed = data_smoothed
        max_pos = np.where(data2d_smoothed==np.nanmax(data2d_smoothed))
        apparent_pix_y, apparent_pix_x = max_pos[0][0], max_pos[1][0] 
    try:
        w = WCS(imagename).celestial
        result = w.array_index_to_world(apparent_pix_y, apparent_pix_x)
        x_cen = result.ra.deg
        y_cen = result.dec.deg
        ra = float(x_cen)
        dec = float(y_cen)
        print(f"Aparent RA DEC: {ra} {dec}")
        print(f"True RA, DEC: {sun_radeg}, {sun_decdeg}")
        if np.sqrt((ra - sun_radeg) ** 2 + (dec - sun_decdeg) ** 2) < cellsize / 3600.0:
            need_shifting = False
        else:
            need_shifting = True
        return 0, need_shifting, ra, dec, sun_radeg, sun_decdeg, apparent_pix_x, apparent_pix_y
    except Exception:
        traceback.print_exc()
        return 1, False, sun_radeg, sun_decdeg, sun_radeg, sun_decdeg, 0, 0


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
    sunra, sundec, shiftsun = cal_solar_phaseshift(imagename, sigma=sigma)
    try:
        if shiftsun:
            w = WCS(imagename).celestial
            pix = w.all_world2pix(np.array([[sunra, sundec]]), 0)
            ra_pix = int(np.round(pix[0][0]))
            dec_pix = int(np.round(pix[0][1]))
            data = fits.getdata(imagename)
            header = fits.getheader(imagename)
            header["CRPIX1"] = float(ra_pix + 1)
            header["CRPIX2"] = float(dec_pix + 1)
            header["CRVAL1"] = float(sunra)
            header["CRVAL2"] = float(sundec)
            if overwrite:
                fits.writeto(imagename, data=data, header=header, overwrite=True)
            else:
                fits.writeto(
                    imagename.split(".fits")[0] + "_centered.fits",
                    data=data,
                    header=header,
                    overwrite=True,
                )
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
    if not os.path.exists(msname + "/.sidereal_cor"):
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
