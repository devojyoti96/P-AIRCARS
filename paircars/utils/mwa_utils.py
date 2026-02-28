import psutil
import numpy as np
import glob
import os
import traceback
import warnings
import astropy.units as u
import requests
import math
from astropy.wcs import FITSFixedWarning
from astropy.io import fits
from astropy.time import Time
from casatools import msmetadata
from .udocker_utils import run_wsclean

warnings.simplefilter("ignore", category=FITSFixedWarning)


def get_MWA_OBSID(msname):
    """
    Get MWA OBSID from ms

    Parameters
    ----------
    msname : str
        Measurement set

    Returns
    -------
    int
        OBSid
    """
    msmd = msmetadata()
    msmd.open(msname)
    start_time = msmd.timerangeforobs(0)["begin"]["m0"]["value"] * 86400
    msmd.close()
    t = Time(start_time * u.s, format="mjd", scale="utc")
    gps = t.gps
    obsid = int((gps // 8) * 8)
    return obsid


def get_ncoarse(msname):
    """
    Get number of coarse channels

    Parameters
    ----------
    msname : str
        Measurement set

    Returns
    -------
    int
        Number of coarse channels
    """
    msmd = msmetadata()
    msmd.open(msname)
    freqs = msmd.chanfreqs(0, unit="MHz")
    bw = max(freqs) - min(freqs)
    ncoarse = max(1, math.ceil(bw / 1.28))
    return ncoarse


def freq_to_MWA_coarse(freq):
    """
    Frequency to MWA coarse channel conversion.

    Parameters
    ----------
    freq : float
        Frequency in MHz

    Returns
    -------
    int
        MWA coarse channel number
    """
    return int(round(freq / 1.28))


def get_MWA_coarse_chan(msname):
    """
    Get MWA coarse channel number

    Parameters
    ----------
    msname : str
        Measurement set

    Returns
    -------
    int
        Coarse channel corresponding to central frequency of the measurement set
    """
    msmd = msmetadata()
    msmd.open(msname)
    meanfreq = msmd.meanfreq(0, unit="MHz")
    msmd.close()
    ncoarse = freq_to_MWA_coarse(meanfreq)
    return ncoarse


def get_MWA_coarse_bands(msname, flag_central_chan=False):
    """
    Get MWA coarse channel bands.

    Parameters
    ----------
    msname : str
        Name of the measurement set
    flag_central_chan : bool, optional
        Flag central channel or not

    Returns
    -------
    list of tuples
        (start_chan, end_chan, good_chan_list)
    """
    bad_spw = get_bad_chans(msname, flag_central_chan=flag_central_chan)
    if bad_spw:
        bad_chans = set(int(i) for i in bad_spw.split("0:")[1].split(";"))
    else:
        bad_chans = set()
    msmd = msmetadata()
    msmd.open(msname)
    freqs = msmd.chanfreqs(0, unit="MHz")
    freqres = msmd.chanres(0, unit="MHz")[0]
    msmd.close()
    msmd.done()
    nchan = len(freqs)
    nchan_coarse = int(round(1.28 / freqres))
    start_ms_freq = np.nanmin(freqs)
    end_ms_freq = np.nanmax(freqs)

    coarse_chans = []
    seen = set()
    for start_chan in range(0, nchan, nchan_coarse):
        end_chan = min(start_chan + nchan_coarse - 1, nchan - 1)
        cent_freq = np.nanmean(freqs[start_chan : end_chan + 1])
        if cent_freq < start_ms_freq or cent_freq > end_ms_freq:
            continue
        good_chunk = [
            ch for ch in range(start_chan, end_chan + 1) if ch not in bad_chans
        ]
        entry = (start_chan, end_chan, tuple(good_chunk))
        if entry not in seen:
            coarse_chans.append((start_chan, end_chan, good_chunk))
            seen.add(entry)
    return coarse_chans


def get_bad_chans(msname, flag_central_chan=False):
    """
    Get bad channels to flag

    Parameters
    ----------
    msname : str
        Name of the ms
    flag_central_chan : bool, optional
        Flag central channel

    Returns
    -------
    str
        SPW string of bad channels
    """
    msmd = msmetadata()
    msmd.open(msname)
    chanres = msmd.chanres(0, unit="MHz")[0]  # MHz
    nchan = msmd.nchan(0)
    msmd.close()
    msmd.done()
    bw = nchan*chanres
    if chanres > 0.16 or bw<1.28:
        return ""
    n_per_coarse = int(round(1.28 / chanres))
    n_edge = max(1, int(round(0.16 / chanres)))
    bad_channels = set()
    for start in range(0, nchan, n_per_coarse):
        coarse_end = min(start + n_per_coarse - 1, nchan - 1)
        # First 160 kHz
        for ch in range(start, min(start + n_edge, coarse_end + 1)):
            bad_channels.add(ch)
        # Last 160 kHz
        for ch in range(max(coarse_end - n_edge + 1, start), coarse_end + 1):
            bad_channels.add(ch)
        if flag_central_chan:
            # Central channel
            central_chan = start + (coarse_end - start) // 2
            bad_channels.add(central_chan)
    if not bad_channels:
        return ""
    # Sort and format
    sorted_chans = sorted(bad_channels)
    chan_string = ";".join(str(ch) for ch in sorted_chans)
    return f"0:{chan_string}"


def get_good_chans(msname):
    """
    Get good channel range of MWA

    Parameters
    ----------
    msname : str
        Name of the ms

    Returns
    -------
    str
        SPW string
    """
    msmd = msmetadata()
    msmd.open(msname)
    nchan = msmd.nchan(0)
    msmd.close()
    msmd.done()
    bad_spw = get_bad_chans(msname)
    if bad_spw == "":
        good_spw = f"0:0~{nchan-1}"
    else:
        bad_chan_list = bad_spw.split("0:")[-1].split(";")
        good_chan_list = []
        start_chan = 0
        for bad_chans in bad_chan_list:
            end_chan = int(bad_chans.split("~")[0])
            if end_chan > start_chan:
                good_chan_list.append(f"{start_chan+1}~{end_chan-1}")
            start_chan = int(bad_chans.split("~")[-1])
        good_chan_list.append(f"{start_chan+1}~{nchan-1}")
        good_spw = f"0:{';'.join(good_chan_list)}"
    return good_spw


def get_mwa_bad_ants(metafits):
    """
    Function to determine non-working MWA tiles for a observation

    Parameters
    ----------
    metafits : str
        Name of the metafits file

    Returns
    -------
    str
        Non-working antenna names
    """
    data = fits.getdata(metafits)
    flags = np.array(data["Flag"])
    tiles = np.array(data["TileName"])
    pos = np.where(flags == 1)
    bad_tiles = tiles[pos]
    bad_tiles = np.unique(bad_tiles)
    bad_antennas = ""
    if len(bad_tiles) > 0:
        for ant in bad_tiles:
            bad_antennas += str(ant) + ","
        bad_antennas = bad_antennas[:-1]
    return bad_antennas


def download_MWA_metafits(OBSID, outdir="."):
    """
    Download MWA metafits file for a given OBSID.

    Parameters
    ----------
    OBSID : int
        MWA observation ID
    outdir : str
        Output directory

    Returns
    -------
    str or None
        Path to metafits file or None if failed
    """
    os.makedirs(outdir, exist_ok=True)
    metafits = os.path.join(outdir, f"{OBSID}.metafits")
    if os.path.isfile(metafits):
        return metafits
    url = f"https://ws.mwatelescope.org/metadata/fits?obs_id={OBSID}"
    for attempt in range(5):
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 200:
                with open(metafits, "wb") as f:
                    f.write(r.content)
                return metafits
        except Exception:
            pass
    print(f"Metafits file could not be downloaded after {max_tries} tries.")
    return None
