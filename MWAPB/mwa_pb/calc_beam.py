import numpy as np,sys,struct,time,os
from . import mwa_sweet_spots
from . import config
from . import primary_beam
from . import beam_tools
from . import mwa_sweet_spots as mss
from casatools import *
from casatasks import *
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import EarthLocation,SkyCoord,AltAz
from astropy import units as u
from numpy.linalg import inv
from optparse import OptionParser
'''
Code is written by Devojyoti Kansabanik, 24 April, 2021
'''
def freq_to_MWA_coarse(freq):
	'''
	Frequency to MWA coarse channel conversion
	Parameters:
	freq = Frequency in MHz
	Return:
	MWA coarse channel number
	'''
	coarse_chans=[[(i*1.28)-0.64,(i*1.28)+0.64] for i in range(300)]
	for i in range(len(coarse_chans)):
		ch0=coarse_chans[i][0]
		ch1=coarse_chans[i][1]
		if freq>=ch0 and freq<ch1:
			return i

def radec_to_altaz(ra,dec,obstime,LAT,LON,ALT):
	'''
	Function to convert radec to altaz for a given Earth location
	Parameters:
	ra = RA either in degree or 'hh:mm:ss' or '%fh%fm%fs' format
	dec = DEC either in degree or 'dd:mm:ss' or '%fd%fm%fs'format
	obstime = Time of the observation in 'yyyy-mm-dd hh:mm:ss' format
	LAT = Latitude of the Earth location in degree
	LON = Longitude of Earth location in degree 
	ALT = Altitude of the Earth location in meter
	Return:
	Elevation, Azimuth in degree
	'''
	LOCATION=EarthLocation.from_geodetic(lat=LAT*u.deg,lon=LON*u.deg,height=ALT*u.m)
	observing_time=Time(obstime)  
	aa=AltAz(location=LOCATION,obstime=observing_time)
	try:
		ra=float(ra)
		dec=float(dec)
		coord=SkyCoord(ra,dec,frame='icrs',unit='deg')
	except:
		try:
			coord=SkyCoord(ra,dec)
		except:
			coord=SkyCoord(ra,dec,unit=(u.hourangle,u.deg))
	altaz_object=coord.transform_to(aa)
	alt=altaz_object.alt.degree
	az=altaz_object.az.degree
	return alt,az


def calc_beamjones_info(ra=0,dec=0,obstime='2014-05-04 02:48:00',gridpoint=0,freq=[150],metafits='',coord='MWA'):
	'''
	Function to calculate the Beam Jones at specific ra, dec and gridpoint at certain epoch for a list of frequencies.
	Parameters:
	ra = RA in degree (default : 0)
	dec = DEC in degree (default : 0)
	obstime = Observing epoch in 'yyyy-mm-dd hh:mm:ss' format
	gridpoint = Grid point number (default : 0, if metafits is given automatically read it from there)
	freq = [], list of frequencies 
	metafits = Metafits file
	coord = MWA or IAU convention (default : MWA)
	Return :
	List of beam jones and Stokes I beam value corresponding to frequency list
	'''
	LAT=-26.703319
	LON=116.67081
	ALT=377.0
	if type(freq)!=list:
		freq=[freq]
	coarse_freqs=np.array([freq_to_MWA_coarse(i)*float(1.28) for i in freq])
	alt,az=radec_to_altaz(ra,dec,obstime,LAT,LON,ALT)
	if metafits!='':
		metaheader=fits.getheader(metafits)
		gridpoint=metaheader['GRIDNUM']
		delays=[int(h) for h in metaheader['DELAYS'].split(',')]
		delays = np.vstack((delays, delays))
	else:
		delays=mss.get_delays(gridpoint)
	beam_jones=[]
	stokesI_beam=[]
	prefreq=0
	for f in freq:
		nearest_coarse_freq=coarse_freqs[np.argmin(abs(coarse_freqs-f))]
		if nearest_coarse_freq!=prefreq:
			f=f*10**6
			Jones_FullEE=primary_beam.MWA_Tile_full_EE(np.deg2rad(90-alt),np.deg2rad(az),f,delays=delays,zenithnorm=True,jones=True,interp=True)
			if coord=='MWA':
				coord_transform=np.matrix([[0.0,-1.0],[-1.0,0.0]])
				Jones_FullEE=np.matmul(coord_transform,np.matrix(Jones_FullEE))
			stokesI=((np.abs(Jones_FullEE[0,0]))**2+(np.abs(Jones_FullEE[0,1]))**2+(np.abs(Jones_FullEE[1,0]))**2+(np.abs(Jones_FullEE[1,1]))**2)/2
			prefreq=nearest_coarse_freq
		else:
			Jones_FullEE=beam_jones[-1]
			stokesI=stokesI_beam[-1]
		beam_jones.append(Jones_FullEE)
		stokesI_beam.append(stokesI)
	return beam_jones,stokesI_beam














