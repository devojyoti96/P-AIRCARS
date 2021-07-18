import numpy as np,os,julian,smtplib,imaplib,datetime as dtt,psutil,json,urllib.request,copy,time
from casatools import *
from . import access_ms as am
from astropy.io import fits
from astropy.coordinates import EarthLocation,SkyCoord,AltAz
from astropy.time import Time
from astropy import units as u
from datetime import datetime
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.application import MIMEApplication
'''
Code is written by Devojyoti Kansabanik, 05 Jan, 2021
'''
############ Basic functions ###############
datadir = os.path.dirname(__file__)

class ImageBasic:
	'''
	Generic class to calculate different imaging related parameters
	Attribute:
	msname = Name of the measurement set
	'''
	def __init__(self,msname):
		self.msname=msname
		self.md=msmetadata()
		self.md.open(msname)
		self.freq=self.md.meanfreq(0)
		self.md.close()
		self.tb=table()
		self.max_baseline=self.calc_max_baseline()

	############################################
	# Imaging related #
	############################################
	def calc_psf(self,freq=0):
		'''
		Function to calculate PSF size in arcsec
		Parameter :
		freq = Frequency in MHz (default : 0, using central frequency of the ms)
		Return:
		PSF size in arcsec
		'''
		if freq==0:
			wavelength = 299792458.0/self.freq
		else:
			wavelength = 299792458.0/(freq*10**6)
		psf	= (1.42*(wavelength/self.max_baseline))*(180/np.pi*3600.0) # In arcsec
		return psf

	def calc_cellsize(self,num_pixel_in_psf,freq=0):
		'''
		Calculate pixel size in arcsec
		Parameters:
		num_pixel_in_psf = Number of pixels in one PSF
		freq = Frequency in MHz (default : 0, using central frequency of the ms)
		Return:
		Pixel size in arcsec
		'''
		psf	=	self.calc_psf(freq=freq)	
		pixel	=	int(psf/num_pixel_in_psf) 
		return pixel

	def choose_scales(self,num_pixel_in_psf,max_size,freq=0):
		'''
		Function to calculate multiscale scales
		Parameters:
		num_pixel_in_psf = Number of pixels in one PSF
		max_size = Maximum source size in arcsec
		freq = Frequency in MHz (default : 0, using central frequency of the ms)
		Return:
		List of multiscale lists in number of pixels
		'''
		psf=self.calc_psf(freq=freq)
		cellsize=self.calc_cellsize(num_pixel_in_psf,freq=freq)
		max_size_rad=max_size/2.0
		psf_pix	=int(psf/cellsize)
		scale=[0,psf_pix,3*psf_pix,int(max_size_rad/cellsize)]  ### Choosing scale to be [0,psf,3*psf,max_size/5,max_size/3,max_size] in pixel
		if int(max_size_rad/5*cellsize)<max(scale):
			scale.append(int(max_size_rad/5*cellsize))
			if int(max_size_rad/3*cellsize)>int(max_size_rad/5*cellsize) and int(max_size_rad/3*cellsize)<max(scale) :
				scale.append(int(max_size_rad/3*cellsize))
		return sorted(scale)

	def field_of_view(self,freq=0):
		'''
		Calculate optimum field of view in arcsec
		Parameter :
		freq = Frequency in MHz (default : 0, using central frequency of the ms)
		Return:
		Field of view in arcsec
		'''
		if freq==0:
			FOV=np.sqrt(610)*150*10**6/self.freq  # 610 deg^2 is the image FoV at 150MHz for MWA. So extrapolating this to central frequency
		else:
			FOV=np.sqrt(610)*150/freq  # 610 deg^2 is the image FoV at 150MHz for MWA. So extrapolating this to central frequency
		return FOV*3600 ### In arcsecs

	def num_pixels(self,num_pixel_in_psf,freq=0):
		'''
		Number of image pixels
		Parameters:
		num_pixel_in_psf = Number of pixels in one PSF
		freq = Frequency in MHz (default : 0, using central frequency of the ms)
		Return:
		Number of pixels in the image
		'''
		FOV=self.field_of_view(freq=freq)
		cellsize=self.calc_cellsize(num_pixel_in_psf,freq=freq)
		num	=	FOV/cellsize
		pow2	=	int(np.log2(num))
		possibility=	np.array([2**(pow2-1)*3,2**(pow2-2)*5,2**(pow2-2)*7,2**(pow2+1)])
		return possibility[getnearpos(possibility,num)[0]]

	def calc_calib_uvrange(self,max_angular_scale):
		'''
		This function calculate the uvrange to be used for calibration. 
		Maximum uv-range beyond which number of baselines is less than 1% of the total number of baselines are excluded. 
		max_angular_scale = Maximum angular scale to exclude short baselines from calibration
		Return:
		uv-range string, "uvmin~uvmax lambda", "uvmin" in meter, "uvmax" in meter, "uvmin" in wavelength unit, "uvmax" in wavelength unit.
		'''
		wavelength=299792458.0/self.freq
		uvmin_lambda=1.22/np.deg2rad(max_angular_scale)
		uvmin=uvmin_lambda*wavelength
		self.tb.open(self.msname)
		uvw=self.tb.getcol('UVW')
		self.tb.close()
		u,v,w=[uvw[i, :] for i in range(3)]
		uvdist=np.sqrt(u**2+v**2)
		uvlambda=uvdist/wavelength
		uvlambda_hist=np.histogram(uvlambda,bins=int(max(uvlambda)/5))
		max_uvpoints=np.max(uvlambda_hist[1])
		cutpos1=np.min(np.where(uvlambda_hist[0]<max_uvpoints*0.1))
		cutpos2=np.min(np.where(uvlambda_hist[0]==0))
		uvlambda1=uvlambda_hist[1][cutpos1]
		uvlambda2=uvlambda_hist[1][cutpos2]
		if uvlambda1>200 and uvlambda2<200:
			uvmax=uvlambda1
		elif uvlambda1<200 and uvlambda2>200:
			uvmax=uvlambda2
		elif uvlambda1<max(uvlambda_hist[1])*0.7 and uvlambda2<max(uvlambda_hist[1])*0.7:
			uvmax=max(uvlambda_hist[1])*0.9
		else:
			uvmax=max(uvlambda1,uvlambda2)
		uvmax=uvmax*wavelength
		uvmax_lambda=uvmax/wavelength
		return str(int(uvmin_lambda))+'~'+str(int(uvmax_lambda))+'lambda',int(uvmin),int(uvmax),int(uvmin_lambda),int(uvmax_lambda)

	def calc_uvtaper(self):
		'''
		Function return uv-taper value
		Return:
		UV taper string "uvtaper lambda"
		'''
		wavelength=299792458.0/self.freq
		self.tb.open(self.msname)
		uvw=self.tb.getcol('UVW')
		self.tb.close()
		u,v,w=[uvw[i, :] for i in range(3)]
		uvdist=np.sqrt(u**2+v**2)
		uvlambda=uvdist/wavelength
		uvlambda_hist=np.histogram(uvlambda,bins=int(max(uvlambda)/5))
		max_uvpoints=np.max(uvlambda_hist[1])
		cutpos1=np.min(np.where(uvlambda_hist[0]<max_uvpoints*0.1))
		cutpos2=np.min(np.where(uvlambda_hist[0]==0))
		uvlambda1=uvlambda_hist[1][cutpos1]
		uvlambda2=uvlambda_hist[1][cutpos2]
		if uvlambda1>200 and uvlambda2<200:
			uvmax=uvlambda1
		elif uvlambda1<200 and uvlambda2>200:
			uvmax=uvlambda2
		elif uvlambda1<max(uvlambda_hist[1])*0.7 and uvlambda2<max(uvlambda_hist[1])*0.7:
			uvmax=max(uvlambda_hist[1])*0.9
		else:
			uvmax=max(uvlambda1,uvlambda2)
		uvmax=uvmax*wavelength
		uvmax_lambda=uvmax/wavelength
		return str(uvmax_lambda)+'lambda'

	def calc_suntaper(self):
		'''
		Function return uv-taper value to treat Sun as a point source of size 16 arcmin
		Return:
		UV taper string "uvtaper lambda"
		'''
		uvlimstring,umin,umax,uvminlambda,uvmaxlambda=self.calc_calib_uvrange(16/60.0)
		return str(uvminlambda)+'lambda'

	def calc_max_baseline(self):
		'''
		Get the maximum baseline in meter
		Return:
		Maximum baseline length in meter
		'''
		self.tb.open(self.msname)
		uvw=self.tb.getcol('UVW')
		self.tb.close()
		u,v,w=[uvw[i, :] for i in range(3)]
		uvdist=np.sqrt(u**2+v**2)
		return int(np.max(uvdist))

#########################################
# Calibration parameter estimation class
#########################################

class CalcParams:
	def __init__(self,msname,quality_factor,safety_factor):
		self.msname=msname
		self.quality_factor=quality_factor
		self.safety_factor=safety_factor
		self.IB=ImageBasic(self.msname)
		self.AM=am.AccessMS(self.msname)
		
	def calc_calibration_params(self):
		'''
		Function to calculate calibration parameters automatically based on given quality factor and safety standard
		For details read the README
		Return:
		start_sigma,sigma_step,residual_frac,min_sigma,gain_minsnr,DR_delta_rms,DR_delta_neg,min_DR,max_DR,min_selfcal_snr,skip_time,skip_freq,uvrange_to_cal
		'''
		if self.quality_factor==0:
			if self.safety_factor==0:
				start_sigma		=	9.0
				sigma_step		=	1.0
				residual_frac	=	0.2
				min_sigma		=	6
				gain_minsnr		=	3.0		
				DR_delta_rms	=	25
				DR_delta_neg	=	20
				min_DR			=	20
				max_DR			=	100
				min_selfcal_snr	=	2.5
				skip_time		=	960
				skip_freq		=	2560
			elif self.safety_factor==1:
				start_sigma		=	9.0
				sigma_step		=	1.0
				residual_frac	=	0.17
				min_sigma		=	7
				gain_minsnr		=	3.0		
				DR_delta_rms	=	22
				DR_delta_neg	=	18
				min_DR			=	25
				max_DR			=	500
				min_selfcal_snr	=	2.5
				skip_time		=	720
				skip_freq		=	2560
			else:
				start_sigma		=	9.0
				sigma_step		=	1.0
				residual_frac	=	0.15
				min_sigma		=	8
				gain_minsnr		=	3.0		
				DR_delta_rms	=	20
				DR_delta_neg	=	15
				min_DR			=	30
				max_DR			=	1000
				min_selfcal_snr	=	2.5
				skip_time		=	480
				skip_freq		=	2560
		elif self.quality_factor==1:
			if self.safety_factor==0:
				start_sigma		=	10.0
				sigma_step		=	0.5
				residual_frac	=	0.17
				min_sigma		=	7
				gain_minsnr		=	4		
				DR_delta_rms	=	20
				DR_delta_neg	=	15
				min_DR			=	30
				max_DR			=	1000
				min_selfcal_snr	=	3
				skip_time		=	240
				skip_freq		=	1280
			elif self.safety_factor==1:
				start_sigma		=	10.0
				sigma_step		=	0.5
				residual_frac	=	0.15
				min_sigma		=	8
				gain_minsnr		=	4		
				DR_delta_rms	=	18
				DR_delta_neg	=	12
				min_DR			=	35
				max_DR			=	5000
				min_selfcal_snr	=	3
				skip_time		=	120
				skip_freq		=	1280
			else:
				start_sigma		=	10.0
				sigma_step		=	0.5
				residual_frac	=	0.12
				min_sigma		=	9
				gain_minsnr		=	4		
				DR_delta_rms	=	15
				DR_delta_neg	=	10
				min_DR			=	40
				max_DR			=	10000
				min_selfcal_snr	=	3
				skip_time		=	60
				skip_freq		=	1280
		else:
			if self.safety_factor==0:
				start_sigma		=	11.0
				sigma_step		=	0.25
				residual_frac	=	0.15
				min_sigma		=	8
				gain_minsnr		=	4.5		
				DR_delta_rms	=	18
				DR_delta_neg	=	12
				min_DR			=	40
				max_DR			=	10000
				min_selfcal_snr	=	3
				skip_time		=	80
				skip_freq		=	640
			elif self.safety_factor==1:
				start_sigma		=	11.0
				sigma_step		=	0.25
				residual_frac	=	0.1
				min_sigma		=	9
				gain_minsnr		=	4.5		
				DR_delta_rms	=	15
				DR_delta_neg	=	10
				min_DR			=	45
				max_DR			=	50000
				min_selfcal_snr	=	3
				skip_time		=	60
				skip_freq		=	640
			else:
				start_sigma		=	11.0
				sigma_step		=	0.25
				residual_frac	=	0.07
				min_sigma		=	10
				gain_minsnr		=	4.5		
				DR_delta_rms	=	12
				DR_delta_neg	=	8
				min_DR			=	50
				max_DR			=	100000
				min_selfcal_snr	=	3
				skip_time		=	30
				skip_freq		=	640
		uvrange_to_cal=''
		if skip_freq<self.AM.calc_freqres():
			skip_freq=self.AM.calc_freqres()
		if skip_time<self.AM.calc_timeres():
			skip_time=self.AM.calc_timeres()
		return start_sigma,sigma_step,residual_frac,min_sigma,gain_minsnr,DR_delta_rms,DR_delta_neg,min_DR,max_DR,min_selfcal_snr,skip_time,skip_freq,uvrange_to_cal


###########################################
# Flagging related #
###########################################
def calc_flag_fraction(msname):
	'''
	Function to calculate the fraction of total data flagged
	Parameters:
	msname = Name of the measurement set
	Return:
	Fraction of the total data flagged
	'''
	tb=table()
	tb.open(msname)
	flag=tb.getcol('FLAG')
	tb.close()
	flagged_data=np.sum(flag)
	total_data=flag.size
	flagged_fraction=flagged_data/float(total_data)
	return flagged_fraction

def calc_flag_fraction_caltable(caltable):
	'''
	Calculate flagg fraction from caltable
	Parameters:
	caltable = Name of the CASA caltable	
	Return:
	Fraction of the total solutions are flagged
	'''
	#TODO: As of now only based on CASA caltable, later CALIBRATE caltables also included
	
	tb=table()
	tb.open(caltable)
	flag=tb.getcol('FLAG')
	tb.close()
	flagged_fraction=np.sum(flag)/float(flag.size)
	return flagged_fraction

##############################################
# General usuage functions #
##############################################

def getnearpos(array,value):
	'''
	Function to return index of two elements nearest to a given number
	Parameters:
	array = Input numpy array
	value = Value to which nearest element search
	Return:
	Index of two elements nearest to the value
	'''
	a = abs(array-value)
	b=np.argsort(a)
	del a
	if len(b)>1:
		return b[0],b[1]
	else:
		return b[0],b[0]

def error_msgs(err_code):
	'''
	Error code to error message
	'''
	if err_code==1:
		return "Split problem."
	elif err_code==2:
		return "Dirty image is not produced."
	elif err_code==3:
		return "Failed to make image in a selfcal iteration."
	elif err_code==4:
		return "Model image is not present."
	elif err_code==5:
		return "No flux is CLEANed in the model."
	elif err_code==6:
		return "Start sigma is too low to start."
	elif err_code==7:
		return "No good solution found."
	elif err_code==8:
		return "Dynamic range decreasing."
	elif err_code==9:
		return "Maximum selfcal iteration is reached."
	elif err_code==10:
		return 'SNR is not sufficient for selfcal.'
	elif err_code==11:
		return 'No flux above for the present sigma threshold for starting selfcal.'
	elif err_code==12:
		return 'No start sigma and threshold information found for reference time channel.'
	elif err_code==13:
		return 'Maximum selfcal iteration is reached and final image has DR less than minimum DR.'
	elif err_code==100:
		return "Reference time and frequency imaging."
	else:
		return "Succeeded."

def radec_con_deg_to_hhmmss(radeg,decdeg):
	'''
	Convert RA-DEC from degree to hh mm ss dd mm ss format
	Parameter:
	radeg = Value of the RA in degree
	decdeg = Value of the DEC in degree
	Return:
	Numpy array ['RA','DEC'] in hh mm ss dd mm ss format
	'''
	radeg_copy=copy.deepcopy(radeg)
	decdeg_copy=copy.deepcopy(decdeg)
	radeg=abs(radeg)
	decdeg=abs(decdeg)
	ra_sign=int(radeg_copy/radeg)
	dec_sign=int(decdeg_copy/decdeg)
	if ra_sign<0:
		radeg=copy.deepcopy(radeg_copy)+360.0
	rahh=radeg/15.0
	ramm=(rahh-int(rahh))*60.0
	rass=(ramm-int(ramm))*60.0
	ra=str(int(rahh))+'h'+str(int(ramm))+'m'+str('%.2f'%rass)+'s'
	decdd=decdeg
	decmm=(decdd-int(decdd))*60.0
	decss=(decmm-int(decmm))*60.0
	dec=str(int(dec_sign*decdd))+'d'+str(abs(int(decmm)))+'m'+str('%.2f'%abs(decss))+'s'
	return np.array([ra,dec])

def mjdsec_to_timestamp(mjdsec,includedate=True,format=0):
	'''
	Convert CASA MJD seceonds to CASA timestamp
	Parameters:
	mjdsec = CASA MJD seconds in float
	includedate = True, include date in timestamp
	format = Datetime string format (0 : 'YYYY/MM/DD/hh:mm:ss', 1: 'YYYY-MM-DDThh:mm:ss', 2: 'YYYY-MM-DD hh:mm:ss', 3: 'YYYY_MM_DD_hh_mm_ss') 
	Return:
	CASA time stamp in UTC at the format 'YYYY/MM/DD/hh:mm:ss' while includedate=True or 'hh:mm:ss' while includedate=False
	'''
	#TODO : Add more formats
	me=measures()
	qa=quanta()
	today = me.epoch('utc','today')
	mjd = np.array(mjdsec)/86400.0
	today['m0']['value'] =  mjd
	hhmmss = qa.time(today['m0'],prec=8)[0]
	date = qa.splitdate(today['m0'])
	qa.done()
	if (includedate):
		if format==0:
			utcstring = "%s/%02d/%02d/%s" % (date['year'],date['month'],date['monthday'],hhmmss)
		elif format==1:
			utcstring = "%s-%02d-%02dT%s" % (date['year'],date['month'],date['monthday'],hhmmss)
		elif format==2:
			utcstring = "%s-%02d-%02d %s" % (date['year'],date['month'],date['monthday'],hhmmss)
		elif format==3:
			utcstring = "%s_%02d_%02d" % (date['year'],date['month'],date['monthday'])+'_'+'_'.join(hhmmss.split(':'))
	else:
		utcstring = hhmmss
	return utcstring

def timestamp_to_mjdsec(timestamp,format=0):
	'''
	Convert timestamp to mjd
	Parameters:
	timestamp = Time stamp to convert 
	format = Datetime string format (0 : 'YYYY/MM/DD/hh:mm:ss', 1: 'YYYY-MM-DDThh:mm:ss', 2: 'YYYY-MM-DD hh:mm:ss') 	
	Return:
	Return correspondong second of the day
	'''
	if format==0:
		try:
			timestamp_datetime=datetime.strptime(timestamp,'%Y/%m/%d/%H:%M:%S.%f')
		except:
			timestamp_datetime=datetime.strptime(timestamp,'%Y/%m/%d/%H:%M:%S')
	elif format==1:
		try:
			timestamp_datetime=datetime.strptime(timestamp,'%Y-%m-%dT%H:%M:%S.%f')
		except:
			timestamp_datetime=datetime.strptime(timestamp,'%Y-%m-%dT%H:%M:%S')
	elif format==2:
		try:
			timestamp_datetime=datetime.strptime(timestamp,'%Y-%m-%d %H:%M:%S.%f')
		except:
			timestamp_datetime=datetime.strptime(timestamp,'%Y-%m-%d %H:%M:%S')
	mjd=float("{:.2f}".format((julian.to_jd(timestamp_datetime)-2400000.5)*(24.*3600.)))
	return mjd

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

def altaz_to_radec(alt,az,obstime,LAT,LON,ALT):
	'''
	Function to convert altaz to radec for a given Earth location
	Parameters:
	alt = Elevation in degree
	az = Azimuth in degree
	obstime = Time of the observation in 'yyyy-mm-dd hh:mm:ss' format
	LAT = Latitude of the Earth location in degree
	LON = Longitude of Earth location in degree 
	ALT = Altitude of the Earth location in meter
	Return:
	Numpy array ['RA','DEC'] in hh mm ss dd mm ss format
	'''
	LOCATION=EarthLocation.from_geodetic(lat=LAT*u.deg,lon=LON*u.deg,height=ALT*u.m)
	observing_time=Time(obstime) 
	AltAz=SkyCoord(alt=alt*u.deg,az=az*u.deg,obstime=observing_time,frame='altaz',location=LOCATION)
	radec=AltAz.icrs
	radeg=radec.ra.deg
	decdeg=radec.dec.deg
	radec_str=radec_con_deg_to_hhmmss(radeg,decdeg)
	return radec_str

def freq_to_MWA_coarse(freq):
	'''
	Frequency to MWA coarse channel conversion
	Parameters:
	freq = Frequency in MHz
	Return:
	MWA coarse channel number
	'''
	freq=float(freq)
	coarse_chans=[[(i*1.28)-0.64,(i*1.28)+0.64] for i in range(300)]
	for i in range(len(coarse_chans)):
		ch0=coarse_chans[i][0]
		ch1=coarse_chans[i][1]
		if freq>=ch0 and freq<ch1:
			return i 

def update_mwa_obsids(obsid_file=''):
	'''
	Function to update MWA OBSIDs 
	Parameter:
	obsid_file = Name of the file to save MWA OBSIDs
	Return:
	OBSID file name, update code 
	'''
	print ('Updating local MWA OBSid file......\n')
	if obsid_file=='':
		obsid_file=datadir+'/MWA_OBSids'
	BASEURL='http://ws.mwatelescope.org/'
	temp_array=np.empty(0,dtype='int')
	if os.path.isfile(obsid_file+'.npy')==True:
		try:
			obsids=np.load(obsid_file+'.npy',allow_pickle=True)
			start_obsid=np.max(obsids)
			temp_array=np.append(temp_array,obsids)
		except:
			start_obsid=972654120
	else:
		start_obsid=972654120
	end_obsid=3786480018  # Till 2100-01-01
	searchurl=BASEURL+'metadata/find?maxtime='+str(end_obsid)+'&page=20000000000000'
	try:
		end_obsid=json.load(urllib.request.urlopen(searchurl,timeout=10))[-1][0]
		print ('Last OBSID in MWA metadata server : '+end_obsid+'\n')
		while True:
			searchurl=BASEURL+'metadata/find?mintime='+str(start_obsid)+'&maxtime='+str(start_obsid+432000)
			try:
				OBSid=json.load(urllib.request.urlopen(searchurl,timeout=150))
				OBSid=np.array(OBSid)[:,0].astype('int')
				start_obsid=np.max(OBSid)+235
			except:
				OBSid=np.empty(0,dtype='int')
				start_obsid=start_obsid+3600
			if len(OBSid)!=0:
				temp_array=np.append(temp_array,OBSid)
			if start_obsid>=end_obsid:
				break
		np.save(obsid_file,temp_array)
		return obsid_file+'.npy',0
	except:
		return obsid_file+'.npy',1


def get_OBSID_from_ms(msname):
	'''
	Function to return OBSID of an MWA observation
	Parameters:
	msname = Name of the measurement set
	Return:
	MWA OBSID
	'''
	obsid_file=datadir+'/MWA_OBSids.npy'
	if os.path.isfile(obsid_file)==True:
		obsids=np.load(obsid_file,allow_pickle=True)
	else:
		obsids=np.empty(0)
	if len(obsids)!=0:
		md=msmetadata()
		md.open(msname)
		obs_mjd_ms=md.timerangeforobs(0)['begin']['m0']['value']
		md.close()
		utc_string=mjdsec_to_timestamp(obs_mjd_ms*24*3600,includedate=True,format=1)	
		GPStime=int(Time(utc_string,format='isot',scale='utc').gps)
		OBSid=obsids[np.argmin(abs(GPStime-obsids))]
		return OBSid
	else:
		c=0
		while c<5:
			try:
				BASEURL='http://ws.mwatelescope.org/'
				
				ms_path=os.path.dirname(os.path.realpath(msname))
				searchurl=BASEURL+'metadata/tconv/?utciso='+utc_string
				GPStime=json.load(urllib.request.urlopen(searchurl,timeout=10))
				searchurl=BASEURL+'metadata/find?maxtime='+str(GPStime)+'&mintime='+str(GPStime-500)+'&page=20000000000000'
				OBSid=json.load(urllib.request.urlopen(searchurl,timeout=15))[-1][0]
				return OBSid
			except Exception as e:
				print ('Error occured : '+str(e)+'\n')
				c+=1
				time.sleep(5.0)
				if c>=5:
					return 0

def get_OBSID_from_metafits(metafits):
	'''
	Function to return OBSID of an MWA observation
	Parameters:
	metafits = Name of the metafits file
	MWA OBSID
	'''
	OBSid=fits.getheader(metafits)['GPSTIME']
	return OBSid 

def download_metafits(msname,outdir):
	'''
	Function to download MWA metafits of a given measurement set.
	Parameters :
	msname : Name of the measurement set
	outdir : Name of the outputdir
	Return :
	Download the metafits file of the given measurement set and return metafits file name
	'''
	BASEURL='http://ws.mwatelescope.org/'
	OBSid=get_OBSID_from_ms(msname)
	if os.path.isfile(outdir+'/'+str(OBSid)+'.metafits')==False:
		try:
			os.system('wget -O '+outdir+'/'+str(OBSid)+'.metafits http://ws.mwatelescope.org/metadata/fits?obs_id='+str(OBSid))
			metafits=outdir+'/'+str(OBSid)+'.metafits'
		except:
			metafits=None
	else:
		metafits=outdir+'/'+str(OBSid)+'.metafits'
	return metafits

def compress_files(filelist,outputfile):
	'''
	Compress a list of numpy files
	Parameters:
	filelist = List of numpy table files
	outputfile = Output compressed file name
	Return:
	Compressed file name (Compressed file have arrays in format [original_filename,data_array])
	'''
	file_array_list=[]
	file_name_list=[]
	for i in filelist:
		a=np.load(i,allow_pickle=True)
		a=np.append(i,a)
		file_array_list.append(a)
	np.savez_compressed(outputfile,a=np.array(file_array_list))
	os.system('mv '+outputfile+'.npz '+outputfile)
	return outputfile








