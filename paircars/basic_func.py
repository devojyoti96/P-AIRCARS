import os
import numpy as np,julian,smtplib,imaplib,datetime as dtt,psutil,json,urllib.request,copy,time,requests
from casatools import msmetadata,table,measures,quanta,agentflagger,image,calibrater,ms
from casatasks import tclean,imhead
from . import access_ms as am
from astropy.io import fits
from astropy.coordinates import EarthLocation,SkyCoord,AltAz
from astropy.time import Time
from astropy import units as u
from datetime import datetime
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.application import MIMEApplication
from pathlib import Path
import matplotlib.pyplot as plt
'''
Code is written by Devojyoti Kansabanik, 05 Jan, 2021
'''
############ Basic functions ###############
datadir = os.path.dirname(__file__)

class ImageBasic:
	'''
	Generic class to calculate different imaging related parameters

	Parameters
	----------
	msname : str 
		Name of the measurement set
	includeflag : bool
			Include flag data or not whilw calculating maximum baseline
	'''
	def __init__(self,msname,includeflag=False):
		self.msname=msname
		self.md=msmetadata()
		self.md.open(msname)
		self.freq=self.md.meanfreq(0)
		self.md.close()
		self.tb=table()
		self.max_baseline=self.calc_max_baseline(includeflag=includeflag)

	############################################
	# Imaging related #
	############################################
	def calc_psf(self,freq=0):
		'''
		Function to calculate PSF size in arcsec

		Parameters
		----------
		freq : float
	 		Frequency in MHz (default : 0, using central frequency of the ms)
		Returns
		-------
		float
			PSF size in arcsec
		'''
		if freq==0:
			wavelength = 299792458.0/self.freq
		else:
			wavelength = 299792458.0/(freq*10**6)
		psf	= (1.2*(wavelength/self.max_baseline))*(180/np.pi*3600.0) # In arcsec
		return psf

	def calc_cellsize(self,num_pixel_in_psf,freq=0):
		'''
		Calculate pixel size in arcsec

		Parameters
		----------
		num_pixel_in_psf : int 
			Number of pixels in one PSF
		freq : float 
			Frequency in MHz (default : 0, using central frequency of the ms)
		Returns
		-------
		float
			Pixel size in arcsec
		'''
		psf	=	self.calc_psf(freq=freq)	
		pixel	=	int(psf/num_pixel_in_psf) 
		return pixel

	def choose_scales(self,num_pixel_in_psf,max_size,freq=0):
		'''
		Function to calculate multiscale scales

		Parameters
		----------
		num_pixel_in_psf : int 
			Number of pixels in one PSF
		max_size : float 
			Maximum source size in arcsec
		freq : float 
			Frequency in MHz (default : 0, using central frequency of the ms)
		Returns
		-------
		list
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

		Parameters
		----------
		freq : float 
			Frequency in MHz (default : 0, using central frequency of the ms)
		Returns
		-------
		float
			Field of view in arcsec
		'''
		if freq==0:
			FOV=np.sqrt(610)*150*10**6/self.freq  # 600 deg^2 is the image FoV at 150MHz for MWA. So extrapolating this to central frequency
		else:
			FOV=np.sqrt(610)*150/freq  # 610 deg^2 is the image FoV at 150MHz for MWA. So extrapolating this to central frequency
		return FOV*3600 ### In arcsecs

	def num_pixels(self,num_pixel_in_psf,freq=0):
		'''
		Number of image pixels
	
		Parameters
		----------
		num_pixel_in_psf : int
			Number of pixels in one PSF
		freq : float 
			Frequency in MHz (default : 0, using central frequency of the ms)
		Returns
		-------
		int
			Number of pixels in the image
		'''
		FOV=self.field_of_view(freq=freq)
		cellsize=self.calc_cellsize(num_pixel_in_psf,freq=freq)
		num=FOV/cellsize
		pow2=round(np.log2(num),0)
		possibility=np.array([2**(pow2),2**(pow2)*2,2**(pow2)*3,2**(pow2)*4,2**(pow2)*5,2**(pow2)*6,2**(pow2)*7,2**(pow2)*8,2**(pow2)*9,2**(pow2)*10])
		argmin=np.argmin(np.abs(possibility-num))
		return int(possibility[argmin])

	def calc_calib_uvrange(self,max_angular_scale,uvbin=10,includeflag=False):
		'''
		This function calculate the uvrange to be used for calibration. 
	
		Parameters
		----------
		max_angular_scale : float 
			Maximum angular scale to exclude short baselines from calibration
		uvbin : float
			Binning in uv-lambda
		includeflag : bool
			Include flag data or not whilw calculating maximum baseline
		Returns
		-------
		str
			uv-range string (uvmin~uvmax lambda)
		float
			uvmin in meter
		float
			uvmax in meter
		float
			uvmin in wavelength unit
		float
			uvmax in wavelength unit
		'''
		wavelength=299792458.0/self.freq
		uvmin_lambda=1.22/np.deg2rad(max_angular_scale)
		uvmin=uvmin_lambda*wavelength
		self.tb.open(self.msname)
		uvw=self.tb.getcol('UVW')
		if includeflag==False:
			flag=self.tb.getcol('FLAG')
			flag=np.prod(np.prod(flag,axis=0,dtype='bool'),axis=0,dtype='bool')
			for i in range(3):
				uvw[i][flag]=np.nan
		self.tb.close()
		u,v,w=[uvw[i, :] for i in range(3)]
		uvdist=np.sqrt(u**2+v**2)
		nanpos=np.where(np.isnan(uvdist)==True)
		uvdist=np.delete(uvdist,nanpos)
		uvlambda=uvdist/wavelength
		uvlambda_hist=np.histogram(uvlambda,bins=int(max(uvlambda)/uvbin))
		max_uvpoints=np.max(uvlambda_hist[1])
		cutpos1=np.min(np.where(uvlambda_hist[0]<max_uvpoints*0.1))
		uvlambda1=uvlambda_hist[1][cutpos1]
		try:
			cutpos2=np.min(np.where(uvlambda_hist[0]==0))
			uvlambda2=uvlambda_hist[1][cutpos2]
		except:
			uvlambda2=max(uvlambda_hist[1])*0.66
		if uvlambda1>200 and uvlambda2<200:
			uvmax=uvlambda1
		elif uvlambda1<200 and uvlambda2>200:
			uvmax=uvlambda2
		elif uvlambda1<max(uvlambda_hist[1])*0.7 and uvlambda2<max(uvlambda_hist[1])*0.7:
			uvmax=max(uvlambda_hist[1])*0.7
		else:
			uvmax=max(uvlambda1,uvlambda2)
		uvmax=uvmax*wavelength
		uvmax_lambda=uvmax/wavelength
		return str(int(uvmin_lambda))+'~'+str(int(uvmax_lambda))+'lambda',int(uvmin),int(uvmax),int(uvmin_lambda),int(uvmax_lambda)

	def calc_uvtaper(self,uvbin=10,includeflag=False):
		'''
		Function return uv-taper value
		Parameters
		----------
		uvbin : float
			Binning in uv-lambda
		includeflag : bool
			Include flag data or not whilw calculating maximum baseline
		Returns
		-------
		str
			UV taper (uvtaper lambda)
		'''
		wavelength=299792458.0/self.freq
		self.tb.open(self.msname)
		uvw=self.tb.getcol('UVW')
		if includeflag==False:
			flag=self.tb.getcol('FLAG')
			flag=np.prod(np.prod(flag,axis=0,dtype='bool'),axis=0,dtype='bool')
			for i in range(3):
				uvw[i][flag]=np.nan
		self.tb.close()
		u,v,w=[uvw[i, :] for i in range(3)]
		uvdist=np.sqrt(u**2+v**2)
		nanpos=np.where(np.isnan(uvdist)==True)
		uvdist=np.delete(uvdist,nanpos)
		uvlambda=uvdist/wavelength
		uvlambda_hist=np.histogram(uvlambda,bins=int(max(uvlambda)/uvbin))
		max_uvpoints=np.max(uvlambda_hist[1])
		cutpos1=np.min(np.where(uvlambda_hist[0]<max_uvpoints*0.1))
		uvlambda1=uvlambda_hist[1][cutpos1]
		try:
			cutpos2=np.min(np.where(uvlambda_hist[0]==0))
			uvlambda2=uvlambda_hist[1][cutpos2]
		except:
			uvlambda2=max(uvlambda_hist[1])*0.66
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

	def calc_suntaper(self,uvbin=10,includeflag=False):
		'''
		Function return uv-taper value to treat Sun as a point source of size 16 arcmin
		Parameters
		----------
		uvbin : float
			Binning in uv-lambda
		includeflag : bool
			Include flag data or not whilw calculating maximum baseline
		Returns
		-------
		str
			UV taper (uvtaper lambda)
		'''
		uvlimstring,umin,umax,uvminlambda,uvmaxlambda=self.calc_calib_uvrange(16/60.0,uvbin=uvbin,includeflag=includeflag)
		return str(uvminlambda)+'lambda'

	def calc_max_baseline(self,includeflag=False):
		'''
		Get the maximum baseline in meter
		Parameters
		----------
		includeflag : bool
			Include flag data or not whilw calculating maximum baseline
		Returns
		-------
		float
			Maximum baseline length in meter
		'''
		self.tb.open(self.msname)
		uvw=self.tb.getcol('UVW')
		if includeflag==False:
			flag=self.tb.getcol('FLAG')
			flag=np.prod(np.prod(flag,axis=0,dtype='bool'),axis=0,dtype='bool')
			for i in range(3):
				uvw[i][flag]=np.nan
		self.tb.close()
		u,v,w=[uvw[i, :] for i in range(3)]
		uvdist=np.sqrt(u**2+v**2)
		return int(np.nanmax(uvdist))

	def calc_max_uvw(self):
		'''
		Get the maximum u,v,w

		Returns
		-------
		float
			Maximum U
		float
			Maximum V
		float 
			Maximum W
		'''
		self.tb.open(self.msname)
		uvw=self.tb.getcol('UVW')
		flag=self.tb.getcol('FLAG')
		flag=np.prod(np.prod(flag,axis=0,dtype='bool'),axis=0,dtype='bool')
		self.tb.close()
		for i in range(3):
			uvw[i][flag]=np.nan
		u,v,w=[uvw[i, :] for i in range(3)]
		maxu=np.nanmax(np.abs(u))
		maxv=np.nanmax(np.abs(v))
		maxw=np.nanmax(np.abs(w))
		return maxu,maxv,maxw

#########################################
# Calibration parameter estimation class
#########################################

class CalcParams:
	'''
	Calculate calibration parameters	

	Parameters
	----------
	msname : str
		Name of the measurement set
	quality_factor : int 
		Imaging quality factor (0,1,2)
	safety_factor : int
		Safety factor for calibration (0,1,2)
	'''
	def __init__(self,msname,quality_factor,safety_factor):
		self.msname=msname
		self.quality_factor=quality_factor
		self.safety_factor=safety_factor
		self.IB=ImageBasic(self.msname)
		self.AM=am.AccessMS(self.msname)
		
	def calc_calibration_params(self):
		'''
		Return calibration parameters automatically based on given quality factor and safety standard
		For details read the documentation of quality_factor and safety_factor

		Returns
		-------
		float
			Starting sigma value for calibration
		float
			Step to reduce sigma value with self-calibration
		float
			Residual flux fraction to stop the calibration
		float
			Minimum allowed sigma for threshold
		float
			Minimum gain SNR
		float
			Increment of DR_rms
		float
			Increment of DR_neg
		float
			Minimum allowed dynamic range
		float
			Maximum dynamic range
		float
			Minimum antenna based SNR for self-calibration
		float
			Time interval for calibration 
		float
			Frequency interval of calibration
		float
			uvrange for calibration
		'''
		if self.quality_factor==0:
			residual_frac	=	0.03
			if self.safety_factor==0:
				start_sigma		=	9.0
				sigma_step		=	1.0
				min_sigma		=	6
				gain_minsnr		=	3.0		
				DR_delta_rms	=	25
				DR_delta_neg	=	20
				min_DR			=	20
				max_DR			=	100
				min_selfcal_snr	=	2.5
				skip_time		=	480
				skip_freq		=	2560
			elif self.safety_factor==1:
				start_sigma		=	9.0
				sigma_step		=	1.0
				min_sigma		=	7
				gain_minsnr		=	3.0		
				DR_delta_rms	=	22
				DR_delta_neg	=	18
				min_DR			=	25
				max_DR			=	500
				min_selfcal_snr	=	2.5
				skip_time		=	240
				skip_freq		=	2560
			else:
				start_sigma		=	9.0
				sigma_step		=	1.0
				min_sigma		=	8
				gain_minsnr		=	3.0		
				DR_delta_rms	=	20
				DR_delta_neg	=	15
				min_DR			=	30
				max_DR			=	1000
				min_selfcal_snr	=	2.5
				skip_time		=	120
				skip_freq		=	2560
		elif self.quality_factor==1:
			residual_frac	=	0.015
			if self.safety_factor==0:
				start_sigma		=	10.0
				sigma_step		=	0.5
				min_sigma		=	7
				gain_minsnr		=	4		
				DR_delta_rms	=	20
				DR_delta_neg	=	15
				min_DR			=	30
				max_DR			=	1000
				min_selfcal_snr	=	3
				skip_time		=	120
				skip_freq		=	1280
			elif self.safety_factor==1:
				start_sigma		=	10.0
				sigma_step		=	0.5
				min_sigma		=	8
				gain_minsnr		=	4		
				DR_delta_rms	=	18
				DR_delta_neg	=	12
				min_DR			=	35
				max_DR			=	5000
				min_selfcal_snr	=	3
				skip_time		=	60
				skip_freq		=	1280
			else:
				start_sigma		=	10.0
				sigma_step		=	0.5
				min_sigma		=	9
				gain_minsnr		=	4		
				DR_delta_rms	=	15
				DR_delta_neg	=	10
				min_DR			=	40
				max_DR			=	10000
				min_selfcal_snr	=	3
				skip_time		=	30
				skip_freq		=	1280
		else:
			residual_frac	=	0.01
			if self.safety_factor==0:
				start_sigma		=	11.0
				sigma_step		=	0.25
				min_sigma		=	8
				gain_minsnr		=	4.5		
				DR_delta_rms	=	18
				DR_delta_neg	=	12
				min_DR			=	40
				max_DR			=	10000
				min_selfcal_snr	=	3
				skip_time		=	15
				skip_freq		=	640
			elif self.safety_factor==1:
				start_sigma		=	11.0
				sigma_step		=	0.25
				min_sigma		=	9
				gain_minsnr		=	4.5		
				DR_delta_rms	=	15
				DR_delta_neg	=	10
				min_DR			=	45
				max_DR			=	50000
				min_selfcal_snr	=	3
				skip_time		=	12
				skip_freq		=	640
			else:
				start_sigma		=	11.0
				sigma_step		=	0.25
				min_sigma		=	10
				gain_minsnr		=	4.5		
				DR_delta_rms	=	12
				DR_delta_neg	=	8
				min_DR			=	50
				max_DR			=	100000
				min_selfcal_snr	=	3
				skip_time		=	10
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

	Parameters
	----------
	msname : str 
		Name of the measurement set
	Returns
	-------
	float
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
	Calculate flaged fraction from caltable
	
	Parameters
	----------
	caltable : str 
		Name of the CASA caltable	
	Returns
	-------
	float
		Fraction of the total solutions are flagged
	'''
	#TODO: As of now only based on CASA caltable, later CALIBRATE caltables also included
	
	tb=table()
	tb.open(caltable)
	flag=tb.getcol('FLAG')
	tb.close()
	flagged_fraction=np.sum(flag)/float(flag.size)
	return flagged_fraction

def calc_flag_chans_caltable(caltable,flag_frac=1.0):
	'''
	Calculate flaged channels from caltable

	Parameters
	----------
	caltable : str 
		Name of the CASA caltable
	flag_frac : float 
		Minimum fraction of data flagged for a single channel	
	Returns
	-------
	list
		Flagged channels list
	float
		flag fraction
	'''
	#TODO: As of now only based on CASA caltable, later CALIBRATE caltables also included
	
	tb=table()
	tb.open(caltable)
	flag=tb.getcol('FLAG')
	tb.close()
	y=(flag==True)
	flag_ants=np.nansum(np.nansum(y,axis=-1),axis=0).astype('float')
	flag1=flag+True
	x=(flag1==True)
	x_tot=np.nansum(np.nansum(x,axis=-1),axis=0).astype('float')
	frac=flag_ants/x_tot
	flagged_chans=np.where(frac>=flag_frac)[0].tolist()
	return flagged_chans,frac
##############################################
# General usuage functions #
##############################################

def getnearpos(array,value):
	'''
	Function to return index of two elements nearest to a given number
	
	Parameters
	----------
	array : numpy.array 
		Input numpy array
	value : float 
		Value to which nearest element search
	Returns
	-------
	float
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
	Function to transform error code to error message
		
	Parameters
	----------
	err_code : int
		Error code
	Returns
	-------
	str
		Error message
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

	Parameters
	----------
	radeg : float 
		Value of the RA in degree
	decdeg : float 
		Value of the DEC in degree
	Returns
	-------
	numpy.array
		['RA','DEC'] in hh mm ss dd mm ss format
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

	Parameters
	----------
	mjdsec : float 
		CASA MJD seconds
	includedate : bool 
		Include date in timestamp
	format : int 
		Datetime string format 
			0: 'YYYY/MM/DD/hh:mm:ss'

			1: 'YYYY-MM-DDThh:mm:ss'

			2: 'YYYY-MM-DD hh:mm:ss' 

			3: 'YYYY_MM_DD_hh_mm_ss'
	Returns
	-------
	str
		CASA time stamp in UTC at the format
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
	Convert timestamp to mjd second

	Parameters
	----------
	timestamp : str 
		Time stamp to convert 
	format : int 
		Datetime string format 
			0: 'YYYY/MM/DD/hh:mm:ss'
 
			1: 'YYYY-MM-DDThh:mm:ss'

			2: 'YYYY-MM-DD hh:mm:ss' 

			3: 'YYYY_MM_DD_hh_mm_ss'	
	Returns
	-------
	float
		Return correspondong MJD second of the day
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
	elif format==3:
		try:
			timestamp_datetime=datetime.strptime(timestamp,'%Y_%m_%d_%H_%M_%S.%f')
		except:
			timestamp_datetime=datetime.strptime(timestamp,'%Y_%m_%d_%H_%M_%S')
	mjd=float("{:.2f}".format((julian.to_jd(timestamp_datetime)-2400000.5)*(24.*3600.)))
	return mjd

def radec_to_altaz(ra,dec,obstime,LAT,LON,ALT):
	'''
	Function to convert radec to altaz for a given Earth location

	Parameters
	----------
	ra : str 
		RA either in degree or 'hh:mm:ss' or '%fh%fm%fs' format
	dec : str 
		DEC either in degree or 'dd:mm:ss' or '%fd%fm%fs'format
	obstime : str 
		Time of the observation in 'yyyy-mm-dd hh:mm:ss' format
	LAT : float
		Latitude of the Earth location in degree
	LON : float 
		Longitude of Earth location in degree 
	ALT : float 
		Altitude of the Earth location in meter
	Returns
	-------
	float
		Elevation
	float
		Azimuth in degree
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

	Parameters
	----------
	alt : float 
		Elevation in degree
	az : float 
		Azimuth in degree
	obstime : str 
		Time of the observation in 'yyyy-mm-dd hh:mm:ss' format
	LAT : float 
		Latitude of the Earth location in degree
	LON : float 
		Longitude of Earth location in degree 
	ALT : float 
		Altitude of the Earth location in meter
	Returns
	-------
	numpy.array
		['RA','DEC'] in hh mm ss dd mm ss format
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

	Parameters
	----------
	freq : float 
		Frequency in MHz
	Returns
	-------
	int
		MWA coarse channel number
	'''
	freq=float(freq)
	coarse_chans=[[(i*1.28)-0.64,(i*1.28)+0.64] for i in range(300)]
	for i in range(len(coarse_chans)):
		ch0=round(coarse_chans[i][0],2)
		ch1=round(coarse_chans[i][1],2)
		if freq>=ch0 and freq<ch1:
			return i 

def update_mwa_obsids(obsid_file='',verbose=False):
	'''
	Function to update MWA OBSIDs 

	Parameters
	----------
	obsid_file : str 
		Name of the file to save MWA OBSIDs
	verbose : bool
		Verbose output
	Returns
	-------
	str
		OBSID file name
	int
		Update success or failure code (0 or 1)
	'''
	if verbose==True:
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
		if verbose==True:
			print ('Last OBSID in MWA metadata server : '+str(end_obsid)+'\n')
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
		if verbose==True:
			print ('Updated successfully.\n')
		os.system('rm -rf casa*log')
		return obsid_file+'.npy',0
	except Exception as e:
		if verbose==True:
			print ('Error in update : '+str(e)+'\n')
			print ('Update not successful.\n')
		os.system('rm -rf casa*log')
		return obsid_file+'.npy',1

def get_OBSID_from_ms(msname):
	'''
	Function to return OBSID of an MWA observation

	Parameters
	----------
	msname : str 
		Name of the measurement set
	Returns
	-------
	int
		MWA Observation ID
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
		os.system('rm -rf casa*log')
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
				os.system('rm -rf casa*log')
				return OBSid
			except Exception as e:
				print ('Error occured : '+str(e)+'\n')
				c+=1
				time.sleep(5.0)
				if c>=5:
					os.system('rm -rf casa*log')
					return 0

def get_OBSID_from_metafits(metafits):
	'''
	Function to return OBSID of an MWA observation

	Parameters
	----------
	metafits : str 
		Name of the metafits file
	Returns
	-------
	int
		MWA Observation ID
	'''
	OBSid=fits.getheader(metafits)['GPSTIME']
	return OBSid 

def download_metafits(msname,outdir):
	'''
	Function to download MWA metafits of a given measurement set (Require internet connection)

	Parameters
	----------
	msname : str 
		Name of the measurement set
	outdir : str
		Name of the outputdir
	Returns
	-------
	str
		Metafits file name
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

	Parameters
	----------
	filelist : list 
		List of numpy table files
	outputfile : str 
		Output compressed file name
	Returns
	-------
	str
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

def multifreq_gaincal_interpolate(gaintables=[],output_gaintable='',outputfreq=0):
	'''
	Function to calculate linearly interpolated gain phase from a set of gaintables (at-least two) at multiple frequencies

	Parameters
	----------
	gaintables : list 
		List of gaintables at multiple frequencies
	output_gaintable : str 
		Name of the output gaintable 
	outputfreq : float
		Output frequency in MHz
	Returns
	-------
	str
		Output gaintable name
	'''
	if outputfreq==0:
		print ('Output frequency is not given.\n')
		return
	elif len(gaintables)<2:
		print ('Minimum two gaintables at different frequencies are required.\n')
		return
	else:
		if output_gaintable=='':
			output_gaintable='Gaintable_interp_'+str(outputfreq)+'MHz.gcal'
		if output_gaintable[-1]=='/':
			output_gaintable=output_gaintable[:-1]
		tb=table()
		freqs=[]
		gains=[]
		for g in gaintables:
			tb.open(g+'/SPECTRAL_WINDOW')
			freqs.append(tb.getcol('REF_FREQUENCY')[0])		
			tb.close()
			tb.open(g)
			gains.append(tb.getcol('CPARAM'))
			tb.close()
		gains=np.array(gains)
		if os.path.exists(output_gaintable)==True:
			print ('Removing previous existing caltable :'+output_gaintable+'\n')
			os.system('rm -rf '+output_gaintable)
		os.system('cp -r '+gaintables[0]+' '+output_gaintable)
		tb.open(output_gaintable)
		interpolated_gains=tb.getcol('CPARAM')
		tb.close()
		for i in range(interpolated_gains.shape[0]):
			for j in range(interpolated_gains.shape[-1]):
				x=np.polyfit(freqs,np.angle(gains[:,i,0,j]),deg=1)
				pq=np.poly1d(x,)
				phase=pq(outputfreq*10**6)
				interpolated_gains[i,0,j]=np.cos(phase)+1j*np.sin(phase)
		tb.open(output_gaintable+'/SPECTRAL_WINDOW',nomodify=False)
		chan_freq=tb.getcol('CHAN_FREQ')
		chan_freq=chan_freq+(outputfreq*10**6-chan_freq[0])
		tb.putcol('CHAN_FREQ',chan_freq)
		ref_freq=tb.getcol('REF_FREQUENCY')
		ref_freq[0]=outputfreq*10**6
		tb.putcol('REF_FREQUENCY',ref_freq)
		tb.flush()
		tb.close()
		tb.open(output_gaintable,nomodify=False)
		tb.putcol('CPARAM',interpolated_gains)
		tb.putkeyword('CALTYPE','p')
		tb.flush()
		tb.close()
		return output_gaintable

def get_caltable_metadata(caltable):
	'''
	Function to get caltable metadata

	Parameters
	----------
	caltable : str 
		Name of the caltable
	Returns
	-------
	dict
		A python dictionary with keywords MSNAME, JonesType, Channel 0 frequency (MHz), Central channel frequency (MHz), Channel width (kHz), Bandwidth (MHz), Start time, End time
	'''
	tb=table()
	tb.open(caltable)
	caltype=tb.getkeywords()['VisCal']
	msname=tb.getkeywords()['MSName']
	tb.close()
	tb.open(caltable+'/SPECTRAL_WINDOW')
	ch0=(tb.getcol('REF_FREQUENCY')[0])/10**6 # In MHz
	chanwidth=(tb.getcol('CHAN_WIDTH')[0]/10**3)[0] # In kHz
	freqlist=tb.getcol('CHAN_FREQ')
	chm=(freqlist[int(len(freqlist)/2)]/10**6)[0] # In MHz
	bw=(tb.getcol('TOTAL_BANDWIDTH')[0]/10**6) # In MHz
	tb.close()
	tb.open(caltable+'/OBSERVATION')
	timerange= tb.getcol('TIME_RANGE')
	start_time=mjdsec_to_timestamp(np.min(timerange),includedate=True,format=0)
	end_time=mjdsec_to_timestamp(np.max(timerange),includedate=True,format=0)
	tb.close()
	result={'MSNAME':msname,'JonesType':caltype,'Channel 0 frequency (MHz)':ch0,'Central channel frequency (MHz)':chm,'Channel width (kHz)':chanwidth,'Bandwidth (MHz)':bw,\
			'Start time':start_time,'End time':end_time}
	return result

def get_MWA_phase(metafits):
	'''
	Function to get MWA phase

	Parameters
	----------
	metafits : str 
		Name of the metafits file
	Returns
	-------
	str
		MWA phase
	'''
	OBSID=get_OBSID_from_metafits(metafits)
	try:
		url='http://ws.mwatelescope.org/metadata/con?obs_id='+str(OBSID)+'&summary'
		config=json.load(urllib.request.urlopen(url,timeout=15))[0]
		if config=='PHASE1':
			config='MWAPhaseI'
		elif config=='LB':
			config='MWAPhaseIILB'
		elif config=='COMPACT':
			config='MWAPhaseIICOMPACT'
	except:
		tilename=fits.getdata(metafits)['TileName']
		LB=['LB' in i for i in tilename]
		Hex=['Hex' in i for i in tilename]
		if np.sum(LB)!=0:
			config='MWAPhaseIILB'
		elif np.sum(Hex)!=0:
			config='MWAPhaseIICOMPACT'
		else:
			config='MWAPhaseI'
	return config

def paircars_instance_runner(cmd,basedir,paircars_basedir,screen_name,finished_touch_file,jobid,prefix_cmds=[]):
	'''
	Function to run a P-AIRCARS instance

	Parameters
	----------
	cmd : str 
		Command to run
	basedir : str
		Base directory of the measurement set
	paircars_basedir : str
		Base directory for a particular P-AIRCARS job
	screen_name : str 
		Name of the screen
	finished_touch_file : str
		Hidden file to generate at the base directory after finishing the instance
	'''
	batch_file=basedir+'/'+screen_name+'.batch'
	cmd_batch=basedir+'/'+screen_name+'_cmd.batch'
	if os.path.isdir(paircars_basedir+'/Logs_and_Errors/Logs')==False:
		os.makedirs(paircars_basedir+'/Logs_and_Errors/Logs')
	if os.path.isdir(paircars_basedir+'/Logs_and_Errors/Errors')==False:
		os.makedirs(paircars_basedir+'/Logs_and_Errors/Errors')
	outputfile=paircars_basedir+'/Logs_and_Errors/Logs/'+screen_name+'.log'
	output_error=paircars_basedir+'/Logs_and_Errors/Errors/'+screen_name+'.error'
	pid_file=paircars_basedir+'/'+str(jobid)+'_pids.log'
	cmd+=';sleep 2 ;if ! ls '+finished_touch_file+'_* ; then  touch '+finished_touch_file+'_error ;  fi'
	cmd='export PYTHONUNBUFFERED=1;echo \"'+cmd+'\" > '+cmd_batch+';sleep 2; chmod a+rwx '+cmd_batch+'; sleep 2; nohup sh '+cmd_batch+' > '+outputfile+' 2>'+output_error+\
		' < /dev/null &; echo $! >> '+pid_file+'; sleep 2'
	cmd=cmd.split(';')
	cmds=[i+'\n' for i in cmd]
	if os.path.exists(cmd_batch):
		os.system('rm -rf '+cmd_batch)
	if os.path.exists(batch_file):
		os.system('rm -rf '+batch_file)
	if os.path.isfile(batch_file):
		fil=open(batch_file,'r+')
	else:
		fil=open(batch_file,'w')
	if len(prefix_cmds)!=0:
		fil.writelines(prefix_cmds)
	fil.writelines(cmds)
	fil.close()
	os.system('chmod a+rwx '+batch_file)
	del cmd,prefix_cmds
	return basedir+'/'+screen_name+'.batch'

def get_available_paircars_instance(basedir,jobid,total_instances):
	'''
	Function to get available P-AIRCARS instances

	Parameters
	----------
	basedir : str
		Name of the P-AIRCARS base directory
	job_id : int
		P-AIRCARS job ID
	total_instances : int
		Total P-AIRCARS instances
	Returns
	-------
	int
		Available P-AIRCARS instance
	'''
	c=0
	while True:
		if os.path.isfile(basedir+'/'+str(jobid)+'_pids.log')==False:
			running_pids=[]
		else:
			pids=np.loadtxt(basedir+'/'+str(jobid)+'_pids.log',unpack=True).astype('int')
			if pids.shape==():
				pids=np.array([int(pids)])			
			running_pids=[psutil.pid_exists(pid) for pid in pids]
		available_paircars_instance=int(total_instances-np.sum(running_pids))
		total_cpus_required=int(available_paircars_instance*1.5)
		available_cpus=psutil.cpu_count()*(1-(psutil.cpu_percent(interval=10)/100.0))
		if total_cpus_required>available_cpus:
			c+=1
			if c>5:
				available_paircars_instance=0
				break
		else:
			break
	return available_paircars_instance

def get_final_psf(msname,imager='wsclean',weight='briggs',robust=0.5):
	'''
	Function to get a common final psf parameters

	Parameters
	----------
	msname : str
		Name of the measurement set
	imager : str
		Imager, 'wsclean' or 'CASA'
	weight : str
		Visibility weighting for imaging
	robust : float
		Robust parameter for 'briggs' weighting
	Returns
	-------
	float 
		Major axis (FWHM) in arcsec
	float
		Minor axis (FWHM) in arcsec
	float 
		Position angle in degree
	'''
	imagename=msname.split('.ms')[0]+'_test_image'
	IB=ImageBasic(msname,includeflag=True)
	cell=IB.calc_cellsize(3) # Assuming 3 pixels in one PSF
	imsize=int(IB.num_pixels(3)/4)
	if imsize>100:
		imsize=100
	scales=IB.choose_scales(3,32*60)
	if imager=='wsclean':
		if weight=='briggs':
			weight=weight+' '+str(robust)
		wsclean_args=['-scale '+str(cell)+'asec','-size '+str(imsize)+' '+str(imsize),'-no-dirty','-quiet','-j 3','-weight '+weight,'-name '+str(imagename),'-nwlayers 1','-niter 1']
		os.system('wsclean '+' '.join(wsclean_args)+' '+msname)
		header=fits.getheader(imagename+'-image.fits')
		maj=(header['BMAJ']*3600)
		minor=(header['BMIN']*3600)
		pa=header['BPA']
		os.system('rm -rf '+imagename+'*')
		return maj,minor,pa
	else:
		robust=robust*2
		tclean(vis=msname,imsize=[imsize],cell=str(cell)+'arcsec',niter=1,imagename=imagename,weighting=weight,robust=robust)
		header=imhead(imagename=imagename+'.image',mode='list')
		maj=header['beammajor']['value']
		minor=header['beamminor']['value']
		pa=header['beampa']['value']
		os.system('rm -rf '+imagename+'* casa*log')
		return maj,minor,pa

def check_internet():
	'''
	Function to check internet connection
	Returns
	-------
	bool
		Whether internet connection is available or not
	'''
	url='https://www.google.com'
	try_count=5
	timeout=5
	c=0
	while c<=try_count:
		try:
			requests.get(url,timeout=timeout)
			return True
		except:
			continue
		c+=1
	return False 

