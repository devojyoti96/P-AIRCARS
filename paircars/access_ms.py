import numpy as np,os,copy,glob
from datetime import datetime as dt, timedelta
from casatools import *
from casatasks import *
from astropy.time import Time
from astropy.coordinates import get_sun
from . import basic_func as B
'''
Code is written by Devojyoti Kansabanik, 06 Jan, 2021
'''

### Functions related to measurement set operations ###
class AccessMS:
	'''
	Generic class to perform measurement set related operation
	Attribute:
	msname = Name of the measurement set
	'''
	def __init__(self,msname):
		self.msname=msname
		self.md=msmetadata()
		self.tb=table()
		self.me=measures()

	def get_phasecenter(self):
		'''
		Get phasecenter of the measurement set
		Return:
		Phasecenter of the measurement set in ['RA','DEC'] in hh mm ss dd mm ss format, RA in degree, DEC in degree 
		'''
		self.md.open(self.msname)
		phasecenter=self.md.phasecenter()
		freq=self.md.meanfreq(0)
		self.md.close()
		radeg=np.rad2deg(phasecenter['m0']['value'])
		decdeg=np.rad2deg(phasecenter['m1']['value'])
		radec_str=B.radec_con_deg_to_hhmmss(radeg,decdeg)
		return radec_str,radeg,decdeg		

	def convert_mwa_to_iau(self):
		'''
		Convert the data from MWA coordinate (X=E->W and Y=N->S) to IAU coordinate (X=S->N,Y=W->E) system.
		Thus X=>-Y, and Y=>-X. So, XX=>YY,YY=>XX,XY=>YX,YX=>XY
		Returns:
		Measurement set in IAU convention and confirms the measurement set convention
		'''
		code=vishead(vis=self.msname,mode='get',hdkey='fld_code')[0][0]
		code_list=code.split(',')
		if 'IAU' not in code_list:
			self.tb.open(self.msname,nomodify=False)
			data=self.tb.getcol('DATA')
			data_iau=copy.deepcopy(data)
			data_iau[0,:,:]=data[3,:,:]
			data_iau[1,:,:]=data[2,:,:]
			data_iau[2,:,:]=data[1,:,:]
			data_iau[3,:,:]=data[0,:,:]
			self.tb.putcol('DATA',data_iau)
			self.tb.flush()
			self.tb.close()
			del data,data_iau
			if len(code_list)==1 and code_list[0]=='':
				code+='IAU'
			else:
				code+=',IAU'
			vishead(vis=self.msname,mode='put',hdkey='fld_code',hdvalue=np.array([code]))
			msg='Measurement set :'+self.msname+' has been converted in IAU convention'
			return msg
		else:
			msg='Measurement set :'+self.msname+' is already in IAU convention'
			return msg
		
	def calc_ncoarse(self):
		'''
		Calculate the number of MWA coarse channels in the measurement set
		Return:
		Number of coarse channel in the measurement set
		'''
		self.md.open(self.msname)
		ncoarse=int(self.md.bandwidths()[0]/(1.28*10**6))
		self.md.close()
		return ncoarse

	def calc_freqres(self):
		'''
		Return frequency resolution of the data
		Return:
		Frequency resolution in kHz
		'''
		self.md.open(self.msname)
		freqres=self.md.chanres(0)[0]/10**3
		self.md.close()
		return freqres
		
	def calc_bandwidth(self):
		'''
		Function to calculate bandwidth of the data
		Return:
		Bandwidth in Hz
		'''
		self.md.open(self.msname)
		bandwidth=self.md.bandwidths()[0]
		self.md.close()
		return bandwidth

	def calc_timeres(self):
		'''
		Function to calculate time resolution of the measurement set
		Return :
		Time resolution in second
		'''
		self.md.open(self.msname)
		times=self.md.timesforfield(0)
		if len(times)==1:
			timeres=self.md.exposuretime(1)['value']
		else:
			timeres=times[1]-times[0]
		self.md.close()
		return timeres

	def get_num_timestamps(self):
		'''
		Function to calculate number of timestamps
		Return:
		Number of timestamps
		'''
		self.md.open(self.msname)
		times=self.md.timesforfield(0)
		self.md.close()
		return len(times)

	def get_num_channels(self):
		'''
		Function to calculate number of channels
		Return:
		Number of channels
		'''
		self.md.open(self.msname)
		nchan=self.md.nchan(0)
		self.md.close()
		return nchan

	def get_timestamps_in_mjdsecs(self):
		'''
		Function to return list of timestamps
		Return:
		List of timestamps in MJD seconds
		'''
		self.md.open(self.msname)
		mjds=self.md.timesforfield(0)
		self.md.close()
		return mjds

	def get_timestamps(self):
		'''
		Function to return list of timestamps
		Return:
		List of timestamps in 'YYYY/MM/DD/hh:mm:ss' format
		'''
		timestamps=[B.mjdsec_to_timestamp(mjdsec,includedate=True,format=0) for mjdsec in self.get_timestamps_in_mjdsecs()]	
		return timestamps	
	
	def get_obs_startend_time(self):
		'''
		Function to get observation start time
		Return:
		Observation start time and end time in 'YYYY/MM/DD/hh:mm:ss' format, start_mjdsec, end_mjdsec
		'''
		self.md.open(self.msname)
		start_mjdsec=self.md.timerangeforobs(0)['begin']['m0']['value']*24*3600.0
		end_mjdsec=self.md.timerangeforobs(0)['end']['m0']['value']*24*3600.0
		self.md.close()
		start=B.mjdsec_to_timestamp(start_mjdsec,includedate=True,format=0)
		end=B.mjdsec_to_timestamp(end_mjdsec,includedate=True,format=0)
		return start,end,start_mjdsec,end_mjdsec

	def get_scan_startend_time(self):
		'''
		Function to get scan start time
		Return:
		Scan start time and end time in 'YYYY/MM/DD/hh:mm:ss' format, start_mjdsec, end_mjdsec
		'''
		self.md.open(self.msname)
		time_list=self.md.timesforfield(0)
		self.md.close()
		start_mjdsec=time_list[0]
		end_mjdsec=time_list[-1]
		start=B.mjdsec_to_timestamp(start_mjdsec,includedate=True,format=0)
		end=B.mjdsec_to_timestamp(end_mjdsec,includedate=True,format=0)
		return start,end,start_mjdsec,end_mjdsec

	def get_num_antenna(self):
		'''
		Function to get total number of antennas in the ms
		Return:
		Number of antennas
		'''
		self.md.open(self.msname)
		nant=self.md.nantennas()
		self.md.close()
		return nant

	def get_freqs(self):
		'''
		Function to return list of channel frequencies
		Return:
		List of frequencies in Hz
		'''
		self.md.open(self.msname)
		freqs=self.md.chanfreqs(0)
		self.md.close()
		return freqs

	def get_unflag_chan(self,flagfrac=1):
		'''
		Function to get the unflagged channels if flag fraction is less than certain value
		Parameter:
		flagfrac = Flag fraction per channel (default : 1)
		Return:
		List of unflagged channels
		'''
		self.tb.open(self.msname)
		flag=self.tb.getcol('FLAG')
		self.tb.close()
		nchan=self.get_num_channels()
		unflagged_chan=[]
		for chan in range(nchan):
			flagged_data=flag[:,chan,:]
			total_data=float(len(flagged_data.flatten()))
			flagged_data=float(np.sum(flagged_data.flatten()))
			if (flagged_data/total_data)<flagfrac:
				unflagged_chan.append(chan)
		return unflagged_chan	

	def get_model_nomodel_chan(self):
		'''
		Function to get the channels wioth no model data.
		Return:
		List of model and nomodel channels
		'''
		self.tb.open(self.msname)
		model=self.tb.getcol('MODEL_DATA')
		self.tb.close()
		nchan=self.get_num_channels()
		nomodel_chan=[]
		model_chan=[]
		for chan in range(nchan):
			model_data=model[:,chan,:]
			if np.abs(np.sum(model_data))==0.0 or np.abs(np.sum(model_data))==len(model_data.flatten())/2:
				nomodel_chan.append(chan)
			else:
				model_chan.append(chan)
		return model_chan,nomodel_chan

	def calc_meanfreq(self):
		'''
		Function to return central frequency of the measurement set (Only MS with single SPW)
		Return:
		Central frequency in Hz
		'''
		self.md.open(self.msname)
		cent_freq=self.md.meanfreq(0)
		self.md.close()
		return cent_freq

	def calc_meanwavelength(self):
		'''
		Function to return central wavelength of the measurement set in metre
		Return:
		Central wavelength in metre
		'''
		freq=self.calc_meanfreq()
		wavelength=299792458.0/freq
		return wavelength

	def radec_sun(self):
		'''
		RA DEC of the Sun at the middle of the measurement set
		Return:
		RA DEC of the Sun in J2000, RA in degree and DEC in degree
		'''
		self.md.open(self.msname)
		start_end_time=self.md.timerangeforobs(0)
		self.md.close()
		start_time_mjd=start_end_time['begin']['m0']['value']
		end_time_mjd=start_end_time['end']['m0']['value']
		mid_time_mjd=(start_time_mjd+end_time_mjd)/2.0
		mid_time_mjd=Time(str(mid_time_mjd),format='mjd')
		sun_coord=get_sun(mid_time_mjd)
		sun_ra_deg=sun_coord.ra.deg
		sun_dec_deg=sun_coord.dec.deg
		sun_radec=B.radec_con_deg_to_hhmmss(sun_ra_deg,sun_dec_deg) 
		sun_radec_string='J2000 '+str(sun_radec[0])+' '+str(sun_radec[1])
		return sun_radec_string,sun_ra_deg,sun_dec_deg
		
	def move_phasecenter_to_sun(self):
		'''
		Move the phasecenter of the measurement set at the center of the Sun.
		Return:
		Measurement set at phase center at the center of the Sun
		'''
		sun_radec_string,sunra_deg,sundec_deg=self.radec_sun()
		radec_str,ra,dec=self.get_phasecenter()
		IB=B.ImageBasic(self.msname)
		psfsize=IB.calc_psf()
		diff=(sunra_deg-ra)**2+(sundec_deg-dec)**2
		code=vishead(vis=self.msname,mode='get',hdkey='fld_code')[0][0]
		code_list=code.split(',')
		if 'FIXVIS' not in code_list:
			if len(code_list)==1 and code_list[0]=='':
				code+='FIXVIS'
			else:
				code+=',FIXVIS'
			vishead(vis=self.msname,mode='put',hdkey='fld_code',hdvalue=np.array([code]))
			if diff<(psfsize/3600.0)**2:
				return 'Phasecenter shift is less than PSF size. No shift is required.\n' 
			elif 'FIXVIS' not in code_list:
				fixvis(vis=self.msname,outputvis=self.msname,phasecenter=sun_radec_string)
				return 'Phasecenter of the observation is moved to Sun center at :'+sun_radec_string+'.\n'
		else:
			return 'Phasecenter is already shifted to the Sun.\n'

	def move_phasecenter_to_source(self,radec=''):
		'''
		Function to move the phasecenter of measurement set at particular RA-DEC
		Parameters:
		radec = RA DEC string, format \'J2000 00h00m00s +00d00m00s\'
		Return:
		Phasecenter of the measurement set is moved to the new phasecenter
		'''
		if radec=='':
			return 'No RA-DEC is given'
		else:
			code=vishead(vis=self.msname,mode='get',hdkey='fld_code')[0][0]
			code_list=code.split(',')
			radec_string='_'.join(radec.split(' ')[1:])
			if 'FIXVIS_'+radec_string not in code_list:
				for i in range(len(code_list)):
					if 'FIXVIS' in code_list[i]:
						code_list.remove(code_list[i])
				code=','.join(code_list)
				fixvis(vis=self.msname,outputvis=self.msname,phasecenter=radec)
				if len(code_list)==1 and code_list[0]=='':
					code+='FIXVIS_'+radec_string
				else:
					code+=',FIXVIS_'+radec_string
				vishead(vis=self.msname,mode='put',hdkey='fld_code',hdvalue=np.array([code]))
				return 'Move phasecenter to :'+radec+'.\n'
			else:
				return 'Phasecenter is already at :'+radec+'.\n'
		return
			
	def get_max_baseline(self):
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

	def get_listobs(self,listobsfile=''):
		'''
		Function to save listobs
		Return:
		listobs file name
		'''
		if listobsfile=='':
			listobsfile=self.msname.split('.ms')[0]+'.listobs'
		if os.path.isfile(listobsfile):
			os.system('rm -rf '+listobsfile)
		listobs(vis=self.msname,listfile=listobsfile)
		return listobsfile

	def make_antenna_list(self,num_bins=5,antenna_list_file=''):
		'''
		Make the antenna addition list
		Parameters:
		num_bin = Number of antenna bins
		antenna_list_file = Name of the file to save antenna list (Keep blank to avoid save)
		Return:
		Lists of antennas to be added in steps
		'''
		listobs_file=self.msname.split('.ms')[0]+'.listobs'
		if os.path.isfile(listobs_file)==False:
			listobs(vis=self.msname,listfile=listobs_file)
		else:
			pass
		fil=open(listobs_file,'r')
		lines=fil.readlines()
		start_line_num=0
		for i in range(len(lines)):
			if 'East' in lines[i] and 'North' in lines[i] and 'Elevation' in lines[i]:
				start_line_num=i
			if 'Observation:' in lines[i]:
				instrument=lines[i].split('Observation:')[-1].split('\n')[0].split(' ')[-1]
		coords=np.loadtxt(listobs_file,skiprows=start_line_num+1,usecols=(-6,-5,0))
		antenna_ids=coords[:,-1]
		dist=np.sqrt(coords[:,0]**2+coords[:,1]**2)	
		pos=np.argsort(dist)
		while num_bins>=1:
			val,bins=np.histogram(dist,num_bins)
			if val[0]<len(coords)*0.3:
				num_bins-=1				
				continue
			else:
				break
		antenna_list=[]
		start=0
		for elem in val:
			antenna_list.append(np.sort(antenna_ids[0:start+elem].astype('int')+1))
			start+=elem
		if antenna_list_file!='':
			antenna_list_file=antenna_list_file.split('.')[0]
			np.save(antenna_list_file,np.array(['Antennafile',len(coords),antenna_list]))
		return antenna_list,len(coords)
	
	def get_observatory_loc(self):
		'''
		Give the observatory geodetic location
		Return:
		Latitude, Longitude in degree and Altitude in meter
		'''
		self.tb.open(self.msname+'::ANTENNA')
		pos = self.tb.getcol('POSITION')
		meanpos = np.mean(pos, axis=1)
		frame = self.tb.getcolkeyword('POSITION','MEASINFO')['Ref']
		units = self.tb.getcolkeyword('POSITION','QuantumUnits')
		mpos  = self.me.position(frame,str(meanpos[0])+units[0],str(meanpos[1])+units[1],str(meanpos[2])+units[2])
		self.me.doframe(mpos)
		self.tb.close()
		loc=self.me.measure(mpos,'WGS84')
		LAT=np.rad2deg(loc['m1']['value'])
		LON=np.rad2deg(loc['m0']['value'])
		ALT=loc['m2']['value']/2
		return LAT,LON,ALT

	def get_altaz(self,source_field=0,source_scan=1):
		'''
		Calculate alt az of the phasecenter
		Parameters:
		source_field = FIELD ID of the source
		source_scan = Scan number 
		Return:
		Alt,Az in radian
		''' 
		LAT,LON,ALT=self.get_observatory_loc()
		radec_str,radeg,decdeg=self.get_phasecenter()
		self.md.open(self.msname)
		scans=self.md.scansforfield(source_field)
		if source_scan not in scans:
			print ('Scan is not in given field ID')
			self.md.close()
			return
		else:
			mjd=np.mean(self.md.timesforscan(source_scan))
			time_string=B.mjdsec_to_timestamp(mjd,includedate=True,format=2)
			alt,az=B.radec_to_altaz(radeg,decdeg,time_string,LAT,LON,ALT)
			alt=np.deg2rad(alt)
			az=np.deg2rad(az)
			self.md.close()
			return alt,az

	def get_phasecenter_parang(self,source_field=0,combine='field'):
		'''
		Calculate parallactic for phasecenter at a given Earth location
		Parameters:
		source_file = FIELD id of the source
		combine = 'field' or 'scan' or '' for no combine
		Note : Parallactic angle is defined as the orientation of the sky in telescope coordinate. All angles are defined positive in IAU defined sky coordiunate (North to East).
			   So, the rotation of the sky with respect to telescope is negative in the sky coordinate. To account this effect parallactic angle is given in 360-parang form.
		Return:
		1. combine = 'field', A single parallactic angle for the entire field in degree
		2. combine = 'scan', A python dictionary {'scan':parang} format
		3. combine = '', A list of parang for all timestamps
		'''
		LAT,LON,ALT=self.get_observatory_loc()
		radec_str,radeg,decdeg=self.get_phasecenter()
		self.md.open(self.msname)
		if combine=='field':
			mjd=np.mean(self.md.timesforfield(source_field))
			time_string=B.mjdsec_to_timestamp(mjd,includedate=True,format=2)
			alt,az=B.radec_to_altaz(radeg,decdeg,time_string,LAT,LON,ALT)
			alt=np.deg2rad(alt)
			az=np.deg2rad(az)
			lat=np.deg2rad(LAT)
			p=-np.arctan2(np.sin(az)*np.cos(lat),np.cos(alt)*np.sin(lat) - np.sin(alt)*np.cos(lat)*np.cos(az))
			self.md.close()
			return 360-np.rad2deg(p)
		elif combine=='scan':
			scans=self.md.scansforfield(source_field)
			mjds=[np.mean(self.md.timesforscan(source_scan)) for source_scan in scans]
			time_strings=[B.mjdsec_to_timestamp(mjdsec,includedate=True,format=2) for mjdsec in mjds]
			parang={}		
			for i in range(len(scans)):
				obstime=time_strings[i]
				scan=scans[i]
				alt,az=B.radec_to_altaz(radeg,decdeg,obstime,LAT,LON,ALT)
				alt=np.deg2rad(alt)
				az=np.deg2rad(az)
				lat=np.deg2rad(LAT)
				p=-np.arctan2(np.sin(az)*np.cos(lat),np.cos(alt)*np.sin(lat) - np.sin(alt)*np.cos(lat)*np.cos(az))
				parang[scan]=360-np.rad2deg(p)
				self.md.close()
			return parang 
		else:
			mjds=self.md.timesforfield(source_field) 
			time_strings=[B.mjdsec_to_timestamp(mjdsec,includedate=True,format=2) for mjdsec in mjds]
			parang=[]		
			for i in range(len(time_strings)):
				obstime=time_strings[i]
				alt,az=B.radec_to_altaz(radeg,decdeg,obstime,LAT,LON,ALT)
				alt=np.deg2rad(alt)
				az=np.deg2rad(az)
				lat=np.deg2rad(LAT)
				p=-np.arctan2(np.sin(az)*np.cos(lat),np.cos(alt)*np.sin(lat) - np.sin(alt)*np.cos(lat)*np.cos(az))
				parang.append(360-np.rad2deg(p))
				self.md.close()
			return parang 

	def get_parang(self,ra,dec):# TODO : Finish this function
		'''
		Calculate parallactic for phasecenter at a given Earth location
		Parameters:
		source_file = FIELD id of the source
		combine = 'field' or 'scan' or '' for no combine
		Note : Parallactic angle is defined as the orientation of the sky in telescope coordinate. All angles are defined positive in IAU defined sky coordiunate (North to East).
			   So, the rotation of the sky with respect to telescope is negative in the sky coordinate. To account this effect parallactic angle is given in 360-parang form.
		Return:
		1. combine = 'field', A single parallactic angle for the entire field in degree
		2. combine = 'scan', A python dictionary {'scan':parang} format
		3. combine = '', A list of parang for all timestamps
		'''
		LAT,LON,ALT=self.get_observatory_loc()
		radec_str,radeg,decdeg=self.get_phasecenter()
		self.md.open(self.msname)
		if combine=='field':
			mjd=np.mean(self.md.timesforfield(source_field))
			time_string=B.mjdsec_to_timestamp(mjd,includedate=True,format=2)
			alt,az=B.radec_to_altaz(radeg,decdeg,time_string,LAT,LON,ALT)
			alt=np.deg2rad(alt)
			az=np.deg2rad(az)
			lat=np.deg2rad(LAT)
			p=-np.arctan2(np.sin(az)*np.cos(lat),np.cos(alt)*np.sin(lat) - np.sin(alt)*np.cos(lat)*np.cos(az))
			self.md.close()
			return 360-np.rad2deg(p)
		elif combine=='scan':
			scans=self.md.scansforfield(source_field)
			mjds=[np.mean(self.md.timesforscan(source_scan)) for source_scan in scans]
			time_strings=[B.mjdsec_to_timestamp(mjdsec,includedate=True,format=2) for mjdsec in mjds]
			parang={}		
			for i in range(len(scans)):
				obstime=time_strings[i]
				scan=scans[i]
				alt,az=B.radec_to_altaz(radeg,decdeg,obstime,LAT,LON,ALT)
				alt=np.deg2rad(alt)
				az=np.deg2rad(az)
				lat=np.deg2rad(LAT)
				p=-np.arctan2(np.sin(az)*np.cos(lat),np.cos(alt)*np.sin(lat) - np.sin(alt)*np.cos(lat)*np.cos(az))
				parang[scan]=360-np.rad2deg(p)
				self.md.close()
			return parang 
		else:
			mjds=self.md.timesforfield(source_field) 
			time_strings=[B.mjdsec_to_timestamp(mjdsec,includedate=True,format=2) for mjdsec in mjds]
			parang=[]		
			for i in range(len(time_strings)):
				obstime=time_strings[i]
				alt,az=B.radec_to_altaz(radeg,decdeg,obstime,LAT,LON,ALT)
				alt=np.deg2rad(alt)
				az=np.deg2rad(az)
				lat=np.deg2rad(LAT)
				p=-np.arctan2(np.sin(az)*np.cos(lat),np.cos(alt)*np.sin(lat) - np.sin(alt)*np.cos(lat)*np.cos(az))				
				parang.append(360-np.rad2deg(p))
				self.md.close()
			return parang 

def splited_ms_rename(msname,ref_time_chan=True,change_msname=True):
	'''
	Convert the name of a splited measurement set according to its frequency and timestamp
	Parameters:
	msname = Name of the measurement set
	ref_time_chan = True, whether the time and frequency slice is refernce
	change_msname = True, change the msname or just return the name
	Return:
	Chnage the name of the measurement set at return its changed name
	'''
	md=msmetadata()
	md.open(msname)
	mjdsecond=md.timesforscan(1)[0]
	mean_freq=md.meanfreq(0)/10**6 # In MHz
	utcstring=B.mjdsec_to_timestamp(mjdsecond,includedate=True,format=0)
	yyyy=utcstring.split('/')[0]
	mm=utcstring.split('/')[1]
	dd=utcstring.split('/')[2]
	timestamp='_'.join(utcstring.split('/')[3].split(':')) 
	ms_path=os.path.dirname(os.path.realpath(msname))
	if ref_time_chan==False:
		filename=ms_path+'/time_'+yyyy+'_'+mm+'_'+dd+'_'+timestamp+'_freq_'+str(mean_freq)+'.ms'
	else:
		filename=ms_path+'/time_'+yyyy+'_'+mm+'_'+dd+'_'+timestamp+'_freq_'+str(mean_freq)+'_ref.ms'
	if change_msname==True and msname!=filename:
		if os.path.exists(filename):
			os.system('rm -rf '+filename)
		os.system('mv '+msname+' ' +filename)
	return filename










