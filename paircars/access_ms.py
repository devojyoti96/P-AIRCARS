import os
import numpy as np,copy,glob,julian
from datetime import datetime as dt, timedelta
from casatools import msmetadata,table,measures,quanta,agentflagger,image,calibrater
from casatasks import *
from astropy.time import Time
from astropy.coordinates import get_sun
from . import basic_func as B
os.system('rm -rf casa*log')
'''
Code is written by Devojyoti Kansabanik, 06 Jan, 2021
'''

### Functions related to measurement set operations ###
class AccessMS:
	'''
	Generic class to perform measurement set related operation
	
	Parameters
	----------
	msname : str 
		Name of the measurement set
	'''
	def __init__(self,msname):
		if msname[-1]=='/':
			msname=msname[:-1]
		self.msname=msname
		self.md=msmetadata()
		self.tb=table()
		self.me=measures()

	def get_nbaseline(self,autocor=True):
		'''		
		Function to get number of baselines

		Parameters
		----------
		autocor : bool 
			Include auto-correlations into account or not
		Returns
		-------
		int
			Number of baselines
		'''
		self.md.open(self.msname)
		nbs=self.md.nbaselines(ac=autocor)
		self.md.close()
		return nbs

	def get_phasecenter(self):
		'''
		Get phasecenter of the measurement set

		Returns
		-------
		str
			Phasecenter of the measurement set in ['RA','DEC'] in hh mm ss dd mm ss format
		float 
			RA in degree
		float 
			DEC in degree 
		'''
		self.md.open(self.msname)
		phasecenter=self.md.phasecenter()
		freq=self.md.meanfreq(0)
		self.md.close()
		radeg=np.rad2deg(phasecenter['m0']['value'])
		decdeg=np.rad2deg(phasecenter['m1']['value'])
		radec_str=B.radec_con_deg_to_hhmmss(radeg,decdeg)
		return radec_str,radeg,decdeg

	def get_pol(self):
		'''
		Function to get number of polarisations 
		Returns
		-------
		int
			Number of correlation products
		'''
		self.md.open(self.msname)
		npol=self.md.ncorrforpol()[0]
		self.md.close()
		return npol		

	def convert_mwa_to_iau(self):
		'''
		Convert the data from MWA coordinate (X=E->W and Y=N->S) to IAU coordinate (X=S->N,Y=W->E) system.
		Thus X=>-Y, and Y=>-X. So, XX=>YY,YY=>XX,XY=>YX,YX=>XY

		Returns
		------
		str
 			Confirmation message of coordinate conversion the measurement set
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

		Returns
		-------
		int
			Number of coarse channel in the measurement set
		'''
		self.md.open(self.msname)
		ncoarse=int(self.md.bandwidths()[0]/(1.28*10**6))
		self.md.close()
		return ncoarse

	def calc_freqres(self):
		'''
		Return frequency resolution of the data
	
		Returns
		-------
		float
			Frequency resolution in kHz
		'''
		self.md.open(self.msname)
		freqres=self.md.chanres(0)[0]/10**3
		self.md.close()
		return freqres
		
	def calc_bandwidth(self):
		'''
		Function to calculate bandwidth of the data
		
		Returns
		-------
		float
			Bandwidth in Hz
		'''
		self.md.open(self.msname)
		bandwidth=self.md.bandwidths()[0]
		self.md.close()
		return bandwidth

	def calc_total_time(self):
		'''
		Function to calculate total time of the measurement set

		Returns
		-------
		float
			Total time in seceond
		'''
		self.md.open(self.msname)
		times=self.md.timesforfield(0)
		if len(times)==1:
			timewidth=self.md.exposuretime(1)['value']
		else:
			timewidth=(times[-1]-times[0])+self.calc_timeres()
		self.md.close()
		return timewidth

	def calc_timeres(self):
		'''
		Function to calculate time resolution of the measurement set

		Returns
		-------
		float
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

		Returns
		-------
		int
			Number of timestamps
		'''
		return len(self.get_timestamps_in_mjdsecs()[0])

	def get_num_channels(self):
		'''
		Function to calculate number of channels

		Returns
		-------
		int
			Number of channels
		'''
		self.md.open(self.msname)
		nchan=self.md.nchan(0)
		self.md.close()
		return nchan

	def get_timestamps_in_mjdsecs(self):
		'''
		Function to return list of timestamps
	
		Returns
		-------
		list
			List of timestamps in MJD seconds
		list
			List of timestamps with errors
		'''
		self.md.open(self.msname)
		mjds=self.md.timesforfield(0)
		time_diff=np.ediff1d(mjds)
		pos=np.where(time_diff<self.calc_timeres()/2.0)
		wrong_timestamps=[]
		if len(pos[0])!=0:
			for p in pos[0]:
				try:
					wrong_timestamps.append(B.mjdsec_to_timestamp(mjds[p],includedate=True,format=0)+'~'+B.mjdsec_to_timestamp(mjds[p+1],includedate=True,format=0))
				except:
					wrong_timestamps.append(B.mjdsec_to_timestamp(mjds[p],includedate=True,format=0))			
		mjds=np.delete(mjds,pos)
		self.md.close()
		return mjds,wrong_timestamps

	def get_timestamps(self,format=0):
		'''
		Function to return list of timestamps

		Parameters
		----------
		format : int 
			Datetime string format 
				(0 : 'YYYY/MM/DD/hh:mm:ss', 1: 'YYYY-MM-DDThh:mm:ss', 2: 'YYYY-MM-DD hh:mm:ss', 3: 'YYYY_MM_DD_hh_mm_ss') 
		Returns
		-------
		list
			List of timestamps in UTC at the format 'YYYY/MM/DD/hh:mm:ss'
		'''
		timestamps=[B.mjdsec_to_timestamp(mjdsec,includedate=True,format=format) for mjdsec in self.get_timestamps_in_mjdsecs()[0]]	
		return timestamps	
	
	def get_obs_startend_time(self):
		'''
		Function to get observation start time

		Returns
		-------
		str
			Observation start time in 'YYYY/MM/DD/hh:mm:ss' format
		str 
			Observation end time in 'YYYY/MM/DD/hh:mm:ss' format
		float
			Observation start time in MJD second
		float
			Observation end time in MJD second
		'''
		self.md.open(self.msname)
		start_mjdsec=self.md.timerangeforobs(0)['begin']['m0']['value']*24*3600.0
		end_mjdsec=self.md.timerangeforobs(0)['end']['m0']['value']*24*3600.0
		self.md.close()
		start=B.mjdsec_to_timestamp(start_mjdsec,includedate=True,format=0)
		end=B.mjdsec_to_timestamp(end_mjdsec,includedate=True,format=0)
		return start,end,start_mjdsec,end_mjdsec

	def get_obs_date(self,format=0):
		'''
		Function to get observation date

		Parameters
		----------
		format : str
			Date string format 
				0: 'YYYY/MM/DD'

				1: 'YYYY-MM-DD'
 
				2: 'YYYY_MM_DD'
		Returns
		-------
		str
			Observation start date
		str
			Observation end date
		'''
		start,end,start_mjdsec,end_mjdsec=self.get_obs_startend_time()
		start_mjd=start_mjdsec/(24*3600.0)
		end_mjd=end_mjdsec/(24*3600.0)
		start_year=julian.from_jd(start_mjd+2400000.5,fmt='jd').date().year
		start_month=julian.from_jd(start_mjd+2400000.5,fmt='jd').date().month
		start_day=julian.from_jd(start_mjd+2400000.5,fmt='jd').date().day
		end_year=julian.from_jd(end_mjd+2400000.5,fmt='jd').date().year
		end_month=julian.from_jd(end_mjd+2400000.5,fmt='jd').date().month
		end_day=julian.from_jd(end_mjd+2400000.5,fmt='jd').date().day
		if format==0:
			start_date=str(start_year)+'/'+str(start_month)+'/'+str(start_day)
			end_date=str(end_year)+'/'+str(end_month)+'/'+str(end_day)
		elif format==1:
			start_date=str(start_year)+'-'+str(start_month)+'-'+str(start_day)
			end_date=str(end_year)+'-'+str(end_month)+'-'+str(end_day)
		else:
			start_date=str(start_year)+'_'+str(start_month)+'_'+str(start_day)
			end_date=str(end_year)+'_'+str(end_month)+'_'+str(end_day)
		return start_date,end_date

	def get_scan_startend_time(self):
		'''
		Function to get scan start time

		Returns
		-------
		str
			Scan start time in 'YYYY/MM/DD/hh:mm:ss' format
		str 
			Scan end time in 'YYYY/MM/DD/hh:mm:ss' format
		float
			Scan start time in MJD second
		float
			Scan end time in MJD second
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
		Function to get total number of antennas in the measurement set

		Returns
		-------
		int
			Number of antennas
		'''
		self.md.open(self.msname)
		nant=int(self.md.nantennas())
		self.md.close()
		return nant

	def get_antenna_string(self):
		'''
		Function to get total antenna string

		Returns
		-------
		str
			Antenna string
		'''
		nant=self.get_num_antenna()
		ant=''
		for i in range(nant):
			ant+=str(i)+','
		ant=ant[:-1]
		return ant

	def get_freqs(self):
		'''
		Function to return list of channel frequencies

		Returns
		-------
		list
			List of frequencies in Hz
		'''
		self.md.open(self.msname)
		freqs=self.md.chanfreqs(0)
		self.md.close()
		return freqs

	def get_unflag_chan(self,flagfrac=1):
		'''
		Function to get the unflagged channels if flag fraction is less than certain value

		Parameters
		----------
		flagfrac : float 
			Flag fraction per channel (default : 1.0)
		Returns
		-------
		list
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

	def get_unflag_times_mjd(self,flagfrac=1):
		'''
		Function to get the unflagged timestamps in MJD seconds if flag fraction is less than certain value

		Parameters
		----------
		flagfrac : float 
			Flag fraction per channel (default : 1.0)
		Returns
		-------
		list
			List of unflagged times in MJD second
		'''
		self.tb.open(self.msname)
		flag=self.tb.getcol('FLAG')
		self.tb.close()
		ntime=self.get_num_timestamps()
		nbaselines_frac=(flag.shape[-1]/(ntime))-int(flag.shape[-1]/(ntime))
		if nbaselines_frac==0:
			nbaselines=int(flag.shape[-1]/ntime)
			ntime=int(ntime)
		else:
			nbaselines=int(flag.shape[-1]/(ntime-1))
			ntime=int(ntime-1)
		flag=flag.reshape(flag.shape[0],flag.shape[1],ntime,nbaselines)
		unflagged_time=[]
		for t in range(flag.shape[2]):
			flagged_data=flag[:,:,t,:]
			total_data=float(len(flagged_data.flatten()))
			flagged_data=float(np.sum(flagged_data.flatten()))
			if (flagged_data/total_data)<flagfrac:
				unflagged_time.append(t)
		timestamps=self.get_timestamps_in_mjdsecs()[0]
		unflagged_times_mjd=[float(timestamps[i]) for i in unflagged_time]
		return unflagged_times_mjd

	def get_unflag_timestamps(self,flagfrac=1):
		'''
		Function to get the unflagged timestamps in MJD seconds if flag fraction is less than certain value

		Parameters
		----------
		flagfrac : float 
			Flag fraction per channel (default : 1.0)
		Returns
		-------
		list
			List of unflagged channels
		'''
		self.tb.open(self.msname)
		flag=self.tb.getcol('FLAG')
		self.tb.close()
		ntime=self.get_num_timestamps()
		nbaselines_frac=(flag.shape[-1]/(ntime))-int(flag.shape[-1]/(ntime))
		if nbaselines_frac==0:
			nbaselines=int(flag.shape[-1]/ntime)
			ntime=int(ntime)
		else:
			nbaselines=int(flag.shape[-1]/(ntime-1))
			ntime=int(ntime-1)
		flag=flag.reshape(flag.shape[0],flag.shape[1],ntime,nbaselines)
		unflagged_time=[]
		for t in range(flag.shape[2]):
			flagged_data=flag[:,:,t,:]
			total_data=float(len(flagged_data.flatten()))
			flagged_data=float(np.sum(flagged_data.flatten()))
			if (flagged_data/total_data)<flagfrac:
				unflagged_time.append(t)
		timestamps=self.get_timestamps()
		unflagged_timestamps=[timestamps[i] for i in unflagged_time]
		return unflagged_timestamps

	def get_flag_timestamps(self,flagfrac=1):
		'''
		Function to get the flagged timestamps in MJD seconds if flag fraction is less than certain value

		Parameters
		----------
		flagfrac : float
 			Flag fraction per channel (default : 1.0)
		Returns
		-------
		list
			List of flagged timestamps
		'''
		unflagged_time=self.get_unflag_timestamps()
		timestamps=self.get_timestamps()
		flagged_times=[]
		for i in timestamps:
			if i not in unflagged_time:
				flagged_times.append(i)
		return flagged_times
	
	def get_flag_times_mjd(self,flagfrac=1):
		'''
		Function to get the flagged timestamps in MJD seconds if flag fraction is less than certain value
	
		Parameters
		----------
		flagfrac : float 
			Flag fraction per channel (default : 1.0)
		Returns
		-------
		list
			List of flagged times in MJD seconds
		'''
		unflagged_time=self.get_unflag_times_mjd()
		timestamps=self.get_timestamps_in_mjdsecs()[0]
		flagged_times=[]
		for i in timestamps:
			if i not in unflagged_time:
				flagged_times.append(float(i))
		return flagged_times

	def get_model_nomodel_chan(self):
		'''
		Function to get the channels wioth no model data

		Returns
		-------
		list
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

	def get_antenna_id(self,antenna_name=''):
		'''
		Function to get antenna id from antenna name

		Parameters
		----------
		antenna_name : str 
			Name of the antenna
		Returns
		-------
		int
			Antenna id
		'''
		if antenna_name=='':
			return None
		self.md.open(self.msname)
		antenna_names=self.md.antennanames()
		antenna_ids=self.md.antennaids()
		self.md.close()
		if antenna_name not in antenna_names:
			print ('Antenna name not in list.\n')
		else:
			pos=antenna_names.index(antenna_name)
			antenna_id=antenna_ids[pos]
			return antenna_id

	def calc_meanfreq(self):
		'''
		Function to return central frequency of the measurement set (Only MS with single SPW)

		Returns
		-------
		float
			Central frequency in Hz
		'''
		self.md.open(self.msname)
		cent_freq=self.md.meanfreq(0)
		self.md.close()
		return cent_freq

	def calc_meanwavelength(self):
		'''
		Function to return central wavelength of the measurement set in meter
	
		Returns
		-------		
		float
			Central wavelength in metre
		'''
		freq=self.calc_meanfreq()
		wavelength=299792458.0/freq
		return wavelength

	def calc_flagfrac(self):
		'''
		Function to calculate flag fraction

		Returns
		-------
		float
			Fraction of the data flagged
		'''
		self.tb.open(self.msname)
		flag=self.tb.getcol('FLAG')
		flag_frac=np.nansum(flag)/len(flag.flatten())
		self.tb.close()
		return flag_frac

	def radec_sun(self):
		'''
		RA DEC of the Sun at the middle of the measurement set

		Returns
		-------
		str
			RA DEC of the Sun in J2000
		float
			RA in degree 
		float
			DEC in degree
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
		Move the phasecenter of the measurement set at the center of the Sun

		Returns
		-------
		str
			Message of the measurement set phase center chaged at the center of the Sun
		'''
		sun_radec_string,sunra_deg,sundec_deg=self.radec_sun()
		radec_str,ra,dec=self.get_phasecenter()
		code=vishead(vis=self.msname,mode='get',hdkey='fld_code')[0][0]
		code_list=code.split(',')
		if 'FIXVIS' not in code_list:
			if len(code_list)==1 and code_list[0]=='':
				code+='FIXVIS'
			else:
				code+=',FIXVIS'
			if os.path.islink(self.msname):
				msdir=os.path.realpath(self.msname)
				os.unlink(self.msname)
				linked_path=True
			else:
				msdir=self.msname
				linked_path=False
			phaseshift(vis=msdir,outputvis=msdir+'.phaseshift',phasecenter=sun_radec_string)
			if os.path.isdir(msdir):
				os.system('rm -rf '+msdir)
			os.system('mv '+msdir+'.phaseshift '+msdir)
			if linked_path:
				os.system('ln -s '+msdir+' '+self.msname)
			vishead(vis=self.msname,mode='put',hdkey='fld_code',hdvalue=np.array([code]))
			return 'Phasecenter of the observation is moved to Sun center at :'+sun_radec_string+'.\n'
		else:
			return 'Phasecenter is already shifted to the Sun.\n'

	def move_phasecenter_to_source(self,radec=''):
		'''
		Function to move the phasecenter of measurement set at particular RA-DEC

		Parameters
		----------
		radec : str 
			RA DEC string, format 'J2000 00h00m00s +00d00m00s'
		Returns
		-------
		str
			Message of the phasecenter of the measurement set is moved to the new phasecenter
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
				if os.path.islink(self.msname):
					msdir=os.path.realpath(self.msname)
					os.unlink(self.msname)
					linked_path=True
				else:
					msdir=self.msname
					linked_path=False
				phaseshift(vis=msdir,outputvis=msdir+'.phaseshift',phasecenter=sun_radec_string)
				if os.path.isdir(msdir):
					os.system('rm -rf '+msdir)
				os.system('mv '+msdir+'.phaseshift '+msdir)
				if linked_path:
					os.system('ln -s '+msdir+' '+self.msname)
				if len(code_list)==1 and code_list[0]=='':
					code+='FIXVIS_'+radec_string
				else:
					code+=',FIXVIS_'+radec_string
				vishead(vis=self.msname,mode='put',hdkey='fld_code',hdvalue=np.array([code]))
				return 'Move phasecenter to :'+radec+'.\n'
			else:
				return 'Phasecenter is already at :'+radec+'.\n'
		return
			
	def get_max_baseline(self,includeflag=False):
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
		nanpos=np.where(np.isnan(uvdist)==True)
		uvdist=np.delete(uvdist,nanpos)
		return int(np.max(uvdist))

	def get_listobs(self,listobsfile=''):
		'''
		Function to save listobs
	
		Parameters
		----------
		listobsfile : str
			File name to save listobs output
		Returns
		-------
		str
			listobs file name
		'''
		if listobsfile=='':
			listobsfile=self.msname.split('.ms')[0]+'.listobs'
		if os.path.isfile(listobsfile):
			os.system('rm -rf '+listobsfile)
		listobs(vis=self.msname,listfile=listobsfile)
		return listobsfile

	def model_imported(self):
		'''
		Fuction to check whether model imported or model column exists or not
		Returns
		-------
		bool
			True if model exists and False if does not
		''' 
		self.tb.open(self.msname)
		try:
			model=self.tb.getcol('MODEL_DATA')
			if model.shape[0]==2:
				if np.nansum(np.real(model))==model.flatten().shape[0]:
					return False
				else:
					return True
			elif model.shape[0]==4:
				if np.nansum(np.real(model))==model.flatten().shape[0]/2:
					return False
				else:
					return True
			else:
				return False
		except:
			self.tb.close()
			return False	

	
	def make_antenna_list(self,num_bins=5,antenna_list_file=''):
		'''
		Make the antenna addition list splited in number of bins

		Parameters
		----------
		num_bin : int 
			Number of antenna bins
		antenna_list_file : str
 			Name of the file to save antenna list. If an antenna list file given, antenna list will be loaded from that file
		Returns
		-------
		list
			Lists of antennas to be added in steps
		int 
			Number of antenna
		'''
		if antenna_list_file!='':
			try:
				tag,num_ant,antenna_list=np.load(antenna_list_file,allow_pickle=True)
			except:
				tag,num_ant,antenna_list=np.load(antenna_list_file+'.npy',allow_pickle=True)
			return antenna_list,num_ant
		listobs_file=self.msname.split('.ms')[0]+'.listobs'
		if os.path.isfile(listobs_file)==False:
			listobs(vis=self.msname,listfile=listobs_file)
		else:
			pass
		fil=open(listobs_file,'r')
		lines=fil.readlines()
		fil.close()
		start_line_num=0
		for i in range(len(lines)):
			if 'East' in lines[i] and 'North' in lines[i] and 'Elevation' in lines[i]:
				start_line_num=i
			if 'Observation:' in lines[i]:
				instrument=lines[i].split('Observation:')[-1].split('\n')[0].split(' ')[-1]
		for i in range(len(lines)):
			if 'Offset from array center' in lines[i]:
				start_line_num=i+2
		coords=np.genfromtxt(listobs_file,skip_header=start_line_num)[:,6:9]
		cm_coords=np.median(coords,axis=0)
		rel_coords=coords-cm_coords
		dist=np.sqrt(rel_coords[:,0]**2+rel_coords[:,1]**2)	
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
			antenna_list.append(np.sort(pos[0:start+elem].astype('int')+1))
			start+=elem
		if antenna_list_file!='':
			antenna_list_file=antenna_list_file.split('.')[0]
			np.save(antenna_list_file,np.array(['Antennafile',len(coords),antenna_list]))
		return antenna_list,len(coords)

	def get_observatory_loc(self):
		'''
		Give the observatory geodetic location

		Returns
		-------
		float
			Latitude in degree
		float 
			Longitude in degree
		float
			Altitude in meter
		'''
		self.md.open(self.msname)
		pos=self.md.observatoryposition()
		LON=round(np.rad2deg(pos['m0']['value']),2)
		LAT=round(np.rad2deg(pos['m1']['value']),2)
		ALT=pos['m2']['value']
		self.md.close()
		return LAT,LON,ALT

	def get_altaz(self,source_field=0,source_scan=1):
		'''
		Calculate alt az of the phasecenter

		Parameters
		----------
		source_field : int 
			FIELD ID of the source
		source_scan : int 
			Scan number 
		Returns
		-------
		float
			Altitude in radian
		float
			Azimuth in radian
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
		Calculate parallactic for phasecenter at a given Earth location.
		Note = Parallactic angle is defined as the orientation of the sky in telescope coordinate. All angles are defined positive in IAU defined sky coordiunate (North to East).
		So, the rotation of the sky with respect to telescope is negative in the sky coordinate. To account this effect parallactic angle is given in 360-parang form.

		Parameters
		----------
		source_field : int 
			FIELD id of the source
		combine : str 
			Combine 'field' or 'scan'  for calculating parallactic angle or '' for no combine
		Returns
		-------
		float	
			1. combine = 'field', A single parallactic angle for the entire field in degree
		dict
			2. combine = 'scan', A python dictionary {'scan':parang} format
		list
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

	def get_parang(self,ra,dec,source_field=0,combine='field'):
		'''
		Calculate parallactic for phasecenter at a given Earth location for a given RA DEC
		Note = Parallactic angle is defined as the orientation of the sky in telescope coordinate. All angles are defined positive in IAU defined sky coordiunate (North to East).
		So, the rotation of the sky with respect to telescope is negative in the sky coordinate. To account this effect parallactic angle is given in 360-parang form.

		Parameters
		----------
		ra : float
			RA in degree
		dec : float
			DEC in degree
		source_field : int 
			FIELD id of the source
		combine : str 
			Combine 'field' or 'scan' for calculating parallactic angle or '' for no combine
		Returns
		-------
		float
			1. combine = 'field', A single parallactic angle for the entire field in degree
		dict
			2. combine = 'scan', A python dictionary {'scan':parang} format
		list
			3. combine = '', A list of parang for all timestamps
		'''
		LAT,LON,ALT=self.get_observatory_loc()
		radeg=ra
		decdeg=dec
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

	Parameters
	----------
	msname : str 
		Name of the measurement set
	ref_time_chan : bool 
		Whether the time and frequency slice is refernce
	change_msname : bool
		Change the msname or just return the name
	Returns
	-------
	str
		Modified name of the measurement set
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










