'''
Code is written by Devojyoti Kansabanik , 28 Jan, 2021
'''

from optparse import OptionParser
if __name__=='__main__':
	usage= ' P-AIRCARS master controller for each day calibration'
	parser = OptionParser(usage=usage)
	parser.add_option('--basedir',dest="basedir",default=None,help="Name of base directory for a given day",metavar="Directory path")
	(options, args) = parser.parse_args()

import os,sys
os.chdir(options.basedir)
sys.path.append(os.getcwd())
a=os.system('validating_paircars_input')
if os.WEXITSTATUS(a)!=0:
	os._exit(1)
from selfcal_inputs import basedir
from casatools import msmetadata,table,measures,quanta,agentflagger,image,calibrater,ms
from casatasks import *
import logging,numpy as np,copy,glob,psutil,time,subprocess,getpass
from paircars.basic_func import *
from paircars.access_ms import *
from paircars.decor import *
from paircars.flagger import *
from astropy.io import fits
from pathlib import Path
from CALIBRATE.access_calibrate import *

def spliting_timechan(msname,metafits,channel,timestamp,caltype='',ref_timechan=False,input_file='',datacolumn='corrected'):
	'''
	Function to split specific time and frequency slice and keep the necessary files in one directory to run the P-AIRCARS
	Parameters:
	msname = Name of the source measurement set
	metafits = Name of the metafits file
	channel = Channel number to split
	timestamp = Time stamp to split
	caltype = 'B' or 'G' or 'P'
	ref_timechan = True, if True the time and frequency slice is the reference
	input_file = Name of the selfcal input file
	datacolumn = 'corrected', datacolumn to split
	Return:
	Splited measurement set
	'''
	cwd=os.getcwd()
	md=msmetadata()
	md.open(msname)
	freqlist=md.chanfreqs(0)
	md.close()
	inttime=float(fits.getheader(metafits)['INTTIME'])/2
	if type(channel)==str:
		if '~' in channel:
			ch0=int(channel.split('~')[0])
			ch1=int(channel.split('~')[1])
			ch=int((ch0+ch1)/2.0)
			spw='0:'+channel
		else:
			ch=int(channel)
			spw='0:'+str(channel)
	else:
		spw='0:'+str(channel)
		ch=channel
	mjd=timestamp_to_mjdsec(timestamp,format=0)
	t0=mjdsec_to_timestamp(mjd-inttime)
	t1=mjdsec_to_timestamp(mjd+inttime)
	timestamp=t0+'~'+t1
	if ref_timechan==True:
		mainlog.info('Reference channel : '+str(channel)+' at '+str(freqlist[ch]/10**6)+' MHz and reference time : '+str(timestamp)+'\n')
		# Spliting ref time chan ms
		###########################
		mainlog.info('Spliting reference time and channel..............\n')	
		if os.path.isdir(cwd+'/reftimechan.ms')==True:
			os.system('rm -rf '+cwd+'/reftimechan.ms* '+cwd+'/reftimechan.ms.flagversions')
		split(vis=msname,outputvis=cwd+'/reftimechan.ms',datacolumn=datacolumn,spw=spw,timerange=timestamp)
		mainlog.info('split(vis=\''+msname+'\',outputvis=\''+cwd+'/reftimechan.ms\',datacolumn=\''+datacolumn+'\',spw=\''+spw+'\',timerange=\''+timestamp+'\')\n')
		ref_timechan_ms=splited_ms_rename(cwd+'/reftimechan.ms',ref_time_chan=True,change_msname=True)
		ref_timechan_dir=cwd+'/'+os.path.basename(ref_timechan_ms).split('.ms')[0]+'_'+caltype
		if os.path.isdir(ref_timechan_dir)==True:
			os.system('rm -rf '+ref_timechan_dir)
		if os.path.isdir(ref_timechan_dir)==False:
			os.makedirs(ref_timechan_dir) # Making ref time chan directory
		os.system('cp -r '+input_file+' '+ref_timechan_dir)
		if os.path.isdir(ref_timechan_dir+'/'+os.path.basename(ref_timechan_ms))==True:
			os.system('rm -rf '+ref_timechan_dir+'/'+os.path.basename(ref_timechan_ms)+'* '+ref_timechan_dir+'/'+os.path.basename(ref_timechan_ms)+'.flagversions')
		os.system('mv '+ref_timechan_ms+' '+ref_timechan_dir+'/'+os.path.basename(ref_timechan_ms))
		return ref_timechan_dir+'/'+os.path.basename(ref_timechan_ms),ref_timechan_dir
	else:
		# Spliting specific time chan ms
		################################
		mainlog.info('Spliting measurement set for time : '+timestamp+' and frequency : '+str(freqlist[ch]/10**6)+' MHz ............\n')
		if os.path.isdir(cwd+'/timechan.ms')==True:
			os.system('rm -rf '+cwd+'/timechan.ms* '+cwd+'/timechan.ms.flagversions')
		split(vis=msname,outputvis=cwd+'/timechan.ms',datacolumn=datacolumn,spw=spw,timerange=timestamp)
		mainlog.info('split(vis=\''+msname+'\',outputvis=\''+cwd+'/timechan.ms\',datacolumn=\''+datacolumn+'\',spw=\''+spw+'\',timerange=\''+timestamp+'\')\n')
		timechan_ms=splited_ms_rename(cwd+'/timechan.ms',ref_time_chan=False,change_msname=True)
		timechan_dir=cwd+'/'+os.path.basename(timechan_ms).split('.ms')[0]+'_'+caltype
		if os.path.isdir(timechan_dir):
			os.system('rm -rf '+timechan_dir)
		if os.path.isdir(timechan_dir)==False:
			os.makedirs(timechan_dir) # Making ref time chan directory
		os.system('cp -r '+input_file+' '+timechan_dir)
		if os.path.isdir(timechan_dir+'/'+os.path.basename(timechan_ms))==True:
			os.system('rm -rf '+timechan_dir+'/'+os.path.basename(timechan_ms)+'* '+timechan_dir+'/'+os.path.basename(timechan_ms)+'.flagversions')
		os.system('mv '+timechan_ms+' '+timechan_dir+'/'+os.path.basename(timechan_ms))
		return timechan_dir+'/'+os.path.basename(timechan_ms),timechan_dir

# cpulimit check
###########
def CPULIMIT_check():
	a=os.system('cpulimit -h > cpulimit_tmp')
	if a==256:
		os.system('rm -rf cpulimit_tmp')
		return 0
	else:
		os.system('rm -rf cpulimit_tmp')
		return 1


# PAIRCARS master controller
############################
basedir=str(options.basedir)
if basedir[-1]=='/':
	basedir=basedir[:-1]

os.chdir(basedir)
sys.path.append(basedir)
import selfcal_inputs as inputs
if os.path.isdir(inputs.basedir+'/data')==False:
	os.makedirs(inputs.basedir+'/data')

if inputs.use_wsclean==True:
	a=os.system('wsclean > wsclean_test')
	if a!=0:
		print('WSClean is not installed. Using CASA for imaging.\n')
		use_wsclean=False
	else:
		use_wsclean=True
	os.system('rm -rf wsclean_test')
else:
	use_wsclean=inputs.use_wsclean

cpulimit_check=CPULIMIT_check()

# Estimating total casa instances
#################################
available_paircars_instance=get_available_paircars_instance(inputs.paircars_dir,inputs.job_id,inputs.instance)

# Removing screen logs
######################
if inputs.basedir[-1]=='/':
	inputs.basedir=inputs.basedir[:-1]
paircars_base=os.path.dirname(inputs.basedir)
if os.path.isdir(paircars_base+'/Logs_and_Errors/')==False:
	os.makedirs(paircars_base+'/Logs_and_Errors/')
else:
	if os.path.isdir(paircars_base+'/Logs_and_Errors/Logs')==False:
		os.makedirs(paircars_base+'/Logs_and_Errors/Logs')
	else:
		log_list=glob.glob(paircars_base+'/Logs_and_Errors/Logs/*')
		if len(log_list)>0:
			for logfile in log_list:
				if 'Control_' in logfile and str(inputs.job_id) in logfile:
					pass
				else:
					os.system('rm -rf '+logfile)
	if os.path.isdir(paircars_base+'/Logs_and_Errors/Errors')==False:
		os.makedirs(paircars_base+'/Logs_and_Errors/Errors')
	else:
		error_list=glob.glob(paircars_base+'/Logs_and_Errors/Errors/*')
		if len(error_list)>0:
			for errorfile in error_list:
				if 'Control_' in errorfile and str(inputs.job_id) in errorfile:
					pass
				else:
					os.system('rm -rf '+errorfile)


# Logger initiating
###################
formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
mainlog = logging.getLogger('paircars_main_log')
mainlog.setLevel(logging.DEBUG)
console=logging.StreamHandler(sys.stdout)
console.setFormatter(formatter)
mainlog.addHandler(console)
filehandle=logging.FileHandler(paircars_base+'/Logs_and_Errors/Logs/P-AIRCARS_mainlog_'+str(os.path.basename(inputs.basedir))+'.log')
filehandle.setFormatter(formatter)
mainlog.addHandler(filehandle)
mainlog.propagate = False	
os.system('touch '+inputs.basedir+'/.paircars_running')
is_internet=check_internet()

# Deciding bandpass interval
############################
if inputs.quality_factor==0:
	bandpass_interval=0
	if inputs.safety_factor==0:
		calibrator_interval=90
	elif inputs.safety_factor==0:
		calibrator_interval=60
	elif inputs.safety_factor==0:
		calibrator_interval=30
elif inputs.quality_factor==1:
	if inputs.safety_factor==0:
		bandpass_interval=7
		calibrator_interval=60
	elif inputs.safety_factor==1:
		bandpass_interval=5
		calibrator_interval=30
	elif inputs.safety_factor==2:
		bandpass_interval=3
		calibrator_interval=10
elif inputs.quality_factor==2:
	if inputs.safety_factor==0:
		bandpass_interval=3
		calibrator_interval=20
	elif inputs.safety_factor==1:
		bandpass_interval=1
		calibrator_interval=10
	elif inputs.safety_factor==2:
		bandpass_interval=0.5
		calibrator_interval=4

# User given timerange
######################
if inputs.timerange!='':
	time_list=[]
	for i in inputs.timerange.split(','):
		l=i.split('~')
		for j in l:
			time_list.append(j)
	timerange_list_mjdsecs=sorted([float("{:.2f}".format(timestamp_to_mjdsec(i,format=0))) for i in time_list])
	start_time_mjd=min(timerange_list_mjdsecs)
	end_time_mjd=max(timerange_list_mjdsecs)
else:
	start_time_mjd=0
	end_time_mjd=0

# Making local cal database folder
##################################
local_caldatabase=inputs.basedir+'/localdatabase'
if os.path.isdir(local_caldatabase)==False:
	os.makedirs(local_caldatabase)
mainlog.info('Local calibration database is at : '+local_caldatabase+'\n')

# Download available cal files from paircars database #TODO
#####################################################





# Organising ms
###############
updated_MWA_obsids=0
mainlog.info('Organising measurement sets .....\n')
measurement_set_list=glob.glob(inputs.basedir+'/data/*.ms')
for ms in measurement_set_list:
	if os.path.isdir(os.path.realpath(ms))==False:
		os.system('unlink '+ms)
		measurement_set_list.remove(ms)
msfreqs=[float(os.path.basename(a).split('.ms')[0].split('_')[-1]) for a in measurement_set_list]
mstimes=[float(timestamp_to_mjdsec('/'.join(os.path.basename(a).split('time_')[-1].split('_freq')[0].split('_')[:3])+'/'+\
			':'.join(os.path.basename(a).split('time_')[-1].split('_freq')[0].split('_')[3:]))) for a in measurement_set_list]
mstimes_iso=[('-'.join(os.path.basename(a).split('.ms')[0].split('_')[1:4])+' '+':'.join(os.path.basename(a).split('.ms')[0].split('_')[4:7])) for a in measurement_set_list]
ms_OBSIDs=[]
metafits_obsids_msdir=glob.glob(inputs.msdir+'/*.metafits')

if len(metafits_obsids_msdir)!=0:
	for metafits in metafits_obsids_msdir:
		os.system('cp -r '+metafits+' '+inputs.basedir+'/data/'+os.path.basename(metafits))
metafits_obsids=[int(os.path.basename(x).split('.metafits')[0]) for x in glob.glob(inputs.basedir+'/data/*.metafits')]

for i in range(len(measurement_set_list)):
	msname=measurement_set_list[i]
	AMtimerange=AccessMS(msname)
	mjdstamps=AMtimerange.get_timestamps_in_mjdsecs()[0]
	start_mjdtimestamp=min(mjdstamps)
	end_mjdtimestamp=max(mjdstamps)
	allow_ms=False
	if start_time_mjd!=0 and end_time_mjd!=0:
		if start_mjdtimestamp>end_time_mjd or end_mjdtimestamp<start_time_mjd:
			allow_ms=False
		else:
			allow_ms=True
	else:
		allow_ms=True
	if allow_ms==True:
		obsid=get_OBSID_from_ms(msname)
		if obsid==0 or len(metafits_obsids)==0:
			try:
				GPStime=int(Time(mstimes_iso[i],format='iso',scale='utc').gps)
				diff_gpstime=[]
				for a in metafits_obsids:
					if abs(a-GPStime)<480:
						diff_gpstime.append(a)
				if len(diff_gpstime)!=0:
					obsid=np.min(np.array(diff_gpstime))
					ms_OBSIDs.append(obsid)
				elif is_internet==False:
					mainlog.info('No local metafits file is present and no internet connection is available to download metafits file.'+\
									' Either supply metafits file in Data Directory or connect to the internet. Exiting P-AIRCARS....\n')
					os._exit(0)
				else:
					mainlog.info('Trying to download unavailable metafits for OBS ID :'+str(obsid)+' at : '+inputs.basedir+'/data/'+str(obsid)+'.metafits.\n')
					metafits=download_metafits(msname,inputs.basedir+'/data')
					if metafits!=None:
						mainlog.info('Metafits file downloaded at : '+metafits+'\n')
						ms_OBSIDs.append(get_OBSID_from_metafits(metafits))
					else:
						mainlog.info('Metafits file could not be downloaded for ms : '+msname+'. Removing ms from list.\n')
						measurement_set_list.remove(msname)
						msfreqs.remove(msfreqs[i])
						mstimes.remove(mstimes[i])
						mstimes_iso.remove(mstimes_iso[i]) 
			except:
				mainlog.info('Could not connect to MWA metadata server. No metafits files are found in local data directory. Exiting P-AIRCARS.....\n')
				os._exit(0)
		else:
			ms_OBSIDs.append(obsid)
	else:
		measurement_set_list.remove(msname)

if inputs.calc_selfcalib_params==True:
	if inputs.quality_factor==0:
		skip_freq_pol=2560
		if inputs.safety_factor==0:
			skip_time		=	960
			skip_freq		=	1280
		elif inputs.safety_factor==1:
			skip_time		=	720
			skip_freq		=	1280
		else:
			skip_time		=	480
			skip_freq		=	1280
	elif inputs.quality_factor==1:
		skip_freq_pol=1280
		if inputs.safety_factor==0:
			skip_time		=	240
			skip_freq		=	1280
		elif inputs.safety_factor==1:
			skip_time		=	120
			skip_freq		=	1280
		else:
			skip_time		=	60
			skip_freq		=	1280
	else:
		skip_freq_pol=1280
		if inputs.safety_factor==0:
			skip_time		=	80
			skip_freq		=	640
		elif inputs.safety_factor==1:
			skip_time		=	60
			skip_freq		=	640
		else:
			skip_time		=	30
			skip_freq		=	640
else:
	skip_time=inputs.skip_time
	skip_freq=inputs.skip_freq
	if skip_freq>1280:
		mainlog.info('Skip frequency greater than 1.28 MHz is not appropriate for bandpass self-calibration. Using 1.28 MHz for bandpass.\n')
		skip_freq=1280
		inputs.skip_freq=1280


# Making metafits dictionary
############################
ms_gridpoints=[]
ms_OBSIDs_gcal=[]
ms_OBSIDs_bcal=[]
ms_OBSIDs_pcal=[]
final_image_dic={}
if len(measurement_set_list)!=0:
	metafits_dic={}
	for i in range(len(measurement_set_list)):
		msname=measurement_set_list[i]
		metafits=str(ms_OBSIDs[i])+'.metafits'
		if os.path.isfile(inputs.basedir+'/data/'+metafits)==False: 
			mainlog.info('Trying to download unavailable metafits for OBS ID :'+str(ms_OBSIDs[i])+' at : '+inputs.basedir+'/data/'+str(ms_OBSIDs[i])+'.metafits.\n')
			downloaded_metafits=download_metafits(msname,inputs.basedir+'/data')
			if downloaded_metafits!=None:
				metafits=os.path.basename(downloaded_metafits)				
				mainlog.info('Metafits file downloaded at : '+downloaded_metafits+'\n')
				metafits_dic[msname.split('.ms')[0].split('_timesliced')[0].split('_ref')[0]]=inputs.basedir+'/data/'+metafits
				ms_gridpoints.append(fits.getheader(inputs.basedir+'/data/'+metafits)['GRIDNUM'])
			else:
				mainlog.info('Metafits file could not be downloaded for ms : '+msname+' and also not found locally at : '+inputs.basedir+'/data/'+metafits+' for ms : '+msname+'\n')	
				measurement_set_list.remove(msname)
				msfreqs.remove(msfreqs[i])
				mstimes.remove(mstimes[i])
				mstimes_iso.remove(mstimes_iso[i])
				ms_OBSIDs.remove(ms_OBSIDs[i])
		else:
			mainlog.info('Metafits file found at : '+inputs.basedir+'/data/'+metafits+'\n')
			metafits_dic[msname.split('.ms')[0].split('_timesliced')[0].split('_ref')[0]]=inputs.basedir+'/data/'+metafits
			ms_gridpoints.append(fits.getheader(inputs.basedir+'/data/'+metafits)['GRIDNUM'])

	ms_list_copy=copy.deepcopy(measurement_set_list)
	ms_OBSIDs_copy=copy.deepcopy(ms_OBSIDs)
	metafits_dic_copy=copy.deepcopy(metafits_dic)
	ref_timechan_success=False

	ms_list_cals=[]
	ms_OBSIDs_cal=[]
	msfreqs_cal=[]
	mstimes_cal=[]
	for i in range(len(ms_list_copy)):
		if i==0:
			ms_list_cals.append(ms_list_copy[i])
			ms_OBSIDs_cal.append(ms_OBSIDs_copy[i])
			msfreqs_cal.append(msfreqs[i])
			mstimes_cal.append(mstimes[i])
			if inputs.do_polcal==True:
				ms=ms_list_copy[i]
				grid0=ms_gridpoints[i]
		else:
			if inputs.do_polcal==True:
				ms=ms_list_copy[i]
				meta=metafits_dic[ms.split('.ms')[0].split('_timesliced')[0].split('_ref')[0]]
				grid1=ms_gridpoints[i]
			time_diff=np.array(abs(mstimes[i]-np.array(mstimes_cal)))
			freq_diff=np.array(abs(msfreqs[i]-np.array(msfreqs_cal)))
			if (((len(time_diff)==0 and len(freq_diff)==0) or grid0!=grid1) and inputs.do_polcal==True)	or (len(time_diff)==0 and len(freq_diff)==0 and inputs.do_polcal==False):
				ms_list_cals.append(ms_list_copy[i])
				ms_OBSIDs_cal.append(ms_OBSIDs_copy[i])
				msfreqs_cal.append(msfreqs[i])
				mstimes_cal.append(mstimes[i])
			elif (((min(time_diff)>skip_time or min(freq_diff)>(skip_freq/10**3)) or grid0!=grid1) and inputs.do_polcal==True)	or\
				 (min(time_diff)>skip_time or min(freq_diff)>(skip_freq/10**3) and inputs.do_polcal==False):
				ms_list_cals.append(ms_list_copy[i])
				ms_OBSIDs_cal.append(ms_OBSIDs_copy[i])
				msfreqs_cal.append(msfreqs[i])
				mstimes_cal.append(mstimes[i])
		
	ms_list_copy_cals=copy.deepcopy(ms_list_cals)
	ms_OBSIDs_copy_cals=copy.deepcopy(ms_OBSIDs_cal)
	mainlog.info('Measurement sets to perform calibrations : \n')
	for i in ms_list_cals:
		mainlog.info(i+'\n')

	screen_list=[os.path.basename(i) for i in glob.glob('/var/run/screen/S-'+str(getpass.getuser())+'/*')]
	delete_screen_list=[]
	for j in screen_list:
		if str(inputs.job_id) in j and 'P-AIRCARS_mainlog' not in j and 'screen_basedir_for' not in j:
			delete_screen_list.append(j) 
	for i in delete_screen_list:
		os.system('screen -S '+i+' -X quit')

	ref_timechan_done_list=glob.glob(inputs.basedir+'/.ref_timechan_done_*')
	for i in ref_timechan_done_list:
		try:
			ref_timechan_done_msg=int(i.split('_')[-1])
			if ref_timechan_done_msg!=0:
				os.system('rm -rf '+i)
				ref_timechan_done_list.remove(i)
		except:
			os.system('rm -rf '+i)
			ref_timechan_done_list.remove(i)
	if len(ref_timechan_done_list)>=2:
		mainlog.info('Reference time frequency calibration has been done already.\n')
		ref_timechan_done_list=glob.glob(inputs.basedir+'/.ref_timechan_done_*')
		ref_freq_time_msname=[inputs.basedir+'/data/time_'+(a.split('time_')[-1].split('.ms_')[0])+'.ms' for a in ref_timechan_done_list]
		for i in ref_freq_time_msname:
			if i not in measurement_set_list:
				ref_freq_time_msname.remove(i)
		mainlog.info('Reference time frequency measurement set list is : '+str(ref_freq_time_msname)+'\n')
		mainlog.info('######################\n')
		ref_ms_OBSID_list=[]
		ref_time_freq_metafits=[]
		for i in ref_freq_time_msname:
			reftimefreq_ms_OBSID=get_OBSID_from_metafits(metafits_dic[i.split('.ms')[0].split('_timesliced')[0].split('_ref')[0]])
			ref_ms_OBSID_list.append(reftimefreq_ms_OBSID)
			ref_time_freq_metafits.append(metafits_dic[i.split('.ms')[0].split('_timesliced')[0].split('_ref')[0]])
		ref_timechan_success=True
	############################################## # TODO : Do not apply, only make list
	# Applying solutions from calibration database




	# If no calibration is found in calibration database
	# Apply calibrator solution
	####################################################


	calibrator_found=False # TODO : change after done this part
	gaincal_list=[]
	bandpass_list=[]
	caltable_list=[]
	spawned_ms_jobs={}
	finished_ms=[]
	spawned_casa_instances=0
	while ref_timechan_success==False:
		if len(msfreqs)==0 or len(mstimes)==0:
			mainlog.info('No measurement set is present for performing reference time frequency calibration.\n')
			os._exit(1)
		if len(msfreqs_cal)>=3:
			ref_freqstamp=[str(msfreqs_cal[0]),str(msfreqs_cal[int(len(msfreqs_cal)/2)]),str(msfreqs_cal[-1])]
		elif len(msfreqs_cal)==2:
			ref_freqstamp=[str(msfreqs_cal[0]),str(msfreqs_cal[-1])]
		else:
			ref_freqstamp=[str(msfreqs_cal[int(len(msfreqs_cal)/2)])]
		ref_timestamp=mjdsec_to_timestamp(mstimes_cal[int(len(mstimes_cal)/2)],format=3)
		ref_timestamp_mjd=mstimes_cal[int(len(mstimes_cal)/2)]
		# Selecting reference time frequnency ms
		######################################
		ref_freq_time_msname=[]
		ref_time_freq_metafits=[]
		ref_ms_OBSID_list=[]

		for i in range(len(measurement_set_list)):
			msname=measurement_set_list[i]
			for j in ref_freqstamp:
				if j in msname and ref_timestamp in msname:
					ref_freq_time_msname.append(msname)
					ref_time_freq_metafits.append(metafits_dic[msname.split('.ms')[0].split('_timesliced')[0].split('_ref')[0]])
		mainlog.info('Reference time frequency measurement set list is : '+str(ref_freq_time_msname)+'\n')
		mainlog.info('######################\n')

		# Self calibration for reference time frequency ms
		################################################## 
		for i in ref_freq_time_msname:
			reftimefreq_ms_OBSID=get_OBSID_from_metafits(metafits_dic[i.split('.ms')[0].split('_timesliced')[0].split('_ref')[0]])
			ref_ms_OBSID_list.append(reftimefreq_ms_OBSID)
			AM=AccessMS(i)
			workdir=inputs.basedir+'/'+os.path.basename(i).split('.ms')[0]
			if os.path.isdir(workdir)==False:
				os.makedirs(workdir)
		
			touch_file_list=glob.glob(inputs.basedir+'/.ref_timechan_done_*'+str(reftimefreq_ms_OBSID)+'*')
			if len(touch_file_list)!=0:
				for t in touch_file_list:
					msg=t.split('_')[-1]
					try:
						msg=int(msg)
					except:
						pass
					if type(msg)==str:
						os.system('rm -rf '+t)
					else:
						msg=int(msg)
						if msg>100:
							msg-=100
						if msg!=0 and msg!=8 and msg!=9:
							os.system('rm -rf '+t)

			touch_file_list1=glob.glob(inputs.basedir+'/.Finished_*cal*'+str(reftimefreq_ms_OBSID)+'*')
			if len(touch_file_list1)!=0:
				for t in touch_file_list1:
					msg=t.split('_')[-1]
					try:
						msg=int(msg)
					except:
						pass
					if type(msg)==str:
						os.system('rm -rf '+t)
					else:
						msg=int(msg)
						if msg>100:
							msg-=100
						if msg!=0 and msg!=8 and msg!=9:
							os.system('rm -rf '+t)
		
		if calibrator_found==True and len(gaincal_list)!=0: # If calibration found from paircars database or calibrator
			calibrator_OBSID=np.min(np.array([x.split('_')[0] for x in gaincal_list])) # Considering naming convention OBSID_COARSECHAN.gcal
			# Not debugged this part
			ms_gaincal_list=[]
			ms_bandpass_list=[]
			for gcal in gaincal_list:
				if int(gcal.split('_')[0])==int(calibrator_OBSID):
					ms_gaincal_list.append(gcal)
			for bcal in bandpass_list:
				if int(bcal.split('_')[0])==int(calibrator_OBSID):
					ms_bandpass_list.append(bcal)
			bandpass_freq_list=[float(get_caltable_metadata(i)['Channel 0 frequency (MHz)']) for i in ms_bandpass_list]
			bandpass_coarse_chan=np.array([freq_to_MWA_coarse(freq) for freq in bandpass_freq_list])
			reftimefreq_ms_OBSID=ref_ms_OBSID_list[int(len(ref_ms_OBSID_list)/2)]
			single_ref_freq_time_metafits=ref_time_freq_metafits[int(len(ref_time_freq_metafits)/2)]
			ref_freqstamp=ref_freqstamp[int(len(ref_freqstamp)/2)]
			single_ref_freq_time_ms=ref_freq_time_msname[int(len(ref_freq_time_msname)/2)]
			AMref=AccessMS(single_ref_freq_time_ms)
			ms_freq=float(AMref.calc_meanfreq()/10**6)
			ms_coarse_chan=freq_to_MWA_coarse(ms_freq)
			ms_bandpass_list=ms_bandpass_list[np.where(bandpass_coarse_chan==ms_coarse_chan)]
			caltable_list=ms_gaincal_list+ms_bandpass_list
			#
			if abs(reftimefreq_ms_OBSID-calibrator_OBSID)>calibrator_interval*60:
				ref_time_freq=True
			else:
				ref_time_freq=False
			workdir=inputs.basedir+'/'+os.path.basename(single_ref_freq_time_ms).split('.ms')[0]
			if os.path.isdir(workdir)==False:
				if os.path.exists(workdir)==True:
					os.system('rm -rf '+workdir)
				os.makedirs(workdir)
			os.system('cp -r selfcal_inputs.py '+workdir+'/selfcal_inputs.py')
			cmd='run_paircars --msname '+single_ref_freq_time_ms+' --metafits '+single_ref_freq_time_metafits+' --basedir '+inputs.basedir+' --workdir '+workdir+\
				' --ref_freq_avg 0 --ref_time_avg 0 '+' --ref_time_freq True --do_bandpass '+str(inputs.do_bandpass)+' --do_polcal '+str(inputs.do_polcal)\
				+' --cal_attenuation '+str(1.0)+' --scratch True --caltables '+str(','.join(caltable_list))+' --wsclean '\
				+str(use_wsclean)+' --cpu_frac '+str(inputs.cpu_frac)  # TODO : Change calibrator attenuator
			screen_name=str(reftimefreq_ms_OBSID)+'_'+str(os.path.basename(ref_freq_time_msname).split('.ms')[0])+'_runpaircars_'+str(inputs.job_id)
			finished_touch_file=inputs.basedir+'/.Finished_runpaircars_'+str(reftimefreq_ms_OBSID)+'_'+str(single_ref_freq_time_ms.split('.ms')[0])
			screen_batch_file=paircars_instance_runner(cmd,inputs.basedir,inputs.paircars_dir,screen_name,finished_touch_file,inputs.job_id)
			screen_cmd='sh '+screen_batch_file
			mainlog.info(screen_cmd+'\n')
			os.system(screen_cmd)	
			while True:
				touch_file_list=glob.glob(inputs.basedir+'/.ref_timechan_done_'+str(reftimefreq_ms_OBSID)+'_'+str(os.path.basename(single_ref_freq_time_ms))+'*')
				if len(touch_file_list)==1:
					other_error=True
					for i in touch_file_list:
						if '_selfcalerror' in i:
							selfcal_fail=True
						elif '_runtimeerror' in i:
							runtime_fail=True
						elif '_error' in i:
							other_error=True
					if selfcal_fail==True:
						mainlog.info('Reference time frequency calibration failed. Trying with new reference time channel.\n')
						msfreqs_cal.remove(float(ref_freqstamp))
						mstimes_cal.remove(float(ref_timestamp_mjd))
						ref_timechan_success=False
						break
					elif runtime_fail==True:
						mainlog.info('Reference time frequency calibration failed during run time. Some error occured during runtime. Contact developer to fix the problem.\n')
						msfreqs_cal.remove(float(ref_freqstamp))
						mstimes_cal.remove(float(ref_timestamp_mjd))
						ref_timechan_success=False
						break
					elif other_error==True:
						mainlog.info('Reference time frequency calibration failed because some error occured. Contact developer to fix the problem.\n')
						msfreqs_cal.remove(float(ref_freqstamp))
						mstimes_cal.remove(float(ref_timestamp_mjd))
						ref_timechan_success=False
						break
					else:
						mainlog.info('Reference time frequency calibration is finished.\n')
						ref_timechan_success=True
						ms_OBSIDs_gcal.append(single_reftimefreq_ms_OBSID)
						if inputs.do_bandpass==True:
							ms_OBSIDs_bcal.append(single_reftimefreq_ms_OBSID)
						if inputs.do_polcal==True:
							ms_OBSIDs_pcal.append(single_reftimefreq_ms_OBSID)
						result_list=np.load(inputs.basedir+'/Ref_time_freq_slice_output.npy',allow_pickle=True)
						for i in result_list:
							if ref_freq_time_msname in i:
								result=i
								break
						return_msg,ref_time,ref_chan,spawned_casa_instances,ref_freq_avg,ref_time_avg,ms=result 
						spawned_ms_jobs[os.path.basename(ref_freq_time_msname)]=[ref_time,int(ref_chan),int(spawned_casa_instances),float(ref_freq_avg),float(ref_time_avg)]
					break
				else:
					time.sleep(10.0)
					if updated_MWA_obsids==0:
						obsid_file,update_msg=update_mwa_obsids()
						if update_msg==0:
							updated_MWA_obsids=1
		elif calibrator_found==False and len(touch_file_list)==0:
			for j in range(len(ref_freq_time_msname)):
				single_ref_freq_time_ms=ref_freq_time_msname[j]
				single_ref_freq_time_metafits=ref_time_freq_metafits[j]
				single_reftimefreq_ms_OBSID=ref_ms_OBSID_list[j]
				single_ref_freqstamp=ref_freqstamp[j]
				workdir=inputs.basedir+'/'+os.path.basename(single_ref_freq_time_ms).split('.ms')[0]
				if os.path.isdir(workdir)==False:
					if os.path.exists(workdir)==True:
						os.system('rm -rf '+workdir)
					os.makedirs(workdir)
				os.system('cp -r selfcal_inputs.py '+workdir+'/selfcal_inputs.py')
				cmd='run_paircars --msname '+single_ref_freq_time_ms+' --metafits '+single_ref_freq_time_metafits+' --basedir '+inputs.basedir+' --workdir '+workdir+\
					' --ref_freq_avg 0 --ref_time_avg 0 --ref_time_freq True --do_bandpass '+str(inputs.do_bandpass)+' --do_polcal '+str(inputs.do_polcal)\
					+' --scratch True --wsclean '+str(use_wsclean)+' --cpu_frac '+str(inputs.cpu_frac)
				mainlog.info('Command : '+cmd+'\n')
				screen_name=str(single_reftimefreq_ms_OBSID)+'_'+str(os.path.basename(single_ref_freq_time_ms).split('.ms')[0])+'_runpaircars_'+str(inputs.job_id)
				finished_touch_file=inputs.basedir+'/.Finished_runpaircars_'+str(reftimefreq_ms_OBSID)+'_'+str(os.path.basename(single_ref_freq_time_ms).split('.ms')[0])
				screen_batch_file=paircars_instance_runner(cmd,inputs.basedir,inputs.paircars_dir,screen_name,finished_touch_file,inputs.job_id)
				screen_cmd='sh '+screen_batch_file
				mainlog.info(screen_cmd+'\n')
				os.system(screen_cmd)
				time.sleep(1.0)
			while True:
				all_touch_files=[]
				all_success_list=[]
				for j in range(len(ref_freq_time_msname)):
					single_reftimefreq_ms_OBSID=ref_ms_OBSID_list[j]	
					single_ref_freq_time_ms=ref_freq_time_msname[j]
					single_ref_freqstamp=ref_freqstamp[j]
					touch_file_list=glob.glob(inputs.basedir+'/.ref_timechan_done_'+str(single_reftimefreq_ms_OBSID)+'_*')
					runtime_fail=False
					selfcal_fail=False
					other_error=False
					if len(touch_file_list)>0:
						for j in touch_file_list:
							if j not in all_touch_files:
								all_touch_files.append(j)
						i=touch_file_list[0]
						if '_selfcalerror' in i:
							selfcal_fail=True
						elif '_runtimeerror' in i:
							runtime_fail=True
						elif '_error' in i:
							other_error=True
						else:
							all_success_list.append(single_ref_freq_time_ms)
						if selfcal_fail==True:
							mainlog.info('Reference time frequency calibration failed. Trying with new reference time channel.\n')
							msfreqs_cal.remove(float(single_ref_freqstamp))
							break
						elif runtime_fail==True:
							mainlog.info('Reference time frequency calibration failed during run time. Some error occured during runtime. Contact developer to fix the problem.\n')
							msfreqs_cal.remove(float(single_ref_freqstamp))
							mstimes_cal.remove(float(ref_timestamp_mjd))
							ref_timechan_success=False
							break
						elif other_error==True:				
							mainlog.info('Reference time frequency calibration failed because some error occured. Contact developer to fix the problem.\n')
							msfreqs_cal.remove(float(single_ref_freqstamp))
							break
				if len(all_touch_files)==len(ref_freq_time_msname):
					if len(all_success_list)>=2 and len(ref_freq_time_msname)!=1:						
						mainlog.info('Reference time frequency calibration is finished for '+str(len(all_success_list))+' reference frequency slices.\n') 
						ref_timechan_success=True
						ms_OBSIDs_gcal.append(single_reftimefreq_ms_OBSID)
						if inputs.do_bandpass==True:
							ms_OBSIDs_bcal.append(single_reftimefreq_ms_OBSID)
						if inputs.do_polcal==True:
							ms_OBSIDs_pcal.append(single_reftimefreq_ms_OBSID)
						break
					elif len(all_success_list)==0:
						mainlog.info('Reference time frequency calibration has been failed for all '+str(len(all_success_list))+' reference frequency slices.\n')
						mainlog.info('Trying with new reference ms.\n') 
						mstimes_cal.remove(float(ref_timestamp_mjd))
						ref_timechan_success=False
						break
					elif len(ref_freq_time_msname)==1 and len(all_success_list)==1:
						mainlog.info('One reference frequency ms was given. Reference time frequency calibration is finished for '\
										+str(len(all_success_list))+' reference frequency slice.\n') 
						ref_timechan_success=True
						ms_OBSIDs_gcal.append(single_reftimefreq_ms_OBSID)
						if inputs.do_bandpass==True:
							ms_OBSIDs_bcal.append(single_reftimefreq_ms_OBSID)
						if inputs.do_polcal==True:
							ms_OBSIDs_pcal.append(single_reftimefreq_ms_OBSID)
						break
				else:
					time.sleep(10.0)
					if updated_MWA_obsids==0:
						obsid_file,update_msg=update_mwa_obsids()
						if update_msg==0:
							updated_MWA_obsids=1
		else: # Checked and debugged
			mainlog.info('Reference time frequency calibration is done already.\n')
			for j in range(len(ref_freq_time_msname)):
				single_ref_freq_time_ms=ref_freq_time_msname[j]
				single_ref_freq_time_metafits=ref_time_freq_metafits[j]
				single_reftimefreq_ms_OBSID=ref_ms_OBSID_list[j]
				single_ref_freqstamp=ref_freqstamp[j]
				msfreqs_cal.remove(float(single_ref_freqstamp))
				mstimes_cal.remove(float(ref_timestamp_mjd))
				result_list=np.load(inputs.basedir+'/Ref_time_freq_slice_output.npy',allow_pickle=True)
				for i in result_list:
					if single_ref_freq_time_ms in i:
						result=i
						break
				return_msg,ref_time,ref_chan,spawned_casa_instances,ref_freq_avg,ref_time_avg,ms=result 
				spawned_ms_jobs[os.path.basename(single_ref_freq_time_ms)]=[ref_time,int(ref_chan),int(spawned_casa_instances),float(ref_freq_avg),float(ref_time_avg)]
			ref_timechan_success=True
			break

	for j in range(len(ref_freq_time_msname)):
		index=measurement_set_list.index(ref_freq_time_msname[j])
		measurement_set_list.remove(ref_freq_time_msname[j]) # Removing ref time freq ms from ms list
		del metafits_dic[ref_freq_time_msname[j].split('.ms')[0].split('_timesliced')[0].split('_ref')[0]] # Removing ref time freq metafits from list 

	for i in measurement_set_list: # Removing non-fininshed ms directories
		a=glob.glob(inputs.basedir+'/.Finished_spawned_*_'+os.path.basename(i)+'*')
		if len(a)>0:
			measurement_set_list.remove(i)
		else:
			if os.path.isdir(inputs.basedir+'/'+os.path.basename(i).split('.ms')[0]):
				os.system('rm -rf '+inputs.basedir+'/'+os.path.basename(i).split('.ms')[0])
				os.system('rm -rf '+inputs.basedir+'/.Finished_spawned_*'+os.path.basename(i).split('.ms')[0]+'*')
				os.system('rm -rf '+inputs.basedir+'/.Finished_runpaircars_*'+os.path.basename(msname).split('.ms')[0]+'*')

	time.sleep(10)
	result_list=np.load(inputs.basedir+'/Ref_time_freq_slice_output.npy',allow_pickle=True)
	keys=spawned_ms_jobs.keys()
	for i in result_list:
		return_msg,ref_time,ref_chan,spawned_instances,ref_freq_avg,ref_time_avg,ms=i
		if os.path.basename(ms) not in keys:
			spawned_ms_jobs[os.path.basename(ms)]=[ref_time,int(ref_chan),int(spawned_instances),float(ref_freq_avg),float(ref_time_avg)]

	mainlog.info('Remaining measurement sets : '+','.join(measurement_set_list)+'\n')
	os.system('rm -rf '+inputs.basedir+'/Nonref_time_freq_slice_output.npy')
	if len(measurement_set_list)>0:
		result=np.load(inputs.basedir+'/Ref_time_freq_slice_output.npy',allow_pickle=True)[0]
		ref_time=result[1]
		ref_chan=result[2]
		ref_freq_avg=result[-3]
		ref_time_avg=result[-2]
		ref_caltable_list=[]
		cal_dirs=glob.glob(inputs.basedir+'/caltables/'+str(reftimefreq_ms_OBSID)+'/*')
		for i in cal_dirs:
			ref_cal=glob.glob(i+'/*ref*.cal')
			if len(ref_cal)!=0:
				ref_caltable_list.append(ref_cal[0]) # Ref time caltables
		caltable_str=','.join(ref_caltable_list)
		mainlog.info('Reference caltables : '+caltable_str+'\n')
		ref_freq_list=[]
		for j in range(len(ref_freq_time_msname)):
			single_reftimefreq_ms=ref_freq_time_msname[j]
			AMref=AccessMS(single_reftimefreq_ms)
			ref_freq_list.append(AMref.calc_meanfreq()/10**6)
		ref_freq_list=np.array(ref_freq_list)
		cent_ref_freq=ref_freq_list[int(len(ref_freq_list)/2)]
		c=0
		num_ms_jobs=0
		for i in range(len(measurement_set_list)):
			available_paircars_instance=get_available_paircars_instance(inputs.paircars_dir,inputs.job_id,inputs.instance)
			msname=measurement_set_list[i]
			obsid=ms_OBSIDs_cal[i]
			AMms=AccessMS(msname)
			ms_freq=AMms.calc_meanfreq()/10**6
			freq_diff=np.min(abs(ms_freq-ref_freq_list))
			time_diff=abs(reftimefreq_ms_OBSID-obsid)
			if time_diff>bandpass_interval*3600 and inputs.quality_factor!=0:
				mainlog.info('Performing bandpass as time interval is greater than : '+str(bandpass_interval)+' hr.\n')
				do_bandpass=True
			elif freq_diff>(skip_freq/10**3) and inputs.quality_factor!=0:
				mainlog.info('Performing bandpass as frequency interval is greater than : '+str(skip_freq)+' kHz.\n')
				do_bandpass=True
			else:
				mainlog.info('Do not perform bandpass, becuase either time interval is smaller than : '+str(bandpass_interval)+' or frequency is smaller than '+str(skip_freq)\
							+' kHz or quality_factor = 0.\n')
				do_bandpass=False
			metafits=metafits_dic[msname.split('.ms')[0].split('_timesliced')[0].split('_ref')[0]]
			if inputs.do_polcal==True:
				ref_time_freq_gridpoint=fits.getheader(ref_time_freq_metafits[0])['GRIDNUM']
				ms_gridpoint=fits.getheader(metafits)['GRIDNUM']
				if ref_time_freq_gridpoint!=ms_gridpoint:
					do_polcal=True
					mainlog.info('Beam pointing changed. Performing polarisation calibration.\n')
				elif freq_diff>(skip_freq_pol/10**3):
					do_polcal=True
					mainlog.info('Frequency difference is more than polarisation frequency interval. Performing polarisation calibration.\n')
				elif time_diff>(bandpass_interval*2)*3600:
					do_polcal=True
					mainlog.info('Time difference is more than polarisation time interval. Performing polarisation calibration.\n')
				else:
					mainlog.info('Beam pointing is same. Do not perform polarisation calibration.\n')
					do_polcal=False
			if (do_bandpass==False and do_polcal==False and time_diff>0 and abs(ms_freq-cent_ref_freq)<=1280) or (do_bandpass==True or do_polcal==True):
				if do_bandpass==False and do_polcal==False and time_diff>0 and abs(ms_freq-cent_ref_freq)<=1280:
					mainlog.info('Performing gaincal for other times at central reference frequency : '+str(cent_ref_freq/1000.0)+' MHz')
				workdir=inputs.basedir+'/'+os.path.basename(msname).split('.ms')[0]
				if os.path.isdir(workdir)==False:
					if os.path.isfile(workdir)==True:
						os.system('rm -rf '+workdir)
					os.makedirs(workdir)
				os.system('cp -r selfcal_inputs.py '+workdir+'/selfcal_inputs.py')
				if calibrator_found==True and (len(gaincal_list)!=0 or len(bandpass_list)!=0): # Not debugged
					calibrator_OBSID=np.min(np.array([x.split('_')[0] for x in gaincal_list])) # Considering naming convention OBSID_COARSECHAN.gcal
					ms_gaincal_list=[]
					ms_bandpass_list=[]
					for gcal in gaincal_list:
						if int(gcal.split('_')[0])==int(calibrator_OBSID):
							ms_gaincal_list.append(gcal)
					for bcal in bandpass_list:
						if int(bcal.split('_')[0])==int(calibrator_OBSID):
							ms_bandpass_list.append(bcal)
					bandpass_freq_list=[float(get_caltable_metadata(i)['Channel 0 frequency (MHz)']) for i in ms_bandpass_list]
					bandpass_coarse_chan=np.array([freq_to_MWA_coarse(freq) for freq in bandpass_freq_list])
					AMref=AccessMS(msname)
					ms_freq=float(AMref.calc_meanfreq()/10**6)
					ms_coarse_chan=freq_to_MWA_coarse(ms_freq)
					ms_bandpass_list=ms_bandpass_list[np.where(bandpass_coarse_chan==ms_coarse_chan)]
					caltable_list=ms_gaincal_list+ms_bandpass_list
				else:
					caltable_list=[]
				mainlog.info('Interpolating gaincal table to : '+inputs.basedir+'/caltables/'+str(reftimefreq_ms_OBSID)+'/Interp_'+str(ms_freq)+'MHz.gcalp.\n')
				interpolated_caltable=multifreq_gaincal_interpolate(gaintables=ref_caltable_list,outputfreq=ms_freq,\
								output_gaintable=inputs.basedir+'/caltables/'+str(reftimefreq_ms_OBSID)+'/Interp_'+str(ms_freq)+'MHz.gcalp')
				if len(caltable_list)!=0:
					interpolated_caltable=interpolated_caltable+caltable_list
				ms_obsid=get_OBSID_from_metafits(metafits)
				ms_OBSIDs_gcal.append(ms_obsid)
				if do_bandpass==True:
					ms_OBSIDs_bcal.append(ms_obsid)
				if do_polcal==True:
					ms_OBSIDs_pcal.append(ms_obsid)
				cmd='run_paircars --msname '+msname+' --metafits '+metafits+' --basedir '+inputs.basedir+' --workdir '+workdir+' --ref_freq_avg '+str(ref_freq_avg)+\
					' --ref_time_avg '+str(ref_time_avg)+' --ref_time_freq False --do_bandpass '+str(do_bandpass)+' --do_polcal '+str(do_polcal)+\
					' --scratch False --caltables '+interpolated_caltable+' --wsclean '+str(use_wsclean)+' --cpu_frac '+str(inputs.cpu_frac)
				mainlog.info('Command : '+cmd+'\n')
				screen_name=str(ms_obsid)+'_'+str(os.path.basename(msname).split('.ms')[0])+'_runpaircars_'+str(inputs.job_id)
				finished_touch_file=inputs.basedir+'/.Finished_runpaircars_'+str(reftimefreq_ms_OBSID)+'_'+str(os.path.basename(msname).split('.ms')[0])
				screen_batch_file=paircars_instance_runner(cmd,inputs.basedir,inputs.paircars_dir,screen_name,finished_touch_file,inputs.job_id)
				screen_cmd='sh '+screen_batch_file
				mainlog.info(screen_cmd+'\n')
				os.system(screen_cmd)
				time.sleep(1.0)
				AMmsname=AccessMS(msname)
				BW=AMmsname.calc_bandwidth()/10**3 # In kHz
				TW=AMmsname.calc_total_time()
				spawned_casa_instances+=int(TW/skip_time)
				num_ms_jobs=+int(TW/skip_time)
				if do_bandpass==True:
					num_ms_jobs=+int(BW/skip_freq)
					spawned_casa_instances+=int(BW/skip_freq)
				if do_polcal==True:
					num_ms_jobs=+int(BW/skip_freq_pol)
					spawned_casa_instances+=int(BW/skip_freq_pol)
				spawned_ms_jobs[os.path.basename(msname)]=[ref_time,int(ref_chan),int(num_ms_jobs),float(ref_freq_avg),float(ref_time_avg)]
				available_paircars_instance=get_available_paircars_instance(inputs.paircars_dir,inputs.job_id,inputs.instance)
				basemsdir=os.path.basename(msname).split('.ms')[0]
				touch_file_list=glob.glob(inputs.basedir+'/.Finished_*cal*'+str(obsid)+'*'+basemsdir+'*')
				c+=1	
				if c>=len(measurement_set_list):
					break
				while True:
					if available_paircars_instance>1:
						print('Available P-AIRCARS instance : '+str(available_paircars_instance)+'\n')
						mainlog.info('P-AIRCARS instance is available. Spawn new job...\n')
						break
					else:
						available_paircars_instance=get_available_paircars_instance(inputs.paircars_dir,inputs.job_id,inputs.instance)
						time.sleep(10.0)
						if updated_MWA_obsids==0:
							obsid_file,update_msg=update_mwa_obsids()
							if update_msg==0:
								updated_MWA_obsids=1	
			else:
				mainlog.info('Passing the ms as it is not in same coarse channel of central reference frequency.\n')			

	mainlog.info('Calibration jobs for all measurement sets have been started. Waiting to spawn all calibration instances....\n')
	ms_list=copy.deepcopy(ms_list_copy)
	ms_OBSIDs=copy.deepcopy(ms_OBSIDs_copy)
	metafits_dic=copy.deepcopy(metafits_dic_copy)

	# Estimating total casa instances
	#################################
	available_paircars_instance=get_available_paircars_instance(inputs.paircars_dir,inputs.job_id,inputs.instance)
	mainlog.info('Total number of available CASA instances : '+str(int(inputs.instance))+'\n')
	mainlog.info('Available P-AIRCARS instances : '+str(int(available_paircars_instance))+'\n')
	while True:
		available_paircars_instance=get_available_paircars_instance(inputs.paircars_dir,inputs.job_id,inputs.instance)
		time.sleep(20)
		print('Running P-AIRCARS instances : '+str(int(inputs.instance)-int(available_paircars_instance))+'\n')
		if available_paircars_instance>1:
			print('Available P-AIRCARS instance : '+str(available_paircars_instance)+'\n')
			break

	# Spawning jobs for making final caltables and images
	#####################################################
	for msname in ms_list: 
		metafits=metafits_dic[msname.split('.ms')[0].split('_timesliced')[0].split('_ref')[0]]
		OBSID=get_OBSID_from_metafits(metafits)
		keys=final_image_dic.keys()
		if OBSID not in keys:
			final_image_dic[OBSID]=[msname]
		else:
			x=final_image_dic[OBSID]
			x.append(msname)
			final_image_dic[OBSID]=x
		for i in glob.glob(inputs.basedir+'/imagemodels/*'):
			if i not in ms_OBSIDs_gcal:
				ms_OBSIDs_gcal.append(int(os.path.basename(i)))
		for i in glob.glob(inputs.basedir+'/bpimagemodels/*'):
			if i not in ms_OBSIDs_bcal:
				ms_OBSIDs_bcal.append(int(os.path.basename(i)))
		for i in glob.glob(inputs.basedir+'/polimagemodels/*'):
			if i not in ms_OBSIDs_pcal:
				ms_OBSIDs_pcal.append(int(os.path.basename(i)))
		gaincal_modeldirs=','.join([inputs.basedir+'/imagemodels/'+str(i) for i in ms_OBSIDs_gcal])
		bandpass_modeldirs=','.join([inputs.basedir+'/bpimagemodels/'+str(i) for i in ms_OBSIDs_bcal])
		polcal_modeldirs=','.join([inputs.basedir+'/polimagemodels/'+str(i) for i in ms_OBSIDs_pcal])
		mainlog.info('Making final calibration tables for : '+msname+'\n')
		num_jobs=spawned_ms_jobs[os.path.basename(msname)][2]
		freq_avg=spawned_ms_jobs[os.path.basename(msname)][3]
		time_avg=spawned_ms_jobs[os.path.basename(msname)][4]
		screen_name=str(OBSID)+'_'+str(os.path.basename(msname).split('.ms')[0])+'_final_'+str(inputs.job_id)
		batch_file=inputs.basedir+'/'+screen_name+'.batch'
		cmd_batch_file=inputs.basedir+'/'+screen_name+'_cmd.batch'
		cmd='manage_database --msname '+msname+' --metafits '+metafits+' --num_jobs '+str(num_jobs)+' --basedir '+basedir+\
			' --gaincal_modeldir '+gaincal_modeldirs+' --bandpass_modeldir '+bandpass_modeldirs+' --polcal_modeldir '+polcal_modeldirs+' --localdatabase '+local_caldatabase+\
			' --freqavg '+str(freq_avg)+' --timeavg '+str(time_avg)+' --inputfile '+basedir+'/selfcal_inputs.py --verbose '+str(inputs.verbose)
		outputfile_log=paircars_base+'/Logs_and_Errors/Logs/'+screen_name+'.log'
		outputfile_error=paircars_base+'/Logs_and_Errors/Errors/'+screen_name+'.error'
		cmd='echo \"'+cmd+'\" > '+cmd_batch_file+';sleep 2; chmod a+rwx '+cmd_batch_file+'; sleep 2; nohup sh '+cmd_batch_file+' > '+outputfile_log+' 2>'+\
			outputfile_error+' < /dev/null &; sleep 2'
		cmd=cmd.split(';')
		cmds=[i+'\n' for i in cmd]
		if os.path.exists(cmd_batch_file):
			os.system('rm -rf '+cmd_batch_file)
		if os.path.exists(batch_file):
			os.system('rm -rf '+batch_file)
		if os.path.isfile(batch_file):
			fil=open(batch_file,'r+')
		else:
			fil=open(batch_file,'w')
		fil.writelines(cmds)
		fil.close()
		os.system('chmod a+rwx '+batch_file)
		screen_cmd='sh '+batch_file
		os.system(screen_cmd)
		mainlog.info(screen_cmd+'\n')
		while True:
			available_paircars_instance=get_available_paircars_instance(inputs.paircars_dir,inputs.job_id,inputs.instance)
			time.sleep(20)
			print('Running P-AIRCARS instances : '+str(int(inputs.instance)-int(available_paircars_instance))+'\n')
			if available_paircars_instance>1:
				print('Available P-AIRCARS instance : '+str(available_paircars_instance)+'\n')
				break
	basemsdir=os.path.basename(msname).split('.ms')[0]
	if updated_MWA_obsids==0:
		obsid_file,update_msg=update_mwa_obsids()
		if update_msg==0:
			updated_MWA_obsids=1
	if updated_MWA_obsids==0:
		obsid_file,update_msg=update_mwa_obsids()
		if update_msg==0:
			updated_MWA_obsids=1		
else:
	if updated_MWA_obsids==0:
		obsid_file,update_msg=update_mwa_obsids()
		if update_msg==0:
			updated_MWA_obsids=1
	if inputs.timerange!='':
		mainlog.info('No measurement set is present in the timerange : '+inputs.timerange+'\n')
	else:
		mainlog.info('No measurment set is found.\n')
if inputs.send_notification:
	keys=final_image_dic.keys()
	for OBSID in keys:
		ms_list=final_image_dic[OBSID]
		cmd='track_final_imaging --basedir '+str(inputs.basedir)+' --savedir '+str(inputs.final_image_dir)+' --num_ms '+str(len(ms_list))+' --OBSID '+str(OBSID)\
				+' --email '+str(inputs.email)
		screen_name=str(OBSID)+'_track_final_imaging'
		outputfile_log=paircars_base+'/Logs_and_Errors/Logs/'+screen_name+'.log'
		outputfile_error=paircars_base+'/Logs_and_Errors/Errors/'+screen_name+'.error'
		batch_file=inputs.basedir+'/'+screen_name+'.batch'
		cmd_batch_file=inputs.basedir+'/'+screen_name+'_cmd.batch'
		cmd='echo \"'+cmd+'\" > '+cmd_batch_file+';sleep 2; chmod a+rwx '+cmd_batch_file+'; sleep 2; nohup sh '+cmd_batch_file+' > '+outputfile_log+' 2>'+\
			outputfile_error+' < /dev/null &; sleep 2'
		cmd=cmd.split(';')
		cmds=[i+'\n' for i in cmd]
		if os.path.exists(cmd_batch_file):
			os.system('rm -rf '+cmd_batch_file)
		if os.path.exists(batch_file):
			os.system('rm -rf '+batch_file)
		if os.path.isfile(batch_file):
			fil=open(batch_file,'r+')
		else:
			fil=open(batch_file,'w')
		fil.writelines(cmds)
		fil.close()
		os.system('chmod a+rwx '+batch_file)
		screen_cmd='sh '+batch_file
		os.system(screen_cmd)
		mainlog.info(screen_cmd+'\n')

# Applying solution to whole ms
############################### 

# Final imaging mode
####################
# TODO : include mwa_hyperbeam instead of mwa_pb
