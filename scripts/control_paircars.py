'''
Code is written by Devojyoti Kansabanik , 28 Jan, 2021
'''

from optparse import OptionParser
if __name__=='__main__':
	usage= ' PAIRCARS master controller for each day calibration'
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
from casatools import *
from casatasks import *
import logging,numpy as np,copy,glob,psutil,time,subprocess,getpass
from paircars.basic_func import *
from paircars.access_ms import *
from paircars.decor import *
from paircars.flagger import *
from astropy.io import fits
from CALIBRATE.access_calibrate import *

def spliting_timechan(msname,metafits,channel,timestamp,caltype='',ref_timechan=False,input_file='',datacolumn='corrected'):
	'''
	Function to split specific time and frequency slice and keep the necessary files in one directory to run the PAIRCARS
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
# MPI check
###########
def MPI_check():
	a=subprocess.getstatusoutput('mpirun -h')[0]
	if a==0:
		return 0
	else:
		return 1

def casa_instance_runner(cmd,screen_name,finished_touch_file,prefix_cmds=[]):
	'''
	Function to run a casa instance
	Parameters:
	cmd = Command to run
	screen_name = Name of the screen
	'''
	batch_file=inputs.basedir+'/'+screen_name+'.batch'
	cmd_batch=inputs.basedir+'/'+screen_name+'_cmd.batch'
	cmd+=';sleep 2 ;if ! ls '+finished_touch_file+'_* ; then  touch '+finished_touch_file+'_error ;  fi'
	cmd='screen -S '+screen_name+' -X quit; sleep 2; screen -mdS '+screen_name+';sleep 2; echo \"'+cmd+'\" > '+cmd_batch+';sleep 2; chmod a+rwx '+cmd_batch+\
			'; sleep 2; screen -S '+screen_name+' -X stuff \"sh '+cmd_batch+' \\n\"; sleep 2'
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
	return inputs.basedir+'/'+screen_name+'.batch'


# PAIRCARS master controller
############################

basedir=str(options.basedir)
if basedir[-1]=='/':
	basedir=basedir[:-1]

if os.path.isdir(basedir+'/data')==False:
	os.makedirs(basedir+'/data')

os.chdir(basedir)
sys.path.append(basedir)
import selfcal_inputs as inputs

mpi=MPI_check()

# Estimating total casa instances
#################################
total_available_cpu=psutil.cpu_count()-(psutil.cpu_count()*psutil.cpu_percent()/100.0)
available_cpu_for_paircars=int(total_available_cpu*inputs.cpu_frac)
total_cpu_frac=int(psutil.cpu_count()*inputs.cpu_frac)
if available_cpu_for_paircars>total_cpu_frac:
	available_cpu_for_paircars=total_cpu_frac
cpu_sockets =  int(subprocess.check_output('cat /proc/cpuinfo | grep "physical id" | sort -u | wc -l', shell=True))

# Logger initiating
###################
formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
mainlog = logging.getLogger('paircars_main_log')
mainlog.setLevel(logging.DEBUG)
console=logging.StreamHandler(sys.stdout)
console.setFormatter(formatter)
mainlog.addHandler(console)
filehandle=logging.FileHandler(inputs.basedir+'/PAIRCARS_mainlog.log')
filehandle.setFormatter(formatter)
mainlog.addHandler(filehandle)
mainlog.propagate = False	
os.system('touch '+inputs.basedir+'/.paircars_running')

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
	timerange_list_mjdsecs=sorted([float("{:.2f}".format(timestamp_to_mjdsec(i,format=2))) for i in time_list])
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
	mjdstamps=AMtimerange.get_timestamps_in_mjdsecs()
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
				mainlog.info('Could not connect to MWA metadata server. No metafits files are found in local data directory. Exiting PAIRCARS.....\n')
				os._exit(0)
		else:
			ms_OBSIDs.append(obsid)
	else:
		measurement_set_list.remove(msname)

if inputs.calc_selfcalib_params==True:
	if inputs.quality_factor==0:
		skip_freq_pol=1280
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


# Making metafits dictionary
############################
ms_gridpoints=[]
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
	for i in ms_OBSIDs_cal:
		for j in screen_list:
			if str(i)+'_' in j:
				delete_screen_list.append(j) 
	for i in delete_screen_list:
		os.system('screen -S '+i+' -X quit')


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

		############################################## # TODO : Do not apply, only make list
		# Applying solutions from calibration database




		# If no calibration is found in calibration database
		# Apply calibrator solution
		####################################################


		calibrator_found=False # TODO : change after done this part
		gaincal_list=[]
		# Self calibration for reference time frequency ms
		################################################## 
		for i in ref_freq_time_msname:
			reftimefreq_ms_OBSID=get_OBSID_from_metafits(metafits_dic[i.split('.ms')[0].split('_timesliced')[0].split('_ref')[0]])
			ref_ms_OBSID_list.append(reftimefreq_ms_OBSID)
			AM=AccessMS(i)
			workdir=inputs.basedir+'/'+os.path.basename(i).split('.ms')[0]
			if os.path.isdir(workdir)==False:
				os.makedirs(workdir)
			
			touch_file_list=glob.glob(inputs.basedir+'/.ref_timechan_done_*'+str(ref_ms_OBSID_list)+'*')
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

			touch_file_list1=glob.glob(inputs.basedir+'/.Finished_*cal*'+str(ref_ms_OBSID_list)+'*')
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
		caltable_list=[]
		spawned_ms_jobs={}
		finished_ms=[]

		if calibrator_found==True and len(gaincal_list)!=0: # If calibration found from paircars database or calibrator
			calibrator_OBSID=np.min(np.array([x.split('_')[0] for x in gaincal_list]))
			caltable_list=gaincal_list+bandpass_list
			reftimefreq_ms_OBSID=ref_ms_OBSID_list[int(len(ref_ms_OBSID_list)/2)]
			single_ref_freq_time_metafits=ref_time_freq_metafits[int(len(ref_time_freq_metafits)/2)]
			ref_freqstamp=ref_freqstamp[int(len(ref_freqstamp)/2)]
			single_ref_freq_time_ms=ref_freq_time_msname[int(len(ref_freq_time_msname)/2)]
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
				+' --num_threads '+str(available_cpu_for_paircars)+' --cal_attenuation '+str(1.0)+' --caltables '+str(','.join(caltable_list)) # TODO : Change calibrator attenuator
			screen_name=str(reftimefreq_ms_OBSID)+'_'+str(os.path.basename(ref_freq_time_msname).split('.ms')[0])+'_runpaircars'
			finished_touch_file=inputs.basedir+'/.Finished_runpaircars_'+str(reftimefreq_ms_OBSID)+'_'+str(single_ref_freq_time_ms.split('.ms')[0])
			screen_batch_file=casa_instance_runner(cmd,screen_name,finished_touch_file)
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
						result_list=np.load(inputs.basedir+'/Ref_time_freq_slice_output.npy',allow_pickle=True)
						for i in result_list:
							if ref_freq_time_msname in i:
								result=i
								break
						return_msg,ref_time,ref_chan,spawned_casa_instances,ref_freq_avg,ref_time_avg,ms=result 
						spawned_ms_jobs[ref_freq_time_msname]=[ref_time,ref_chan,spawned_casa_instances,ref_freq_avg,ref_time_avg]
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
					+' --num_threads '+str(available_cpu_for_paircars)
				screen_name=str(single_reftimefreq_ms_OBSID)+'_'+str(os.path.basename(single_ref_freq_time_ms).split('.ms')[0])+'_runpaircars'
				finished_touch_file=inputs.basedir+'/.Finished_runpaircars_'+str(reftimefreq_ms_OBSID)+'_'+str(os.path.basename(single_ref_freq_time_ms).split('.ms')[0])
				screen_batch_file=casa_instance_runner(cmd,screen_name,finished_touch_file)
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
					touch_file_list=glob.glob(inputs.basedir+'/.ref_timechan_done_'+str(single_reftimefreq_ms_OBSID)+'_'+str(os.path.basename(single_ref_freq_time_ms))+'*')
					if len(touch_file_list)>0:
						if touch_file_list[0] not in all_touch_files:
							all_touch_files.append(touch_file_list[0])
						runtime_fail=False
						selfcal_fail=False
						other_error=True
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
					elif len(all_touch_files)==len(ref_freq_time_msname):
						if len(all_success_list)>=2 and len(ref_freq_time_msname)!=1:						
							mainlog.info('Reference time frequency calibration is finished for '+str(len(all_success_list))+' reference frequency slices.\n') 
							ref_timechan_success=True
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
							break
				else:
					time.sleep(10.0)
					if updated_MWA_obsids==0:
						obsid_file,update_msg=update_mwa_obsids()
						if update_msg==0:
							updated_MWA_obsids=1
		else:
			mainlog.info('Reference time frequency calibration is done already.\n')
			msfreqs_cal.remove(float(ref_freqstamp))
			mstimes_cal.remove(float(ref_timestamp_mjd))
			spawned_ms_jobs[ref_freq_time_msname]=[ref_timestamp,0,len(touch_file_list1),160,2.0]
			ref_timechan_success=True
			break

	for j in range(len(ref_freq_time_msname)):
		index=measurement_set_list.index(ref_freq_time_msname[j])
		measurement_set_list.remove(ref_freq_time_msname[j]) # Removing ref time freq ms from ms list
		del metafits_dic[ref_freq_time_msname[j].split('.ms')[0].split('_timesliced')[0].split('_ref')[0]] # Removing ref time freq metafits from list 

	if len(measurement_set_list)>0:
		caltable_list=glob.glob(inputs.basedir+'/caltables/'+str(reftimefreq_ms_OBSID)+'/'+str(os.path.basename(ref_freq_time_msname).split('.ms')[0])+'/*ref*.cal') # Ref time caltables
		caltable_str=','.join(caltable_list)
		mainlog.info('Caltables to be applied : '+caltable_str+'\n')
		for i in range(len(measurement_set_list)):
			casa_instance=int(available_cpu_for_paircars/2)
			if mpi==1:
				casa_instance/=2
			msname=measurement_set_list[i]
			obsid=ms_OBSIDs_cal[i]
			time_diff=abs(reftimefreq_ms_OBSID-obsid)
			AMref=AccessMS(ref_freq_time_msname)
			ref_freq=AMref.calc_meanfreq()/10**6
			AMms=AccessMS(msname)
			ms_freq=AMms.calc_meanfreq()/10**6
			if time_diff>bandpass_interval*3600 and inputs.quality_factor!=0:
				mainlog.info('Performing bandpass as time interval is greater than : '+str(bandpass_interval)+' hr.\n')
				do_bandpass=True
			elif abs(ms_freq-ref_freq)>(skip_freq/10**3) and inputs.quality_factor!=0:
				mainlog.info('Performing bandpass as frequency interval is greater than : '+str(skip_freq)+' kHz.\n')
				do_bandpass=True
			else:
				mainlog.info('Do not perform bandpass, becuase either time interval is smaller than : '+str(bandpass_interval)+' or frequency is smaller than '+str(skip_freq)\
							+' kHz orquality_factor = 0.\n')
				do_bandpass=False
			metafits=metafits_dic[msname.split('.ms')[0].split('_timesliced')[0].split('_ref')[0]]
			if inputs.do_polcal==True:
				ref_time_freq_gridpoint=fits.getheader(ref_time_freq_metafits)['GRIDNUM']
				ms_gridpoint=fits.getheader(metafits)['GRIDNUM']
				if ref_time_freq_gridpoint!=ms_gridpoint:
					do_polcal=True
					mainlog.info('Beam pointing changed. Performing polarisation calibration.\n')
				elif abs(ms_freq-ref_freq)>(skip_freq_pol/10**3):
					do_polcal=True
					mainlog.info('Frequency difference is more than polarisation frequency interval. Performing polarisation calibration.\n')
				elif time_diff>(bandpass_interval*2)*3600:
					do_polcal=True
					mainlog.info('Time difference is more than polarisation time interval. Performing polarisation calibration.\n')
				else:
					mainlog.info('Beam pointing is same. Do not perform polarisation calibration.\n')
					do_polcal=False
			workdir=inputs.basedir+'/'+os.path.basename(msname).split('.ms')[0]
			if os.path.isdir(workdir)==False:
				if os.path.isfile(workdir)==True:
					os.system('rm -rf '+workdir)
				os.makedirs(workdir)
			os.system('cp -r selfcal_inputs.py '+workdir+'/selfcal_inputs.py')
			return_msg,ref_time,ref_chan,spawned_casa_instances,ref_freq_avg,ref_time_avg=\
				run_paircars_ms(msname,metafits,workdir,ref_freq_avg=ref_freq_avg,ref_time_avg=ref_time_avg,ref_time_freq=False,do_bandpass=do_bandpass,\
								do_polcal=do_polcal,num_threads=available_cpu_for_paircars,calibrator_caltable=caltable_list)	
			spawned_ms_jobs[msname]=[ref_time,ref_chan,spawned_casa_instances,ref_freq_avg,ref_time_avg]
			available_casa_instance=casa_instance-spawned_casa_instances
			basemsdir=os.path.basename(msname).split('.ms')[0]
			touch_file_list=glob.glob(inputs.basedir+'/.Finished_*cal*'+str(obsid)+'*'+basemsdir+'*')
			while True:
				if available_casa_instance>1:
					mainlog.info('At least 1 casa instance is available. Spawn new job...\n')
					break
				else:
					available_casa_instance=len(touch_file_list)-len(glob.glob(inputs.basedir+'/.Finished_*cal*'+str(obsid)+'*'+basemsdir+'*'))
					time.sleep(10.0)
					if updated_MWA_obsids==0:
						obsid_file,update_msg=update_mwa_obsids()
						if update_msg==0:
							updated_MWA_obsids=1				
	else:
		mainlog.info('Calibration jobs for all measurement sets are spawned.\n')
			
	ms_list=copy.deepcopy(ms_list_copy)
	ms_OBSIDs=copy.deepcopy(ms_OBSIDs_copy)
	metafits_dic=copy.deepcopy(metafits_dic_copy)

	# Estimating total casa instances
	#################################
	casa_instance=int(available_cpu_for_paircars/2)
	if mpi==1:
		casa_instance/=2
	open_casa_instance=0
	touch_count=0
	mainlog.info('Available cpus for P-AIRCARS: '+str(available_cpu_for_paircars)+'\n')
	mainlog.info('Total number of available CASA instances : '+str(casa_instance)+'\n')

	# Spawning jobs for making final caltables and images
	#####################################################
	mpicmd=[]
	ms_OBSIDs_cal=np.array(copy.deepcopy(ms_OBSIDs_copy_cals))
	for msname in ms_list: 
		metafits=metafits_dic[msname.split('.ms')[0].split('_timesliced')[0].split('_ref')[0]]
		OBSID=get_OBSID_from_metafits(metafits)
		if OBSID not in ms_OBSIDs_cal:
			nearest_OBSID=ms_OBSIDs_cal[np.argmin(abs(OBSID-ms_OBSIDs_cal))]
		else:
			nearest_OBSID=OBSID
		gaincal_modeldir=inputs.basedir+'/imagemodels/'+str(nearest_OBSID)
		bandpass_modeldir=inputs.basedir+'/bpimagemodels/'+str(nearest_OBSID)
		polcal_modeldir=inputs.basedir+'/polimagemodels/'+str(nearest_OBSID)
		mainlog.info('Making final calibration tables for : '+msname+'\n')
		num_jobs=spawned_ms_jobs[msname][2]
		freq_avg=spawned_ms_jobs[msname][3]
		time_avg=spawned_ms_jobs[msname][4]
		screen_name=str(OBSID)+'_screen_'+str(nearest_OBSID)+'_manage_database'
		batch_file=inputs.basedir+'/'+screen_name+'.batch'
		cmd_batch_file=inputs.basedir+'/'+screen_name+'_cmd.batch'
		cmd='manage_database --msname '+msname+' --metafits '+metafits+' --num_jobs '+str(num_jobs)+' --basedir '+basedir+\
			' --gaincal_modeldir '+gaincal_modeldir+' --bandpass_modeldir '+bandpass_modeldir+' --polcal_modeldir '+polcal_modeldir+' --localdatabase '+local_caldatabase+\
			' --freqavg '+str(freq_avg)+' --timeavg '+str(time_avg)+' --cal_obsid '+str(nearest_OBSID)+' --inputfile '+basedir+'/selfcal_inputs.py --verbose '+str(inputs.verbose)
		cmd='screen -S '+screen_name+' -X quit;sleep 2; screen -mdS '+screen_name+'; sleep 2;echo \"'+cmd+'\" > '+cmd_batch_file+';sleep 2; chmod a+rwx '+cmd_batch_file+\
			';sleep 2; screen -S '+screen_name+' -X stuff \"'+cmd_batch_file+'\\n\"; sleep 2'
		if os.path.isfile(batch_file):
			fil=open(batch_file,'r+')
		else:
			fil=open(batch_file,'w')
		fil.write(cmd)
		fil.close()
		os.system('chmod a+rwx '+batch_file)
		if mpi==1:
			screen_cmd='sh '+batch_file
			os.system('screen -S '+screen_name+' -X quit')	
			time.sleep(0.5)
			os.system('screen -mdS '+screen_name)
			time.sleep(0.5)
			mainlog.info('########################\n')
			mainlog.info('Made Screen : '+screen_name+'\n')
			mainlog.info('Command : '+cmd+'\n')
			os.system('screen -S '+screen_name+' -X stuff \"'+screen_cmd+'\n"')	
		elif mpi==0:
			mpicmd.append('-np 1 --map-by ppr:'+str(int((available_cpu_for_paircars/3)/cpu_sockets))+':core:pe=2 -x OMP_NUM_THREADS='+\
				str(available_cpu_for_paircars)+' sh '+batch_file+'\n')
			mpicmd.append('-np 1 sleep 1\n')
	basemsdir=os.path.basename(msname).split('.ms')[0]
	if updated_MWA_obsids==0:
		obsid_file,update_msg=update_mwa_obsids()
		if update_msg==0:
			updated_MWA_obsids=1
	if mpi==0:
		mpicmd_file=inputs.basedir+'/'+basemsdir+'.manage_database_mpicmd'
		if os.path.exists(mpicmd_file):
			os.system('rm -rf '+mpicmd_file)
		mpifil=open(mpicmd_file,'w')
		mainlog.info('MPI commands .....\n')
		for i in mpicmd:
			mainlog.info(i)
		mpifil.writelines(mpicmd)
		mpifil.close()
		os.system('chmod a+rwx '+mpicmd_file)
		screen_cmd='mpirun --app '+mpicmd_file
		os.system('screen -S '+basemsdir+'_manage_database -X quit')	
		time.sleep(0.5)
		os.system('screen -mdS '+basemsdir+'_manage_database')
		time.sleep(0.5)
		mainlog.info('########################\n')
		mainlog.info('Made Screen : '+basemsdir+'_manage_database\n')
		mainlog.info('Command : '+cmd+'\n')
		os.system('screen -S '+basemsdir+'_manage_database -X stuff \"'+screen_cmd+'\n"')	
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
# Applying solution to whole ms
############################### 

# Final imaging mode
####################
# TODO : include mwa_hyperbeam instead of mwa_pb
