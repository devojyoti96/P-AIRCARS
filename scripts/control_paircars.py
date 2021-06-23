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

def run_paircars_ms(msname,metafits,workdir,ref_freq_avg=0,ref_time_avg=0,ref_time_freq=False,do_bandpass=False,do_polcal=False,calatten=1.0,num_threads=0,calibrator_caltable=[]): #TODO: XY phasecal
	'''
	Function to run paircars on a measurement set
	Parameters:
	msname = Name of the measurement set
	metafits = Name of the metafits file
	workdir = Name of the working directory
	ref_freq_avg = Frequency average in kHz used in reference time. (default : 0 to take frequency resolution of the data)
	ref_time_avg = Temporal average in secind used in reference time. (default : 0 to take time resolution of the data)
	ref_time_freq = False, reference time frequency ms
	do_bandpass = False, perform bandpass or not
	do_polcal = False, perform polcal or not
	'''
	if workdir[-1]=='/':
		workdir=workdir[:-1]

	if os.path.isdir(workdir)==False:
		if os.path.exists(workdir)==True:
			os.system('rm -rf '+workdir)
		os.makedirs(workdir)

	cwd=os.getcwd()
	os.chdir(workdir)
	
	inttime=float(fits.getheader(metafits)['INTTIME'])/2
	AM=AccessMS(msname)
	antenna=AM.get_antenna_string()
	unflagchan,flagchan=flag_MWA_coarse(msname,edgewidth=160,do_flag=False,force=False)
	#mainlog.info('flagdata(vis=\''+msname+'\',mode=\'unflag\',spw=\''+unflagchan+'\',antenna=\''+antenna+'\')\n')
	#flagdata(vis=msname,mode='unflag',spw=unflagchan,antenna=antenna)
	cpu_sockets =  int(subprocess.check_output('cat /proc/cpuinfo | grep "physical id" | sort -u | wc -l', shell=True))

	open_casa_instance=0
	ms_obsid=get_OBSID_from_metafits(metafits)
	obs_atten=float(fits.getheader(metafits)['ATTEN_DB'])
	basemsdir=os.path.basename(msname).split('.ms')[0]

	last_selfcal_msg=1
	ref_time_freq_copy=copy.deepcopy(ref_time_freq)
	try:
		last_selfcal_msg,ref_time,ref_freq,ref_chan,spawned_casa_instances,ref_freq_avg,ref_time_avg=np.load(inputs.basedir+'/Ref_time_cal_record.npy')
		last_selfcal_msg=int(last_selfcal_msg)
		ref_time=str(ref_time)
		ref_chan=int(ref_chan)
		ref_freq=float(ref_freq)
		spawned_casa_instances=int(spawned_casa_instances)
		ref_freq_avg=float(ref_freq_avg)
		ref_time_avg=float(ref_time_avg)
		if last_selfcal_msg==1:
			mainlog.info('Reference time frequency calibration was failed. Try ms : '+msname+' as reference ms.\n')
			ref_time_freq=True
	except:
		ref_time_freq=copy.deepcopy(ref_time_freq_copy)
		ref_time_freq==True
		pass
	del ref_time_freq_copy
		

	if len(glob.glob(inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'*'))>0:  # Checking whether calibration already done or not 
		# Removing .Finished files if error occured
		###########################################
		touch_file_list=glob.glob(inputs.basedir+'/.Finished_*cal*'+str(ms_obsid)+'*'+basemsdir+'*')
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
		touch_file_list=glob.glob(inputs.basedir+'/.Finished_*cal*'+str(ms_obsid)+'*'+basemsdir+'*')
		spawned_jobs=int(glob.glob(inputs.basedir+'/.Finished_spawned*'+str(ms_obsid)+'*'+basemsdir+'*')[0].split('_')[-1])
		if len(touch_file_list)==spawned_jobs:
			mainlog.info('Calibration has already been done for ms : '+msname+'\n')
			mainlog.info('#########################\n')
		if os.path.exists(inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(spawned_jobs))==True:
			os.system('rm -rf '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(spawned_jobs))
			os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(len(touch_file_list)))
		else:
			os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(len(touch_file_list)))

		if last_selfcal_msg==0 and ref_time_freq==True:
			os.system('touch '+inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_0')
		elif ref_time_freq==True and os.path.exists(inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_0')==False:
			os.system('touch '+inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_error')
		try:
			return last_selfcal_msg,ref_time,ref_chan,spawned_jobs,ref_freq_avg,ref_time_avg
		except:
			return 0,0,0,spawned_jobs,160,2.0
	
	perform_leakcor=True

	mpi=MPI_check()
	# Decorrelation correction and convention correction
	####################################################

	AM=AccessMS(msname)
	output=AM.move_phasecenter_to_sun()  # Moving the phasecenter to the Sun
	mainlog.info(output)
	if inputs.do_decor_correction: # Performing decorrelation correction and IAU convention change
		mainlog.info('Performing de-correlation correction and IAU convention correction for ms : '+msname+'\n')
		decor(msname,metafits,10,False)
	else: # If user do not want decorrelation correction perform only IAU convention
		mainlog.info('Correcting to IAU convention......\n')
		AM.convert_mwa_to_iau()

	# Validating timerange
	######################
	AM=AccessMS(msname)
	unflag_times_mjd=AM.get_unflag_times_mjd()
	unflag_start_time=min(unflag_times_mjd)
	unflag_end_time=max(unflag_times_mjd)
	if inputs.timerange!='':
		try:
			ms_mjds=AM.get_timestamps_in_mjdsecs()
			ms_start_mjd=min(ms_mjds)
			ms_end_mjd=max(ms_mjds)
			time_list=[]
			for i in inputs.timerange.split(','):
				l=i.split('~')
				for j in l:
					time_list.append(j)
			time_list_mjds=sorted([float(timestamp_to_mjdsec(i,format=2)) for i in time_list])
			for i in range(len(time_list_mjds)):
				if time_list_mjds[i]<=ms_start_mjd:
					time_list_mjds[i]=ms_start_mjd
				elif time_list_mjds[i]>=ms_end_mjd:
					time_list_mjds[i]=ms_end_mjd
			time_list_mjds=sorted(time_list_mjds)
			if (max(time_list_mjds)-min(time_list_mjds))>=10:
				if min(time_list_mjds)<=unflag_start_time:
					start_time=unflag_start_time
				else:
					start_time=min(time_list_mjds)
				if max(time_list_mjds)>unflag_end_time:
					end_time=unflag_end_time
				else:
					end_time=max(time_list_mjds)
				new_timerange=mjdsec_to_timestamp(start_time,format=0)+'~'+mjdsec_to_timestamp(end_time,format=0)
			elif (max(time_list_mjds)-min(time_list_mjds))<10 and max(time_list_mjds)!=min(time_list_mjds):
				mean_mjd=(max(time_list_mjds)+min(time_list_mjds))/2.0
				if (mean_mjd-5.0)<unflag_start_time:
					start_time=unflag_start_time
				else:
					start_time=(mean_mjd-5.0)
				if (mean_mjd+5.0)>unflag_end_time:
					end_time=unflag_end_time
				else:
					end_time=(mean_mjd+5.0)
				new_timerange=mjdsec_to_timestamp(start_time,format=0)+'~'+mjdsec_to_timestamp(end_time+5.0,format=0)
			else:
				mean_mjd=max(time_list_mjds)
				if (mean_mjd-5.0)<unflag_start_time:
					start_time=unflag_start_time
				else:
					start_time=(mean_mjd-5.0)
				if (mean_mjd+5.0)>unflag_end_time:
					end_time=unflag_end_time
				else:
					end_time=(mean_mjd+5.0)
				new_timerange=mjdsec_to_timestamp(start_time,format=0)+'~'+mjdsec_to_timestamp(end_time+5.0,format=0)
		except:
			mainlog.info('Error in user given time range. Choosing full ms.\n')
			new_timerange=mjdsec_to_timestamp(unflag_start_time,format=0)+'~'+mjdsec_to_timestamp(unflag_end_time,format=0)
	else:
		new_timerange=mjdsec_to_timestamp(unflag_start_time,format=0)+'~'+mjdsec_to_timestamp(unflag_end_time,format=0)

	ms_times=AM.get_timestamps(format=0)
	ms_timerange=ms_times[0]+'~'+ms_times[-1]

	# Validating channel ranges
	###########################
	try:
		if inputs.chanrange!='':
			chans=inputs.chanrange.split(',')
			nchan=AM.get_num_channels()
			new_chan_list=[]
			for i in range(len(chans)):
				c=chans[i].split('~')
				chan_list=[]
				for chan in c:
					try:
						a=int(chan)
					except:
						mainlog.error('Channel number is wrong. Removing chanstamp '+str(i))
						break
					chan_list.append(a)
				chan_list=sorted(chan_list)
				if chan_list[0]>nchan or chan_list[1]<0:
					mainlog.info('Channel range is not in ms, removing chanstamp '+str(i))
				elif chan_list[0]<0:
					mainlog.info('Start channel is less than 0. Shifted it to channel 0.\n')
					chan_list[0]=0
				elif chan_list[1]>nchan:
					mainlog.info('End chan is greater than total number of channels. Shifted it to total number of channels.\n')
					chan_list[1]=nchan
				for x in range(chan_list[0],chan_list[1]): 
					new_chan_list.append(str(x))
			
			if len(new_chan_list)!=0:
				new_chan_list=','.join(new_chan_list)
			else:
				mainlog.info('Given channels are not in ms. Continue with all channels.\n')
				new_chan_list=','.join([str(i) for i in range(nchan)])
		else:
			nchan=AM.get_num_channels()
			new_chan_list=','.join([str(i) for i in range(nchan)])
	except:
		nchan=AM.get_num_channels()
		new_chan_list=','.join([str(i) for i in range(nchan)])

	# Auto calculate calibration parameters
	#######################################
	if inputs.calc_selfcalib_params==True:
		CP=CalcParams(msname,inputs.quality_factor,inputs.safety_factor)
		start_sigma,sigma_step,residual_frac,min_sigma,gain_minsnr,DR_delta_rms,DR_delta_neg,min_DR,max_DR,min_selfcal_snr,skip_time,skip_freq,uvrange_to_cal=CP.calc_calibration_params()
		inputs.skip_freq=skip_freq
		inputs.skip_time=skip_time
		inpfile=open(workdir+'/selfcal_inputs.py','r+')
		lines=inpfile.readlines()
		for i in range(len(lines)):
			if 'start_sigma' in lines[i]:
				lines[i]='start_sigma\t\t=\t'+str(start_sigma)+'\n'
			if 'sigma_step' in lines[i]:
				lines[i]='sigma_step\t\t=\t'+str(sigma_step)+'\n'
			if 'residual_frac' in lines[i]:
				lines[i]='residual_frac\t\t=\t'+str(residual_frac)+'\n'
			if 'min_sigma' in lines[i]:
				lines[i]='min_sigma\t\t=\t'+str(min_sigma)+'\n'
			if 'uvrange_to_cal' in lines[i]:
				lines[i]='uvrange_to_cal\t\t=\t\''+str(uvrange_to_cal)+'\'\n'	
			if 'gain_minsnr' in lines[i]:
				lines[i]='gain_minsnr\t\t=\t'+str(gain_minsnr)+'\n'
			if 'DR_delta_rms' in lines[i]:
				lines[i]='DR_delta_rms\t\t=\t'+str(DR_delta_rms)+'\n'
			if 'DR_delta_neg' in lines[i]:
				lines[i]='DR_delta_neg\t\t\t=\t'+str(DR_delta_neg)+'\n'
			if 'min_DR' in lines[i]:
				lines[i]='min_DR\t\t\t=\t'+str(min_DR)+'\n'
			if 'max_DR' in lines[i]:
				lines[i]='max_DR\t\t\t=\t'+str(max_DR)+'\n'
			if 'min_selfcal_snr' in lines[i]:
				lines[i]='min_selfcal_snr\t=\t'+str(min_selfcal_snr)+'\n'
			if 'skip_time' in lines[i]:
				lines[i]='skip_time\t\t=\t'+str(skip_time)+'\n'
			if 'skip_freq' in lines[i]:
				lines[i]='skip_freq\t\t=\t'+str(skip_freq)+'\n'
		inpfile.seek(0)
		inpfile.writelines(lines)
		inpfile.close()

	# Validating image_delta_freq and image_delta_time
	##################################################
	freqres=AM.calc_freqres()
	timeres=AM.calc_timeres()
	if inputs.image_delta_freq<freqres:
		mainlog.info('Frequency resolution of the data is less than the intended imaging frequency resolution. Setting imaging frequency resolution to frequency of the data\n')
		inputs.image_delta_freq=freqres
	if inputs.image_delta_time<timeres:
		mainlog.info('Time resolution of the data is less than the intended imaging time resolution. Setting imaging time resolution to time resolution of the data\n')
		inputs.image_delta_time=timeres
	if inputs.skip_freq<freqres:
		mainlog.info('Frequency resolution of the data is less than the skip frequency. Setting skip frequency resolution to frequency of the data.\n')
		inputs.skip_freq=freqres
	if inputs.skip_time<timeres:
		mainlog.info('Time resolution of the data is less than the skip time. Setting skip time to the time resolution of the data.\n')
		inputs.skip_time=timeres
	if inputs.image_freq<freqres:
		mainlog.info('Frequency resolution of the data is less than the image bandwidth. Setting image bandwidth resolution to frequency of the data.\n')
		inputs.image_freq=freqres
	elif inputs.image_freq>inputs.image_delta_freq:
		mainlog.info('Image bandwidth is greater than image frequency interval. Setting image bandwidth to image frequency interval.\n')
		inputs.image_freq=inputs.image_delta_freq
	if inputs.image_time>inputs.image_delta_time:
		mainlog.info('Time resolution of the image is greater than the time interval. Setting time resolution to the time interval.\n')
		inputs.image_time=inputs.image_delta_time

	mainlog.info('###############################################\n')
	mainlog.info('Skip frequency : '+str(inputs.skip_freq)+' kHz\n')
	mainlog.info('Skip time : '+str(inputs.skip_time)+' s\n')
	mainlog.info('Image frequency interval : '+str(inputs.image_delta_freq)+' kHz\n')
	mainlog.info('Image time interval : '+str(inputs.image_delta_time)+' s\n')
	mainlog.info('Image bandwidth : '+str(inputs.image_freq)+' kHz\n')
	mainlog.info('Image time resolution : '+str(inputs.image_time)+' s\n')
	mainlog.info('Channel range to image : '+str(new_chan_list)+'\n')
	mainlog.info('Time range to image : '+str(new_timerange)+'\n')
	mainlog.info('###############################################\n')

	# Flagging coarse channel edges and center
	##########################################
	mainlog.info('flag_MWA_coarse(\''+msname+'\',edgewidth=160,do_flag=True)\n')
	good_channels,channel_per_coarse=flag_MWA_coarse(msname,edgewidth=160,do_flag=True)
	unflag_channels=np.array(AM.get_unflag_chan(flagfrac=1))
	new_chan_list=np.array(new_chan_list.split(','),dtype='int')
	new_unflag_chan_list=np.intersect1d(unflag_channels,new_chan_list)
	start_chan=np.min(new_unflag_chan_list)
	end_chan=np.max(new_unflag_chan_list)

	# Spliting timeranges for reference times
	#########################################
	mainlog.info('Measurement set timeslice range : '+new_timerange+'\n')

	if new_timerange!=ms_timerange:
		mainlog.info('Spliting reference timerange : '+str(new_timerange)+'\n')	
		if os.path.isdir(workdir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftimesliced.ms')==True:
			os.system('rm -rf '+workdir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftimesliced.ms* '+\
						workdir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftimesliced.ms.flagversions')
		mainlog.info('split(vis=\''+msname+'\',outputvis=\''+workdir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftimesliced.ms\',timerange=\''+new_timerange\
					+'\',datacolumn=\'data\')\n')
		split(vis=msname,outputvis=workdir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftimesliced.ms',\
				timerange=new_timerange,datacolumn='data')
		ref_timesliced_measurement_set=workdir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftimesliced.ms'
	else:
		mainlog.info('Reference timerange is similar to MS timerange. Linking the MS.....\n')
		ref_timesliced_measurement_set=workdir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftimesliced.ms'
		if os.path.exists(ref_timesliced_measurement_set):
			if os.path.islink(ref_timesliced_measurement_set):
				os.unlink(ref_timesliced_measurement_set)
			else:
				os.system('rm -rf '+ref_timesliced_measurement_set)
		elif os.path.islink(ref_timesliced_measurement_set):
			os.unlink(ref_timesliced_measurement_set)
		os.system('ln -s '+os.path.realpath(msname)+' '+ref_timesliced_measurement_set)

	mainlog.info('Linking reference timeslice ms to timesliced ms....\n')
	timesliced_measurement_set=workdir+'/'+os.path.basename(msname).split('.ms')[0]+'_timesliced.ms'
	if os.path.exists(timesliced_measurement_set):
		if os.path.islink(timesliced_measurement_set):
			os.unlink(timesliced_measurement_set)
		else:
			os.system('rm -rf '+timesliced_measurement_set)
	elif os.path.islink(timesliced_measurement_set):
		os.unlink(timesliced_measurement_set)
	os.system('ln -s '+ref_timesliced_measurement_set+' '+timesliced_measurement_set)
	os.system('rm -rf casa*.log')	# Time slices list
	
	# Setting time and frequency averaging parameters
	#################################################
	AMtimesliced=AccessMS(timesliced_measurement_set)
	timestamps=AMtimesliced.get_timestamps()
	if freqres<160:
		freq_avg=160
	else:
		freq_avg=freqres
	if timeres<2:
		time_avg=2
	else:
		time_avg=timeres
	freq_averaging_count=0
	time_averaging_count=0
	
	touch_file_list=glob.glob(inputs.basedir+'/.Finished_*cal*'+str(ms_obsid)+'*'+basemsdir+'*')
	# Removing .Finished files if error occured
	###########################################
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

	ref_time_chan_loop_count=0 # Counter for total number of referece time channel image

	spawned_casa_instances=len(glob.glob(basedir+'/.Finished*cal*'+str(ms_obsid)+'*'+basemsdir+'*/*'+os.path.basename(msname).split('.ms')[0]+'*'))
																				 # Total casa instances already spawned or done
	cur_spawned_casa_instances=spawned_casa_instances
	pass_flag=False # Pass the while loop if True and failed to selfcal with all time and frequency averaging
	try_reduce_flag=True # Try to reduce flag solutions by increasing time averaging
	reduce_moreflag=True # Try to reduce more number of solutions flag
	reduce_flag_count=0	#
	previous_selfcal_record=''
	previous_ms=''
	previous_caltable=''
	preworkdi=''
	touch_file_list=glob.glob(inputs.basedir+'/.Finished_gcal*'+str(ms_obsid)+'*'+basemsdir+'*ref*')
	if len(touch_file_list)>0:
		touch_file=touch_file_list[0]
	else:
		touch_file=''
	ref_time=''
	ref_chan=''
	ref_freq=''
	ref_freq_avg=freq_avg
	ref_time_avg=time_avg
	# In this while loop we are checking whether the present time and frequency averaging is enough to start the self calibration. If it has it will leave the loop and go for selfcal
	##################################################################################################################################################################################
	while os.path.isfile('.ref_timechan_done')==False and ref_time_freq==True and os.path.exists(touch_file)==False:  	
		mainlog.info('Choosing averaging frequency width : '+str(ref_freq_avg)+' kHz, averaging temporal width : '+str(ref_time_avg)+' s, Skip frequency : '\
				+str(skip_freq)+' kHz, Skip time : '+str(skip_time)+' s\n')
		
		# Averaging reference timesliced measurement set
		################################################
		do_averaging=False
		AMtimesliced=AccessMS(ref_timesliced_measurement_set)
		if ref_freq_avg>AMtimesliced.calc_freqres() or ref_time_avg>AMtimesliced.calc_timeres():
			chan_width=int(ref_freq_avg/AMtimesliced.calc_freqres())
			previous_averaged_ms=glob.glob(ref_timesliced_measurement_set.split('.ms')[0]+'*s_ref_averaged.ms')
			for pams in previous_averaged_ms:
				os.system('rm -rf '+pams+'* '+pams+'.flagversions')
			ref_averaged_msname=ref_timesliced_measurement_set.split('.ms')[0]+'_'+str(float(ref_freq_avg))+'kHz_'+str(ref_time_avg)+'s_ref_averaged.ms'
			if os.path.isdir(ref_averaged_msname)==True:
				os.system('rm -rf '+ref_averaged_msname+'* '+ref_averaged_msname+'.flagversions')	
			mainlog.info('Avearging reference time measurement width frequency average :'+str(ref_freq_avg)+' kHz, temporal average :'+str(ref_time_avg)+'s\n')
			mainlog.info('split(vis=\''+ref_timesliced_measurement_set+'\',outputvis=\''+ref_averaged_msname+'\',width='+\
						str(chan_width)+',timerange=\''+new_timerange+'\',timebin=\''+str(ref_time_avg)+'s\',datacolumn=\'data\')\n')
			split(vis=ref_timesliced_measurement_set,outputvis=ref_averaged_msname,width=chan_width,timerange=new_timerange,timebin=str(ref_time_avg)+'s',datacolumn='data')
			if ref_time_freq==True:
				np.save(inputs.basedir+'/Ref_time_cal_record',np.array([0,ref_time,ref_freq,ref_chan,cur_spawned_casa_instances,ref_freq_avg,ref_time_avg]))
			do_averaging=True
		else:
			mainlog.info('No averaging is required.\n')
			previous_averaged_ms=glob.glob(ref_timesliced_measurement_set.split('.ms')[0]+'*s_ref_averaged.ms')
			for pams in previous_averaged_ms:
				os.system('rm -rf '+pams+'* '+pams+'.flagversions')
			ref_averaged_msname=ref_timesliced_measurement_set.split('.ms')[0]+'_'+str(float(ref_freq_avg))+'kHz_'+str(ref_time_avg)+'s_ref_averaged.ms'
			if os.path.exists(ref_averaged_msname):
				if os.path.islink(ref_averaged_msname):
					os.unlink(ref_averaged_msname)
				else:
					os.system('rm -rf '+ref_averaged_msname)
			elif os.path.islink(ref_averaged_msname):
				os.unlink(ref_averaged_msname)
			os.system('ln -s '+ref_timesliced_measurement_set+' '+ref_averaged_msname)
			if ref_time_freq==True:
				np.save(inputs.basedir+'/Ref_time_cal_record',np.array([0,ref_time,ref_freq,ref_chan,cur_spawned_casa_instances,ref_freq_avg,ref_time_avg]))

		AMref=AccessMS(ref_averaged_msname)
		# Making reference time and channel and time frequency grid
		##########################################################
		unflagged_channels=AMref.get_unflag_chan(flagfrac=1) # Unflagged averaged channels
		if len(unflagged_channels)==0:
			mainlog.info('No unflagged channel is present.\n')
			if ref_time_freq==True:	
				np.save(inputs.basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))	
				os.system('touch '+inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_selfcalerror')	
			if os.path.exists(inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
				os.system('rm -rf '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
			else:
				os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
			return 1,0,0,0,0,0
		if do_averaging==True or ref_time_chan_loop_count==0:
			timestamps=AMref.get_timestamps()
			total_time=AMref.calc_total_time()
			timeres=AMref.calc_timeres()
			mjdlist=AMref.get_timestamps_in_mjdsecs()
			ref_time_grid=copy.deepcopy(timestamps)
			while True:
				ref_time=ref_time_grid[int(len(ref_time_grid)/2)]
				ref_time_mjd=timestamp_to_mjdsec(ref_time)
				ref_time_mjd_frac_sec=ref_time_mjd-int(ref_time_mjd)
				frac_multiple=ref_time_mjd_frac_sec/inttime
				if (frac_multiple-int(frac_multiple))==0:
					break
				else:
					ref_time_grid.remove(ref_time)
			ref_time_grid_copy=copy.deepcopy(ref_time_grid)
		elif ref_time_chan_loop_count==0:
			timestamps=AMref.get_timestamps()
			ref_time_grid=timestamps
			ref_time_grid_copy=copy.deepcopy(ref_time_grid)
		else:
			ref_time_grid=copy.deepcopy(ref_time_grid_copy)

		ref_time_copy=copy.deepcopy(ref_time) # Copy this timestamp to remove from ref time grid if failed
		AMref=AccessMS(ref_averaged_msname)
		# Making reference time and channel and time frequency grid
		##########################################################
		unflagged_channels=AMref.get_unflag_chan(flagfrac=1) # Unflagged averaged channels
		if len(unflagged_channels)==0:
			mainlog.info('No unflagged channel is present.\n')
			if ref_time_freq==True:	
				os.system('touch '+inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_selfcalerror')
				np.save(inputs.basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))	
			if os.path.exists(inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
				os.system('rm -rf '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
			else:
				os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))	
			return 1,0,0,0,0,0

		timestamps=AMref.get_timestamps()
		total_time=AMref.calc_total_time()
		ref_freqs=AMref.get_freqs()/10**6
		if do_averaging==True or ref_time_chan_loop_count==0:		
			mjdlist=AMref.get_timestamps_in_mjdsecs()
			skip_channel=int(skip_freq/AMref.calc_freqres())
			ref_channel_grid=[]		
			if skip_channel<len(unflagged_channels):
				for i in range(0,len(unflagged_channels),skip_channel):
					ref_channel_grid.append(unflagged_channels[i])
			else:
				ref_channel_grid=unflagged_channels

			ref_chan=ref_channel_grid[int(len(ref_channel_grid)/2)]
			ref_freq=float(AMref.get_freqs()[ref_chan]/10**6)

		mainlog.info('Reference frequency : '+str(ref_freqs[ref_chan])+' MHz.\n')
		mainlog.info('Reference time : '+str(ref_time)+'\n')

		# Estimating total casa instances
		#################################
		if num_threads==0:
			total_available_cpu=psutil.cpu_count()-(psutil.cpu_count()*psutil.cpu_percent()/100.0)
			available_cpu_for_paircars=int(total_available_cpu*inputs.cpu_frac)
			total_cpu_frac=int(psutil.cpu_count()*inputs.cpu_frac)
			if available_cpu_for_paircars>total_cpu_frac:
				available_cpu_for_paircars=total_cpu_frac
		else:
			available_cpu_for_paircars=num_threads
		casa_instance=int(available_cpu_for_paircars/3)
		if mpi==1: # If mpi is not available could not restrict cpu usuage, thus using half of the available casa instances
			casa_instance/=2
		if casa_instance<1:
			casa_instance=1
		touch_count=0
		mainlog.info('Available cpus for P-AIRCARS: '+str(available_cpu_for_paircars)+'\n')
		mainlog.info('Total number of available CASA instances : '+str(casa_instance)+'\n')

		# Spliting ref time chan ms
		###########################
		ref_timechan_ms,ref_timechan_dir=spliting_timechan(ref_averaged_msname,metafits,ref_chan,ref_time,caltype='G',ref_timechan=True,\
						input_file=workdir+'/selfcal_inputs.py',datacolumn='data')
		cur_workdir=ref_timechan_dir
		mpicmd_file=inputs.basedir+'/'+basemsdir+'.ref_mpicmd'
		# Run selfcal
		while True:
			ref_time_chan_loop_count+=1
			if mpi==0: # MPI command file
				if os.path.isfile(mpicmd_file):
					os.system('rm -rf '+mpicmd_file)
			try:
				touch_file_list=glob.glob(inputs.basedir+'/.Finished_gcal*'+str(ms_obsid)+'*'+basemsdir+'*'+os.path.basename(ref_timechan_ms)+'*')
				if reduce_flag_count==1: # Restarting calibration with more time averaging
					if previous_caltable!='' and previous_record!='':
						try:
							mainlog.info('applycal(vis=\''+ref_timechan_ms+'\',gaintable=\''+previous_caltable+'\',applymode=\'calflag\',flagbackup=True)\n')
							applycal(vis=ref_timechan_ms,gaintable=previous_caltable,applymode='calflag',flagbackup=True)
							if os.path.isdir(cur_workdir+'/junk1.ms')==False:
								os.system('cp -r '+ref_timechan_ms+' '+cur_workdir+'/junk1.ms')
							if os.path.isdir(cur_workdir+'/junk1.cal')==False:
								os.system('cp -r '+previous_caltable+' '+cur_workdir+'/junk1.cal')
							if len(glob.glob(inputs.basedir+'/freq_*datetime*'))>0:
								mainlog.info('Copying previous round backup directories : '+str(glob.glob(inputs.basedir+'/freq_*datetime*'))+'\n')
								os.system('cp -r '+inputs.basedir+'/freq_*datetime* '+cur_workdir)
							num_iter,DR1,DR3,DR5,DR2,DR4,DR6,rms_list,calmode,scratch,antenna_list_index,start_sigma,antenna_added,num_ant_current_iteration,\
								num_iter_fixed_sigma,num_iter_fixed_ant,num_iteration_after_ap,stokes,phasecenter_changed,startmodel,startmask,uvsub_flag_count,\
								ra,dec,num_iter_after_phasecenter_change,phasecenter_change_done,solmode,start_time=np.load(previous_record,allow_pickle=True)
							startmodel=''
							startmask=''
							if os.path.isfile(workdir+'/Intensity_selfcal_record.npy'):
								os.system('rm -rf '+workdir+'/Intensity_selfcal_record.npy')
							selfcal_record=np.array([num_iter,DR1,DR3,DR5,DR2,DR4,DR6,rms_list,calmode,scratch,antenna_list_index,start_sigma,antenna_added,num_ant_current_iteration,\
									num_iter_fixed_sigma,num_iter_fixed_ant,num_iteration_after_ap,stokes,phasecenter_changed,startmodel,startmask,uvsub_flag_count,ra,dec,\
									num_iter_after_phasecenter_change,phasecenter_change_done,solmode,start_time],dtype='object')
							np.save(cur_workdir+'/Intensity_selfcal_record',selfcal_record)
							if prelog!='':
								mainlog.info('Copying previous log....\n')
								os.system('cp -r '+prelog+' '+cur_workdir+'/Intensity_Selfcal.log')
							if preverboselog!='':
								mainlog.info('Copying previous verbose log......\n')
								os.system('cp -r '+preverboselog+' '+cur_workdir+'/Intensity_Selfcal_verbose.log')
							if prerms!='':
								mainlog.info('Copying previous DR_rms record......\n')
								os.system('cp -r '+prerms+' '+cur_workdir+'/DR_rms.npy')
							if prerneg!='':
								mainlog.info('Copying previous DR_neg record......\n')
								os.system('cp -r '+preneg+' '+cur_workdir+'/DR_neg.npy')
							os.system('rm -rf '+previous_caltable+' '+previous_selfcal_record+' '+previous_ms+' '+prelog+' '+preverboselog+' '+prerms+' '+preneg+' '+\
										inputs.basedir+'/freq_*datetime*')
							fresh=False
						except:
							fresh=True
					else:
						fresh=True
				else:
					fresh=True	
				mainlog.info('Perforimg self-calibration using fresh = '+str(fresh)+'\n')
				if try_reduce_flag==True and reduce_flag_count<1 and reduce_moreflag==False:
					reduce_moreflag=True	
				cur_workdir=ref_timechan_dir
				if len(calibrator_caltable)!=0:
					calstring=','.join(calibrator_caltable)
					cmd='run_intensity_selfcal --msname '+ref_timechan_ms+' --metafits '+metafits+' --workdir '+cur_workdir+\
						' --dopoint True --verbose '+str(inputs.verbose)+' --interactive '+str(inputs.interactive)+' --reduce_flags '+str(reduce_moreflag)\
							+' --caltables '+calstring+' --fresh '+str(fresh)
				else:
					cmd='run_intensity_selfcal --msname '+ref_timechan_ms+' --metafits '+metafits+' --workdir '+cur_workdir+\
						' --dopoint True --verbose '+str(inputs.verbose)+' --interactive '+str(inputs.interactive)+' --reduce_flags '+str(reduce_moreflag)+' --fresh '+str(fresh)
				screen_name=str(ms_obsid)+'_'+os.path.basename(ref_timechan_ms).split('.ms')[0]+'_screen_refG'
				finished_touch_file=inputs.basedir+'/.Finished_gcal_'+str(ms_obsid)+'_'+basemsdir+'_'+os.path.basename(ref_timechan_ms)
				screen_batch_file=casa_instance_runner(cmd,screen_name,finished_touch_file)
				if mpi==0: # MPI command file
					if os.path.isfile(mpicmd_file):
						mpifil=open(mpicmd_file,'a+')
					else:
						mpifil=open(mpicmd_file,'w')
				if mpi==1: # If MPI is not available spawned screens serially
					screen_cmd='sh '+screen_batch_file
					os.system('screen -S '+screen_name+' -X quit')	
					time.sleep(0.5)
					os.system('screen -mdS '+screen_name)
					time.sleep(0.5)
					mainlog.info('########################\n')
					mainlog.info('Made Screen : '+screen_name+'\n')
					mainlog.info('Command : '+cmd+'\n')
					os.system('screen -S '+screen_name+' -X stuff \"'+screen_cmd+'\n"')	
				elif mpi==0:
					mpicmd=['-np 1 --map-by ppr:'+str(int((available_cpu_for_paircars/3)/cpu_sockets))+':core:pe=2 -x OMP_NUM_THREADS='+\
						str(available_cpu_for_paircars)+' sh '+screen_batch_file+'\n']
					mpicmd.append('-np 1 sleep 1\n')
					mainlog.info('MPI commands .....\n')
					for i in mpicmd:
						mainlog.info(i)
					mpifil.writelines(mpicmd)
					mpifil.close()
					os.system('chmod a+rwx '+mpicmd_file)
					screen_cmd='mpirun --app '+mpicmd_file
					os.system('screen -S '+str(ms_obsid)+'_'+basemsdir+'_refcal -X quit')	
					time.sleep(0.5)
					os.system('screen -mdS '+str(ms_obsid)+'_'+basemsdir+'_refcal')
					time.sleep(0.5)
					mainlog.info('########################\n')
					mainlog.info('Made Screen : '+str(ms_obsid)+'_'+basemsdir+'_refcal\n')
					mainlog.info('Command : '+cmd+'\n')
					os.system('screen -S '+str(ms_obsid)+'_'+basemsdir+'_refcal -X stuff \"'+screen_cmd+'\n"')	
					time.sleep(2.0)	
				mainlog.info('Self calibration for ms : '+ref_timechan_ms+' is spawned in screen : '+screen_name+'\n')
				mainlog.info('Waiting to finish self calibration for reference time frequency ms :'+ref_timechan_ms+'................\n') 
				while True:
					time.sleep(2)
					touch_file_list=glob.glob(inputs.basedir+'/.Finished_gcal*'+str(ms_obsid)+'*'+basemsdir+'*'+os.path.basename(ref_timechan_ms)+'*')
					if len(touch_file_list)!=0:
						msg=touch_file_list[0].split('_')[-1]
						break	
			except Exception as e: # If runtime error occured
				mainlog.error('Error occured :'+str(e)+'\n')
				mainlog.error('Error in running selfcal.\n')
				os.system('rm -rf '+basedir+'/.paircars_running')
				os.system('touch '+basedir+'/.paircars_failed')
				os.system('rm -rf '+cur_workdir)
				if ref_time_freq==True:
					os.system('touch '+inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_runtimeerror')		
					np.save(inputs.basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))
				if os.path.exists(inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
					os.system('rm -rf '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
					os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				else:
					os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				return 1,0,0,0,0,0
			if msg=='error': # If error occured in run_intensity_selfcal
				mainlog.info('Runtime error occured.\n')
				os.system('rm -rf '+cur_workdir)	
				if ref_time_freq==True:
					os.system('touch '+inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_runtimeerror')	
					np.save(inputs.basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))
				if os.path.exists(inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
					os.system('rm -rf '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
					os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				else:
					os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				return 1,0,0,0,0,0
			elif msg=='noms': # If ms is not present
				mainlog.info('Runtime error occured. No measurement set found.\n')
				os.system('rm -rf '+cur_workdir)	
				if ref_time_freq==True:
					os.system('touch '+inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_runtimeerror')	
					np.save(inputs.basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))
				if os.path.exists(inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
					os.system('rm -rf '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
					os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				else:
					os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				return 1,0,0,0,0,0
			elif msg=='nometa': # If metafits not present
				mainlog.info('Runtime error occured. No metafits file found.\n')
				os.system('rm -rf '+cur_workdir)	
				if ref_time_freq==True:
					os.system('touch '+inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_runtimeerror')	
					np.save(inputs.basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))
				if os.path.exists(inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
					os.system('rm -rf '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
					os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				else:
					os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				return 1,0,0,0,0,0
			elif msg=='moreflag' and try_reduce_flag==True: # Try to reduce more flags with more time averaging
				time.sleep(2)
				mainlog.info('More than 5 % solutions are flagged. Increasing time averaging.\n')
				new_time_avg=int(ref_time_avg+4.0) # Averaging extra 4 seconds
				if try_reduce_flag==True and reduce_flag_count<1:
					try_reduce_flag=False
					reduce_moreflag=False
					reduce_flag_count+=1
					os.system('cp -r '+cur_workdir+'/Intensity_selfcal_record.npy '+workdir+'/prerecord.npy')					
					os.system('cp -r '+cur_workdir+'/junk.precal '+workdir+'/precal.cal')
					os.system('cp -r '+cur_workdir+'/Intensity_Selfcal.log '+workdir+'/pre.log')
					os.system('cp -r '+cur_workdir+'/DR_neg.npy '+workdir+'/pre_DR_neg.npy')
					os.system('cp -r '+cur_workdir+'/DR_rms.npy '+workdir+'/pre_DR_rms.npy')
					prerms=workdir+'/pre_DR_rms.npy'
					preneg=workdir+'/pre_DR_neg.npy'
					prelog=workdir+'/pre.log'
					if os.path.exists(cur_workdir+'/Intensity_Selfcal_verbose.log')==True:
						os.system('cp -r '+cur_workdir+'/Intensity_Selfcal_verbose.log '+workdir+'/preverbose.log')
						preverboselog=workdir+'/preverbose.log'
					else:
						preverboselog=''
					os.system('rm -rf '+cur_workdir)
					previous_caltable=workdir+'/precal.cal'
					previous_record=workdir+'/prerecord.npy'
					preworkdir=cur_workdir
				for i in touch_file_list:
					msg=i.split('_')[-1]
					if type(msg)==str and msg=='moreflag':						
						os.system('rm -rf '+i)
				if new_time_avg<=10.0:
					ref_time_avg=new_time_avg
					mainlog.info('Increasing time averaging to '+str(ref_time_avg)+'s\n')
					break
			elif int(msg)>=100:
				msg=int(msg)-100
			if int(msg)==10 and pass_flag==False:  # Checking for selfcal SNR and increasing the time and frequency averaging if required
				mainlog.error('SNR for self calibration is not sufficent.\n')
				selfcal_snr=float(np.load(inputs.basedir+'/selfcal_minsnr.npy'))
				mainlog.info('Selfcal SNR : '+str(selfcal_snr)+'\n')
				new_time_avg=int(ref_time_avg*np.sqrt(selfcal_snr/min_selfcal_snr)*3)
				if new_time_avg==ref_time_avg:
					new_time_avg+=AMref.calc_timeres()
				if new_time_avg<=skip_time and new_time_avg>ref_time_avg and new_time_avg<=total_time:
					ref_time_avg=new_time_avg
					mainlog.info('Increasing time averaging to '+str(ref_time_avg)+'s\n')
					os.system('rm -rf '+touch_file_list[0])
					break
				new_freq_avg=int(freq_avg*np.sqrt(selfcal_snr/min_selfcal_snr)*3)
				if new_freq_avg<=skip_freq and new_freq_avg>freq_avg:
					mainlog.info('Time averaging reached skip time limit. Increasing frequency averaging.')
					if new_freq_avg<=skip_freq and new_freq_avg>freq_avg:
						freq_avg=new_freq_avg
						mainlog.info('Increasing frequency averaging to '+str(freq_avg)+'kHz\n')
						os.system('rm -rf '+touch_file_list[0])
						break
				else:
					if selfcal_snr>2: # If selfcal_snr is not improved after all trials but it is greater than 2, then contine, otherwise exit
						mainlog.info('Both time and frequency averaging has been tried. Still SNR is not sufficient but is greater than 2. Thus continuing with present averaging.\n')
						pass_flag=True						
						continue
					else:
						mainlog.info('Both time and frequency averaging has been tried. Still SNR is not sufficient. Trying with other time frequency.\n')
						os.system('rm -rf '+cur_workdir)
						if ref_time_freq==True:
							os.system('touch '+inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_selfcalerror')	
							np.save(inputs.basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))	
						if os.path.exists(inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
							os.system('rm -rf '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
							os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
						else:
							os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))	
						return 1,0,0,0,0,0
			elif int(msg)!=0 and int(msg)!=9 and int(msg)!=8: # If not succeeded, or max iteration reached or DR decreased but more than min DR, removing the ref time
				mainlog.info('Message : '+error_msgs(100)+' : '+error_msgs(int(msg))+'\n')
				if inputs.verbose==False:
					mainlog.info('Removing the directory : '+ref_timechan_dir)
					os.system('rm -rf '+ref_timechan_dir)
				if len(ref_time_grid)!=0:
					mainlog.info('Removing timestamp : '+str(ref_time)+' from time grid.\n')
					ref_time_grid.remove(ref_time_copy)
					ref_index=int(len(ref_time_grid)/2)
					ref_time=ref_time_grid[ref_index]
					mainlog.info('Trying for new timestamp :'+str(ref_time)+'\n')
					if inputs.verbose==False:
						os.system('rm -rf '+ref_timechan_dir)
					ref_timechan_ms,ref_timechan_dir=spliting_timechan(ref_averaged_msname,metafits,ref_chan,ref_time,caltype='G',ref_timechan=True,input_file=workfir\
									+'/selfcal_inputs.py',datacolumn='data')
					cur_workdir=ref_timechan_dir
					continue
				elif len(ref_channel_grid)!=0:
					ref_time_grid=copy.deepcopy(ref_time_grid_copy)
					ref_channel_grid.remove(ref_chan)
					ref_chan_index=int(len(ref_channel_grid)/2)
					ref_chan=ref_channel_grid[ref_chan_index]
					ref_freq=float(AMref.get_freqs()[ref_chan]/10**6)
					mainlog.info('Trying with new channel : '+str(ref_chan)+'\n')
					if inputs.verbose==False:
						os.system('rm -rf '+ref_timechan_dir)
					ref_timechan_ms,ref_timechan_dir=spliting_timechan(ref_averaged_msname,metafits,ref_chan,ref_time,caltype='G',ref_timechan=True,input_file=\
								workfir+'/selfcal_inputs.py',datacolumn='data')
					cur_workdir=ref_timechan_dir
					continue
				else:
					mainlog.info('Reference imaging has been tried over full measurement set. No good starting point is found. Exiting PAIRCARS...\n')	
					os.system('rm -rf '+basedir+'/.paircars_running')
					os.system('touch '+basedir+'/.paircars_failed')
					os.system('rm -rf '+cur_workdir)	
					if ref_time_freq==True:	
						os.system('touch '+inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_selfcalerror')		
						np.save(inputs.basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))
					if os.path.exists(inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
						os.system('rm -rf '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
						os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
					else:
						os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
					return 1,0,0,0,0,0
			elif int(msg)==0 or int(msg)==8 or int(msg)==9: # if succeeded or max iteration reached or DR decreases but more than min DR
				os.system('touch .ref_timechan_done')
				os.system('touch '+inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(msg))		
				mainlog.info('Reference time frequency calibration done.\n')
				cur_spawned_casa_instances+=1	
				if ref_time_freq==True:
					np.save(inputs.basedir+'/Ref_time_cal_record',np.array([0,ref_time,ref_freq,ref_chan,cur_spawned_casa_instances,ref_freq_avg,ref_time_avg]))
				ref_time_grid.remove(ref_time_copy)
				ref_channel_grid.remove(ref_chan)
				break

	try:
		del ref_channel_grid
		del ref_time_grid
	except:
		pass

	try:
		last_selfcal_msg,ref_time,ref_freq,ref_chan,spawned_casa_instances,ref_freq_avg,ref_time_avg=np.load(inputs.basedir+'/Ref_time_cal_record.npy')
		last_selfcal_msg=int(last_selfcal_msg)
		ref_time=str(ref_time)
		ref_chan=int(ref_chan)
		ref_freq=float(ref_freq)
		spawned_casa_instances=int(spawned_casa_instances)
		if spawned_casa_instances>cur_spawned_casa_instances and ref_time_freq==True:
			cur_spawned_casa_instances=spawned_casa_instances
		ref_freq_avg=float(ref_freq_avg)
		ref_time_avg=float(ref_time_avg)
	except:
		ref_freq_avg=freq_avg
		ref_time_avg=2
	if ref_freq_avg==0:
		ref_freq_avg=freq_avg
	if ref_time_avg==0:
		ref_time_avg=2
	# If not the reference time frequency ms, making time and channel and time frequency grid or reference time frequency ms but touch file does not exist
	######################################################################################################################################################
	if ref_time_freq==False or (len(glob.glob(inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_*'))==0 and ref_time_freq==True and ref_time_chan_loop_count==0):
		if ref_time_freq==False:
			try:
				last_selfcal_msg,ref_time,ref_freq,ref_chan,spawned_casa_instances,ref_freq_avg,ref_time_avg=np.load(inputs.basedir+'/Ref_time_cal_record.npy')
				last_selfcal_msg=int(last_selfcal_msg)
				ref_time=str(ref_time)
				ref_chan=int(ref_chan)
				ref_freq=float(ref_freq)
				spawned_casa_instances=int(spawned_casa_instances)
				if spawned_casa_instances>cur_spawned_casa_instances and ref_time_freq==True:
					cur_spawned_casa_instances=spawned_casa_instances
				ref_freq_avg=float(ref_freq_avg)
				ref_time_avg=float(ref_time_avg)
			except:
				ref_freq_avg=freq_avg
				ref_time_avg=2
			if ref_freq_avg==0:
				ref_freq_avg=freq_avg
			if ref_time_avg==0:
				ref_time_avg=2
			mainlog.info('Choosing averaging frequency width : '+str(ref_freq_avg)+' kHz, averaging temporal width : '+str(ref_time_avg)+' s, Skip frequency : '\
				+str(skip_freq)+' kHz, Skip time : '+str(skip_time)+' s\n')

			# Averaging reference measurement set
			#####################################
			AMtimesliced=AccessMS(ref_timesliced_measurement_set)
			do_averaging=False
			if ref_freq_avg>AMtimesliced.calc_freqres() or ref_time_avg>AMtimesliced.calc_timeres():
				chan_width=int(ref_freq_avg/AMtimesliced.calc_freqres())
				previous_averaged_ms=glob.glob(ref_timesliced_measurement_set.split('.ms')[0]+'*s_ref_averaged.ms')
				for pams in previous_averaged_ms:
					os.system('rm -rf '+pams+'* '+pams+'.flagversions')
				ref_averaged_msname=ref_timesliced_measurement_set.split('.ms')[0]+'_'+str(float(ref_freq_avg))+'kHz_'+str(2)+'s_ref_averaged.ms'
				if os.path.isdir(ref_averaged_msname)==True:
					os.system('rm -rf '+ref_averaged_msname+'* '+ref_averaged_msname+'.flagversions')	
				mainlog.info('Avearging reference time  measurement width frequency average :'+str(ref_freq_avg)+' kHz, temporal average :'+str(2)+'s\n')
				mainlog.info('split(vis=\''+ref_timesliced_measurement_set+'\',outputvis=\''+ref_averaged_msname+'\',width='+\
							str(chan_width)+',timerange=\''+new_timerange+'\',timebin=\''+str(ref_time_avg)+'s\',datacolumn=\'data\')\n')
				split(vis=ref_timesliced_measurement_set,outputvis=ref_averaged_msname,width=chan_width,timerange=new_timerange,timebin=str(ref_time_avg)+'s',datacolumn='data')
				do_averaging=True
			else:
				mainlog.info('No averaging is required.\n')
				previous_averaged_ms=glob.glob(ref_timesliced_measurement_set.split('.ms')[0]+'*s_ref_averaged.ms')
				for pams in previous_averaged_ms:
					os.system('rm -rf '+pams+'* '+pams+'.flagversions')
				ref_averaged_msname=ref_timesliced_measurement_set.split('.ms')[0]+'_'+str(float(ref_freq_avg))+'kHz_'+str(ref_time_avg)+'s_ref_averaged.ms'
				if os.path.exists(ref_averaged_msname):
					if os.path.islink(ref_averaged_msname):
						os.unlink(ref_averaged_msname)
					else:
						os.system('rm -rf '+ref_averaged_msname)
				elif os.path.islink(ref_averaged_msname):
					os.unlink(ref_averaged_msname)
				os.system('ln -s '+ref_timesliced_measurement_set+' '+ref_averaged_msname)
			
			AMref=AccessMS(ref_averaged_msname)
			# Making reference time and channel and time frequency grid
			##########################################################
			unflagged_channels=AMref.get_unflag_chan(flagfrac=1) # Unflagged averaged channels
			if len(unflagged_channels)==0:
				mainlog.info('No unflagged channel is present.\n')
				if ref_time_freq==True:	
					os.system('touch '+inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_selfcalerror')	
					np.save(inputs.basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))
				if os.path.exists(inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
					os.system('rm -rf '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
					os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				else:
					os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))	
				return 1,0,0,0,0,0

			timestamps=AMref.get_timestamps()
			total_time=AMref.calc_total_time()
			timeres=AMref.calc_timeres()
			mjdlist=AMref.get_timestamps_in_mjdsecs()
			ref_time_grid=copy.deepcopy(timestamps)
			while True:
				ref_time=ref_time_grid[int(len(ref_time_grid)/2)]
				ref_time_mjd=timestamp_to_mjdsec(ref_time)
				ref_time_mjd_frac_sec=ref_time_mjd-int(ref_time_mjd)
				frac_multiple=ref_time_mjd_frac_sec/inttime
				if (frac_multiple-int(frac_multiple))==0:
					break
				else:
					ref_time_grid.remove(ref_time)
			ref_time_grid_copy=copy.deepcopy(ref_time_grid)
			ref_time_copy=copy.deepcopy(ref_time) # Copy this timestamp to remove from ref time grid if failed
			try:
				del ref_time_grid
			except:
				pass			

		elif ref_time_freq==True:
			try:
				last_selfcal_msg,ref_time,ref_freq,ref_chan,spawned_casa_instances,ref_freq_avg,ref_time_avg=np.load(inputs.basedir+'/Ref_time_cal_record.npy')
				last_selfcal_msg=int(last_selfcal_msg)
				ref_time=str(ref_time)
				ref_chan=int(ref_chan)
				ref_freq=float(ref_freq)
				spawned_casa_instances=int(spawned_casa_instances)
				if spawned_casa_instances>cur_spawned_casa_instances and ref_time_freq==True:
					cur_spawned_casa_instances=spawned_casa_instances
				ref_freq_avg=float(ref_freq_avg)
				ref_time_avg=float(ref_time_avg)
			except:
				ref_freq_avg=freq_avg
				ref_time_avg=2
			if ref_freq_avg==0:
				ref_freq_avg=freq_avg
			if ref_time_avg==0:
				ref_time_avg=2

			previous_averaged_ms=glob.glob(ref_timesliced_measurement_set.split('.ms')[0]+'*s_ref_averaged.ms')
			for pams in previous_averaged_ms:
				os.system('rm -rf '+pams+'* '+pams+'.flagversions')
			ref_averaged_msname_list=glob.glob(ref_timesliced_measurement_set.split('.ms')[0]+'_*s_ref_averaged.ms')
			if len(ref_averaged_msname_list)==0:
				AMtimesliced=AccessMS(ref_timesliced_measurement_set)
				if ref_freq_avg>AMtimesliced.calc_freqres() or ref_time_avg>AMtimesliced.calc_timeres():
					chan_width=int(ref_freq_avg/AMtimesliced.calc_freqres())
					previous_averaged_ms=glob.glob(ref_timesliced_measurement_set.split('.ms')[0]+'*s_ref_averaged.ms')
					for pams in previous_averaged_ms:
						os.system('rm -rf '+pams+'* '+pams+'.flagversions')
					ref_averaged_msname=ref_timesliced_measurement_set.split('.ms')[0]+'_'+str(float(ref_freq_avg))+'kHz_'+str(ref_time_avg)+'s_ref_averaged.ms'
					if os.path.isdir(ref_averaged_msname)==True:
						os.system('rm -rf '+ref_averaged_msname+'* '+ref_averaged_msname+'.flagversions')	
					mainlog.info('Avearging reference time measurement width frequency average :'+str(ref_freq_avg)+' kHz, temporal average :'+str(ref_time_avg)+'s\n')
					mainlog.info('split(vis=\''+ref_timesliced_measurement_set+'\',outputvis=\''+ref_averaged_msname+'\',width='+\
								str(chan_width)+',timerange=\''+new_timerange+'\',timebin=\''+str(ref_time_avg)+'s\',datacolumn=\'data\')\n')
					split(vis=ref_timesliced_measurement_set,outputvis=ref_averaged_msname,width=chan_width,timerange=new_timerange,timebin=str(ref_time_avg)+'s',datacolumn='data')
				else:
					mainlog.info('No averaging is required.\n')
					previous_averaged_ms=glob.glob(ref_timesliced_measurement_set.split('.ms')[0]+'*s_ref_averaged.ms')
					for pams in previous_averaged_ms:
						os.system('rm -rf '+pams+'* '+pams+'.flagversions')
					ref_averaged_msname=ref_timesliced_measurement_set.split('.ms')[0]+'_'+str(float(ref_freq_avg))+'kHz_'+str(ref_time_avg)+'s_ref_averaged.ms'
					if os.path.exists(ref_averaged_msname):
						if os.path.islink(ref_averaged_msname):
							os.unlink(ref_averaged_msname)
						else:
							os.system('rm -rf '+ref_averaged_msname)
					elif os.path.islink(ref_averaged_msname):
						os.unlink(ref_averaged_msname)
					os.system('ln -s '+ref_timesliced_measurement_set+' '+ref_averaged_msname)
			else:
				ref_averaged_msname=ref_averaged_msname_list[0]

	# Averaging timesliced measurement set
	######################################
	AM=AccessMS(timesliced_measurement_set)	
	try:
		last_selfcal_msg,ref_time,ref_freq,ref_chan,spawned_casa_instances,ref_freq_avg,ref_time_avg=np.load(inputs.basedir+'/Ref_time_cal_record.npy')
		last_selfcal_msg=int(last_selfcal_msg)
		ref_time=str(ref_time)
		ref_chan=int(ref_chan)
		ref_freq=float(ref_freq)
		spawned_casa_instances=int(spawned_casa_instances)
		if spawned_casa_instances>cur_spawned_casa_instances and ref_time_freq==True:
			cur_spawned_casa_instances=spawned_casa_instances
		ref_freq_avg=float(ref_freq_avg)
		ref_time_avg=float(ref_time_avg)
	except:
		ref_freq_avg=freq_avg
		ref_time_avg=2
	if ref_freq_avg==0:
		ref_freq_avg=freq_avg
	if ref_time_avg==0:
		ref_time_avg=2


	if ref_freq_avg>AM.calc_freqres() or ref_time_avg>AM.calc_timeres():
		mainlog.info('Linking averaged reference timesliced ms to timesliced ms.....\n')
		averaged_msname=timesliced_measurement_set.split('.ms')[0]+'_'+str(float(ref_freq_avg))+'kHz_'+str(ref_time_avg)+'s_averaged.ms'
		if os.path.exists(averaged_msname):
			if os.path.islink(averaged_msname):
				os.unlink(averaged_msname)
			else:
				os.system('rm -rf '+averaged_msname)
		elif os.path.islink(averaged_msname):
			os.unlink(averaged_msname)
		os.system('ln -s '+ref_averaged_msname+' '+averaged_msname)
	else:
		mainlog.info('No averaging is required.\n')
		averaged_msname=timesliced_measurement_set.split('.ms')[0]+'_'+str(float(ref_freq_avg))+'kHz_'+str(ref_time_avg)+'s_averaged.ms'
		if os.path.exists(averaged_msname):
			if os.path.islink(averaged_msname):
				os.unlink(averaged_msname)
			else:
				os.system('rm -rf '+averaged_msname)
		elif os.path.islink(averaged_msname):
			os.unlink(averaged_msname)
		os.system('ln -s '+timesliced_measurement_set+' '+averaged_msname)


	AMref=AccessMS(ref_averaged_msname)
	# Making reference time and channel and time frequency grid
	##########################################################
	unflagged_channels=AMref.get_unflag_chan(flagfrac=1) # Unflagged averaged channels
	if len(unflagged_channels)==0:
		mainlog.info('No unflagged channel is present.\n')
		if ref_time_freq==True:
			os.system('touch '+inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_selfcalerror')		
			np.save(inputs.basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))	
		if os.path.exists(inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
			os.system('rm -rf '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
			os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
		else:
			os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
		return 1,0,0,0,0,0

	timestamps=AMref.get_timestamps()
	total_time=AMref.calc_total_time()
	ref_freqs=AMref.get_freqs()/10**6
	unflagged_channels=AMref.get_unflag_chan(flagfrac=1) # Unflagged averaged channels
	if len(unflagged_channels)==0:
		mainlog.info('No unflagged channel is present.\n')
		if ref_time_freq==True:
			os.system('touch '+inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_selfcalerror')		
			np.save(inputs.basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))	
		if os.path.exists(inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
			os.system('rm -rf '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
			os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
		else:
			os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
		return 1,0,0,0,0,0

	mainlog.info('Reference frequency : '+str(ref_freqs[ref_chan])+' MHz.\n')
	mainlog.info('Reference time : '+str(ref_time)+'\n')
	try:
		del ref_channel_grid
		del ref_time_grid
	except:
		pass

	AM=AccessMS(averaged_msname)
	unflagged_channels=AM.get_unflag_chan(flagfrac=1) # Unflagged averaged channels
	if len(unflagged_channels)==0:	
		mainlog.info('No unflagged channel is present.\n')
		if ref_freq_time==True:
			os.system('touch '+inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_selfcalerror')		
			np.save(inputs.basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))	
		if os.path.exists(inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
			os.system('rm -rf '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
			os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
		else:
			os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
		return 1,0,0,0,0,0
				
	timestamps=AM.get_timestamps()
	mjd_timestamps=AM.get_timestamps_in_mjdsecs()

	skip_channel=int(skip_freq/AM.calc_freqres())
	skip_timestamp=int(skip_time/AM.calc_timeres())

	channel_grid=[]
	time_grid=[]

	if skip_channel<len(unflagged_channels): # Final channel grid
		for i in range(min(unflagged_channels),max(unflagged_channels),skip_channel):
			if (i!=ref_chan and ref_time_freq==True) or ref_time_freq==False:
				channel_grid.append(i)
	else:
		channel_grid=unflagged_channels

	if skip_timestamp<len(timestamps): # Final time grid
		for i in range(0,len(mjd_timestamps),skip_timestamp):
			if (timestamps[i]!=ref_time and ref_time_freq==True) or ref_time_freq==False:
				time_grid.append(timestamps[i])
	else:
		time_grid=timestamps

	if ref_chan=='' or ref_time_freq==False:
		ref_chan=channel_grid[int(len(channel_grid)/2)]
		ref_freq=float(AM.get_freqs()[ref_chan]/10**6)
	if ref_time=='' or ref_time_freq==False:
		ref_time=time_grid[int(len(time_grid)/2)]

	if ref_time_freq==True:
		mjd_time_grid=np.array([timestamp_to_mjdsec(i,format=0) for i in time_grid])
		ref_mjd=timestamp_to_mjdsec(ref_time,format=0)
		time_diff_list=np.abs(np.array(mjd_time_grid)-ref_mjd)
		pos=np.where(time_diff_list<=skip_time)
		mjd_time_grid=np.delete(mjd_time_grid,pos)
		time_grid=[mjdsec_to_timestamp(i,format=0) for i in mjd_time_grid]

	time_grid_copy=copy.deepcopy(time_grid)

	mainlog.info('Channel grid list : '+str(channel_grid)+'\n')
	mainlog.info('Timestamp grid list : '+str(time_grid)+'\n')
	mainlog.info('Reference time : '+str(ref_time)+'\n')
	mainlog.info('Reference frequency : '+str(ref_freqs[ref_chan])+' MHz.\n')

	# Making ref time freq gaintable list
	#####################################
	if ref_time_freq==True:
		ref_timechan_caltable=glob.glob(inputs.basedir+'/caltables/'+str(ms_obsid)+'/'+basemsdir+'/*ref*.cal')
		ref_gaintable=copy.deepcopy(calibrator_caltable)
		ref_gaintable+=ref_timechan_caltable
	else:
		ref_gaintable=copy.deepcopy(calibrator_caltable)

	# Estimating total casa instances
	#################################
	if num_threads==0:
		total_available_cpu=psutil.cpu_count()-(psutil.cpu_count()*psutil.cpu_percent()/100.0)
		available_cpu_for_paircars=int(total_available_cpu*inputs.cpu_frac)
		total_cpu_frac=int(psutil.cpu_count()*inputs.cpu_frac)
		if available_cpu_for_paircars>total_cpu_frac:
			available_cpu_for_paircars=total_cpu_frac
	else:
		available_cpu_for_paircars=num_threads
	casa_instance=int(available_cpu_for_paircars/2)
	if mpi==1:
		casa_instance/=2
	if casa_instance<1:
		casa_instance=1
	touch_count=0
	mainlog.info('Available cpus for P-AIRCARS: '+str(available_cpu_for_paircars)+'\n')
	mainlog.info('Total number of available CASA instances : '+str(casa_instance)+'\n')

	# Spliting gaincal measurement set
	##################################
	gaincal_cmd_list=[]
	gaincal_screen_list=[]
	gaincal_finished_file_list=[]
	batch_file_list=[]
	calstring=','.join(ref_gaintable)
	temp_ref_ms=''
	if ref_time_freq==True:
		ref_touch_list=glob.glob(inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_*')
		if len(ref_touch_list)==0:
			os.system('touch '+inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_0')		
	if len(time_grid)!=0:
		mainlog.info('Spliting reference chan data for performing gaincal.\n')
		for timestamp in time_grid:
			splited_msname,splited_msdir=spliting_timechan(averaged_msname,metafits,str(ref_chan),timestamp,caltype='G',ref_timechan=False,\
											input_file=workdir+'/selfcal_inputs.py',datacolumn='data')
			if timestamp==ref_time:
				temp_ref_ms=splited_msname
			touch_file=glob.glob(inputs.basedir+'/.Finished_gcal_'+str(ms_obsid)+'_'+basemsdir+'_'+os.path.basename(splited_msname)+'*')
			if len(touch_file)==0:
				cmd='run_intensity_selfcal --msname '+splited_msname+' --metafits '+metafits+' --workdir '+splited_msdir+\
					' --dopoint True --verbose '+str(inputs.verbose)+' --interactive '+str(inputs.interactive)+' --fresh True --reduce_flags True --caltables '+calstring
				gaincal_cmd_list.append(cmd)
				gaincal_screen_list.append(str(ms_obsid)+'_'+os.path.basename(splited_msname).split('.ms')[0]+'_screen_G')
				gaincal_finished_file_list.append(inputs.basedir+'/.Finished_gcal_'+str(ms_obsid)+'_'+basemsdir+'_'+os.path.basename(splited_msname))
		mpicmd=[]
		mpicount=0
		while len(gaincal_cmd_list)!=0:  # Loop while all gaincal cmds are spawned
			touch_file_list=glob.glob(inputs.basedir+'/.Finished_*cal*'+str(ms_obsid)+'*'+basemsdir+'*')
			if (len(touch_file_list)-touch_count)>0:
				if open_casa_instance>1:
					open_casa_instance-=(len(touch_file_list)-touch_count)
					mainlog.info('New CASA instance available : '+str((len(touch_file_list)-touch_count))+'\n')
			while len(gaincal_screen_list)!=0 and (casa_instance-open_casa_instance)>=1:
				screen_name=gaincal_screen_list[0]
				cmd=gaincal_cmd_list[0]
				finished_file=gaincal_finished_file_list[0]
				screen_batch_file=casa_instance_runner(cmd,screen_name,finished_file)
				batch_file_list.append(screen_batch_file)
				touch_count+=1
				if mpi==1:
					screen_cmd='sh '+screen_batch_file
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
							str(available_cpu_for_paircars)+' sh '+screen_batch_file+'\n')
					mpicmd.append('-np 1 sleep 1\n')
				open_casa_instance+=1
				cur_spawned_casa_instances+=1
				gaincal_screen_list.remove(screen_name)
				gaincal_cmd_list.remove(cmd)
				gaincal_finished_file_list.remove(finished_file)
				if open_casa_instance>=casa_instance or len(gaincal_screen_list)==0:
					mainlog.info('Maximum casa instances spawned. Waiting for complete those jobs.\n')
					if mpi==0:
						mpicmd_file=inputs.basedir+'/'+basemsdir+'.gcal_mpicmd_'+str(mpicount)
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
						os.system('screen -S '+str(ms_obsid)+'_'+basemsdir+'_gcal_'+str(mpicount)+' -X quit')	
						time.sleep(0.5)
						os.system('screen -mdS '+str(ms_obsid)+'_'+basemsdir+'_gcal_'+str(mpicount))
						time.sleep(0.5)
						mainlog.info('########################\n')
						mainlog.info('Made Screen : '+str(ms_obsid)+'_'+basemsdir+'_gcal_'+str(mpicount)+'\n')
						mainlog.info('Command : '+screen_cmd+'\n')
						os.system('screen -S '+str(ms_obsid)+'_'+basemsdir+'_gcal_'+str(mpicount)+' -X stuff \"'+screen_cmd+'\n"')
						if casa_instance*2.0>10:
							sleep_time=casa_instance*2.0
						else:
							sleep_time=10	
						time.sleep(float(sleep_time))
						mpicount+=1
						mpicmd=[]
					break
			time.sleep(2.0)			
		mainlog.info('All gaincal jobs are spawned.\n')
	else:
		mainlog.info('No timestamp is left for calibration.\n')

	# Waiting for gaincal to finish # TODO : modify this part
	######################################
	if ref_time_freq==False:
		while True:
			gtables=glob.glob(inputs.basedir+'/caltables/'+str(ms_obsid)+'/'+basemsdir+'/*'+os.path.basename(temp_ref_ms).split('.ms')[0]+'*.cal')
			if len(gtables)>0:
				mainlog.info('Gaintable found : '+gtables[0]+'\n')
				break
			else:
				time.sleep(2.0)
		ref_gaintable=[gtables[0]]
		
	# Applying ref time solution
	############################
	mainlog.info('Applying gain solution for reference time in all times and all channels.......\n')
	mainlog.info('applycal(vis=\''+ref_averaged_msname+'\',gaintable='+str(ref_gaintable)+',applymode=\'calflag\',flagbackup=True)\n')
	applycal(vis=ref_averaged_msname,gaintable=ref_gaintable,applymode='calflag',flagbackup=True)
	mainlog.info('applycal(vis=\''+ref_averaged_msname+'\',gaintable='+str(ref_gaintable)+',applymode=\'calflag\',flagbackup=True)\n')
	applycal(vis=ref_averaged_msname,gaintable=ref_gaintable,applymode='calflag',flagbackup=True)
	flaglist=flagmanager(vis=ref_averaged_msname,mode='list')
	flaglist_keys=list(flaglist.keys())
	flaglist_keys.remove('MS')
	if len(flaglist_keys)>0:
		for i in flaglist_keys:
			last_flagversion=flaglist[i]['name']
			# Restore the flag and delete the present flag version
			mainlog.info('flagmanager(vis=\''+ref_averaged_msname+'\',mode=\'restore\',versionname=\''+str(last_flagversion)+'\',merge=\'replace\')\n')
			mainlog.info('flagmanager(vis=\''+ref_averaged_msname+'\',mode=\'delete\',versionname=\''+str(last_flagversion)+'\')\n')
			flagmanager(vis=ref_averaged_msname,mode='restore',versionname=last_flagversion,merge='replace')
			flagmanager(vis=ref_averaged_msname,mode='delete',versionname=last_flagversion)

	#Deciding bandpass selfcal conditions
	#####################################
	AM=AccessMS(ref_averaged_msname)
	unflagged_channels=AM.get_unflag_chan(flagfrac=1) # Unflagged averaged channels
	if len(unflagged_channels)<=1:
		mainlog.info('Only 1 unflagged channel is present. No bandpass is required.\n')
		do_bandpass==False
	elif inputs.quality_factor==0:
		mainlog.info('Quality factor is 0. Skipping bandpass self calibration.\n')
		do_bandpass==False
	elif do_bandpass==True:
		mainlog.info('Proceed for bandpass self-calibration considering same source model for '+str(skip_freq)+' kHz\n')
		if inputs.interactive==True:
			want_change=input('Want to change bandpass bandwidth? If yes type frequency bandwidth in kHz or press enter\n')
			if want_change!='':
				skip_freq=float(want_change)
				mainlog.info('Now proceed for bandpass self-calibration considering same source model for modified bandwidth '+str(skip_freq)+' kHz\n')

	# Performing bandpass selfcal
	#############################
	num_bp=0
	if do_bandpass==True:
		if ref_time_freq==True:
			ref_touch_list=glob.glob(inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_*')
			if len(ref_touch_list)==0:
				os.system('touch '+inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_0')

		AM=AccessMS(ref_averaged_msname)
		nchan_per_bandpass=int(skip_freq/AM.calc_freqres())
		nchan=AM.get_num_channels()
		bandpass_cmd_list=[]
		bandpass_screen_list=[]
		bandpass_finished_file_list=[]
		mainlog.info('Spliting reference time data for performing bandpass........\n')
		for i in range(0,nchan,nchan_per_bandpass):
			start_chan=i
			end_chan=i+nchan_per_bandpass-1
			if end_chan>nchan:
				end_chan=nchan-1
			mainlog.info('Spliting ms of channel range : '+str(start_chan)+'~'+str(end_chan)+'\n')
			if ref_time_freq==True:
				splited_msname,splited_msdir=spliting_timechan(ref_averaged_msname,metafits,str(start_chan)+'~'+str(end_chan),ref_time,caltype='B',\
						ref_timechan=True,input_file=workdir+'/selfcal_inputs.py',datacolumn='corrected')
			else:
				splited_msname,splited_msdir=spliting_timechan(ref_averaged_msname,metafits,str(start_chan)+'~'+str(end_chan),ref_time,caltype='B',\
						ref_timechan=False,input_file=workdir+'/selfcal_inputs.py',datacolumn='corrected')
			touch_file=glob.glob(inputs.basedir+'/.Finished_bcal_'+str(ms_obsid)+'_'+basemsdir+'_'+os.path.basename(splited_msname)+'*')
			if len(touch_file)==0:
				cmd='run_bandpass_selfcal --msname '+splited_msname+' --metafits '+metafits+' --workdir '+splited_msdir+\
				' --verbose '+str(inputs.verbose)+' --interactive '+str(inputs.interactive)+' --fresh True'
				bandpass_cmd_list.append(cmd)
				bandpass_screen_list.append(str(ms_obsid)+'_'+os.path.basename(splited_msname).split('.ms')[0]+'_screen_B')
				bandpass_finished_file_list.append(inputs.basedir+'/.Finished_bcal_'+str(ms_obsid)+'_'+basemsdir+'_'+os.path.basename(splited_msname))
			if end_chan>=nchan:
				break
		finished_bandpass=False
		num_bp=len(bandpass_screen_list)
		bp_finish_list=glob.glob(inputs.basedir+'/.Finished_*bcal*'+str(ms_obsid)+'*'+basemsdir+'*')
		if len(bp_finish_list)==num_bp:
			mainlog.info('Bandpass for channel blocks are finished.\n')
			bandpass_screen_list=[]
			bandpass_cmd_list=[]
			bandpass_finished_file_list=[]
			finished_bandpass=True
		else:
			mainlog.info('Waiting for available casa instance.....\n')
		mpicmd=[]	
		mpicount=0
		while finished_bandpass==False:
			touch_file_list=glob.glob(inputs.basedir+'/.Finished_*cal*'+str(ms_obsid)+'*'+basemsdir+'*')
			if (len(touch_file_list)-touch_count)>0:
				if open_casa_instance>1:
					open_casa_instance-=(len(touch_file_list)-touch_count)
					mainlog.info('New CASA instance available : '+str((len(touch_file_list)-touch_count))+'\n')		
			while len(bandpass_screen_list)!=0 and (casa_instance-open_casa_instance)>=1:

				screen_name=bandpass_screen_list[0]
				cmd=bandpass_cmd_list[0]
				finished_file=bandpass_finished_file_list[0]
				screen_batch_file=casa_instance_runner(cmd,screen_name,finished_file)
				touch_count+=1
				if mpi==1:
					screen_cmd='sh '+screen_batch_file
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
							str(available_cpu_for_paircars)+' sh '+screen_batch_file+'\n')
					mpicmd.append('-np 1 sleep 1\n')
				open_casa_instance+=1
				cur_spawned_casa_instances+=1
				bandpass_screen_list.remove(screen_name)
				bandpass_cmd_list.remove(cmd)
				bandpass_finished_file_list.remove(finished_file)
				if open_casa_instance>=casa_instance or len(bandpass_screen_list)==0:
					mainlog.info('Maximum casa instances spawned. Waiting for complete those jobs.\n')
					if mpi==0:
						mpicmd_file=inputs.basedir+'/'+basemsdir+'.bcal_mpicmd_'+str(mpicount)
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
						os.system('screen -S '+str(ms_obsid)+'_'+basemsdir+'_bcal='+str(mpicount)+' -X quit')	
						time.sleep(0.5)
						os.system('screen -mdS '+str(ms_obsid)+'_'+basemsdir+'_bcal_'+str(mpicount))
						time.sleep(0.5)
						mainlog.info('########################\n')
						mainlog.info('Made Screen : '+str(ms_obsid)+'_'+basemsdir+'_bcal_'+str(mpicount)+'\n')
						mainlog.info('Command : '+screen_cmd+'\n')
						os.system('screen -S '+str(ms_obsid)+'_'+basemsdir+'_bcal_'+str(mpicount)+' -X stuff \"'+screen_cmd+'\n"')	
						if casa_instance*2.0>10:
							sleep_time=casa_instance*2.0
						else:
							sleep_time=10	
						time.sleep(float(sleep_time))
						mpicount+=1
						mpicmd=[]
					break
			time.sleep(2.0)
			bp_finish_list=glob.glob(inputs.basedir+'/.Finished_*bcal*'+str(ms_obsid)+'*'+basemsdir+'*')
			if len(bandpass_cmd_list)==0:
				mainlog.info('All bandpass tasks are spawned for all spectral slices.\n')
				bandpass_screen_list=[]
				bandpass_cmd_list=[]
				bandpass_finished_file_list=[]
				finished_bandpass=True
				break
	
	# Waiting for bandpass bash script
	##################################
	prefix_cmds=[]
	prefix_cmds.append('#!/bin/bash\n')
	prefix_cmds.append('shopt -s nullglob\n')
	prefix_cmds.append('logfiles=('+inputs.basedir+'/.Finished_bcal*'+str(ms_obsid)+'*'+str(basemsdir)+'*)\n')
	prefix_cmds.append('C=${#logfiles[@]}\n')
	prefix_cmds.append('\techo "Waiting : "$C\n')
	prefix_cmds.append('until [ $C -ge '+str(num_bp)+' ]\n')
	prefix_cmds.append('do\n') 
	prefix_cmds.append('\tsleep 5\n')
	prefix_cmds.append('\tshopt -s nullglob\n')
	prefix_cmds.append('\tlogfiles=('+inputs.basedir+'/.Finished_bcal*'+str(ms_obsid)+'*'+str(basemsdir)+'*)\n')
	prefix_cmds.append('\tC=${#logfiles[@]}\n')
	prefix_cmds.append('done\n')
	prefix_cmds.append('shopt -s nullglob\n')
	prefix_cmds.append('logfiles=('+inputs.basedir+'/bpcaltables/'+str(ms_obsid)+'/'+basemsdir+'/*.bcal)\n')
	prefix_cmds.append('x=\'\'\n')
	prefix_cmds.append('for i in ${logfiles[@]}\n')
	prefix_cmds.append('do\n')
	prefix_cmds.append('x=$x$i\',\'\n')
	prefix_cmds.append('done\n')
	prefix_cmds.append('y=${x%?}\n')
	
	# Spliting gain calibrated reference time channnel measurement set for polarisation calibration
	###############################################################################################
	if do_polcal==True:# and ((do_bandpass==True and finished_bandpass==True) or do_bandpass==False):
		if ref_time_freq==True:
			ref_touch_list=glob.glob(inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_*')
			if len(ref_touch_list)==0:
				os.system('touch '+inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_0')
	
		AM=AccessMS(ref_averaged_msname)
		unflagged_channels=AM.get_unflag_chan(flagfrac=1) # Unflagged averaged channels
		if len(unflagged_channels)==0:
			mainlog.info('No unflagged channel is present.\n')
			if os.path.exists(inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
				os.system('rm -rf '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
			else:
				os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
			if ref_time_freq==True:
				np.save(inputs.basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))	
				os.system('touch '+inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_selfcalerror')	
			return 1,0,0,0,0,0
				
		# Deciding polcal bandwidth
		###########################
		if inputs.quality_factor==0:
			skip_freq_pol=2560
		elif inputs.quality_factor==1:
			skip_freq_pol=1280
		elif inputs.quality_factor==2:
			skip_freq_pol=1280

		skip_channel_pol=int(skip_freq_pol/AM.calc_freqres())
		pol_channel_grid=[]
		for i in range(min(unflagged_channels),max(unflagged_channels),skip_channel_pol):
			pol_channel_grid.append(i)
		
		mainlog.info('Polcal channel grid list : '+str(pol_channel_grid)+' for calibration per '+str(skip_freq_pol)+' MHz.\n')
		polcal_cmd_list=[]
		polcal_screen_list=[]
		polcal_finished_file_list=[]
		for i in pol_channel_grid:
			mainlog.info('Spliting data for performing polarisation calibration of channel : '+str(i)+' and timerange : '+ref_time+'\n')
			if ref_time_freq==True:
				splited_msname,splited_msdir=spliting_timechan(ref_averaged_msname,metafits,str(i),ref_time,caltype='P',ref_timechan=True,\
											input_file=workdir+'/selfcal_inputs.py',datacolumn='data')
			else:
				splited_msname,splited_msdir=spliting_timechan(ref_averaged_msname,metafits,str(i),ref_time,caltype='P',ref_timechan=False,\
											input_file=workdir+'/selfcal_inputs.py',datacolumn='data')
			calstring=','.join(ref_gaintable)
			touch_file=glob.glob(inputs.basedir+'/.Finished_pcal_'+str(ms_obsid)+'_'+basemsdir+'_'+os.path.basename(splited_msname)+'*')
			if len(touch_file)==0:
				if len(ref_gaintable)!=0:
					cmd='run_pol_selfcal --msname '+splited_msname+' --metafits '+metafits+' --workdir '+splited_msdir+' --verbose '+str(inputs.verbose)+\
					' --interactive '+str(inputs.interactive)+' --fresh True --gaincal '+str(perform_leakcor)+' --caltables '+calstring+',\"$y\"'
					polcal_cmd_list.append(cmd)
					polcal_screen_list.append(str(ms_obsid)+'_'+os.path.basename(splited_msname).split('.ms')[0]+'_screen_P')
					polcal_finished_file_list.append(inputs.basedir+'/.Finished_pcal_'+str(ms_obsid)+'_'+basemsdir+'_'+os.path.basename(splited_msname))
		num_pol=len(polcal_screen_list)
		finished_polcal=False
		polcal_finish_list=glob.glob(inputs.basedir+'/.Finished_*pcal*'+str(ms_obsid)+'*'+basemsdir+'*')
		if len(polcal_finish_list)==num_pol:
			mainlog.info('Polcal for all coarse channels have been finished.\b')
			polcal_screen_list=[]
			polcal_cmd_list=[]
			polcal_finished_file_list=[]
			finished_polcal=True
		mpicmd_file=inputs.basedir+'/'+basemsdir+'.polcal_mpicmd'
		mainlog.info('Waiting for available casa instance......\n')
		mpicmd=[]	
		mpicount=0
		while finished_polcal==False:	
			touch_file_list=glob.glob(inputs.basedir+'/.Finished_*cal*'+str(ms_obsid)+'*'+basemsdir+'*')
			if (len(touch_file_list)-touch_count)>0:
				if open_casa_instance>1:
					open_casa_instance-=(len(touch_file_list)-touch_count)
					mainlog.info('New CASA instance available : '+str((len(touch_file_list)-touch_count))+'\n')
			while len(polcal_screen_list)!=0 and (casa_instance-open_casa_instance)>=1:
				screen_name=polcal_screen_list[0]
				cmd=polcal_cmd_list[0]
				finished_file=polcal_finished_file_list[0]
				screen_batch_file=casa_instance_runner(cmd,screen_name,finished_file,prefix_cmds=prefix_cmds)
				touch_count+=1
				if mpi==1:
					screen_cmd='sh '+screen_batch_file
					os.system('screen -S '+screen_name+' -X quit')	
					time.sleep(0.5)
					os.system('screen -mdS '+screen_name)
					time.sleep(0.5)
					mainlog.info('########################\n')
					mainlog.info('Made Screen : '+screen_name+'\n')
					mainlog.info('Command : '+cmd+'\n')
					os.system('screen -S '+screen_name+' -X stuff \"'+screen_cmd+'\n"')	
				elif mpi==0:
					mpicmd.append('-np 1 --map-by ppr:'+str(int((available_cpu_for_paircars/4)/cpu_sockets))+':core:pe=4 -x OMP_NUM_THREADS='+\
							str(available_cpu_for_paircars)+' sh '+screen_batch_file+'\n')
					mpicmd.append('-np 1 sleep 1\n')
				open_casa_instance+=1
				cur_spawned_casa_instances+=1
				polcal_screen_list.remove(screen_name)
				polcal_cmd_list.remove(cmd)
				polcal_finished_file_list.remove(finished_file)
				if open_casa_instance>=casa_instance or len(polcal_screen_list)==0:
					mainlog.info('Maximum casa instances spawned. Waiting for complete those jobs.\n')	
					if mpi==0:
						mpicmd_file=inputs.basedir+'/'+basemsdir+'.pcal_mpicmd_'+str(mpicount)
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
						os.system('screen -S '+str(ms_obsid)+'_'+basemsdir+'_polcal_'+str(mpicount)+' -X quit')	
						time.sleep(0.5)
						os.system('screen -mdS '+str(ms_obsid)+'_'+basemsdir+'_polcal_'+str(mpicount))
						time.sleep(0.5)
						mainlog.info('########################\n')
						mainlog.info('Made Screen : '+str(ms_obsid)+'_'+basemsdir+'_polcal_'+str(mpicount)+'\n')
						mainlog.info('Command : '+screen_cmd+'\n')
						os.system('screen -S '+str(ms_obsid)+'_'+basemsdir+'_polcal_'+str(mpicount)+' -X stuff \"'+screen_cmd+'\n"')	
						if casa_instance*2.0>10:
							sleep_time=casa_instance*2.0
						else:
							sleep_time=10	
						time.sleep(float(sleep_time))
						mpicount+=1
						mpicmd=[]
					break
			time.sleep(2.0)
			if len(polcal_cmd_list)==0:
				polcal_screen_list=[]
				polcal_cmd_list=[]
				polcal_finished_file_list=[]
				finished_polcal=True
				time.sleep(5)
				break	
		
		mainlog.info('All calibration job spawned for ms : '+msname+'\n')
		mainlog.info('#########################\n')
	#	os.system('rm -rf '+ref_timesliced_measurement_set+'* '+ref_timesliced_measurement_set+'.flagversions')
	#	os.system('rm -rf '+ref_averaged_msname+'* '+ref_timesliced_measurement_set+'.flagversions')
	#	os.system('rm -rf '+timesliced_measurement_set+'* '+ref_timesliced_measurement_set+'.flagversions')
	#	os.system('rm -rf '+averaged_msname+'* '+ref_timesliced_measurement_set+'.flagversions')
		del pol_channel_grid
		if ref_time_freq==True:
			np.save(inputs.basedir+'/Ref_time_cal_record',np.array([0,ref_time,ref_freq,ref_chan,cur_spawned_casa_instances,ref_freq_avg,ref_time_avg]))
			os.system('touch '+inputs.basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_0')
		if os.path.exists(inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(cur_spawned_casa_instances))==True:
			os.system('rm -rf '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(cur_spawned_casa_instances))
			os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(cur_spawned_casa_instances))
		else:
			os.system('touch '+inputs.basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(cur_spawned_casa_instances))
		return 0,ref_time,ref_chan,cur_spawned_casa_instances,ref_freq_avg,ref_time_avg


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
		skip_freq_pol=2560
		if inputs.safety_factor==0:
			skip_time		=	960
			skip_freq		=	2560
		elif inputs.safety_factor==1:
			skip_time		=	720
			skip_freq		=	2560
		else:
			skip_time		=	480
			skip_freq		=	2560
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
		ref_freqstamp=str(msfreqs_cal[int(len(msfreqs_cal)/2)])
		ref_timestamp=mjdsec_to_timestamp(mstimes_cal[int(len(mstimes_cal)/2)],format=3)
		ref_timestamp_mjd=mstimes_cal[int(len(mstimes_cal)/2)]
		# Selecting reference time frequnency ms
		######################################
		for i in range(len(measurement_set_list)):
			msname=measurement_set_list[i]
			if ref_freqstamp in msname and ref_timestamp in msname:
				ref_freq_time_msname=msname
				ref_time_freq_metafits=metafits_dic[ref_freq_time_msname.split('.ms')[0].split('_timesliced')[0].split('_ref')[0]]
		mainlog.info('Reference time frequency measurement set : '+ref_freq_time_msname+'\n')
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
		reftimefreq_ms_OBSID=get_OBSID_from_metafits(metafits_dic[ref_freq_time_msname.split('.ms')[0].split('_timesliced')[0].split('_ref')[0]])
		AM=AccessMS(ref_freq_time_msname)
		workdir=inputs.basedir+'/'+os.path.basename(ref_freq_time_msname).split('.ms')[0]
		if os.path.isdir(workdir)==False:
			os.makedirs(workdir)

		spawned_ms_jobs={}
		finished_ms=[]
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
		caltable_list=[]
		
		if calibrator_found==True and len(gaincal_list)!=0: # If calibration found from paircars database or calibrator
			calibrator_OBSID=np.min(np.array([x.split('_')[0] for x in gaincal_list]))
			caltable_list=gaincal_list+bandpass_list
			if abs(reftimefreq_ms_OBSID-calibrator_OBSID)>calibrator_interval*60:
				ref_time_freq=True
			else:
				ref_time_freq=False
			if os.path.isdir(workdir)==False:
				if os.path.exists(workdir)==True:
					os.system('rm -rf '+workdir)
				os.makedirs(workdir)
			os.system('cp -r selfcal_inputs.py '+workdir+'/selfcal_inputs.py')
			return_msg,ref_time,ref_chan,spawned_casa_instances,ref_freq_avg,ref_time_avg=run_paircars_ms(ref_freq_time_msname,ref_time_freq_metafits,workdir,\
					ref_freq_avg=0,ref_time_avg=0,ref_time_freq=ref_time_freq,do_bandpass=inputs.do_bandpass,do_polcal=inputs.do_polcal,\
					num_threads=available_cpu_for_paircars,calibrator_caltable=caltable_list)
			spawned_ms_jobs[ref_freq_time_msname]=[ref_time,ref_chan,spawned_casa_instances,ref_freq_avg,ref_time_avg]
			while True:
				touch_file_list=glob.glob(inputs.basedir+'/.ref_timechan_done_'+str(reftimefreq_ms_OBSID)+'*')
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
					break
				else:
					time.sleep(2.0)
		elif calibrator_found==False and len(touch_file_list)==0:
			if os.path.isdir(workdir)==False:
				if os.path.isfile(workdir)==True:
					os.system('rm -rf '+workdir)
				os.makedirs(workdir)
			os.system('cp -r selfcal_inputs.py '+workdir+'/selfcal_inputs.py')
			result=run_paircars_ms(ref_freq_time_msname,ref_time_freq_metafits,workdir,ref_freq_avg=0,ref_time_avg=0,ref_time_freq=True,do_bandpass=inputs.do_bandpass,\
					do_polcal=inputs.do_polcal,num_threads=available_cpu_for_paircars,calibrator_caltable=[])
			return_msg,ref_time,ref_chan,spawned_casa_instances,ref_freq_avg,ref_time_avg=result
			spawned_ms_jobs[ref_freq_time_msname]=[ref_time,ref_chan,spawned_casa_instances,ref_freq_avg,ref_time_avg]
			while True:
				touch_file_list=glob.glob(inputs.basedir+'/.ref_timechan_done_'+str(reftimefreq_ms_OBSID)+'*')
				if len(touch_file_list)==1:
					runtime_fail=False
					selfcal_fail=False
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
						mainlog.info('Reference time frequency calibration failed because some error occured. Contact developer to fix the problem.\n')
						msfreqs_cal.remove(float(ref_freqstamp))
						mstimes_cal.remove(float(ref_timestamp_mjd))
						ref_timechan_success=False
						break
					else:
						mainlog.info('Reference time frequency calibration is finished.\n') 
						ref_timechan_success=True
					break
				else:
					time.sleep(2.0)
		else:
			mainlog.info('Reference time frequency calibration is done already.\n')
			msfreqs_cal.remove(float(ref_freqstamp))
			mstimes_cal.remove(float(ref_timestamp_mjd))
			spawned_ms_jobs[ref_freq_time_msname]=[ref_timestamp,0,len(touch_file_list1),160,2.0]
			ref_timechan_success=True
			break

	index=measurement_set_list.index(ref_freq_time_msname)
	measurement_set_list.remove(ref_freq_time_msname) # Removing ref time freq ms from ms list
	del metafits_dic[ref_freq_time_msname.split('.ms')[0].split('_timesliced')[0].split('_ref')[0]] # Removing ref time freq metafits from list 
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
					time.sleep(2.0)				
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
else:
	if inputs.timerange!='':
		mainlog.info('No measurement set is present in the timerange : '+inputs.timerange+'\n')
	else:
		mainlog.info('No measurment set is found.\n')
# Applying solution to whole ms
############################### 

# Final imaging mode
####################










#final_finished_file='/'.join(inputs.basedir.split('/')[:-1])+'/.Finished_'+inputs.basedir.split('/')[-1]
#os.system('touch '+final_finished_file)
# TODO : Include nasedir name in screen name so that multiple mode of run can be done simultaneusly on same dataset on different basedir
# TODO : Increase integration time or frequency based on the integration at reference time frequency and scaled sensitivity at observed band
# TODO : Frequency dependent solar size for convolving and phase aligning (Using the expression from Rohit's paper)
# TODO : Cleaning data directory except ms and metafits
# TODO : Make a flag whether internet connection is availbale or not and pass that flag to functions which depends on internet connection
# TODO : Final leakage correction using background sources (using the times where Sun is almost quiet)
# TODO : Calibration database
# TODO : Diagnostic plots
# TODO : Calibrator solution
# TODO : Start from the failed part or stopped part
# End of PAIRCARS
#################
