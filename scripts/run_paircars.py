'''
Code is written by Devojyoti Kansabanik , 28 Jan, 2021
'''

from optparse import OptionParser
import os,sys
if __name__=='__main__':
	usage= ' Perform self calibration of a single time and frequency slice'
	parser = OptionParser(usage=usage)
	parser.add_option('--msname',dest="chantime_msname",default=None,help="Name of measurement set of a single time and frequency slice",metavar="Measurement Set")
	parser.add_option('--metafits',dest="metafits",default=None,help="Name of metafits file of the observation",metavar="Metafits file")
	parser.add_option('--basedir',dest='basedir',default=None,help='Name of the base directory',metavar='Directory path')
	parser.add_option('--workdir',dest='workdir',default=None,help='Name of the working directory',metavar='Directory path')
	parser.add_option('--ref_freq_avg',dest='ref_freq_avg',default=0,help='Frequency averaging for reference ms',metavar="Float")
	parser.add_option('--ref_time_avg',dest="ref_time_avg",default=0,help="Time averaging for reference ms",metavar="Float")
	parser.add_option('--ref_time_freq',dest="ref_time_freq",default=False,help="Reference measurement set or not",metavar="Boolean")
	parser.add_option('--do_bandpass',dest="do_bandpass",default=True,help="Perform bandpass calibration or not",metavar="Boolean")
	parser.add_option('--do_polcal',dest="do_polcal",default=True,help="Perform polarisation calibration or not",metavar="Boolean")
	parser.add_option('--cal_attenuation',dest="calatten",default=1.0,help="Attenuation in dB for calibrator observation",metavar="Float")
	parser.add_option('--num_threads',dest="num_threads",default=0,help="Number of processing threads to use",metavar="Integer")
	parser.add_option('--caltables',dest="caltables",default=None,help="Previous calibration tables",metavar="Comma separated string")
	(options, args) = parser.parse_args()

os.chdir(options.basedir)
sys.path.append(os.getcwd())
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
		ms_mainlog.info('Reference channel : '+str(channel)+' at '+str(freqlist[ch]/10**6)+' MHz and reference time : '+str(timestamp)+'\n')
		# Spliting ref time chan ms
		###########################
		ms_mainlog.info('Spliting reference time and channel..............\n')	
		if os.path.isdir(cwd+'/reftimechan.ms')==True:
			os.system('rm -rf '+cwd+'/reftimechan.ms* '+cwd+'/reftimechan.ms.flagversions')
		split(vis=msname,outputvis=cwd+'/reftimechan.ms',datacolumn=datacolumn,spw=spw,timerange=timestamp)
		ms_mainlog.info('split(vis=\''+msname+'\',outputvis=\''+cwd+'/reftimechan.ms\',datacolumn=\''+datacolumn+'\',spw=\''+spw+'\',timerange=\''+timestamp+'\')\n')
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
		ms_mainlog.info('Spliting measurement set for time : '+timestamp+' and frequency : '+str(freqlist[ch]/10**6)+' MHz ............\n')
		if os.path.isdir(cwd+'/timechan.ms')==True:
			os.system('rm -rf '+cwd+'/timechan.ms* '+cwd+'/timechan.ms.flagversions')
		split(vis=msname,outputvis=cwd+'/timechan.ms',datacolumn=datacolumn,spw=spw,timerange=timestamp)
		ms_mainlog.info('split(vis=\''+msname+'\',outputvis=\''+cwd+'/timechan.ms\',datacolumn=\''+datacolumn+'\',spw=\''+spw+'\',timerange=\''+timestamp+'\')\n')
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

def casa_instance_runner(cmd,basedir,screen_name,finished_touch_file,prefix_cmds=[]):
	'''
	Function to run a casa instance
	Parameters:
	cmd = Command to run
	screen_name = Name of the screen
	'''
	batch_file=basedir+'/'+screen_name+'.batch'
	cmd_batch=basedir+'/'+screen_name+'_cmd.batch'
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
	return basedir+'/'+screen_name+'.batch'

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
	#ms_mainlog.info('flagdata(vis=\''+msname+'\',mode=\'unflag\',spw=\''+unflagchan+'\',antenna=\''+antenna+'\')\n')
	#flagdata(vis=msname,mode='unflag',spw=unflagchan,antenna=antenna)
	cpu_sockets =  int(subprocess.check_output('cat /proc/cpuinfo | grep "physical id" | sort -u | wc -l', shell=True))

	open_casa_instance=0
	ms_obsid=get_OBSID_from_metafits(metafits)
	obs_atten=float(fits.getheader(metafits)['ATTEN_DB'])
	basemsdir=os.path.basename(msname).split('.ms')[0]

	last_selfcal_msg=1
	ref_time_freq_copy=copy.deepcopy(ref_time_freq)
	try:
		last_selfcal_msg,ref_time,ref_freq,ref_chan,spawned_casa_instances,ref_freq_avg,ref_time_avg=np.load(basedir+'/Ref_time_cal_record.npy')
		last_selfcal_msg=int(last_selfcal_msg)
		ref_time=str(ref_time)
		ref_chan=int(ref_chan)
		ref_freq=float(ref_freq)
		spawned_casa_instances=int(spawned_casa_instances)
		ref_freq_avg=float(ref_freq_avg)
		ref_time_avg=float(ref_time_avg)
		if last_selfcal_msg==1:
			ms_mainlog.info('Reference time frequency calibration was failed. Try ms : '+msname+' as reference ms.\n')
			ref_time_freq=True
	except:
		ref_time_freq=copy.deepcopy(ref_time_freq_copy)
		ref_time_freq==True
		pass
	del ref_time_freq_copy
		

	if len(glob.glob(basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'*'))>0:  # Checking whether calibration already done or not 
		# Removing .Finished files if error occured
		###########################################
		touch_file_list=glob.glob(basedir+'/.Finished_*cal*'+str(ms_obsid)+'*'+basemsdir+'*')
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
		touch_file_list=glob.glob(basedir+'/.Finished_*cal*'+str(ms_obsid)+'*'+basemsdir+'*')
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
					touch_file_list.remove(t)
				else:
					msg=int(msg)
					if msg>100:
						msg-=100
					if msg!=0 and msg!=8 and msg!=9:
						os.system('rm -rf '+t)
						touch_file_list.remove(t)
	
		spawned_jobs=int(glob.glob(basedir+'/.Finished_spawned*'+str(ms_obsid)+'*'+basemsdir+'*')[0].split('_')[-1])
		if len(touch_file_list)==spawned_jobs:
			ms_mainlog.info('Calibration has already been done for ms : '+msname+'\n')
			ms_mainlog.info('#########################\n')
		if os.path.exists(basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(spawned_jobs))==True:
			os.system('rm -rf '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(spawned_jobs))
			os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(len(touch_file_list)))
		else:
			os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(len(touch_file_list)))

		if last_selfcal_msg==0 and ref_time_freq==True:
			os.system('touch '+basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_0')
		elif ref_time_freq==True and os.path.exists(basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_0')==False:
			os.system('touch '+basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_error')
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
	ms_mainlog.info(output)
	if inputs.do_decor_correction: # Performing decorrelation correction and IAU convention change
		ms_mainlog.info('Performing de-correlation correction and IAU convention correction for ms : '+msname+'\n')
		decor(msname,metafits,10,False)
	else: # If user do not want decorrelation correction perform only IAU convention
		ms_mainlog.info('Correcting to IAU convention......\n')
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
			ms_mainlog.info('Error in user given time range. Choosing full ms.\n')
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
						ms_mainlog.error('Channel number is wrong. Removing chanstamp '+str(i))
						break
					chan_list.append(a)
				chan_list=sorted(chan_list)
				if chan_list[0]>nchan or chan_list[1]<0:
					ms_mainlog.info('Channel range is not in ms, removing chanstamp '+str(i))
				elif chan_list[0]<0:
					ms_mainlog.info('Start channel is less than 0. Shifted it to channel 0.\n')
					chan_list[0]=0
				elif chan_list[1]>nchan:
					ms_mainlog.info('End chan is greater than total number of channels. Shifted it to total number of channels.\n')
					chan_list[1]=nchan
				for x in range(chan_list[0],chan_list[1]): 
					new_chan_list.append(str(x))
			
			if len(new_chan_list)!=0:
				new_chan_list=','.join(new_chan_list)
			else:
				ms_mainlog.info('Given channels are not in ms. Continue with all channels.\n')
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
		ms_mainlog.info('Frequency resolution of the data is less than the intended imaging frequency resolution. Setting imaging frequency resolution to frequency of the data\n')
		inputs.image_delta_freq=freqres
	if inputs.image_delta_time<timeres:
		ms_mainlog.info('Time resolution of the data is less than the intended imaging time resolution. Setting imaging time resolution to time resolution of the data\n')
		inputs.image_delta_time=timeres
	if inputs.skip_freq<freqres:
		ms_mainlog.info('Frequency resolution of the data is less than the skip frequency. Setting skip frequency resolution to frequency of the data.\n')
		inputs.skip_freq=freqres
	if inputs.skip_time<timeres:
		ms_mainlog.info('Time resolution of the data is less than the skip time. Setting skip time to the time resolution of the data.\n')
		inputs.skip_time=timeres
	if inputs.image_freq<freqres:
		ms_mainlog.info('Frequency resolution of the data is less than the image bandwidth. Setting image bandwidth resolution to frequency of the data.\n')
		inputs.image_freq=freqres
	elif inputs.image_freq>inputs.image_delta_freq:
		ms_mainlog.info('Image bandwidth is greater than image frequency interval. Setting image bandwidth to image frequency interval.\n')
		inputs.image_freq=inputs.image_delta_freq
	if inputs.image_time>inputs.image_delta_time:
		ms_mainlog.info('Time resolution of the image is greater than the time interval. Setting time resolution to the time interval.\n')
		inputs.image_time=inputs.image_delta_time

	ms_mainlog.info('###############################################\n')
	ms_mainlog.info('Skip frequency : '+str(inputs.skip_freq)+' kHz\n')
	ms_mainlog.info('Skip time : '+str(inputs.skip_time)+' s\n')
	ms_mainlog.info('Image frequency interval : '+str(inputs.image_delta_freq)+' kHz\n')
	ms_mainlog.info('Image time interval : '+str(inputs.image_delta_time)+' s\n')
	ms_mainlog.info('Image bandwidth : '+str(inputs.image_freq)+' kHz\n')
	ms_mainlog.info('Image time resolution : '+str(inputs.image_time)+' s\n')
	ms_mainlog.info('Channel range to image : '+str(new_chan_list)+'\n')
	ms_mainlog.info('Time range to image : '+str(new_timerange)+'\n')
	ms_mainlog.info('###############################################\n')

	# Flagging coarse channel edges and center
	##########################################
	ms_mainlog.info('flag_MWA_coarse(\''+msname+'\',edgewidth=160,do_flag=True)\n')
	good_channels,channel_per_coarse=flag_MWA_coarse(msname,edgewidth=160,do_flag=True)
	unflag_channels=np.array(AM.get_unflag_chan(flagfrac=1))
	new_chan_list=np.array(new_chan_list.split(','),dtype='int')
	new_unflag_chan_list=np.intersect1d(unflag_channels,new_chan_list)
	start_chan=np.min(new_unflag_chan_list)
	end_chan=np.max(new_unflag_chan_list)

	# Spliting timeranges for reference times
	#########################################
	ms_mainlog.info('Measurement set timeslice range : '+new_timerange+'\n')

	if new_timerange!=ms_timerange:
		ms_mainlog.info('Spliting reference timerange : '+str(new_timerange)+'\n')	
		if os.path.isdir(workdir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftimesliced.ms')==True:
			os.system('rm -rf '+workdir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftimesliced.ms* '+\
						workdir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftimesliced.ms.flagversions')
		ms_mainlog.info('split(vis=\''+msname+'\',outputvis=\''+workdir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftimesliced.ms\',timerange=\''+new_timerange\
					+'\',datacolumn=\'data\')\n')
		split(vis=msname,outputvis=workdir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftimesliced.ms',\
				timerange=new_timerange,datacolumn='data')
		ref_timesliced_measurement_set=workdir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftimesliced.ms'
	else:
		ms_mainlog.info('Reference timerange is similar to MS timerange. Linking the MS.....\n')
		ref_timesliced_measurement_set=workdir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftimesliced.ms'
		if os.path.exists(ref_timesliced_measurement_set):
			if os.path.islink(ref_timesliced_measurement_set):
				os.unlink(ref_timesliced_measurement_set)
			else:
				os.system('rm -rf '+ref_timesliced_measurement_set)
		elif os.path.islink(ref_timesliced_measurement_set):
			os.unlink(ref_timesliced_measurement_set)
		os.system('ln -s '+os.path.realpath(msname)+' '+ref_timesliced_measurement_set)

	ms_mainlog.info('Linking reference timeslice ms to timesliced ms....\n')
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
	
	touch_file_list=glob.glob(basedir+'/.Finished_*cal*'+str(ms_obsid)+'*'+basemsdir+'*')
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
	touch_file_list=glob.glob(basedir+'/.Finished_gcal*'+str(ms_obsid)+'*'+basemsdir+'*ref*')
	if len(touch_file_list)>0:
		touch_file=touch_file_list[0]
	else:
		touch_file=''
	ref_time=''
	ref_chan=''
	ref_freq=''
	ref_freq_avg=freq_avg
	ref_time_avg=time_avg
	ref_timechan_done=False
	# In this while loop we are checking whether the present time and frequency averaging is enough to start the self calibration. If it has it will leave the loop and go for selfcal
	##################################################################################################################################################################################
	while ref_timechan_done==False and ref_time_freq==True and os.path.exists(touch_file)==False:  	
		ms_mainlog.info('Choosing averaging frequency width : '+str(ref_freq_avg)+' kHz, averaging temporal width : '+str(ref_time_avg)+' s, Skip frequency : '\
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
			ms_mainlog.info('Avearging reference time measurement width frequency average :'+str(ref_freq_avg)+' kHz, temporal average :'+str(ref_time_avg)+'s\n')
			ms_mainlog.info('split(vis=\''+ref_timesliced_measurement_set+'\',outputvis=\''+ref_averaged_msname+'\',width='+\
						str(chan_width)+',timerange=\''+new_timerange+'\',timebin=\''+str(ref_time_avg)+'s\',datacolumn=\'data\')\n')
			split(vis=ref_timesliced_measurement_set,outputvis=ref_averaged_msname,width=chan_width,timerange=new_timerange,timebin=str(ref_time_avg)+'s',datacolumn='data')
			if ref_time_freq==True:
				np.save(basedir+'/Ref_time_cal_record',np.array([0,ref_time,ref_freq,ref_chan,cur_spawned_casa_instances,ref_freq_avg,ref_time_avg]))
			do_averaging=True
		else:
			ms_mainlog.info('No averaging is required.\n')
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
				np.save(basedir+'/Ref_time_cal_record',np.array([0,ref_time,ref_freq,ref_chan,cur_spawned_casa_instances,ref_freq_avg,ref_time_avg]))

		AMref=AccessMS(ref_averaged_msname)
		# Making reference time and channel and time frequency grid
		##########################################################
		unflagged_channels=AMref.get_unflag_chan(flagfrac=1) # Unflagged averaged channels
		if len(unflagged_channels)==0:
			ms_mainlog.info('No unflagged channel is present.\n')
			if ref_time_freq==True:	
				np.save(basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))	
				os.system('touch '+basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(msname)+'_selfcalerror')	
			if os.path.exists(basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
				os.system('rm -rf '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
			else:
				os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
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

	#	ref_time_copy=copy.deepcopy(ref_time) # Copy this timestamp to remove from ref time grid if failed
		AMref=AccessMS(ref_averaged_msname)
		# Making reference time and channel and time frequency grid
		##########################################################
		unflagged_channels=AMref.get_unflag_chan(flagfrac=1) # Unflagged averaged channels
		if len(unflagged_channels)==0:
			ms_mainlog.info('No unflagged channel is present.\n')
			if ref_time_freq==True:	
				os.system('touch '+basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(msname)+'_selfcalerror')
				np.save(basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))	
			if os.path.exists(basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
				os.system('rm -rf '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
			else:
				os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))	
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

		ms_mainlog.info('Reference frequency : '+str(ref_freqs[ref_chan])+' MHz.\n')
		ms_mainlog.info('Reference time : '+str(ref_time)+'\n')

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
		ms_mainlog.info('Available cpus for P-AIRCARS: '+str(available_cpu_for_paircars)+'\n')
		ms_mainlog.info('Total number of available CASA instances : '+str(casa_instance)+'\n')

		# Spliting ref time chan ms
		###########################
		ref_timechan_ms,ref_timechan_dir=spliting_timechan(ref_averaged_msname,metafits,ref_chan,ref_time,caltype='G',ref_timechan=True,\
						input_file=workdir+'/selfcal_inputs.py',datacolumn='data')
		cur_workdir=ref_timechan_dir
		mpicmd_file=basedir+'/'+basemsdir+'.ref_mpicmd'
		# Run selfcal
		while True:
			ref_time_chan_loop_count+=1
			if mpi==0: # MPI command file
				if os.path.isfile(mpicmd_file):
					os.system('rm -rf '+mpicmd_file)
			try:
				touch_file_list=glob.glob(basedir+'/.Finished_gcal*'+str(ms_obsid)+'*'+basemsdir+'*'+os.path.basename(ref_timechan_ms)+'*')
				if reduce_flag_count==1: # Restarting calibration with more time averaging
					if previous_caltable!='' and previous_record!='':
						try:
							ms_mainlog.info('applycal(vis=\''+ref_timechan_ms+'\',gaintable=\''+previous_caltable+'\',applymode=\'calflag\',flagbackup=True)\n')
							applycal(vis=ref_timechan_ms,gaintable=previous_caltable,applymode='calflag',flagbackup=True)
							if os.path.isdir(cur_workdir+'/junk1.ms')==True:
								os.system('rm -rf '+cur_workdir+'/junk1.ms')
							ms_mainlog.info('cp -r '+ref_timechan_ms+' '+cur_workdir+'/junk1.ms\n')
							os.system('cp -r '+ref_timechan_ms+' '+cur_workdir+'/junk1.ms')
							if os.path.isdir(cur_workdir+'/junk1.cal')==True:
								os.system('rm -rf '+cur_workdir+'/junk1.cal')
							ms_mainlog.info('cp -r '+previous_caltable+' '+cur_workdir+'/junk1.cal\n')
							os.system('cp -r '+previous_caltable+' '+cur_workdir+'/junk1.cal')
							print (glob.glob(workdir+'/presession_backup/freq_*datetime*'))
							print (workdir+'/presession_backup/freq_*datetime*')
							if len(glob.glob(workdir+'/presession_backup/freq_*datetime*'))>0:
								ms_mainlog.info('Copying previous round backup directories : '+str(glob.glob(workdir+'/presession_backup/freq_*datetime*'))+'\n')
								if os.path.isdir(cur_workdir+'/pre_backup')==False:
									os.makedirs(cur_workdir+'/pre_backup')
								pre_backups=glob.glob(workdir+'/freq_*datetime*')
								total_prebackups=len(glob.glob(cur_workdir+'/pre_backup/'))
								for i in range(len(pre_backups)):
									j=pre_backups[i]
									ms_mainlog.info('cp -r '+j+' '+cur_workdir+'/pre_backup/'+j+'_'+str(total_prebackups+1)+'\n')
									os.system('cp -r '+j+' '+cur_workdir+'/pre_backup/'+j+'_'+str(total_prebackups+1))
									total_prebackups+=1
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
								ms_mainlog.info('Copying previous log....\n')
								ms_mainlog.info('cp -r '+prelog+' '+cur_workdir+'/Intensity_Selfcal.log\n')
								os.system('cp -r '+prelog+' '+cur_workdir+'/Intensity_Selfcal.log')
							if preverboselog!='':
								ms_mainlog.info('Copying previous verbose log......\n')
								ms_mainlog.info('cp -r '+preverboselog+' '+cur_workdir+'/Intensity_Selfcal_verbose.log\n')
								os.system('cp -r '+preverboselog+' '+cur_workdir+'/Intensity_Selfcal_verbose.log')
							if prerms!='':
								ms_mainlog.info('Copying previous DR_rms record......\n')
								ms_mainlog.info('cp -r '+prerms+' '+cur_workdir+'/DR_rms.npy\n')
								os.system('cp -r '+prerms+' '+cur_workdir+'/DR_rms.npy')
							if preneg!='':
								ms_mainlog.info('Copying previous DR_neg record......\n')
								ms_mainlog.info('cp -r '+preneg+' '+cur_workdir+'/DR_neg.npy\n')
								os.system('cp -r '+preneg+' '+cur_workdir+'/DR_neg.npy')
							os.system('rm -rf '+previous_caltable+' '+previous_selfcal_record+' '+previous_ms+' '+prelog+' '+preverboselog+' '+prerms+' '+previous_record+' '+preneg)
				#+' '+\									workdir+'/freq_*datetime*')
							fresh=False
						except:
							fresh=True
					else:
						fresh=True
				else:
					fresh=True	
				ms_mainlog.info('Perforimg self-calibration using fresh = '+str(fresh)+'\n')
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
				finished_touch_file=basedir+'/.Finished_gcal_'+str(ms_obsid)+'_'+basemsdir+'_'+os.path.basename(ref_timechan_ms)
				screen_batch_file=casa_instance_runner(cmd,basedir,screen_name,finished_touch_file)
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
					ms_mainlog.info('########################\n')
					ms_mainlog.info('Made Screen : '+screen_name+'\n')
					ms_mainlog.info('Command : '+cmd+'\n')
					os.system('screen -S '+screen_name+' -X stuff \"'+screen_cmd+'\n"')	
				elif mpi==0:
					mpicmd=['-np 1 --bind-to core --map-by ppr:'+str(int(available_cpu_for_paircars))+':node:pe=4 -x OMP_NUM_THREADS='+\
						str(available_cpu_for_paircars)+' sh '+screen_batch_file+'\n']
					mpicmd.append('-np 1 sleep 1\n')
					ms_mainlog.info('MPI commands .....\n')
					for i in mpicmd:
						ms_mainlog.info(i)
					mpifil.writelines(mpicmd)
					mpifil.close()
					os.system('chmod a+rwx '+mpicmd_file)
					screen_cmd='mpirun --app '+mpicmd_file
					os.system('screen -S '+str(ms_obsid)+'_'+basemsdir+'_refcal -X quit')	
					time.sleep(0.5)
					os.system('screen -mdS '+str(ms_obsid)+'_'+basemsdir+'_refcal')
					time.sleep(0.5)
					ms_mainlog.info('########################\n')
					ms_mainlog.info('Made Screen : '+str(ms_obsid)+'_'+basemsdir+'_refcal\n')
					ms_mainlog.info('Command : '+cmd+'\n')
					os.system('screen -S '+str(ms_obsid)+'_'+basemsdir+'_refcal -X stuff \"'+screen_cmd+'\n"')	
					time.sleep(2.0)	
				ms_mainlog.info('Self calibration for ms : '+ref_timechan_ms+' is spawned in screen : '+screen_name+'\n')
				ms_mainlog.info('Waiting to finish self calibration for reference time frequency ms :'+ref_timechan_ms+'................\n') 
				while True:
					time.sleep(2)
					touch_file_list=glob.glob(basedir+'/.Finished_gcal*'+str(ms_obsid)+'*'+basemsdir+'*'+os.path.basename(ref_timechan_ms)+'*')
					if len(touch_file_list)!=0:
						msg=touch_file_list[0].split('_')[-1]
						break	
			except Exception as e: # If runtime error occured
				ms_mainlog.error('Error occured :'+str(e)+'\n')
				ms_mainlog.error('Error in running selfcal.\n')
				os.system('rm -rf '+basedir+'/.paircars_running')
				os.system('touch '+basedir+'/.paircars_failed')
				os.system('rm -rf '+cur_workdir)
				if ref_time_freq==True:
					os.system('touch '+basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_runtimeerror')		
					np.save(basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))
				if os.path.exists(basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
					os.system('rm -rf '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
					os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				else:
					os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				return 1,0,0,0,0,0
			if msg=='error': # If error occured in run_intensity_selfcal
				ms_mainlog.info('Runtime error occured.\n')
				os.system('rm -rf '+cur_workdir)	
				if ref_time_freq==True:
					os.system('touch '+basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_runtimeerror')	
					np.save(basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))
				if os.path.exists(basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
					os.system('rm -rf '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
					os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				else:
					os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				return 1,0,0,0,0,0
			elif msg=='noms': # If ms is not present
				ms_mainlog.info('Runtime error occured. No measurement set found.\n')
				os.system('rm -rf '+cur_workdir)	
				if ref_time_freq==True:
					os.system('touch '+basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_runtimeerror')	
					np.save(basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))
				if os.path.exists(basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
					os.system('rm -rf '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
					os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				else:
					os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				return 1,0,0,0,0,0
			elif msg=='nometa': # If metafits not present
				ms_mainlog.info('Runtime error occured. No metafits file found.\n')
				os.system('rm -rf '+cur_workdir)	
				if ref_time_freq==True:
					os.system('touch '+basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_runtimeerror')	
					np.save(basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))
				if os.path.exists(basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
					os.system('rm -rf '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
					os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				else:
					os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				return 1,0,0,0,0,0
			elif msg=='moreflag' and try_reduce_flag==True: # Try to reduce more flags with more time averaging
				time.sleep(2)
				ms_mainlog.info('More than 5 % solutions are flagged. Increasing time averaging.\n')
				new_time_avg=int(ref_time_avg+4.0) # Averaging extra 4 seconds
				if try_reduce_flag==True and reduce_flag_count<1:
					try_reduce_flag=False
					reduce_moreflag=False
					reduce_flag_count+=1
					if os.path.exists(workdir+'/presession_backup')==False:
						os.makedirs(workdir+'/presession_backup')
					if os.path.exists(workdir+'/presession_backup/prerecord.npy')==True:
						os.system('rm -rf '+workdir+'/presession_backup/prerecord.npy')
					os.system('cp -r '+cur_workdir+'/Intensity_selfcal_record.npy '+workdir+'/presession_backup/prerecord.npy')		
					if os.path.exists(workdir+'/presession_backup/precal.cal')==True:
						os.system('rm -rf '+workdir+'/presession_backup/precal.cal')			
					os.system('cp -r '+cur_workdir+'/junk.precal '+workdir+'/presession_backup/precal.cal')
					if os.path.exists(workdir+'/presession_backup/pre.log')==True:
						os.system('rm -rf '+workdir+'/presession_backup/pre.log')
					os.system('cp -r '+cur_workdir+'/Intensity_Selfcal.log '+workdir+'/presession_backup/pre.log')
					if os.path.exists(workdir+'/presession_backup/pre_DR_neg.npy')==True:
						os.system('rm -rf '+workdir+'/presession_backup/pre_DR_neg.npy')
					os.system('cp -r '+cur_workdir+'/DR_neg.npy '+workdir+'/presession_backup/pre_DR_neg.npy')
					if os.path.exists(workdir+'/presession_backup/pre_DR_rms.npy')==True:
						os.system('rm -rf '+workdir+'/presession_backup/pre_DR_rms.npy')
					os.system('cp -r '+cur_workdir+'/DR_rms.npy '+workdir+'/presession_backup/pre_DR_rms.npy')
					prerms=workdir+'/presession_backup/pre_DR_rms.npy'
					preneg=workdir+'/presession_backup/pre_DR_neg.npy'
					prelog=workdir+'/presession_backup/pre.log'
					if os.path.exists(cur_workdir+'/Intensity_Selfcal_verbose.log')==True:
						if os.path.exists(workdir+'/presession_backup/preverbose.log')==True:
							os.system('rm -rf '+workdir+'/presession_backup/preverbose.log')
						os.system('cp -r '+cur_workdir+'/Intensity_Selfcal_verbose.log '+workdir+'/presession_backup/preverbose.log')
						preverboselog=workdir+'/presession_backup/preverbose.log'
					else:
						preverboselog=''
					os.system('rm -rf '+cur_workdir)
					previous_caltable=workdir+'/presession_backup/precal.cal'
					previous_record=workdir+'/presession_backup/prerecord.npy'
					preworkdir=cur_workdir
				for i in touch_file_list:
					msg=i.split('_')[-1]
					if type(msg)==str and msg=='moreflag':						
						os.system('rm -rf '+i)
				if new_time_avg<=10.0:
					ref_time_avg=new_time_avg
					ms_mainlog.info('Increasing time averaging to '+str(ref_time_avg)+'s\n')
					break
			elif int(msg)>=100:
				msg=int(msg)-100
			if int(msg)==10 and pass_flag==False:  # Checking for selfcal SNR and increasing the time and frequency averaging if required
				ms_mainlog.error('SNR for self calibration is not sufficent.\n')
				selfcal_snr=float(np.load(basedir+'/selfcal_minsnr.npy'))
				ms_mainlog.info('Selfcal SNR : '+str(selfcal_snr)+'\n')
				new_time_avg=int(ref_time_avg*np.sqrt(selfcal_snr/min_selfcal_snr)*3)
				if new_time_avg==ref_time_avg:
					new_time_avg+=AMref.calc_timeres()
				if new_time_avg<=skip_time and new_time_avg>ref_time_avg and new_time_avg<=total_time:
					ref_time_avg=new_time_avg
					ms_mainlog.info('Increasing time averaging to '+str(ref_time_avg)+'s\n')
					os.system('rm -rf '+touch_file_list[0])
					break
				new_freq_avg=int(freq_avg*np.sqrt(selfcal_snr/min_selfcal_snr)*3)
				if new_freq_avg<=skip_freq and new_freq_avg>freq_avg:
					ms_mainlog.info('Time averaging reached skip time limit. Increasing frequency averaging.')
					if new_freq_avg<=skip_freq and new_freq_avg>freq_avg:
						freq_avg=new_freq_avg
						ms_mainlog.info('Increasing frequency averaging to '+str(freq_avg)+'kHz\n')
						os.system('rm -rf '+touch_file_list[0])
						break
				else:
					if selfcal_snr>2: # If selfcal_snr is not improved after all trials but it is greater than 2, then contine, otherwise exit
						ms_mainlog.info('Both time and frequency averaging has been tried. Still SNR is not sufficient but is greater than 2. Thus continuing with present averaging.\n')
						pass_flag=True						
						continue
					else:
						ms_mainlog.info('Both time and frequency averaging has been tried. Still SNR is not sufficient. Trying with other time frequency.\n')
						os.system('rm -rf '+cur_workdir)
						if ref_time_freq==True:
							os.system('touch '+basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_selfcalerror')	
							np.save(basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))	
						if os.path.exists(basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
							os.system('rm -rf '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
							os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
						else:
							os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))	
						return 1,0,0,0,0,0
			elif int(msg)!=0 and int(msg)!=9 and int(msg)!=8: # If not succeeded, or max iteration reached or DR decreased but more than min DR, removing the ref time
				ms_mainlog.info('Message : '+error_msgs(100)+' : '+error_msgs(int(msg))+'\n')
				if inputs.verbose==False:
					ms_mainlog.info('Removing the directory : '+ref_timechan_dir)
					os.system('rm -rf '+ref_timechan_dir)
				if len(ref_time_grid)!=0:
					ms_mainlog.info('Removing timestamp : '+str(ref_time)+' from time grid.\n')
					ref_time_grid.remove(ref_time)
					ref_index=int(len(ref_time_grid)/2)
					ref_time=ref_time_grid[ref_index]
					ms_mainlog.info('Trying for new timestamp :'+str(ref_time)+'\n')
					if inputs.verbose==False:
						os.system('rm -rf '+ref_timechan_dir)
					ref_timechan_ms,ref_timechan_dir=spliting_timechan(ref_averaged_msname,metafits,ref_chan,ref_time,caltype='G',ref_timechan=True,input_file=workdir\
									+'/selfcal_inputs.py',datacolumn='data')
					cur_workdir=ref_timechan_dir
					continue
				elif len(ref_channel_grid)!=0:
					ref_time_grid=copy.deepcopy(ref_time_grid_copy)
					ref_channel_grid.remove(ref_chan)
					ref_chan_index=int(len(ref_channel_grid)/2)
					ref_chan=ref_channel_grid[ref_chan_index]
					ref_freq=float(AMref.get_freqs()[ref_chan]/10**6)
					ms_mainlog.info('Trying with new channel : '+str(ref_chan)+'\n')
					if inputs.verbose==False:
						os.system('rm -rf '+ref_timechan_dir)
					ref_timechan_ms,ref_timechan_dir=spliting_timechan(ref_averaged_msname,metafits,ref_chan,ref_time,caltype='G',ref_timechan=True,input_file=\
								workdir+'/selfcal_inputs.py',datacolumn='data')
					cur_workdir=ref_timechan_dir
					continue
				else:
					ms_mainlog.info('Reference imaging has been tried over full measurement set. No good starting point is found. Exiting PAIRCARS...\n')	
					os.system('rm -rf '+basedir+'/.paircars_running')
					os.system('touch '+basedir+'/.paircars_failed')
					os.system('rm -rf '+cur_workdir)	
					if ref_time_freq==True:	
						os.system('touch '+basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_selfcalerror')		
						np.save(basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))
					if os.path.exists(basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
						os.system('rm -rf '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
						os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
					else:
						os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
					return 1,0,0,0,0,0
			elif int(msg)==0 or int(msg)==8 or int(msg)==9: # if succeeded or max iteration reached or DR decreases but more than min DR
				os.system('touch '+basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_'+str(msg))	
				ref_timechan_done=True	
				ms_mainlog.info('Reference time frequency calibration done.\n')
				cur_spawned_casa_instances+=1	
				if ref_time_freq==True:
					np.save(basedir+'/Ref_time_cal_record',np.array([0,ref_time,ref_freq,ref_chan,cur_spawned_casa_instances,ref_freq_avg,ref_time_avg]))
				ref_time_grid.remove(ref_time)
				ref_channel_grid.remove(ref_chan)
				break

	try:
		del ref_channel_grid
		del ref_time_grid
	except:
		pass

	try:
		last_selfcal_msg,ref_time,ref_freq,ref_chan,spawned_casa_instances,ref_freq_avg,ref_time_avg=np.load(basedir+'/Ref_time_cal_record.npy')
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
	if ref_time_freq==False or (len(glob.glob(basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_*'))==0 and\
			 ref_time_freq==True and ref_time_chan_loop_count==0):
		if ref_time_freq==False:
			try:
				last_selfcal_msg,ref_time,ref_freq,ref_chan,spawned_casa_instances,ref_freq_avg,ref_time_avg=np.load(basedir+'/Ref_time_cal_record.npy')
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
			ms_mainlog.info('Choosing averaging frequency width : '+str(ref_freq_avg)+' kHz, averaging temporal width : '+str(ref_time_avg)+' s, Skip frequency : '\
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
				ms_mainlog.info('Avearging reference time  measurement width frequency average :'+str(ref_freq_avg)+' kHz, temporal average :'+str(2)+'s\n')
				ms_mainlog.info('split(vis=\''+ref_timesliced_measurement_set+'\',outputvis=\''+ref_averaged_msname+'\',width='+\
							str(chan_width)+',timerange=\''+new_timerange+'\',timebin=\''+str(ref_time_avg)+'s\',datacolumn=\'data\')\n')
				split(vis=ref_timesliced_measurement_set,outputvis=ref_averaged_msname,width=chan_width,timerange=new_timerange,timebin=str(ref_time_avg)+'s',datacolumn='data')
				do_averaging=True
			else:
				ms_mainlog.info('No averaging is required.\n')
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
				ms_mainlog.info('No unflagged channel is present.\n')
				if ref_time_freq==True:	
					os.system('touch '+basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_selfcalerror')	
					np.save(basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))
				if os.path.exists(basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
					os.system('rm -rf '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
					os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				else:
					os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))	
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
			#ref_time_copy=copy.deepcopy(ref_time) # Copy this timestamp to remove from ref time grid if failed
			try:
				del ref_time_grid
			except:
				pass			

		elif ref_time_freq==True:
			try:
				last_selfcal_msg,ref_time,ref_freq,ref_chan,spawned_casa_instances,ref_freq_avg,ref_time_avg=np.load(basedir+'/Ref_time_cal_record.npy')
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
					ms_mainlog.info('Avearging reference time measurement width frequency average :'+str(ref_freq_avg)+' kHz, temporal average :'+str(ref_time_avg)+'s\n')
					ms_mainlog.info('split(vis=\''+ref_timesliced_measurement_set+'\',outputvis=\''+ref_averaged_msname+'\',width='+\
								str(chan_width)+',timerange=\''+new_timerange+'\',timebin=\''+str(ref_time_avg)+'s\',datacolumn=\'data\')\n')
					split(vis=ref_timesliced_measurement_set,outputvis=ref_averaged_msname,width=chan_width,timerange=new_timerange,timebin=str(ref_time_avg)+'s',datacolumn='data')
				else:
					ms_mainlog.info('No averaging is required.\n')
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
		last_selfcal_msg,ref_time,ref_freq,ref_chan,spawned_casa_instances,ref_freq_avg,ref_time_avg=np.load(basedir+'/Ref_time_cal_record.npy')
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
		ms_mainlog.info('Linking averaged reference timesliced ms to timesliced ms.....\n')
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
		ms_mainlog.info('No averaging is required.\n')
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
		ms_mainlog.info('No unflagged channel is present.\n')
		if ref_time_freq==True:
			os.system('touch '+basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_selfcalerror')		
			np.save(basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))	
		if os.path.exists(basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
			os.system('rm -rf '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
			os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
		else:
			os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
		return 1,0,0,0,0,0

	timestamps=AMref.get_timestamps()
	total_time=AMref.calc_total_time()
	ref_freqs=AMref.get_freqs()/10**6
	unflagged_channels=AMref.get_unflag_chan(flagfrac=1) # Unflagged averaged channels
	if len(unflagged_channels)==0:
		ms_mainlog.info('No unflagged channel is present.\n')
		if ref_time_freq==True:
			os.system('touch '+basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_selfcalerror')		
			np.save(basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))	
		if os.path.exists(basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
			os.system('rm -rf '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
			os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
		else:
			os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
		return 1,0,0,0,0,0

	ms_mainlog.info('Reference frequency : '+str(ref_freqs[ref_chan])+' MHz.\n')
	ms_mainlog.info('Reference time : '+str(ref_time)+'\n')
	try:
		del ref_channel_grid
		del ref_time_grid
	except:
		pass

	AM=AccessMS(averaged_msname)
	unflagged_channels=AM.get_unflag_chan(flagfrac=1) # Unflagged averaged channels
	if len(unflagged_channels)==0:	
		ms_mainlog.info('No unflagged channel is present.\n')
		if ref_freq_time==True:
			os.system('touch '+basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_selfcalerror')		
			np.save(basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))	
		if os.path.exists(basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
			os.system('rm -rf '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
			os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
		else:
			os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
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

	ms_mainlog.info('Channel grid list : '+str(channel_grid)+'\n')
	ms_mainlog.info('Timestamp grid list : '+str(time_grid)+'\n')
	ms_mainlog.info('Reference time : '+str(ref_time)+'\n')
	ms_mainlog.info('Reference frequency : '+str(ref_freqs[ref_chan])+' MHz.\n')

	# Making ref time freq gaintable list
	#####################################
	if ref_time_freq==True:
		ref_timechan_caltable=glob.glob(basedir+'/caltables/'+str(ms_obsid)+'/'+basemsdir+'/*ref*.cal')
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
	ms_mainlog.info('Available cpus for P-AIRCARS: '+str(available_cpu_for_paircars)+'\n')
	ms_mainlog.info('Total number of available CASA instances : '+str(casa_instance)+'\n')

	# Spliting gaincal measurement set
	##################################
	gaincal_cmd_list=[]
	gaincal_screen_list=[]
	gaincal_finished_file_list=[]
	batch_file_list=[]
	calstring=','.join(ref_gaintable)
	temp_ref_ms=''
	if ref_time_freq==True:
		ref_touch_list=glob.glob(basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_*')
		if len(ref_touch_list)==0:
			os.system('touch '+basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_0')		
	if len(time_grid)!=0:
		ms_mainlog.info('Spliting reference chan data for performing gaincal.\n')
		for timestamp in time_grid:
			splited_msname,splited_msdir=spliting_timechan(averaged_msname,metafits,str(ref_chan),timestamp,caltype='G',ref_timechan=False,\
											input_file=workdir+'/selfcal_inputs.py',datacolumn='data')
			if timestamp==ref_time:
				temp_ref_ms=splited_msname
			touch_file=glob.glob(basedir+'/.Finished_gcal_'+str(ms_obsid)+'_'+basemsdir+'_'+os.path.basename(splited_msname)+'*')
			if len(touch_file)==0:
				cmd='run_intensity_selfcal --msname '+splited_msname+' --metafits '+metafits+' --workdir '+splited_msdir+\
					' --dopoint True --verbose '+str(inputs.verbose)+' --interactive '+str(inputs.interactive)+' --fresh True --reduce_flags True --caltables '+calstring
				gaincal_cmd_list.append(cmd)
				gaincal_screen_list.append(str(ms_obsid)+'_'+os.path.basename(splited_msname).split('.ms')[0]+'_screen_G')
				gaincal_finished_file_list.append(basedir+'/.Finished_gcal_'+str(ms_obsid)+'_'+basemsdir+'_'+os.path.basename(splited_msname))
		mpicmd=[]
		mpicount=0
		while len(gaincal_cmd_list)!=0:  # Loop while all gaincal cmds are spawned
			touch_file_list=glob.glob(basedir+'/.Finished_*cal*'+str(ms_obsid)+'*'+basemsdir+'*')
			if (len(touch_file_list)-touch_count)>0:
				if open_casa_instance>1:
					open_casa_instance-=(len(touch_file_list)-touch_count)
					ms_mainlog.info('New CASA instance available : '+str((len(touch_file_list)-touch_count))+'\n')
			while len(gaincal_screen_list)!=0 and (casa_instance-open_casa_instance)>=1:
				screen_name=gaincal_screen_list[0]
				cmd=gaincal_cmd_list[0]
				finished_file=gaincal_finished_file_list[0]
				screen_batch_file=casa_instance_runner(cmd,basedir,screen_name,finished_file)
				batch_file_list.append(screen_batch_file)
				touch_count+=1
				if mpi==1:
					screen_cmd='sh '+screen_batch_file
					os.system('screen -S '+screen_name+' -X quit')	
					time.sleep(0.5)
					os.system('screen -mdS '+screen_name)
					time.sleep(0.5)
					ms_mainlog.info('########################\n')
					ms_mainlog.info('Made Screen : '+screen_name+'\n')
					ms_mainlog.info('Command : '+cmd+'\n')
					os.system('screen -S '+screen_name+' -X stuff \"'+screen_cmd+'\n"')	
				elif mpi==0:
					mpicmd.append('-np 1 --bind-to core --map-by ppr:'+str(int(available_cpu_for_paircars))+':node:pe=4 -x OMP_NUM_THREADS='+\
							str(available_cpu_for_paircars)+' sh '+screen_batch_file+'\n')
					mpicmd.append('-np 1 sleep 1\n')
				open_casa_instance+=1
				cur_spawned_casa_instances+=1
				gaincal_screen_list.remove(screen_name)
				gaincal_cmd_list.remove(cmd)
				gaincal_finished_file_list.remove(finished_file)
				if open_casa_instance>=casa_instance or len(gaincal_screen_list)==0:
					ms_mainlog.info('Maximum casa instances spawned. Waiting for complete those jobs.\n')
					if mpi==0:
						mpicmd_file=basedir+'/'+basemsdir+'.gcal_mpicmd_'+str(mpicount)
						if os.path.exists(mpicmd_file):
							os.system('rm -rf '+mpicmd_file)
						mpifil=open(mpicmd_file,'w')
						ms_mainlog.info('MPI commands .....\n')
						for i in mpicmd:
							ms_mainlog.info(i)
						mpifil.writelines(mpicmd)
						mpifil.close()
						os.system('chmod a+rwx '+mpicmd_file)
						screen_cmd='mpirun --app '+mpicmd_file
						os.system('screen -S '+str(ms_obsid)+'_'+basemsdir+'_gcal_'+str(mpicount)+' -X quit')	
						time.sleep(0.5)
						os.system('screen -mdS '+str(ms_obsid)+'_'+basemsdir+'_gcal_'+str(mpicount))
						time.sleep(0.5)
						ms_mainlog.info('########################\n')
						ms_mainlog.info('Made Screen : '+str(ms_obsid)+'_'+basemsdir+'_gcal_'+str(mpicount)+'\n')
						ms_mainlog.info('Command : '+screen_cmd+'\n')
						os.system('screen -S '+str(ms_obsid)+'_'+basemsdir+'_gcal_'+str(mpicount)+' -X stuff \"'+screen_cmd+'\n"')
						if casa_instance*2.0>10:
							sleep_time=casa_instance*2.0
						else:
							sleep_time=10	
						time.sleep(float(sleep_time))
						mpicount+=1
						mpicmd=[]
				#	break
			time.sleep(2.0)			
		ms_mainlog.info('All gaincal jobs are spawned.\n')
	else:
		ms_mainlog.info('No timestamp is left for calibration.\n')

	# Waiting for gaincal to finish # TODO : modify this part
	######################################
	if ref_time_freq==False:
		while True:
			gtables=glob.glob(basedir+'/caltables/'+str(ms_obsid)+'/'+basemsdir+'/*'+os.path.basename(temp_ref_ms).split('.ms')[0]+'*.cal')
			if len(gtables)>0:
				ms_mainlog.info('Gaintable found : '+gtables[0]+'\n')
				break
			else:
				time.sleep(2.0)
		ref_gaintable=[gtables[0]]
		#if len(calibrator_caltable)!=0:
		#	ref_gaintable=ref_gaintable+calibrator_caltable
		
	# Applying ref time solution
	############################
	ms_mainlog.info('Applying gain solution for reference time in all times and all channels.......\n')
	ms_mainlog.info('applycal(vis=\''+ref_averaged_msname+'\',gaintable='+str(ref_gaintable)+',applymode=\'calflag\',flagbackup=True)\n')
	applycal(vis=ref_averaged_msname,gaintable=ref_gaintable,applymode='calflag',flagbackup=True)
	ms_mainlog.info('applycal(vis=\''+ref_averaged_msname+'\',gaintable='+str(ref_gaintable)+',applymode=\'calflag\',flagbackup=True)\n')
	applycal(vis=ref_averaged_msname,gaintable=ref_gaintable,applymode='calflag',flagbackup=True)
	flaglist=flagmanager(vis=ref_averaged_msname,mode='list')
	flaglist_keys=list(flaglist.keys())
	flaglist_keys.remove('MS')
	if len(flaglist_keys)>0:
		for i in flaglist_keys:
			last_flagversion=flaglist[i]['name']
			# Restore the flag and delete the present flag version
			ms_mainlog.info('flagmanager(vis=\''+ref_averaged_msname+'\',mode=\'restore\',versionname=\''+str(last_flagversion)+'\',merge=\'replace\')\n')
			ms_mainlog.info('flagmanager(vis=\''+ref_averaged_msname+'\',mode=\'delete\',versionname=\''+str(last_flagversion)+'\')\n')
			flagmanager(vis=ref_averaged_msname,mode='restore',versionname=last_flagversion,merge='replace')
			flagmanager(vis=ref_averaged_msname,mode='delete',versionname=last_flagversion)

	#Deciding bandpass selfcal conditions
	#####################################
	AM=AccessMS(ref_averaged_msname)
	unflagged_channels=AM.get_unflag_chan(flagfrac=1) # Unflagged averaged channels
	if len(unflagged_channels)<=1:
		ms_mainlog.info('Only 1 unflagged channel is present. No bandpass is required.\n')
		do_bandpass==False
	elif inputs.quality_factor==0:
		ms_mainlog.info('Quality factor is 0. Skipping bandpass self calibration.\n')
		do_bandpass==False
	elif do_bandpass==True:
		ms_mainlog.info('Proceed for bandpass self-calibration considering same source model for '+str(skip_freq)+' kHz\n')
		if inputs.interactive==True:
			want_change=input('Want to change bandpass bandwidth? If yes type frequency bandwidth in kHz or press enter\n')
			if want_change!='':
				skip_freq=float(want_change)
				ms_mainlog.info('Now proceed for bandpass self-calibration considering same source model for modified bandwidth '+str(skip_freq)+' kHz\n')

	# Performing bandpass selfcal
	#############################
	num_bp=0
	if do_bandpass==True:
		if ref_time_freq==True:
			ref_touch_list=glob.glob(basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_*')
			if len(ref_touch_list)==0:
				os.system('touch '+basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_0')

		AM=AccessMS(ref_averaged_msname)
		nchan_per_bandpass=int(skip_freq/AM.calc_freqres())
		nchan=AM.get_num_channels()
		bandpass_cmd_list=[]
		bandpass_screen_list=[]
		bandpass_finished_file_list=[]
		ms_mainlog.info('Spliting reference time data for performing bandpass........\n')
		for i in range(0,nchan,nchan_per_bandpass):
			start_chan=i
			end_chan=i+nchan_per_bandpass-1
			if end_chan>nchan:
				end_chan=nchan-1
			ms_mainlog.info('Spliting ms of channel range : '+str(start_chan)+'~'+str(end_chan)+'\n')
			if ref_time_freq==True:
				splited_msname,splited_msdir=spliting_timechan(ref_averaged_msname,metafits,str(start_chan)+'~'+str(end_chan),ref_time,caltype='B',\
						ref_timechan=True,input_file=workdir+'/selfcal_inputs.py',datacolumn='corrected')
			else:
				splited_msname,splited_msdir=spliting_timechan(ref_averaged_msname,metafits,str(start_chan)+'~'+str(end_chan),ref_time,caltype='B',\
						ref_timechan=False,input_file=workdir+'/selfcal_inputs.py',datacolumn='corrected')
			touch_file=glob.glob(basedir+'/.Finished_bcal_'+str(ms_obsid)+'_'+basemsdir+'_'+os.path.basename(splited_msname)+'*')
			if len(touch_file)==0:
				cmd='run_bandpass_selfcal --msname '+splited_msname+' --metafits '+metafits+' --workdir '+splited_msdir+\
				' --verbose '+str(inputs.verbose)+' --interactive '+str(inputs.interactive)+' --fresh True'
				bandpass_cmd_list.append(cmd)
				bandpass_screen_list.append(str(ms_obsid)+'_'+os.path.basename(splited_msname).split('.ms')[0]+'_screen_B')
				bandpass_finished_file_list.append(basedir+'/.Finished_bcal_'+str(ms_obsid)+'_'+basemsdir+'_'+os.path.basename(splited_msname))
			if end_chan>=nchan:
				break
		finished_bandpass=False
		num_bp=len(bandpass_screen_list)
		bp_finish_list=glob.glob(basedir+'/.Finished_*bcal*'+str(ms_obsid)+'*'+basemsdir+'*')
		if len(bp_finish_list)==num_bp:
			ms_mainlog.info('Bandpass for channel blocks are finished.\n')
			bandpass_screen_list=[]
			bandpass_cmd_list=[]
			bandpass_finished_file_list=[]
			finished_bandpass=True
		else:
			ms_mainlog.info('Waiting for available casa instance.....\n')
		mpicmd=[]	
		mpicount=0
		while finished_bandpass==False:
			touch_file_list=glob.glob(basedir+'/.Finished_*cal*'+str(ms_obsid)+'*'+basemsdir+'*')
			if (len(touch_file_list)-touch_count)>0:
				if open_casa_instance>1:
					open_casa_instance-=(len(touch_file_list)-touch_count)
					ms_mainlog.info('New CASA instance available : '+str((len(touch_file_list)-touch_count))+'\n')		
			while len(bandpass_screen_list)!=0 and (casa_instance-open_casa_instance)>=1:
				screen_name=bandpass_screen_list[0]
				cmd=bandpass_cmd_list[0]
				finished_file=bandpass_finished_file_list[0]
				screen_batch_file=casa_instance_runner(cmd,basedir,screen_name,finished_file)
				touch_count+=1
				if mpi==1:
					screen_cmd='sh '+screen_batch_file
					os.system('screen -S '+screen_name+' -X quit')	
					time.sleep(0.5)
					os.system('screen -mdS '+screen_name)
					time.sleep(0.5)
					ms_mainlog.info('########################\n')
					ms_mainlog.info('Made Screen : '+screen_name+'\n')
					ms_mainlog.info('Command : '+cmd+'\n')
					os.system('screen -S '+screen_name+' -X stuff \"'+screen_cmd+'\n"')	
				elif mpi==0:
					mpicmd.append('-np 1 --bind-to core --map-by ppr:'+str(int(available_cpu_for_paircars))+':node:pe=4 -x OMP_NUM_THREADS='+\
							str(available_cpu_for_paircars)+' sh '+screen_batch_file+'\n')
					mpicmd.append('-np 1 sleep 1\n')
				open_casa_instance+=1
				cur_spawned_casa_instances+=1
				bandpass_screen_list.remove(screen_name)
				bandpass_cmd_list.remove(cmd)
				bandpass_finished_file_list.remove(finished_file)
				if open_casa_instance>=casa_instance or len(bandpass_screen_list)==0:
					ms_mainlog.info('Maximum casa instances spawned. Waiting for complete those jobs.\n')
					if mpi==0:
						mpicmd_file=basedir+'/'+basemsdir+'.bcal_mpicmd_'+str(mpicount)
						if os.path.exists(mpicmd_file):
							os.system('rm -rf '+mpicmd_file)
						mpifil=open(mpicmd_file,'w')
						ms_mainlog.info('MPI commands .....\n')
						for i in mpicmd:
							ms_mainlog.info(i)
						mpifil.writelines(mpicmd)
						mpifil.close()
						os.system('chmod a+rwx '+mpicmd_file)
						screen_cmd='mpirun --app '+mpicmd_file
						os.system('screen -S '+str(ms_obsid)+'_'+basemsdir+'_bcal'+str(mpicount)+' -X quit')	
						time.sleep(0.5)
						os.system('screen -mdS '+str(ms_obsid)+'_'+basemsdir+'_bcal_'+str(mpicount))
						time.sleep(0.5)
						ms_mainlog.info('########################\n')
						ms_mainlog.info('Made Screen : '+str(ms_obsid)+'_'+basemsdir+'_bcal_'+str(mpicount)+'\n')
						ms_mainlog.info('Command : '+screen_cmd+'\n')
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
			bp_finish_list=glob.glob(basedir+'/.Finished_*bcal*'+str(ms_obsid)+'*'+basemsdir+'*')
			if len(bandpass_cmd_list)==0:
				ms_mainlog.info('All bandpass tasks are spawned for all spectral slices.\n')
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
	prefix_cmds.append('logfiles=('+basedir+'/.Finished_bcal*'+str(ms_obsid)+'*'+str(basemsdir)+'*)\n')
	prefix_cmds.append('C=${#logfiles[@]}\n')
	prefix_cmds.append('\techo "Waiting : "$C\n')
	prefix_cmds.append('until [ $C -ge '+str(num_bp)+' ]\n')
	prefix_cmds.append('do\n') 
	prefix_cmds.append('\tsleep 5\n')
	prefix_cmds.append('\tshopt -s nullglob\n')
	prefix_cmds.append('\tlogfiles=('+basedir+'/.Finished_bcal*'+str(ms_obsid)+'*'+str(basemsdir)+'*)\n')
	prefix_cmds.append('\tC=${#logfiles[@]}\n')
	prefix_cmds.append('done\n')
	prefix_cmds.append('shopt -s nullglob\n')
	prefix_cmds.append('logfiles=('+basedir+'/bpcaltables/'+str(ms_obsid)+'/'+basemsdir+'/*.bcal)\n')
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
			ref_touch_list=glob.glob(basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_*')
			if len(ref_touch_list)==0:
				os.system('touch '+basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_0')
	
		AM=AccessMS(ref_averaged_msname)
		unflagged_channels=AM.get_unflag_chan(flagfrac=1) # Unflagged averaged channels
		if len(unflagged_channels)==0:
			ms_mainlog.info('No unflagged channel is present.\n')
			if os.path.exists(basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))==True:
				os.system('rm -rf '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
				os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
			else:
				os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
			if ref_time_freq==True:
				np.save(basedir+'/Ref_time_cal_record',np.array([1,0,0,0,0,0,0]))	
				os.system('touch '+basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_selfcalerror')	
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
		
		ms_mainlog.info('Polcal channel grid list : '+str(pol_channel_grid)+' for calibration per '+str(skip_freq_pol)+' MHz.\n')
		polcal_cmd_list=[]
		polcal_screen_list=[]
		polcal_finished_file_list=[]
		for i in pol_channel_grid:
			ms_mainlog.info('Spliting data for performing polarisation calibration of channel : '+str(i)+' and timerange : '+ref_time+'\n')
			if ref_time_freq==True:
				splited_msname,splited_msdir=spliting_timechan(ref_averaged_msname,metafits,str(i),ref_time,caltype='P',ref_timechan=True,\
											input_file=workdir+'/selfcal_inputs.py',datacolumn='data')
			else:
				splited_msname,splited_msdir=spliting_timechan(ref_averaged_msname,metafits,str(i),ref_time,caltype='P',ref_timechan=False,\
											input_file=workdir+'/selfcal_inputs.py',datacolumn='data')
			calstring=','.join(ref_gaintable)
			touch_file=glob.glob(basedir+'/.Finished_pcal_'+str(ms_obsid)+'_'+basemsdir+'_'+os.path.basename(splited_msname)+'*')
			if len(touch_file)==0:
				if len(ref_gaintable)!=0:
					cmd='run_pol_selfcal --msname '+splited_msname+' --metafits '+metafits+' --workdir '+splited_msdir+' --verbose '+str(inputs.verbose)+\
					' --interactive '+str(inputs.interactive)+' --fresh True --gaincal '+str(perform_leakcor)+' --caltables '+calstring+',\"$y\"'
					polcal_cmd_list.append(cmd)
					polcal_screen_list.append(str(ms_obsid)+'_'+os.path.basename(splited_msname).split('.ms')[0]+'_screen_P')
					polcal_finished_file_list.append(basedir+'/.Finished_pcal_'+str(ms_obsid)+'_'+basemsdir+'_'+os.path.basename(splited_msname))
		num_pol=len(polcal_screen_list)
		finished_polcal=False
		polcal_finish_list=glob.glob(basedir+'/.Finished_*pcal*'+str(ms_obsid)+'*'+basemsdir+'*')
		if len(polcal_finish_list)==num_pol:
			ms_mainlog.info('Polcal for all coarse channels have been finished.\b')
			polcal_screen_list=[]
			polcal_cmd_list=[]
			polcal_finished_file_list=[]
			finished_polcal=True
		mpicmd_file=basedir+'/'+basemsdir+'.polcal_mpicmd'
		ms_mainlog.info('Waiting for available casa instance......\n')
		mpicmd=[]	
		mpicount=0
		while finished_polcal==False:	
			touch_file_list=glob.glob(basedir+'/.Finished_*cal*'+str(ms_obsid)+'*'+basemsdir+'*')
			if (len(touch_file_list)-touch_count)>0:
				if open_casa_instance>1:
					open_casa_instance-=(len(touch_file_list)-touch_count)
					ms_mainlog.info('New CASA instance available : '+str((len(touch_file_list)-touch_count))+'\n')
			while len(polcal_screen_list)!=0 and (casa_instance-open_casa_instance)>=1:
				screen_name=polcal_screen_list[0]
				cmd=polcal_cmd_list[0]
				finished_file=polcal_finished_file_list[0]
				screen_batch_file=casa_instance_runner(cmd,basedir,screen_name,finished_file,prefix_cmds=prefix_cmds)
				touch_count+=1
				if mpi==1:
					screen_cmd='sh '+screen_batch_file
					os.system('screen -S '+screen_name+' -X quit')	
					time.sleep(0.5)
					os.system('screen -mdS '+screen_name)
					time.sleep(0.5)
					ms_mainlog.info('########################\n')
					ms_mainlog.info('Made Screen : '+screen_name+'\n')
					ms_mainlog.info('Command : '+cmd+'\n')
					os.system('screen -S '+screen_name+' -X stuff \"'+screen_cmd+'\n"')	
				elif mpi==0:
					mpicmd.append('-np 1 --bind-to core --map-by ppr:'+str(int(available_cpu_for_paircars))+':node:pe=4 -x OMP_NUM_THREADS='+\
							str(available_cpu_for_paircars)+' sh '+screen_batch_file+'\n')
					mpicmd.append('-np 1 sleep 1\n')
				open_casa_instance+=1
				cur_spawned_casa_instances+=1
				polcal_screen_list.remove(screen_name)
				polcal_cmd_list.remove(cmd)
				polcal_finished_file_list.remove(finished_file)
				if open_casa_instance>=casa_instance or len(polcal_screen_list)==0:
					ms_mainlog.info('Maximum casa instances spawned. Waiting for complete those jobs.\n')	
					if mpi==0:
						mpicmd_file=basedir+'/'+basemsdir+'.pcal_mpicmd_'+str(mpicount)
						if os.path.exists(mpicmd_file):
							os.system('rm -rf '+mpicmd_file)
						mpifil=open(mpicmd_file,'w')
						ms_mainlog.info('MPI commands .....\n')
						for i in mpicmd:
							ms_mainlog.info(i)
						mpifil.writelines(mpicmd)
						mpifil.close()
						os.system('chmod a+rwx '+mpicmd_file)
						screen_cmd='mpirun --app '+mpicmd_file
						os.system('screen -S '+str(ms_obsid)+'_'+basemsdir+'_polcal_'+str(mpicount)+' -X quit')	
						time.sleep(0.5)
						os.system('screen -mdS '+str(ms_obsid)+'_'+basemsdir+'_polcal_'+str(mpicount))
						time.sleep(0.5)
						ms_mainlog.info('########################\n')
						ms_mainlog.info('Made Screen : '+str(ms_obsid)+'_'+basemsdir+'_polcal_'+str(mpicount)+'\n')
						ms_mainlog.info('Command : '+screen_cmd+'\n')
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
		
		ms_mainlog.info('All calibration job spawned for ms : '+msname+'\n')
		ms_mainlog.info('#########################\n')
		os.system('rm -rf '+ref_timesliced_measurement_set+'* '+ref_timesliced_measurement_set+'.flagversions')
		os.system('rm -rf '+ref_averaged_msname+'* '+ref_timesliced_measurement_set+'.flagversions')
		os.system('rm -rf '+timesliced_measurement_set+'* '+ref_timesliced_measurement_set+'.flagversions')
		os.system('rm -rf '+averaged_msname+'* '+ref_timesliced_measurement_set+'.flagversions')
		del pol_channel_grid
		if ref_time_freq==True:
			np.save(basedir+'/Ref_time_cal_record',np.array([0,ref_time,ref_freq,ref_chan,cur_spawned_casa_instances,ref_freq_avg,ref_time_avg]))
			os.system('touch '+basedir+'/.ref_timechan_done_'+str(ms_obsid)+'_'+str(os.path.basename(msname))+'_0')
		if os.path.exists(basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(cur_spawned_casa_instances))==True:
			os.system('rm -rf '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(cur_spawned_casa_instances))
			os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(cur_spawned_casa_instances))
		else:
			os.system('touch '+basedir+'/.Finished_spawned_'+str(ms_obsid)+'_'+basemsdir+'_'+str(cur_spawned_casa_instances))
		return 0,ref_time,ref_chan,cur_spawned_casa_instances,ref_freq_avg,ref_time_avg

# Function to run the script stand alone from command line
if __name__=='__main__':
	import selfcal_inputs as inputs
	start_time=time.time()
	if options.chantime_msname==None:
		print ('No Measurement set is given.\n')
		os._exit(1)
	
	if options.metafits==None:
		print ('No metafits file is given.\n')
		os._exit()

	if options.caltables!=None:
		calibrator_caltable=str(options.caltables).split(',')
	else:
		calibrator_caltable=[]			
	ms_obsid=get_OBSID_from_metafits(options.metafits)
	# Logger initiating
	###################
	formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
	ms_mainlog = logging.getLogger('paircars_log_'+os.path.basename(str(options.chantime_msname)).split('.ms')[0])
	ms_mainlog.setLevel(logging.DEBUG)
	console=logging.StreamHandler(sys.stdout)
	console.setFormatter(formatter)
	ms_mainlog.addHandler(console)
	filehandle=logging.FileHandler(str(options.basedir)+'/PAIRCARS_mainlog_'+os.path.basename(str(options.chantime_msname)).split('.ms')[0]+'.log')
	filehandle.setFormatter(formatter)
	ms_mainlog.addHandler(filehandle)
	ms_mainlog.propagate = False	
	basemsdir=os.path.basename(options.chantime_msname).split('.ms')[0]
	try:
		ms_mainlog.info('Starting calibration for ms : '+str(options.chantime_msname)+'\n')		
		ms_mainlog.info('run_paircars_ms(\''+str(options.chantime_msname)+'\',\''+str(options.metafits)+'\',\''+str(options.workdir)+'\',ref_freq_avg='+str(options.ref_freq_avg)+\
						',ref_time_avg='+str(options.ref_time_avg)+',ref_time_freq='+str(options.ref_time_freq)+',do_bandpass='+str(options.do_bandpass)+\
						',do_polcal='+str(options.do_polcal)+',calatten='+str(options.calatten)+',num_threads='+str(options.num_threads)+\
						',calibrator_caltable='+str(calibrator_caltable)+')\n')
		result=run_paircars_ms(str(options.chantime_msname),str(options.metafits),str(options.workdir),ref_freq_avg=float(options.ref_freq_avg),ref_time_avg=float(options.ref_time_avg),\
				ref_time_freq=eval(str(options.ref_time_freq)),do_bandpass=eval(str(options.do_bandpass)),do_polcal=eval(str(options.do_polcal)),\
				calatten=float(options.calatten),num_threads=int(options.num_threads),calibrator_caltable=calibrator_caltable)
		result=list(result)
		result[0]=str(result[0])
		result[2]=int(result[2])
		result[3]=int(result[3])
		result[4]=float(result[4])
		result[5]=float(result[5])
		if os.path.isfile(str(options.basedir)+'/Ref_time_freq_slice_output.npy')==False and eval(str(options.ref_time_freq))==True:
			result.append(str(options.chantime_msname))
			result=np.array([result],dtype='object')
			np.save(str(options.basedir)+'/Ref_time_freq_slice_output',result)
		elif eval(str(options.ref_time_freq))==True:
			pre_result=np.load(str(options.basedir)+'/Ref_time_freq_slice_output.npy',allow_pickle=True)
			result.append(str(options.chantime_msname))
			result=np.array([result],dtype='object')
			result=np.append(pre_result,result,axis=0)
			np.save(str(options.basedir)+'/Ref_time_freq_slice_output',result)
		if os.path.isfile(str(options.basedir)+'/Nonref_time_freq_slice_output.npy')==False and eval(str(options.ref_time_freq))==False:
			result.append(str(options.chantime_msname))
			result=np.array([result],dtype='object')
			np.save(str(options.basedir)+'/Nonref_time_freq_slice_output',result)
		elif eval(str(options.ref_time_freq))==False:
			pre_result=np.load(str(options.basedir)+'/Nonref_time_freq_slice_output.npy',allow_pickle=True)
			result.append(str(options.chantime_msname))
			result=np.array([result],dtype='object')
			result=np.append(pre_result,result,axis=0)
			np.save(str(options.basedir)+'/Nonref_time_freq_slice_output',result)
		os.system('touch '+basedir+'/.Finished_runpaircars_'+str(ms_obsid)+'_'+basemsdir+'_'+str(0))
	except:
		os.system('touch '+basedir+'/.Finished_runpaircars_'+str(ms_obsid)+'_'+basemsdir+'_error')

