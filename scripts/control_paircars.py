'''
Code is written by Devojyoti Kansabanik , 28 Jan, 2021
'''
import os,sys
a=os.system('python3 validating_paircars_input.py\n')
if os.WEXITSTATUS(a)!=0:
	os._exit(1)
from paircars_inputs import basedir
from casatools import *
from casatasks import *
import logging,numpy as np,copy,glob,psutil,time
from paircars.basic_func import *
from paircars.access_ms import *
from paircars.decor import *
from paircars.flagger import *
from astropy.io import fits
from validating_paircars_input import download_metafits,get_OBSID

def spliting_timechan(msname,channel,timestamp,caltype='',ref_timechan=False,input_file='',datacolumn='corrected'):
	'''
	Function to split specific time and frequency slice and keep the necessary files in one directory to run the PAIRCARS
	Parameters:
	msname = Name of the source measurement set
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
	mjd_timestamps=md.timesforfield(0)
	freqlist=md.chanfreqs(0)
	md.close()
	if type(channel)==str:
		if '~' in channel:
			ch0=int(channel.split('~')[0])
			ch1=int(channel.split('~')[1])
			ch=int((ch0+ch1)/2.0)
			spw='0:'+channel
	else:
		spw='0:'+str(channel)
		ch=channel
	if ref_timechan==True:
		mainlog.info('Reference channel : '+str(channel)+' at '+str(freqlist[ch]/10**6)+' MHz and reference time : '+str(timestamp))
		# Spliting ref time chan ms
		###########################
		mainlog.info('Spliting reference time and channel..............\n')	
		if os.path.isdir(cwd+'/reftimechan.ms')==True:
			os.system('rm -rf '+cwd+'/reftimechan.ms')
		split(vis=msname,outputvis=cwd+'/reftimechan.ms',datacolumn=datacolumn,spw=spw,timerange=timestamp)
		mainlog.info('split(vis=\''+msname+'\',outputvis=\''+cwd+'/reftimechan.ms\',datacolumn=\''+datacolumn+'\',spw=\''+spw+'\',timerange=\''+timestamp+'\')\n')
		ref_timechan_ms=splited_ms_rename(cwd+'/reftimechan.ms',ref_time_chan=True,change_msname=True)
		ref_timechan_dir=cwd+'/'+os.path.basename(ref_timechan_ms).split('.ms')[0]+'_'+caltype
		if os.path.isdir(ref_timechan_dir)==False:
			os.makedirs(ref_timechan_dir) # Making ref time chan directory
		os.system('cp -r '+input_file+' '+ref_timechan_dir)
		os.system('mv '+ref_timechan_ms+' '+ref_timechan_dir+'/'+os.path.basename(ref_timechan_ms))
		return ref_timechan_dir+'/'+os.path.basename(ref_timechan_ms),ref_timechan_dir
	else:
		# Spliting specific time chan ms
		################################
		mainlog.info('Spliting measurement set for time : '+timestamp+' and frequency : '+str(freqlist[ch]/10**6)+' MHz ............\n')
		if os.path.isdir(cwd+'/timechan.ms')==True:
			os.system('rm -rf '+cwd+'/timechan.ms')
		split(vis=msname,outputvis=cwd+'/timechan.ms',datacolumn=datacolumn,spw=spw,timerange=timestamp)
		mainlog.info('split(vis=\''+msname+'\',outputvis=\''+cwd+'/timechan.ms\',datacolumn=\''+datacolumn+'\',spw=\''+spw+'\',timerange=\''+timestamp+'\')\n')
		timechan_ms=splited_ms_rename(cwd+'/timechan.ms',ref_time_chan=False,change_msname=True)
		timechan_dir=cwd+'/'+os.path.basename(timechan_ms).split('.ms')[0]+'_'+caltype
		if os.path.isdir(timechan_dir)==False:
			os.makedirs(timechan_dir) # Making ref time chan directory
		os.system('cp -r '+input_file+' '+timechan_dir)
		os.system('mv '+timechan_ms+' '+timechan_dir+'/'+os.path.basename(timechan_ms))
		return timechan_dir+'/'+os.path.basename(timechan_ms),timechan_dir

def casa_instance_runner(cmd,screen_name):
	'''
	Function to run a casa instance
	Parameters:
	cmd = Command to run
	screen_name = Name of the screen
	'''
	if os.path.isfile(inputs.basedir+'/.mpi_enabled'):
		cmd='mpirun -np 1 -x OMP_NUM_THREADS=3 '+cmd
	os.system('screen -S '+screen_name+' -X quit')	
	time.sleep(0.5)
	os.system('screen -mdS '+screen_name)
	time.sleep(0.5)
	mainlog.info('########################\n')
	mainlog.info('Made Screen : '+screen_name+'\n')
	mainlog.info('Command : '+cmd+'\n')
	os.system('screen -S '+screen_name+' -X stuff \"'+cmd+'\n"')	
	return screen_name

def run_paircars_ms(msname,metafits,workdir,ref_time_freq=False,do_bandpass=False,do_polcal=False,calibrator_caltable=[]): #TODO: XY phasecal
	'''
	Function to run paircars on a measurement set
	Parameters:
	msname = Name of the measurement set
	metafits = Name of the metafits file
	workdir = Name of the working directory
	ref_time_freq = False, reference time frequency ms
	do_bandpass = False, perform bandpass or not
	do_polcal = False, perform polcal or not
	'''
	if workdir[-1]=='/':
		workdir=workdir[:-1]

	if os.path.isdir(workdir)==False:
		os.makedirs(workdir)

	cwd=os.getcwd()
	os.chdir(workdir)

	# Validating timerange
	######################
	AM=AccessMS(msname)
	if len(inputs.timerange)!=0:
		start,end,start_mjdsec,end_mjdsec=AM.get_scan_startend_time()
		new_timerange=[]
		times=inputs.timerange.split(',') # Spliting timeblocks separated by ','
		for i in range(len(times)):
			t=times[i].split('~')
			mjd_list=[]
			for timestamp in t:
				try:
					a=timestamp_to_mjdsec(timestamp,format=2)
					mjd_list.append(a)
				except:
					mainlog.error('Timestamp format is wrong. Removing timestamp '+str(i)+'\n')
					break
			if len(mjd_list)!=0:
				mjd_list=sorted(mjd_list)
				if len(mjd_list)==1:
					if mjd_list[0]>end_mjdsec or mjd_list[0]<start_mjdsec:
						mainlog.info('Time range is not in ms, removing timestamp '+str(i))
					elif mjd_list[0]<start_mjdsec:
						mainlog.info('Start time is less than observation start time. Shifted it to observation start time.\n')
						mjd_list[0]=start_mjdsec
					elif mjd_list[0]>end_mjdsec:
						mainlog.info('End time is greater than observation start time. Shifted it to observation end time.\n')
						mjd_list[0]=end_mjdsec
					new_timerange.append(mjdsec_to_timestamp(mjd_list[0],includedate=True,format=0))
				else:
					if mjd_list[0]>end_mjdsec or mjd_list[1]<start_mjdsec:
						mainlog.info('Time range is not in ms, removing timestamp '+str(i))
					elif mjd_list[0]<start_mjdsec:
						mainlog.info('Start time is less than observation start time. Shifted it to observation start time.\n')
						mjd_list[0]=start_mjdsec
					elif mjd_list[1]>end_mjdsec:
						mainlog.info('End time is greater than observation start time. Shifted it to observation end time.\n')
						mjd_list[1]=end_mjdsec
					new_timerange.append(mjdsec_to_timestamp(mjd_list[0],includedate=True,format=0)+'~'+mjdsec_to_timestamp(mjd_list[1],includedate=True,format=0))

		if len(new_timerange)!=0:
			new_timerange=','.join(new_timerange)
			mainlog.info('New time range : '+str(new_timerange)+'\n')
		else:
			mainlog.info('Given time stamps are not in ms. Continue with all times.\n')
			new_timerange=''
	else:
		new_timerange=''

	# Validating channel ranges
	###########################
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
			mainlog.info('New channel range : '+str(new_chan_list)+'\n')
		else:
			mainlog.info('Given channels are not in ms. Continue with all channels.\n')
			new_chan_list=''
	else:
		new_chan_list=''

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
	mainlog.info('###############################################\n')

	# Flagging coarse channel edges and center
	mainlog.info('Flagging edges and center of coarese chnnels\n')
	mainlog.info('flag_MWA_coarse(\''+msname+'\',edgewidth=160,do_flag=True)\n')
	good_channels,channel_per_coarse=flag_MWA_coarse(msname,edgewidth=160,do_flag=True)
	unflag_channels=np.array(AM.get_unflag_chan(flagfrac=1))
	new_chan_list=np.array(new_chan_list.split(','),dtype='int')
	new_unflag_chan_list=np.intersect1d(unflag_channels,new_chan_list)
	start_chan=np.min(new_unflag_chan_list)
	end_chan=np.max(new_unflag_chan_list)

	# Spliting chan and timeranges
	##############################
	mainlog.info('Spliting channel range : '+str(start_chan)+'~'+str(end_chan)+' and timerange : '+str(new_timerange)+'\n')	
	if os.path.isdir(workdir+'/'+os.path.basename(msname).split('.ms')[0]+'_chantimesliced.ms')==True:
		os.system('rm -rf '+workdir+'/'+os.path.basename(msname).split('.ms')[0]+'_chantimesliced.ms')
	mainlog.info('split(vis=\''+msname+'\',outputvis=\''+workdir+'/'+os.path.basename(msname).split('.ms')[0]+'_chantimesliced.ms\',timerange=\''+new_timerange+\
						',spw=\'0:'+str(start_chan)+'~'+str(end_chan)+'\',datacolumn=\'data\')\n')
	split(vis=msname,outputvis=workdir+'/'+os.path.basename(msname).split('.ms')[0]+'_chantimesliced.ms',timerange=new_timerange,\
						spw='0:'+str(start_chan)+'~'+str(end_chan),datacolumn='data')
	timesliced_measurement_set=workdir+'/'+os.path.basename(msname).split('.ms')[0]+'_chantimesliced.ms'
	os.system('rm -rf casa*.log')	# Time slices list
	
	##################
	md=msmetadata()
	md.open(timesliced_measurement_set)
	mjd_timestamps=md.timesforfield(0)
	md.close()
	timestamps=[mjdsec_to_timestamp(mjdsec,includedate=True,format=0) for mjdsec in mjd_timestamps]
	
	if freqres<160:
		freq_avg=160
	else:
		freq_avg=freqres
	time_avg=timeres
	freq_averaging_count=0
	time_averaging_count=0
	
	# Removing .Finished files if error occured
	###########################################
	touch_file_list=glob.glob(inputs.basedir+'/.Finished_*cal_*')
	if len(touch_file_list)!=0:
		for t in touch_file_list:
			msg=t.split('_')[-1]
			if msg=='error':
				os.system('rm -rf '+t)
			else:
				msg=int(msg)
				if msg>100:
					msg-=100
				if msg!=0 and mgs!=8 and msg!=9:
					os.system('rm -rf '+t)

	# In this while loop we are checking whether the present time and frequency averaging is enough to start the self calibration. If it has it will leave the loop and go for selfcal
	while os.path.isfile('.ref_timechan_done')==False and ref_time_freq==True:  	
		mainlog.info('Choosing averaging frequency width : '+str(freq_avg)+' kHz, averaging temporal width : '+str(time_avg)+' s, Skip frequency : '\
				+str(skip_freq)+' kHz, Skip time : '+str(skip_time)+' s\n')

		# Averaging measurement set
		###########################
		if freq_avg>AM.calc_freqres() or time_avg>AM.calc_timeres():
			chan_width=int(freq_avg/AM.calc_freqres())
			previous_averaged_ms=glob.glob(timesliced_measurement_set.split('.ms')[0]+'*averaged.ms')
			for pams in previous_averaged_ms:
				os.system('rm -rf '+pams)
			averaged_msname=timesliced_measurement_set.split('.ms')[0]+'_'+str(freq_avg)+'kHz_'+str(time_avg)+'s_averaged.ms'
			if os.path.isdir(averaged_msname)==True:
				os.system('rm -rf '+averaged_msname)	
			mainlog.info('Avearging measurement width frequency average :'+str(freq_avg)+' kHz, temporal average :'+str(time_avg)+'s\n')
			split(vis=timesliced_measurement_set,outputvis=averaged_msname,width=chan_width,timebin=str(time_avg)+'s',datacolumn='data')
			mainlog.info('split(vis=\''+timesliced_measurement_set+'\',outputvis=\''+averaged_msname+'\',width='+str(chan_width)+',timebin=\''+str(time_avg)+'s\',datacolumn=\'data\')\n')
		else:
			mainlog.info('No averaging is required.\n')
			averaged_msname=timesliced_measurement_set

		AM=AccessMS(averaged_msname)
		# Making reference time and channel and time frequency grid
		##########################################################
		unflagged_channels=AM.get_unflag_chan(flagfrac=1) # Unflagged averaged channels
		if unflagged_channels==0:
			mainlog.info('No unflagged channel is present.\n')
			return 0,0,0

		timestamps=[mjdsec_to_timestamp(mjdsec,includedate=True,format=0) for mjdsec in mjd_timestamps]

		skip_channel=int(skip_freq/AM.calc_freqres())
		skip_timestamp=int(skip_time/AM.calc_timeres())

		channel_grid=[]
		time_grid=[]

		if skip_channel>len(unflagged_channels):
			for i in range(0,len(unflagged_channels),skip_channel):
				channel_grid.append(unflagged_channels[i])
		else:
			channels_grid=unflagged_channels
		if skip_timestamp>len(timestamps):
			for i in range(0,len(timestamps),skip_timestamp):
				time_grid.append(timestamps[i])
		else:
			time_grid=timestamps

		channel_grid_copy=copy.deepcopy(channel_grid)
		time_grid_copy=copy.deepcopy(time_grid)

		mainlog.info('Channel grid list : '+str(channel_grid)+'\n')
		mainlog.info('Timestamp grid list : '+str(time_grid)+'\n')
		
		ref_chan=channel_grid[int(len(channel_grid)/2)]
		ref_time=time_grid[int(len(time_grid)/2)]

		# Spliting ref time chan ms
		###########################
		ref_timechan_ms,ref_timechan_dir=spliting_timechan(averaged_msname,ref_chan,ref_time,caltype='G',ref_timechan=True,input_file=workdir+'/selfcal_inputs.py',datacolumn='data')
		cur_workdir=ref_timechan_dir
		# Run selfcal
		while True:
			try:
				touch_file_list=glob.glob(inputs.basedir+'/.Finished_gcal_'+os.path.basename(ref_timechan_ms)+'_*')
				if len(touch_file_list)!=0:
					for t in touch_file_list:
						msg=t.split('_')[-1]
						if msg=='error':
							os.system('rm -rf '+t)
						else:
							msg=int(msg)
							if msg>100:
								msg-=100
							if msg!=0 and mgs!=8 and msg!=9:
								os.system('rm -rf '+t)
				if len(calibrator_caltable)!=0:
					calstring=','.join(calibrator_caltable)
					cmd='run_intensity_selfcal --msname '+ref_timechan_ms+' --metafits '+metafits+' --workdir '+cur_workdir+' --dopoint True --verbose '+str(inputs.verbose)\
						+' --interactive '+str(inputs.interactive)+' --fresh True --caltables '+calstring 
				else:
					cmd='run_intensity_selfcal --msname '+ref_timechan_ms+' --metafits '+metafits+' --workdir '+cur_workdir+' --dopoint True --verbose '+str(inputs.verbose)\
						+' --interactive '+str(inputs.interactive)+' --fresh True'
				screen_name=os.path.basename(ref_timechan_ms).split('.ms')[0]+'_screen_refG'
				result=casa_instance_runner(cmd,screen_name)
				mainlog.info('Self calibration for ms : '+ref_timechan_ms+' is spawned in screen : '+result+'\n')
				mainlog.info('Waiting to finish self calibration for reference time frequency ms :'+ref_timechan_ms+'................\n') 
				while True:
					time.sleep(2)
					touch_file_list=glob.glob(inputs.basedir+'/.Finished_gcal_'+os.path.basename(ref_timechan_ms)+'_*')
					msg=int(touch_file_list[0].split('_')[-1])
					break	
			except Exception as e:
				mainlog.error('Error occured :'+str(e)+'\n')
				mainlog.error('Error in running selfcal. Exiting PAIRCARS.....\n')
				os.system('rm -rf '+basedir+'/.paircars_running')
				os.system('touch '+basedir+'/.paircars_failed')
				return 0,0,0
			if msg>=100:
				msg=msg-100
			if msg==10 and time_averaging_count<1 and freq_averaging_count<1:  # Checking for selfcal SNR and increasing the time and frequency averaging if required
				mainlog.error('SNR for self calibration is not sufficent.\n')
				selfcal_snr=float(np.load(inputs.basedir+'/selfcal_snr.npy'))
				new_time_avg=time_avg*np.sqrt(selfcal_snr/min_selfcal_snr)*1.5
				if new_time_avg<=skip_time and time_averaging_count<1:
					time_avg=new_time_avg
					mainlog.info('Increasing time averaging to '+str(time_avg)+'s\n')
					time_averaging_count+=1
					break
				elif freq_averaging_count<1:
					mainlog.info('Time averaging reached skip time limit. Increasing frequency averaging.')
					new_freq_avg=freq_avg*np.sqrt(selfcal_snr/min_selfcal_snr)
					if new_freq_avg<=skip_freq:
						freq_avg=new_freq_avg
						mainlog.info('Increasing frequency averaging to '+str(freq_avg)+'kHz\n')
						freq_averaging_count+=1
						break
				else:
					mainlog.info('Both time and frequency averaging has been tried. Still SNR is not sufficient. Trying with other time frequency.\n')
					continue
			elif msg!=0:
				mainlog.info('Message : '+error_msgs(100)+','+error_msgs(msg)+'\n')
				if inputs.verbose==False:
					mainlog.info('Removing the directory : '+ref_timechan_dir)
					os.system('rm -rf '+ref_timechan_dir)
				if len(time_grid)!=0:
					mainlog.info('Removing timestamp : '+str(ref_time)+' from time grid.\n')
					time_grid.remove(ref_time)
					ref_time=time_grid[int(len(time_grid)/2)]
					mainlog.info('Trying for new timestamp :'+str(ref_time)+'\n')
					if inputs.verbose==False:
						os.system('rm -rf '+ref_timechan_dir)
					ref_timechan_ms,ref_timechan_dir=spliting_timechan(averaged_msname,ref_chan,ref_time,caltype='G',ref_timechan=True,input_file=workfir\
									+'/selfcal_inputs.py',datacolumn='data')
					cur_workdir=ref_timechan_dir
					continue
				elif len(channel_grid)!=0:
					time_grid=copy.deepcopy(time_grid_copy)
					channel_grid.remove(ref_chan)
					ref_chan=channel_grid[int(len(channel_grid)/2)]
					ref_time=time_grid[int(len(time_grid)/2)]
					mainlog.info('Trying with new channel : '+str(ref_chan)+'\n')
					if inputs.verbose==False:
						os.system('rm -rf '+ref_timechan_dir)
					ref_timechan_ms,ref_timechan_dir=spliting_timechan(averaged_msname,ref_chan,ref_time,caltype='G',ref_timechan=True,input_file=\
								workfir+'/selfcal_inputs.py',datacolumn='data')
					cur_workdir=ref_timechan_dir
					continue
				else:
					mainlog.info('Reference imaging has been tried over full measurement set. No good starting point is found. Exiting PAIRCARS...\n')	
					os.system('rm -rf '+basedir+'/.paircars_running')
					os.system('touch '+basedir+'/.paircars_failed')
					return 0,0,0
			elif msg==0:
				os.system('touch .ref_timechan_done')	
				mainlog.info('Reference time frequency calibration done.\n')
				time_grid.remove(ref_time)
				channel_grid.remove(ref_chan)	
				break

	unflagged_channels=AM.get_unflag_chan(flagfrac=1) # Unflagged averaged channels
	if unflagged_channels==0:
		mainlog.info('No unflagged channel is present.\n')
		return 0,0,0
	
	# If not the reference time frequency ms, making time and channel and time frequency grid
	#########################################################################################
	if ref_time_freq==False or (os.path.isfile('.ref_timechan_done')==False):
		timestamps=[mjdsec_to_timestamp(mjdsec,includedate=True,format=0) for mjdsec in mjd_timestamps]

		skip_channel=int(skip_freq/AM.calc_freqres())
		skip_timestamp=int(skip_time/AM.calc_timeres())

		channel_grid=[]
		time_grid=[]

		if skip_channel>len(unflagged_channels):
			for i in range(0,len(unflagged_channels),skip_channel):
				channel_grid.append(unflagged_channels[i])
		else:
			channels_grid=unflagged_channels
		if skip_timestamp>len(timestamps):
			for i in range(0,len(timestamps),skip_timestamp):
				time_grid.append(timestamps[i])
		else:
			time_grid=timestamps

		channel_grid_copy=copy.deepcopy(channel_grid)
		time_grid_copy=copy.deepcopy(time_grid)

		mainlog.info('Channel grid list : '+str(channel_grid)+'\n')
		mainlog.info('Timestamp grid list : '+str(time_grid)+'\n')
		
		ref_chan=channel_grid[int(len(channel_grid)/2)]
		ref_time=time_grid[int(len(time_grid)/2)]


	#Deciding bandpass selfcal conditions
	#####################################
	if len(unflagged_channels)<=1:
		mainlog.info('Only 1 unflagged channel is present. No bandpass is required.\n')
	elif inputs.quality_factor==0:
		mainlog.info('Quality factor is 0. Skipping bandpass self calibration.\n')
		do_bandpass==False
	elif len(calibrator_caltable)!=0 and do_bandpass==True:
		mainlog.warning('Bandpass calibration is already applied using calibrator observation.\n')
		if inpurs.interactive==True:
			want_bandpass=input('Do you still want to perform bandpass selfcal? Y/N\n')
			if want_bandpass=='N':
				do_bandpass=False
	elif do_bandpass==True:
		mainlog.info('Proceed for bandpass self-calibration considering same source model for '+str(skip_freq)+' kHz\n')
		if inputs.interactive==True:
			want_change=input('Want to change bandpass bandwidth? If yes type frequency bandwidth in kHz or press enter\n')
			if want_change!='':
				skip_freq=float(want_change)
				mainlog.info('Now proceed for bandpass self-calibration considering same source model for modified bandwidth '+str(skip_freq)+' kHz\n')

	ms_obsid=get_OBSID(metafits)
	# Making ref time freq gaintable list
	#####################################
	if ref_time_freq==True:
		ref_timechan_caltable=glob.glob(inputs.basedir+'/caltables/'+str(ms_obsid)+'/*ref*.cal')
		ref_gaintable=copy.deepcopy(calibrator_caltable)
		ref_gaintable.append(ref_timechan_caltable)
	else:
		ref_gaintable=copy.deepcopy(calibrator_caltable)
		averaged_msname=timesliced_measurement_set

	# Estimating total casa instances
	#################################
	total_available_cpu=psutil.cpu_count()-(psutil.cpu_count()*psutil.cpu_percent()/100.0)
	available_cpu_for_paircars=int(total_available_cpu*inputs.cpu_frac)
	casa_instance=int(available_cpu_for_paircars/3)
	open_casa_instance=0
	spawned_casa_instances=0
	touch_count=0
	mainlog.info('Available cpus for P-AIRCARS: '+str(available_cpu_for_paircars)+'\n')
	mainlog.info('Total number of CASA instances : '+str(casa_instance)+'\n')

	# Spliting gaincal measurement set
	##################################
	mainlog.info('Spliting reference chan data for performing gaincal.\n')
	gaincal_cmd_list=[]
	gaincal_screen_list=[]
	for timestamp in time_grid:
		msname,msdir=spliting_timechan(averaged_msname,str(ref_chan),timestamp,caltype='G',ref_timechan=False,\
										input_file=workdir+'/selfcal_inputs.py',datacolumn='data')

		calstring=','.join(ref_gaintable)
		cmd='run_intensity_selfcal --msname '+msname+' --metafits '+metafits+' --workdir '+msdir+' --dopoint True --verbose '+str(inputs.verbose)\
				+' --interactive '+str(inputs.interactive)+' --fresh True --caltables '+calstring
		gaincal_cmd_list.append(cmd)
		gaincal_screen_list.append(os.path.basename(msname).split('.ms')[0]+'_screen_G')

	while len(gaincal_cmd_list)!=0:  # Loop while all gaincal cmds are spawned
		touch_file_list=glob.glob(inputs.basedir+'/.Finished_*cal_*'+os.path.basename(msname)+'*')
		if (len(touch_file_list)-touch_count)!=0:
			if open_casa_instance!=0:
				open_casa_instance-=(len(touch_file_list)-touch_count)
			touch_count=len(touch_file_list)
		if open_casa_instance<casa_instance: # If casa instance available
			while True:
				screen_name=gaincal_screen_list[0]
				cmd=gaincal_cmd_list[0]
				result=casa_instance_runner(cmd,screen_name)
				open_casa_instance+=1
				spawned_casa_instances+=1
				gaincal_screen_list.remove(screen_name)
				gaincal_cmd_list.remove(cmd)
				if open_casa_instance>casa_instance:
					mainlog.info('Maximum casa instances spawned. Waiting for complete those jobs.\n')
					break	
		time.sleep(2.0)	
	mainlog.info('All gaincal jobs are spawned.\n')

	# Applying ref time solution
	############################

	index=time_grid.index(ref_time)
	if len(time_grid)>=10:
		timerange=','.join(time_grid[index-5:index+5])
	else:
		timerange=','.join(time_grid)
	mainlog.info('Applying gain solution for reference time in the timerange : '+str(timerange)+' in all channels.......\n')
	mainlog.info('applycal(vis=\''+averaged_msname+'\',gaintable='+str(ref_gaintable)+',timerange=\''+timerange+'\',applymode=\'calflag\',flagbackup=False)\n')
	applycal(vis=averaged_msname,gaintable=ref_gaintable,timerange=timerange,applymode='calflag',flagbackup=False)

	# Performing bandpass selfcal
	#############################
	if do_bandpass==True:
		AM=AccessMS(averaged_msname)
		nchan_per_bandpass=int(skip_freq/AM.calc_freq())
		nchan=AM.get_num_channels()
		bandpass_cmd_list=[]
		bandpass_screen_list=[]
		mainlog.info('Spliting reference time data for performing bandpass........\n')
		index=time_grid.index(ref_time)
		if len(time_grid)>=10:
			timerange=','.join(time_grid[index-5:index+5])
		else:
			timerange=','.join(time_grid)
		for i in range(0,nchan,nchan_per_bandpass):
			start_chan=i
			end_chan=i+nchan_per_bandpass
			if (nchan-end_chan)<int(nchan_per_bandpass/2):
				end_chan=nchan
			mainlog.info('Spliting ms of channel range : '+str(start_chan)+'~'+str(end_chan)+'\n')
			if ref_time_freq==True:
				msname,msdir=spliting_timechan(averaged_msname,str(start_chan)+'~'+str(end_chan),timerange,caltype='B',\
						ref_timechan=True,input_file=workdir+'/selfcal_inputs.py',datacolumn='corrected')
			else:
				msname,msdir=spliting_timechan(averaged_msname,str(start_chan)+'~'+str(end_chan),timerange,caltype='B',\
						ref_timechan=False,input_file=workdir+'/selfcal_inputs.py',datacolumn='corrected')
			cmd='run_bandpass_selfcal --msname '+msname+' --metafits '+metafits+' --workdir '+msdir+' --verbose '+str(inputs.verbose)+' --interactive '\
					+str(inputs.interactive)+' --fresh True'
			bandpass_cmd_list.append(cmd)
			bandpass_screen_list.append(os.path.basename(averaged_msname).split('.ms')[0]+'_screen_B')
			if end_chan>=nchan:
				break
		finished_bandpass=False
		while finished_bandpass==False:
			touch_file_list=glob.glob(inputs.basedir+'/.Finished_*cal_*'+os.path.basename(msname)+'*')
			if (len(touch_file_list)-touch_count)!=0:
				if open_casa_instance!=0:
					open_casa_instance-=(len(touch_file_list)-touch_count)
				touch_count=len(touch_file_list)
			if open_casa_instance<casa_instance:
				while True:
					screen_name=bandpass_screen_list[0]
					cmd=bandpass_cmd_list[0]
					result=casa_instance_runner(cmd,screen_name)
					open_casa_instance+=1
					spwaned_casa_instances+=1
					bandpass_screen_list.remove(screen_name)
					bandpass_cmd_list.remove(cmd)
					if open_casa_instance>casa_instance:
						mainlog.info('Maximum casa instances spawned. Waiting for complete those jobs.\n')
						break
			time.sleep(2.0)
			bp_finish_list=glob.glob(inputs.basedir+'/.Finished_*bcal_*'+os.path.basename(msname)+'*')
			if len(bp_finish_list)==0:
				finished_bandpass=True
				break	

		bpcaltable=glob.glob(inputs.basedir+'/bpcaltables/'+str(ms_obsid)+'/*ref*.cal')+ref_gaintable
		index=time_grid.index(ref_time)
		if len(time_grid)>=10:
			timerange=','.join(time_grid[index-5:index+5])
		else:
			timerange=','.join(time_grid)
		mainlog.info('Applying bandpass solution for timerange : '+str(timerange)+' in all channels.......\n')
		mainlog.info('applycal(vis=\''+averaged_msname+'\',gaintable='+str(bpcaltable)+',timerange=\''+timerange+'\',applymode=\'calflag\',flagbackup=False)\n')
		applycal(vis=averaged_msname,gaintable=bpcaltable,timerange=timerange,applymode='calflag',flagbackup=False)
				
	# Spliting gain calibrated reference time channnel measurement set for polarisation calibration
	###############################################################################################
	if do_polcal==True and ((do_bandpass==True and finished_bandpass==True) or do_bandpass==False):
		good_avg_channels,channel_per_coarse=flag_MWA_coarse(averaged_msname,edgewidth=160,do_flag=False)
		polcal_channels=channel_per_coarse.values()
		polcal_cmd_list=[]
		polcal_screen_list=[]
		index=time_grid.index(ref_time)
		if len(time_grid)>=10:
			timerange=','.join(time_grid[index-5:index+5])
		else:
			timerange=','.join(time_grid)
		for i in polcal_channels:
			mainlog.info('Spliting data for performing polarisation calibration of channel : '+str(i)+' and timerange : '+timerange+'\n')
			if ref_time_freq==True:
				msname,msdir=spliting_timechan(averaged_msname,str(i),timerange,caltype='P',ref_timechan=True,\
											input_file=workdir+'/selfcal_inputs.py',datacolumn='corrected')
			else:
				msname,msdir=spliting_timechan(averaged_msname,str(i),timerange,caltype='P',ref_timechan=False,\
											input_file=workdir+'/selfcal_inputs.py',datacolumn='corrected')
			if len(calibrator_caltable)!=0:
				cmd='run_pol_selfcal --msname '+msname+' --metafits '+metafits+' --workdir '+msdir+' --verbose '+str(inputs.verbose)+\
				' --interactive '+str(inputs.interactive)+' --fresh True --gaincal False'
			else:
				cmd='run_pol_selfcal --msname '+msname+' --metafits '+metafits+' --workdir '+msdir+' --verbose '+str(inputs.verbose)+\
				' --interactive '+str(inputs.interactive)+' --fresh True --gaincal True'
			polcal_cmd_list.append(cmd)
			polcal_screen_list.append(os.path.basename(msname).split('.ms')[0]+'_screen_P')
		finished_polcal=False
		while finished_polcal==False:	
			touch_file_list=glob.glob(inputs.basedir+'/.Finished_*cal_*'+os.path.basename(msname)+'*')
			if (len(touch_file_list)-touch_count)!=0:
				if open_casa_instance!=0:
					open_casa_instance-=(len(touch_file_list)-touch_count)
				touch_count=len(touch_file_list)
			if open_casa_instance<casa_instance: # If casa instance available
				while True:
					screen_name=polcal_screen_list[0]
					cmd=polcal_cmd_list[0]
					result=casa_instance_runner(cmd,screen_name)
					open_casa_instance+=1
					spawned_casa_instances+=1
					polcal_screen_list.remove(screen_name)
					polcal_cmd_list.remove(cmd)
					if open_casa_instance>casa_instance:
						mainlog.info('Maximum casa instances spawned. Waiting for complete those jobs.\n')
						break			
			polcal_finish_list=glob.glob(inputs.basedir+'/.Finished_*pcal_*'+os.path.basename(msname)+'*')
			if len(polcal_finish_list)==0:
				finished_polcal=True
				break
	
		if ref_time_freq==False:
			mainlog.info('All calibration job spawned for ms : '+msname+'\n')
			mainlog.info('#########################\n')
			return ref_time,ref_chan,spawned_casa_instances
		else:
			if (do_bandpass==True and do_polcal==True and len(bandpass_cmd_list)==0 and len(polcal_cmd_list)==0 and len(gaincal_cmd_list)==0) or \
			(do_bandpass==True and do_polcal==False and len(bandpass_cmd_list)==0 and len(gaincal_cmd_list)==0) or (do_bandpass==False and do_polcal==False and len(gaincal_cmd_list)==0):
				mainlog.info('All calibration finished for reference time frequency ms : '+msname+'\n')
				mainlog.info('#########################\n')
				if clear_screen==True:
					os.system('screen -ls | tail -n +2 | head -n -2 | awk \'{print $1}\'| xargs -I{} screen -S {} -X quit')
				return ref_time,ref_chan,open_casa_instances

def managing_caldatabase(msname,ref_time,gaincal_modedir,bandpass_modeldir,polcal_caldir,localdatabase):
	'''
	Function to manager calibration database
	msname : Averaged msname
	ref_time : Reference timestamp
	gaincal_modeldir = Model directory for gaincal
	bandpass_modeldir = Model directory for bandpass
	polcal_caldir = Polarisation caltable directory
	localdatabase = Local database directory
	'''
	OBSID=get_OBSID(msname)
	if localdatabase=='' or os.path.isdir(localdatabase)==False:
		mainlog.error('Local data base not found. Making local database at basedir.\n')
		localdatabase=inputs.basedir+'/localdatabase/'+str(OBSID)
		if os.path.isdir(localdatbase)==False:
			os.makedirs(localdatbase)
	else:
		mainlog.info('Local data base is at : '+localdatbase+'\n')
		localdatabase=localdatabase+'/'+str(OBSID)
		if os.path.isdir(localdatbase)==False:
			os.makedirs(localdatbase)
	if gaincal_modeldir=='' or os.path.isdir(gaincal_modeldir)==False or len(glob.glob(gaincal_modeldir+'/*.model'))==0:
		mainlog.info('No models available.\n')
		return 1
	else:
		AM=AccessMS(msname)
		freqs=AM.get_freqs()/10**6
		coarse_chan_0=int(freqs[0]/1.28)
		coarse_chan_1=int(freqs[-1]/1.28)
		caltable_name=str(OBSID)+'_'+str(coarse_chan_0)+'_'+str(coarse_chan_1)+'.gcal'
		model_list=glob.glob(gaincal_modeldir+'/*.model')
		for i in range(len(model_list)):
			modelname=model_list[i]
			modelbasename=os.path.basename(modelname)
			timestamp='/'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[:3])+'/'+':'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[3:])
			delmod(vis=msname,scr=True)
			ft(vis=msname,model=modelname,usescratch=True)
			IB=ImageBasic(msname)
			uvrange_to_cal=IB.calc_calib_uvrange(4)[0]
			gaincal(vis=msname,caltable=caltable_name,timerange=timestamp,append=True,uvrange=uvrange_to_cal,solnorm=True,rmsthresh=[10,8,6],\
					refant=str(inputs.ref_ant),minsnr=inputs.gain_minsnr)
		os.system('cp -r '+caltable_name+' '+localdatabase)
		os.system('rm -rf '+caltable_name)
		if os.path.isdir(bandpass_modeldir)==False or len(glob.glob(bandpass_modeldir+'/*.model'))==0: # Backup bandpass
			mainlog.info('No bandpass models are available.\n')
		else:
			bp_caltable_name=str(OBSID)+'_'+str(coarse_chan_0)+'_'+str(coarse_chan_1)+'.bcal'
			model_list=glob.glob(bandpass_modeldir+'/*.model')
			modelname=model_list[0]
			modelbasename=os.path.basename(modelname)
			timestamp='/'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[:3])+'/'+':'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[3:])
			if len(timestamp)<10:
				applycal(vis=msname,gaintable=caltable_name,applymode='calflag',flagbackup=False)
				split(vis=msname,outputvis=msname.split('.ms')[0]+'_reftime.ms',datacolumn='corrected')
			else:
				index=timestamp.index(ref_time)
				timerange=','.join(timestamp[index-5:index+5])
				applycal(vis=msname,gaintable=caltable_name,timerange=timerange,applymode='calflag',flagbackup=False)
				split(vis=msname,outputvis=msname.split('.ms')[0]+'_reftime.ms',timerange=timerange,datacolumn='corrected')
			msname=msname.split('.ms')[0]+'_reftime.ms'
			freqlist=(AM.get_freqs()/10**6)
			for i in range(len(model_list)):
				modelname=model_list[i]
				f=imhead(imagename=modelname,mode='list')['crval4']/10**6
				df=(imhead(imagename=modelname,mode='list')['cdelt4']/10**6)/2.0
				spw='0:'+str(f-df)+'~'+str(f+df)+'MHz'
				ft(vis=msname,model=modelname,spw=spw,usescratch=True)
			bandpass(vis=msname,caltable=bp_caltable_name,solnorm=True,refant=str(ref_ant),minsnr=inputs.gain_minsnr)
			os.system('cp -r '+bp_caltable_name+' '+localdatabase)
			os.system('rm -rf '+bp_caltable_name)
		if os.path.isdir(polcal_caldir)==False or len(glob.glob(polcal_caldir+'/*.bin'))==0: # Backup polcal
			mainlog.info('No polarisation caltables are available.\n')
		else:
			polcaltable_list=glob.glob(polcal_caldir+'/*.bin')
			for polcal in polcaltable_list:
				freq=float(polcal.split('.bin')[0].split('freq_')[-1].split('_')[0]) # In MHz
				coarse_chan=int(freq/1.28)
				polcaltable_name=str(OBSID)+'_'+str(coarse_chan)+'.bin'
				os.system('mv '+polcal+' '+localdatabase+'/'+polcaltable_name)		
	return 0


# PAIRCARS master controller
############################

if basedir[-1]=='/':
	basedir=basedir[:-1]

if os.path.isdir(basedir+'/data')==False:
	os.makedirs(basedir+'/data')

os.system('mv selfcal_inputs_temp.py '+basedir+'/selfcal_inputs.py')
#os.system('cp -r *.py '+basedir) #TODO : remove after packaging
os.chdir(basedir)
sys.path.append(basedir)
import selfcal_inputs as inputs


# MPI check
###########
if os.path.isfile('.mpi_enabled'):
	os.system('rm -rf .mpi_enabled')
os.system('echo "import os\nos.system(\'touch .mpi_enabled\')" > test_mpi.py')
os.system('mpirun -np 1 -x OMP_NUM_THREADS=3 python3 test_mpi.py')
os.system('rm -rf test_mpi.py') 

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
	if inputs.safety_factor==0:
		calibrator_interval=90
	elif inputs.safety_factor==0:
		calibrator_interval=60
	elif inputs.safety_factor==0:
		calibrator_interval=30
elif inputs.quality_factor==1:
	if inputs.safety_factor==0:
		bandpass_interval=10
		calibrator_interval=60
	elif inputs.safety_factor==1:
		bandpass_interval=7
		calibrator_interval=30
	elif inputs.safety_factor==2:
		bandpass_interval=4
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
# Making local cal database folder
##################################
if os.path.isdir(inputs.local_caldatabase)==False:
	os.makedirs(inputs.local_caldatabase)
mainlog.info('Local calibration database is at : '+inputs.local_caldatabase+'\n')

# Download available cal files from paircars database #TODO
#####################################################





# Organising ms
###############
mainlog.info('Organising measurement sets .....\n')
measurement_set_list=glob.glob(inputs.basedir+'/data/*.ms')
msfreqs=[float(a.split('.ms')[0].split('_')[-1]) for a in measurement_set_list]
mstimes=[float(timestamp_to_mjdsec('/'.join(a.split('.ms')[0].split('_')[1:4])+'/'+':'.join(a.split('.ms')[0].split('_')[4:7]))) for a in measurement_set_list]
mstimes_iso=[('-'.join(a.split('.ms')[0].split('_')[1:4])+' '+':'.join(a.split('.ms')[0].split('_')[4:7])) for a in measurement_set_list]
ms_OBSIDs=[]
metafits_obsids_msdir=glob.glob(inputs.msdir+'/*.metafits')
if len(metafits_obsids_msdir)!=0:
	for metafits in metafits_obsids_msdir:
		os.system('cp -r '+metafits+' '+inputs.basedir+'/data/'+os.path.basename(metafits))
metafits_obsids=[int(os.path.basename(x).split('.metafits')[0]) for x in glob.glob(inputs.basedir+'/data/*.metafits')]
for i in range(len(measurement_set_list)):
	msname=measurement_set_list[i]
	obsid=get_OBSID(msname)
	if obsid==0 and len(metafits_obsids)!=0:
		mainlog.info('Could not connect to MWA metadata server. Searching for metafits in local data directory....\n')
		GPStime=int(Time(mstimes_iso[i],format='iso',scale='utc').gps)
		diff_gpstime=[]
		for a in metafits_obsids:
			if abs(a-GPStime)<480:
				diff_gpstime.append(a)
		if len(diff_gpstime)!=0:
			obsid=np.min(np.array(diff_gpstime))
			ms_OBSIDs.append(obsid)
		else:
			mainlog.info('Metafits file is not found for ms : '+msname+'. Removing ms from list.\n')
			measurement_set_list.remove(msname)
			msfreqs.remove(msfreqs[i])
			mstimes.remove(mstimes[i])
			mstimes_iso.remove(mstimes_iso[i])
	elif obsid==0 and len(metafits_obsids)==0:
		mainlog.info('Could not connect to MWA metadata server. No metafits files are found in local data directory. Exiting PAIRCARS.....\n')
		os._exit(0)
	else:
		ms_OBSIDs.append(obsid)

# Download metafits
###################
metafits_list=[]
mainlog.info('Downloading metafits files if does not exist.\n')
for i in range(len(measurement_set_list)):
	msname=measurement_set_list[i]
	metafits=str(ms_OBSIDs[i])+'.metafits'
	if os.path.isfile(inputs.basedir+'/data/'+metafits)==False: 
		metafits=download_metafits(msname,inputs.basedir+'/data')
		mainlog.info('Metafits file downloaded at : '+inputs.basedir+'/data/'+metafits+'\n')	
		metafits_list.append(inputs.basedir+'/data/'+metafits)
	else:
		mainlog.info('Metafits file found at : '+inputs.basedir+'/data/'+metafits+'\n')
		metafits_list.append(inputs.basedir+'/data/'+metafits)

ref_freq=str(msfreqs[int(len(msfreqs)/2)])
ref_time=mjdsec_to_timestamp(mstimes[int(len(mstimes)/2)],format=3)

# Phase center correction and decor correction
##############################################
for i in range(len(measurement_set_list)):
	msname=measurement_set_list[i]
	metafits=metafits_list[i]
	AM=AccessMS(msname)
	output=AM.move_phasecenter_to_sun()  # Moving the phasecenter to the Sun
	mainlog.info(output)
	if inputs.do_decor_correction: # Performing decorrelation correction and IAU convention change
		mainlog.info('Performing de-correlation correction and IAU convention correction for ms : '+msname+'\n')
		decor(msname,metafits,10,False)
	else: # If user do not want decorrelation correction perform only IAU convention
		mainlog.info('Correcting to IAU convention')
		AM.convert_mwa_to_iau()

# Selecting reference time frequnency ms
######################################

for i in range(len(measurement_set_list)):
	msname=measurement_set_list[i]
	if ref_freq in msname and ref_time in msname:
		ref_freq_time_msname=msname
		ref_time_freq_metafits=metafits_list[i]
mainlog.info('Reference time frequency measurement set : '+ref_freq_time_msname+'\n')
mainlog.info('######################\n')

############################################## # TODO : Do not apply, only make list
# Applying solutions from calibration database




# If no calibration is found in calibration database
# Apply calibrator solution
####################################################



calibrator_found=False
gaincal_list=[]
# Self calibration for reference time frequency ms
################################################## 
reftimefreq_ms_OBSID=get_OBSID(ref_freq_time_msname)
AM=AccessMS(ref_freq_time_msname)
workdir=inputs.basedir+'/'+os.path.basename(ref_freq_time_msname).split('.ms')[0]
if os.path.isdir(workdir)==False:
	os.makedirs(workdir)

spawned_ms_jobs={}
finished_ms=[]
if calibrator_found==True and len(gaincal_list)!=0: # If calibration found from paircars database or calibrator
	calibrator_OBSID=np.min(np.array([x.split('_')[0] for x in gaincal_list]))
	caltable_list=gaincal_list+bandpass_list
	if abs(reftimefreq_ms_OBSID-calibrator_OBSID)>calibrator_interval*60:
		ref_time_freq=True
	else:
		ref_time_freq=False
	os.system('cp -r selfcal_inputs.py '+workdir)
	ref_time,ref_chan,spawned_casa_instances=run_paircars_ms(ref_freq_time_msname,ref_time_freq_metafits,workdir,\
			ref_time_freq=ref_time_freq,do_bandpass=inputs.do_bandpass,do_polcal=inputs.do_polcal,calibrator_caltable=caltable_list)
	spawned_ms_jobs[ref_freq_time_msname]=[ref_time,ref_chan,spawned_casa_instances]
	while True:
		touch_file_list=glob.glob(inputs.basedir+'/.Finished_*cal_*'+os.path.basename(ref_freq_time_msname)+'*')
		if len(touch_file_list)==spawned_casa_instances:
			mainlog.info('Reference time frequency calibration is finished.\n')
			break
		else:
			time.sleep(2.0)
else:
	os.system('cp -r selfcal_inputs.py '+workdir)
	ref_time,ref_chan,spawned_casa_instances=run_paircars_ms(ref_freq_time_msname,ref_time_freq_metafits,workdir,\
			ref_time_freq=True,do_bandpass=inputs.do_bandpass,do_polcal=inputs.do_polcal,calibrator_caltable=[])
	spawned_ms_jobs[ref_freq_time_msname]=[ref_time,ref_chan,spawned_casa_instances]
	while True:
		touch_file_list=glob.glob(inputs.basedir+'/.Finished_*cal_*'+os.path.basename(ref_freq_time_msname)+'*')
		if len(touch_file_list)==spawned_casa_instances:
			mainlog.info('Reference time frequency calibration is finished.\n')
			break
		else:
			time.sleep(2.0)

index=measurement_set_list.index(ref_freq_time_msname)
measurement_set_list.remove(ref_freq_time_msname) # Removing ref time freq ms from ms list
ms_OBSIDs.remove(ms_OBSIDs[index]) # Removing ref time freq ms from ms list
metafits_list.remove(ref_time_freq_metafits) # Removing ref time freq metafit from list 
if len(measurement_set_list)>0:
	gtable=glob.glob(inputs.basedir+'/caltables/*'+str(os.path.basename(ref_freq_time_msname).split('.ms')[0])+'*.cal') # Ref time caltables
	bptable=glob.glob(inputs.basedir+'/bpcaltables/*'+str(os.path.basename(ref_freq_time_msname).split('.ms')[0])+'*.cal')
	caltable_list.append(gtable)
	caltable_list.append(bptable)
	mainlog.info('Caltables to be applied : '+','.join(caltable_list)+'\n')
	for i in range(len(measurement_set_list)):
		# Estimating total casa instances
		#################################
		total_available_cpu=psutil.cpu_count()-(psutil.cpu_count()*psutil.cpu_percent())
		available_cpu_for_paircars=int(total_available_cpu*inputs.cpu_frac)
		casa_instance=int(available_cpu_for_paircars/3)
		msname=measurement_set_list[i]
		obsid=ms_OBSIDs[i]
		time_diff=abs(reftimefreq_ms_OBSID-obsid)
		if time_diff>bandpass_interval*3600 and inputs.quality_factor!=0:
			mainlog.info('Performing bandpass as time interval is greater than : '+str(bandpass_interval)+' hr.\n')
			do_bandpass=True
		else:
			mainlog.info('Do not perform bandpass, becuase either time interval is smaller than : '+str(bandpass_interval)+' or quality_factor = 0.\n')
			do_bandpass=False
		metafits=metafits_list[i]
		if inputs.do_polcal==True:
			ref_time_freq_gridpoint=fits.getheader(ref_time_freq_metafits)['GRIDNUM']
			ms_gridpoint=fits.getheader(metafits)['GRIDNUM']
			if ref_time_freq_gridpoint!=ms_gridpoint:
				do_polcal=True
				mainlog.info('Beam pointing changed. Performing polarisation calibration.\n')
			else:
				mainlog.info('Beam pointing is same. Do not perform polarisation calibration.\n')
				do_polcal=False
		workdir=inputs.basedir+'/'+ps.path.basename(msname).split('.ms')[0]
		os.system('cp -r selfcal_inputs.py '+workdir)
	#	os.system('cp -r *.py '+workdir) #TODO : remove after packaging
		ref_time,ref_chan,spawned_casa_instances=\
			run_paircars_ms(msname,workdir,ref_time_freq=False,do_bandpass=do_bandpass,do_polcal=do_polcal,calibrator_caltable=caltable_list)	
		spawned_ms_jobs[msname]=[ref_time,ref_chan,spawned_casa_instances]
		available_casa_instance=casa_instance-spawned_casa_instance
		while True:
			if available_casa_instance>1:
				mainlog.info('At least 1 casa instance is available. Spawn new job...\n')
				break
			else:
				touch_file_list=glob.glob(inputs.basedir+'/.Finished_*cal_*'+os.path.basename(msname)+'*')
				available_casa_instance=len(touch_file_list)-spawned_casa_instances
				time.sleep(2.0)				

job_spawned_msname=spawned_ms_jobs.keys()
total_spawned_jobs=0
for ms in job_spawned_msname:
	total_spawned_job+=spawned_ms_jobs[ms][-1]

mainlog.info('Waiting for finishing all calibrations.\n')
while True:
	touch_files=len(glob.glob(inputs.basedir+'/Finished*cal*'))
	if touch_files==total_spawned_job:
		ms_list=glob.glob('time*')
		for ms in ms_list:
			msname=glob.glob(ms+'/*averaged.ms')[0]
			ref_time=spwaned_ms_jobs[ms][0]
			OBSID=get_OBSID(msname)
			gaincal_modeldir=inputs.basedir+'/imagemodels/'+str(OBSID)
			gaincal_modeldir=inputs.basedir+'/bpimagemodels/'+str(OBSID)
			polcal_caldir=inputs.basedir+'/polcaltables/'+str(OBSID)
			mainlog.info('Making final calibration tables for : '+msname+'\n')
			managing_caldatabase(msname,ref_time,gaincal_modedir,bandpass_modeldir,polcal_caldir,inputs.localdatabase)	
	else:
		time.sleep(2.0)
		

# TODO : Final leakage correction using background sources
# TODO : Calibration database
# TODO : Diagnostic plots
# TODO : Calibrator solution
# TODO : Start from the failed part or stopped part
# End of PAIRCARS
#################

















