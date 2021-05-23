'''
Code is written by Devojyoti Kansabanik , 2 May, 2021
'''
from casatools import *
from casatasks import *
import os,sys,logging,numpy as np,copy,glob,psutil,time,multiprocessing as mp,subprocess,getpass
from paircars.basic_func import *
from paircars.access_ms import *
from paircars.decor import *
from paircars.flagger import *
from paircars.fullpol_selfcal_LTS import *
from optparse import OptionParser
from astropy.io import fits
from CALIBRATE.access_calibrate import *
from paircars.libpaircars import send_paircars_notification,send_to_database

def fill_models(msname):
	'''
	Function to fill models from nearest freuencies in the ms
	Parameters:
	msname = Name of the meausreent set
	Return:
	Model filled ms
	'''
	AM=AccessMS(msname)
	model,nomodel=AM.get_model_nomodel_chan()
	model=np.array(model)
	nomodel=np.array(nomodel)
	tb=table()
	tb.open(msname,nomodify=False)
	modeldata=tb.getcol('MODEL_DATA')
	for i in nomodel:
		nearest_model_chan=model[np.argmin(abs(model-i))]
		modeldata[:,i,:]=modeldata[:,nearest_model_chan,:]
	tb.putcol('MODEL_DATA',modeldata)
	tb.flush()
	tb.close()
	print ('Fill models are done.\n')
	return

# MPI check
###########
def MPI_check():
	a=subprocess.getstatusoutput('mpirun -h')[0]
	if a==0:
		return 0
	else:
		return 1

def casa_instance_runner(cmd,screen_name,basedir,finished_touch_file,num_thread,num_casa_instance):
	'''
	Function to run a casa instance
	Parameters:
	cmd = Command to run
	screen_name = Name of the screen
	'''
	mpi_check=MPI_check()
	if basedir[-1]=='/':
		basedir=basedir[:-1]
	if mpi_check==0:
		cmd='mpirun -n 1 -x OMP_NUM_THREADS='+str(num_thread)+' -cpus-per-proc '+str(num_thread)+' '+cmd
	cmd+=';wait; if ! ls '+finished_touch_file+'_* ; then  touch '+finished_touch_file+'_error ;  fi; rm -rf '+basedir+'/'+screen_name+'.batch'
	os.system('echo "'+cmd+'" > '+basedir+'/'+screen_name+'.batch')
	screen_cmd='sh '+basedir+'/'+screen_name+'.batch'
	os.system('screen -S '+screen_name+' -X quit')	
	time.sleep(0.5)
	os.system('screen -mdS '+screen_name)
	time.sleep(0.5)
	obslog.info('########################\n')
	obslog.info('Made Screen : '+screen_name+'\n')
	obslog.info('Command : '+cmd+'\n')
	os.system('screen -S '+screen_name+' -X stuff \"'+screen_cmd+'\n"')	
	return screen_name


def make_final_gaincal(msname,workdir,caltable_name_prefix,freqavg=160,timeavg=2.0,ref_ant='1',gain_minsnr=3,modellist=[]):
	'''
	Function to make final gain caltable
	Parameters:
	msname = Name of the measurement set
	workdir = Name of the working directory
	caltable_name_prefix = Caltable name prefix with full path. For local caltable .gcal extension will be added, for global caltable .gcal.bin extension will be added
	freqavg = Frequency averaging in kHz
	timeavg = Time averaging in s 
	ref_ant = Reference antenna (default : 1)
	gain_minsnr = Minimum gain SNR for calibration (default : 3)
	modellist = [], list of gain calibration models
	Return:
	0,1 (For success or failure), List of local caltables, List of global caltables
	'''
	caltable_for_local_database=[]
	caltable_for_global_database=[]
	modellist=sorted(modellist)
	cal=CALIBRATE()
	original_msname=copy.deepcopy(msname)
	if msname[-1]=='/':
		msname=msname[:-1]
	if os.path.exists(workdir+'/'+os.path.basename(msname)+'.temp'):
		os.system('rm -rf '+workdir+'/'+os.path.basename(msname)+'.temp')
	os.system('cp -r '+msname+' '+workdir+'/'+os.path.basename(msname)+'.temp')
	msname=workdir+'/'+os.path.basename(msname)+'.temp'
	if len(modellist)==0:
		obslog.info('No models are present for gain calibration.\n')
		return 1,caltable_for_global_database,caltable_for_local_database
	else:
		local_caltable=caltable_name_prefix+'.gcal'
		global_caltable=caltable_name_prefix+'.gcal.bin'
		if os.path.isdir(local_caltable):
			os.system('rm -rf '+local_caltable)
		if os.path.isdir(global_caltable):
			os.system('rm -rf '+global_caltable)
		obslog.info('#######################################\n')
		obslog.info('Making final gaincal table for ms : '+original_msname+'\n')
		msname_copy=copy.deepcopy(msname)
		# Determining frequency and time stamps
		#######################################
		timerange_list=[]
		ref_timestamp=''
		ref_model=''
		for i in range(len(modellist)):
			modelname=modellist[i]
			modelbasename=os.path.basename(modelname)
			timestamp='/'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[:3])+'/'+':'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[3:])
			timerange_list.append(timestamp)
			if 'ref' in modelname:
				ref_model=modelname
				ref_timestamp+=timestamp+','	
		if ref_model=='':
			ref_index=int(len(modellist)/2)
			ref_model=modellist[ref_index]
			ref_timestamp=timerange_list[ref_index]
		else:
			ref_timestamp=ref_timestamp[:-1] # Reference time
		timerange_list=sorted(timerange_list)
		timerange=','.join(timerange_list)
		f=imhead(imagename=ref_model,mode='list')['crval4']/10**6
		df=(imhead(imagename=ref_model,mode='list')['cdelt4']/10**6)/2.0
		spw=str(f-df)+'~'+str(f+df)+'MHz' # Spectral window range for reference model
		AM=AccessMS(msname)
		freqres=AM.calc_freqres()
		chanwidth=int((float(freqavg))/freqres)
		timeres=AM.calc_timeres()

		ref_time_start=mjdsec_to_timestamp(timestamp_to_mjdsec(ref_timestamp,format=0)-float(timeavg/2),includedate=True,format=0)
		ref_time_end=mjdsec_to_timestamp(timestamp_to_mjdsec(ref_timestamp,format=0)+float(timeavg/2),includedate=True,format=0)
		ref_timestamp=ref_time_start+'~'+ref_time_end

		time_start=mjdsec_to_timestamp(timestamp_to_mjdsec(timerange_list[0],format=0)-float(timeavg/2),includedate=True,format=0)
		time_end=mjdsec_to_timestamp(timestamp_to_mjdsec(timerange_list[-1],format=0)+float(timeavg/2),includedate=True,format=0)
		timerange=time_start+'~'+time_end
		timebin=str(timeavg)+'s'
	
		# Spliting ms for local and global caltable, for global caltable only reference time solution will be there
		###########################################################################################################
		obslog.info('Spliting ms for gain calibration.....\n')
		if os.path.isdir(workdir+'/'+os.path.basename(msname)+'.gcalms'):
			os.system('rm -rf '+workdir+'/'+os.path.basename(msname)+'.gcalms '+workdir+'/'+os.path.basename(msname)+'.gcalms.flagversions')
		if os.path.isdir(workdir+'/'+os.path.basename(msname)+'.gcalms.ref'):
			os.system('rm -rf '+workdir+'/'+os.path.basename(msname)+'.gcalms.ref '+workdir+'/'+os.path.basename(msname)+'.gcalms.ref.flagversions')
		obslog.info('split(vis=\''+msname+'\',outputvis=\''+workdir+'/'+os.path.basename(msname)+'.gcalms\',datacolumn=\'DATA\',spw=\'0:'+spw+'\',timerange=\''+\
					timerange+'\',width='+str(chanwidth)+',timebin=\''+timebin+'\')\n')
		split(vis=msname,outputvis=workdir+'/'+os.path.basename(msname)+'.gcalms',datacolumn='DATA',spw='0:'+spw,timerange=timerange,width=chanwidth,timebin=timebin)
																																			# Reference channel ms
		local_cal_ms=workdir+'/'+os.path.basename(msname)+'.gcalms'
		obslog.info('split(vis=\''+msname+'\',outputvis=\''+workdir+'/'+os.path.basename(msname)+'.gcalms.ref\',datacolumn=\'DATA\',spw=\'0:'+spw+'\',timerange=\''+ref_timestamp\
					+'\',width='+str(chanwidth)+',timebin=\''+timebin+'\')\n')
		split(vis=msname,outputvis=workdir+'/'+os.path.basename(msname)+'.gcalms.ref',datacolumn='DATA',spw='0:'+spw,timerange=ref_timestamp,width=chanwidth,timebin=timebin)
																																				# Only reference channel and time ms
		global_cal_ms=workdir+'/'+os.path.basename(msname)+'.gcalms.ref'
		# Deleteing previous models and importing reference time model in reference time ms
		###################################################################################
		obslog.info('delmod(vis=\''+global_cal_ms+'\',scr=True)\n') # Deleting any previous models in ms and importing refernece time chan model
		delmod(vis=global_cal_ms,scr=True)
		AMgcal=AccessMS(global_cal_ms)
		freqs=AMgcal.get_freqs()
		lowfreq=(f-df)*10**6
		highfreq=(f+df)*10**6
		lowchan=np.argmin(abs(freqs-lowfreq))
		highchan=np.argmin(abs(freqs-highfreq))
		if highchan<lowchan:
			if highchan==0:
				highchan=len(freqs)-1
			else:
				highchan=lowchan
		spw=str(lowchan)+'~'+str(highchan) # Spectral window range for reference model
		obslog.info('ft(vis=\''+global_cal_ms+'\',model=\''+ref_model+'\',spw=\'0:'+str(spw)+'\',usescratch=True)\n')
		ft(vis=global_cal_ms,model=ref_model,spw='0:'+str(spw),usescratch=True)

		# Iteratively add models, perform gaincal, append them in a single caltable for all times 
		#########################################################################################
		IB=ImageBasic(local_cal_ms)
		uvrange_to_cal=IB.calc_calib_uvrange(12)[0]
		for i in range(len(modellist)): # Importing model for all times dequentially for reference chan and calibrating 
			modelname=modellist[i]
			modelbasename=os.path.basename(modelname)
			timestamp='/'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[:3])+'/'+':'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[3:])
			timerange_list.append(timestamp)
			obslog.info('delmod(vis=\''+local_cal_ms+'\',scr=True)\n')
			delmod(vis=local_cal_ms,scr=True)
			obslog.info('ft(vis=\''+local_cal_ms+'\',model=\''+modelname+'\',spw=\'0:'+str(spw)+'\',usescratch=True)\n')
			ft(vis=local_cal_ms,model=modelname,spw='0:'+str(spw),usescratch=True)
			obslog.info('gaincal(vis=\''+local_cal_ms+'\',caltable=\''+local_caltable+'\',spw=\'0:'+str(spw)+'\',timerange=\''+timestamp+\
						'\',append=True,uvrange=\''+uvrange_to_cal+'\',solnorm=True,rmsthresh=[10,8,6],refant=\''+str(ref_ant)+'\',minsnr='+str(gain_minsnr)+')\n')
			gaincal(vis=local_cal_ms,caltable=local_caltable,spw='0:'+str(spw),timerange=timestamp,append=True,\
					uvrange=uvrange_to_cal,solnorm=True,rmsthresh=[],refant=str(ref_ant),minsnr=gain_minsnr) # Appedning solutions into a same gain table

		# Making glocal caltable in CALIBRATE numpy format
		##################################################
		AMgcal=AccessMS(global_cal_ms)
		ntimes=AMgcal.get_num_timestamps()
		nchan=AMgcal.get_num_channels()
		IB1=ImageBasic(global_cal_ms)
		calib_uvrange_min=IB1.calc_calib_uvrange(12)[1]
		calib_uvrange_max=IB1.calc_calib_uvrange(12)[2]
		obslog.info('cal.calibrate(msname=\''+global_cal_ms+'\',caltable=\''+global_caltable+'\',calmode=\'diag\',minuv='+\
					str(calib_uvrange_min)+',maxuv='+str(calib_uvrange_max)+',quiet=True,a=\'0.001,0.0001\'t='+\
					str(int(ntimes/len(modellist)))+',j=3,ch='+str(int(AMgcal.get_num_channels()))+')\n') # Making gaincal table for reference time and chan
		cal.calibrate(msname=global_cal_ms,caltable=global_caltable,calmode='diag',minuv=calib_uvrange_min,maxuv=calib_uvrange_max,\
			quiet=True,a='0.001,0.0001',t=int(ntimes/len(modellist)),j=3,ch=int(AMgcal.get_num_channels()))
		os.system('rm -rf '+local_cal_ms+' '+local_cal_ms+'.flagversions')
		os.system('rm -rf '+global_cal_ms+' '+global_cal_ms+'.flagversions')
		if os.path.islink(msname_copy):
			os.system('unlink '+msname_copy)
		else:
			os.system('rm -rf '+msname_copy+'*')
		caltable_for_global_database.append(global_caltable)
		caltable_for_local_database.append(local_caltable)
		return 0,caltable_for_global_database,caltable_for_local_database


def make_final_bandpass(msname,workdir,caltable_name_prefix,gaintable=[],freqavg=160,timeavg=2.0,ref_ant='1',gain_minsnr=3,modellist=[]):
	'''
	Function to make final bandpass caltable
	Parameters:
	msname = Name of the measurement set
	workdir = Name of the working directory
	caltable_name_prefix = Caltable name prefix with full path. For local caltable .bcal extension will be added, for global caltable .bcal.bin extension will be added
	gaintable = [] , previous gaincal table list	
	ref_ant = Reference antenna (default : 1)
	gain_minsnr = Minimum gain SNR for calibration (default : 3)
	modellist = [], list of gain calibration models
	Return:
	0,1 (For success or failure), List of local caltables, List of global caltables
	'''
	caltable_for_local_database=[]
	caltable_for_global_database=[]
	modellist=sorted(modellist)
	original_msname=copy.deepcopy(msname)
	if msname[-1]=='/':
		msname=msname[:-1]
	if os.path.exists(workdir+'/'+os.path.basename(msname)+'.temp'):
		os.system('rm -rf '+workdir+'/'+os.path.basename(msname)+'.temp')
	os.system('cp -r '+msname+' '+workdir+'/'+os.path.basename(msname)+'.temp')
	msname=workdir+'/'+os.path.basename(msname)+'.temp'
	cal=CALIBRATE()
	if len(modellist)==0:
		obslog.info('No models are present for bandpass calibration.\n')
		return 1,caltable_for_global_database,caltable_for_local_database
	else:
		local_caltable=caltable_name_prefix+'.bcal'
		global_caltable=caltable_name_prefix+'.bcal.bin'
		if os.path.isdir(local_caltable):
			os.system('rm -rf '+local_caltable)
		if os.path.isdir(global_caltable):
			os.system('rm -rf '+global_caltable)
		obslog.info('#######################################\n')
		obslog.info('Making final bandpass table for ms : '+original_msname+'\n')
		msname_copy=copy.deepcopy(msname)

		# Determining frequency and time stamps
		#######################################
		ref_timestamp=''
		ref_models=[]
		spw=''
		for i in range(len(modellist)):
			modelname=modellist[i]
			modelbasename=os.path.basename(modelname)
			timestamp='/'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[:3])+'/'+':'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[3:])
			f=imhead(imagename=modelname,mode='list')['crval4']/10**6
			df=(imhead(imagename=modelname,mode='list')['cdelt4']/10**6)/2.0
			spw+=str(f-df)+'~'+str(f+df)+'MHz;' # Spectral window range for reference model
			ref_models.append(modelname)
			ref_timestamp=timestamp	
		spw=spw[:-1]
		if len(ref_models)==0:
			obslog.info('No models are present.\n')
			return 1,caltable_for_global_database,caltable_for_local_database
		AM=AccessMS(msname)
		freqres=AM.calc_freqres()
		chanwidth=int((float(freqavg))/freqres)
		timeres=AM.calc_timeres()
		timestamps=AM.get_timestamps()
		timeres=AM.calc_timeres()
		ref_time_mjd=timestamp_to_mjdsec(ref_timestamp,format=0)
		ref_time_start=mjdsec_to_timestamp(ref_time_mjd-float(timeavg/2),includedate=True,format=0)
		ref_time_end=mjdsec_to_timestamp(ref_time_mjd+float(timeavg/2),includedate=True,format=0)
		timerange=ref_time_start+'~'+ref_time_end
		timebin=str(timeavg)+'s'

		# Applying previous gaincal solutions
		#####################################
		if len(gaintable)!=0:
			obslog.info('applycal(vis=\''+msname+'\',gaintable='+str(gaintable)+',timerange=\''+timerange+\
								'\',applymode=\'calflag\',flagbackup=False)\n')			
			applycal(vis=msname,gaintable=gaintable,timerange=timerange,applymode='calflag',flagbackup=False)
			datacolumn_to_split='corrected'
		else:
			datacolumn_to_split='data'

		# Spliting ms for local and global caltable only at reference time 
		##################################################################
		obslog.info('Spliting ms for bandpass calibration.....\n')
		if os.path.isdir(workdir+'/'+os.path.basename(msname)+'.bcalms'):
			os.system('rm -rf '+workdir+'/'+os.path.basename(msname)+'.bcalms '+workdir+'/'+os.path.basename(msname)+'.bcalms.flagversions')
		obslog.info('split(vis=\''+msname+'\',outputvis=\''+workdir+'/'+os.path.basename(msname)+'.bcalms\',datacolumn=\''+str(datacolumn_to_split)+'\',spw=\'0:'+spw+'\',timerange=\''+\
					timerange+'\',width='+str(chanwidth)+',timebin=\''+timebin+'\')\n')
		split(vis=msname,outputvis=workdir+'/'+os.path.basename(msname)+'.bcalms',datacolumn=datacolumn_to_split,spw='0:'+spw,timerange=timerange,width=chanwidth,timebin=timebin)
																																			# Bandpass calibration ms
		local_cal_ms=workdir+'/'+os.path.basename(msname)+'.bcalms'
		global_cal_ms=workdir+'/'+os.path.basename(msname)+'.bcalms'

		# Deleteing previous models and importing reference time model
		##############################################################
		obslog.info('delmod(vis=\''+local_cal_ms+'\',scr=True)\n') # Deleting any previous models in ms and importing models
		delmod(vis=local_cal_ms,scr=True)
		AMbcal=AccessMS(local_cal_ms)
		freqs=AMbcal.get_freqs()
		for i in range(len(modellist)):
			modelname=modellist[i]
			modelbasename=os.path.basename(modelname)
			timestamp='/'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[:3])+'/'+':'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[3:])
			f=imhead(imagename=modelname,mode='list')['crval4']/10**6
			df=(imhead(imagename=modelname,mode='list')['cdelt4']/10**6)/2.0
			lowfreq=(f-df)*10**6
			highfreq=(f+df)*10**6
			lowchan=np.argmin(abs(freqs-lowfreq))
			highchan=np.argmin(abs(freqs-highfreq))
			if highchan<lowchan:
				if highchan==0:
					highchan=len(freqs)-1
				else:
					highchan=lowchan
			spw=str(lowchan)+'~'+str(highchan) # Spectral window range for model
			obslog.info('ft(vis=\''+local_cal_ms+'\',model=\''+modelname+'\',spw=\'0:'+str(spw)+'\',usescratch=True)\n')
			ft(vis=local_cal_ms,model=modelname,spw='0:'+str(spw),usescratch=True)

		# Perform bandpass in a single caltable
		#######################################
		IB=ImageBasic(local_cal_ms)
		uvrange_to_cal=IB.calc_calib_uvrange(12)[0]
		obslog.info('bandpass(vis=\''+local_cal_ms+'\',caltable='+str(local_caltable)+',solnorm=True,refant=\''\
					+str(ref_ant)+'\',minsnr='+str(gain_minsnr)+',uvrange=\''+uvrange_to_cal+'\')')
		bandpass(vis=local_cal_ms,caltable=local_caltable,solnorm=True,refant=str(ref_ant),minsnr=gain_minsnr,uvrange=uvrange_to_cal) # Performing bandpass
			
		# Making glocal caltable in CALIBRATE numpy format
		##################################################
		AMgcal=AccessMS(global_cal_ms)
		ntimes=AMgcal.get_num_timestamps()
		nchan=AMgcal.get_num_channels()
		IB1=ImageBasic(global_cal_ms)
		calib_uvrange_min=IB1.calc_calib_uvrange(12)[1]
		calib_uvrange_max=IB1.calc_calib_uvrange(12)[2]
		model_chan,nomodel_chan=AMgcal.get_model_nomodel_chan()
		if len(nomodel_chan)!=0:
			flagchans=[str(i) for i in nomodel_chan]
			if len(flagchans)!=0:
				flag_chan_str=';'.join(flagchans)
				obslog.info('flagdata(vis=\''+global_cal_ms+'\',spw=\'0:'+flag_chan_str+'\',flagbackup=False)\n')
				flagdata(vis=global_cal_ms,spw='0:'+flag_chan_str,flagbackup=False)
		obslog.info('cal.calibrate(msname=\''+global_cal_ms+'\',caltable=\''+global_caltable+'\',calmode=\'diag\',minuv='+\
					str(calib_uvrange_min)+',maxuv='+str(calib_uvrange_max)+',quiet=True,a=\'0.001,0.0001\'t='+\
					str(int(ntimes))+',j=3,ch=1)\n') # Making bandpass table for reference time
		cal.calibrate(msname=global_cal_ms,caltable=global_caltable,calmode='diag',minuv=calib_uvrange_min,maxuv=calib_uvrange_max,\
			quiet=True,a='0.001,0.0001',t=int(ntimes),j=3,ch=1)
		os.system('rm -rf '+local_cal_ms+' '+local_cal_ms+'.flagversions')
		os.system('rm -rf '+global_cal_ms+' '+global_cal_ms+'.flagversions')
		if os.path.islink(msname_copy):
			os.system('unlink '+msname_copy)
		else:
			os.system('rm -rf '+msname_copy+'*')
		caltable_for_global_database.append(global_caltable)
		caltable_for_local_database.append(local_caltable)
		return 0,caltable_for_global_database,caltable_for_local_database

def make_final_leakcal(msname,workdir,caltable_name_prefix,gaintable=[],freqavg=160,timeavg=2.0,ref_ant='1',gain_minsnr=3,modellist=[]):
	'''
	Function to make final leakage corrected differential gain caltable
	Parameters:
	msname = Name of the measurement set
	workdir = Name of the working directory
	caltable_name_prefix = Caltable name prefix with full path. For local caltable .gcal extension will be added, for global caltable .gcal.bin extension will be added
	gaintable = [] , previous gaincal table list	
	ref_ant = Reference antenna (default : 1)
	gain_minsnr = Minimum gain SNR for calibration (default : 3)
	modellist = [], list of gain calibration models
	Return:
	0,1 (For success or failure), List of local caltables, List of global caltables
	'''
	caltable_for_local_database=[]
	caltable_for_global_database=[]
	modellist=sorted(modellist)
	original_msname=copy.deepcopy(msname)
	if msname[-1]=='/':
		msname=msname[:-1]
	if os.path.exists(workdir+'/'+os.path.basename(msname)+'.temp'):
		os.system('rm -rf '+workdir+'/'+os.path.basename(msname)+'.temp')
	os.system('cp -r '+msname+' '+workdir+'/'+os.path.basename(msname)+'.temp')
	msname=workdir+'/'+os.path.basename(msname)+'.temp'
	cal=CALIBRATE()
	if len(modellist)==0:
		obslog.info('No models are present for bandpass calibration.\n')
		return 1,caltable_for_global_database,caltable_for_local_database
	else:
		local_caltable=caltable_name_prefix+'.lcal'
		global_caltable=caltable_name_prefix+'.lcal.bin'
		if os.path.isdir(local_caltable):
			os.system('rm -rf '+local_caltable)
		if os.path.isdir(global_caltable):
			os.system('rm -rf '+global_caltable)
		obslog.info('#######################################\n')
		obslog.info('Making final leakage corrected table for ms : '+original_msname+'\n')
		msname_copy=copy.deepcopy(msname)

		# Determining frequency and time stamps
		#######################################
		ref_timestamp=''
		ref_models=[]
		spw=''
		for i in range(len(modellist)):
			modelname=modellist[i]
			modelbasename=os.path.basename(modelname)
			timestamp='/'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[:3])+'/'+':'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[3:])
			f=imhead(imagename=modelname,mode='list')['crval4']/10**6
			df=(imhead(imagename=modelname,mode='list')['cdelt4']/10**6)/2.0
			spw+=str(f-df)+'~'+str(f+df)+'MHz;' # Spectral window range for model
			ref_models.append(modelname)
			ref_timestamp=timestamp	
		spw=spw[:-1]
		if len(ref_models)==0:
			obslog.info('No models are present.\n')
			return 1,caltable_for_global_database,caltable_for_local_database
		AM=AccessMS(msname)
		freqres=AM.calc_freqres()
		chanwidth=int((float(freqavg))/freqres)
		timeres=AM.calc_timeres()
		timestamps=AM.get_timestamps()
		timeres=AM.calc_timeres()
		ref_time_mjd=timestamp_to_mjdsec(ref_timestamp,format=0)
		ref_time_start=mjdsec_to_timestamp(ref_time_mjd-float(timeavg/2),includedate=True,format=0)
		ref_time_end=mjdsec_to_timestamp(ref_time_mjd+float(timeavg/2),includedate=True,format=0)
		timerange=ref_time_start+'~'+ref_time_end
		timebin=str(timeavg)+'s'

		# Applying previous gaincal solutions
		#####################################
		if len(gaintable)!=0:
			obslog.info('applycal(vis=\''+msname+'\',gaintable='+str(gaintable)+',timerange=\''+timerange+\
								'\',applymode=\'calflag\',flagbackup=False)\n')			
			applycal(vis=msname,gaintable=gaintable,timerange=timerange,applymode='calflag',flagbackup=False)
			datacolumn_to_split='corrected'
		else:
			datacolumn_to_split='data'
		
		# Spliting ms for local and global caltable for reference time only
		###################################################################
		obslog.info('Spliting ms for leakage corrected calibration.....\n')
		if os.path.isdir(workdir+'/'+os.path.basename(msname)+'.lcalms'):
			os.system('rm -rf '+workdir+'/'+os.path.basename(msname)+'.lcalms '+workdir+'/'+os.path.basename(msname)+'.lcalms.flagversions')
		obslog.info('split(vis=\''+msname+'\',outputvis=\''+workdir+'/'+os.path.basename(msname)+'.lcalms\',datacolumn=\''+str(datacolumn_to_split)+'\',spw=\'0:'+spw+'\',timerange=\''+\
					timerange+'\',width='+str(chanwidth)+',timebin=\''+timebin+'\')\n')
		split(vis=msname,outputvis=workdir+'/'+os.path.basename(msname)+'.lcalms',datacolumn=datacolumn_to_split,spw='0:'+spw,timerange=timerange,width=chanwidth,timebin=timebin)
																																				 # Leakge correction ms
		local_cal_ms=workdir+'/'+os.path.basename(msname)+'.lcalms'
		global_cal_ms=workdir+'/'+os.path.basename(msname)+'.lcalms'
	
		# Deleteing previous models and importing reference time model 
		##############################################################
		obslog.info('delmod(vis=\''+local_cal_ms+'\',scr=True)\n') # Deleting any previous models in ms and importing refernece time model
		delmod(vis=local_cal_ms,scr=True)
		AMlcal=AccessMS(local_cal_ms)
		freqs=AMlcal.get_freqs()
		for i in range(len(modellist)):
			modelname=modellist[i]
			modelbasename=os.path.basename(modelname)
			timestamp='/'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[:3])+'/'+':'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[3:])
			f=imhead(imagename=modelname,mode='list')['crval4']/10**6
			df=(imhead(imagename=modelname,mode='list')['cdelt4']/10**6)/2.0
			lowfreq=(f-df)*10**6
			highfreq=(f+df)*10**6
			lowchan=np.argmin(abs(freqs-lowfreq))
			highchan=np.argmin(abs(freqs-highfreq))
			if highchan<lowchan:
				if highchan==0:
					highchan=len(freqs)-1
				else:
					highchan=lowchan
			spw=str(lowchan)+'~'+str(highchan) # Spectral window range for model
			obslog.info('ft(vis=\''+local_cal_ms+'\',model=\''+modelname+'\',spw=\'0:'+str(spw)+'\',usescratch=True)\n')
			ft(vis=local_cal_ms,model=modelname,spw='0:'+str(spw),usescratch=True)

		# Perform leakage corrected differential bandpass in a single caltable
		######################################################################
		IB=ImageBasic(local_cal_ms)
		uvrange_to_cal=IB.calc_calib_uvrange(12)[0]
		obslog.info('bandpass(vis=\''+local_cal_ms+'\',caltable='+str(local_caltable)+',solnorm=True,refant=\''\
					+str(ref_ant)+'\',minsnr='+str(gain_minsnr)+',uvrange=\''+uvrange_to_cal+'\')')
		bandpass(vis=local_cal_ms,caltable=local_caltable,solnorm=True,refant=str(ref_ant),minsnr=gain_minsnr,uvrange=uvrange_to_cal) # Performing leakage corrected differential bandpass
			
		# Making glocal caltable in CALIBRATE numpy format
		##################################################
		AMgcal=AccessMS(global_cal_ms)
		ntimes=AMgcal.get_num_timestamps()
		nchan=AMgcal.get_num_channels()
		IB1=ImageBasic(global_cal_ms)
		model_chan,nomodel_chan=AMgcal.get_model_nomodel_chan()
		if len(nomodel_chan)!=0:
			flagchans=[str(i) for i in nomodel_chan]
			if len(flagchans)!=0:
				flag_chan_str=';'.join(flagchans)
				obslog.info('flagdata(vis=\''+global_cal_ms+'\',spw=\'0:'+flag_chan_str+'\',flagbackup=False)\n')
				flagdata(vis=global_cal_ms,spw='0:'+flag_chan_str,flagbackup=False)
		calib_uvrange_min=IB1.calc_calib_uvrange(12)[1]
		calib_uvrange_max=IB1.calc_calib_uvrange(12)[2]
		obslog.info('cal.calibrate(msname=\''+global_cal_ms+'\',caltable=\''+global_caltable+'\',calmode=\'diag\',minuv='+\
					str(calib_uvrange_min)+',maxuv='+str(calib_uvrange_max)+',quiet=True,a=\'0.001,0.0001\'t='+\
					str(int(ntimes))+',j=3,ch=1)\n') # Making leakage corrected bandpass table for reference time and chan
		cal.calibrate(msname=global_cal_ms,caltable=global_caltable,calmode='diag',minuv=calib_uvrange_min,maxuv=calib_uvrange_max,\
			quiet=True,a='0.001,0.0001',t=int(ntimes),j=3,ch=1)
		os.system('rm -rf '+local_cal_ms+' '+local_cal_ms+'.flagversions')
		os.system('rm -rf '+global_cal_ms+' '+global_cal_ms+'.flagversions')
		if os.path.islink(msname_copy):
			os.system('unlink '+msname_copy)
		else:
			os.system('rm -rf '+msname_copy+'*')
		caltable_for_global_database.append(global_caltable)
		caltable_for_local_database.append(local_caltable)
		return 0,caltable_for_global_database,caltable_for_local_database

def make_final_polcal(msname,metafits,workdir,caltable_name_prefix,gaintable=[],freqavg=160,timeavg=2.0,pol_skip_freq=1280,modellist=[]):
	'''
	Function to make final gain caltable
	Parameters:
	msname = Name of the measurement set
	metafits = Name of the metafits file
	workdir = Name of the working directory
	caltable_name_prefix = Caltable name prefix with full path. For local caltable .gcal extension will be added, for global caltable .gcal.bin extension will be added
	gaintable = [] , previous gaincal table list	
	pol_skip_freq = Frequency to skip polarisation calibration in kHz
	modellist = [], list of gain calibration models
	Return:
	0,1 (For success or failure), List of local caltables, List of global caltables
	'''
	caltable_for_local_database=[]
	caltable_for_global_database=[]
	modellist=sorted(modellist)
	original_msname=copy.deepcopy(msname)
	if msname[-1]=='/':
		msname=msname[:-1]
	if os.path.exists(workdir+'/'+os.path.basename(msname)+'.temp'):
		os.system('rm -rf '+workdir+'/'+os.path.basename(msname)+'.temp')
	os.system('cp -r '+msname+' '+workdir+'/'+os.path.basename(msname)+'.temp')
	msname=workdir+'/'+os.path.basename(msname)+'.temp'		
	cal=CALIBRATE()
	if len(modellist)==0:
		obslog.info('No models are present for polarisation calibration.\n')
		return 1,caltable_for_global_database,caltable_for_local_database
	else:
		beam_caltable=caltable_name_prefix+'.beam.bin'
		polcal_caltable=caltable_name_prefix+'.pcal.bin'
		if os.path.isdir(beam_caltable):
			os.system('rm -rf '+beam_caltable)
		if os.path.isdir(polcal_caltable):
			os.system('rm -rf '+polcal_caltable)
		obslog.info('#######################################\n')
		obslog.info('Making final polarisation caltable for ms : '+original_msname+'\n')
		msname_copy=copy.deepcopy(msname)

		# Determining frequency and time stamps
		#######################################
		ref_timestamp=''
		ref_models=[]
		spw=''
		for i in range(len(modellist)):
			modelname=modellist[i]
			modelbasename=os.path.basename(modelname)
			timestamp='/'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[:3])+'/'+':'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[3:])
			f=imhead(imagename=modelname,mode='list')['crval4']/10**6
			df=(imhead(imagename=modelname,mode='list')['cdelt4']/10**6)/2.0
			spw+=str(f-df)+'~'+str(f+df)+'MHz;' # Spectral window range for reference model
			ref_models.append(modelname)
			ref_timestamp=timestamp	
		spw=spw[:-1]
		if len(ref_models)==0:
			obslog.info('No models are present.\n')
			return 1,caltable_for_global_database,caltable_for_local_database
		AM=AccessMS(msname)
		freqres=AM.calc_freqres()
		chanwidth=int((float(freqavg))/freqres)
		timeres=AM.calc_timeres()
		timestamps=AM.get_timestamps()
		timeres=AM.calc_timeres()
		ref_time_mjd=timestamp_to_mjdsec(ref_timestamp,format=0)
		ref_time_start=mjdsec_to_timestamp(ref_time_mjd-float(timeavg/2),includedate=True,format=0)
		ref_time_end=mjdsec_to_timestamp(ref_time_mjd+float(timeavg/2),includedate=True,format=0)
		timerange=ref_time_start+'~'+ref_time_end
		timebin=str(timeavg)+'s'

		# Applying previous gaincal solutions
		#####################################
		if len(gaintable)!=0:
			obslog.info('applycal(vis=\''+msname+'\',gaintable='+str(gaintable)+',timerange=\''+timerange+\
								'\',applymode=\'calflag\',flagbackup=False)\n')			
			applycal(vis=msname,gaintable=gaintable,timerange=timerange,applymode='calflag',flagbackup=False)
			datacolumn_to_split='corrected'
		else:
			datacolumn_to_split='data'

		# Spliting ms for local and global caltable, for global caltable only reference time solution will be there
		###########################################################################################################
		obslog.info('Spliting ms for polarisation calibration.....\n')
		obslog.info('###########################################\n')
		if os.path.isdir(workdir+'/'+os.path.basename(msname)+'.beamcalms'):
			os.system('rm -rf '+workdir+'/'+os.path.basename(msname)+'.beamcalms '+workdir+'/'+os.path.basename(msname)+'.beamcalms.flagversions')
		obslog.info('split(vis=\''+msname+'\',outputvis=\''+workdir+'/'+os.path.basename(msname)+'.beamcalms\',datacolumn=\''+\
				str(datacolumn_to_split)+'\',spw=\'0:'+spw+'\',timerange=\''+timerange+'\',width='+str(chanwidth)+',timebin=\''+timebin+'\')\n')
		split(vis=msname,outputvis=workdir+'/'+os.path.basename(msname)+'.beamcalms',datacolumn=datacolumn_to_split,spw='0:'+spw,timerange=timerange,width=chanwidth,timebin=timebin)
																																				 # Reference time ms
		beamcal_ms=workdir+'/'+os.path.basename(msname)+'.beamcalms'
		polcal_ms=workdir+'/'+os.path.basename(msname)+'.pcalms'
	
		# Perform ideal beam correction
		###############################
		obslog.info('Performing ideal beam correction......\n')
		obslog.info('###########################################\n')
		PSC=PolSelfcal(beamcal_ms,metafits,32*60,verbose=False,interactive=False) # Performing ideal beam correction at phase center for every coarse channels
		obslog.info('PSC.correct_visibility_single_beam_jones(modify_datacolumn=False,skip_freq='+str(pol_skip_freq)+',save_beamfile=\''+str(beam_caltable)+'\')\n')
		PSC.correct_visibility_single_beam_jones(modify_datacolumn=False,skip_freq=float(pol_skip_freq),save_beamfile=beam_caltable)
		caltable_for_global_database.append(beam_caltable)
		caltable_for_local_database.append(beam_caltable)

		# Spliting beam corrected ms
		############################
		obslog.info('Spliting beam corrected measurement set ........\n')
		if os.path.exists(polcal_ms)==True:
			os.system('rm -rf '+polcal_ms+' '+polcal_ms+'.flagversions')
		obslog.info('split(vis=\''+beamcal_ms+'\',outputvis=\''+polcal_ms+'\',datacolumn=\'corrected\')\n')		
		split(vis=beamcal_ms,outputvis=polcal_ms,datacolumn='corrected')
		# Deleteing previous models and importing reference time model in reference time ms
		###################################################################################
		obslog.info('delmod(vis=\''+polcal_ms+'\',scr=True)\n') # Deleting any previous models in ms and importing refernece time chan model
		delmod(vis=polcal_ms,scr=True)
		AMpcal=AccessMS(polcal_ms)
		freqs=AMpcal.get_freqs()

		for i in modellist:
			if 'leakage' in i:
				modellist.remove(i)

		for i in range(len(modellist)):
			modelname=modellist[i]
			modelbasename=os.path.basename(modelname)
			timestamp='/'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[:3])+'/'+':'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[3:])
			f=imhead(imagename=modelname,mode='list')['crval4']/10**6
			df=(imhead(imagename=modelname,mode='list')['cdelt4']/10**6)/2.0
			lowfreq=(f-df)*10**6
			highfreq=(f+df)*10**6
			lowchan=np.argmin(abs(freqs-lowfreq))
			highchan=np.argmin(abs(freqs-highfreq))
			if highchan<lowchan:
				if highchan==0:
					highchan=len(freqs)-1
				else:
					highchan=lowchan
			spw=str(lowchan)+'~'+str(highchan) # Spectral window range for reference model
			obslog.info('ft(vis=\''+polcal_ms+'\',model=\''+modelname+'\',spw=\'0:'+str(spw)+'\',usescratch=True)\n')
			ft(vis=polcal_ms,model=modelname,spw='0:'+str(spw),usescratch=True)
		
		# Making polcal caltable in CALIBRATE numpy format
		##################################################
		AMgcal=AccessMS(polcal_ms)
		ntimes=AMgcal.get_num_timestamps()
		nchan=AMgcal.get_num_channels()
		IB1=ImageBasic(polcal_ms)
		model_chan,nomodel_chan=AMgcal.get_model_nomodel_chan()
		if len(nomodel_chan)!=0:
			flagchans=[str(i) for i in nomodel_chan]
			if len(flagchans)!=0:
				flag_chan_str=';'.join(flagchans)
				obslog.info('flagdata(vis=\''+polcal_ms+'\',spw=\'0:'+flag_chan_str+'\',flagbackup=False)\n')
				flagdata(vis=polcal_ms,spw='0:'+flag_chan_str,flagbackup=False)
		calib_uvrange_min=IB1.calc_calib_uvrange(12)[1]
		calib_uvrange_max=IB1.calc_calib_uvrange(12)[2]
		obslog.info('cal.calibrate(msname=\''+polcal_ms+'\',caltable=\''+polcal_caltable+'\',minuv='+\
					str(calib_uvrange_min)+',maxuv='+str(calib_uvrange_max)+',quiet=True,a=\'0.001,0.0001\'t='+\
					str(int(ntimes))+',j=3,ch=1)\n') # Making gaincal table for reference time and chan
		cal.calibrate(msname=polcal_ms,caltable=polcal_caltable,minuv=calib_uvrange_min,maxuv=calib_uvrange_max,\
			quiet=True,t=int(ntimes),j=3,ch=1)
		obslog.info('cal.applycal(msname=\''+polcal_ms+'\',gaintable=\''+polcal_caltable+'\',applymode=\'calflag\')\n')
		cal.applycal(msname=polcal_ms,gaintable=polcal_caltable,applymode='calflag')
		os.system('rm -rf '+beamcal_ms+' '+beamcal_ms+'.flagversions')
		os.system('rm -rf '+polcal_ms+' '+polcal_ms+'.flagversions')
		if os.path.islink(msname_copy):
			os.system('unlink '+msname_copy)
		else:
			os.system('rm -rf '+msname_copy+'*')
		caltable_for_global_database.append(polcal_caltable)
		caltable_for_local_database.append(polcal_caltable)
		return 0,caltable_for_global_database,caltable_for_local_database

def managing_caldatabase(msname,metafits,total_spawned_jobs,basedir,gaincal_modeldir,bandpass_modeldir,polcal_modeldir,localdatabase,gain_minsnr,\
						ref_ant,freq_avg=160.0,time_avg=2.0,pol_skip_freq=1.28):
	'''
	Function to manager calibration database
	Parameters:
	msname = Name of the measurement set for a set of contiguous coarse channels with all timestamps
	metafits = Name of the metafits file
	total_spawned_jobs = Total calibration jobs spawned for the ms
	gaincal_modeldir = Model directory for gaincal
	bandpass_modeldir = Model directory for bandpass
	polcal_modeldir = Model directory for polarisation
	localdatabase = Local database directory
	gain_minsnr = Minimum SNR for gain calibration
	ref_ant = Reference antenna
	freq_avg = Frequency averaging done during calibration
	time_avg = Time averaging done during calibration
	pol_skip_freq = Frequency interval to perform polarisation calibration
	Return:
	Message code, global database caltable list, local database caltable list
	'''
	# Setting up logs and waiting for calibrations jobs to finish
	#############################################################
	cal=CALIBRATE()
	OBSID=int(fits.getheader(metafits)['GPSTIME'])
	obslog.info('Waiting for finishing '+str(total_spawned_jobs)+' calibration jobs for ms : '+msname+'\n')
	caltable_for_global_database=[]
	caltable_for_local_database=[]
	# Setting up different directory names
	######################################
	if basedir[-1]=='/':
		basedir=basedir[:-1]
	if msname[-1]=='/':
		msname=msname[:-1]
	if gaincal_modeldir[-1]=='/':
		gaincal_modeldir=gaincal_modeldir[:-1]
	if bandpass_modeldir[-1]=='/':
		bandpass_modeldir=bandpass_modeldir[:-1]
	if polcal_modeldir[-1]=='/':
		polcal_modeldir=polcal_modeldir[:-1]
	if localdatabase[-1]=='/':
		localdatabase=localdatabase[:-1]

	basemsdir=os.path.basename(msname).split('.ms')[0] # Base directory for the ms inside model directories
	
	while True:  # Waiting for all jobs to finish
		touch_files=len(glob.glob(basedir+'/.Finished*cal*'+str(OBSID)+'*'+basemsdir+'*'))
		if touch_files>=total_spawned_jobs:
			obslog.info('All calibration jobs for ms : '+msname+' is completed.\n')
			break	
		else:
			time.sleep(2.0)

	batch_files=glob.glob(basedir+'/'+basemsdir+'/*.batch')
	mpiapp_files=glob.glob(basedir+'/'+basemsdir+'/*cal_mpicmd*')
	if len(batch_files)!=0:
		for i in batch_files:
			os.system('rm -rf '+i)
	if len(mpiapp_files)!=0:
		for j in mpiapp_files:
			os.system('rm -rf '+j)	

	if inputs.clear_screen==True:
		screen_list=[os.path.basename(i) for i in glob.glob('/var/run/screen/S-'+str(getpass.getuser())+'/*')]
		for i in screen_list:
			if str(OBSID)+'_' not in i:
				screen_list.remove(i) 

		for i in screen_list:
			os.system('screen -S '+i+' -X quit')

	if localdatabase=='' or os.path.isdir(localdatabase)==False: # If localdabase is available, making the localdatabase directory
		obslog.error('Local data base not found. Making local database at basedir.\n')
		localdatabase=basedir+'/localdatabase/'+str(OBSID)
		if os.path.isdir(localdatabase)==False:
			os.makedirs(localdatabase)
	else:
		obslog.info('Local data base is at : '+localdatabase+'\n')
		localdatabase=localdatabase+'/'+str(OBSID)
		if os.path.isdir(localdatabase)==False:
			os.makedirs(localdatabase)
	os.system('rm -rf '+basedir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftime*')	
	
	if len(glob.glob(gaincal_modeldir+'/*.model'))==0:
		gaincal_modeldir=gaincal_modeldir+'/'+basemsdir
	if len(glob.glob(bandpass_modeldir+'/*.model'))==0:
		bandpass_modeldir=bandpass_modeldir+'/'+basemsdir
	if len(glob.glob(polcal_modeldir+'/*.model'))==0:
		polcal_modeldir=polcal_modeldir+'/'+basemsdir
		
	os.system('rm -rf '+basedir+'/*.ms.temp*')
	# Final gain calibrations
	#########################
	if gaincal_modeldir=='' or os.path.isdir(gaincal_modeldir)==False or len(glob.glob(gaincal_modeldir+'/*.model'))==0:
		obslog.info('No gaincal models are available.\n')  # Checking if any gaincal model is present or not. If not exit at this stage
		gcal_msg=2
		return 1,caltable_for_global_database,caltable_for_local_database
	else:
		AM=AccessMS(msname)
		antenna=AM.get_antenna_string()
		unflagchan,flagchan=flag_MWA_coarse(msname,edgewidth=160,do_flag=False,force=False)
	#	obslog.info('flagdata(vis=\''+msname+'\',mode=\'unflag\',antenna=\''+antenna+'\',spw=\''+unflagchan+'\')\n')
	#	flagdata(vis=msname,mode='unflag',antenna=antenna,spw=unflagchan)
		freqs=AM.get_freqs()/10**6
		nchan_avg=int(freq_avg/AM.calc_freqres())
		coarse_chan_0=freq_to_MWA_coarse(freqs[0])
		coarse_chan_1=freq_to_MWA_coarse(freqs[-1])
		caltable_name_prefix=basedir+'/'+str(OBSID)+'_'+str(coarse_chan_0)+'_'+str(coarse_chan_1) # Caltable name prefix; OBSID_startcoarsechan_endcoarsechan format
		obslog.info('Searching gaincal models in model directory : '+gaincal_modeldir+'\n')
		model_list=glob.glob(gaincal_modeldir+'/*.model') # Gaincal model list
		gcal_msg,global_database_list,local_database_list=make_final_gaincal(msname,basedir,caltable_name_prefix,freqavg=freq_avg,timeavg=time_avg,ref_ant=ref_ant,\
																				gain_minsnr=gain_minsnr,modellist=model_list)
		caltable_for_local_database+=local_database_list
		caltable_for_global_database+=global_database_list	
	
		if gcal_msg==0:
			obslog.info('Searching bandpass models in model directory : '+bandpass_modeldir+'\n')
			bpmodel_list=glob.glob(bandpass_modeldir+'/*.model')
			if len(bpmodel_list)==0:
				obslog.info('No bandpass models are available.\n')  # Checking if any gaincal model is present or not. If not exit at this stage
				bp_msg=2
			else:
				bp_msg,global_database_list,local_database_list=make_final_bandpass(msname,basedir,caltable_name_prefix,gaintable=caltable_for_local_database,\
								freqavg=freq_avg,timeavg=time_avg,ref_ant=ref_ant,gain_minsnr=gain_minsnr,modellist=bpmodel_list)
				caltable_for_local_database+=local_database_list
				caltable_for_global_database+=global_database_list	

		if gcal_msg==0:
			obslog.info('Searching leakage corrected models in model directory : '+polcal_modeldir+'\n')
			polcalmodel_list=glob.glob(polcal_modeldir+'/*.model')
			lcal_model_list=[]
			for i in polcalmodel_list:
				if 'leakage' in i:
					lcal_model_list.append(i)
			if len(lcal_model_list)==0:
				obslog.info('No leakage corrected models are available.\n')  # Checking if any gaincal model is present or not. If not exit at this stage
				lcal_msg=2
			else:
				lcal_msg,global_database_list,local_database_list=make_final_leakcal(msname,basedir,caltable_name_prefix,gaintable=caltable_for_local_database,\
								freqavg=freq_avg,timeavg=time_avg,ref_ant=ref_ant,gain_minsnr=gain_minsnr,modellist=lcal_model_list)
				caltable_for_local_database+=local_database_list
				caltable_for_global_database+=global_database_list	

		if gcal_msg==0 and lcal_msg==0:
			obslog.info('Searching polarisation calibration models in model directory : '+polcal_modeldir+'\n')
			polcalmodel_list=glob.glob(polcal_modeldir+'/*.model')
			for i in polcalmodel_list:
				if 'leakage' in i:
					polcalmodel_list.remove(i)
			if len(polcalmodel_list)==0:
				obslog.info('No polarisation calibration models are available.\n')
				pcal_msg=2
			else:
				pcal_msg,global_database_list,local_database_list=make_final_polcal(msname,metafits,basedir,caltable_name_prefix,\
									gaintable=caltable_for_local_database,freqavg=freq_avg,timeavg=time_avg,pol_skip_freq=pol_skip_freq,modellist=polcalmodel_list)
				caltable_for_local_database+=local_database_list
				caltable_for_global_database+=global_database_list
		else:
			if gcal_msg!=0:
				if gcal_msg==2:
					obslog.info('No gaincal models are present.\n')
				else:
					obslog.info('Error in gaincal.\n')	
			elif lcal_msg!=0:
				if lcal_msg==2:
					obslog.info('No leakage corrected gaincal models are present.\n')
				else:
					obslog.info('Error in leakage corrected gaincal.\n')
			else:
				obslog.info('Other error occured.\n')

		if len(caltable_for_local_database)!=0:
			for i in caltable_for_local_database:				
				os.system('cp -r '+i+' '+localdatabase) # Copying to local database
		if len(caltable_for_global_database)!=0:
			for j in caltable_for_global_database:
				os.system('cp -r '+j+' '+localdatabase) # Copying to local database for further copying to global database
		if len(caltable_for_local_database)!=0:
			for i in caltable_for_local_database:				
				os.system('rm -rf '+i) # Deleting to local database
		if len(caltable_for_global_database)!=0:
			for j in caltable_for_global_database:
				os.system('rm -rf '+j) # Deleting to local database for further copying to global database	
		return gcal_msg,caltable_for_global_database,caltable_for_local_database

def final_imaging_for_database(msname,metafits,basedir,casacals=[],calibratecals=[],residual_frac=0.1,\
		quality_factor=1,inputfile='',localdatabase='',savedir='',cutoutbox='',want_automask=False,\
		savemodel=False,saveres=False,use_ankflag=False,do_pol=False,mask='',freq_interval=160,\
		time_interval=10.0,freq_width=160,time_width=10,sigma=10,thresh=[0.1]):
	'''
	Function to make final images for database
	Parameters:
	msname = Name of the measurement set
	workdir = Name of the working directory
	casacals = CASA caltables to apply
	calibratecals = CALIBRATE caltables to apply
	Return:
	0 on success, 1 on failure
	'''
	if len(casacals)!=0:
		casa_caltables=','.join(casacals)
	else:
		casa_caltables=''
	if len(calibratecals)!=0:
		calibrate_caltables=','.join(calibratecals)
	else:
		calibrate_caltables=''
	if basedir[-1]=='/':
		basedir=basedir[:-1]
	if casa_caltables=='' and calibrate_caltables=='':
		obslog.info('parallel_ms_split --msname '+msname+' --savedir '+basedir+' --freq_interval '+str(freq_interval)+' --time_interval '+str(time_interval)+' --freq_width '+\
						str(freq_width)+' --time_width '+str(time_width)+' --datacolumn corrected\n')
		a=os.system('parallel_ms_split --msname '+msname+' --savedir '+basedir+' --freq_interval '+str(freq_interval)+' --time_interval '+str(time_interval)+' --freq_width '+\
						str(freq_width)+' --time_width '+str(time_width)+' --datacolumn corrected')
	elif casa_caltables!='' and calibrate_caltables=='':
		obslog.info('parallel_ms_split --msname '+msname+' --savedir '+basedir+' --freq_interval '+str(freq_interval)+' --time_interval '+str(time_interval)+' --freq_width '+\
						str(freq_width)+' --time_width '+str(time_width)+' --datacolumn corrected --casa_caltables '+casa_caltables)
		a=os.system('parallel_ms_split --msname '+msname+' --savedir '+basedir+' --freq_interval '+str(freq_interval)+' --time_interval '+str(time_interval)+' --freq_width '+\
						str(freq_width)+' --time_width '+str(time_width)+' --datacolumn corrected --casa_caltables '+casa_caltables)
	elif casa_caltables=='' and calibrate_caltables!='':
		obslog.info('parallel_ms_split --msname '+msname+' --savedir '+basedir+' --freq_interval '+str(freq_interval)+' --time_interval '+str(time_interval)+' --freq_width '+\
						str(freq_width)+' --time_width '+str(time_width)+' --datacolumn corrected --calibrate_caltables '+calibrate_caltables)
		a=os.system('parallel_ms_split --msname '+msname+' --savedir '+basedir+' --freq_interval '+str(freq_interval)+' --time_interval '+str(time_interval)+' --freq_width '+\
						str(freq_width)+' --time_width '+str(time_width)+' --datacolumn corrected --calibrate_caltables '+calibrate_caltables)
	elif casa_caltables!='' and calibrate_caltables!='':
		obslog.info('parallel_ms_split --msname '+msname+' --savedir '+basedir+' --freq_interval '+str(freq_interval)+' --time_interval '+str(time_interval)+' --freq_width '+\
						str(freq_width)+' --time_width '+str(time_width)+' --datacolumn corrected --casa_caltables '+casa_caltables+' --calibrate_caltables '+calibrate_caltables)
		a=os.system('parallel_ms_split --msname '+msname+' --savedir '+basedir+' --freq_interval '+str(freq_interval)+' --time_interval '+str(time_interval)+' --freq_width '+\
						str(freq_width)+' --time_width '+str(time_width)+' --datacolumn corrected --casa_caltables '+casa_caltables+' --calibrate_caltables '+calibrate_caltables)
	if os.WEXITSTATUS(a)!=0:
		return 1
	else:
		obslog.info('Spliting of ms is done.\n')
		splited_ms=glob.glob(basedir+'/splited_ms/*.ms')
		splited_ms_copy=copy.deepcopy(splited_ms)
		thresh=[str(i) for i in thresh]
		if do_pol==True:
			stokes='IQUV'
			if len(thresh)!=4:
				thresh=thresh*4
		else:
			if len(thresh)>1:
				thresh=thresh[0]
			stokes='psudoI'
		threshold=','.join(thresh)
		touch_count=0
		# Estimating total casa instances
		#################################
		total_available_cpu=psutil.cpu_count()-(psutil.cpu_count()*psutil.cpu_percent()/100.0)
		available_cpu_for_paircars=int(total_available_cpu*inputs.cpu_frac)
		total_cpu_frac=int(psutil.cpu_count()*inputs.cpu_frac)
		if available_cpu_for_paircars>total_cpu_frac:
			available_cpu_for_paircars=total_cpu_frac
		casa_instance=int(available_cpu_for_paircars/2)
		finished_list=glob.glob(basedir+'/.Finished_final_imaging_'+os.path.basename(msname)+'_*')
		for i in finished_list:
			if 'nometa' in i or 'noms' in i:	
				os.system('rm -rf '+i)
		finished_num=0
		while True:	
			for i in range(casa_instance):
				ms=splited_ms[0]
				if ms[-1]=='/':
					ms=ms[:-1]
				workdir=os.path.dirname(os.path.abspath(ms))+'/'+os.path.basename(ms).split('.ms')[0]
				cmd='final_imaging --msname '+ms+' --metafits '+metafits+' --basedir '+basedir+' --workdir '+workdir+' --savedir '+savedir+' --savemodel '+str(savemodel)\
					+' --saveres '+str(saveres)+' --stokes '+stokes+' --cutoutbox '+cutoutbox+' --sigma '+str(sigma)+' --threshold '+str(threshold)\
					+' --want_automask '+str(want_automask)+' --maskfile '+str(mask)+' --quality_factor '+str(quality_factor)+' --inputfile '+str(inputfile)\
						+' --use_ankflag '+str(use_ankflag)+' --residual_frac '+str(residual_frac)
				screen_name='Final_imaging_'+os.path.basename(ms).split('.ms')[0]
				finished_touch_file=basedir+'/.Finished_final_imaging_'+os.path.basename(ms)
				if finished_touch_file in finished_list:
					obslog.info('Imaging is already completed for ms : '+os.path.basename(ms))
				else:
					casa_instance_runner(cmd,screen_name,basedir,finished_touch_file,2,casa_instance)
					touch_count+=1
				splited_ms.remove(ms)
				if len(splited_ms)==0:
					break
			obslog.info('Waiting for final imaging spawned to be finished for '+str(touch_count-finished_num)+' jobs......\n')
			while True:
				time.sleep(2.0)
				touch_files=len(glob.glob(basedir+'/.Finished_final_imaging_'+os.path.basename(msname)+'_*'))
				if (touch_files-finished_num)>=1:
					finished_num+=touch_files
					break
			if len(glob.glob(basedir+'/.Finished_final_imaging_'+os.path.basename(msname)+'_*'))==len(splited_ms_copy):
				break
		return 0


from optparse import OptionParser
if __name__=='__main__':
	usage= ' PAIRCARS database manager'
	parser = OptionParser(usage=usage)
	parser.add_option('--msname',dest="msname",default=None,help="Name of the measurement set",metavar="Measurement set")
	parser.add_option('--metafits',dest="metafits",default=None,help="Name of the metafits file",metavar="Metafits file")
	parser.add_option('--num_jobs',dest="num_jobs",default=0,help="Total number of calibration jobs spawned for this measurement set",metavar="Integer")
	parser.add_option('--basedir',dest="basedir",default=None,help="Name of base directory for a given day",metavar="Directory path")
	parser.add_option('--gaincal_modeldir',dest="gaincal_modeldir",default=None,help="Name of gaincal model directory",metavar="Directory path")
	parser.add_option('--bandpass_modeldir',dest="bandpass_modeldir",default=None,help="Name of bandpass model directory",metavar="Directory path")
	parser.add_option('--polcal_modeldir',dest="polcal_modeldir",default=None,help="Name of polarisation model directory",metavar="Directory path")
	parser.add_option('--localdatabase',dest="localdatabase",default=None,help="Name of local database",metavar="Directory path")
	parser.add_option('--freqavg',dest="freqavg",default=160,help="Frequency averaging during calibration in kHz",metavar="Float")
	parser.add_option('--timeavg',dest="timeavg",default=2.0,help="Time averaging during calibration in second",metavar="Float")
	parser.add_option('--cal_obsid',dest="cal_obsid",default=None,help="Caltable Observation ID",metavar="Integer")
	parser.add_option('--inputfile',dest='inputfile',default=None,help='Path of the P-AIRCARS input file',metavar="File path")
	
	(options, args) = parser.parse_args()
	
if os.path.isdir(str(options.msname))==False:
	print ('Measurement set is not present.\n')
	os._exit(1)
	
if options.basedir==None:
	print ('No base directory path is given. Please give base directory path to continue.\n')
	os._exit(1)
elif os.path.isdir(str(options.basedir))==False:
	os.makedirs(str(options.basedir))

if options.inputfile==None or os.path.isfile(str(options.inputfile))==False:
	print ('P-AIRCARS input file is not present. Please give the input file and re run.\n')
	os._exit(1)

os.chdir(str(options.basedir))
inputfile=str(options.inputfile)
if inputfile[-1]=='/':
	inputfile=inputfile[:-1]
sys.path.append(os.path.dirname(os.path.abspath(inputfile)))
import selfcal_inputs as inputs
if inputs.quality_factor==0:
	polskip=5.12
elif inputs.quality_factor==1:
	polskip=2.56
else:
	polskip=1.28

AM=AccessMS(str(options.msname))
freqs=AM.get_freqs()/10**6
coarse_chan_0=freq_to_MWA_coarse(freqs[0])
coarse_chan_1=freq_to_MWA_coarse(freqs[-1])
OBSID=int(fits.getheader(str(options.metafits))['GPSTIME'])
formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
obslog = logging.getLogger(str(OBSID)+'_log')
obslog.setLevel(logging.DEBUG)
if inputs.verbose==True:
	console=logging.StreamHandler(sys.stdout)
	console.setFormatter(formatter)
	obslog.addHandler(console)
if os.path.exists(str(options.basedir)+'/'+str(OBSID)+'_obslog.log'):
	os.system('rm -rf '+str(options.basedir)+'/'+str(OBSID)+'_obslog.log')
filehandle=logging.FileHandler(str(options.basedir)+'/'+str(OBSID)+'_obslog.log')
filehandle.setFormatter(formatter)
obslog.addHandler(filehandle)
obslog.propagate = False

if int(options.cal_obsid)==OBSID:
	msg,caltable_for_global_database,caltable_for_local_database=managing_caldatabase(str(options.msname),str(options.metafits),int(options.num_jobs),str(options.basedir),\
	str(options.gaincal_modeldir),str(options.bandpass_modeldir),str(options.polcal_modeldir),str(options.localdatabase),float(inputs.gain_minsnr),str(inputs.ref_ant),\
	freq_avg=float(options.freqavg),time_avg=float(options.timeavg),pol_skip_freq=float(polskip))
else:
	obslog.info('Calibration was not performed for this MS. Applying solutions from nearest caltables.\n')
	caltable_for_local_database=glob.glob(inputs.basedir+'/localdatabase/'+str(int(options.cal_obsid))

'''
if len(caltable_for_global_database)>0:
	calarrays=[]
	calmodes=[]
	for i in caltable_for_global_database:
		calmode=os.path.basename(i).split('.bin')[0].split('.')[-1]
		calarrays.append(np.load(i,allow_pickle=True))	
		calmodes.append(calmode)
	if str(options.localdatabase)[-1]=='/':
		localdatabase=str(options.localdatabase)[:-1]
	else:
		localdatabase=str(options.localdatabase)
	final_caltable=str(localdatabase)+'/'+str(OBSID)+'/'+str(OBSID)+'_'+str(coarse_chan_0)+'_'+str(coarse_chan_1)+'.bin'
	if os.path.exists(final_caltable):	
		os.system('rm -rf '+final_caltable)
	caltable_list=','.join(caltable_for_global_database)
	obslog.info('Compressing caltables........\n')
	os.system('compress_caltables --caltables '+caltable_list+' --compressed_file '+final_caltable)
	if os.path.exists(final_caltable):
		attachments=[final_caltable]
		msg='Dear PAIRCARS developers,\n\nCaltables for : OBSID : '+str(OBSID)+' Coarse channels : '+str(coarse_chan_0)+'-'+\
								str(coarse_chan_1)+'\n\nBest,\nPAIRCARS developing team.\n'
		print ('Sending caltables to database....\n')
		send_msg_code,send_msg=send_to_database('Caltables for : OBSID : '+str(OBSID)+' Coarse channels : '+str(coarse_chan_0)+'-'+\
								str(coarse_chan_1),msg,attachments=attachments)
		if send_msg_code==0:
			for i in caltable_for_global_database:
				if 'pcal' not in i and 'beam' not in i: 
					os.system('rm -rf '+i)
			os.system('rm -rf '+final_caltable)	

# Appying solution to main ms
#############################
casa_gaintable=[]
for i in caltable_for_local_database:
	if 'gcal' in i:
		casa_gaintable.append(i)
	if 'bcal' in i:
		casa_gaintable.append(i)
	if 'lcal' in i:
		casa_gaintable.append(i)
calibrate_gaintable=[]
for i in caltable_for_local_database:
	if 'beam' in i:
		calibrate_gaintable.append(i)
	elif 'pcal' in i:
		calibrate_gaintable.append(i)

if inputs.calc_selfcalib_params==True:
	calparams=CalcParams(str(options.msname),inputs.quality_factor,inputs.safety_factor)
	residual_frac=float(calparams.calc_calibration_params()[2])
if inputs.maskfile!='':
	mask=inputs.maskfile
elif inputs.maskstr!='':
	mask=inputs.maskstr
else:
	mask=''

start_sigma=np.load(str(options.basedir)+'/Ref_time_chan_sigma.npy',allow_pickle=True)[0]
rms_list=np.load(str(options.basedir)+'/Ref_time_chan_sigma.npy',allow_pickle=True)[1]

result=final_imaging_for_database(str(options.msname),str(options.metafits),str(options.basedir),casacals=casa_gaintable,calibratecals=calibrate_gaintable,residual_frac=residual_frac,\
		quality_factor=inputs.quality_factor,inputfile=inputfile,localdatabase=localdatabase,savedir=localdatabase+'/'+str(OBSID)+'/images',cutoutbox='3,3',\
		want_automask=False,savemodel=False,saveres=False,use_ankflag=False,do_pol=inputs.do_polcal,mask=mask,\
		freq_interval=1280,time_interval=30,freq_width=1280,time_width=30,sigma=start_sigma,thresh=rms_list)
'''
'''
result=final_imaging_for_database(str(options.msname),str(options.metafits),basedir,casacals=casa_gaintable,calibratecals=calibrate_gaintable,residual_frac=residual_frac,\
		quality_factor=inputs.quality_factor,inputfile=inputfile,localdatabase=localdatabase,savedir=inputs.savedir,cutoutbox=inputs.cutoutbox,want_automask=inputs.want_auto_masking,\
		savemodel=inputs.savemodel,saveres=inputs.saveresidual,use_ankflag=inputs.use_ankflagger,do_pol=inputs.do_polcal,mask=mask,freq_interval=inputs.image_delta_freq,\
		time_interval=image_delta_time,freq_width=inputs.image_freq,time_width=inputs.image_time,sigma=start_sigma,threshold=rms_thresh)

'''













