'''
Code is written by Devojyoti Kansabanik , 2 May, 2021
'''
import os
from casatools import msmetadata,table,measures,quanta,agentflagger,image,calibrater,ms
from casatasks import *
import sys,logging,numpy as np,copy,glob,psutil,time,multiprocessing as mp,subprocess,getpass
from paircars.basic_func import *
from paircars.access_ms import *
from paircars.decor import *
from paircars.flagger import *
from paircars.fullpol_selfcal_LTS import *
from paircars.libpaircars import send_paircars_notification
from optparse import OptionParser
from astropy.io import fits
from CALIBRATE.access_calibrate import *
from paircars.libpaircars import send_paircars_notification,send_to_database
from multiprocessing import Process
from astropy import wcs

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


def make_final_gaincal(msname,workdir,caltable_name_prefix,freqavg=10,timeavg=0.5,ref_ant='1',gain_minsnr=4,modellist=[],verbose=False):
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
	verbose = False, verbose output or not
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
		freq_list=[]
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
		freq_list.append(f-df)
		freq_list.append(f+df)
		AM=AccessMS(msname)
		freqres=AM.calc_freqres()
		if freqavg<freqres:
			freqavg=freqres
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
		spw=str(min(freq_list))+'~'+str(max(freq_list))+'MHz' # Spectral window range for reference model
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
		t=int(ntimes/len(modellist))
		if t<=0:
			t=1
		obslog.info('cal.calibrate(msname=\''+global_cal_ms+'\',caltable=\''+global_caltable+'\',calmode=\'diag\',minuv='+\
					str(calib_uvrange_min)+',maxuv='+str(calib_uvrange_max)+',quiet=True,t='+\
					str(t)+',j=3,ch=1,verbose='+str(verbose)+')\n') # Making gaincal table for reference time and chan
		cal.calibrate(msname=global_cal_ms,caltable=global_caltable,calmode='diag',minuv=calib_uvrange_min,maxuv=calib_uvrange_max,\
			quiet=True,t=t,j=3,ch=1,verbose=verbose)
		os.system('rm -rf '+local_cal_ms+' '+local_cal_ms+'.flagversions')
		os.system('rm -rf '+global_cal_ms+' '+global_cal_ms+'.flagversions')
		if os.path.islink(msname_copy):
			os.system('unlink '+msname_copy)
		else:
			os.system('rm -rf '+msname_copy+'*')
		caltable_for_global_database.append(global_caltable)
		caltable_for_local_database.append(local_caltable)
		return 0,caltable_for_global_database,caltable_for_local_database


def make_final_bandpass(msname,workdir,caltable_name_prefix,gaintable=[],freqavg=10,timeavg=0.5,ref_ant='1',gain_minsnr=4,modellist=[],verbose=False):
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
	verbose = False, verbose output or not
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
		freq_list=[]
		for i in range(len(modellist)):
			modelname=modellist[i]
			modelbasename=os.path.basename(modelname)
			timestamp='/'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[:3])+'/'+':'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[3:])
			f=imhead(imagename=modelname,mode='list')['crval4']/10**6
			df=(imhead(imagename=modelname,mode='list')['cdelt4']/10**6)/2.0
			freq_list.append(f-df)
			freq_list.append(f+df)
			ref_models.append(modelname)
			ref_timestamp=timestamp	
		if len(ref_models)==0:
			obslog.info('No models are present.\n')
			return 1,caltable_for_global_database,caltable_for_local_database
		AM=AccessMS(msname)
		freqres=AM.calc_freqres()
		if freqavg<freqres:
			freqavg=freqres
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
								'\',applymode=\'calflag\',calwt=[False],flagbackup=False)\n')			
			applycal(vis=msname,gaintable=gaintable,timerange=timerange,applymode='calflag',calwt=[False],flagbackup=False)
			datacolumn_to_split='corrected'
		else:
			datacolumn_to_split='data'

		# Spliting ms for local and global caltable only at reference time 
		##################################################################
		spw=str(min(freq_list))+'~'+str(max(freq_list))+'MHz'
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
		freqs_coarse=np.array([freq_to_MWA_coarse(freq/10**6) for freq in freqs])
		for i in range(len(modellist)):
			modelname=modellist[i]
			modelbasename=os.path.basename(modelname)
			imhead(imagename=modelname,mode='put',hdkey='cdelt4',hdvalue=str(2.56*10**6))
			timestamp='/'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[:3])+'/'+':'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[3:])
			f=imhead(imagename=modelname,mode='list')['crval4']/10**6
			model_coarse=freq_to_MWA_coarse(f)
			pos=np.where(freqs_coarse==model_coarse)[0]
			spw=str(np.min(pos))+'~'+str(np.max(pos))
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
					str(calib_uvrange_min)+',maxuv='+str(calib_uvrange_max)+',quiet=True,t='+\
					str(int(ntimes))+',j=3,ch=1,verbose='+str(verbose)+')\n') # Making bandpass table for reference time
		cal.calibrate(msname=global_cal_ms,caltable=global_caltable,calmode='diag',minuv=calib_uvrange_min,maxuv=calib_uvrange_max,\
			quiet=True,t=int(ntimes),j=3,ch=1,verbose=verbose)
		os.system('rm -rf '+local_cal_ms+' '+local_cal_ms+'.flagversions')
		os.system('rm -rf '+global_cal_ms+' '+global_cal_ms+'.flagversions')
		if os.path.islink(msname_copy):
			os.system('unlink '+msname_copy)
		else:
			os.system('rm -rf '+msname_copy+'*')
		caltable_for_global_database.append(global_caltable)
		caltable_for_local_database.append(local_caltable)
		return 0,caltable_for_global_database,caltable_for_local_database

def make_final_leakcal(msname,workdir,caltable_name_prefix,gaintable=[],freqavg=10,timeavg=0.5,ref_ant='1',gain_minsnr=4,modellist=[],verbose=False):
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
	verbose = False, verbose output or not
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
		freq_list=[]
		for i in range(len(modellist)):
			modelname=modellist[i]
			modelbasename=os.path.basename(modelname)
			timestamp='/'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[:3])+'/'+':'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[3:])
			f=imhead(imagename=modelname,mode='list')['crval4']/10**6
			df=(imhead(imagename=modelname,mode='list')['cdelt4']/10**6)/2.0
			freq_list.append(f-df)
			freq_list.append(f+df)
			ref_models.append(modelname)
			ref_timestamp=timestamp	
		if len(ref_models)==0:
			obslog.info('No models are present.\n')
			return 1,caltable_for_global_database,caltable_for_local_database
		AM=AccessMS(msname)
		freqres=AM.calc_freqres()
		if freqavg<freqres:
			freqavg=freqres
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
								'\',applymode=\'calflag\',calwt=[False],flagbackup=False)\n')			
			applycal(vis=msname,gaintable=gaintable,timerange=timerange,applymode='calflag',calwt=[False],flagbackup=False)
			datacolumn_to_split='corrected'
		else:
			datacolumn_to_split='data'
		
		# Spliting ms for local and global caltable for reference time only
		###################################################################
		spw=str(min(freq_list))+'~'+str(max(freq_list))+'MHz'
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
		freqs_coarse=np.array([freq_to_MWA_coarse(freq/10**6) for freq in freqs])
		for i in range(len(modellist)):
			modelname=modellist[i]
			modelbasename=os.path.basename(modelname)
			imhead(imagename=modelname,mode='put',hdkey='cdelt4',hdvalue=str(2.56*10**6))
			timestamp='/'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[:3])+'/'+':'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[3:])
			f=imhead(imagename=modelname,mode='list')['crval4']/10**6
			model_coarse=freq_to_MWA_coarse(f)
			pos=np.where(freqs_coarse==model_coarse)[0]
			spw=str(np.min(pos))+'~'+str(np.max(pos))
			obslog.info('ft(vis=\''+local_cal_ms+'\',model=\''+modelname+'\',spw=\'0:'+str(spw)+'\',usescratch=True)\n')
			ft(vis=local_cal_ms,model=modelname,spw='0:'+str(spw),usescratch=True)

		AM=AccessMS(local_cal_ms)
		model_chan,nomodel_chan=AM.get_model_nomodel_chan()		
		model_spw='0:'
		for i in model_chan:
			model_spw+=str(i)+';'
		model_spw=model_spw[:-1]

		# Perform leakage corrected differential bandpass in a single caltable
		######################################################################
		IB=ImageBasic(local_cal_ms)
		uvrange_to_cal=IB.calc_calib_uvrange(12)[0]
	
		obslog.info('bandpass(vis=\''+local_cal_ms+'\',caltable='+str(local_caltable)+'\',refant=\''\
				+str(ref_ant)+'\',minsnr='+str(gain_minsnr)+',uvrange=\''+uvrange_to_cal+'\')')
		bandpass(vis=local_cal_ms,caltable=local_caltable,refant=str(ref_ant),minsnr=gain_minsnr,uvrange=uvrange_to_cal)
																							 # Performing leakage corrected differential bandpass
		obslog.info('applycal(vis=\''+local_cal_ms+'\',gaintable=[\''+local_caltable+'\'],applymode=\'calflag\',flagbackup=False,calwt=[False])\n')
		applycal(vis=local_cal_ms,gaintable=[local_caltable],applymode='calflag',flagbackup=False,calwt=[False])
		obslog.info('flagdata(vis=\''+local_cal_ms+'\',mode=\'rflag\',datacolumn=\'corrected\',timedevscale=10.0,freqdevscale=7.0,flagbackup=False)\n')
		flagdata(vis=local_cal_ms,mode='rflag',datacolumn='corrected',timedevscale=10.0,freqdevscale=7.0,flagbackup=False)
		if os.path.isdir(local_caltable):
			os.system('rm -rf '+local_caltable)
		obslog.info('bandpass(vis=\''+local_cal_ms+'\',caltable='+str(local_caltable)+'\',refant=\''\
				+str(ref_ant)+'\',minsnr='+str(gain_minsnr)+',uvrange=\''+uvrange_to_cal+'\')')
		bandpass(vis=local_cal_ms,caltable=local_caltable,refant=str(ref_ant),minsnr=gain_minsnr,uvrange=uvrange_to_cal)
																							 # Performing leakage corrected differential bandpass
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
					str(calib_uvrange_min)+',maxuv='+str(calib_uvrange_max)+',quiet=True,j=3,ch=1,verbose='+str(verbose)+')\n') 
							# Making leakage corrected bandpass table for reference time and chan
		cal.calibrate(msname=global_cal_ms,caltable=global_caltable,calmode='diag',minuv=calib_uvrange_min,maxuv=calib_uvrange_max,quiet=True,j=3,ch=1,verbose=verbose)
		#os.system('rm -rf '+local_cal_ms+' '+local_cal_ms+'.flagversions')
		#os.system('rm -rf '+global_cal_ms+' '+global_cal_ms+'.flagversions')
		if os.path.islink(msname_copy):
			os.system('unlink '+msname_copy)
		else:
			os.system('rm -rf '+msname_copy+'*')
		caltable_for_global_database.append(global_caltable)
		caltable_for_local_database.append(local_caltable)
		return 0,caltable_for_global_database,caltable_for_local_database

def make_final_polcal(msname,metafits,workdir,caltable_name_prefix,gaintable=[],freqavg=10,timeavg=0.5,pol_skip_freq=1280,modellist=[],verbose=False):
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
	verbose = False, verbose output or not
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
		freq_list=[]
		for i in range(len(modellist)):
			modelname=modellist[i]
			modelbasename=os.path.basename(modelname)
			timestamp='/'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[:3])+'/'+':'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[3:])
			f=imhead(imagename=modelname,mode='list')['crval4']/10**6
			df=(imhead(imagename=modelname,mode='list')['cdelt4']/10**6)/2.0
			freq_list.append(f-df)
			freq_list.append(f+df)
			ref_models.append(modelname)
			ref_timestamp=timestamp	
		if len(ref_models)==0:
			obslog.info('No models are present.\n')
			return 1,caltable_for_global_database,caltable_for_local_database
		AM=AccessMS(msname)
		freqres=AM.calc_freqres()
		if freqavg<freqres:
			freqavg=freqres
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
								'\',applymode=\'calflag\',calwt=[False],flagbackup=False)\n')			
			applycal(vis=msname,gaintable=gaintable,timerange=timerange,applymode='calflag',calwt=[False],flagbackup=False)
			datacolumn_to_split='corrected'
		else:
			datacolumn_to_split='data'

		# Spliting ms for local and global caltable, for global caltable only reference time solution will be there
		###########################################################################################################
		spw=str(min(freq_list))+'~'+str(max(freq_list))+'MHz'
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

		# Applying cross-hand phase correction
		######################################
		mwa_config=get_MWA_phase(metafits) # TODO : Include from cross phase cal solutions
		PSC=PolSelfcal(beamcal_ms,metafits,32*60,verbose=False,interactive=False,savelog=False) # Performing ideal beam correction at phase center for every coarse channels
		if mwa_config=='MWAPhaseI':
			crossphase=15
		elif mwa_config=='MWAPhaseIILB' or mwa_config=='MWAPhaseIICOMPACT':
			crossphase=135
		obslog.info('Applying cross hand phase solution. Cross hand phase : '+str(crossphase)+' deg.\n')
		PSC.apply_cross_hand_phase(cross_phase=crossphase,caltable='',polbasis='Linear',modify_datacolumn=True,datacolumn='DATA')
	
		# Perform ideal beam correction
		###############################
		obslog.info('Performing ideal beam correction......\n')
		obslog.info('###########################################\n')
		obslog.info('PolSelfcal(\''+beamcal_ms+'\',\''+metafits+'\','+str(32*60)+',verbose=False,interactive=False)\n')
		obslog.info('PSC.correct_visibility_single_beam_jones(modify_datacolumn=False,skip_freq='+str(pol_skip_freq)+',save_beamfile=\''+str(beam_caltable)+'\')\n')
		PSC.correct_visibility_single_beam_jones(modify_datacolumn=False,skip_freq=float(pol_skip_freq),save_beamfile=beam_caltable)
		caltable_for_global_database.append(beam_caltable)
		caltable_for_local_database.append(beam_caltable)
		obslog.info('flagdata(vis=\''+beamcal_ms+'\',mode=\'rflag\',datacolumn=\'corrected\',timedevscale=7,freqdevscale=7,flagbackup=False)\n')
		flagdata(vis=beamcal_ms,mode='rflag',datacolumn='corrected',timedevscale=7,freqdevscale=7,flagbackup=False)

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
		freqs_coarse=np.array([freq_to_MWA_coarse(freq/10**6) for freq in freqs])
	
		for i in modellist:
			if 'leakage' in i:
				modellist.remove(i)

		for i in range(len(modellist)):
			modelname=modellist[i]
			modelbasename=os.path.basename(modelname)
			imhead(imagename=modelname,mode='put',hdkey='cdelt4',hdvalue=str(2.56*10**6))
			timestamp='/'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[:3])+'/'+':'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[3:])
			f=imhead(imagename=modelname,mode='list')['crval4']/10**6
			model_coarse=freq_to_MWA_coarse(f)
			pos=np.where(freqs_coarse==model_coarse)[0]
			spw=str(np.min(pos))+'~'+str(np.max(pos))
			obslog.info('ft(vis=\''+polcal_ms+'\',model=\''+modelname+'\',spw=\'0:'+str(spw)+'\',usescratch=True)\n')
			ft(vis=polcal_ms,model=modelname,spw='0:'+str(spw),usescratch=True)
		
		# Making polcal caltable in CALIBRATE numpy format
		##################################################
		AMgcal=AccessMS(polcal_ms)
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
					str(calib_uvrange_min)+',maxuv='+str(calib_uvrange_max)+',quiet=True,j=3,ch=1,verbose='+str(verbose)+')\n') # Making gaincal table for reference time and chan
		cal.calibrate(msname=polcal_ms,caltable=polcal_caltable,minuv=calib_uvrange_min,maxuv=calib_uvrange_max,quiet=True,j=3,ch=1,verbose=verbose)
		#os.system('rm -rf '+beamcal_ms+' '+beamcal_ms+'.flagversions')
		#os.system('rm -rf '+polcal_ms+' '+polcal_ms+'.flagversions')
		if os.path.islink(msname_copy):
			os.system('unlink '+msname_copy)
		else:
			os.system('rm -rf '+msname_copy+'*')
		caltable_for_global_database.append(polcal_caltable)
		caltable_for_local_database.append(polcal_caltable)
		return 0,caltable_for_global_database,caltable_for_local_database

def managing_caldatabase(msname,metafits,total_spawned_jobs,basedir,gaincal_modeldir,bandpass_modeldir,polcal_modeldir,localdatabase,gain_minsnr,\
						ref_ant,freq_avg=40.0,time_avg=0.5,pol_skip_freq=1.28,verbose=False):
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
	verbose = False, print verbose output
	Return:
	Message code, global database caltable list, local database caltable list
	'''
	# Setting up logs and waiting for calibrations jobs to finish
	#############################################################
	cal=CALIBRATE()
	OBSID=int(fits.getheader(metafits)['GPSTIME'])
	if verbose==False:
		print ('Waiting for spawning jobs for ms : '+msname+'\n')
	obslog.info('Waiting for spawning jobs for ms : '+msname+'\n')
	caltable_for_global_database=[]
	caltable_for_local_database=[]

	if basedir[-1]=='/':
		basedir=basedir[:-1]
	if msname[-1]=='/':
		msname=msname[:-1]

	basemsdir=os.path.basename(msname).split('.ms')[0] # Base directory for the ms inside model directories
	while True:
		spawned_file=glob.glob(basedir+'/.Finished_spawned_*'+str(OBSID)+'*_*'+os.path.basename(msname).split('.ms')[0]+'*')
		if len(spawned_file)>0:
			spawned_file=spawned_file[0]
			total_spawned_jobs=int(spawned_file.split(os.path.basename(msname).split('.ms')[0]+'_')[-1])
			if verbose==False:
				print ('Waiting for finishing '+str(total_spawned_jobs)+' calibration jobs for ms : '+msname+'\n')
			obslog.info('Waiting for finishing '+str(total_spawned_jobs)+' calibration jobs for ms : '+msname+'\n')
			break

	while True:  # Waiting for all jobs to finish
		touch_files=len(glob.glob(basedir+'/.Finished*cal*'+str(OBSID)+'*'+basemsdir+'*'))
		if touch_files>=total_spawned_jobs:
			if verbose==False:
				print('All calibration jobs for ms : '+msname+' is completed.\n')
			obslog.info('All calibration jobs for ms : '+msname+' is completed.\n')
			break	
		else:
			time.sleep(2.0)

	# Searching for nearest cal directories
	#######################################
	gcal_dirs=str(gaincal_modeldir).split(',')
	bcal_dirs=str(bandpass_modeldir).split(',')
	pcal_dirs=str(polcal_modeldir).split(',')
	
	for i in gcal_dirs:
		if os.path.isdir(i)==False:
			gcal_dirs.remove(i)
	for i in bcal_dirs:
		if os.path.isdir(i)==False:
			bcal_dirs.remove(i)
	for i in pcal_dirs:
		if os.path.isdir(i)==False:
			pcal_dirs.remove(i)
	gcal_obsids=np.array([int(os.path.basename(i)) for i in gcal_dirs])
	bcal_obsids=np.array([int(os.path.basename(i)) for i in bcal_dirs])
	pcal_obsids=np.array([int(os.path.basename(i)) for i in pcal_dirs])

	if len(gcal_dirs)==0 or len(gcal_obsids)==0:
		obslog.info('No gain calibration tables present.\n')
		gaincal_modeldir=''
	else:
		gcal_obsid=gcal_obsids[np.argmin(abs(OBSID-gcal_obsids))]
		gaincal_modeldir=gcal_dirs[np.argmin(abs(OBSID-gcal_obsids))]

	if len(bcal_obsids)!=0:
		bcal_obsid=bcal_obsids[np.argmin(abs(OBSID-bcal_obsids))]
		bandpass_modeldir=bcal_dirs[np.argmin(abs(OBSID-pcal_obsids))]

	if len(pcal_obsids)!=0:
		pcal_obsid=pcal_obsids[np.argmin(abs(OBSID-pcal_obsids))]
		polcal_modeldir=pcal_dirs[np.argmin(abs(OBSID-pcal_obsids))]

	# Setting up different directory names
	######################################
	if gaincal_modeldir[-1]=='/':
		gaincal_modeldir=gaincal_modeldir[:-1]
	if bandpass_modeldir[-1]=='/':
		bandpass_modeldir=bandpass_modeldir[:-1]
	if polcal_modeldir[-1]=='/':
		polcal_modeldir=polcal_modeldir[:-1]
	if localdatabase[-1]=='/':
		localdatabase=localdatabase[:-1]

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
		
	# Final gain calibrations
	#########################
	if verbose==False:
		print ('Making final calibration tables........\n')
	if gaincal_modeldir=='' or os.path.isdir(gaincal_modeldir)==False or len(glob.glob(gaincal_modeldir+'/*.model'))==0:
		if verbose==False:
			print ('No gaincal models are available.\n')
		obslog.info('No gaincal models are available.\n')  # Checking if any gaincal model is present or not. If not exit at this stage
		gcal_msg=2
		return 1,caltable_for_global_database,caltable_for_local_database
	else:
		AM=AccessMS(msname)
		antenna=AM.get_antenna_string()
		unflagchan,flagchan=flag_MWA_coarse(msname,edgewidth=160,do_flag=False,force=False)
		freqs=AM.get_freqs()/10**6
		nchan_avg=int(freq_avg/AM.calc_freqres())
		coarse_chan_0=freq_to_MWA_coarse(freqs[0])
		coarse_chan_1=freq_to_MWA_coarse(freqs[-1])
		if coarse_chan_0!=None and coarse_chan_1!=None:
			pass
		elif coarse_chan_0==None and coarse_chan_1!=None:
			coarse_chan_0=coarse_chan_1
		elif coarse_chan_1==None and coarse_chan_0!=None:
			coarse_chan_1=coarse_chan_0
		else:
			print ('Could not find coarse channels for frequency : '+str(freqs[0])+' and '+str(freqs[-1])+'\n')
			return
		caltable_name_prefix=basedir+'/'+str(OBSID)+'_'+str(coarse_chan_0)+'_'+str(coarse_chan_1) # Caltable name prefix; OBSID_startcoarsechan_endcoarsechan format
		obslog.info('Searching gaincal models in model directory : '+gaincal_modeldir+'\n')
		model_list=glob.glob(gaincal_modeldir+'/*.model') # Gaincal model list
		gcal_msg,global_database_list,local_database_list=make_final_gaincal(msname,basedir,caltable_name_prefix,freqavg=freq_avg,timeavg=time_avg,ref_ant=ref_ant,\
																				gain_minsnr=gain_minsnr,modellist=model_list,verbose=verbose)
		caltable_for_local_database+=local_database_list
		caltable_for_global_database+=global_database_list	
	
		if gcal_msg==0:
			obslog.info('Searching bandpass models in model directory : '+bandpass_modeldir+'\n')
			bpmodel_list=glob.glob(bandpass_modeldir+'/*.model')
			if len(bpmodel_list)==0:
				if verbose==False:
					print ('No bandpass models are available.\n')
				obslog.info('No bandpass models are available.\n')  # Checking if any gaincal model is present or not. If not exit at this stage
				bp_msg=2
			else:
				bp_msg,global_database_list,local_database_list=make_final_bandpass(msname,basedir,caltable_name_prefix,gaintable=caltable_for_local_database,\
								freqavg=freq_avg,timeavg=time_avg,ref_ant=ref_ant,gain_minsnr=gain_minsnr,modellist=bpmodel_list,verbose=verbose)
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
				if verbose==False:
					print ('No leakage corrected models are available.\n')
				obslog.info('No leakage corrected models are available.\n')  # Checking if any gaincal model is present or not. If not exit at this stage
				lcal_msg=2
			else:
				lcal_msg,global_database_list,local_database_list=make_final_leakcal(msname,basedir,caltable_name_prefix,gaintable=caltable_for_local_database,\
								freqavg=freq_avg,timeavg=time_avg,ref_ant=ref_ant,gain_minsnr=gain_minsnr,modellist=lcal_model_list,verbose=verbose)
				caltable_for_local_database+=local_database_list
				caltable_for_global_database+=global_database_list	

		if gcal_msg==0 and lcal_msg==0:
			obslog.info('Searching polarisation calibration models in model directory : '+polcal_modeldir+'\n')
			polcalmodel_list=glob.glob(polcal_modeldir+'/*.model')
			for i in polcalmodel_list:
				if 'leakage' in i:
					polcalmodel_list.remove(i)
			if len(polcalmodel_list)==0:
				if verbose==False:
					print ('No polarisation calibration models are available.\n')
				obslog.info('No polarisation calibration models are available.\n')
				pcal_msg=2
			else:
				pcal_msg,global_database_list,local_database_list=make_final_polcal(msname,metafits,basedir,caltable_name_prefix,\
									gaintable=caltable_for_local_database,freqavg=freq_avg,timeavg=time_avg,pol_skip_freq=pol_skip_freq,modellist=polcalmodel_list,verbose=verbose)
				caltable_for_local_database+=local_database_list
				caltable_for_global_database+=global_database_list
		else:
			if gcal_msg!=0:
				if gcal_msg==2:
					if verbose==False:
						print ('No gaincal models are present.\n')
					obslog.info('No gaincal models are present.\n')
				else:
					if verbose==False:
						print ('Error in gaincal.\n')
					obslog.info('Error in gaincal.\n')	
			elif lcal_msg!=0:
				if lcal_msg==2:
					if verbose==False:
						print ('No leakage corrected gaincal models are present.\n')
					obslog.info('No leakage corrected gaincal models are present.\n')
				else:
					if verbose==False:
						print ('Error in leakage corrected gaincal.\n')
					obslog.info('Error in leakage corrected gaincal.\n')
			else:
				if verbose==False:
					print('Other error occured.\n')
				obslog.info('Other error occured.\n')

		final_local_caltables=[]
		final_global_caltables=[]
		if len(caltable_for_local_database)!=0:
			for i in caltable_for_local_database:				
				os.system('cp -r '+i+' '+localdatabase) # Copying to local database
				final_local_caltables.append(localdatabase+'/'+os.path.basename(i))
		if len(caltable_for_global_database)!=0:
			for j in caltable_for_global_database:
				os.system('cp -r '+j+' '+localdatabase) # Copying to local database for further copying to global database
				final_global_caltables.append(localdatabase+'/'+os.path.basename(j))
		if len(caltable_for_local_database)!=0:
			for i in caltable_for_local_database:				
				os.system('rm -rf '+i) # Deleting to local database
		if len(caltable_for_global_database)!=0:
			for j in caltable_for_global_database:
				os.system('rm -rf '+j) # Deleting to local database for further copying to global database	
		return gcal_msg,final_global_caltables,final_local_caltables

def final_imaging(msname,metafits,basedir,mode,casacals=[],calibratecals=[],residual_frac=0.1,do_diffcal=False,\
		quality_factor=1,inputfile='',localdatabase='',savedir='',cutoutbox='',want_automask=False,\
		savemodel=False,saveres=False,use_ankflag=False,do_pol=False,mask='',freq_interval=160,\
		time_interval=10.0,freq_width=40,time_width=0.5,time_avg=0.5,freq_avg=40,sigma=10,thresh=[0.1],use_wsclean=True):
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
	c=0
	if use_wsclean:
		obslog.info('Calculating final PSF size.....\n')
		maj,minor,pa=get_final_psf(msname,imager='wsclean',weight=inputs.weight,robust=inputs.robust)
	else:
		obslog.info('Calculating final PSF size.....\n')
		maj,minor,pa=get_final_psf(msname,imager='CASA',weight=inputs.weight,robust=inputs.robust)
	while True:
		available_paircars_instance=get_available_paircars_instance(inputs.paircars_dir,inputs.job_id,inputs.instance)
		if c==0:
			obslog.info('Waiting for free cores......\n')
		time.sleep(5)
		c+=1
		if available_paircars_instance>0:
			break	
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
	if savedir=='':
		savedir=basedir+'/Final_images'
	if mode=='Database_Imaging':
		imaging_mode='database'
	elif mode=='Final_Imaging':
		imaging_mode='final'
	AM=AccessMS(msname)
	freqres=AM.calc_freqres()
	timeres=AM.calc_timeres()
	
	if time_avg<=timeres:
		time_avg=timeres
	else:
		time_avg=timeres

	if freq_avg<=freqres:
		freq_avg=freqres
	else:
		freq_avg=freqres
	freqs=AM.get_freqs()/10**6

	ms_cmd_args=['--msname='+msname,'--savedir='+savedir,'--freq_interval='+str(freq_interval),'--time_interval='+str(time_interval),'--freq_width='+str(freq_width),\
		'--time_width='+str(time_width),'--cpu_frac='+str(inputs.cpu_frac)]
	if imaging_mode!='database':
		if inputs.timerange!='':
			ms_cmd_args.append('--time_list='+str(inputs.timerange))
		if inputs.chanrange!='':
			freqrange=[]	
			chanrange=inputs.chanrange.split(',')
			for chan in chanrange:
				s_chan=int(chan.split('~')[0])
				e_chan=int(chan.split('~')[-1])
				freqrange.append(str(freqs[s_chan])+'~'+str(freqs[e_chan]))
			freqrange=','.join(freqrange)
			ms_cmd_args.append('--freq_list='+str(freqrange))
	if time_avg!='':
		ms_cmd_args.append('--time_avg='+str(time_avg))
	if freq_avg!='':
		ms_cmd_args.append('--freq_avg='+str(freq_avg))
	if casa_caltables=='' and calibrate_caltables=='':
		ms_cmd_args.append('--datacolumn=corrected')
	else:
		ms_cmd_args.append('--datacolumn=DATA')

	AM=AccessMS(msname)
	total_bandwidth=AM.calc_bandwidth()/1000.0
	total_time=AM.calc_total_time()
	total_chunks=int((total_bandwidth/freq_interval)*(total_time/time_interval))
	total_blocks=int(total_chunks/20)
	cur_block=0
	pre_casa_instance=0
	casa_instance=0
	if mode=='Final_Imaging':
		while True:
			if cur_block<total_blocks:
				ms_cmd_args.append('--total_block='+str(total_blocks))
				ms_cmd_args.append('--cur_block='+str(cur_block))
				obslog.info('parallel_ms_split '+' '.join(ms_cmd_args)+'\n')
				a=os.system('parallel_ms_split '+' '.join(ms_cmd_args))
				ms_cmd_args=ms_cmd_args[:-2]	
				if os.WEXITSTATUS(a)!=0:
					obslog.info('Error in spliting block '+str(cur_block)+'.\n')
				else:
					splited_ms=glob.glob(savedir+'/splited_ms/'+os.path.basename(msname).split('.ms')[0]+'/*.ms')
					finished_touch_files_for_ms=[basedir+'/.Finished_final_imaging_'+str(imaging_mode)+'_'+os.path.basename(ms) for ms in splited_ms]
					if len(splited_ms)>0:
						obslog.info('Spliting of ms is done for block '+str(cur_block)+'.\n')
						break
					else:
						obslog.info('No ms in block '+str(cur_block)+'\n')
				cur_block+=1
	else:
		obslog.info('parallel_ms_split '+' '.join(ms_cmd_args)+'\n')
		a=os.system('parallel_ms_split '+' '.join(ms_cmd_args))
	if os.WEXITSTATUS(a)!=0:
		return 1
	else:
		obslog.info('Spliting of ms is done.\n')
		splited_ms=glob.glob(savedir+'/splited_ms/'+os.path.basename(msname).split('.ms')[0]+'/*.ms')
		finished_touch_files_for_ms=[basedir+'/.Finished_final_imaging_'+str(imaging_mode)+'_'+os.path.basename(ms) for ms in splited_ms]
		if len(thresh)==2:
			thresh=[np.mean(np.array(thresh))]
		thresh=[str(i) for i in thresh]
		if do_pol==True:
			stokes='IQUV'
			if len(thresh)!=4:
				thresh=thresh*4
		else:
			if len(thresh)>1:
				thresh=thresh[0]
			stokes='pseudoI'
		threshold=','.join(thresh)
		touch_count=0
		ms_done=[]
		# Estimating total casa instances
		#################################
		c=0
		while True:
			available_paircars_instance=get_available_paircars_instance(inputs.paircars_dir,inputs.job_id,inputs.instance)
			if c==0:
				obslog.info('Waiting for free P-AIRCARS instance......\n')
			time.sleep(5)
			c+=1
			if available_paircars_instance>1:
				break
		finished_list=glob.glob(basedir+'/.Finished_final_imaging_'+str(imaging_mode)+'_*success')
		for i in finished_list:
			if 'nometa' in i or 'noms' in i:	
				os.system('rm -rf '+i)
		touch_files=glob.glob(basedir+'/.Finished_final_imaging_'+str(imaging_mode)+'*success')
		error_files=glob.glob(basedir+'/.Finished_final_imaging_'+str(imaging_mode)+'*')
		for e in error_files:
			if e not in touch_files:
				os.system('rm -rf '+e)
		touch_files=[]
		imaging_list=[]
		while True:	
			splited_ms=glob.glob(savedir+'/splited_ms/'+os.path.basename(msname).split('.ms')[0]+'/*.ms')
			splited_ms_copy=copy.deepcopy(splited_ms)
			finished_touch_files_for_ms=[basedir+'/.Finished_final_imaging_'+str(imaging_mode)+'_'+os.path.basename(ms) for ms in splited_ms]
			count=0
			for i in range(len(splited_ms)):
				if count>min(pre_casa_instance,casa_instance):
					break
				ms=splited_ms[0]
				if ms[-1]=='/':
					ms=ms[:-1]
				workdir=os.path.dirname(os.path.abspath(ms))+'/'+os.path.basename(ms).split('.ms')[0]
				cmd_args=['--msname '+ms,' --metafits '+metafits,' --basedir '+basedir,' --workdir '+workdir,' --savedir '+savedir,' --savemodel '+str(savemodel),\
						' --saveres '+str(saveres),' --stokes '+stokes,' --sigma '+str(sigma),' --threshold '+str(threshold),' --wsclean '+str(use_wsclean)\
						,' --want_automask '+str(want_automask),' --quality_factor '+str(quality_factor),' --use_ankflag '+str(use_ankflag),' --imaging_mode '+str(imaging_mode),\
						' --residual_frac '+str(residual_frac),'--cpu_frac '+str(inputs.cpu_frac),'--major_axis '+str(maj),'--minor_axis '+str(minor),'--pa '+str(pa)]
				if do_diffcal==True:
					cmd_args.append(' --do_diffcal '+str(do_diffcal))
				if casa_caltables!='':
					cmd_args.append(' --casa_caltables '+casa_caltables)
				if calibrate_caltables!='':
					cmd_args.append(' --calibrate_caltables '+calibrate_caltables)
				if mask!='':
					cmd_args.append(' --maskfile '+str(mask))
				if cutoutbox!='':
					cmd_args.append(' --cutoutbox '+cutoutbox)
				cmd='final_imaging '+' '.join(cmd_args)
				screen_name=mode+'_'+os.path.basename(ms).split('.ms')[0]+'_'+str(inputs.job_id)
				finished_touch_file=basedir+'/.Finished_final_imaging_'+str(imaging_mode)+'_'+os.path.basename(ms)
				if len(glob.glob(finished_touch_file+'_*'))>0:
					if glob.glob(finished_touch_file+'_*')[0] in finished_touch_files_for_ms:
						touch_files.append(glob.glob(finished_touch_file+'_*')[0])
				imaging_list.append(finished_touch_file)
				if finished_touch_file+'_success' in finished_list:
					obslog.info('Imaging is already completed for ms : '+os.path.basename(ms)+'\n')
					ms_done.append(ms)
					if ms in splited_ms:
						splited_ms.remove(ms)
				elif ms in ms_done:
					obslog.info('Imaging for ms : '+ms+' is going...\n')
					if ms in splited_ms:
						splited_ms.remove(ms)
				elif os.path.exists(finished_touch_file+'_moreflag') or os.path.exists(finished_touch_file+'_error') or len(glob.glob(finished_touch_file+'_*'))>0:
					obslog.info('Imaging already attempted. Either more data flagged or error occured.\n')
					if ms in splited_ms:
						splited_ms.remove(ms)
				else:
					if 'temp_aocal' in ms or 'chan' in ms or 'avg' in ms:
						splited_ms.remove(ms)
						continue
					else:
						obslog.info(cmd+'\n')
						batch_file=paircars_instance_runner(cmd,basedir,inputs.paircars_dir,screen_name,finished_touch_file,inputs.job_id)
						if CPULIMIT_check()==0 and inputs.use_wsclean==False:
							obslog.info('cpulimit --limit 500 -z sh '+batch_file)
							os.system('cpulimit --limit 500 -z sh '+batch_file)
						else:
							obslog.info('sh '+batch_file)
							os.system('sh '+batch_file) 
						time.sleep(0.5)
						splited_ms.remove(ms)
						ms_done.append(ms)
						casa_instance+=1
						count+=1
				touch_count=len(touch_files)
			pre_touch_files=copy.deepcopy(touch_files)
			pre_touch_count=touch_count
			if (len(splited_ms)==0 and mode!='Final_Imaging') or (len(splited_ms)==0 and mode=='Final_Imaging' and cur_block>=total_blocks):
				break
			c=0
			while True:
				for ms in splited_ms_copy:
					finished_touch_file=basedir+'/.Finished_final_imaging_'+str(imaging_mode)+'_'+os.path.basename(ms)
					if len(glob.glob(finished_touch_file+'_*'))>0:
						if glob.glob(finished_touch_file+'_*')[0] not in pre_touch_files and glob.glob(finished_touch_file+'_*')[0] in finished_touch_files_for_ms:
							pre_touch_files.append(glob.glob(finished_touch_file+'_*')[0])
				pre_touch_count=len(pre_touch_files)
				available_paircars_instance=get_available_paircars_instance(inputs.paircars_dir,inputs.job_id,inputs.instance)
				if c==0:
					obslog.info('Waiting for free P-AIRCARS instance......\n')
				if available_paircars_instance>1:
					pre_casa_instance=(pre_touch_count-touch_count)
					if mode=='Final_Imaging' and cur_block<total_blocks:
						obslog.info('Splilting block no : '+str(cur_block)+'\n')
						cmds=['parallel_ms_split']+ms_cmd_args+['--total_block='+str(total_blocks)]+['--cur_block='+str(cur_block)]
						os.system(' '.join(cmds))
						cur_block+=1
					pre_touch_files=[]
					for i in imaging_list:
						if len(glob.glob(i+'_*')):
							if glob.glob(i+'_*')[0] in finished_touch_files_for_ms:
								pre_touch_files.append(glob.glob(i+'_*')[0])
							if mode=='Final_Imaging':
								fil=open(savedir+'/Final_Image_logs.log','a')
								x=os.path.basename(glob.glob(i+'_*')[0]).split('.Finished_final_imaging_'+str(imaging_mode)+'_')[-1].split('.ms')[0]
								y=os.path.basename(glob.glob(i+'_*')[0]).split('.Finished_final_imaging_'+str(imaging_mode)+'_')[-1].split('.ms_')[1]
								fil.write(x+' : '+y+'\n')
								fil.seek(0)
								fil.close()
					break		
				time.sleep(2.0)
				available_paircars_instance=get_available_paircars_instance(inputs.paircars_dir,inputs.job_id,inputs.instance)
				c+=1
			if (len(splited_ms)==0 and mode!='Final_Imaging') or (len(splited_ms)==0 and mode=='Final_Imaging' and cur_block>=total_blocks):
				break
		return 0

def get_quicklook_image(imagename,outfile,freq,timestamp,field_of_view=2): 
	'''
	Function to get a quick look image
	Parameters:
	imagename = Name of the CASA image
	outfile = Output file name
	freq = Frequency in MHz
	timestamp = Timestamp string
	field_of_view = Field of view to cut the image in degree (default : 2)
	Return:
	Outfile name
	'''
	os.system('cp -r '+imagename+' '+'quick_look_'+os.path.basename(imagename))
	imagename='quick_look_'+os.path.basename(imagename)
	org_image=copy.deepcopy(imagename)
	header=imhead(imagename=imagename,mode='list')
	xcent=int(header['shape'][0]/2)
	ycent=int(header['shape'][1]/2)
	cell=np.rad2deg(abs(header['cdelt2'])) # In degree
	freq="{:.2f}".format(float(freq))
	xwidth=ywidth=int((field_of_view)/cell)
	box=str(xcent-int(xwidth/2))+','+str(ycent-int(ywidth/2))+','+str(xcent+int(xwidth/2))+','+str(ycent+int(ywidth/2))
	try:
		header=fits.getheader(imagename)
		if header['NAXIS']==4:
			if header['CTYPE3']=='STOKES':
				stokes_length=header['NAXIS3']
			elif header['CTYPE4']=='STOKES':
				stokes_length=header['NAXIS4']
		else:
			stokes_length=1
		importfits(fitsimage=imagename,imagename=imagename.split('.fits')[0]+'.image')
		imagename=imagename.split('.fits')[0]+'.image'
	except:
		header=imhead(imagename=imagename)
		stokes_axis=np.where(header['axisnames']=='Stokes')[0][0]
		stokes_length=header['shape'][stokes_axis]
	if stokes_length==1:
		stokes_list=['I']
	elif stokes_length==4:
		stokes_list=['I','Q','U','V']
	else:
		print ('Stokes axes are not I or IQUV.\n')
	fig = plt.figure(figsize=(8,8))
	plt.subplots_adjust(wspace=0.45, hspace=0.1)
	for i in range(len(stokes_list)):
		stokes=stokes_list[i]
		try:
			imsubimage(imagename=imagename,outfile='temp_'+stokes+'_'+os.path.basename(imagename)+'.image',box=box,stokes=stokes)
		except:
			return 
		exportfits(imagename='temp_'+stokes+'_'+os.path.basename(imagename)+'.image',fitsimage='temp_'+stokes+'_'+os.path.basename(imagename)+'.fits',dropdeg=True,dropstokes=True)
		data=fits.getdata('temp_'+stokes+'_'+os.path.basename(imagename)+'.fits')
		wlist=fits.getheader('temp_'+stokes+'_'+os.path.basename(imagename)+'.fits')
		w = wcs.WCS(wlist)
		if i==0:
			ax1 = fig.add_subplot(221, projection = w)
			im=ax1.imshow(data,cmap='seismic',origin='lower')
			fig.colorbar(im, ax=ax1)
			ax1.set_title('Stokes : '+stokes)
			ax1.set_xlabel('RA')
			ax1.set_ylabel('DEC')
		elif i==1:		
			ax2 = fig.add_subplot(222, projection = w)
			im=ax2.imshow(data,cmap='seismic',origin='lower')
			fig.colorbar(im, ax=ax2)
			ax2.set_title('Stokes : '+stokes)
			ax2.set_xlabel('RA')
			ax2.set_ylabel('DEC')
		elif i==2:
			ax3 = fig.add_subplot(223, projection = w)
			im=ax3.imshow(data,cmap='seismic',origin='lower')
			fig.colorbar(im, ax=ax3)
			ax3.set_title('Stokes : '+stokes)
			ax3.set_xlabel('RA')
			ax3.set_ylabel('DEC')
		elif i==3:
			ax4 = fig.add_subplot(224, projection = w)
			im=ax4.imshow(data,cmap='seismic',origin='lower')
			fig.colorbar(im, ax=ax4)
			ax4.set_title('Stokes : '+stokes)
			ax4.set_xlabel('RA')
			ax4.set_ylabel('DEC')
	title='Frequency : '+str(freq)+' MHz, Timestamp : '+str(timestamp)+' UTC'
	plt.suptitle(title,fontsize=12)	
	cwd=os.getcwd()
	outfile_dir=os.path.dirname(outfile)
	if outfile_dir=='':
		outfile=cwd+'/'+outfile
	plt.savefig(outfile)
	os.system('rm -rf casa*log temp_*'+os.path.basename(imagename)+'* '+imagename+' '+org_image+' '+imagename.split('.fits')[0]+'.image')
	return outfile

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
	parser.add_option('--freqavg',dest="freqavg",default=10,help="Frequency averaging during calibration in kHz",metavar="Float")
	parser.add_option('--timeavg',dest="timeavg",default=0.5,help="Time averaging during calibration in second",metavar="Float")
	parser.add_option('--inputfile',dest='inputfile',default=None,help='Path of the P-AIRCARS input file',metavar="File path")
	parser.add_option('--verbose',dest='verbose',default=False,help='Verbose output',metavar="Boolean")
	parser.add_option('--wsclean',dest="use_wsclean",default=True,help="Use WSClean for imaging or not",metavar="Boolean")
	(options, args) = parser.parse_args()
#try:	
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
polskip=1.28

AM=AccessMS(str(options.msname))
freqs=AM.get_freqs()/10**6
coarse_chan_0=freq_to_MWA_coarse(freqs[0])
coarse_chan_1=freq_to_MWA_coarse(freqs[-1])
OBSID=int(fits.getheader(str(options.metafits))['GPSTIME'])
is_internet=check_internet()
formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
obslog = logging.getLogger(str(OBSID)+'_log')
obslog.setLevel(logging.DEBUG)
if eval(str(options.verbose))==True:
	console=logging.StreamHandler(sys.stdout)
	console.setFormatter(formatter)
	obslog.addHandler(console)
if os.path.exists(str(options.basedir)+'/logs/'+str(OBSID)+'/'+str(OBSID)+'_obslog.log'):
	os.system('rm -rf '+str(options.basedir)+'/logs/'+str(OBSID)+'/'+str(OBSID)+'_obslog.log')
if inputs.keep_logger==True:
	filehandle=logging.FileHandler(str(options.basedir)+'/logs/'+str(OBSID)+'/'+str(OBSID)+'_obslog.log')
	filehandle.setFormatter(formatter)
	obslog.addHandler(filehandle)
obslog.propagate = False

results=managing_caldatabase(str(options.msname),str(options.metafits),int(options.num_jobs),str(options.basedir),\
str(options.gaincal_modeldir),str(options.bandpass_modeldir),str(options.polcal_modeldir),str(options.localdatabase),float(inputs.gain_minsnr),str(inputs.ref_ant),\
freq_avg=float(options.freqavg),time_avg=float(options.timeavg),pol_skip_freq=float(polskip),verbose=eval(str(options.verbose)))
if results==None:
	print ('Errors in making caltables.\n')
else:
	msg,caltable_for_global_database,caltable_for_local_database=results

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
		if is_internet:
			obslog.info('Compressing caltables........\n')
			os.system('compress_caltables --caltables '+caltable_list+' --compressed_file '+final_caltable)
			if os.path.exists(final_caltable):
				attachments=[final_caltable]
				msg='Dear PAIRCARS developers,\n\nCaltables for : OBSID : '+str(OBSID)+' Coarse channels : '+str(coarse_chan_0)+'-'+\
										str(coarse_chan_1)+'\n\nBest,\nPAIRCARS developing team.\n'
				obslog.info('Sending caltables to database....\n')
				send_msg_code,send_msg=send_to_database('Caltables for : OBSID : '+str(OBSID)+' Coarse channels : '+str(coarse_chan_0)+'-'+\
										str(coarse_chan_1),msg,attachments=attachments)
				if send_msg_code==0:
					for i in caltable_for_global_database:
						if 'pcal' not in i and 'beam' not in i: 
							os.system('rm -rf '+i)
					os.system('rm -rf '+final_caltable)	
	os.system('touch '+inputs.basedir+'/.Finished_final_cal_'+os.path.basename(str(options.msname)))

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
		residual_frac=float(calparams.calc_calibration_params()[2])/1.5
	if inputs.maskfile!='':
		mask=inputs.maskfile
	elif inputs.maskstr!='':
		mask=inputs.maskstr
	else:
		mask=''

	start_sigma=np.load(str(options.basedir)+'/Ref_time_chan_sigma.npy',allow_pickle=True)[0]
	rms_list=np.load(str(options.basedir)+'/Ref_time_chan_sigma.npy',allow_pickle=True)[1]
	
	obslog.info('Flagging MWA coarse channel edges.\n')
	flag_MWA_coarse(str(options.msname),edgewidth=160,do_flag=True)
	obslog.info('Flagging QUACK times.\n')
	quacktime=flag_MWA_quack(str(options.msname),str(options.metafits))
	obslog.info('Flagged '+str(quacktime)+' s at beginning and end.\n')

	if is_internet:
		obslog.info('Start making images for global MWA solar database.\n')
		c=0
		while True:
			available_paircars_instance=get_available_paircars_instance(inputs.paircars_dir,inputs.job_id,inputs.instance)
			if c==0:
				obslog.info('Waiting for free P-AIRCARS instance......\n')
			c+=1
			if available_paircars_instance>1:
				break	
		
		a=final_imaging(str(options.msname),str(options.metafits),str(options.basedir),'Database_Imaging',casacals=casa_gaintable,calibratecals=calibrate_gaintable,\
				residual_frac=residual_frac,quality_factor=inputs.quality_factor,inputfile=inputfile,localdatabase=localdatabase,time_avg=float(options.timeavg),\
				freq_avg=float(options.freqavg),savedir=localdatabase+'/'+str(OBSID)+'/images',cutoutbox='3,3',want_automask=False,savemodel=False,saveres=False,\
				use_ankflag=inputs.use_ankflagger,do_pol=inputs.do_polcal,mask=mask,do_diffcal=True,freq_interval=1280,time_interval=30,freq_width=1280,\
				time_width=30,sigma=start_sigma,thresh=rms_list,use_wsclean=eval(str(options.use_wsclean)))

		AM=AccessMS(str(options.msname))
		ms_freq=AM.calc_meanfreq()/10**6	
		if inputs.send_notification==True:
			while True:
				database_images=glob.glob(localdatabase+'/'+str(OBSID)+'/images/All_final_images/'+str(OBSID)+'/*')
				if len(database_images)>0:
					database_images_freq=np.array([float(database_image.split('freq_')[-1].split('_image')[0]) for database_image in database_images])
					pos=np.argmin(np.abs(database_images_freq-ms_freq))
					database_image=np.array(database_images)[pos]
					break
			freqstr=database_image.split('freq_')[-1].split('_image')[0]
			datestrfile=database_image.split('time_')[-1].split('_freq')[0].split('_')
			datetimestr='/'.join(datestrfile[:3])+'/'+':'.join(datestrfile[3:])
			quickimage=get_quicklook_image(database_image,'sample_image_freq_'+freqstr+'_time_'+str('_'.join(datestrfile))+'.png',freqstr,datetimestr,field_of_view=2)
			msg_str='Dear PAIRCARS User,\n\nDatabase imaging for ms : '+str(os.path.basename(options.msname))+' is done\n\nBest regards,\nPAIRCARS developing team'
			msg_subject='Notification from PAIRCARS : Database imaging : OBSID = '+str(OBSID)
			send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
			os.system('rm -rf '+quickimage)
		os.system('rm -rf '+localdatabase+'/'+str(OBSID)+'/images/splited_ms/'+os.path.basename(options.msname).split('.ms')[0])

	try:
		cutoutbox=inputs.cutoutbox
	except:
		cutoutbox='3,3'

	c=0
	while True:
		available_paircars_instance=get_available_paircars_instance(inputs.paircars_dir,inputs.job_id,inputs.instance)
		if c==0:
			obslog.info('Waiting for free cores......\n')
		time.sleep(5)
		c+=1
		if available_paircars_instance>1:
			break

	obslog.info('Start making final images.\n')

	obslog.info('Waiting for finishing final imaging.......\n')
	a=final_imaging(str(options.msname),str(options.metafits),str(options.basedir),'Final_Imaging',casacals=casa_gaintable,calibratecals=calibrate_gaintable,\
			residual_frac=residual_frac,quality_factor=inputs.quality_factor,inputfile=inputfile,localdatabase=localdatabase,savedir=inputs.savedir,cutoutbox=inputs.cutoutbox,\
			want_automask=inputs.want_auto_masking,savemodel=inputs.savemodel,saveres=inputs.saveresidual,use_ankflag=inputs.use_ankflagger,do_pol=inputs.do_polcal,mask=mask,\
			freq_interval=inputs.image_delta_freq,time_interval=inputs.image_delta_time,freq_width=inputs.image_freq,time_width=inputs.image_time,sigma=start_sigma,thresh=rms_list,\
			use_wsclean=eval(str(options.use_wsclean)),time_avg=float(options.timeavg),freq_avg=float(options.freqavg))
	obslog.info('\n###########################\nFinal imaging finished.\n###########################\n')
	os.system('touch '+str(options.basedir)+'/.Imaging_done_'+str(os.path.basename(options.msname))+'_success')
#except Exception as e:
#		print ('Error occured : '+ str(e)+'\n')
#		os.system('touch '+str(options.basedir)+'/.Imaging_done_'+str(os.path.basename(options.msname))+'_error')
