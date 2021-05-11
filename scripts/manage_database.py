'''
Code is written by Devojyoti Kansabanik , 2 May, 2021
'''
from casatools import *
from casatasks import *
import os,sys,logging,numpy as np,copy,glob,psutil,time,multiprocessing as mp
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
	try:
		os.system('echo "import os\nprint (\'CHECK MPI\')" > test_mpi.py')
		os.system('mpirun -np 1 -x OMP_NUM_THREADS=3 python3 test_mpi.py')
		os.system('rm -rf test_mpi.py') 
		return 0
	except:	
		return 1

def casa_instance_runner(cmd,screen_name,finished_touch_file,num_thread):
	'''
	Function to run a casa instance
	Parameters:
	cmd = Command to run
	screen_name = Name of the screen
	'''
	mpi_check=MPI_check()
	print (mpi_check)
	if mpi_check==0:
		cmd='mpirun -np 1 -x OMP_NUM_THREADS='+str(num_thread)+' '+cmd
	cmd+=';wait; if ! ls '+finished_touch_file+'_* ; then  touch '+finished_touch_file+'_error ;  fi; rm -rf '+inputs.basedir+'/'+screen_name+'.batch'
	os.system('echo "'+cmd+'" > '+inputs.basedir+'/'+screen_name+'.batch')
	screen_cmd='sh '+inputs.basedir+'/'+screen_name+'.batch'
	os.system('screen -S '+screen_name+' -X quit')	
	time.sleep(0.5)
	os.system('screen -mdS '+screen_name)
	time.sleep(0.5)
	mainlog.info('########################\n')
	mainlog.info('Made Screen : '+screen_name+'\n')
	mainlog.info('Command : '+cmd+'\n')
	os.system('screen -S '+screen_name+' -X stuff \"'+screen_cmd+'\n"')	
	return screen_name


def managing_caldatabase(msname,metafits,total_spawned_jobs,basedir,gaincal_modeldir,bandpass_modeldir,polcal_modeldir,localdatabase,gain_minsnr,\
						ref_ant,freq_avg=160.0,pol_skip_freq=1.28):
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
	gaincal_modeldir=gaincal_modeldir+'/'+basemsdir
	bandpass_modeldir=bandpass_modeldir+'/'+basemsdir
	polcal_modeldir=polcal_modeldir+'/'+basemsdir

	while True:  # Waiting for all jobs to finish
		touch_files=len(glob.glob(basedir+'/.Finished*cal*'+str(OBSID)+'*'+basemsdir+'*'))
		if touch_files==total_spawned_jobs:
			obslog.info('All calibration jobs for ms : '+msname+' is completed.\n')
			break	
		else:
			time.sleep(2.0)
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

	# Final gain calibrations
	#########################
	if gaincal_modeldir=='' or os.path.isdir(gaincal_modeldir)==False or len(glob.glob(gaincal_modeldir+'/*.model'))==0:
		obslog.info('No models available.\n')  # Checking if any gaincal model is present or not. If not exit at this stage
		return 1,caltable_for_global_database,caltable_for_local_database
	else:
		AM=AccessMS(msname)
		freqs=AM.get_freqs()/10**6
		nchan_avg=int(freq_avg/AM.calc_freqres())
		coarse_chan_0=freq_to_MWA_coarse(freqs[0])
		coarse_chan_1=freq_to_MWA_coarse(freqs[-1])
		caltable_name_prefix=str(OBSID)+'_'+str(coarse_chan_0)+'_'+str(coarse_chan_1) # Caltable name prefix; OBSID_startcoarsechan_endcoarsechan format
		caltable_name=basedir+'/'+caltable_name_prefix+'.gcal' # Final gaincal table name
		model_list=glob.glob(gaincal_modeldir+'/*.model') # Gaincal model list
		if os.path.isdir(caltable_name):
			os.system('rm -rf '+caltable_name)
		obslog.info('#######################################\n')
		obslog.info('Making final gaincal table for ms : '+msname+'\n')
		timerange_list=[]
		ref_timestamp=''
		ref_model=''
		for i in range(len(model_list)):
			modelname=model_list[i]
			modelbasename=os.path.basename(modelname)
			timestamp='/'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[:3])+'/'+':'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[3:])
			timerange_list.append(timestamp)
			if 'ref' in modelname:
				ref_model=modelname
				ref_timestamp+=timestamp+','
		ref_timestamp=ref_timestamp[:-1] # Reference time
		timerange=','.join(timerange_list)
		f=imhead(imagename=ref_model,mode='list')['crval4']/10**6
		df=(imhead(imagename=ref_model,mode='list')['cdelt4']/10**6)/2.0
		spw=str(f-df)+'~'+str(f+df)+'MHz' # Spectral window range for reference model
		obslog.info('Spliting ms for gain calibration.....\n')
		if os.path.isdir(basedir+'/'+os.path.basename(msname)+'.gcalms'):
			os.system('rm -rf '+basedir+'/'+os.path.basename(msname)+'.gcalms '+basedir+'/'+os.path.basename(msname)+'.gcalms.flagversions')
		if os.path.isdir(basedir+'/'+os.path.basename(msname)+'.gcalms.ref'):
			os.system('rm -rf '+basedir+'/'+os.path.basename(msname)+'.gcalms.ref '+basedir+'/'+os.path.basename(msname)+'.gcalms.ref.flagversions')
		obslog.info('split(vis=\''+msname+'\',outputvis=\''+basedir+'/'+os.path.basename(msname)+'.gcalms\',datacolumn=\'DATA\',spw=\'0:'+spw+'\',timerange=\''+timerange+'\')\n')
		split(vis=msname,outputvis=basedir+'/'+os.path.basename(msname)+'.gcalms',datacolumn='DATA',spw='0:'+spw,timerange=timerange) # Reference channel ms
		obslog.info('split(vis=\''+msname+'\',outputvis=\''+basedir+'/'+os.path.basename(msname)+'.gcalms.ref\',datacolumn=\'DATA\',spw=\'0:'+spw+'\',timerange=\''+ref_timestamp+'\')\n')
		split(vis=msname,outputvis=basedir+'/'+os.path.basename(msname)+'.gcalms.ref',datacolumn='DATA',spw='0:'+spw,timerange=ref_timestamp) # Only reference channel and time ms
		AMgcal=AccessMS(basedir+'/'+os.path.basename(msname)+'.gcalms.ref')
		ntimes=AMgcal.get_num_timestamps()
		nchan=AMgcal.get_num_channels()
		obslog.info('delmod(vis=\''+basedir+'/'+os.path.basename(msname)+'.gcalms.ref\',scr=True)\n') # Deleting any previous models in ms and importing refernece time chan model
		delmod(vis=basedir+'/'+os.path.basename(msname)+'.gcalms.ref',scr=True)
		obslog.info('ft(vis=\''+basedir+'/'+os.path.basename(msname)+'.gcalms.ref\',model=\''+ref_model+'\',spw=\'0:'+str(spw)+'\',usescratch=True)\n')
		ft(vis=basedir+'/'+os.path.basename(msname)+'.gcalms.ref',model=ref_model,spw='0:'+str(spw),usescratch=True)
		for i in range(len(model_list)): # Importing model for all times dequentially for reference chan and calibrating 
			modelname=model_list[i]
			modelbasename=os.path.basename(modelname)
			timestamp='/'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[:3])+'/'+':'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[3:])
			timerange_list.append(timestamp)
			obslog.info('delmod(vis=\''+basedir+'/'+os.path.basename(msname)+'.gcalms\',scr=True)\n')
			delmod(vis=basedir+'/'+os.path.basename(msname)+'.gcalms',scr=True)
			obslog.info('ft(vis=\''+basedir+'/'+os.path.basename(msname)+'.gcalms\',model=\''+modelname+'\',spw=\'0:'+str(spw)+'\',usescratch=True)\n')
			ft(vis=basedir+'/'+os.path.basename(msname)+'.gcalms',model=modelname,spw='0:'+str(spw),usescratch=True)
			IB=ImageBasic(msname)
			uvrange_to_cal=IB.calc_calib_uvrange(4)[0]
			obslog.info('gaincal(vis=\''+basedir+'/'+os.path.basename(msname)+'.gcalms\',caltable=\''+caltable_name+'\',spw=\'0:'+str(spw)+'\',timerange=\''+timestamp+\
						'\',append=True,uvrange=\''+uvrange_to_cal+'\',solnorm=True,rmsthresh=[10,8,6],refant=\''+str(ref_ant)+'\',minsnr='+str(gain_minsnr)+')\n')
			gaincal(vis=basedir+'/'+os.path.basename(msname)+'.gcalms',caltable=caltable_name,spw='0:'+str(spw),timerange=timestamp,append=True,\
					uvrange=uvrange_to_cal,solnorm=True,rmsthresh=[10,8,6],refant=str(ref_ant),minsnr=gain_minsnr) # Appedning solutions into a same gain table
		IB1=ImageBasic(msname)
		calib_uvrange_min=IB1.calc_calib_uvrange(4)[1]
		calib_uvrange_max=IB1.calc_calib_uvrange(4)[2]
		obslog.info('cal.calibrate(msname=\''+basedir+'/'+os.path.basename(msname)+'.gcalms.ref\',caltable=\''+caltable_name+'.bin\',calmode=\'diag\',minuv='+\
					str(calib_uvrange_min)+',maxuv='+str(calib_uvrange_max)+',quiet=True,a=\'0.001,0.0001\'t='+\
					str(int(ntimes/len(model_list)))+',j=2,ch='+str(int(AMgcal.get_num_channels()))+')\n') # Making gaincal table for reference time and chan
		cal.calibrate(msname=basedir+'/'+os.path.basename(msname)+'.gcalms.ref',caltable=caltable_name+'.bin',calmode='diag',minuv=calib_uvrange_min,maxuv=calib_uvrange_max,\
			quiet=True,a='0.001,0.0001',t=int(ntimes/len(model_list)),j=2,ch=int(AMgcal.get_num_channels()))
		os.system('rm -rf '+basedir+'/'+os.path.basename(msname)+'.gcalms '+basedir+'/'+os.path.basename(msname)+'.gcalms.flagversions')
		os.system('rm -rf '+basedir+'/'+os.path.basename(msname)+'.gcalms.ref '+basedir+'/'+os.path.basename(msname)+'.gcalms.ref.flagversions')
		os.system('cp -r '+caltable_name+' '+localdatabase) # Copying to local database
		os.system('cp -r '+caltable_name+'.bin '+localdatabase) # Copying to local database for further copying to global database
		caltable_for_global_database.append(localdatabase+'/'+os.path.basename(caltable_name)+'.bin')
		caltable_for_local_database.append(localdatabase+'/'+os.path.basename(caltable_name))
		os.system('rm -rf '+caltable_name+'*')
		# Applying gain solutions and spliting for bandpass or polarisation calibration 
		###############################################################################
		AM=AccessMS(msname)
		timestamps=AM.get_timestamps()
		timeres=AM.calc_timeres()
		index=timestamps.index(timerange)
		if len(timestamps)<10 or timeres>2: # If time resolution is greater than 2 sec, only split the bandpass model timerange
			timerange=timestamps[index]
			obslog.info('applycal(vis=\''+msname+'\',gaintable=[\''+localdatabase+'/'+os.path.basename(caltable_name)+\
						'\'],timerange=\''+timerange+'\',applymode=\'calonly\',flagbackup=False)\n')
			applycal(vis=msname,gaintable=[localdatabase+'/'+os.path.basename(caltable_name)],timerange=timerange,applymode='calonly',flagbackup=False)
			if os.path.isdir(basedir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftime.ms'):
				os.system('rm -rf '+basedir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftime.ms '+basedir+'/'\
						+os.path.basename(msname).split('.ms')[0]+'_reftime.ms.flagversions')
			obslog.info('split(vis=\''+msname+'\',outputvis=\''+basedir+'/'+\
						os.path.basename(msname).split('.ms')[0]+'_reftime.ms\',timerange=\''+timerange+'\',datacolumn=\'corrected\')\n')
			split(vis=msname,outputvis=basedir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftime.ms',width=nchan_avg,timerange=timerange,datacolumn='corrected')
		elif timeres<2: # If time resolution is less than 2 sec, split at 2 sec around bandpass model time and use same model for calibration
			if (index-5)<0:
				timerange=','.join(timestamps[0:10])
			else:
				timerange=','.join(timestamps[index-5:index+5])
			obslog.info('applycal(vis=\''+msname+'\',gaintable=[\''+localdatabase+'/'+os.path.basename(caltable_name)+\
						'\'],timerange=\''+timerange+'\',applymode=\'calonly\',flagbackup=False)\n')
			applycal(vis=msname,gaintable=[localdatabase+'/'+os.path.basename(caltable_name)],timerange=timerange,applymode='calonly',flagbackup=False)
			if os.path.isdir(basedir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftime.ms'):
				os.system('rm -rf '+basedir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftime.ms '+basedir+'/'+\
					os.path.basename(msname).split('.ms')[0]+'_reftime.ms.flagversions')
			obslog.info('split(vis=\''+msname+'\',outputvis=\''+basedir+'/'\
						+os.path.basename(msname).split('.ms')[0]+'_reftime.ms\',timerange=\''+timerange+'\',datacolumn=\'corrected\')\n')
			split(vis=msname,outputvis=basedir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftime.ms',width=nchan_avg,timerange=timerange,datacolumn='corrected')
		bpmsname=basedir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftime.ms' # Reftime ms for bandpass calibration
	
		# Performing bandpass
		#####################
		done_bandpass=False
		if os.path.isdir(bandpass_modeldir)==False or len(glob.glob(bandpass_modeldir+'/*.model'))==0: # Final bandpass
			obslog.info('No bandpass models are available.\n') # If no bandpass table avaiable do not proceed
			done_bandpass=False
		else:
			obslog.info('#######################################\n')
			obslog.info('Making final bandpass table for ms : '+msname+'\n')
			bp_caltable_name=basedir+'/'+caltable_name_prefix+'.bcal'
			if os.path.isdir(bp_caltable_name):
				os.system('rm -rf '+bp_caltable_name)
			model_list=glob.glob(bandpass_modeldir+'/*.model')
			modelname=model_list[0]
			modelbasename=os.path.basename(modelname)
			timerange='/'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[:3])+'/'+':'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[3:])
			AM2=AccessMS(bpmsname)
			freqlist=(AM2.get_freqs()/10**6)
			obslog.info('delmod(vis=\''+bpmsname+'\',scr=True)\n')
			delmod(vis=bpmsname,scr=True)
			for i in range(len(model_list)): # Importing models to different channels
				modelname=model_list[i]
				f=imhead(imagename=modelname,mode='list')['crval4']/10**6
				df=(imhead(imagename=modelname,mode='list')['cdelt4']/10**6)/2.0
				spw='0:'+str(f-df)+'~'+str(f+df)+'MHz'
				obslog.info('ft(vis=\''+bpmsname+'\',model=\''+modelname+'\',spw=\''+spw+'\',usescratch=True)\n')
				ft(vis=bpmsname,model=modelname,spw=spw,usescratch=True)
			AMbp=AccessMS(bpmsname)
			model,nomodel=AMbp.get_model_nomodel_chan()
			unflag_chan=AMbp.get_unflag_chan()
			nomodel_str_list=[str(i) for i in nomodel]
			flagchans='0:'+';'.join(nomodel_str_list)
			obslog.info('flagdata(vis=\''+bpmsname+'\',mode=\'manual\',spw=\''+flagchans+'\')\n')
			flagdata(vis=bpmsname,mode='manual',spw=flagchans) # Flagging channels which do not have models, such that no bad solutions will be stored.
			obslog.info('bandpass(vis=\''+bpmsname+'\',caltable='+str(bp_caltable_name)+',solnorm=True,refant=\''\
					+str(ref_ant)+'\',minsnr='+str(gain_minsnr)+',uvrange=\''+uvrange_to_cal+'\',bandtype=\'B\',fillgaps=3)')
			bandpass(vis=bpmsname,caltable=bp_caltable_name,solnorm=True,refant=str(ref_ant),minsnr=gain_minsnr,uvrange=uvrange_to_cal,bandtype='B',fillgaps=3) # Performing bandpass
			if model==unflag_chan: # If models are available for all unflagged channels perform bandpass for global database
				obslog.info('Calibrating for global database......\n')
				obslog.info('cal.calibrate(msname=\''+bpmsname+'\',caltable=\''+bp_caltable_name+'.bin\',calmode=\'diag\',minuv='+str(calib_uvrange_min)+\
							',maxuv='+str(calib_uvrange_max)+',quiet=True,a=\'0.001,0.0001\',i=1000,j=2,ch=1)\n')
				cal.calibrate(msname=bpmsname,caltable=bp_caltable_name+'.bin',calmode='diag',minuv=calib_uvrange_min,maxuv=calib_uvrange_max,quiet=True,a='0.001,0.0001',i=1000,j=2,ch=1)
				os.system('cp -r '+bp_caltable_name+'.bin '+localdatabase)
				caltable_for_global_database.append(localdatabase+'/'+os.path.basename(bp_caltable_name)+'.bin')
			os.system('cp -r '+bp_caltable_name+' '+localdatabase)
			os.system('rm -rf '+bp_caltable_name+'*')
			caltable_for_local_database.append(localdatabase+'/'+os.path.basename(bp_caltable_name))
			obslog.info('flagmanager(vis=\''+bpmsname+'\',mode=\'restore\',versionname=\'flagdata_1\')\n')
			flagmanager(vis=bpmsname,mode='restore',versionname='flagdata_1') # Restoring flags
			obslog.info('flagmanager(vis=\''+bpmsname+'\',mode=\'delete\',versionname=\'flagdata_1\')\n')
			flagmanager(vis=bpmsname,mode='delete',versionname='flagdata_1')
			done_bandpass=True
		# Performing polarisation calibration
		#####################################
		if os.path.isdir(polcal_modeldir)==False or len(glob.glob(polcal_modeldir+'/*.model'))==0: # Final polarisation calibration
			obslog.info('No polarisation caltables are available.\n') # If no polarisation model available, do not continue
		else:
			obslog.info('#######################################\n')
			obslog.info('Making final polcal table for ms : '+msname+'\n')
			if done_bandpass==True:
				bptables=[localdatabase+'/'+os.path.basename(bp_caltable_name)]
				obslog.info('applycal(vis=\''+bpmsname+'\',gaintable='+str(bptables)+',applymode=\'calflag\',flagbackup=False)\n')			
				applycal(vis=bpmsname,gaintable=bptables,applymode='calflag',flagbackup=False)
				if os.path.isdir(bpmsname+'.polleakcal')==True:
					os.system('rm -rf '+bpmsname+'.polleakcal '+bpmsname+'.polleakcal.flagversions')
				obslog.info('Spliting bandpass calibrated ms.\n')
				obslog.info('split(vis=\''+bpmsname+'\',outputvis=\''+bpmsname+'.polleakcal\',datacolumn=\'corrected\')\n')
				split(vis=bpmsname,outputvis=bpmsname+'.polleakcal',datacolumn='corrected')
			else:
				if os.path.isdir(bpmsname+'.polleakcal')==True:
					os.system('rm -rf '+bpmsname+'.polleakcal '+bpmsname+'.polleakcal.flagversions')
				obslog.info('Spliting gain calibrated ms.\n')
				obslog.info('split(vis=\''+bpmsname+'\',outputvis=\''+bpmsname+'.polleakcal\',datacolumn=\'data\')\n')
				split(vis=bpmsname,outputvis=bpmsname+'.polleakcal',datacolumn='data')
			leakcor_gain_models=glob.glob(polcal_modeldir+'/*_leakage.model') # Leakage calibrated source models
			obslog.info('delmod(vis=\''+bpmsname+'.polleakcal\',scr=True)\n')
			delmod(vis=bpmsname+'.polleakcal',scr=True)
			for i in range(len(leakcor_gain_models)):
				modelname=leakcor_gain_models[i]
				f=imhead(imagename=modelname,mode='list')['crval4']/10**6
				df=(imhead(imagename=modelname,mode='list')['cdelt4']/10**6)/2.0
				spw='0:'+str(f-df)+'~'+str(f+df)+'MHz'
				obslog.info('ft(vis=\''+bpmsname+'.polleakcal\',model=\''+modelname+'\',spw=\''+spw+'\',usescratch=True)\n')
				ft(vis=bpmsname+'.polleakcal',model=modelname,spw=spw,usescratch=True)
			AM1=AccessMS(bpmsname+'.polleakcal')
			unflagged_chan=AM1.get_unflag_chan()
			model_chan,nomodel_chan=AM1.get_model_nomodel_chan()
			nomodel_chan_str_list=[str(i) for i in nomodel_chan]
			flagchans='0:'+';'.join(nomodel_chan_str_list)
			obslog.info('flagdata(vis=\''+bpmsname+'.polleakcal\',mode=\'manual\',spw=\''+flagchans+'\')\n')
			flagdata(vis=bpmsname+'.polleakcal',mode='manual',spw=flagchans)
			freqs=AM1.get_freqs()/10**6
			polleak_caltable=basedir+'/'+caltable_name_prefix+'.lcal' # Leakage corrected differential gain correction gaintable
			beamfile=basedir+'/'+caltable_name_prefix+'.beam' # Ideal beam jones table
			pol_caltable=basedir+'/'+caltable_name_prefix+'.pcal.bin' # Polarisation calibration table
			obslog.info('bandpass(vis=\''+bpmsname+'.polleakcal\',caltable=\''+polleak_caltable+'\',solnorm=True,refant=\''+str(ref_ant)+\
						'\',minsnr='+str(gain_minsnr)+',uvrange=\''+uvrange_to_cal+'\',bandtype=\'B\',fillgaps=3)\n') # Performing leakage corrected gain calibration
			bandpass(vis=bpmsname+'.polleakcal',caltable=polleak_caltable,solnorm=True,refant=str(ref_ant),minsnr=gain_minsnr,uvrange=uvrange_to_cal,bandtype='B',fillgaps=3)
			obslog.info('Calibrating for global database......\n')
			obslog.info('cal.calibrate(msname=\''+bpmsname+'.polleakcal\',caltable=\''+polleak_caltable+'.bin\',calmode=\'diag\',minuv='+str(calib_uvrange_min)\
						+',maxuv='+str(calib_uvrange_max)+',quiet=True,a=\'0.001,0.0001\',i=1000,j=2,ch=1)\n')
			cal.calibrate(msname=bpmsname+'.polleakcal',caltable=polleak_caltable+'.bin',calmode='diag',minuv=calib_uvrange_min,\
							maxuv=calib_uvrange_max,quiet=True,a='0.001,0.0001',i=1000,j=2,ch=1)
			os.system('cp -r '+polleak_caltable+'.bin '+localdatabase)
			caltable_for_global_database.append(localdatabase+'/'+os.path.basename(polleak_caltable)+'.bin')
			os.system('cp -r '+polleak_caltable+' '+localdatabase)
			os.system('rm -rf '+polleak_caltable+'*')
			caltable_for_local_database.append(localdatabase+'/'+os.path.basename(polleak_caltable))
			obslog.info('flagmanager(vis=\''+bpmsname+'.polleakcal\',mode=\'restore\',versionname=\'flagdata_1\')\n')
			flagmanager(vis=bpmsname+'.polleakcal',mode='restore',versionname='flagdata_1')
			obslog.info('flagmanager(vis=\''+bpmsname+'.polleakcal\',mode=\'delete\',versionname=\'flagdata_1\')\n')
			flagmanager(vis=bpmsname+'.polleakcal',mode='delete',versionname='flagdata_1')
			obslog.info('applycal(vis=\''+bpmsname+'.polleakcal\',gaintable='+str(localdatabase+'/'+os.path.basename(polleak_caltable))+\
						',applymode=\'calflag\',flagbackup=False)\n')		
			applycal(vis=bpmsname+'.polleakcal',gaintable=[localdatabase+'/'+os.path.basename(polleak_caltable)],applymode='calflag',flagbackup=False) # Applying leakage corrected gains
			polcal_ms=bpmsname+'.polcal'
			if os.path.isdir(polcal_ms)==True:
				os.system('rm -rf '+polcal_ms+' '+polcal_ms+'.flagversions')
			obslog.info('split(vis=\''+bpmsname+'.polleakcal\',outputvis=\''+polcal_ms+'\',datacolumn=\'corrected\')\n')
			split(vis=bpmsname+'.polleakcal',outputvis=polcal_ms,datacolumn='corrected')
			polcal_models=glob.glob(polcal_modeldir+'/*.model')
			for i in leakcor_gain_models:
				polcal_models.remove(i)
			os.system('rm -rf '+localdatabase+'/'+bp_caltable_name+'_temp*cal')
			obslog.info('delmod(vis=\''+polcal_ms+'\',scr=True)\n')
			delmod(vis=polcal_ms,scr=True)
			for i in range(len(polcal_models)):
				modelname=polcal_models[i]
				f=imhead(imagename=modelname,mode='list')['crval4']/10**6
				df=(imhead(imagename=modelname,mode='list')['cdelt4']/10**6)/2.0
				spw='0:'+str(f-df)+'~'+str(f+df)+'MHz'
				obslog.info('ft(vis=\''+polcal_ms+'\',model=\''+modelname+'\',spw=\''+spw+'\',usescratch=True)\n')
				ft(vis=polcal_ms,model=modelname,spw=spw,usescratch=True)
			obslog.info('PolSelfcal(\''+polcal_ms+'\',\''+metafits+'\','+str(32*60)+',verbose=False,interactive=False)\n')
			PSC=PolSelfcal(polcal_ms,metafits,32*60,verbose=False,interactive=False) # Performing ideal beam correction at phase center for every coarse channels
			obslog.info('PSC.correct_visibility_single_beam_jones(modify_datacolumn=True)\n')
			PSC.correct_visibility_single_beam_jones(modify_datacolumn=True,skip_freq=float(pol_skip_freq),save_beamfile=beamfile)
			obslog.info('cal.calibrate(msname=\''+polcal_ms+'\',caltable=\''+pol_caltable+'\',minuv='+str(calib_uvrange_min)\
						+',quiet=True,maxuv='+str(calib_uvrange_max)+',j=2)\n')
			cal.calibrate(msname=polcal_ms,caltable=pol_caltable,minuv=calib_uvrange_min,quiet=True,maxuv=calib_uvrange_max,j=2)
			os.system('cp -r '+pol_caltable+' '+localdatabase)
			os.system('cp -r '+beamfile+' '+localdatabase)	
			os.system('rm -rf '+pol_caltable+'* '+beamfile+'*')		
			caltable_for_global_database.append(localdatabase+'/'+os.path.basename(pol_caltable))	
			caltable_for_local_database.append(localdatabase+'/'+os.path.basename(beamfile))			
			caltable_for_local_database.append(localdatabase+'/'+os.path.basename(pol_caltable))			
	os.system('rm -rf '+basedir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftime* '+basedir+'/'+os.path.basename(msname).split('.ms')[0]+'.temp*')
	obslog.info('All final calibrations are finished.\n')
	return 0,caltable_for_global_database,caltable_for_local_database


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
		splited_ms=glob.glob(workdir+'/*.ms')
		splited_ms_copy=copy.deepcopy(splited_ms)
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
		while len(splited_ms)!=0:
			# Estimating total casa instances
			#################################
			total_available_cpu=psutil.cpu_count()-(psutil.cpu_count()*psutil.cpu_percent()/100.0)
			available_cpu_for_paircars=int(total_available_cpu*inputs.cpu_frac)
			total_cpu_frac=int(psutil.cpu_count()*inputs.cpu_frac)
			if available_cpu_for_paircars>total_cpu_frac:
				available_cpu_for_paircars=total_cpu_frac
			casa_instance=int(available_cpu_for_paircars/2)
			open_casa_instance=0
			for i in range(casa_instance):
				ms=splited_ms[0]
				if ms[-1]=='/':
					ms=ms[:-1]
				workdir=os.path.abspath(ms)+'/'+os.path.basename(ms)
				cmd='final_imaging --msname '+ms+' --metafits '+metafits+' --basedir '+basedir+' --workdir '+workdir+' --savedir '+savedir+' --savemodel '+str(savemodel)\
					+' --saveres '+str(saveres)+' --stokes '+stokes+' --cutoutbox '+cutoutbox+' --sigma '+str(sigma)+' --threshold '+str(threshold)\
					+' --want_automask '+str(want_automask)+'--maskfile '+str(mask)+' --quality_factor '+str(quality_factor)+' --inputfile '+str(inputfile)\
						+' --use_ankflag '+str(use_ankflag)+' --residual_frac '+str(residual_frac)
				screen_name='Final_imaging_'+os.path.basename(ms).split('.ms')[0]
				finished_touch_file=basedir+'/.Finished_final_imaging_'+os.path.basename(msname)+'_'
				casa_instance_runner(cmd,screen_name,finished_touch_file,2)
				splited_ms.remove(ms)
				touch_count+=1
				open_casa_instance+=1
			while True:
				obslog.info('Waiting for final imaging spawned to be finished......\n')
				time.sleep(2.0)
				if len(glob.glob(basedir+'/.Finished_final_imaging_'+os.path.basename(msname)+'_*'))>=(touch_count-open_casa_instance):
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
	parser.add_option('--inputfile',dest='inputfile',default='',help='Path of the P-AIRCARS input file',metavar="File path")
	(options, args) = parser.parse_args()
	
if os.path.isdir(str(options.msname))==False:
	print ('Measurement set is not present.\n')
	os._exit(0)
	
if options.basedir==None:
	print ('No base directory path is given. Please give base directory path to continue.\n')
	os._exit(0)
elif os.path.isdir(str(options.basedir))==False:
	os.makedirs(str(options.basedir))
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

msg,caltable_for_global_database,caltable_for_local_database=managing_caldatabase(str(options.msname),\
	str(options.metafits),int(options.num_jobs),str(options.basedir),str(options.gaincal_modeldir),str(options.bandpass_modeldir),\
	str(options.polcal_modeldir),str(options.localdatabase),float(inputs.gain_minsnr),str(inputs.ref_ant),freq_avg=float(options.freqavg),pol_skip_freq=float(polskip))
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

result=final_imaging_for_database(str(options.msname),str(options.metafits),basedir,casacals=casa_gaintable,calibratecals=calibrate_gaintable,residual_frac=residual_frac,\
		quality_factor=inputs.quality_factor,inputfile=inputfile,localdatabase=localdatabase,savedir=localdatabase+'/'+str(OBSID)+'/images',cutoutbox='3,3',\
		want_automask=False,savemodel=False,saveres=False,use_ankflag=False,do_pol=inputs.do_polcal,mask=mask,\
		freq_interval=1280,time_interval=30,freq_width=1280,time_width=30,sigma=start_sigma,threshold=rms_thresh)

'''
result=final_imaging_for_database(str(options.msname),str(options.metafits),basedir,casacals=casa_gaintable,calibratecals=calibrate_gaintable,residual_frac=residual_frac,\
		quality_factor=inputs.quality_factor,inputfile=inputfile,localdatabase=localdatabase,savedir=inputs.savedir,cutoutbox=inputs.cutoutbox,want_automask=inputs.want_auto_masking,\
		savemodel=inputs.savemodel,saveres=inputs.saveresidual,use_ankflag=inputs.use_ankflagger,do_pol=inputs.do_polcal,mask=mask,freq_interval=inputs.image_delta_freq,\
		time_interval=image_delta_time,freq_width=inputs.image_freq,time_width=inputs.image_time,sigma=start_sigma,threshold=rms_thresh)

'''













