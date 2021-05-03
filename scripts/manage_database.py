'''
Code is written by Devojyoti Kansabanik , 2 May, 2021
'''
from casatools import *
from casatasks import *
import os,sys,logging,numpy as np,copy,glob,psutil,time
from paircars.basic_func import *
from paircars.access_ms import *
from paircars.decor import *
from paircars.flagger import *
from astropy.io import fits
from CALIBRATE.access_calibrate import *

def managing_caldatabase(msname,OBSID,total_spawned_jobs,basedir,gaincal_modeldir,bandpass_modeldir,polcal_caldir,localdatabase,gain_minsnr,ref_ant,ncoarse_pol):
	'''
	Function to manager calibration database
	msname : Averaged msname
	gaincal_modeldir = Model directory for gaincal
	bandpass_modeldir = Model directory for bandpass
	polcal_caldir = Polarisation caltable directory
	localdatabase = Local database directory
	'''
	formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
	obslog = logging.getLogger(str(OBSID)+'_log')
	obslog.setLevel(logging.DEBUG)
	console=logging.StreamHandler(sys.stdout)
	console.setFormatter(formatter)
	obslog.addHandler(console)
	filehandle=logging.FileHandler(basedir+'/'+str(OBSID)+'_obslog.log')
	filehandle.setFormatter(formatter)
	obslog.addHandler(filehandle)
	obslog.propagate = False
	while True:
		obslog.info('Waiting for finishing '+str(total_spawned_jobs)+' calibration jobs for ms : '+msname+'\n')
		touch_files=len(glob.glob(basedir+'/.Finished*cal*'+str(OBSID)+'*'))
		if touch_files==total_spawned_jobs:
			obslog.info('All calibration jobs for ms : '+msname+' is completed.\n')
			break	
		else:
			time.sleep(2.0)
	if localdatabase=='' or os.path.isdir(localdatabase)==False:
		obslog.error('Local data base not found. Making local database at basedir.\n')
		localdatabase=basedir+'/localdatabase/'+str(OBSID)
		if os.path.isdir(localdatabase)==False:
			os.makedirs(localdatabase)
	else:
		obslog.info('Local data base is at : '+localdatabase+'\n')
		localdatabase=localdatabase+'/'+str(OBSID)
		if os.path.isdir(localdatabase)==False:
			os.makedirs(localdatabase)
	if gaincal_modeldir=='' or os.path.isdir(gaincal_modeldir)==False or len(glob.glob(gaincal_modeldir+'/*.model'))==0:
		obslog.info('No models available.\n')
		return 1
	else:
		AM=AccessMS(msname)
		freqs=AM.get_freqs()/10**6
		coarse_chan_0=freq_to_MWA_coarse(freqs[0])
		coarse_chan_1=freq_to_MWA_coarse(freqs[-1])
		caltable_name=str(OBSID)+'_'+str(coarse_chan_0)+'_'+str(coarse_chan_1)+'.gcal'
		model_list=glob.glob(gaincal_modeldir+'/*.model')
		if os.path.isdir(caltable_name):
			os.system('rm -rf '+caltable_name)
		obslog.info('#######################################\n')
		obslog.info('Making final gaincal table for ms : '+msname+'\n')
		for i in range(len(model_list)):
			modelname=model_list[i]
			modelbasename=os.path.basename(modelname)
			timestamp='/'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[:3])+'/'+':'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[3:])
			f=imhead(imagename=modelname,mode='list')['crval4']/10**6
			df=(imhead(imagename=modelname,mode='list')['cdelt4']/10**6)/2.0
			spw=str(f-df)+'~'+str(f+df)+'MHz'
			obslog.info('delmod(vis=\''+msname+'\',scr=True)\n')
			delmod(vis=msname,scr=True)
			obslog.info('ft(vis=\''+msname+'\',model=\''+modelname+'\',spw=\'0:'+str(spw)+'\',usescratch=True)\n')
			ft(vis=msname,model=modelname,spw='0:'+str(spw),usescratch=True)
			IB=ImageBasic(msname)
			uvrange_to_cal=IB.calc_calib_uvrange(4)[0]
			obslog.info('gaincal(vis=\''+msname+'\',caltable=\''+caltable_name+'\',spw=\'0:'+str(spw)+'\'timerange=\''+timestamp+'\',append=True,uvrange=\''+uvrange_to_cal\
				+'\',solnorm=True,rmsthresh=[10,8,6],refant=\''+str(ref_ant)+'\',minsnr='+str(gain_minsnr)+')\n')
			gaincal(vis=msname,caltable=caltable_name,spw='0:'+str(spw),timerange=timestamp,append=True,uvrange=uvrange_to_cal,solnorm=True,rmsthresh=[10,8,6],\
					refant=str(ref_ant),minsnr=gain_minsnr)
		os.system('cp -r '+caltable_name+' '+localdatabase)
		os.system('rm -rf '+caltable_name)
		if os.path.isdir(bandpass_modeldir)==False or len(glob.glob(bandpass_modeldir+'/*.model'))==0: # Backup bandpass
			obslog.info('No bandpass models are available.\n')
		else:
			obslog.info('#######################################\n')
			obslog.info('Making final bandpass table for ms : '+msname+'\n')
			bp_caltable_name=str(OBSID)+'_'+str(coarse_chan_0)+'_'+str(coarse_chan_1)+'.bcal'
			if os.path.isdir(bp_caltable_name):
				os.system('rm -rf '+bp_caltable_name)
			model_list=glob.glob(bandpass_modeldir+'/*.model')
			modelname=model_list[0]
			modelbasename=os.path.basename(modelname)
			timerange='/'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[:3])+'/'+':'.join(modelbasename.split('time_')[1].split('_freq')[0].split('_')[3:])
			timestamps=AM.get_timestamps()
			timeres=AM.calc_timeres()
			index=timestamps.index(timerange)
			if len(timestamps)<10 or timeres>2:
				timerange=timestamps[index]
				obslog.info('applycal(vis=\''+msname+'\',gaintable=[\''+localdatabase+'/'+caltable_name+'\'],timerange=\''+timerange+'\',applymode=\'calonly\',flagbackup=False)\n')
				applycal(vis=msname,gaintable=[localdatabase+'/'+caltable_name],timerange=timerange,applymode='calonly',flagbackup=False)
				if os.path.isdir(basedir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftime.ms'):
					os.system('rm -rf '+basedir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftime.ms')
				obslog.info('split(vis=\''+msname+'\',outputvis=\''+basedir+'/'+\
							os.path.basename(msname).split('.ms')[0]+'_reftime.ms\',timerange=\''+timerange+'\',datacolumn=\'corrected\')\n')
				split(vis=msname,outputvis=basedir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftime.ms',timerange=timerange,datacolumn='corrected')
			elif timeres<2:
				if (index-5)<0:
					timerange=','.join(timestamps[0:10])
				else:
					timerange=','.join(timestamps[index-5:index+5])
				obslog.info('applycal(vis=\''+msname+'\',gaintable=[\''+localdatabase+'/'+caltable_name+'\'],timerange=\''+timerange+'\',applymode=\'calonly\',flagbackup=False)\n')
				applycal(vis=msname,gaintable=[localdatabase+'/'+caltable_name],timerange=timerange,applymode='calonly',flagbackup=False)
				if os.path.isdir(basedir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftime.ms'):
					os.system('rm -rf '+basedir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftime.ms')
				obslog.info('split(vis=\''+msname+'\',outputvis=\''+basedir+'/'\
							+os.path.basename(msname).split('.ms')[0]+'_reftime.ms\',timerange=\''+timerange+'\',datacolumn=\'corrected\')\n')
				split(vis=msname,outputvis=basedir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftime.ms',timerange=timerange,datacolumn='corrected')
			bpmsname=basedir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftime.ms'
			AM2=AccessMS(bpmsname)
			freqlist=(AM2.get_freqs()/10**6)
			obslog.info('delmod(vis=\''+bpmsname+'\',scr=True)\n')
			delmod(vis=bpmsname,scr=True)
			for i in range(len(model_list)):
				modelname=model_list[i]
				f=imhead(imagename=modelname,mode='list')['crval4']/10**6
				df=(imhead(imagename=modelname,mode='list')['cdelt4']/10**6)/2.0
				spw='0:'+str(f-df)+'~'+str(f+df)+'MHz'
				obslog.info('ft(vis=\''+bpmsname+'\',model=\''+modelname+'\',spw=\''+spw+'\',usescratch=True)\n')
				ft(vis=bpmsname,model=modelname,spw=spw,usescratch=True)
			AM1=AccessMS(bpmsname)
			unflagged_chan=AM.get_unflag_chan()
			model_chan,nomdel_chan=AM.get_model_nomodel_chan()
			obslog.info('bandpass(vis=\''+bpmsname+'\',caltable=\''+bp_caltable_name+'\',spw=\''+spw+'\',solnorm=True,refant=\''+str(ref_ant)+'\',minsnr='+str(gain_minsnr)+')\n')
			bandpass(vis=bpmsname,caltable=bp_caltable_name,spw=spw,solnorm=True,refant=str(ref_ant),minsnr=gain_minsnr)
			os.system('cp -r '+bp_caltable_name+' '+localdatabase)
			#if model_chan==unflagged_chan:
				#TODO: Updating global database
			os.system('rm -rf '+bp_caltable_name)
		if os.path.isdir(polcal_caldir)==False or len(glob.glob(polcal_caldir+'/*.bin'))==0: # Backup polcal
			obslog.info('No polarisation caltables are available.\n')
		else:
			obslog.info('#######################################\n')
			obslog.info('Making final polcal table for ms : '+msname+'\n')
			polcaltable_list=glob.glob(polcal_caldir+'/*.bin')
			for polcal in polcaltable_list:
				freq=float(polcal.split('.bin')[0].split('freq_')[-1].split('_')[0]) # In MHz
				coarse_chan_0=freq_to_MWA_coarse(freq)
				coarse_chan_1+=ncoarse_pol
				polcaltable_name=str(OBSID)+'_'+str(coarse_chan_0)+'_'+str(coarse_chan_1)+'.bin'
				os.system('cp -r '+polcal+' '+localdatabase+'/'+polcaltable_name)	
	os.system('rm -rf '+basedir+'/time*reftime.ms')
	return 0


#def run_final_imaging(msname,OBSID,channel_grid,time_grid):
	




from optparse import OptionParser
if __name__=='__main__':
	usage= ' PAIRCARS database manager'
	parser = OptionParser(usage=usage)
	parser.add_option('--msname',dest="msname",default=None,help="Name of the measurement set",metavar="Measurement set")
	parser.add_option('--OBSID',dest="OBSID",default=None,help="Observation ID of MWA observation",metavar="Integer")
	parser.add_option('--num_jobs',dest="num_jobs",default=0,help="Toatal number of calibration jobs spawned for this measurement set",metavar="Integer")
	parser.add_option('--basedir',dest="basedir",default=None,help="Name of base directory for a given day",metavar="Directory path")
	parser.add_option('--gaincal_modeldir',dest="gaincal_modeldir",default=None,help="Name of gaincal model directory",metavar="Directory path")
	parser.add_option('--bandpass_modeldir',dest="bandpass_modeldir",default=None,help="Name of bandpass model directory",metavar="Directory path")
	parser.add_option('--polcal_caldir',dest="polcal_caldir",default=None,help="Name of polarisation caltable directory",metavar="Directory path")
	parser.add_option('--localdatabase',dest="localdatabase",default=None,help="Name of local database",metavar="Directory path")
	(options, args) = parser.parse_args()
	
os.chdir(options.basedir)
sys.path.append(os.getcwd())
import selfcal_inputs as inputs
if inputs.quality_factor==0:
	ncoarse_pol=3
elif inputs.quality_factor==1:
	ncoarse_pol=1
elif inputs.quality_factor==2:
	ncoase_pol=0
managing_caldatabase(str(options.msname),str(options.OBSID),int(options.num_jobs),str(options.basedir),str(options.gaincal_modeldir),\
						str(options.bandpass_modeldir),str(options.polcal_caldir),str(options.localdatabase),int(inputs.gain_minsnr),str(inputs.ref_ant),int(ncoarse_pol))

