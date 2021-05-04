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
from paircars.fullpol_selfcal_LTS import *
from astropy.io import fits
from CALIBRATE.access_calibrate import *


def fill_models(msname):
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

def copy_gcals_to_bcal(gcals=[],bcal='',outputname=''):
	if len(gcals)==0:
		obslog.info('No gaincal table is provided to copy.\n')
		return
	elif bcal=='':
		obslog.info('No template bandpass table is not provided to copy.\n')
		return
	else:
		


def managing_caldatabase(msname,metafits,OBSID,total_spawned_jobs,basedir,gaincal_modeldir,bandpass_modeldir,polcal_modeldir,localdatabase,gain_minsnr,ref_ant,freq_avg=160.0):
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
	obslog.info('Waiting for finishing '+str(total_spawned_jobs)+' calibration jobs for ms : '+msname+'\n')

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

	basemsdir=os.path.basename(msname).split('.ms')[0]
	gaincal_modeldir=gaincal_modeldir+'/'+basemsdir
	bandpass_modeldir=bandpass_modeldir+'/'+basemsdir
	polcal_modeldir=polcal_modeldir+'/'+basemsdir

	while True:
		touch_files=len(glob.glob(basedir+'/.Finished*cal*'+str(OBSID)+'*'+basemsdir+'*'))
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
		nchan_avg=int(freq_avg/AM.calc_freqres())
		coarse_chan_0=freq_to_MWA_coarse(freqs[0])
		coarse_chan_1=freq_to_MWA_coarse(freqs[-1])
		caltable_name_prefix=str(OBSID)+'_'+str(coarse_chan_0)+'_'+str(coarse_chan_1)
		caltable_name=caltable_name_prefix+'.gcal'
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
			bp_caltable_name=caltable_name_prefix+'.bcal'
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
				split(vis=msname,outputvis=basedir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftime.ms',width=nchan_avg,timerange=timerange,datacolumn='corrected')
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
				split(vis=msname,outputvis=basedir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftime.ms',width=nchan_avg,timerange=timerange,datacolumn='corrected')
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
			obslog.info('fill_models(\''+bpmsname+'\')\n')
			fill_models(bpmsname)
			obslog.info('bandpass(vis=\''+bpmsname+'\',caltable=\''+bp_caltable_name+'\',solnorm=True,refant=\''+str(ref_ant)+'\',minsnr='+str(gain_minsnr)+')\n')
			bandpass(vis=bpmsname,caltable=bp_caltable_name,solnorm=True,refant=str(ref_ant),minsnr=gain_minsnr)
			os.system('cp -r '+bp_caltable_name+' '+localdatabase)
			#if model_chan==unflagged_chan:
				#TODO: Updating global database
			os.system('rm -rf '+bp_caltable_name)
		if os.path.isdir(polcal_modeldir)==False or len(glob.glob(polcal_modeldir+'/*.model'))==0: # Backup polcal
			obslog.info('No polarisation caltables are available.\n')
		else:
			obslog.info('#######################################\n')
			obslog.info('Making final polcal table for ms : '+msname+'\n')
			bptables=[localdatabase+'/'+bp_caltable_name]
			obslog.info('applycal(vis=\''+bpmsname+'\',gaintable='+str(bptables)+',applymode=\'calflag\',flagbackup=False)\n')			
			applycal(vis=bpmsname,gaintable=bptables,applymode='calflag',flagbackup=False)
			if os.path.isdir(bpmsname+'.polleakcal')==True:
				os.system('rm -rf '+bpmsname+'.polleakcal')
			obslog.info('split(vis=\''+bpmsname+'\',outputvis=\''+bpmsname+'.polleakcal\',datacolumn=\'corrected\')\n')
			split(vis=bpmsname,outputvis=bpmsname+'.polleakcal',datacolumn='corrected')
			leakcor_gain_models=glob.glob(polcal_modeldir+'/*_leakage.model')
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
			model_chan,nomdel_chan=AM1.get_model_nomodel_chan()
			freqs=AM1.get_freqs()/10**6
			polleak_caltable=caltable_name_prefix+'.lcal'
			pol_caltable=caltable_name_prefix+'.pcal'
			gcals=[]
			for i in model_chan:
				obslog.info('bandpass(vis=\''+bpmsname+'.polleakcal\',caltable=\'temp_'+str(freqs[i])+'.cal\',spw=\''+spw+'\',solnorm=True,refant=\''+\
						str(ref_ant)+'\',minsnr='+str(gain_minsnr)+')\n')
				gaincal(vis=bpmsname+'.polleakcal',caltable='temp_'+str(freqs[i])+'.cal',spw=spw,solnorm=True,refant=str(ref_ant),minsnr=gain_minsnr)
				gcals.append('temp_'+str(freqs[i])+'.cal')
			os.system('cp -r '+bp_caltable_name+' temp.bcal')
			copy_gcals_to_bcal(gcals=gcals,bcal='temp.bcal',outputname=polleak_caltable)
			os.system('cp -r '+polleak_caltable+' '+localdatabase)	
			obslog.info('applycal(vis=\''+bpmsname+'.polleakcal\',gaintable='+str(polleak_caltable)+',applymode=\'calflag\',flagbackup=False)\n')		
			applycal(vis=bpmsname+'.polleakcal',gaintable=[polleak_caltable],applymode='calflag',flagbackup=False)
			polcal_ms=bpmsname+'.polcal'
			if os.path.isdir(polcal_ms)==True:
				os.system('rm -rf '+polcal_ms)
			obslog.info('split(vis=\''+bpmsname+'.polleakcal\',outputvis=\''+polcal_ms+'\',datacolumn=\'corrected\')\n')
			split(vis=bpmsname+'.polleakcal',outputvis=polcal_ms,datacolumn='corrected')
			polcal_models=glob.glob(polcal_modeldir+'/*.model')
			for i in leakcor_gain_models:
				polcal_models.remove(i)
			obslog.info('delmod(vis=\''+polcal_ms+'\',scr=True)\n')
			delmod(vis=polcal_ms,scr=True)
			for i in range(len(polcal_models)):
				modelname=polcal_models[i]
				f=imhead(imagename=modelname,mode='list')['crval4']/10**6
				df=(imhead(imagename=modelname,mode='list')['cdelt4']/10**6)/2.0
				spw='0:'+str(f-df)+'~'+str(f+df)+'MHz'
				obslog.info('ft(vis=\''+polcal_ms+'\',model=\''+modelname+'\',spw=\''+spw+'\',usescratch=True)\n')
				ft(vis=polcal_ms,model=modelname,spw=spw,usescratch=True)
			
			PSC=PolSelfcal(polcal_ms,metafits,32*60,verbose=False,interactive=False)
			obslog.info('PSC.correct_visibility_single_beam_jones(modify_datacolumn=True)\n')
			PSC.correct_visibility_single_beam_jones(modify_datacolumn=True)
			IB=ImageBasic(polcal_ms)
			calib_uvrange_min=IB1.calc_calib_uvrange(4)[1]
			calib_uvrange_max=IB1.calc_calib_uvrange(4)[2]
			cal=CALIBRATE()
			obslog.info('cal.calibrate(msname=\''+polcal_ms+'\',caltable=\''+pol_caltable+'\',minuv='+str(calib_uvrange_min)\
						+',quiet=True,maxuv='+str(calib_uvrange_max)+',j=1,absmem=1,solmode=\'\')\n')
			cal.calibrate(msname=polcal_ms,caltable=pol_caltable,minuv=calib_uvrange_min,quiet=True,maxuv=calib_uvrange_max,j=1,absmem=1,solmode='')			
			os.system('cp -r '+pol_caltable+' '+localdatabase)				
	os.system('rm -rf '+basedir+'/'+os.path.basename(msname).split('.ms')[0]+'_reftime*')
	return 0


#def run_final_imaging(msname,OBSID,channel_grid,time_grid):
	

from optparse import OptionParser
if __name__=='__main__':
	usage= ' PAIRCARS database manager'
	parser = OptionParser(usage=usage)
	parser.add_option('--msname',dest="msname",default=None,help="Name of the measurement set",metavar="Measurement set")
	parser.add_option('--metafits',dest="metafits",default=None,help="Name of the metafits file",metavar="Metafits file")
	parser.add_option('--OBSID',dest="OBSID",default=None,help="Observation ID of MWA observation",metavar="Integer")
	parser.add_option('--num_jobs',dest="num_jobs",default=0,help="Toatal number of calibration jobs spawned for this measurement set",metavar="Integer")
	parser.add_option('--basedir',dest="basedir",default=None,help="Name of base directory for a given day",metavar="Directory path")
	parser.add_option('--gaincal_modeldir',dest="gaincal_modeldir",default=None,help="Name of gaincal model directory",metavar="Directory path")
	parser.add_option('--bandpass_modeldir',dest="bandpass_modeldir",default=None,help="Name of bandpass model directory",metavar="Directory path")
	parser.add_option('--polcal_modeldir',dest="polcal_modeldir",default=None,help="Name of polarisation model directory",metavar="Directory path")
	parser.add_option('--localdatabase',dest="localdatabase",default=None,help="Name of local database",metavar="Directory path")
	(options, args) = parser.parse_args()
	
os.chdir(options.basedir)
sys.path.append(os.getcwd())
import selfcal_inputs as inputs

managing_caldatabase(str(options.msname),str(options.metafits),int(options.OBSID),int(options.num_jobs),str(options.basedir),str(options.gaincal_modeldir),\
str(options.bandpass_modeldir),str(options.polcal_modeldir),str(options.localdatabase),float(inputs.gain_minsnr),str(inputs.ref_ant))

