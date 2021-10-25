import os,psutil
import numpy as np,sys,matplotlib.pyplot as plt,time,logging,matplotlib,json,urllib.request,glob
from casatools import msmetadata,table,measures,quanta,agentflagger,image,calibrater,ms
from casatasks import *
from paircars.access_ms import *
from paircars.basic_func import *
from paircars.fullpol_selfcal_LTS import *
from paircars.flagger import *
from optparse import OptionParser
from astropy.io import fits
from astropy import wcs
from CALIBRATE.access_calibrate import *
from paircars.libpaircars import send_paircars_notification
from mwa_pb.mwapb import *
matplotlib.use('Agg')

if __name__!='__main__':
	cwd=os.getcwd()
	sys.path.append(cwd)
	if os.path.isfile(cwd+'/selfcal_inputs.py')==False:
		print ('Input file does not exist.\n')
		os._exit(0)
	else:
		import selfcal_inputs as inputs
		from selfcal_inputs import *	

'''
Code is written by Devojyoti Kansabanik, March 6, 2021
'''
def get_OBSID(metafits):
	'''
	Function to return OBSID of an MWA observation
	Parameters:
	metafits = Name of the metafits file
	MWA OBSID
	'''
	OBSid=fits.getheader(metafits)['GPSTIME']
	return OBSid  

def get_quicklook_image(imagename,outfile,freq,timestamp,DR_rms,DR_neg,field_of_view=2): 
	'''
	Function to get a quick look image
	Parameters:
	imagename = Name of the CASA image
	outfile = Output file name
	freq = Frequency in MHz
	timestamp = Timestamp string
	DR_rms = Dynamic range of the image based on rms
	DR_neg = Dynamic range of the image based on negative flux near the source
	field_of_view = Field of view to cut the image in degree (default : 2)
	Return:
	Outfile name
	'''
	header=imhead(imagename=imagename,mode='list')
	xcent=int(header['shape'][0]/2)
	ycent=int(header['shape'][1]/2)
	cell=np.rad2deg(abs(header['cdelt2'])) # In degree
	freq="{:.2f}".format(float(freq))
	xwidth=ywidth=int((field_of_view)/cell)
	box=str(xcent-int(xwidth/2))+','+str(ycent-int(ywidth/2))+','+str(xcent+int(xwidth/2))+','+str(ycent+int(ywidth/2))
	os.system('rm -rf temp*')
	stokes_list=['I','Q','U','V']
	fig = plt.figure(figsize=(8,8))
	plt.subplots_adjust(wspace=0.45, hspace=0.1)
	for i in range(len(stokes_list)):
		stokes=stokes_list[i]
		try:
			imsubimage(imagename=imagename,outfile='temp_'+stokes+'.image',box=box,stokes=stokes)
		except:
			return 
		exportfits(imagename='temp_'+stokes+'.image',fitsimage='temp_'+stokes+'.fits',dropdeg=True,dropstokes=True)
		data=fits.getdata('temp_'+stokes+'.fits')
		wlist=fits.getheader('temp_'+stokes+'.fits')
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
	title='Frequency : '+str(freq)+' MHz, Timestamp : '+str(timestamp)+' UTC\n Dynamic range I (rms) : '+str(int(DR_rms))+', Dynamic range I (negative) : '+str(int(DR_neg))
	plt.suptitle(title,fontsize=12)	
	cwd=os.getcwd()
	outfile_dir=os.path.dirname(outfile)
	if outfile_dir=='':
		outfile=cwd+'/'+outfile
	plt.savefig(outfile)
	os.system('rm -rf temp* casa*log')
	return outfile

def run_pol_selfcal(msname,metafits,working_dir,verbose=False,interactive=False,start_fresh=True,perform_gaincal=False,caltables='',use_wsclean=True):
	'''
	Heart of the polarisation selfcal part of the PAIRCARS
	This function performs the polarisation selfcal for PAIRCARS
	This script can be run directly from python IDE or can be imported to other scripts.
	Here some functionalities are chosen specific to MWA solar imaging. Which may not be valid for other instruments and imaging of other sources.
	Use this module only for Solar Imaging with MWA.
	Parameters:
	msname = Name of the measurement set
	metafits = Name of observation metafits file
	working_dir = Name of the working directory
	verbose = False, If True keep all intermediate selfcal records
	interactive = False, If True perform interactive selfcal
	start_fresh = True, start fresh selfcal rounds from scratch or start from last round
	perform_gaincal = False, perform gaincal using leakage corrected model (Only do when no calibrator observation is present)
	caltables = Previous caltables, comma separated
	use_wsclean = Use WSClean for imaging or not
	Return:
	Meassages about the selfcal success or errors
	'''
	if __name__!='__main__':
		start_time=time.time()
	if working_dir[-1]=='/':
		working_dir=working_dir[:-1]

	cwd=os.getcwd()
	os.chdir(msname)
	mspath=os.path.dirname(os.path.realpath(os.getcwd()))
	os.chdir(cwd)
	DR_delta_rms=inputs.DR_delta_rms
	if mspath[-1]=='/':
		mspath=mspath[:-1] 

	if mspath!=working_dir:
		if os.path.isdir(working_dir+'/'+os.path.basename(msname)):
			os.system('rm -rf '+working_dir+'/'+os.path.basename(msname))
		os.system('mv '+msname+' '+working_dir)
		msname=working_dir+'/'+os.path.basename(msname)

	os.chdir(working_dir)
	if __name__!='__main__':
		if (os.path.isfile(working_dir+'/Pol_Selfcal.log') and start_fresh==True) or \
				(os.path.isfile(working_dir+'/Pol_Selfcal.log') and os.path.isdir(working_dir+'/junk1.ms')==False and start_fresh==False):
			os.system('rm -rf '+working_dir+'/Pol_Selfcal.log')
		if os.path.isfile(working_dir+'/Pol_Selfcal_verbose.log') and verbose==True:
				os.system('rm -rf '+working_dir+'/Pol_Selfcal_verbose.log')
		formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
	logger = logging.getLogger('pol_selfcal_log')
	if __name__!='__main__':
		logger.setLevel(logging.DEBUG)
		if verbose==True:
			console=logging.StreamHandler(sys.stdout)
			console.setFormatter(formatter)
			logger.addHandler(console)
		filehandle=logging.FileHandler(working_dir+'/Pol_Selfcal.log')
		filehandle.setFormatter(formatter)
		logger.addHandler(filehandle)
		logger.propagate = False
	print('\n')


	if start_fresh==False and os.path.isdir('junk1.ms'):
		os.system('rm -rf '+msname)
		os.system('cp -r junk1.ms '+msname)
	else:
		start_fresh=True
		if os.path.isfile('IMGSTAT_pol.npy')==True:
			os.system('rm -rf IMSTAT_pol.npy')   # TODO : check IMSTAT_pol.npy not recording all iteration
		if os.path.isfile(msname+'/.usedby_paircars'):
			try:
				tb=table()
				tb.open(msname,nomodify=False)
				flag=tb.getcol('FLAG')
				flag*=False
				tb.putcol('FLAG',flag)
				tb.flush()
				tb.close()
				os.system('rm -rf '+msname+'.flagversions')	
			except:
				pass
		if caltables!='':
			AM=AccessMS(msname)
			msfreq=AM.calc_meanfreq()/10**6
			caltable_list=caltables.split(',')
			bpcaltables=[]
			for i in caltable_list:
				if '.bcal' in i:
					bpcaltables.append(i)
			bp_freqs=[float(os.path.basename(i).split('freq_')[-1].split('_ref')[0].split('.bcal')[0]) for i in bpcaltables]
			bptable=bpcaltables[np.argmin(abs(msfreq-np.array(bp_freqs)))]
			for i in caltable_list:
				if '.bcal' in i:
					if bptable!=i:
						caltable_list.remove(i)
			caltables=','.join(caltable_list)
			logger.info('Applying solutions from previous calibrations : '+str(caltables)+'\n')
			logger.info('applycal(vis=\''+msname+'\',gaintable='+str(caltable_list)+',applymode=\'calflag\',calwt=[False],flagbackup=False)\n')
			applycal(vis=msname,gaintable=caltable_list,applymode='calflag',calwt=[False],flagbackup=True)
		flagged_chans,flag_frac=calc_flag_chans_caltable(bptable,flag_frac=1.0)
		min_flagged_chan=np.argmin(flag_frac)
		logger.info('split(vis=\''+msname+'\',outputvis=\''+msname+'.temp\',datacolumn=\'corrected\',spw=\'0:'+str(min_flagged_chan)+'\')\n')
		try:
			split(vis=msname,outputvis=msname+'.temp',datacolumn='corrected',spw='0:'+str(min_flagged_chan))
		except Exception as e:
			logger.info('Split error : '+str(e)+'\n')	
		os.system('rm -rf '+msname+' '+msname+'.flagversions')
		if 'ref' in msname:
			ref_time_chan=True
		else:
			ref_time_chan=False
		msname=splited_ms_rename(msname+'.temp',ref_time_chan=ref_time_chan,change_msname=True)
		logger.info('splited_ms_rename(\''+msname+'.temp\',ref_time_chan='+str(ref_time_chan)+',change_msname=True)\n')
		logger.info('Polarisation calibration ms : '+msname+'\n')		

	if msname[-1]=='/':
		msname=msname[:-1]
	if inputs.basedir[-1]=='/':
		basedir=inputs.basedir[:-1]
	else:
		basedir=inputs.basedir

	if __name__!='__main__' and inputs.send_notification==True:
		OBSID=get_OBSID(metafits)

	if 'ref' in msname:
		msname_str=os.path.basename(splited_ms_rename(msname,ref_time_chan=True,change_msname=False))
	else:
		msname_str=os.path.basename(splited_ms_rename(msname,ref_time_chan=False,change_msname=False)) 
	freqstr=msname_str.split('.ms')[0].split('_freq_')[1].split('_')[0]  # Frequency string in MHz
	datestr_list=msname_str.split('.ms')[0].split('_freq_')[0].split('time_')[1].split('_')
	datestr='/'.join(datestr_list[:3])+'/'+':'.join(datestr_list[3:]) # Datetime string s
	datestrfile='_'.join(datestr_list[:3])+'_'+'_'.join(datestr_list[3:]) # Datetime string 
	file_str_prefix='freq_'+freqstr+'_datetime_'+datestrfile+'_pol'

	OBSID=get_OBSID(metafits)
	basemsdir=os.path.dirname(working_dir).split('/')[-1]
	c=0
	while c<=10:
		c+=1
		try:
			if os.path.isdir(basedir+'/polcaltables/'+str(OBSID)+'/'+basemsdir)==False: # Directory to keep caltables
				os.makedirs(basedir+'/polcaltables/'+str(OBSID)+'/'+basemsdir)
			if os.path.isdir(basedir+'/polimagemodels/'+str(OBSID)+'/'+basemsdir)==False: # Directory to keep models
				os.makedirs(basedir+'/polimagemodels/'+str(OBSID)+'/'+basemsdir)
			if os.path.isdir(basedir+'/polms/'+str(OBSID)+'/'+basemsdir)==False: # Directory to keep models
				os.makedirs(basedir+'/polms/'+str(OBSID)+'/'+basemsdir)
			if os.path.isdir(basedir+'/logs/'+str(OBSID)+'/'+basemsdir)==False: # Directory to keep models
				os.makedirs(basedir+'/logs/'+str(OBSID)+'/'+basemsdir)
			if os.path.isdir(basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir)==False and inputs.keep_logger==True and verbose==True:
				os.makedirs(basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir)
			break
		except:
			time.sleep(2.0)
			pass

	if 'ref' in msname:
		refcals=glob.glob(basedir+'/polcaltables/'+str(OBSID)+'/'+basemsdir+'/*ref*')
		refimages=glob.glob(basedir+'/polimagemodels/'+str(OBSID)+'/'+basemsdir+'/*ref*')
		refms=glob.glob(basedir+'/polms/'+str(OBSID)+'/'+basemsdir+'/*ref*')
		if len(refcals)!=0:
			for i in refcals:
				os.system('rm -rf '+i)
		if len(refimages)!=0:
			for j in refimages:
				os.system('rm -rf '+j)
		if len(refms)!=0:
			for k in refms:
				os.system('rm -rf '+k)

	#################
	# Searching for previous solutions in not ref time 
	#################
	nearest_freq_leakcal=[]
	nearest_freq_polcal=[]
	
	if start_fresh==False: # Reading selfcal record
		num_iter,DR1,DR3,DR5,DR2,DR4,DR6,FX3_I,FX3_Q,FX3_U,FX3_V,FX3_T,FX3_P,FX2_I,FX2_Q,FX2_U,FX2_V,FX2_T,FX2_P,FX1_I,FX1_Q,FX1_U,FX1_V,FX1_T,FX1_P,\
		rms_list,scratch,start_sigma,num_iteration_after_poldist,num_iter_after_qucor,\
		num_iter_fixed_sigma,startmodel,startmask,uvsub_flag_count,do_solarqu_cor,do_poldist,do_pbcor,do_gaincal,gaincal_count,done_qucor,pre_res\
			=np.load('Pol_selfcal_record.npy',allow_pickle=True)
	else:
		num_iter=0
	
	if 'ref' in msname:
		if start_fresh:
			scratch=True # For reference time and frequency scratch = True
		if start_fresh==False:
			logger.info('Restarting selfcal from selfcal round : '+str(num_iter)+'\n')
		logger.info('Starting imaging for Reference time : '+str(datestr)+' and frequency : '+str(freqstr)+' MHz\n')
		logger.info('####################################\n')
		logger.info('Scratch = '+str(scratch)+'\n')
		if verbose==False:
			if start_fresh==False:
				print('Restarting selfcal from selfcal round : '+str(num_iter)+'\n')
			print ('Starting imaging for Reference time : '+str(datestr)+' and frequency : '+str(freqstr)+' MHz\n')
			print('####################################\n')
			print ('Scratch = '+str(scratch)+'\n')
	else: # For other time and frequency
		if start_fresh: 
			scratch=False
		if scratch==False:
			OBSID_list=[int(os.path.basename(i)) for i in glob.glob(basedir+'/polcaltables/*')]
			nearest_OBSID=OBSID_list[np.argmin(np.abs(OBSID-np.array(OBSID_list)))]				
			polcal=np.array(glob.glob(basedir+'/polcaltables/'+str(nearest_OBSID)+'/*'))
			freqs_polcals=np.array([float(os.path.basename(i).split('freq_')[-1]) for i in polcal])
			basems_freq=float(os.path.basename(basemsdir).split('freq_')[-1])
			if basems_freq in freqs_polcals:
				pos=np.where(freqs_polcals==basems_freq)
				ref_time_chan_leakcal=glob.glob(str(polcal[pos][0])+'/*.lcal')
				ref_time_chan_polcal=glob.glob(str(polcal[pos][0])+'/*.bin')
				ref_time_chan_polcal_copy=copy.deepcopy(ref_time_chan_polcal)
				for i in ref_time_chan_polcal:
					if 'beam' in i:
						ref_time_chan_polcal_copy.remove(i)
				ms_freq=str(freqstr)
				ms_freq_coarse=freq_to_MWA_coarse(ms_freq)
				for i in ref_time_chan_leakcal:
					f1=float(os.path.basename(i).split('freq_')[-1].split('_ref')[0].split('.lcal')[0])
					f1_coarse=freq_to_MWA_coarse(f1)
					if f1_coarse==ms_freq_coarse:
						nearest_freq_leakcal.append(i)
				ref_time_chan_polcal=copy.deepcopy(ref_time_chan_polcal_copy)
				for i in ref_time_chan_polcal:
					f1=float(os.path.basename(i).split('freq_')[-1].split('_ref')[0].split('.bin')[0])
					f1_coarse=freq_to_MWA_coarse(f1)
					if f1_coarse==ms_freq_coarse:
						nearest_freq_polcal.append(i)
			
		if len(nearest_freq_leakcal)==0 or len(nearest_freq_polcal)==0 and num_iter==0:
			scratch=True
			logger.info('Change to scratch = True\n')
		if start_fresh==False:
			logger.info('Restarting selfcal from selfcal round : '+str(num_iter)+'\n')
		logger.info('Reference time frequency slice imaging has been done. Starting imaging for time : '+str(datestr)+' and frequency : '+str(freqstr)+' MHz\n')
		logger.info('####################################\n')
		logger.info('Scratch = '+str(scratch)+'\n')
		if verbose==False:
			if start_fresh==False:
				print('Restarting selfcal from selfcal round : '+str(num_iter)+'\n')
			print ('Reference time frequency slice imaging has been done. Starting imaging for time : '+str(datestr)+' and frequency : '+str(freqstr)+' MHz\n')
			print('####################################\n')
			print ('Scratch = '+str(scratch)+'\n')
	
	if inputs.quality_factor==0:
		num_pixel_in_psf=3
	elif inputs.quality_factor==1:
		num_pixel_in_psf=5
	else:
		num_pixel_in_psf=7
			
	PSC=PolSelfcal(msname,metafits,32*60,num_pixel_in_psf=num_pixel_in_psf,largest_scale=12,verbose=verbose,interactive=interactive,use_wsclean=use_wsclean,savelog=inputs.keep_logger) 
						# Creating selfcal object 32 arcmin maximum scale size
	AM=AccessMS(msname)
	
	###################
	# Putting user defined inputs if exisis or go with default values
	###################

	if calc_image_parameters==False:
		if cellsize!='' and cellsize!='nan':
			PSC.cellsize=inputs.cellsize
		if len(imsize)!=0 and imsize[0]!='nan':
			PSC.imsize=inputs.imsize
		if len(multiscale_scales)!=0 and multiscale_scales[0]!='nan':
			PSC.multiscale_scales=inputs.multiscale_scales
		if uvtaper!='':
			PSC.uvtaper=inputs.uvtaper
		if inputs.weight=='' or (inputs.weight!='uniform' and inputs.weight!='natural' and inputs.weight!='briggs'):
			weight='briggs'
		else:
			weight=inputs.weight
		if weight=='briggs':
			if inputs.robust<-1.0 or inputs.robust>1.0:
				robust=1.0
			else:
				robust=inputs.robust
	else:
		weight='briggs'
		robust=0.8
		
	if calc_selfcalib_params==False:
		if uvrange_to_cal!='':
			if '~' in uvrange_to_cal:
				uvrange=uvrange_to_cal.split('~')
				uvmin=float(uvrange[0])*AM.calc_meanwavelength()
				uvmax=float(uvrange[1])*AM.calc_meanwavelength()
				PSC.calib_uvrange_min=uvmin
				PSC.calib_uvrange_max=uvmax
			elif '>' in uvrange_to_cal:
				uvrange=uvrange_to_cal.split('>')
				uvmin=float(uvrange[1])*AM.calc_meanwavelength()
				uvmax=AM.get_max_baseline()
				PSC.calib_uvrange_min=uvmin
				PSC.calib_uvrange_max=uvmax
			elif '<' in uvrange_to_cal:
				uvrange=uvrange_to_cal.split('<')
				uvmin=0
				uvmax=float(uvrange[1])*AM.calc_meanwavelength()
				PSC.calib_uvrange_min=uvmin
				PSC.calib_uvrange_max=uvmax
	
	end_selfcal=False
	if start_fresh:
		do_gaincal=False
		gaincal_count=0
		num_iter=0
		done_qucor=False
		pre_res=0.0
	else:
		do_gaincal=do_gaincal
		gaincal_count=gaincal_count
		num_iter=num_iter
		done_qucor=done_qucor
		pre_res=pre_res

	try:
		start_sigma=np.load(basedir+'/Ref_time_chan_sigma.npy',allow_pickle=True)[0] # Starting with last gaincal start_sigma and threshold
		rms_list=np.load(basedir+'/Ref_time_chan_sigma.npy',allow_pickle=True)[1]
		rms=np.mean(np.array(rms_list)) # TODO : calculate Stokes I beam and divide
		rms_list=[rms]*4			
	except:
		logger.info('Start sigma and threshold information for last intensity selfcal round for reference time channel is not found.\n')
		if verbose==False:
			print('Start sigma and threshold information for last intensity selfcal round for reference time channel is not found.\n')
		os.chdir(cwd)
		if __name__!='__main__':
			touch_file=basedir+'/.Finished_pcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_12'
			msg_str='Dear PAIRCARS user,\n\nPolarisation self-calibration for : '+\
					os.path.basename(msname)+'\nMessage :'+error_msgs(12)+'\n\nBest regards,\nPAIRCARS developing team'
			msg_subject='Notification from PAIRCARS : Polarisation Selfcal : OBSID = '+str(OBSID)
			if inputs.send_notification==True:
				send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[])
			os.system('touch '+touch_file)
			os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_gaincaled.ms')
			end_time=time.time()
			run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
			logger.info('Total runtime : '+str(run_time))
			os.system('cp -r '+working_dir+'/Pol_Selfcal.log '+basedir+'/logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
			if inputs.keep_logger and verbose==True:
				os.system('cp -r '+working_dir+'/Pol_Selfcal_verbose.log '+basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
			os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice* '+working_dir+'/I*image')
		return 12

	background_mask_rad=int((200*60)/PSC.cellsize) # Creating a mask with 3deg radius centered on the image
	background_mask_str='circle[['+str(PSC.imsize/2)+'pix,'+str(PSC.imsize/2)+'pix],'+str(background_mask_rad)+'pix]'

	subms,subtracted=PSC.subtract_background_sources('IQUV',rms_list,start_sigma+1.5,start_sigma+1.5,maskregion=background_mask_str,includeregion=False,\
													overwrite=True,modify_datacolumn=True,weight=weight,robust=robust)
	cal=CALIBRATE()
	mwa_config=get_MWA_phase(metafits) # TODO : Include from cross phase cal solutions
	if mwa_config=='MWAPhaseI':
		crossphase=15
	elif mwa_config=='MWAPhaseIILB' or mwa_config=='MWAPhaseIICOMPACT':
		crossphase=135
	logger.info('Applying cross hand phase solution. Cross hand phase : '+str(crossphase)+' deg.\n')
	PSC.apply_cross_hand_phase(cross_phase=crossphase,caltable='',polbasis='Linear',modify_datacolumn=True)
	if os.path.isdir(working_dir+'/Backup_gaincaled.ms')==True:
		os.system('rm -rf '+working_dir+'/Backup_gaincaled.ms')
	logger.info('split(vis=\''+msname+'\',outputvis=\''+working_dir+'/Backup_gaincaled.ms\',datacolumn=\'data\')\n')
	try:
		split(vis=msname,outputvis=working_dir+'/Backup_gaincaled.ms',datacolumn='data') # Backup of gain calibrated ms
	except Exception as e:
		logger.info('Split error : '+str(e)+'\n')

	while end_selfcal==False:
		if os.path.isfile(msname+'/.usedby_paircars')==False:
			os.system('touch '+msname+'/.usedby_paircars')

		###############################################################
		# Calculating minimum number of iterations and antenna bin size
		###############################################################
		min_num_iter_fixed_sigma,min_iteration,max_iteration,antenna_bin,frac_flux_change,pol_frac_change=PSC.calc_iter_num(inputs.safety_factor,inputs.quality_factor,scratch=scratch)
		if verbose==False:
			print ('########################\nEstimating the number of selfcal iterations\n########################\n')
		logger.info('Estimating the minimum number of selfcal iterations\n')
		logger.info('Minimum number of iteration at fixed sigma : '+str(min_num_iter_fixed_sigma)\
			+', Minimum number of total selfcal iterations : '+str(min_iteration)+', Maximum number of selfcal iterations : '+str(max_iteration)+\
			', Antenna bins : '+str(antenna_bin)+'\n')
		logger.info('#####################################\n')

		antenna_list,num_ant=AM.make_antenna_list(num_bins=antenna_bin)  # Making the antenna list
		
		###########################
		# Initiating loop variables
		###########################
		quvcor_stokes='QU'
		if start_fresh:
			do_selfcal=True
			stokes='IQUV'
			do_poldist=False
			if (do_gaincal==False) and ((num_iter==0 and perform_gaincal==False) or (perform_gaincal==True)):	
				do_pbcor=True	
				do_solarqu_cor=False
				num_iteration_after_poldist=0	
				uvsub_flag_count=0
				num_iter_fixed_sigma=0
				num_iter_after_qucor=0				
				startmodel=''
				startmask=''
		else:
			do_selfcal=True
			stokes='IQUV'
			do_poldist=do_poldist
			do_pbcor=do_pbcor	
			do_solarqu_cor=do_solarqu_cor
			num_iteration_after_poldist=num_iteration_after_poldist	
			uvsub_flag_count=uvsub_flag_count
			num_iter_fixed_sigma=num_iter_fixed_sigma
			num_iter_after_qucor=num_iter_after_qucor			
			startmodel=startmodel
			startmask=startmask

		previous_image=''
		previous_model=''
		file_str=os.path.basename(msname).split('.ms')[0]
	
		if inputs.maskfile=='' and inputs.maskstr=='':
			mask_rad=int((60*60)/PSC.cellsize) # Creating a mask with 60 arcmin radius centered on the image
			mask_str='circle[['+str(PSC.imsize/2)+'pix,'+str(PSC.imsize/2)+'pix],'+str(mask_rad)+'pix]'
		elif inputs.maskstr!='':
			mask_str=inputs.maskstr

		###################
		# Performing gaincal again using the new leakage corrected source model
		###################
		if do_gaincal==True and gaincal_count<1:
			os.system('rm -rf '+msname+' '+msname+'.flagversions')
			logger.info('Spliting from gaincaled backup ms.....\n')
			logger.info('split(vis=\'Backup_gaincaled.ms\',outputvis=\''+msname+'\',datacolumn=\'data\')\n')
			try:
				split(vis='Backup_gaincaled.ms',outputvis=msname,datacolumn='data')
			except Exception as e:
				logger.info('Split error : '+str(e)+'\n')
			os.system('cp -r '+msname+' '+basedir+'/polms/'+str(OBSID)+'/'+basemsdir+'/Gaincaled.ms')
			if os.path.isfile(msname+'/.usedby_paircars')==False:
				os.system('touch '+msname+'/.usedby_paircars')
			logger.info('Performing gain calibration using the leakage corrected source model.\n')
			if verbose==False:
				print('Performing gain calibration using the leakage corrected source model.\n')
			logger.info('mwapb=MWA_PrimaryBeam(\''+msname+'\',\''+metafits+'\',\'inverse_beam=True)\n')
			mwapb=MWA_PrimaryBeam(msname,metafits,inverse_beam=True)  # Inverse beam jones
			logger.info('mwapb.calc_beamjones_phasecenter(outputfile=\'\')\n')
			inv_beam_jones=mwapb.calc_beamjones_phasecenter(outputfile='')[0]
			os.system('cp -r '+working_dir+'/quvcor.image '+working_dir+'/junk1.image')
			os.system('cp -r '+working_dir+'/quvcor.model '+working_dir+'/junk1.model')
			logger.info('PSC.uncorrect_for_single_beam_jones(\''+working_dir+'/quvcor.model\',\''+working_dir+'/quvcor_pbuncor.model\','+\
						'inv_beam_jones,imagetype=\'CASA\',outtype=\'CASA\',pol_basis=\'Linear\')\n')
			PSC.uncorrect_for_single_beam_jones(working_dir+'/quvcor.model',working_dir+'/quvcor_pbuncor.model',inv_beam_jones,imagetype='CASA',outtype='CASA',pol_basis='Linear')
			modelname=working_dir+'/quvcor_pbuncor.model'
			os.system('cp -r '+modelname+' '+basedir+'/polimagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'_leakage.model') # Keeping new leakage model		
			logger.info('delmod(vis=\''+msname+'\',scr=True,otf=True)\n')
			delmod(vis=msname,scr=True,otf=True)
			logger.info('ft(vis=\''+msname+'\',model=\''+modelname+'\',usescratch=True)\n')
			ft(vis=msname,model=modelname,usescratch=True)
			IB=ImageBasic(msname)
			uvrange=IB.calc_calib_uvrange(12)[0]
			logger.info('bandpass(vis=\''+msname+'\',caltable=\''+working_dir+'/Leakage_cor_gaincal.cal\',refant='+str(inputs.ref_ant)+',minsnr='+str(inputs.gain_minsnr)+','+\
					',uvrange=\''+uvrange+'\',bandtype=\'B\')\n')
			bandpass(vis=msname,caltable=working_dir+'/Leakage_cor_gaincal.cal',refant=str(inputs.ref_ant),minsnr=inputs.gain_minsnr,uvrange=uvrange,bandtype='B')
			applycal_caltable=[working_dir+'/Leakage_cor_gaincal.cal']
			logger.info('applycal(vis=\''+msname+'\',gaintable='+str(applycal_caltable)+',applymode=\'calflag\',calwt=[False],flagbackup=False)\n')
			applycal(vis=msname,gaintable=applycal_caltable,applymode='calflag',calwt=[False],flagbackup=False)
			logger.info('split(vis=\''+msname+'\',outputvis=\''+msname+'.temp\',datacolumn=\'corrected\')\n')
			try:
				split(vis=msname,outputvis=msname+'.temp',datacolumn='corrected')
			except Exception as e:
				logger.info('Split error : '+str(e)+'\n')
			os.system('rm -rf '+msname+' '+msname+'.flagversions')
			logger.info('os.system(\'mv '+msname+'.temp '+msname+')\n')
			os.system('mv '+msname+'.temp '+msname)
			os.system('cp -r '+msname+' '+basedir+'/polms/'+str(OBSID)+'/'+basemsdir+'/Leakcor_gaincaled.ms')
			os.system('cp -r '+working_dir+'/Leakage_cor_gaincal.cal '+basedir+'/polcaltables/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.lcal') 
																															# Keeping new leakage corrected gaintable caltable		
			logger.info('Gaincal using leakage corrected source model is done.\n')
			if verbose==False:
				os.system('rm -rf '+working_dir+'/Leakage_cor_gaincal.cal')
			else:
				os.system('mv '+working_dir+'/Leakage_cor_gaincal.cal '+working_dir+'/'+file_str_prefix)
			do_pbcor=True
			gaincal_count+=1
		elif done_qucor==True and num_iter_after_qucor<1:
			do_pbcor=True
			gaincal_count+=1
		
		if (do_pbcor==True and gaincal_count==1) or (do_pbcor==True):
			if verbose==False:
				print ('Performing ideal beam correction.\n')
			logger.info('Performing ideal beam correction.\n')
			if gaincal_count==1:
				force=True
			else:
				force=False
			if os.path.isdir(working_dir+'/Backup_beamcorrected.ms')==False:
				force=True
			else:
				force=False
			save_beamfile=basedir+'/polcaltables/'+str(OBSID)+'/'+basemsdir+'/'+os.path.basename(msname).split('.ms')[0]+'_beam.bin'
			if os.path.exists(working_dir+'/temp_beamcorrected.ms'):
				os.system('rm -rf '+working_dir+'/temp_beamcorrected.ms')
				force=True
			logger.info('PSC.correct_visibility_single_beam_jones(modify_datacolumn=False,force='+str(force)+',skip_freq=1.28,save_beamfile=\''+str(save_beamfile)+'\')\n')
			PSC.correct_visibility_single_beam_jones(modify_datacolumn=False,force=force,skip_freq=1.28,save_beamfile=save_beamfile) # Single pointing beam correction on visibility data
			if os.path.isdir(working_dir+'/temp_beamcorrected.ms'):
				os.system('rm -rf '+working_dir+'/temp_beamcorrected.ms')			
			logger.info('split(vis=\''+msname+'\',outputvis=\''+working_dir+'/temp_beamcorrected.ms\',datacolumn=\'corrected\')\n')
			try:
				split(vis=msname,outputvis=working_dir+'/temp_beamcorrected.ms',datacolumn='corrected') # Backup beam corrected visibility
				if os.path.isdir(working_dir+'/Backup_beamcorrected.ms'):
					os.system('rm -rf '+working_dir+'/Backup_beamcorrected.ms')
				os.system('mv '+working_dir+'/temp_beamcorrected.ms '+working_dir+'/Backup_beamcorrected.ms')
			except Exception as e:
				logger.error('Error occured in spliting : '+str(e)+'\n')
				os.system('rm -rf '+working_dir+'/temp_beamcorrected.ms')
			tb=table()
			tb.open('Backup_beamcorrected.ms')
			beamcor_data=tb.getcol('DATA')
			tb.close()
			tb.open(msname,nomodify=False)
			try:
				tb.putcol('DATA',beamcor_data)
			except:
				pass
			try:
				tb.putcol('CORRECTED_DATA',beamcor_data)
			except:
				pass
			tb.flush()
			tb.close()
			os.system('cp -r '+msname+' '+basedir+'/polms/'+str(OBSID)+'/'+basemsdir+'/Beamcorrected.ms')
			if gaincal_count==1:
				logger.info('Re-calibrating using Stokes Q,U corrected model.\n')
				IB1=ImageBasic(msname)
				calib_uvrange_min=IB1.calc_calib_uvrange(12)[1]
				calib_uvrange_max=IB1.calc_calib_uvrange(12)[2]
				clearcal(vis=msname,addmodel=True)
				logger.info('delmod(vis=\''+msname+'\',scr=True,otf=True)\n')
				delmod(vis=msname,scr=True,otf=True)
				logger.info('ft(vis=\''+msname+'\',model=\''+working_dir+'/quvcor.model\',usescratch=True)\n')
				ft(vis=msname,model=working_dir+'/quvcor.model',usescratch=True)
				logger.info('cal.calibrate(msname=\''+msname+'\',caltable=\''+working_dir+'/Leakage_corrected.bin\',minuv='+str(calib_uvrange_min)+',quiet='+str(verbose)+\
						',maxuv='+str(calib_uvrange_max)+',j=1,absmem=1)\n')
				cal.calibrate(msname=msname,caltable=working_dir+'/Leakage_corrected.bin',minuv=calib_uvrange_min,quiet=verbose,\
						maxuv=calib_uvrange_max,j=1,absmem=1)
				corrected_gaintable=working_dir+'/Leakage_corrected.bin'
				logger.info('cal.applycal(msname=\''+msname+'\',gaintable=\''+corrected_gaintable+'\',applymode=\'calflag\',flagbackup=True)\n')
				cal.applycal(msname=msname,gaintable=corrected_gaintable,applymode='calflag',flagbackup=True) # Applying the solution
				tb=table()
				tb.open(msname)
				leakage_cor_data=tb.getcol('CORRECTED_DATA')
				tb.close()
				tb.open(msname,nomodify=False)
				tb.putcol('DATA',leakage_cor_data)
				tb.flush()
				tb.close()
				if verbose:
					os.system('mv '+working_dir+'/Leakage_corrected.bin '+working_dir+'/'+file_str_prefix)
					os.system('mv '+working_dir+'/quvcor* '+working_dir+'/'+file_str_prefix)
				else:
					os.system('rm -rf '+working_dir+'/Leakage_corrected.bin')
					os.system('rm -rf '+working_dir+'/quvcor*')
			do_pbcor=False

		####################
		# scratch management
		####################
		IB=ImageBasic(msname)
		os.system('cp -r '+working_dir+'/Backup_beamcorrected.ms '+' '+working_dir+'/beamcor_backup.ms')
		if num_iter==0:
			if os.path.isdir(working_dir+'/beamcor_backup.ms') and gaincal_count<1:
				os.system('rm -rf '+msname+' '+msname+'.flagversions')
				os.system('mv '+working_dir+'/beamcor_backup.ms '+msname)
		
		##############
		# Selfcal loop
		##############
		antenna_to_use=PSC.antenna_string(antenna_list,-1)
		while do_selfcal==True:
			if num_iter==min(min_iteration,10):
				do_poldist=True
			else:
				do_poldist=False
			poldistortion_type='poldistortion'	
			do_flag=True
			if subtracted==True and num_iter<=1:
				do_poldist=False
				do_flag=True
			
			if verbose==False:
				print ('#####################\nPolarisation Selfcal iteration:'+str(num_iter)+'\n#####################\n')
			logger.info('#####################\n')
			logger.info('Polarisation Selfcal iteration : '+str(num_iter)+'\n')
			logger.info('#####################\n')

			if os.path.isdir(startmodel)==False:
				startmodel=''
			if os.path.isdir(startmask)==False:
				startmask=''
			if num_iter_after_qucor<2 and (do_solarqu_cor==True or done_qucor==True):
				startmodel=''
				startmask=''
			if num_iter>min(min_iteration,10):
				if num_iter==1:
					previous_image='junk0.image'
					previous_model='junk0.model'
				elif num_iter>1:
					previous_image='junk1.image'
					previous_model='junk1.model'
				else:
					previous_image=''
					previous_model=''		
			
			if done_qucor==False:
				if num_iter<=min(min_iteration,10):
					polmodel_threshold=start_sigma*1.5
				else:
					polmodel_threshold=start_sigma*1.2
			else:
				polmodel_threshold=start_sigma

			if num_iter==min(min_iteration,10):
				do_solarqu_cor=True
				quvcor_stokes='QUV'
				previous_image=''
				previous_model=''		
			
			if inputs.maskfile!='': # Use user defined mask
				mask_str=''
				output_PSC=PSC.polselfcal_iteration(num_iter,rms_list,mask_str,start_sigma,maskfile,antenna_to_use,startmodel,startmask,want_auto_masking=False,\
				stokes=stokes,interactive=interactive,use_ankflagger=inputs.use_ankflagger,do_flag=do_flag,poldistortion_correction=do_poldist,poldistortion_type=poldistortion_type,\
						poldistortion_matrix='UH',calibrator_caltable=[],box_width=3,previous_image=previous_image,previous_model=previous_model,do_solarqu_cor=do_solarqu_cor,\
						polmodel_threshold=polmodel_threshold,quvcor_stokes=quvcor_stokes,weight=weight,robust=robust)  		
			elif inputs.maskstr!='': # If mask user defined string is given
				mask_str=inputs.maskstr
				output_PSC=PSC.polselfcal_iteration(num_iter,rms_list,mask_str,start_sigma,maskfile,antenna_to_use,startmodel,startmask,want_auto_masking=False,\
					stokes=stokes,interactive=interactive,use_ankflagger=inputs.use_ankflagger,do_flag=do_flag,poldistortion_correction=do_poldist,poldistortion_type=poldistortion_type,\
					poldistortion_matrix='UH',calibrator_caltable=[],box_width=3,previous_image=previous_image,previous_model=previous_model,do_solarqu_cor=do_solarqu_cor,\
					polmodel_threshold=polmodel_threshold,quvcor_stokes=quvcor_stokes,weight=weight,robust=robust)
			elif inputs.maskfile=='' and inputs.maskstr=='' and inputs.want_auto_masking==False: # If no mask is given and auto maksing is off, use default central mask
				output_PSC=PSC.polselfcal_iteration(num_iter,rms_list,mask_str,start_sigma,maskfile,antenna_to_use,startmodel,startmask,want_auto_masking=inputs.want_auto_masking,\
					stokes=stokes,interactive=interactive,use_ankflagger=inputs.use_ankflagger,do_flag=do_flag,poldistortion_correction=do_poldist,poldistortion_type=poldistortion_type,\
					poldistortion_matrix='UH',calibrator_caltable=[],box_width=3,previous_image=previous_image,previous_model=previous_model,do_solarqu_cor=do_solarqu_cor,\
					polmodel_threshold=polmodel_threshold,quvcor_stokes=quvcor_stokes,weight=weight,robust=robust)
			elif inputs.want_auto_masking==True:
				maskregion=mask_str
				output_PSC=PSC.polselfcal_iteration(num_iter,rms_list,'',start_sigma,maskfile,antenna_to_use,startmodel,startmask,want_auto_masking=inputs.want_auto_masking,\
					stokes=stokes,interactive=interactive,use_ankflagger=inputs.use_ankflagger,do_flag=do_flag,poldistortion_correction=do_poldist,poldistortion_type=poldistortion_type,\
					poldistortion_matrix='UH',calibrator_caltable=[],box_width=3,previous_image=previous_image,previous_model=previous_model,do_solarqu_cor=do_solarqu_cor,\
					maskregion=maskregion,polmodel_threshold=polmodel_threshold,quvcor_stokes=quvcor_stokes,weight=weight,robust=robust)
			
			if num_iter>min(min_iteration,10) and do_solarqu_cor==True:
				do_solarqu_cor=False
				quvcor_stokes='QU'
			
			if type(output_PSC)==tuple:				
				msg_code,out_dict,negative_dyn_range=output_PSC
			else:
				msg_code=output_PSC
				
			logger.info('do_poldist = '+str(do_poldist)+'\n')
			logger.info('do_solarqu_cor = '+str(do_solarqu_cor)+'\n')
			logger.info('perform_gaincal = '+str(perform_gaincal)+'\n')
			if perform_gaincal==True:
				logger.info('do_gaincal = '+str(do_gaincal)+'\n')
			if do_poldist==True:			
				logger.info('poldistortion_type = '+str(poldistortion_type)+'\n')
			logger.info('scratch = '+str(scratch)+'\n')

			if 'ref' in msname:		
				PSC.file_remover_and_keeper(num_iter,msg_code,ref_time_chan=True)  # Removing files and keeping the required ones
			else:
				PSC.file_remover_and_keeper(num_iter,msg_code,ref_time_chan=False)

			if type(output_PSC)!=tuple: # If some error occured only error message will be returned
				if scratch==True:
					if verbose==False:
						print (error_msgs(msg_code))
					logger.error(error_msgs(msg_code))
					end_selfcal=True
					if 'ref' in msname:
						os.chdir(cwd)
						if __name__!='__main__':
							touch_file=basedir+'/.Finished_pcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_'+str(msg_code+100)
							msg_str='Dear PAIRCARS user,\n\nPolarisation self-calibration for : '+os.path.basename(msname)+'\nMessage : '+error_msgs(100)+', '+error_msgs(msg_code)\
											+'\n\nBest regards,\nPAIRCARS developing team'
							msg_subject='Notification from PAIRCARS : Polarisation Selfcal : OBSID = '+str(OBSID)
							if inputs.send_notification==True:
								send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[])
							os.system('touch '+touch_file)
							os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_*.ms')
							end_time=time.time()
							run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
							logger.info('Total runtime : '+str(run_time))
							os.system('cp -r '+working_dir+'/Pol_Selfcal.log '+basedir+'/logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
							if inputs.keep_logger and verbose==True:
								os.system('cp -r '+working_dir+'/Pol_Selfcal_verbose.log '+basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
							os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice* '+working_dir+'/I*image')
						return msg_code+100
					else:
						os.chdir(cwd)
						if __name__!='__main__':
							touch_file=basedir+'/.Finished_pcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_'+str(msg_code)
							msg_str='Dear PAIRCARS user,\n\nPolarisation self-calibration for : '+\
									os.path.basename(msname)+'\nMessage : '+error_msgs(6)+'\n\nBest regards,\nPAIRCARS developing team'
							msg_subject='Notification from PAIRCARS : Polarisation Selfcal : OBSID = '+str(OBSID)
							if inputs.send_notification==True:
								send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[])
							os.system('touch '+touch_file)
							os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_*.ms')
							end_time=time.time()
							run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
							logger.info('Total runtime : '+str(run_time))
							os.system('cp -r '+working_dir+'/Pol_Selfcal.log '+basedir+'/logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
							if inputs.keep_logger and verbose==True:
								os.system('cp -r '+working_dir+'/Pol_Selfcal_verbose.log '+basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
							os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice* '+working_dir+'/I*image')
						return msg_code
				else:
					scratch=True
					perform_gaincal=True
					if verbose==False:
						print ('######################\nGoing for a full selfcal from scratch = True\n######################\n')
					logger.info('######################\n')
					logger.info('Going for a full selfcal from scratch = True\n')
					logger.info('######################\n')
					os.system('rm -rf junk1*')
					do_selfcal=False
					break
			else:  # If no error occured continuing pol selfcal
				DR_neg=negative_dyn_range
				result_I=out_dict['I'] # Stokes I results
				result_Q=out_dict['Q'] # Stokes Q results
				result_U=out_dict['U'] # Stokes U results
				result_V=out_dict['V'] # Stokes V results
				dyn_I=result_I[0]
				rms_I=result_I[1]
				flux_I=result_I[2]
				dyn_Q=result_Q[0]
				rms_Q=result_Q[1]
				flux_Q=result_Q[2]
				dyn_U=result_U[0]
				rms_U=result_U[1]
				flux_U=result_U[2]
				dyn_V=result_V[0]
				rms_V=result_V[1]
				flux_V=result_V[2]	
				rms_list=[rms_I,rms_Q,rms_U,rms_V]
				flux_T=np.sqrt(flux_Q**2+flux_U**2+flux_V**2)
				flux_P=np.sqrt(flux_Q**2+flux_U**2)	
			
				if num_iter==0:
					DR5=DR3=DR1=dyn_I
					DR6=DR4=DR2=DR_neg
					FX3_I=FX2_I=FX1_I=flux_I/rms_I
					FX3_Q=FX2_Q=FX1_Q=flux_Q/rms_Q
					FX3_U=FX2_U=FX1_U=flux_U/rms_U
					FX3_V=FX2_V=FX1_V=flux_V/rms_V
					FX3_T=FX2_T=FX1_T=flux_T/np.sqrt(rms_Q**2+rms_U**2+rms_V**2)
					FX3_P=FX2_P=FX1_P=flux_P/np.sqrt(rms_Q**2+rms_U**2)
					if os.path.isfile('IMGSTAT_pol.npy')==False:
						PSC.IMSTAT_record(DR5,DR6,FX3_I,FX3_Q,FX3_U,FX3_V,FX3_T,FX3_P,'IMGSTAT_pol',init=True) # Keeping image statistics
					else:
						PSC.IMSTAT_record(DR5,DR6,FX3_I,FX3_Q,FX3_U,FX3_V,FX3_T,FX3_P,'IMGSTAT_pol',init=False) # Keeping image statistics
				elif num_iter==1:
					DR5=dyn_I
					DR6=DR_neg
					FX3_I=flux_I/rms_I
					FX3_Q=flux_Q/rms_Q
					FX3_U=flux_U/rms_U
					FX3_V=flux_V/rms_V
					FX3_T=flux_T/np.sqrt(rms_Q**2+rms_U**2+rms_V**2)
					FX3_P=flux_P/np.sqrt(rms_Q**2+rms_U**2)
					if os.path.isfile('IMGSTAT_pol.npy')==False:
						PSC.IMSTAT_record(DR5,DR6,FX3_I,FX3_Q,FX3_U,FX3_V,FX3_T,FX3_P,'IMGSTAT_pol',init=True) # Keeping image statistics
					else:
						PSC.IMSTAT_record(DR5,DR6,FX3_I,FX3_Q,FX3_U,FX3_V,FX3_T,FX3_P,'IMGSTAT_pol',init=False) # Keeping image statistics
				else:
					DR1=DR3
					DR3=DR5
					DR5=dyn_I
					DR2=DR4
					DR4=DR6
					DR6=DR_neg			
					FX1_I=FX2_I
					FX2_I=FX3_I
					FX3_I=flux_I/rms_I
					FX1_Q=FX2_Q
					FX2_Q=FX3_Q
					FX3_Q=flux_Q/rms_Q
					FX1_U=FX2_U
					FX2_U=FX3_U
					FX3_U=flux_U/rms_U
					FX1_V=FX2_V
					FX2_V=FX3_V
					FX3_V=flux_V/rms_V
					FX1_T=FX2_T
					FX2_T=FX3_T
					FX3_T=flux_T/np.sqrt(rms_Q**2+rms_U**2+rms_V**2)
					FX1_P=FX2_P
					FX2_P=FX3_P
					FX3_P=flux_P/np.sqrt(rms_Q**2+rms_U**2)
					if os.path.isfile('IMGSTAT_pol.npy')==False:
						PSC.IMSTAT_record(DR5,DR6,FX3_I,FX3_Q,FX3_U,FX3_V,FX3_T,FX3_P,'IMGSTAT_pol',init=True) # Keeping image statistics
					else:
						PSC.IMSTAT_record(DR5,DR6,FX3_I,FX3_Q,FX3_U,FX3_V,FX3_T,FX3_P,'IMGSTAT_pol',init=False) # Keeping image statistics
				
				if do_solarqu_cor==True and done_qucor==False:
					done_qucor=True

				if done_qucor==True:
					previous_image='junk1.image'
					previous_model='junk1.model'

				if os.path.isfile('Pol_selfcal_record.npy'):
					os.system('rm -rf Pol_selfcal_record.npy')
				selfcal_record=np.array([num_iter,DR1,DR3,DR5,DR2,DR4,DR6,FX3_I,FX3_Q,FX3_U,FX3_V,FX3_T,FX3_P,FX2_I,FX2_Q,FX2_U,FX2_V,FX2_T,FX2_P,FX1_I,FX1_Q,FX1_U,FX1_V,FX1_T,FX1_P,\
								rms_list,scratch,start_sigma,num_iteration_after_poldist,num_iter_after_qucor,\
								num_iter_fixed_sigma,startmodel,startmask,uvsub_flag_count,do_solarqu_cor,do_poldist,\
								do_pbcor,do_gaincal,gaincal_count,done_qucor,pre_res],dtype='object')
				np.save('Pol_selfcal_record',selfcal_record)

				if verbose==False:
					print ('RMS based dynamic ranges : \n'+str(DR1)+','+str(DR3)+','+str(DR5)+'\n')
					print ('Negative based dynamic ranges : \n'+str(DR2)+','+str(DR4)+','+str(DR6)+'\n')
					print ('Total flux based dynamic ranges :.\n')
					print ('Stokes I : '+str(FX1_I)+', '+str(FX2_I)+', '+str(FX3_I)+'\n')
					print ('Stokes Q : '+str(FX1_Q)+', '+str(FX2_Q)+', '+str(FX3_Q)+'\n')
					print ('Stokes U : '+str(FX1_U)+', '+str(FX2_U)+', '+str(FX3_U)+'\n')
					print ('Stokes V : '+str(FX1_V)+', '+str(FX2_V)+', '+str(FX3_V)+'\n')
					print ('Stokes T : '+str(FX1_T)+', '+str(FX2_T)+', '+str(FX3_T)+'\n')
					print ('Stokes P : '+str(FX1_P)+', '+str(FX2_P)+', '+str(FX3_P)+'\n')

				logger.info('RMS based dynamic ranges:\n')
				logger.info(str(DR1)+','+str(DR3)+','+str(DR5)+'\n')
				logger.info('Negative based dynamic ranges:\n')
				logger.info(str(DR2)+','+str(DR4)+','+str(DR6)+'\n')
				logger.info('Total flux based dynamic ranges :.\n')
				logger.info('Stokes I : '+str(FX1_I)+', '+str(FX2_I)+', '+str(FX3_I)+'\n')
				logger.info('Stokes Q : '+str(FX1_Q)+', '+str(FX2_Q)+', '+str(FX3_Q)+'\n')
				logger.info('Stokes U : '+str(FX1_U)+', '+str(FX2_U)+', '+str(FX3_U)+'\n')
				logger.info('Stokes V : '+str(FX1_V)+', '+str(FX2_V)+', '+str(FX3_V)+'\n')
				logger.info('Stokes T : '+str(FX1_T)+', '+str(FX2_T)+', '+str(FX3_T)+'\n')
				logger.info('Stokes P : '+str(FX1_P)+', '+str(FX2_P)+', '+str(FX3_P)+'\n')

				if do_solarqu_cor==True:
					do_solarqu_cor=False
					
				if done_qucor==True and gaincal_count<1:
					if perform_gaincal==True:
						do_gaincal=True
						if verbose==False:
							print ('Going for a gaincal with leakage corrected models.\n###########################\n')
						logger.info('Going for a gaincal with leakage corrected models.\n')
						logger.info('###########################\n')
					else:
						do_gaincal=False
						if done_qucor==False:
							done_qucor=True
						if verbose==False:
							print ('Not going for a gaincal with leakage corrected models.\n###########################\n')
						logger.info('Not going for a gaincal with leakage corrected models.\n')
						logger.info('###########################\n')
					do_selfcal=False
					os.system('rm -rf junk1.image junk1.model junk1.ms')
					break
				
				############## If statement 1 (DR decrease)

				if (((DR5<0.85*DR3 and DR5<0.9*DR1 and DR3>DR1) and (DR6<0.85*DR4 and DR6<0.9*DR2 and DR4>DR2))or ((DR5<0.9*DR3 and DR1>1.5*DR3) and (DR6<0.9*DR4 and DR2>1.5*DR4))\
					 and ((num_iteration_after_poldist>min_iteration and done_qucor==False)	or (num_iter_after_qucor>min_iteration and done_qucor==True))):
					# If DR decreases.
					# Case 1: If DR decreases less than 90% and 85% of previous two rounds and all antennas are added. This is a check if the rms is diverging. 
					# Case 2: If DR decreases less than 90% of previous round but DR increases more than 1.5 times in last two rounds and no new antennas are added
					if uvsub_flag_count<1 and want_uvsub_flag==True:
						DR3=DR1
						DR4=DR2
						flaglist=flagmanager(vis=msname,mode='list')
						logger.info('Present flag versions : \n')
						logger.info(str(flaglist)+'\n')
						flaglist_keys=list(flaglist.keys())
						flaglist_keys.remove('MS')
						for key in flaglist_keys:	
							flagversion=flaglist[key]['name']
							logger.info('flagmanager(vis=\''+msname+'\',mode=\'delete\',versionname=\''+flagversion+'\')')
							flagmanager(vis=msname,mode='delete',versionname=flagversion)
						if inputs.use_ankflagger:
							os.system('cp -r '+msname+' '+msname+'.backup')
							try:
								logger.info('Performing uvsub flagging using aNKflagger due to DR decrease.\n')
								logger.info('do_uvsub_ankflag(\''+msname+'\',model=\'junk0.model\',nthread=1,verbose='+str(verbose)+',flagbackup=False)\n')
								do_uvsub_ankflag(msname,model='junk0.model',nthread=1,verbose=verbose,flagbackup=False)
								os.system('rm -rf '+msname+'.backup')
							except Exception as e:
								os.system('rm -rf '+msname)
								os.system('mv '+msname+'.backup '+msname)
								logger.error('Error in aNKflagger : '+str(e)+'\n')
								logger.info('Error in running aNKflagger. Using rms threshold flagging.\n')
								logger.info('Performing uvsub flagging due to DR decrease.\n')
								logger.info('do_uvsub_flagger(\''+msname+'\',model=\'junk0.model\',mode=\'uvsub_flag\',rmsthresh=[15,10,8],flagbackup=False)\n')
								do_uvsub_flagger(msname,model='junk0.model',mode='uvsub_flag',rmsthresh=[15,10,8],flagbackup=False)
						else:
							logger.info('Performing uvsub flagging due to DR decrease.\n')
							logger.info('do_uvsub_flagger(\''+msname+'\',model=\'junk0.model\',mode=\'uvsub_flag\',rmsthresh=[15,10,8],flagbackup=False)\n')
							do_uvsub_flagger(msname,model='junk0.model',mode='uvsub_flag',rmsthresh=[15,10,8],flagbackup=False)
						uvsub_flag_count+=1
						os.system('rm -rf junk1.model')
						os.system('cp -r junk0.model junk1.model')
						continue
					if scratch==True:
						if (num_iter_after_qucor>min_iteration):
							# Image based Stokes I to Q,U correction is necessary. If it is not done go for it.
							# If scratch is False then it has failed in spite of a good starting point and QU correction. Hence relaxation time is not needed.
							if DR5>min_DR: # Only considered the rms based DR here
								os.system('rm -rf beamcor_backup.ms')
								end_selfcal=True
								logger.info('----------------------------------\n')
								logger.info('Making final calibration table.\n')
								logger.info('delmod(vis=\''+working_dir+'/Backup_beamcorrected.ms\',scr=True,otf=True)\n')
								delmod(vis=working_dir+'/Backup_beamcorrected.ms',scr=True,otf=True)
								logger.info('ft(vis=\''+working_dir+'/Backup_beamcorrected.ms\',model=\'junk0.model\',usescratch=True)\n')
								ft(vis=working_dir+'/Backup_beamcorrected.ms',model=working_dir+'/junk0.model',usescratch=True)
								calib_uvrange_min=IB.calc_calib_uvrange(12)[1]
								calib_uvrange_max=IB.calc_calib_uvrange(12)[2]
								logger.info('cal.calibrate(msname=\''+working_dir+'/Backup_beamcorrected.ms\',caltable=\''+working_dir+'/junk_final.bin\',minuv=\''+\
										str(calib_uvrange_min)+',maxuv=\''+str(calib_uvrange_max)+',j=1,absmem=1,solmode=\'R\',rmsthresh=[15,10,8],quiet='+str(verbose)+')\n')
								cal.calibrate(msname=working_dir+'/Backup_beamcorrected.ms',caltable=working_dir+'/junk_final.bin',minuv=calib_uvrange_min,\
									maxuv=calib_uvrange_max,j=1,absmem=1,solmode='R',rmsthresh=[15,10,8],quiet=verbose)
								logger.info('cal.applycal(msname=\''+working_dir+'/Backup_beamcorrected.ms\',gaintable=\''+working_dir+'/junk_final.bin\',applymode=\'calflag\')\n')
								cal.applycal(msname=working_dir+'/Backup_beamcorrected.ms',gaintable=working_dir+'/junk_final.bin',applymode='calflag')
								os.system('cp -r '+working_dir+'/junk_final.bin '+basedir+'/polcaltables/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.bin')
								os.system('cp -r '+working_dir+'/junk0.model '+basedir+'/polimagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.model')
								os.system('cp -r '+working_dir+'/junk0.image '+basedir+'/polimagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.image')
								os.system('rm -rf Backup_*.ms*')
								if inputs.send_notification==True:
									quickimage=get_quicklook_image(working_dir+'/junk1.image','quick_image_freq_'+freqstr+'_time_'+datestrfile+'.png',\
									freqstr,datestr,DR5,DR6,field_of_view=2)
								os.system('rm -rf '+working_dir+'/junk*')
								if 'ref' in msname:
									os.chdir(cwd)
									if __name__!='__main__':
										touch_file=basedir+'/.Finished_pcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_108'
										msg_str='Dear PAIRCARS user,\n\nPolarisation self-calibration for : '+os.path.basename(msname)+'\nMessage : '+error_msgs(100)+', '+error_msgs(8)\
											+'\n\nBest regards,\nPAIRCARS developing team'
										msg_subject='Notification from PAIRCARS : Polarisation Selfcal : OBSID = '+str(OBSID)
										if inputs.send_notification==True:
											send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
											os.system('rm -rf '+quickimage)
										os.system('touch '+touch_file)
										os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_*.ms')
										end_time=time.time()
										run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
										logger.info('Total runtime : '+str(run_time))
										os.system('cp -r '+working_dir+'/Pol_Selfcal.log '+basedir+'/logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
										if inputs.keep_logger and verbose==True:
											os.system('cp -r '+working_dir+'/Pol_Selfcal_verbose.log '+basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
										os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice* '+working_dir+'/I*image')
									return 108
								else:
									os.chdir(cwd)
									if __name__!='__main__':
										touch_file=basedir+'/.Finished_pcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_8'
										msg_str='Dear PAIRCARS user,\n\nPolarisation self-calibration for : '+\
											os.path.basename(msname)+'\nMessage : '+error_msgs(8)+'\n\nBest regards,\nPAIRCARS developing team'
										msg_subject='Notification from PAIRCARS : Polarisation Selfcal : OBSID = '+str(OBSID)
										if inputs.send_notification==True:
											send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
											os.system('rm -rf '+quickimage)
										os.system('touch '+touch_file)
										os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_*.ms')
										end_time=time.time()
										run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
										logger.info('Total runtime : '+str(run_time))
										os.system('cp -r '+working_dir+'/Pol_Selfcal.log '+basedir+'/logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
										if inputs.keep_logger and verbose==True:
											os.system('cp -r '+working_dir+'/Pol_Selfcal_verbose.log '+basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
										os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice* '+working_dir+'/I*image')
									return 8
						else:
							do_solarqu_cor=True
							if verbose==False:
								print ('####################\nGoing for a image based Stokes I to Q,U leakage correction because Stokes I max DR reached.\n####################\n')
							logger.info('####################\n')
							logger.info('Going for a image based Stokes I to Q,U leakage correction because Stokes I max DR reached.\n')
							logger.info('####################\n')
							continue
					else:
						scratch=True
						perform_gaincal=True
						if verbose==False:
							print ('######################\nGoing for a full selfcal from scratch = True\n #####################################\n')
						logger.info('######################\n')
						logger.info('Going for a full selfcal from scratch = True\n')
						logger.info('######################\n')
						os.system('rm -rf '+working_dir+'/junk*')
						do_selfcal=False
						break

				#######################################################

				# If statement 2 (Exiting selfcal conditions)
				if ((DR5>=inputs.max_DR and abs(DR5/DR3-1)<frac_flux_change and abs(DR5/DR1-1)<frac_flux_change and abs(FX1_Q/FX1_I-FX3_Q/FX3_I)<pol_frac_change \
					and abs(FX2_Q/FX2_I-FX1_Q/FX1_I)<pol_frac_change and abs(FX1_U/FX1_I-FX3_U/FX3_I)<pol_frac_change and abs(FX2_U/FX2_I-FX1_U/FX1_I)<pol_frac_change\
					 and abs(FX1_V/FX1_I-FX3_V/FX3_I)<pol_frac_change and abs(FX2_V/FX2_I-FX1_V/FX1_I)<pol_frac_change) and\
					 (num_iteration_after_poldist>min_iteration or num_iter_after_qucor>min_iteration)):
					# Stokes I DR reached maximum limit and polarised flux converged
					if gaincal_count==1 and done_qucor==True: # If QU correction has been done and new gaincal using leakage is done.
						if verbose==False:
							print ('Reached limiting dynamic range and polarised flux converged\n')
						logger.info('Reached limiting dynamic range and polarised flux converged\n')
						os.system('rm -rf beamcor_backup.ms')
						end_selfcal=True
						logger.info('----------------------------------\n')
						logger.info('Making final calibration table.\n')
						logger.info('delmod(vis=\''+working_dir+'/Backup_beamcorrected.ms\',scr=True,otf=True)\n')
						delmod(vis=working_dir+'/Backup_beamcorrected.ms',scr=True,otf=True)
						logger.info('ft(vis=\''+working_dir+'/Backup_beamcorrected.ms\',model=\'junk1.model\',usescratch=True)\n')
						ft(vis=working_dir+'/Backup_beamcorrected.ms',model=working_dir+'/junk1.model',usescratch=True)
						calib_uvrange_min=IB.calc_calib_uvrange(12)[1]
						calib_uvrange_max=IB.calc_calib_uvrange(12)[2]
						logger.info('cal.calibrate(msname=\''+working_dir+'/Backup_beamcorrected.ms\',caltable=\''+working_dir+'/junk_final.bin\',minuv=\''+\
								str(calib_uvrange_min)+',maxuv=\''+str(calib_uvrange_max)+',j=1,absmem=1,solmode=\'R\',rmsthresh=[15,10,8],quiet='+str(verbose)+')\n')
						cal.calibrate(msname=working_dir+'/Backup_beamcorrected.ms',caltable=working_dir+'/junk_final.bin',minuv=calib_uvrange_min,\
							maxuv=calib_uvrange_max,j=1,absmem=1,solmode='R',rmsthresh=[15,10,8],quiet=verbose)
						logger.info('cal.applycal(msname=\''+working_dir+'/Backup_beamcorrected.ms\',gaintable=\''+working_dir+'/junk_final.bin\',applymode=\'calflag\')\n')
						cal.applycal(msname=working_dir+'/Backup_beamcorrected.ms',gaintable=working_dir+'/junk_final.bin',applymode='calflag')
						os.system('cp -r '+working_dir+'/junk_final.bin '+basedir+'/polcaltables/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.bin')
						os.system('cp -r '+working_dir+'/junk1.model '+basedir+'/polimagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.model')
						os.system('cp -r '+working_dir+'/junk1.image '+basedir+'/polimagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.image')
						os.system('rm -rf Backup_*.ms*')
						if inputs.send_notification==True:
							quickimage=get_quicklook_image(working_dir+'/junk1.image','quick_image_freq_'+freqstr+'_time_'+datestrfile+'.png',\
							freqstr,datestr,DR5,DR6,field_of_view=2)
						os.system('rm -rf '+working_dir+'/junk*')
						if __name__!='__main__':
							touch_file=basedir+'/.Finished_pcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_0'
							msg_str='Dear PAIRCARS user,\n\nPolarisation self-calibration for : '+\
								os.path.basename(msname)+'\nMessage : '+error_msgs(0)+'\n\nBest regards,\nPAIRCARS developing team'
							msg_subject='Notification from PAIRCARS : Polarisation Selfcal : OBSID = '+str(OBSID)
							if inputs.send_notification==True:
								send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
								os.system('rm -rf '+quickimage)
							os.system('touch '+touch_file)
							os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_*.ms')
							end_time=time.time()
							run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
							logger.info('Total runtime : '+str(run_time))
							os.system('cp -r '+working_dir+'/Pol_Selfcal.log '+basedir+'/logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
							if inputs.keep_logger and verbose==True:
								os.system('cp -r '+working_dir+'/Pol_Selfcal_verbose.log '+basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
							os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice* '+working_dir+'/I*image')
						return 0
					elif gaincal_count<1 and done_qucor==False:
						if num_iter_after_qucor<1 and do_solarqu_cor==False:
							do_solarqu_cor=True
							if verbose==False:
								print ('####################\nGoing for a image based Stokes I to Q,U leakage correction because Stokes I max DR reached and polarised flux '+\
										'converged.\n####################\n')
							logger.info('####################\n')
							logger.info('Going for a image based Stokes I to Q,U leakage correction because Stokes I max DR reached and polarised flux converged.\n')
							logger.info('####################\n')
							continue				
				elif (abs(DR5-DR3)<DR_delta_rms and abs(DR5-DR1)<DR_delta_rms and abs(FX1_Q/FX1_I-FX3_Q/FX3_I)<pol_frac_change and abs(FX2_Q/FX2_I-FX1_Q/FX1_I)<pol_frac_change and \
					abs(FX1_U/FX1_I-FX3_U/FX3_I)<pol_frac_change and abs(FX2_U/FX2_I-FX1_U/FX1_I)<pol_frac_change and abs(FX1_V/FX1_I-FX3_V/FX3_I)<pol_frac_change and \
					abs(FX2_V/FX2_I-FX1_V/FX1_I)<pol_frac_change): # If polarised flux converged
					if num_iter_fixed_sigma>min_num_iter_fixed_sigma and (num_iteration_after_poldist>min_iteration or num_iter_after_qucor>max(min_iteration,5)):
						if gaincal_count<1 and done_qucor==False: # If QU correction and leakage corrected gaincal not done
							sigma,pre_res=PSC.reduce_sigma('junk1.image',start_sigma,inputs.sigma_step,inputs.min_sigma,pre_residual=pre_res,residual_frac=frac_flux_change,\
									stokes_list=['I','Q','U','V'])
							if sigma<start_sigma:						
								start_sigma=sigma	
								num_iter_fixed_sigma=0
							else:
								if num_iter_after_qucor<1 and do_solarqu_cor==False:
									do_solarqu_cor=True
									if verbose==False:
										print ('####################\nGoing for a image based Stokes I to Q,U leakage correction because '+\
												'polarised flux converged.\n####################\n')
									logger.info('####################\n')
									logger.info('Going for a image based Stokes I to Q,U leakage correction because polarised flux converged.\n')
									logger.info('####################\n')
									continue
						else:
							if verbose==False:
								print ('Selfcal converged. Residual flux inside the mask is less than :'+str(frac_flux_change*100)+'%. Stopped sigma :'+str(start_sigma)+'\n') 	
							logger.info('Selfcal converged. Residual flux inside the mask is less than :'+str(frac_flux_change*100)+'%. Stopped sigma :'+str(start_sigma)+'\n')	
							os.system('rm -rf beamcor_backup.ms')
							end_selfcal=True
							logger.info('----------------------------------\n')
							logger.info('Making final calibration table.\n')
							logger.info('delmod(vis=\''+working_dir+'/Backup_beamcorrected.ms\',scr=True,otf=True)\n')
							delmod(vis=working_dir+'/Backup_beamcorrected.ms',scr=True,otf=True)
							logger.info('ft(vis=\''+working_dir+'/Backup_beamcorrected.ms\',model=\'junk1.model\',usescratch=True)\n')
							ft(vis=working_dir+'/Backup_beamcorrected.ms',model=working_dir+'/junk1.model',usescratch=True)
							calib_uvrange_min=IB.calc_calib_uvrange(12)[1]
							calib_uvrange_max=IB.calc_calib_uvrange(12)[2]
							logger.info('cal.calibrate(msname=\''+working_dir+'/Backup_beamcorrected.ms\',caltable=\''+working_dir+'/junk_final.bin\',minuv=\''+\
									str(calib_uvrange_min)+',maxuv=\''+str(calib_uvrange_max)+',j=1,absmem=1,solmode=\'R\',rmsthresh=[15,10,8],quiet='+str(verbose)+')\n')
							cal.calibrate(msname=working_dir+'/Backup_beamcorrected.ms',caltable=working_dir+'/junk_final.bin',minuv=calib_uvrange_min,\
								maxuv=calib_uvrange_max,j=1,absmem=1,solmode='R',rmsthresh=[15,10,8],quiet=verbose)
							logger.info('cal.applycal(msname=\''+working_dir+'/Backup_beamcorrected.ms\',gaintable=\''+working_dir+'/junk_final.bin\',applymode=\'calflag\')\n')
							cal.applycal(msname=working_dir+'/Backup_beamcorrected.ms',gaintable=working_dir+'/junk_final.bin',applymode='calflag')
							os.system('cp -r '+working_dir+'/junk_final.bin '+basedir+'/polcaltables/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.bin')
							os.system('cp -r '+working_dir+'/junk1.model '+basedir+'/polimagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.model')
							os.system('cp -r '+working_dir+'/junk1.image '+basedir+'/polimagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.image')
							os.system('rm -rf Backup_*.ms*')
							if inputs.send_notification==True:
								quickimage=get_quicklook_image(working_dir+'/junk1.image','quick_image_freq_'+freqstr+'_time_'+datestrfile+'.png',\
								freqstr,datestr,DR5,DR6,field_of_view=2)
							os.system('rm -rf '+working_dir+'/junk*')
							if __name__!='__main__':
								touch_file=basedir+'/.Finished_pcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_0'
								msg_str='Dear PAIRCARS user,\n\nPolarisation self-calibration for : '+\
									os.path.basename(msname)+'\nMessage : '+error_msgs(0)+'\n\nBest regards,\nPAIRCARS developing team'
								msg_subject='Notification from PAIRCARS : Polarisation Selfcal : OBSID = '+str(OBSID)
								if inputs.send_notification==True:
									send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
									os.system('rm -rf '+quickimage)
								os.system('touch '+touch_file)
								os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_*.ms')
								end_time=time.time()
								run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
								logger.info('Total runtime : '+str(run_time))
								os.system('cp -r '+working_dir+'/Pol_Selfcal.log '+basedir+'/logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
								if inputs.keep_logger and verbose==True:
									os.system('cp -r '+working_dir+'/Pol_Selfcal_verbose.log '+basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
								os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice* '+working_dir+'/I*image')
							return 0
								
									
				#############################################################			
				# If statement 3 (Using last round model) 
				#(If DR increases at least DR_delta and all antennas are added and number of iteration at fixed antenna is greater than 5)
				
				if (abs(DR5/DR3-1)<frac_flux_change and abs(DR5/DR1-1)<frac_flux_change and abs(FX1_Q/FX1_I-FX3_Q/FX3_I)<pol_frac_change and abs(FX2_Q/FX2_I-FX1_Q/FX1_I)<pol_frac_change\
					 and abs(FX1_U/FX1_I-FX3_U/FX3_I)<pol_frac_change and abs(FX2_U/FX2_I-FX1_U/FX1_I)<pol_frac_change and abs(FX1_V/FX1_I-FX3_V/FX3_I)<pol_frac_change\
					 and abs(FX2_V/FX2_I-FX1_V/FX1_I)<pol_frac_change and (num_iteration_after_poldist>min_iteration or num_iter_after_qucor>min_iteration)):
					startmodel='junk1.model'
				else:
					startmodel=''
				if num_iter>1:
					startmask='junk1.mask'
				else:
					startmask=''
				num_iter+=1
				num_iter_fixed_sigma+=1
				if do_poldist==True:
					num_iteration_after_poldist+=1
				if done_qucor==True:
					num_iter_after_qucor+=1

				###############################################################
				# If statement 4 (Reached maximum selfcal rounds)
				if (num_iter>max_iteration and num_iter_after_qucor>1) or (num_iter>int(max_iteration/3) and num_iter_after_qucor<1 and do_solarqu_cor==False): 
					if scratch==True:
						if (gaincal_count==1):
							if DR5>min_DR:
								os.system('rm -rf beamcor_backup.ms')
								end_selfcal=True
								logger.info('----------------------------------\n')
								logger.info('Making final calibration table.\n')
								logger.info('delmod(vis=\''+working_dir+'/Backup_beamcorrected.ms\',scr=True,otf=True)\n')
								delmod(vis=working_dir+'/Backup_beamcorrected.ms',scr=True,otf=True)
								logger.info('ft(vis=\''+working_dir+'/Backup_beamcorrected.ms\',model=\'junk1.model\',usescratch=True)\n')
								ft(vis=working_dir+'/Backup_beamcorrected.ms',model=working_dir+'/junk1.model',usescratch=True)
								calib_uvrange_min=IB.calc_calib_uvrange(12)[1]
								calib_uvrange_max=IB.calc_calib_uvrange(12)[2]
								logger.info('cal.calibrate(msname=\''+working_dir+'/Backup_beamcorrected.ms\',caltable=\''+working_dir+'/junk_final.bin\',minuv=\''+\
										str(calib_uvrange_min)+',maxuv=\''+str(calib_uvrange_max)+',j=1,absmem=1,solmode=\'R\',rmsthresh=[15,10,8],quiet='+str(verbose)+')\n')
								cal.calibrate(msname=working_dir+'/Backup_beamcorrected.ms',caltable=working_dir+'/junk_final.bin',minuv=calib_uvrange_min,\
									maxuv=calib_uvrange_max,j=1,absmem=1,solmode='R',rmsthresh=[15,10,8],quiet=verbose)
								logger.info('cal.applycal(msname=\''+working_dir+'/Backup_beamcorrected.ms\',gaintable=\''+working_dir+'/junk_final.bin\',applymode=\'calflag\')\n')
								cal.applycal(msname=working_dir+'/Backup_beamcorrected.ms',gaintable=working_dir+'/junk_final.bin',applymode='calflag')
								os.system('cp -r '+working_dir+'/junk_final.bin '+basedir+'/polcaltables/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.bin')
								os.system('cp -r '+working_dir+'/junk1.model '+basedir+'/polimagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.model')
								os.system('cp -r '+working_dir+'/junk1.image '+basedir+'/polimagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.image')
								os.system('rm -rf Backup_*.ms*')
								if inputs.send_notification==True:
									quickimage=get_quicklook_image(working_dir+'/junk1.image','quick_image_freq_'+freqstr+'_time_'+datestrfile+'.png',\
									freqstr,datestr,DR5,DR6,field_of_view=2)
								os.system('rm -rf '+working_dir+'/junk*')
								if 'ref' in msname:
									os.chdir(cwd)
									if __name__!='__main__':
										touch_file=basedir+'/.Finished_pcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_109'
										msg_str='Dear PAIRCARS User,\n\nPolarisation self-calibration for : '+os.path.basename(msname)+'\nMessage : '+error_msgs(100)+', '+error_msgs(9)\
													+'\n\nBest regards,\nPAIRCARS developing team'
										msg_subject='Notification from PAIRCARS : Polarisation Selfcal : OBSID = '+str(OBSID)
										if inputs.send_notification==True:
											send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
											os.system('rm -rf '+quickimage)
										os.system('touch '+touch_file)
										os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_*.ms')
										end_time=time.time()
										run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
										logger.info('Total runtime : '+str(run_time))
										os.system('cp -r '+working_dir+'/Pol_Selfcal.log '+basedir+'/logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
										if inputs.keep_logger and verbose==True:
											os.system('cp -r '+working_dir+'/Pol_Selfcal_verbose.log '+basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
										os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice* '+working_dir+'/I*image')
									return 109		
								else:
									os.chdir(cwd)
									if __name__!='__main__':
										touch_file=basedir+'/.Finished_pcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_9'
										msg_str='Dear PAIRCARS user,\n\nPolarisation self-calibration for : '+\
											os.path.basename(msname)+'\nMessage : '+error_msgs(9)+'\n\nBest regards,\nPAIRCARS developing team'
										msg_subject='Notification from PAIRCARS : Polarisation Selfcal : OBSID = '+str(OBSID)
										if inputs.send_notification==True:
											send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
											os.system('rm -rf '+quickimage)
										os.system('touch '+touch_file)
										os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_*.ms')
										end_time=time.time()
										run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
										logger.info('Total runtime : '+str(run_time))
										os.system('cp -r '+working_dir+'/Pol_Selfcal.log '+basedir+'/logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
										if inputs.keep_logger and verbose==True:
											os.system('cp -r '+working_dir+'/Pol_Selfcal_verbose.log '+basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
										os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice* '+working_dir+'/I*image')
									return 9
							else:
								if verbose==False:
									print (error_msgs(13))
								logger.error(error_msgs(13))
								end_selfcal=True
								os.system('rm -rf '+working_dir+'/junk*')
								os.chdir(cwd)
								if __name__!='__main__':
									touch_file=basedir+'/.Finished_pcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_13'
									msg_str='Dear PAIRCARS user,\n\nPolarisation self-calibration for : '+\
										os.path.basename(msname)+'\nMessage : '+error_msgs(13)+'\n\nBest regards,\nPAIRCARS developing team'
									msg_subject='Notification from PAIRCARS : Intensity Selfcal : OBSID = '+str(OBSID)
									if inputs.send_notification==True:
										send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
										os.system('rm -rf '+quickimage)
									os.system('touch '+touch_file)
									os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_*.ms')
									end_time=time.time()
									run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
									logger.info('Total runtime : '+str(run_time))
									os.system('cp -r '+working_dir+'/Pol_Selfcal.log '+basedir+'/logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
									if inputs.keep_logger and verbose==True:
										os.system('cp -r '+working_dir+'/Pol_Selfcal_verbose.log '+basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
									os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice* '+working_dir+'/I*image')
								return 13
						else:
							if num_iter_after_qucor<1 and do_solarqu_cor==False:
								do_solarqu_cor=True
								if verbose==False:
									print ('#################\nGoing for a image based Stokes I to Q,U leakage correction because maximum iterations reached.\n#################\n')
								logger.info('#################\n')
								logger.info('Going for a image based Stokes I to Q,U leakage correction because maximum iterations reached.\n')
								logger.info('#################\n')
								continue
					else:
						scratch=True
						perform_gaincal=True
						if verbose==False:
							print ('######################\nGoing for a full selfcal from scratch = True\n ##########################\n')
						logger.info('######################\n')
						logger.info('Going for a full selfcal from scratch = True\n')
						logger.info('######################\n')
						os.system('rm -rf '+working_dir+'/junk*')
						do_selfcal=False
						break
					

# Function to run the script stand alone from command line
if __name__=='__main__':
	start_time=time.time()
	usage= ' Perform polarisation self calibration of a single time and frequency slice'
	parser = OptionParser(usage=usage)
	parser.add_option('--msname',dest="chantime_msname",default=None,help="Name of measurement set of a single time anf frequency slice",metavar="Measurement Set")
	parser.add_option('--metafits',dest="metafits",default=None,help="Name of metafits file",metavar="Metafits file")
	parser.add_option('--workdir',dest='workdir',default=None,help='Name of the working directory',metavar='Directory path')
	parser.add_option('--verbose',dest="verbose",default=False,help="Verbose mode",metavar="Boolean")
	parser.add_option('--interactive',dest="interactive",default=False,help="Interactive mode",metavar="Boolean")
	parser.add_option('--fresh',dest="fresh",default=True,help="Start fresh self calibration loop",metavar="Boolean")
	parser.add_option('--gaincal',dest="gaincal",default=False,help="Perform gaincal using leakage corrected model (Only do when no calibrator observation is present)",metavar="Boolean")
	parser.add_option('--caltables',dest="caltables",default='',help="Previous caltables",metavar="String, comma separated")
	parser.add_option('--wsclean',dest="use_wsclean",default=True,help="Use WSClean for imaging or not",metavar="Boolean")
	(options, args) = parser.parse_args()

	if (os.path.isfile(str(options.workdir)+'/Pol_Selfcal.log') and eval(str(options.fresh))==True) or \
			(os.path.isfile(str(options.workdir)+'/Pol_Selfcal.log') and os.path.isdir(str(options.workdir)+'/junk1.ms')==False and eval(str(options.fresh))==False):
		os.system('rm -rf '+str(options.workdir)+'/Pol_Selfcal.log')
	if os.path.isfile(str(options.workdir)+'/Pol_Selfcal_verbose.log') and eval(str(options.verbose))==True:
			os.system('rm -rf '+str(options.workdir)+'/Pol_Selfcal_verbose.log')
	formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
	logger = logging.getLogger('pol_selfcal_log')
	logger.setLevel(logging.DEBUG)
	if eval(str(options.verbose))==True:
		console=logging.StreamHandler(sys.stdout)
		console.setFormatter(formatter)
		logger.addHandler(console)
	filehandle=logging.FileHandler(str(options.workdir)+'/Pol_Selfcal.log')
	filehandle.setFormatter(formatter)
	logger.addHandler(filehandle)
	logger.propagate = False

	cwd=os.getcwd()
	sys.path.append(cwd)
	if os.path.isfile(cwd+'/selfcal_inputs.py')==False:
		print ('Input file does not exist.\n')
		os._exit(0)
	else:
		import selfcal_inputs as inputs
		from selfcal_inputs import *

	if options.chantime_msname[-1]=='/':
		options.chantime_msname=options.chantime_msname[:-1]

	msbasename=os.path.basename(options.chantime_msname)
	OBSID=get_OBSID(options.metafits)
	basemsdir=os.path.dirname(options.workdir).split('/')[-1]
	if options.chantime_msname==None or os.path.isdir(options.chantime_msname)==False:
		logger.info('Measurement set does not exist. Exititing...\n')
		touch_file=inputs.basedir+'/.Finished_pcal_'+str(OBSID)+'_'+basemsdir+'_'+msbasename+'_'+str('noms')
		end_time=time.time()
		run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
		logger.info('#############################\n')
		logger.info('Polarisation selfcal failed for ms : '+options.chantime_msname+'\n')
		logger.info('Total runtime : '+str(run_time)+'\n')
		logger.info('##############################\n')
		msg_str='Dear PAIRCARS user,\n\nPolarisatioin self-calibration for : '+msbasename+'\nMessage : No measurement set is present\nTotal runtime : '+\
					str(run_time)+'\n\nBest regards,\nPAIRCARS developing team'
		msg_subject='Notification from PAIRCARS : Polarisation Selfcal : OBSID = '+str(OBSID)
		if inputs.send_notification==True:
			send_paircars_notification(inputs.email,msg_subject,msg_str)
		os.system('touch '+touch_file)
		os.system('rm -rf '+options.workdir+'/TempLattice*')
		file_str=msbasename.split('.ms')[0]
		if os.path.isdir(basedir+'/logs/'+str(OBSID)+'/'+basemsdir)==False:
			os.makedirs(basedir+'/logs/'+str(OBSID)+'/'+basemsdir)
		if os.path.isdir(basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir)==False and inputs.keep_logger==True and eval(str(options.verbose)):
			os.makedirs(basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir)
		os.system('cp -r '+options.workdir+'/Pol_Selfcal.log '+basedir+'/logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
		if inputs.keep_logger and eval(str(options.verbose))==True:
			os.system('cp -r '+options.workdir+'/Pol_Selfcal_verbose.log '+basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
		os.system('rm -rf '+options.workdir+'/'+file_str+'*')
		os.system('rm -rf '+options.workdir+'/*.log')
		os._exit(0)
	
	if options.metafits==None or os.path.isfile(options.metafits)==False:
		logger.info('Metafits file does not exist. Exititing...\n')
		touch_file=inputs.basedir+'/.Finished_pcal_'+str(OBSID)+'_'+basemsdir+'_'+msbasename+'_'+str('nometa')
		end_time=time.time()
		run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
		logger.info('#############################\n')
		logger.info('Gain selfcal failed for ms : '+options.chantime_msname+'\n')
		logger.info('Total runtime : '+str(run_time)+'\n')
		logger.info('##############################\n')
		msg_str='Dear PAIRCARS user,\n\nPolarisation self-calibration for : '+msbasename+'\nMessage : No metafits file is present\nTotal runtime : '+\
					str(run_time)+'\n\nBest regards,\nPAIRCARS developing team'
		msg_subject='Notification from PAIRCARS : Polarisation Selfcal : OBSID = '+str(OBSID)
		if inputs.send_notification==True:
			send_paircars_notification(inputs.email,msg_subject,msg_str)
		os.system('touch '+touch_file)
		os.system('rm -rf '+options.workdir+'/TempLattice*')
		file_str=msbasename.split('.ms')[0]
		if os.path.isdir(basedir+'/logs/'+str(OBSID)+'/'+basemsdir)==False:
			os.makedirs(basedir+'/logs/'+str(OBSID)+'/'+basemsdir)
		if os.path.isdir(basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir)==False and inputs.keep_logger==True and eval(str(options.verbose)):
			os.makedirs(basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir)
		os.system('cp -r '+options.workdir+'/Pol_Selfcal.log '+basedir+'/logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
		if inputs.keep_logger and eval(str(options.verbose))==True:
			os.system('cp -r '+options.workdir+'/Pol_Selfcal_verbose.log '+basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
		os.system('rm -rf '+options.workdir+'/'+file_str+'*')
		os.system('rm -rf '+options.workdir+'/*.log')
		os._exit(0)

	try:
		previous_touch_list=glob.glob(inputs.basedir+'/.Finished_pcal_'+str(OBSID)+'_'+basemsdir+'_'+msbasename+'_*')
		if len(previous_touch_list)!=0:
			os.system('rm -rf '+inputs.basedir+'/.Finished_pcal_'+str(OBSID)+'_'+basemsdir+'_'+msbasename+'_*')
		print ('\n\t##########################\n\tStarting Polarisation self-calibration.....\n\t##########################\n')
		print ('run_pol_selfcal(\''+options.chantime_msname+'\',\''+options.metafits+'\',\''+options.workdir+'\',verbose='+str(options.verbose)+',interactive='+str(options.interactive)+\
				',start_fresh='+str(options.fresh)+',perform_gaincal='+str(options.gaincal)+',caltables=\''+str(options.caltables)+'\',use_wsclean='+str(use_wsclean)+')\n')
		msg=run_pol_selfcal(options.chantime_msname,options.metafits,options.workdir,verbose=eval(str(options.verbose)),interactive=eval(str(options.interactive)),\
				start_fresh=eval(str(options.fresh)),perform_gaincal=eval(str(options.gaincal)),caltables=str(options.caltables),use_wsclean=eval(str(use_wsclean)))
		if type(msg)==int:
			if msg>100:
				msg1=msg-100
				msg_str='Message : '+error_msgs(100)+', '+error_msgs(msg1)+'\n'
				if msg1==10:
					send_notification=False
				else:
					send_notification=True
				if options.verbose==False:
					print ('Message : '+error_msgs(100)+', '+error_msgs(msg1)+'\n')
				logger.info('Message : '+error_msgs(100)+', '+error_msgs(msg1)+'\n')
			else:
				msg_str='Message : '+error_msgs(msg)+'\n'
				if msg==10:
					send_notification=False
				else:
					send_notification=True
				if options.verbose==False:
					print ('Message : '+error_msgs(msg)+'\n')
				logger.info('Message : '+error_msgs(msg)+'\n')
		touch_file=inputs.basedir+'/.Finished_pcal_'+str(OBSID)+'_'+basemsdir+'_'+msbasename+'_'+str(msg)
		end_time=time.time()
		run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
		logger.info('Total runtime : '+str(run_time))
		msg_str='Dear PAIRCARS user,\n\nPolarisation self-calibration for : '+msbasename+'\n'+msg_str+'\nTotal runtime : '+str(run_time)+'\n\nBest regards,\nPAIRCARS developing team'
		msg_subject='Notification from PAIRCARS : Polarisation Selfcal : OBSID = '+str(OBSID)
		if type(msg)==int:
			if send_notification==True:
				attachments=glob.glob(options.workdir+'/quick_image_*.png')
				send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=attachments)
				os.system('rm -rf '+options.workdir+'/quick_image_*.png')
		os.system('touch '+touch_file)
		file_str=msbasename.split('.ms')[0]
		if os.path.isdir(basedir+'/logs/'+str(OBSID)+'/'+basemsdir)==False:
			os.makedirs(basedir+'/logs/'+str(OBSID)+'/'+basemsdir)
		if os.path.isdir(basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir)==False and inputs.keep_logger==True and eval(str(options.verbose)):
			os.makedirs(basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir)
		os.system('cp -r '+options.workdir+'/Pol_Selfcal.log '+basedir+'/logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
		if inputs.keep_logger and eval(str(options.verbose))==True:
			os.system('cp -r '+options.workdir+'/Pol_Selfcal_verbose.log '+basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
		os.system('rm -rf '+options.workdir+'/*.log '+options.workdir+'/TempLattice* '+options.workdir+'/I*image')
		os.system('rm -rf '+options.workdir+'/'+file_str+'* '+options.workdir+'/Backup_*.ms')	
	except Exception as e:
		touch_file=inputs.basedir+'/.Finished_pcal_'+str(OBSID)+'_'+basemsdir+'_'+msbasename+'_'+str('error')
		end_time=time.time()
		run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
		logger.info('#############################\n')
		logger.info('Polarisation selfcal failed for ms : '+options.chantime_msname+'\n')
		logger.info('Error occured : '+str(e)+'\n')
		logger.info('Total runtime : '+str(run_time)+'\n')
		logger.info('##############################\n')
		msg_str='Dear PAIRCARS user,\n\nPolarisation self-calibration for : '+msbasename+'\nMessage : Error in runtime : '+str(e)+'\nTotal runtime : '+\
					str(run_time)+'\n\nBest regards,\nPAIRCARS developing team'
		msg_subject='Notification from PAIRCARS : Polarisation Selfcal : OBSID = '+str(OBSID)
		if inputs.send_notification==True:
			attachments=glob.glob(options.workdir+'/quick_image_*.png')
			send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=attachments)
			os.system('rm -rf '+options.workdir+'/quick_image_*.png')
		os.system('touch '+touch_file)
		file_str=msbasename.split('.ms')[0]
		if os.path.isdir(basedir+'/logs/'+str(OBSID)+'/'+basemsdir)==False:
			os.makedirs(basedir+'/logs/'+str(OBSID)+'/'+basemsdir)
		if os.path.isdir(basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir)==False and inputs.keep_logger==True and eval(str(options.verbose)):
			os.makedirs(basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir)
		os.system('cp -r '+options.workdir+'/Pol_Selfcal.log '+basedir+'/logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
		if inputs.keep_logger and eval(str(options.verbose))==True:
			os.system('cp -r '+options.workdir+'/Pol_Selfcal_verbose.log '+basedir+'/verbose_logs/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.pollog')
		os.system('rm -rf '+options.workdir+'/*.log '+options.workdir+'/TempLattice* '+options.workdir+'/I*image')
		os.system('rm -rf '+options.workdir+'/'+file_str+'* '+options.workdir+'/Backup_*.ms')
		pass
