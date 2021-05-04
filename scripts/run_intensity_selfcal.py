import numpy as np,os,sys,matplotlib.pyplot as plt,time,logging,matplotlib,json,urllib.request,glob
from casatools import *
from casatasks import *
from paircars.access_ms import *
from paircars.basic_func import *
from paircars.intensity_selfcal_LTS import *
from paircars.flagger import *
from optparse import OptionParser
from astropy.io import fits
from astropy import wcs
from paircars.libpaircars import send_paircars_notification
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
Code is written by Devojyoti Kansabanik, 26 Jan, 2021
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
	try:
		imsubimage(imagename=imagename,outfile='temp.image',box=box)
	except:
		os.system('cp -r '+imagename+' temp.image')	
	imcollapse(imagename='temp.image',outfile='temp1.image',axes=2,function='mean')
	exportfits(imagename='temp1.image',fitsimage='temp.fits',dropdeg=True)
	data=fits.getdata('temp.fits')
	wlist=fits.getheader('temp.fits')
	w = wcs.WCS(wlist)
	fig = plt.figure(figsize=(8,8))
	ax = fig.add_subplot(111, projection = w)
	cax = fig.add_axes([0.91, 0.11, 0.02, 0.77])
	im=ax.imshow(data,cmap='hot',origin='lower')
	fig.colorbar(im,cax=cax,orientation='vertical')
	ax.set_xlabel('Right Ascension',fontsize=12)
	ax.set_ylabel('Declination',fontsize=12)
	ax.set_title('Frequency : '+str(freq)+' MHz, Timestamp : '+str(timestamp)+' UTC\n Dynamic range (rms) : '+str(int(DR_rms))+', Dynamic range (negative) : '+str(int(DR_neg)))
	ax.xaxis.set_tick_params(labelsize=10)
	ax.yaxis.set_tick_params(labelsize=10)
	cwd=os.getcwd()
	outfile_dir=os.path.dirname(outfile)
	if outfile_dir=='':
		outfile=cwd+'/'+outfile
	plt.savefig(outfile)
	os.system('rm -rf temp* casa*log')
	return outfile


# This part will run the self calibration loops. If the code is imported in some other python code, this part will not be executed

def run_intensity_selfcal(msname,metafits,working_dir,do_point_source=False,verbose=False,interactive=False,start_fresh=True,caltables=''):
	'''
	Heart of the intensity selfcal part of the PAIRCARS
	This function performs the intensity selfcal for PAIRCARS
	This script can be run directly from python IDE or can be imported to other scripts.
	Here some functionalities are chosen specific to MWA solar imaging. Which may not be valid for other instruments and imaging of other sources.
	Use this module only for Solar Imaging with MWA.
	Parameters:
	msname = Name of the measurement set
	metafits = Name of the metafits file
	working_dir = Name of the working directory
	do_point_source = False, Try with point source model
	verbose = False, If True keep all intermediate selfcal records
	interactive = False, If True perform interactive selfcal
	start_fresh = True, start fresh selfcal rounds from scratch or start from last round
	caltables = Previous caltables, comma separated
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
	
	if mspath[-1]=='/':
		mspath=mspath[:-1] 

	if mspath!=working_dir:
		if os.path.isdir(working_dir+'/'+os.path.basename(msname)):
			os.system('rm -rf '+working_dir+'/'+os.path.basename(msname))
		os.system('mv '+msname+' '+working_dir)
		msname=working_dir+'/'+os.path.basename(msname)

	os.chdir(working_dir)

	if __name__!='__main__':
		if (os.path.isfile(working_dir+'/Intensity_Selfcal.log') and start_fresh==True) or \
				(os.path.isfile(working_dir+'/Intensity_Selfcal.log') and os.path.isdir(working_dir+'/junk1.ms')==False and start_fresh==False):
			os.system('rm -rf '+working_dir+'/Intensity_Selfcal.log')
		if os.path.isfile(working_dir+'/Intensity_Selfcal_verbose.log') and verbose==True:
				os.system('rm -rf '+working_dir+'/Intensity_Selfcal_verbose.log')
		formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
	logger = logging.getLogger('intensity_selfcal_log')
	if __name__!='__main__':
		logger.setLevel(logging.DEBUG)
		if verbose==True:
			console=logging.StreamHandler(sys.stdout)
			console.setFormatter(formatter)
			logger.addHandler(console)
		filehandle=logging.FileHandler(working_dir+'/Intensity_Selfcal.log')
		filehandle.setFormatter(formatter)
		logger.addHandler(filehandle)
		logger.propagate = False
	print('\n')

	if start_fresh==False and os.path.isdir('junk1.ms'):
		os.system('rm -rf '+msname)
		os.system('cp -r junk1.ms '+msname)
	else:
		start_fresh=True
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
		else:
			if os.path.isdir(working_dir+'/Backup_uncalib.ms')==True:
				os.system('rm -rf '+working_dir+'/Backup_uncalib.ms')
			logger.info('split(vis=\''+msname+'\',outputvis=\''+working_dir+'/Backup_uncalib.ms\',datacolumn=\'data\')\n')
			split(vis=msname,outputvis=working_dir+'/Backup_uncalib.ms',datacolumn='data') # Backup of uncalibrated ms
		if caltables!='':
			caltable_list=caltables.split(',')
			logger.info('Applying solutions from previous calibrations : '+str(caltables)+'\n')
			logger.info('applycal(vis=\''+msname+'\',gaintable='+str(caltable_list)+',applymode=\'calonly\')\n')
			applycal(vis=msname,gaintable=caltable_list,applymode='calonly')
		
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
	datestr_list=msname.split('.ms')[0].split('_freq_')[0].split('time_')[1].split('_')
	datestr='/'.join(datestr_list[:3])+'/'+':'.join(datestr_list[3:]) # Datetime string 
	datestrfile='_'.join(datestr_list[:3])+'_'+'_'.join(datestr_list[3:]) # Datetime string for name 

	OBSID=get_OBSID(metafits)
	basemsdir=os.path.dirname(working_dir).split('/')[-1]
	if os.path.isdir(basedir+'/caltables/'+str(OBSID)+'/'+basemsdir)==False: # Directory to keep caltables
		os.makedirs(basedir+'/caltables/'+str(OBSID)+'/'+basemsdir)
	if os.path.isdir(basedir+'/imagemodels/'+str(OBSID)+'/'+basemsdir)==False: # Directory to keep models
		os.makedirs(basedir+'/imagemodels/'+str(OBSID)+'/'+basemsdir)


	if start_fresh==False:
		num_iter,DR1,DR3,DR5,DR2,DR4,DR6,rms_list,calmode,scratch,antenna_list_index,start_sigma,antenna_added,num_ant_current_iteration,\
					num_iter_fixed_sigma,num_iter_fixed_ant,num_iteration_after_ap,stokes,phasecenter_changed,startmodel,\
					startmask,uvsub_flag_count,ra,dec,num_iter_after_phasecenter_change,phasecenter_change_done=np.load('Intensity_selfcal_record.npy',allow_pickle=True)		
	if 'ref' in msname:
		if start_fresh:
			scratch=True # For reference time and frequency scratch = True
		if start_fresh==False:
			logger.info('Restarting selfcal from selfcal round : '+str(num_iter)+'\n')
		logger.info('Starting imaging for Reference time : '+str(datestr)+' and frequency : '+str(freqstr)+' MHz\n')
		logger.info('Scratch = '+str(scratch)+'\n')
		if verbose==False:
			if start_fresh==False:
				print('Restarting selfcal from selfcal round : '+str(num_iter)+'\n')
			print ('Starting imaging for Reference time : '+str(datestr)+' and frequency : '+str(freqstr)+' MHz\n')
			print ('Scratch = '+str(scratch)+'\n')
	else: # For other time and frequency
		if start_fresh: 
			scratch=False
		if start_fresh==False:
			logger.info('Restarting selfcal from selfcal round : '+str(num_iter)+'\n')
		logger.info('Reference time frequency slice imaging has been done. Starting imaging for time : '+str(datestr)+' and frequency : '+str(freqstr)+' MHz\n')
		logger.info('Scratch = '+str(scratch)+'\n')
		if verbose==False:
			if start_fresh==False:
				print('Restarting selfcal from selfcal round : '+str(num_iter)+'\n')
			print ('Reference time frequency slice imaging has been done. Starting imaging for time : '+str(datestr)+' and frequency : '+str(freqstr)+' MHz\n')
			print ('Scratch = '+str(scratch)+'\n')
	
	ISC=IntensitySelfcal(msname,metafits,32*60,verbose=verbose,interactive=interactive) # Creating selfcal object 32 arcmin maximum scale size
	AM=AccessMS(msname)

	###################
	# Putting user defined inputs if exisis or go with default values
	###################

	if calc_image_parameters==False:
		if cellsize!='' and cellsize!='nan':
			ISC.cellsize=inputs.cellsize
		if len(imsize)!=0 and imsize[0]!='nan':
			ISC.imsize=inputs.imsize
		if len(multiscale_scales)!=0 and multiscale_scales[0]!='nan':
			ISC.multiscale_scales=inputs.multiscale_scales
		if uvtaper!='':
			ISC.uvtaper=inputs.uvtaper
	
	if calc_selfcalib_params==False:
		if inputs.uvrange_to_cal!='':
			ISC.calib_uvrange=inputs.uvrange_to_cal

	end_selfcal=False
	while end_selfcal==False:
		if os.path.isfile(msname+'/.usedby_paircars')==False:
			os.system('touch '+msname+'/.usedby_paircars')

		###############################################################
		# Calculating minimum number of iterations and antenna bin size
		###############################################################
		
		min_num_iter_fixed_sigma,min_iteration,max_iteration,antenna_bin=ISC.calc_iter_num(inputs.safety_factor,inputs.quality_factor,scratch=scratch,bandpass_selfcal=False)

		logger.info('########################\n')
		logger.info('Estimating the number of selfcal iterations\n')
		logger.info('Minimum number of iteration at fixed sigma : '+str(min_num_iter_fixed_sigma)\
			+'; Minimum number of total selfcal iterations : '+str(min_iteration)+'; Maximum number of selfcal iterations : '+str(max_iteration)\
			+'; Antenna bins : '+str(antenna_bin)+'\n')

		if scratch==True and 'ref' not in msname: # If scratch=True due to failure, restore the flag and original data to start
			tb=table()
			tb.open(msname,nomodify=False)
			data=tb.getcol('DATA')
			flag=tb.getcol('FLAG')
			flag*=False
			try:
				tb.putcol('CORRECTED_DATA',data)		
			except:
				pass
			tb.putcol('FLAG',flag)
			tb.flush()
			tb.close()
			os.system('rm -rf '+msname+'.flagversions')	

		antenna_list,num_ant=AM.make_antenna_list(num_bins=antenna_bin)  # Making the antenna list

		###########################
		# Initiating loop variables
		###########################
		if start_fresh:
			do_selfcal=True
			if scratch==True:						
				do_ap=False
				stokes='I'
				calmode='p'
				antenna_added=True
				antenna_list_index=0
			else:
				do_ap=True
				stokes='XXYY'
				calmode='ap'
				antenna_added=False
				antenna_list_index=-1
			num_iter=0
			startmodel=''
			startmask=''
			uvsub_flag_count=0
			num_iter_fixed_sigma=0
			num_iter_fixed_ant=0	
			num_iter_after_phasecenter_change=0
			point_source_trial_count=0
			num_iteration_after_ap = 0
			nomask_try_count = 0
			try_nomask=False
			phasecenter_changed=False
			phasecenter_change_done=False
			ra=0
			dec=0
		else:
			do_selfcal=True
			antenna_list_index=antenna_list_index
			if calmode=='p':						
				do_ap=False
			else:
				do_ap=True
			num_iter=num_iter
			if os.path.isdir(startmodel):
				startmodel=startmodel
			else:
				startmodel=''
			if os.path.isdir(startmask):
				startmask=startmask
			else:
				startmask=''
			calmode=calmode
			antenna_added=antenna_added
			uvsub_flag_count=uvsub_flag_count
			num_iter_fixed_sigma=num_iter_fixed_sigma
			num_iter_fixed_ant=num_iter_fixed_ant
			point_source_trial_count=0
			num_iteration_after_ap = num_iteration_after_ap
			num_iter_after_phasecenter_change=num_iter_after_phasecenter_change
			phasecenter_change_done=phasecenter_change_done
			nomask_try_count = 0
			try_nomask=False
			phasecenter_changed=phasecenter_changed
			ra=ra
			dec=dec
			stokes=stokes

		if inputs.maskfile=='' and inputs.maskstr=='':
			mask_rad=int((32*60)/ISC.cellsize) # Creating a mask with 32 arcmin radius centered on the image
			mask_str='circle[['+str(ISC.imsize/2)+'pix,'+str(ISC.imsize/2)+'pix],'+str(mask_rad)+'pix]'
		elif inputs.maskstr!='':
			mask_str=inputs.maskstr

		file_str=os.path.basename(msname).split('.ms')[0]

		if start_fresh:
			if 'ref' not in msname:
				try:
					start_sigma=np.load(basedir+'/Ref_time_chan_sigma.npy',allow_pickle=True)[0]
					rms_list=np.load(basedir+'/Ref_time_chan_sigma.npy',allow_pickle=True)[1]
				except:
					logger.info('Start sigma and threshold information for last intensity selfcal round for reference time channel is not found.\n')
					if verbose==False:
						print('Start sigma and threshold information for last intensity selfcal round for reference time channel is not found.\n')
						os.chdir(cwd)
						if __name__!='__main__':
							touch_file=basedir+'/.Finished_gcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_12'
							msg_str='Dear PAIRCARS user,\n\nIntensity self-calibration for : '+\
									os.path.basename(msname)+'\nMessage :'+error_msgs(12)+'\n\nBest regards,\nPAIRCARS developing team'
							msg_subject='Notification from PAIRCARS : Intensity Selfcal : OBSID = '+str(OBSID)
							if inputs.send_notification==True:
								send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[])
							os.system('touch '+touch_file)
							if inputs.keep_logger==False:
								os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
							os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_uncalib.ms')
							end_time=time.time()
							run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
							logger.info('Total runtime : '+str(run_time))
						return 12
			else:
				start_sigma=inputs.start_sigma
		else:
			start_sigma=start_sigma
			rms_list=rms_list
						
		######################
		# Making dirty image for ref time
		######################
		if start_fresh:
			if 'ref' in msname or scratch==True:
				while True:
					logger.info('Start sigma : '+str(start_sigma)+'\n')
					if verbose==False:
						print('Start sigma : '+str(start_sigma)+'\n')
					msg_code,out_dict,negative_dyn_range,selfcal_snr=ISC.dirty_image(start_sigma,antenna_to_use=ISC.antenna_string(antenna_list,antenna_list_index))
					if 'ref' in msname:					
						ISC.file_remover_and_keeper('dirty',msg_code,do_bandpass=False,ref_time_chan=True)
					else:
						ISC.file_remover_and_keeper('dirty',msg_code,do_bandpass=False,ref_time_chan=False)
					if msg_code==0:
						logger.info('Initial selfcal SNR : '+str(selfcal_snr)+'\n')
						np.save(inputs.basedir+'/selfcal_minsnr',selfcal_snr)
						DR1=out_dict['I'][0]
						DR2=negative_dyn_range
						rms_list=[out_dict['I'][1]]
						ISC.DR_record(DR1,'DR_rms',init=True)
						ISC.DR_record(DR2,'DR_neg',init=True)
					if np.isnan(selfcal_snr):
						logger.info('No flux above for the present sigma threshold. Lowering threshold.\n')
						start_sigma=ISC.change_start_sigma(start_sigma,inputs.sigma_step,inputs.min_sigma)
						if np.isnan(start_sigma):
							if verbose==False:							
								print('Start sigma is below the minimum allowed sigma.\n')			
							logger.info('Start sigma is below the minimum allowed sigma.\n')
							os.chdir(cwd)
							if __name__!='__main__':
								touch_file=basedir+'/.Finished_gcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_11'
								msg_str='Dear PAIRCARS user,\n\nIntensity self-calibration for : '+\
										os.path.basename(msname)+'\nMessage :'+error_msgs(11)+'\n\nBest regards,\nPAIRCARS developing team'
								msg_subject='Notification from PAIRCARS : Intensity Selfcal : OBSID = '+str(OBSID)
								if inputs.send_notification==True:
									send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[])
								os.system('touch '+touch_file)
								if inputs.keep_logger==False:
									os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
								os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_uncalib.ms')
								end_time=time.time()
								run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
								logger.info('Total runtime : '+str(run_time))
							return 11
						else:
							logger.info('Trying with reducing start sigma.\n')
							if verbose==False:
								print('Trying with reducing start sigma.\n')							
							continue
					
					elif msg_code==0 and selfcal_snr<min_selfcal_snr:
						msg_code=110
						if __name__!='__main__':
							touch_file=basedir+'/.Finished_gcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_'+str(msg_code)
							msg_str='Dear PAIRCARS use,\n\nIntensity self-calibration for : '+os.path.basename(msname)+'\nMessage : '+error_msgs(100)+', '+error_msgs(msg_code-100)\
										+'\n\nBest regards,\nPAIRCARS developing team'
							msg_subject='Notification from PAIRCARS : Intensity Selfcal : OBSID = '+str(OBSID)
							if inputs.send_notification==True:
								send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[])
							os.system('touch '+touch_file)
							if inputs.keep_logger==False:
								os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
							if verbose==False:
								print('Message :'+error_msgs(10)) # If selfcal SNR is not sufficient
							logger.error('Message :'+error_msgs(10))
							os.chdir(cwd)
							os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_uncalib.ms')
							end_time=time.time()
							run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
							logger.info('Total runtime : '+str(run_time))
						return msg_code
					elif msg_code!=0:
						logger.error('Message :'+error_msgs(msg_code))
						if __name__!='__main__':
							touch_file=basedir+'/.Finished_gcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_'+str(msg_code)
							if 'ref' in msname:
								msg_str='Dear PAIRCARS user,\n\nIntensity self-calibration for : '+\
									os.path.basename(msname)+'\nMessage : '+error_msgs(100)+', '+error_msgs(msg_code)+'\n\nBest regards,\nPAIRCARS developing team'
							else:
								msg_str='Dear PAIRCARS user,\n\nIntensity self-calibration for : '+\
									os.path.basename(msname)+'\nMessage : '+error_msgs(msg_code)+'\n\nBest regards,\nPAIRCARS developing team'
							msg_subject='Notification from PAIRCARS : Intensity Selfcal : OBSID = '+str(OBSID)
							if inputs.send_notification==True:
								send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[])
							os.system('touch '+touch_file)
							if inputs.keep_logger==False:
								os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
							os.chdir(cwd)
							os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_uncalib.ms')
							end_time=time.time()
							run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
							logger.info('Total runtime : '+str(run_time))
						if 'ref' in msname:
							msg_code+=100
						return msg_code	# If dirty image is not made no point in continuing, return error
					else:
						logger.info('Dirty image is made successfully.\n')
						if verbose==False:
							print('Dirty image is made successfully.\n')
						break
			elif save_true_loc_image==True: # Save source true location images with respect to the reference time and channel
				logger.info('Source true location imaging is being done.\n')
				if verbose==False:
					print('Source true location imaging is being done.\n')
				if savedir=='':
					outdir=basedir+'/All_gaincal_source_true_loc_images'
				else:
					outdir=savedir+'/All_gaincal_source_true_loc_images'
				ISC.image_source_true_loc(outdir,do_bandpass=False)

		##############
		# Selfcal loop
		##############

		while do_selfcal==True:
			num_ant_current_iteration=len(antenna_list[antenna_list_index])
			if verbose==False:
				print ('#####################\nIntensity Selfcal iteration:'+str(num_iter)+'\n#####################\n')
			logger.info('####################\n')
			logger.info('Intensity Selfcal iteration:'+str(num_iter)+'\n')
			logger.info('#####################\n')

			if os.path.isdir(startmodel)==False:
				startmodel=''
			if os.path.isdir(startmask)==False:
				startmask=''
			if phasecenter_change_done==True and phasecenter_changed==True:
				phasecenter_changed=False

			if (num_iter<10 and nomask_try_count<1): 
					# Use a circular mask of the size of the Sun if calmode=='p' and no mask is provided by user. This is to keep th phasecenter fixed
				output_ISC=ISC.selfcal_iteration(num_iter,rms_list,start_sigma,mask_str,ISC.antenna_string(antenna_list,antenna_list_index),\
					startmodel,startmask,inputs.ref_ant,inputs.gain_minsnr,calmode,maskfile='',want_auto_masking=False,stokes=stokes,interactive=interactive,\
					do_bandpass=False,correct_phasecenter=phasecenter_changed,ra=ra,dec=dec,box_width=3,calibrator_caltable=[])  # Performing selfcal iterations
			elif inputs.maskfile!='' and try_nomask==False: # Use user defined mask
				mask_str=''
				output_ISC=ISC.selfcal_iteration(num_iter,rms_list,start_sigma,mask_str,ISC.antenna_string(antenna_list,antenna_list_index),\
					startmodel,startmask,inputs.ref_ant,inputs.gain_minsnr,calmode,maskfile=inputs.maskfile,want_auto_masking=False,stokes=stokes,interactive=interactive,\
					do_bandpass=False,correct_phasecenter=phasecenter_changed,ra=ra,dec=dec,box_width=3,calibrator_caltable=[]) # Performing selfcal iterations	
			elif inputs.maskstr!='' and try_nomask==False: # Use user defined mask string
				mask_str=inputs.maskstr
				output_ISC=ISC.selfcal_iteration(num_iter,rms_list,start_sigma,mask_str,ISC.antenna_string(antenna_list,antenna_list_index),\
					startmodel,startmask,inputs.ref_ant,inputs.gain_minsnr,calmode,maskfile=inputs.maskfile,want_auto_masking=False,stokes=stokes,interactive=interactive,\
					do_bandpass=False,correct_phasecenter=phasecenter_changed,ra=ra,dec=dec,box_width=3,calibrator_caltable=[]) # Performing selfcal iterations	
			elif try_nomask==True: # Using no mask to pickup flux from over the field
				mask_str=''
				output_ISC=ISC.selfcal_iteration(num_iter,rms_list,start_sigma,mask_str,ISC.antenna_string(antenna_list,antenna_list_index),\
					startmodel,startmask,inputs.ref_ant,inputs.gain_minsnr,calmode,maskfile='',want_auto_masking=False,stokes=stokes,interactive=interactive,\
					do_bandpass=False,correct_phasecenter=phasecenter_changed,ra=ra,dec=dec,box_width=3,calibrator_caltable=[])  # Performing selfcal iterations
				try_nomask=False
				nomask_try_count+=1
			elif inputs.maskfile=='' and inputs.maskstr=='' and inputs.want_auto_masking==False: # If no mask is given and auto masking off use a circular central mask
				output_ISC=ISC.selfcal_iteration(num_iter,rms_list,start_sigma,mask_str,ISC.antenna_string(antenna_list,antenna_list_index),\
					startmodel,startmask,inputs.ref_ant,inputs.gain_minsnr,calmode,maskfile='',want_auto_masking=inputs.want_auto_masking,stokes=stokes,interactive=interactive,\
					do_bandpass=False,correct_phasecenter=phasecenter_changed,ra=ra,dec=dec,box_width=3,calibrator_caltable=[])
			elif inputs.want_auto_masking==True:
				mask_str=''
				output_ISC=ISC.selfcal_iteration(num_iter,rms_list,start_sigma,mask_str,ISC.antenna_string(antenna_list,antenna_list_index),\
					startmodel,startmask,inputs.ref_ant,inputs.gain_minsnr,calmode,maskfile='',want_auto_masking=inputs.want_auto_masking,stokes=stokes,interactive=interactive,\
					do_bandpass=False,correct_phasecenter=phasecenter_changed,ra=ra,dec=dec,box_width=3,calibrator_caltable=[])

			if type(output_ISC)==tuple:				
				msg_code,out_dict,negative_dyn_range=output_ISC
			else:
				msg_code=output_ISC

			if phasecenter_changed==True and phasecenter_change_done==False:
				phasecenter_changed=False
				phasecenter_change_done=True

			if phasecenter_change_done==True:
				num_iter_after_phasecenter_change+=1

			if 'ref' in msname:		
				ISC.file_remover_and_keeper(num_iter,msg_code,do_bandpass=False,ref_time_chan=True)  # Removing files and keeping the required ones
			else:
				ISC.file_remover_and_keeper(num_iter,msg_code,do_bandpass=False,ref_time_chan=False)
			if num_iter==0 and msg_code==5:
				if nomask_try_count<1:
					logger.info('No flux is picked up in the model inside the mask. Trying with no mask.\n')
					if verbose==False:
						print('No flux is picked up in the model inside the mask. Trying with no mask.\n')
					nomask_try_count+=1		
					try_nomask=True
					do_selfcal=True
					continue
				logger.info('No flux is picked up in the model. Lowering threshold.\n')
				if verbose==False:	
					print('No flux is picked up in the model. Lowering threshold.\n')
				start_sigma=ISC.change_start_sigma(start_sigma,inputs.sigma_step,inputs.min_sigma)
				if np.isnan(start_sigma):
					if scratch==True:
						if do_point_source==True and point_source_trial_count<2: # Trying with point source model
							if point_source_trial_count<1:
								IB=ImageBasic(msname)
								uvrange='<'+IB.calc_suntaper()
								if verbose==False:
									print('Start sigma is below the minimum allowed sigma. Thus trying with a point source model for the uvrange : '+str(uvrange)+'\n')
								logger.info('Start sigma is below the minimum allowed sigma. Thus trying with a point source model for the uvrange : '+str(uvrange)+'\n')
							else:
								uvrange=''	
								if verbose==False:
									print('Start sigma is below the minimum allowed sigma. Thus trying with a point source model without any limit on uvrange \n')
								logger.info('Start sigma is below the minimum allowed sigma. Thus trying with a point source model without any limit on uvrange \n')
							start_sigma=inputs.min_sigma
							logger.info('clearcal(vis=\''+msname+'\',addmodel=True)\n')
							clearcal(vis=msname,addmodel=True)		 
							logger.info('gaincal(vis=\''+msname+'\',caltable=\'point_source_try.cal\',refant=\''+str(ref_ant)+'\',calmode=\'p\','+\
										'minsnr='+str(gain_minsnr)+',uvrange=\''+uvrange+'\',solint=\'inf\')\n')
							gaincal(vis=msname,caltable='point_source_try.cal',refant=str(ref_ant),calmode='p',minsnr=gain_minsnr,uvrange=uvrange,solint='inf')
							logger.info('applycal(vis=\''+msname+'\',gaintable=[\'point_source_try.cal\'],applymode=\'calflag\',flagbackup=True)\n')
							applycal(vis=msname,gaintable=['point_source_try.cal'],applymode='calflag',flagbackup=True)
							point_source_trial_count+=1
							continue
						else:
							end_selfcal=True
							if 'ref' in msname:
								os.chdir(cwd)
								if __name__!='__main__':
									touch_file=basedir+'/.Finished_gcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_106'
									msg_str='Dear PAIRCARS user,\n\nIntensity self-calibration for : '+os.path.basename(msname)+'\nMessage : '+error_msgs(100)+', '+error_msgs(6)\
											+'\n\nBest regards,\nPAIRCARS developing team'
									msg_subject='Notification from PAIRCARS : Intensity Selfcal : OBSID = '+str(OBSID)
									if inputs.send_notification==True:
										send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[])									
									os.system('touch '+touch_file)
									if inputs.keep_logger==False:
										os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
									os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_uncalib.ms')
									end_time=time.time()
									run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
									logger.info('Total runtime : '+str(run_time))
								return 106
							else:
								os.chdir(cwd)
								if __name__!='__main__':
									touch_file=basedir+'/.Finished_gcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_6'
									msg_str='Dear PAIRCARS user,\n\nIntensity self-calibration for : '+\
											os.path.basename(msname)+'\nMessage :'+error_msgs(6)+'\n\nBest regards,\nPAIRCARS developing team'
									msg_subject='Notification from PAIRCARS : Intensity Selfcal : OBSID = '+str(OBSID)
									if inputs.send_notification==True:
										send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[])
									os.system('touch '+touch_file)
									if inputs.keep_logger==False:
										os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
									os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_uncalib.ms')
									end_time=time.time()
									run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
									logger.info('Total runtime : '+str(run_time))
								return 6
					else:
						scratch=True
						if verbose==False:
							print ('######################\nGoing for a full selfcal from scratch = True\n #####################################\n')
						logger.info('######################\n')
						logger.info('Going for a full selfcal from scratch = True\n')
						logger.info('#####################################\n')
						os.system('rm -rf '+working_dir+'/junk*')
						do_selfcal=False
						break
				else:				
					logger.info('Trying with reduced start sigma : '+str(start_sigma)+'\n')
					if verbose==False:
						print('Trying with reducing start sigma : '+str(start_sigma)+'\n')							
					continue
			elif msg_code!=0:
				if verbose==False:
					print (error_msgs(msg_code))
				logger.error(error_msgs(msg_code))
				if scratch==True:
					end_selfcal=True
					if 'ref' in msname:
						os.chdir(cwd)
						if __name__!='__main__':
							touch_file=basedir+'/.Finished_gcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_'+str(msg_code+100)
							msg_str='Dear PAIRCARS user,\n\nIntensity self-calibration for : '+os.path.basename(msname)+'\nMessage : '+error_msgs(100)+', '+error_msgs(msg_code)\
											+'\n\nBest regards,\nPAIRCARS developing team'
							msg_subject='Notification from PAIRCARS : Intensity Selfcal : OBSID = '+str(OBSID)
							if inputs.send_notification==True:
								send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[])
							os.system('touch '+touch_file)
							if inputs.keep_logger==False:
								os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
							os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_uncalib.ms')
							end_time=time.time()
							run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
							logger.info('Total runtime : '+str(run_time))
						return msg_code+100
					else:
						os.chdir(cwd)
						if __name__!='__main__':
							touch_file=basedir+'/.Finished_gcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_'+str(msg_code)
							msg_str='Dear PAIRCARS user,\n\nIntensity self-calibration for : '+\
									os.path.basename(msname)+'\nMessage : '+error_msgs(6)+'\n\nBest regards,\nPAIRCARS developing team'
							msg_subject='Notification from PAIRCARS : Intensity Selfcal : OBSID = '+str(OBSID)
							if inputs.send_notification==True:
								send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[])
							os.system('touch '+touch_file)
							if inputs.keep_logger==False:
								os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
							os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_uncalib.ms')
							end_time=time.time()
							run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
							logger.info('Total runtime : '+str(run_time))
						return msg_code
				else:
					scratch=True
					if verbose==False:
						print ('######################\nGoing for a full selfcal from scratch = True\n #####################################\n')
					logger.info('######################\n')
					logger.info('Going for a full selfcal from scratch = True\n')
					logger.info('#####################################\n')
					os.system('rm -rf '+working_dir+'/junk*')
					do_selfcal=False
					break
			else:	
				if do_ap==False:
					dyn1=out_dict['I'][0]
					dyn2=negative_dyn_range
					rms_list=[out_dict['I'][1]]
				else:
					dyn1=(out_dict['XX'][0]+out_dict['YY'][0])/2.0
					dyn2=negative_dyn_range
					rms_list=[out_dict['XX'][1],out_dict['YY'][1]]
				if 'ref' not in msname and num_iter==0:
					DR1=dyn1
					DR2=dyn2				
					ISC.DR_record(DR1,'DR_rms',init=True)
					ISC.DR_record(DR2,'DR_neg',init=True)
				if num_iter==0:
					DR5=DR3=dyn1
					DR6=DR4=dyn2
				elif num_iter==1:
					DR5=dyn1
					DR6=dyn2
					ISC.DR_record(dyn1,'DR_rms',init=False)
					ISC.DR_record(dyn2,'DR_neg',init=False)
				else:
					DR1=DR3
					DR3=DR5
					DR5=dyn1
					DR2=DR4
					DR4=DR6 
					DR6=dyn2
					ISC.DR_record(dyn1,'DR_rms',init=False)
					ISC.DR_record(dyn2,'DR_neg',init=False)

				if os.path.isfile('Intensity_selfcal_record.npy'):
					os.system('rm -rf Intensity_selfcal_record.npy')
				selfcal_record=np.array([num_iter,DR1,DR3,DR5,DR2,DR4,DR6,rms_list,calmode,scratch,antenna_list_index,start_sigma,antenna_added,num_ant_current_iteration,\
					num_iter_fixed_sigma,num_iter_fixed_ant,num_iteration_after_ap,stokes,phasecenter_changed,startmodel,startmask,uvsub_flag_count,ra,dec,\
					num_iter_after_phasecenter_change,phasecenter_change_done],dtype='object')
				np.save('Intensity_selfcal_record',selfcal_record)

				if verbose==False:
					print ('RMS based dynamic ranges:\n')
					print(str(DR1)+','+str(DR3)+','+str(DR5)+'\n')
					print ('Negative based dynamic ranges:\n')
					print(str(DR2)+','+str(DR4)+','+str(DR6)+'\n')
					print ('Antenna added : '+str(antenna_added)+'\n')
					print ('Number of antennas in use : '+str(num_ant_current_iteration)+'\n')
					print ('Calmode : '+calmode+'\n')
					print ('Scartch = '+str(scratch)+'\n')
					print ('Sigma = '+str(start_sigma)+'\n')

				logger.info('RMS based dynamic ranges:\n')
				logger.info(str(DR1)+','+str(DR3)+','+str(DR5)+'\n')
				logger.info('Negative based dynamic ranges:\n')
				logger.info(str(DR2)+','+str(DR4)+','+str(DR6)+'\n')
				logger.info('Antenna added : '+str(antenna_added)+'\n')
				logger.info('Number of antennas in use : '+str(num_ant_current_iteration)+'\n')
				logger.info('Calmode : '+calmode+'\n')
				logger.info('Scratch = '+str(scratch)+'\n')
				logger.info('Sigma = '+str(start_sigma)+'\n')
			
				############## 
				# If statement 1 (DR decrease)
			
				if (((DR5<0.85*DR3 and DR5<0.9*DR1 and DR3>DR1) or (DR6<0.85*DR4 and DR6<0.9*DR2 and DR4>DR2)) and antenna_added==False and num_ant_current_iteration==num_ant)\
					or (((DR5<0.8*DR3 and DR5<0.85*DR1 and DR3>DR1) or (DR6<0.8*DR4 and DR6<0.85*DR2 and DR4>DR2)) and antenna_added==True and num_ant_current_iteration==num_ant)\
					or (((DR5<0.9*DR3 and DR1>1.5*DR3) or (DR6<0.9*DR4 and DR2>1.5*DR4)) and antenna_added==False and num_ant_current_iteration==num_ant):
					# If DR decreases.
					# Case 1: If DR decreases less than 90% and 85% of previous two rounds and all antennas are added. This is a check if the rms is diverging. 
					# Case 2: If DR decreases less than 85% and 80% of previous two rounds and all antennas are addded and last set of antennas are added in the last round.
					# Case 3: If DR decreases less than 90% of previous round but DR increases more than 1.5 times in last two rounds and no new antennas are added
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
						if use_ankflagger:
							logger.info('Performing uvsub flagging using aNKflagger due to DR decrease.\n')
							logger.info('do_uvsub_ankflag(\''+msname+'\',model=\'junk0.model\',nthread=1,verbose='+str(verbose)+',flagbackup=False)\n')
							do_uvsub_ankflag(msname,model='junk0.model',nthread=1,verbose=verbose,flagbackup=False)
						else:
							logger.info('Performing uvsub flagging due to DR decrease.\n')
							logger.info('do_uvsub_flagger(\''+msname+'\',model=\'junk0.model\',mode=\'uvsub_flag\',rmsthresh=[10,7,5,3.5],flagbackup=False)\n')
							do_uvsub_flagger(msname,model='junk0.model',mode='uvsub_flag',rmsthresh=[10,7,5,3.5],flagbackup=False)
						uvsub_flag_count+=1
						os.system('rm -rf junk1.model')
						os.system('cp -r junk0.model junk1.model')
						continue
					if scratch==True:
						if (do_ap==False or num_iteration_after_ap>min_iteration+5) :	
							# Doing a ap calibration may unsettle things. This gives the calibration some relaxation time to find its new stable position.
							# If scratch is False then it has failed in spite of a good starting point. Hence relaxation time is not needed.
							if DR5>min_DR: # Only considered the rms based DR here
								os.system('rm -rf '+file_str+'.cal')
								os.system('rm -rf '+file_str+'_'+str(num_iter)+'.model')
								ft(vis=working_dir+'/Backup_uncalib.ms',model='junk0.model',usescratch=True)
								logger.info('ft(vis=\''+working_dir+'/Backup_uncalib.ms\',model=\'junk0.model\',usescratch=True)\n')
								if inputs.uvrange_to_cal!='':
									IB=ImageBsic(working_dir+'/Backup_uncalib.ms')	
									uvrange_to_cal=IB.calc_calib_uvrange(4)[0]
								else:
									uvrange_to_cal=inputs.uvrange_to_cal
								logger.info('gaincal(vis=\''+working_dir+'/Backup_uncalib.ms\',caltable=\'temp.cal\',solmode=\'R\',rmsthresh=[10,8,6],calmode=\'ap\',uvrange=\''\
										+uvrange_to_cal+'\')\n')
								gaincal(vis=working_dir+'/Backup_uncalib.ms',caltable='temp.cal',solmode='R',rmsthresh=[10,8,6],calmode='ap',uvrange=uvrange_to_cal)
								os.system('mv temp.cal '+basedir+'/caltables/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.cal')  # Keeping the last good caltable
								os.system('cp -r junk0.model '+basedir+'/imagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.model') # Keeping last good model
								os.system('cp -r junk0.image '+basedir+'/imagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.image') # Keeping last good model
								if verbose==False:
									print (error_msgs(8))
								logger.error(error_msgs(8))
								end_selfcal=True
								if inputs.send_notification==True:
									quickimage=get_quicklook_image('junk0.image','quick_image_freq_'+freqstr+'_time_'+datestrfile+'.png',freqstr,datestr,DR3,DR4,field_of_view=2)
								os.system('rm -rf '+working_dir+'/junk*')
								if 'ref' in msname:
									os.chdir(cwd)
									if __name__!='__main__':
										touch_file=basedir+'/.Finished_gcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_108'
										msg_str='Dear PAIRCARS user,\n\nIntensity self-calibration for : '+os.path.basename(msname)+'\nMessage : '+error_msgs(100)+', '+error_msgs(8)\
											+'\n\nBest regards,\nPAIRCARS developing team'
										msg_subject='Notification from PAIRCARS : Intensity Selfcal : OBSID = '+str(OBSID)
										if inputs.send_notification==True:
											send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
											os.system('rm -rf '+quickimage)
										os.system('touch '+touch_file)
										if inputs.keep_logger==False:
											os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
										os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_uncalib.ms')
										end_time=time.time()
										run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
										logger.info('Total runtime : '+str(run_time))
									return 108
								else:
									os.chdir(cwd)
									if __name__!='__main__':
										touch_file=basedir+'/.Finished_gcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_8'
										msg_str='Dear PAIRCARS user,\n\nIntensity self-calibration for : '+\
											os.path.basename(msname)+'\nMessage : '+error_msgs(8)+'\n\nBest regards,\nPAIRCARS developing team'
										msg_subject='Notification from PAIRCARS : Intensity Selfcal : OBSID = '+str(OBSID)
										if inputs.send_notification==True:
											send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
											os.system('rm -rf '+quickimage)
										os.system('touch '+touch_file)
										if inputs.keep_logger==False:
											os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
										os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_uncalib.ms')
										end_time=time.time()
										run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
										logger.info('Total runtime : '+str(run_time))
									return 8
					else:
						scratch=True
						if verbose==False:
							print ('######################\nGoing for a full selfcal from scratch = True\n #####################################\n')
						logger.info('######################\n')
						logger.info('Going for a full selfcal from scratch = True\n')
						logger.info('#####################################\n')
						os.system('rm -rf '+working_dir+'/junk*')
						do_selfcal=False
						break

				#######################################################
				# If statement 2 (Exiting selfcal conditions)

				antenna_added=False
				if (DR5>=inputs.max_DR and num_ant_current_iteration==num_ant):
					if num_iteration_after_ap>min_iteration+5:
						if verbose==False:
							print ('Reached limiting dynamic range\n')
						logger.info('Reached limiting dynamic range\n')
						end_selfcal=True
						if 'ref' in msname:
							np.save(basedir+'/Ref_time_chan_sigma',np.array([start_sigma,rms_list],dtype='object'))
						logger.info('ft(vis=\''+working_dir+'/Backup_uncalib.ms\',model=\'junk1.model\',usescratch=True)\n')
						ft(vis=working_dir+'/Backup_uncalib.ms',model='junk1.model',usescratch=True)
						if inputs.uvrange_to_cal!='':
							IB=ImageBsic(working_dir+'/Backup_uncalib.ms')	
							uvrange_to_cal=IB.calc_calib_uvrange(4)[0]
						else:
							uvrange_to_cal=inputs.uvrange_to_cal
						logger.info('gaincal(vis=\''+working_dir+'/Backup_uncalib.ms\',caltable=\'temp.cal\',solmode=\'R\',rmsthresh=[10,8,6],calmode=\'ap\',uvrange=\''\
										+uvrange_to_cal+'\')\n')
						gaincal(vis=working_dir+'/Backup_uncalib.ms',caltable='temp.cal',solmode='R',rmsthresh=[10,8,6],calmode='ap',uvrange=uvrange_to_cal)
						os.system('mv temp.cal '+basedir+'/caltables/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.cal')  # Keeping the last good caltable
						os.system('cp -r junk1.model '+basedir+'/imagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.model')
						os.system('cp -r junk1.image '+basedir+'/imagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.image')
						if inputs.send_notification==True:
							quickimage=get_quicklook_image('junk1.image','quick_image_freq_'+freqstr+'_time_'+datestrfile+'.png',freqstr,datestr,DR5,DR6,field_of_view=2)
						os.system('rm -rf '+working_dir+'/junk*')
						os.chdir(cwd)
						if __name__!='__main__':
							touch_file=basedir+'/.Finished_gcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_0'
							msg_str='Dear PAIRCARS user,\n\nIntensity self-calibration for : '+\
								os.path.basename(msname)+'\nMessage : '+error_msgs(0)+'\n\nBest regards,\nPAIRCARS developing team'
							msg_subject='Notification from PAIRCARS : Intensity Selfcal : OBSID = '+str(OBSID)
							if inputs.send_notification==True:
								send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
								os.system('rm -rf '+quickimage)
							os.system('touch '+touch_file)
							if inputs.keep_logger==False:
								os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
							os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_uncalib.ms')
							end_time=time.time()
							run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
							logger.info('Total runtime : '+str(run_time))
						return 0
				elif (abs(DR5-DR3)<DR_delta_rms and abs(DR5-DR1)<DR_delta_rms and do_ap==True and abs(DR5/DR3-1)<0.08) and\
					 (abs(DR6-DR4)<DR_delta_neg and abs(DR6-DR2)<DR_delta_neg and do_ap==True and abs(DR6/DR4-1)<0.05):
				#  If DR does not increas more the DR delta in last two steps and DR does not increase 8% for rms based and 5% for negative based => Converge
					if num_iter_fixed_sigma>min_num_iter_fixed_sigma and num_iter>min_iteration:
						sigma=ISC.reduce_sigma('junk1.image',start_sigma,inputs.sigma_step,inputs.min_sigma,residual_frac=inputs.residual_frac,stokes_list=['XX','YY'])
						if sigma<start_sigma: # If the next sigma is less than the present sigma
							start_sigma=sigma	
							num_iter_fixed_sigma=0
						else:
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
								if use_ankflagger:
									logger.info('Perforing final uvsub flag using aNKflagger.\n')
									logger.info('do_uvsub_ankflag(\''+msname+'\',model=\'junk1.model\',nthread=1,verbose='+str(verbose)+',flagbackup=False)\n')
									do_uvsub_ankflag(msname,model='junk1.model',nthread=1,verbose=verbose,flagbackup=False)
								else:
									logger.info('Performing final uvsub flag.\n')
									logger.info('do_uvsub_flagger(\''+msname+'\',model=\'junk1.model\',mode=\'uvsub_flag\',rmsthresh=[10,7,5,3.5],flagbackup=False)\n')
									do_uvsub_flagger(msname,model='junk1.model',mode='uvsub_flag',rmsthresh=[10,7,5,3.5],flagbackup=False)
								uvsub_flag_count+=1
								continue
							else:
								if verbose==False:
									print ('#################\nSelfcal converged. Residual flux inside the mask is less than : '+\
											str(residual_frac*100)+'%. Stopped sigma : '+str(start_sigma)+'\n##################\n') 	
								logger.info('########################\n')							
								logger.info('Selfcal converged. Residual flux inside the mask is less than : '+str(residual_frac*100)+'%. Stopped sigma : '+str(start_sigma)+'\n')	
								logger.info('########################\n')								
								end_selfcal=True
								if 'ref' in msname:
									np.save(basedir+'/Ref_time_chan_sigma',np.array([start_sigma,rms_list],dtype='object'))
								logger.info('ft(vis=\''+working_dir+'/Backup_uncalib.ms\',model=\'junk1.model\',usescratch=True)\n')
								ft(vis=working_dir+'/Backup_uncalib.ms',model='junk1.model',usescratch=True)
								if inputs.uvrange_to_cal!='':
									IB=ImageBsic(working_dir+'/Backup_uncalib.ms')	
									uvrange_to_cal=IB.calc_calib_uvrange(4)[0]
								else:
									uvrange_to_cal=inputs.uvrange_to_cal
								logger.info('gaincal(vis=\''+working_dir+'/Backup_uncalib.ms\',caltable=\'temp.cal\',solmode=\'R\',rmsthresh=[10,8,6],calmode=\'ap\',uvrange=\''\
										+uvrange_to_cal+'\')\n')
								gaincal(vis=working_dir+'/Backup_uncalib.ms',caltable='temp.cal',solmode='R',rmsthresh=[10,8,6],calmode='ap',uvrange=uvrange_to_cal)
								os.system('mv temp.cal '+basedir+'/caltables/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.cal')  # Keeping the last good caltable
								os.system('cp -r junk1.model '+basedir+'/imagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.model')
								os.system('cp -r junk1.image '+basedir+'/imagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.image')
								if inputs.send_notification==True:
									quickimage=get_quicklook_image('junk1.image','quick_image_freq_'+freqstr+'_time_'+datestrfile+'.png',freqstr,datestr,DR5,DR6,field_of_view=2)
								os.system('rm -rf '+working_dir+'/junk*') 
								os.chdir(cwd)
								if __name__!='__main__':
									touch_file=basedir+'/.Finished_gcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_0'
									msg_str='Dear PAIRCARS user,\n\nIntensity self-calibration for : '+\
											os.path.basename(msname)+'\nMessage : '+error_msgs(0)+'\n\nBest regards,\nPAIRCARS developing team'
									msg_subject='Notification from PAIRCARS : Intensity Selfcal : OBSID = '+str(OBSID)
									if inputs.send_notification==True:
										sent=send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
										if sent==0:
											logger.info('Notification sent successfully.\n')
										else:
											logger.info('Notification could not be sent.\n')
										os.system('rm -rf '+quickimage)
									os.system('touch '+touch_file)
									if inputs.keep_logger==False:
										os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
									os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_uncalib.ms')
									end_time=time.time()
									run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
									logger.info('Total runtime : '+str(run_time))
								return 0
				elif (abs(DR5/DR3-1)<0.08 and abs(DR5/DR1-1)<0.08 and num_iter_fixed_ant>=5): # New antenna addition or calmode change
				# If fractional change of DR is less than 8% in last two steps and number of iterations at fixed antenna is greater than 5.
					if (num_ant_current_iteration<num_ant):
						antenna_list_index+=1
						num_iter_fixed_ant=0
						antenna_added=True	
						if nomask_try_count>=1 and scratch==True and antenna_list_index==1:
							phasecenter_changed=True
						if verbose==False:
							print ('New antenna added at iteration : '+str(num_iter)+'\n')
						logger.info('New antenna added at iteration : '+str(num_iter)+'\n')		
					else:
						if do_ap==False and phasecenter_change_done==False:
							ra,dec,phasecenter_changed=ISC.cal_solar_phaseshift('junk1.image',sigma=start_sigma)
							logger.info('Phase center changed required : '+str(phasecenter_changed)+'\n')
							if phasecenter_changed==True:							
								logger.info('New phasecenter : RA = '+str(ra)+' deg, DEC = '+str(dec)+' deg.\n')
						elif do_ap==False and num_iter_after_phasecenter_change>min_iteration:
							if verbose==False:
								print ('Change calmode to \'ap\' at iteration : '+str(num_iter)+'\n')
							logger.info('Change calmode to \'ap\' at iteration : '+str(num_iter)+'\n')
							do_ap=True
							calmode='ap'
							stokes='XXYY'
							rms_list=[out_dict['I'][1],out_dict['I'][1]]
		
				#############################################################			
				# If statement 3 (Using last round model) 
				#(If DR increases at least DR_delta and all antennas are added and number of iteration at fixed antenna is greater than 5)
				
				if ((DR5-DR3)>DR_delta_rms and (DR5-DR1)>DR_delta_rms) and (num_iter_fixed_ant>=5 or (num_ant_current_iteration==num_ant and num_iter_fixed_ant>=5))\
					 and phasecenter_changed==False:
					startmodel='junk1.model'
				else:
					startmodel=''
				if ((DR5-DR3)>DR_delta_rms and (DR5-DR1)>DR_delta_rms) and phasecenter_changed==False:
					startmask='junk1.mask'
				else:
					startmask=''
				num_iter+=1
				num_iter_fixed_sigma+=1
				num_iter_fixed_ant+=1
				if do_ap==True:
					num_iteration_after_ap+=1
				
				###############################################################
				# If statement 4 (Reached maximum selfcal rounds)
				if num_iter>max_iteration:
					if scratch==True:
						if DR5>min_DR:
							os.system('rm -rf '+file_str+'.cal')
							os.system('cp -r junk0.cal junk1.cal')  # Keeping the last good caltable
							os.system('cp -r junk0.ms junk1.ms')  # Keeping the last good calibrated ms
							if verbose==False:
								print (error_msgs(9))
							logger.error(error_msgs(9))
							end_selfcal=True
							logger.info('ft(vis=\''+working_dir+'/Backup_uncalib.ms\',model=\'junk1.model\',usescratch=True)\n')
							ft(vis=working_dir+'/Backup_uncalib.ms',model='junk1.model',usescratch=True)
							if inputs.uvrange_to_cal!='':
								IB=ImageBsic(working_dir+'/Backup_uncalib.ms')	
								uvrange_to_cal=IB.calc_calib_uvrange(4)[0]
							else:
								uvrange_to_cal=inputs.uvrange_to_cal
							logger.info('gaincal(vis=\''+working_dir+'/Backup_uncalib.ms\',caltable=\'temp.cal\',solmode=\'R\',rmsthresh=[10,8,6],calmode=\'ap\',uvrange=\''\
										+uvrange_to_cal+'\')\n')
							gaincal(vis=working_dir+'/Backup_uncalib.ms',caltable='temp.cal',solmode='R',rmsthresh=[10,8,6],calmode='ap',uvrange=uvrange_to_cal)
							os.system('mv temp.cal '+basedir+'/caltables/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.cal')  # Keeping the last good caltable
							os.system('cp -r junk1.model '+basedir+'/imagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.model')
							os.system('cp -r junk1.image '+basedir+'/imagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.image')
							if inputs.send_notification==True:
								quickimage=get_quicklook_image('junk1.image','quick_image_freq_'+freqstr+'_time_'+datestrfile+'.png',freqstr,datestr,DR5,DR6,field_of_view=2)
							os.system('rm -rf '+working_dir+'/junk*')
							if 'ref' in msname:
								np.save(basedir+'/Ref_time_chan_sigma',np.array([start_sigma,rms_list],dtype='object'))	
								os.chdir(cwd)
								if __name__!='__main__':
									touch_file=basedir+'/.Finished_gcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_109'
									msg_str='Dear PAIRCARS User,\n\nIntensity self-calibration for : '+os.path.basename(msname)+'\nMessage : '+error_msgs(100)+', '+error_msgs(9)\
												+'\n\nBest regards,\nPAIRCARS developing team'
									msg_subject='Notification from PAIRCARS : Intensity Selfcal : OBSID = '+str(OBSID)
									if inputs.send_notification==True:
										send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
										os.system('rm -rf '+quickimage)
									os.system('touch '+touch_file)
									if inputs.keep_logger==False:
										os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
									os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_uncalib.ms')
									end_time=time.time()
									run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
									logger.info('Total runtime : '+str(run_time))
								return 109		
							else:
								os.chdir(cwd)
								if __name__!='__main__':
									touch_file=basedir+'/.Finished_gcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_9'
									msg_str='Dear PAIRCARS user,\n\nIntensity self-calibration for : '+\
										os.path.basename(msname)+'\nMessage : '+error_msgs(9)+'\n\nBest regards,\nPAIRCARS developing team'
									msg_subject='Notification from PAIRCARS : Intensity Selfcal : OBSID = '+str(OBSID)
									if inputs.send_notification==True:
										send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
										os.system('rm -rf '+quickimage)
									os.system('touch '+touch_file)
									if inputs.keep_logger==False:
										os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
									os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_uncalib.ms')
									end_time=time.time()
									run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
									logger.info('Total runtime : '+str(run_time))
								return 9
						else:
							if verbose==False:
								print (error_msgs(13))
							logger.error(error_msgs(13))
							end_selfcal=True
							os.system('rm -rf '+working_dir+'/junk*')
							os.chdir(cwd)
							if __name__!='__main__':
								touch_file=basedir+'/.Finished_gcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_13'
								msg_str='Dear PAIRCARS user,\n\nIntensity self-calibration for : '+\
									os.path.basename(msname)+'\nMessage : '+error_msgs(13)+'\n\nBest regards,\nPAIRCARS developing team'
								msg_subject='Notification from PAIRCARS : Intensity Selfcal : OBSID = '+str(OBSID)
								if inputs.send_notification==True:
									send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
									os.system('rm -rf '+quickimage)
								os.system('touch '+touch_file)
								if inputs.keep_logger==False:
									os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
								os.system('rm -rf '+working_dir+'/'+file_str+'* '+working_dir+'/Backup_uncalib.ms')
								end_time=time.time()
								run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
								logger.info('Total runtime : '+str(run_time))
							return 13
					else:
						scratch=True
						if verbose==False:
							print ('######################\nGoing for a full selfcal from scratch = True\n #####################################\n')
						logger.info('######################\n')
						logger.info('Going for a full selfcal from scratch = True\n')
						logger.info('#####################################\n')
						os.system('rm -rf '+working_dir+'/junk*')
						do_selfcal=False
						break

# Function to run the script stand alone from command line
if __name__=='__main__':
	start_time=time.time()
	usage= ' Perform self calibration of a single time and frequency slice'
	parser = OptionParser(usage=usage)
	parser.add_option('--msname',dest="chantime_msname",default=None,help="Name of measurement set of a single time anf frequency slice",metavar="Measurement Set")
	parser.add_option('--metafits',dest="metafits",default=None,help="Name of metafits file of the observation",metavar="Metafits file")
	parser.add_option('--workdir',dest='workdir',default=None,help='Name of the working directory',metavar='Directory path')
	parser.add_option('--dopoint',dest='do_point_source',default=False,help='Want to try with point source model',metavar="Boolean")
	parser.add_option('--verbose',dest="verbose",default=False,help="Verbose mode",metavar="Boolean")
	parser.add_option('--interactive',dest="interactive",default=False,help="Interactive mode",metavar="Boolean")
	parser.add_option('--fresh',dest="fresh",default=True,help="Start fresh self calibration loop",metavar="Boolean")
	parser.add_option('--caltables',dest="caltables",default='',help="Previous caltables",metavar="String, comma separated")
	(options, args) = parser.parse_args()
	if (os.path.isfile(str(options.workdir)+'/Intensity_Selfcal.log') and eval(str(options.fresh))==True) or \
			(os.path.isfile(str(options.workdir)+'/Intensity_Selfcal.log') and os.path.isdir(str(options.workdir)+'/junk1.ms')==False and eval(str(options.fresh))==False):
		os.system('rm -rf '+str(options.workdir)+'/Intensity_Selfcal.log')
	if os.path.isfile(str(options.workdir)+'/Intensity_Selfcal_verbose.log') and eval(str(options.verbose))==True:
			os.system('rm -rf '+str(options.workdir)+'/Intensity_Selfcal_verbose.log')
	formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
	logger = logging.getLogger('intensity_selfcal_log')
	logger.setLevel(logging.DEBUG)
	if eval(str(options.verbose))==True:
		console=logging.StreamHandler(sys.stdout)
		console.setFormatter(formatter)
		logger.addHandler(console)
	filehandle=logging.FileHandler(str(options.workdir)+'/Intensity_Selfcal.log')
	filehandle.setFormatter(formatter)
	logger.addHandler(filehandle)
	logger.propagate = False

	cwd=os.getcwd()
	sys.path.append(cwd)
	if os.path.isfile(cwd+'/selfcal_inputs.py')==False:
		print('Input file does not exist.\n')
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
		touch_file=inputs.basedir+'/.Finished_gcal_'+str(OBSID)+'_'+basemsdir+'_'+msbasename+'_'+str('noms')
		end_time=time.time()
		run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
		logger.info('#############################\n')
		logger.info('Gain selfcal failed for ms : '+options.chantime_msname+'\n')
		logger.info('Total runtime : '+str(run_time)+'\n')
		logger.info('##############################\n')
		msg_str='Dear PAIRCARS user,\n\nIntensity self-calibration for : '+msbasename+'\nMessage : No measurement set is present\nTotal runtime : '+\
					str(run_time)+'\n\nBest regards,\nPAIRCARS developing team'
		msg_subject='Notification from PAIRCARS : Intensity Selfcal : OBSID = '+str(OBSID)
		if inputs.send_notification==True:
			send_paircars_notification(inputs.email,msg_subject,msg_str)
		os.system('touch '+touch_file)
		if inputs.keep_logger==False:
			os.system('rm -rf '+options.workdir+'/*.log')
		os.system('rm -rf '+options.workdir+'/TempLattice*')
		file_str=msbasename.split('.ms')[0]
		os.system('rm -rf '+options.workdir+'/'+file_str+'*')
		os._exit(0)
	
	if options.metafits==None or os.path.isfile(options.metafits)==False:
		logger.info('Metafits file does not exist. Exititing...\n')
		touch_file=inputs.basedir+'/.Finished_gcal_'+str(OBSID)+'_'+basemsdir+'_'+msbasename+'_'+str('nometa')
		end_time=time.time()
		run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
		logger.info('#############################\n')
		logger.info('Gain selfcal failed for ms : '+options.chantime_msname+'\n')
		logger.info('Total runtime : '+str(run_time)+'\n')
		logger.info('##############################\n')
		msg_str='Dear PAIRCARS user,\n\nIntensity self-calibration for : '+msbasename+'\nMessage : No metafits file is present\nTotal runtime : '+\
					str(run_time)+'\n\nBest regards,\nPAIRCARS developing team'
		msg_subject='Notification from PAIRCARS : Intensity Selfcal : OBSID = '+str(OBSID)
		if inputs.send_notification==True:
			send_paircars_notification(inputs.email,msg_subject,msg_str)
		os.system('touch '+touch_file)
		if inputs.keep_logger==False:
			os.system('rm -rf '+options.workdir+'/*.log')
		os.system('rm -rf '+options.workdir+'/TempLattice*')
		file_str=msbasename.split('.ms')[0]
		os.system('rm -rf '+options.workdir+'/'+file_str+'*')
		os._exit(0)

	try:
		print ('\n\t##########################\n\tStarting Intensity self-calibration.....\n\t##########################\n')
		print ('run_intensity_selfcal(\''+options.chantime_msname+'\',\''+options.metafits+'\',\''+options.workdir+'\',do_point_source='+str(options.do_point_source)+\
				',verbose='+str(options.verbose)+',interactive='+str(options.interactive)+',start_fresh='+str(options.fresh)+',caltables=\''+str(options.caltables)+'\')\n')
		msg=run_intensity_selfcal(options.chantime_msname,options.metafits,options.workdir,do_point_source=eval(str(options.do_point_source)),verbose=eval(str(options.verbose)),\
				interactive=eval(str(options.interactive)),start_fresh=eval(str(options.fresh)),caltables=str(options.caltables))
		if msg>100:
			msg1=msg-100
			msg_str='Message : '+error_msgs(100)+', '+error_msgs(msg1)+'\n'
			if options.verbose==False:
				print ('Message : '+error_msgs(100)+', '+error_msgs(msg1)+'\n')
			logger.info('Message : '+error_msgs(100)+', '+error_msgs(msg1)+'\n')
		else:
			msg_str='Message : '+error_msgs(msg)+'\n'
			if options.verbose==False:
				print ('Message : '+error_msgs(msg)+'\n')
			logger.info('Message : '+error_msgs(msg)+'\n')
		touch_file=inputs.basedir+'/.Finished_gcal_'+str(OBSID)+'_'+basemsdir+'_'+msbasename+'_'+str(msg)
		end_time=time.time()
		run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
		logger.info('#############################\n')
		logger.info('Gain selfcal finished for ms : '+options.chantime_msname+'\n')
		logger.info('Total runtime : '+str(run_time)+'\n')
		logger.info('##############################\n')
		msg_str='Dear PAIRCARS user,\n\nIntensity self-calibration for : '+msbasename+'\n'+msg_str+'\nTotal runtime : '+str(run_time)+'\n\nBest regards,\nPAIRCARS developing team'
		msg_subject='Notification from PAIRCARS : Intensity Selfcal : OBSID = '+str(OBSID)
		if inputs.send_notification==True:
			attachments=glob.glob(options.workdir+'/quick_image_*.png')
			send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=attachments)
			os.system('rm -rf '+options.workdir+'/quick_image_*.png')
		os.system('touch '+touch_file)
		if inputs.keep_logger==False:
			os.system('rm -rf '+options.workdir+'/*.log')
		os.system('rm -rf '+options.workdir+'/TempLattice*')
		file_str=msbasename.split('.ms')[0]
		os.system('rm -rf '+options.workdir+'/'+file_str+'* '+options.workdir+'/Backup_uncalib.ms')
	except Exception as e:
		touch_file=inputs.basedir+'/.Finished_gcal_'+str(OBSID)+'_'+basemsdir+'_'+msbasename+'_'+str('error')
		end_time=time.time()
		run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
		logger.info('#############################\n')
		logger.info('Gain selfcal failed for ms : '+options.chantime_msname+'\n')
		logger.info('Error occured : '+str(e)+'\n')
		logger.info('Total runtime : '+str(run_time)+'\n')
		logger.info('##############################\n')
		msg_str='Dear PAIRCARS user,\n\nIntensity self-calibration for : '+msbasename+'\nMessage : Error in runtime : '+str(e)+'\nTotal runtime : '+\
					str(run_time)+'\n\nBest regards,\nPAIRCARS developing team'
		msg_subject='Notification from PAIRCARS : Intensity Selfcal : OBSID = '+str(OBSID)
		if inputs.send_notification==True:
			send_paircars_notification(inputs.email,msg_subject,msg_str)
		os.system('touch '+touch_file)
		if inputs.keep_logger==False:
			os.system('rm -rf '+options.workdir+'/*.log')
		os.system('rm -rf '+options.workdir+'/TempLattice*')
		file_str=msbasename.split('.ms')[0]
		os.system('rm -rf '+options.workdir+'/'+file_str+'* '+options.workdir+'/Backup_uncalib.ms')
		pass
