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
Code is written by Devojyoti Kansabanik, 07 Mar, 2021
'''

def get_OBSID(metafits):
	'''
	Function to return OBSID of an MWA observation
	Parameters:
	metafits = Name of the metafits file
	Return:
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

def run_bandpass_selfcal(msname,metafits,working_dir,verbose=False,interactive=False,start_fresh=True):
	'''
	Heart of the bandpass selfcal part of the PAIRCARS
	This function performs the bandpass selfcal for PAIRCARS
	This script can be run directly from python IDE or can be imported to other scripts.
	Here some functionalities are chosen specific to MWA solar imaging. Which may not be valid for other instruments and imaging of other sources.
	Use this module only for Solar Imaging with MWA.
	Parameters:
	msname = Name of the measurement set
	metafits = Name of the metafits fil
	working_dir = Name of the working directory
	verbose = False, If True keep all intermediate selfcal records
	interactive = False, If True perform interactive selfcal
	start_fresh = True, start fresh selfcal rounds from scratch or start from last round
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
		if (os.path.isfile(working_dir+'/Bandpass_Selfcal.log') and start_fresh==True) or \
				(os.path.isfile(working_dir+'/Bandpass_Selfcal.log') and os.path.isdir(working_dir+'/junk1.ms')==False and start_fresh==False):
			os.system('rm -rf '+working_dir+'/Bandpass_Selfcal.log')
		if os.path.isfile(working_dir+'/Bandpass_Selfcal_verbose.log') and verbose==True:
				os.system('rm -rf '+working_dir+'/Bandpass_Selfcal_verbose.log')
		formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
	logger = logging.getLogger('bandpass_selfcal_log')
	if __name__!='__main__':
		logger.setLevel(logging.DEBUG)
		if verbose==True:
			console=logging.StreamHandler(sys.stdout)
			console.setFormatter(formatter)
			logger.addHandler(console)
		filehandle=logging.FileHandler(working_dir+'/Bandpass_Selfcal.log')
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
	if os.path.isdir(basedir+'/bpcaltables/'+str(OBSID)+'/'+basemsdir)==False: # Directory to keep caltables
		os.makedirs(basedir+'/bpcaltables/'+str(OBSID)+'/'+basemsdir)
	if os.path.isdir(basedir+'/bpimagemodels/'+str(OBSID)+'/'+basemsdir)==False: # Directory to keep models
		os.makedirs(basedir+'/bpimagemodels/'+str(OBSID)+'/'+basemsdir)
	
	if start_fresh==False:
		num_iter,DR1,DR3,DR5,DR2,DR4,DR6,rms_list,calmode,scratch,antenna_list_index,start_sigma,num_iter_fixed_sigma,num_iter_fixed_ant,stokes,startmodel,\
					startmask,uvsub_flag_count=np.load('Bandpass_selfcal_record.npy',allow_pickle=True)		
	
	scratch=False # For bandpass selfcal scratch is always False, since gain calibration is already applied
	if 'ref' in msname:
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
		if uvrange_to_cal!='':
			ISC.calib_uvrange=uvrange_to_cal

	end_selfcal=False

	if os.path.isfile(msname+'/.usedby_paircars')==False:
		os.system('touch '+msname+'/.usedby_paircars')

	###############################################################
	# Calculating minimum number of iterations and antenna bin size
	###############################################################
	
	min_num_iter_fixed_sigma,min_iteration,max_iteration,antenna_bin=ISC.calc_iter_num(inputs.safety_factor,inputs.quality_factor,scratch=scratch,bandpass_selfcal=True)

	logger.info('########################\n')
	logger.info('Estimating the number of selfcal iterations\n')
	logger.info('Minimum number of iteration at fixed sigma : '+str(min_num_iter_fixed_sigma)\
		+'; Minimum number of total selfcal iterations : '+str(min_iteration)+'; Maximum number of selfcal iterations : '+str(max_iteration)\
		+'; Antenna bins : '+str(antenna_bin)+'\n')

	antenna_list,num_ant=AM.make_antenna_list(num_bins=antenna_bin)  # Making the antenna list
	
	###########################
	# Initiating loop variables
	###########################

	if start_fresh:
		stokes='XXYY'
		calmode='ap'
		antenna_list_index=-1
		num_iter=0
		startmodel=''
		startmask=''
		uvsub_flag_count=0
		num_iter_fixed_sigma=0
		num_iter_fixed_ant=0				
	else:
		antenna_list_index=antenna_list_index
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
		uvsub_flag_count=uvsub_flag_count
		num_iter_fixed_sigma=num_iter_fixed_sigma
		num_iter_fixed_ant=num_iter_fixed_ant
		stokes=stokes

	if inputs.maskfile=='' and inputs.maskstr=='':
		mask_rad=int((32*60)/ISC.cellsize) # Creating a mask with 32 arcmin radius centered on the image
		mask_str='circle[['+str(ISC.imsize/2)+'pix,'+str(ISC.imsize/2)+'pix],'+str(mask_rad)+'pix]'
	elif inputs.maskstr!='':
		mask_str=inputs.maskstr

	file_str=os.path.basename(msname).split('.ms')[0]

	try:
		start_sigma=np.load(basedir+'/Ref_time_chan_sigma.npy',allow_pickle=True)[0] # Starting with last gaincal start_sigma and threshold
		rms_list=np.load(basedir+'/Ref_time_chan_sigma.npy',allow_pickle=True)[1]
		if len(rms_list)!=2:
			rms_list.append(rms_list[0])
	except:
		logger.info('Start sigma and threshold information for last intensity selfcal round for reference time channel is not found.\n')
		if verbose==False:
			print('Start sigma and threshold information for last intensity selfcal round for reference time channel is not found.\n')
		os.chdir(cwd)
		if __name__!='__main__':
			touch_file=basedir+'/.Finished_bcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_12'
			msg_str='Dear PAIRCARS user,\n\nBandpass self-calibration for : '+\
					os.path.basename(msname)+'\nMessage :'+error_msgs(12)+'\n\nBest regards,\nPAIRCARS developing team'
			msg_subject='Notification from PAIRCARS : Bandpass Selfcal : OBSID = '+str(OBSID)
			if inputs.send_notification==True:
				send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[])
			os.system('touch '+touch_file)
			if inputs.keep_logger==False:
				os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
			os.system('rm -rf '+working_dir+'/'+file_str+'*')
			end_time=time.time()
			run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
			logger.info('Total runtime : '+str(run_time))
		return 12

	if save_true_loc_image==True: # Save source true location images with respect to the reference time and channel
		logger.info('Source true location imaging is being done.\n')
		if verbose==False:
			print('Source true location imaging is being done.\n')
		if savedir=='':
			outdir=basedir+'/All_bandpass_source_true_loc_images'
		else:
			outdir=savedir+'/All_bandpass_source_true_loc_images'
		ISC.image_source_true_loc(outdir,do_bandpass=True)
	
	while end_selfcal==False:
			
		##############
		# Selfcal loop
		##############
		if verbose==False:
			print ('#####################\nBandpass Selfcal iteration:'+str(num_iter)+'\n#####################\n')
		logger.info('####################\n')
		logger.info('Bandpass Selfcal iteration:'+str(num_iter)+'\n')
		logger.info('#####################\n')

		if os.path.isdir(startmodel)==False:
			startmodel=''
		if os.path.isdir(startmask)==False:
			startmask=''

		if inputs.maskfile!='': # Use user defined mask
			mask_str=''
			output_ISC=ISC.selfcal_iteration(num_iter,rms_list,start_sigma,mask_str,ISC.antenna_string(antenna_list,antenna_list_index),\
				startmodel,startmask,inputs.ref_ant,inputs.gain_minsnr,calmode,maskfile=inputs.maskfile,want_auto_masking=False,stokes=stokes,interactive=interactive,\
							do_bandpass=True,correct_phasecenter=False,box_width=3,calibrator_caltable=[])  # Performing selfcal iterations	
		elif inputs.maskstr!='': # Use user defined mask string
			mask_str=inputs.maskstr
			output_ISC=ISC.selfcal_iteration(num_iter,rms_list,start_sigma,mask_str,ISC.antenna_string(antenna_list,antenna_list_index),\
				startmodel,startmask,inputs.ref_ant,inputs.gain_minsnr,calmode,maskfile=inputs.maskfile,want_auto_masking=False,stokes=stokes,interactive=interactive,\
				do_bandpass=True,correct_phasecenter=False,box_width=3,calibrator_caltable=[])  # Performing selfcal iterations		
		elif inputs.maskfile=='' and inputs.maskstr=='' and inputs.want_auto_masking==False: # If no mask is given and auto masking off use a circular central mask
			output_ISC=ISC.selfcal_iteration(num_iter,rms_list,start_sigma,mask_str,ISC.antenna_string(antenna_list,antenna_list_index),\
				startmodel,startmask,inputs.ref_ant,inputs.gain_minsnr,calmode,maskfile='',want_auto_masking=inputs.want_auto_masking,stokes=stokes,interactive=interactive,\
							do_bandpass=True,correct_phasecenter=False,box_width=3,calibrator_caltable=[])
		elif inputs.want_auto_masking==True:
			maskstr=''
			output_ISC=ISC.selfcal_iteration(num_iter,rms_list,start_sigma,mask_str,ISC.antenna_string(antenna_list,antenna_list_index),\
				startmodel,startmask,inputs.ref_ant,inputs.gain_minsnr,calmode,maskfile='',want_auto_masking=inputs.want_auto_masking,stokes=stokes,interactive=interactive,\
							do_bandpass=True,correct_phasecenter=False,box_width=3,calibrator_caltable=[])

		if type(output_ISC)==tuple:				
			msg_code,out_dict,negative_dyn_range=output_ISC
		else:
			msg_code=output_ISC

		if 'ref' in msname:		
			ISC.file_remover_and_keeper(num_iter,msg_code,do_bandpass=True,ref_time_chan=True)  # Removing files and keeping the required ones
		else:
			ISC.file_remover_and_keeper(num_iter,msg_code,do_bandpass=True,ref_time_chan=False)
		if msg_code!=0:
			if verbose==False:
				print (error_msgs(msg_code))
			logger.error(error_msgs(msg_code))
			end_selfcal=True
			if 'ref' in msname:
				os.chdir(cwd)
				if __name__!='__main__':
					touch_file=basedir+'/.Finished_bcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_'+str(msg_code+100)
					msg_str='Dear PAIRCARS user,\n\nBandpass self-calibration for : '+os.path.basename(msname)+'\nMessage : '+error_msgs(100)+', '+error_msgs(msg_code)\
									+'\n\nBest regards,\nPAIRCARS developing team'
					msg_subject='Notification from PAIRCARS : Bandpass Selfcal : OBSID = '+str(OBSID)
					if inputs.send_notification==True:
						send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[])
					os.system('touch '+touch_file)
					if inputs.keep_logger==False:
						os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
					os.system('rm -rf '+working_dir+'/'+file_str+'*')
					end_time=time.time()
					run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
					logger.info('Total runtime : '+str(run_time))
				return msg_code+100
			else:
				os.chdir(cwd)
				if __name__!='__main__':
					touch_file=basedir+'/.Finished_bcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_'+str(msg_code)
					msg_str='Dear PAIRCARS user,\n\nBandpass self-calibration for : '+\
							os.path.basename(msname)+'\nMessage : '+error_msgs(6)+'\n\nBest regards,\nPAIRCARS developing team'
					msg_subject='Notification from PAIRCARS : Bandpass Selfcal : OBSID = '+str(OBSID)
					if inputs.send_notification==True:
						send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[])
					os.system('touch '+touch_file)
					if inputs.keep_logger==False:
						os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
					os.system('rm -rf '+working_dir+'/'+file_str+'*')
					end_time=time.time()
					run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
					logger.info('Total runtime : '+str(run_time))
				return msg_code
		else:	
			dyn1=(out_dict['XX'][0]+out_dict['YY'][0])/2.0
			dyn2=negative_dyn_range
			rms_list=[out_dict['XX'][1],out_dict['YY'][1]]				
			if num_iter==0:
				ISC.DR_record(dyn1,'DR_rms',init=True)
				ISC.DR_record(dyn2,'DR_neg',init=True)
				DR5=DR3=DR1=dyn1
				DR6=DR4=DR2=dyn2
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
			
			if os.path.isfile('Bandpass_selfcal_record.npy'):
				os.system('rm -rf Bandpass_selfcal_record.npy')
			selfcal_record=np.array([num_iter,DR1,DR3,DR5,DR2,DR4,DR6,rms_list,calmode,scratch,antenna_list_index,start_sigma,\
				num_iter_fixed_sigma,num_iter_fixed_ant,stokes,startmodel,startmask,uvsub_flag_count],dtype='object')
			np.save('Bandpass_selfcal_record',selfcal_record)

			if verbose==False:
				print ('RMS based dynamic ranges:\n')
				print ('Negative based dynamic ranges:\n')
				print(str(DR2)+','+str(DR4)+','+str(DR6)+'\n')
				print ('Calmode : '+calmode+'\n')
				print ('Scartch = '+str(scratch)+'\n')
				print ('Sigma = '+str(start_sigma)+'\n')
			logger.info('RMS based dynamic ranges:\n')
			logger.info(str(DR1)+','+str(DR3)+','+str(DR5)+'\n')
			logger.info('Negative based dynamic ranges:\n')
			logger.info(str(DR2)+','+str(DR4)+','+str(DR6)+'\n')
			logger.info('Calmode : '+calmode+'\n')
			logger.info('Scratch = '+str(scratch)+'\n')
			logger.info('Sigma = '+str(start_sigma)+'\n')
		
			############## 
			# If statement 1 (DR decrease)

			if (((DR5<0.85*DR3 and DR5<0.9*DR1 and DR3>DR1) or (DR6<0.85*DR4 and DR6<0.9*DR2 and DR4>DR2)) and num_iter>min_iteration)\
				or (((DR5<0.9*DR3 and DR1>1.5*DR3) or (DR6<0.9*DR4 and DR2>1.5*DR4)) and num_iter>min_iteration):
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
				if DR5>min_DR and num_iter>min_iteration: # If minimum number of iteration is completed (Only considered the rms based DR here)
					os.system('rm -rf '+file_str+'.cal')
					os.system('rm -rf '+file_str+'_'+str(num_iter)+'.model')
					os.system('cp -r junk0.cal '+basedir+'/bpcaltables/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.cal')  # Keeping the last good caltable
					os.system('cp -r junk0.model '+basedir+'/bpimagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.model') # Keeping last good model
					os.system('cp -r junk0.image '+basedir+'/bpimagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.image') # Keeping last good model
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
							touch_file=basedir+'/.Finished_bcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_108'
							msg_str='Dear PAIRCARS user,\n\nBandpass self-calibration for : '+os.path.basename(msname)+'\nMessage : '+error_msgs(100)+', '+error_msgs(8)\
								+'\n\nBest regards,\nPAIRCARS developing team'
							msg_subject='Notification from PAIRCARS : Bandpass Selfcal : OBSID = '+str(OBSID)
							if inputs.send_notification==True:
								send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
								os.system('rm -rf '+quickimage)
							os.system('touch '+touch_file)
							if inputs.keep_logger==False:
								os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
							os.system('rm -rf '+working_dir+'/'+file_str+'*')
							end_time=time.time()
							run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
							logger.info('Total runtime : '+str(run_time))
						return 108
					else:
						os.chdir(cwd)
						if __name__!='__main__':
							touch_file=basedir+'/.Finished_bcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_8'
							msg_str='Dear PAIRCARS user,\n\nBandpass self-calibration for : '+\
								os.path.basename(msname)+'\nMessage : '+error_msgs(8)+'\n\nBest regards,\nPAIRCARS developing team'
							msg_subject='Notification from PAIRCARS : Bandpass Selfcal : OBSID = '+str(OBSID)
							if inputs.send_notification==True:
								send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
								os.system('rm -rf '+quickimage)
							os.system('touch '+touch_file)
							if inputs.keep_logger==False:
								os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
							os.system('rm -rf '+working_dir+'/'+file_str+'*')
							end_time=time.time()
							run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
							logger.info('Total runtime : '+str(run_time))
						return 8

			#######################################################
			# If statement 2 (Exiting selfcal conditions)

			if (DR5>=inputs.max_DR and num_iter>min_iteration):
				if verbose==False:
					print ('Reached limiting dynamic range\n')
				logger.info('Reached limiting dynamic range\n')
				end_selfcal=True
				os.system('cp -r junk1.cal '+basedir+'/bpcaltables/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.cal')
				os.system('cp -r junk1.model '+basedir+'/bpimagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.model')
				os.system('cp -r junk1.image '+basedir+'/bpimagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.image')
				if inputs.send_notification==True:
					quickimage=get_quicklook_image('junk1.image','quick_image_freq_'+freqstr+'_time_'+datestrfile+'.png',freqstr,datestr,DR5,DR6,field_of_view=2)
				os.system('rm -rf '+working_dir+'/junk*')
				os.chdir(cwd)
				if __name__!='__main__':
					touch_file=basedir+'/.Finished_bcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_0'
					msg_str='Dear PAIRCARS user,\n\nBandpass self-calibration for : '+\
						os.path.basename(msname)+'\nMessage : '+error_msgs(0)+'\n\nBest regards,\nPAIRCARS developing team'
					msg_subject='Notification from PAIRCARS : Bandpass Selfcal : OBSID = '+str(OBSID)
					if inputs.send_notification==True:
						send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
						os.system('rm -rf '+quickimage)
					os.system('touch '+touch_file)
					if inputs.keep_logger==False:
						os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
					os.system('rm -rf '+working_dir+'/'+file_str+'*')
					end_time=time.time()
					run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
					logger.info('Total runtime : '+str(run_time))
				return 0
			elif (abs(DR5-DR3)<DR_delta_rms and abs(DR5-DR1)<DR_delta_rms and abs(DR5/DR3-1)<0.08) and\
				 (abs(DR6-DR4)<DR_delta_neg and abs(DR6-DR2)<DR_delta_neg and abs(DR6/DR4-1)<0.05):
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
								loggr.info('do_uvsub_ankflag(\''+msname+'\',model=\'junk1.model\',nthread=1,verbose='+str(verbose)+',flagbackup=False)\n')
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
							os.system('cp -r junk1.cal '+basedir+'/bpcaltables/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.cal')
							os.system('cp -r junk1.model '+basedir+'/bpimagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.model')
							os.system('cp -r junk1.image '+basedir+'/bpimagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.image')
							if inputs.send_notification==True:
								quickimage=get_quicklook_image('junk1.image','quick_image_freq_'+freqstr+'_time_'+datestrfile+'.png',freqstr,datestr,DR5,DR6,field_of_view=2)
							os.system('rm -rf '+working_dir+'/junk*') 
							os.chdir(cwd)
							if __name__!='__main__':
								touch_file=basedir+'/.Finished_bcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_0'
								msg_str='Dear PAIRCARS user,\n\nBandpass self-calibration for : '+\
										os.path.basename(msname)+'\nMessage : '+error_msgs(0)+'\n\nBest regards,\nPAIRCARS developing team'
								msg_subject='Notification from PAIRCARS : Bandpass Selfcal : OBSID = '+str(OBSID)
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
								os.system('rm -rf '+working_dir+'/'+file_str+'*')
								end_time=time.time()
								run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
								logger.info('Total runtime : '+str(run_time))
							return 0
			
			#############################################################			
			# If statement 3 (Using last round model) 
			#(If DR increases at least DR_delta and all antennas are added and number of iteration at fixed antenna is greater than 5)
			
			if (DR5>DR3 and DR5>DR1) and (num_iter>min_iteration):
				startmodel='junk1.model'
			else:
				startmodel=''
			if (DR5>DR3 and DR5>DR1):
				startmask='junk1.mask'
			else:
				startmask=''
			num_iter+=1
			num_iter_fixed_sigma+=1
			num_iter_fixed_ant+=1
			
			###############################################################
			# If statement 4 (Reached maximum selfcal rounds)
	
			if num_iter>max_iteration:
				if DR5>min_DR:
					os.system('rm -rf '+file_str+'.cal')
					os.system('cp -r junk0.cal junk1.cal')  # Keeping the last good caltable
					os.system('cp -r junk0.ms junk1.ms')  # Keeping the last good calibrated ms
					if verbose==False:
						print (error_msgs(8))
					logger.error(error_msgs(8))
					end_selfcal=True
					os.system('cp -r junk1.cal '+basedir+'/bpcaltables/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.cal')
					os.system('cp -r junk1.model '+basedir+'/bpimagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.model')
					os.system('cp -r junk1.image '+basedir+'/bpimagemodels/'+str(OBSID)+'/'+basemsdir+'/'+file_str+'.image')
					if inputs.send_notification==True:
						quickimage=get_quicklook_image('junk1.image','quick_image_freq_'+freqstr+'_time_'+datestrfile+'.png',freqstr,datestr,DR5,DR6,field_of_view=2)
					os.system('rm -rf '+working_dir+'/junk*')
					if 'ref' in msname:
						os.chdir(cwd)
						if __name__!='__main__':
							touch_file=basedir+'/.Finished_bcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_109'
							msg_str='Dear PAIRCARS User,\n\nBandpass self-calibration for : '+os.path.basename(msname)+'\nMessage : '+error_msgs(100)+', '+error_msgs(9)\
										+'\n\nBest regards,\nPAIRCARS developing team'
							msg_subject='Notification from PAIRCARS : Bandpass Selfcal : OBSID = '+str(OBSID)
							if inputs.send_notification==True:
								send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
								os.system('rm -rf '+quickimage)
							os.system('touch '+touch_file)
							if inputs.keep_logger==False:
								os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
							os.system('rm -rf '+working_dir+'/'+file_str+'*')
							end_time=time.time()
							run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
							logger.info('Total runtime : '+str(run_time))
						return 109		
					else:
						os.chdir(cwd)
						if __name__!='__main__':
							touch_file=basedir+'/.Finished_bcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_9'
							msg_str='Dear PAIRCARS user,\n\nBandpass self-calibration for : '+\
								os.path.basename(msname)+'\nMessage : '+error_msgs(9)+'\n\nBest regards,\nPAIRCARS developing team'
							msg_subject='Notification from PAIRCARS : Bandpass Selfcal : OBSID = '+str(OBSID)
							if inputs.send_notification==True:
								send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
								os.system('rm -rf '+quickimage)
							os.system('touch '+touch_file)
							if inputs.keep_logger==False:
								os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
							os.system('rm -rf '+working_dir+'/'+file_str+'*')
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
						touch_file=basedir+'/.Finished_bcal_'+str(OBSID)+'_'+basemsdir+'_'+os.path.basename(msname)+'_13'
						msg_str='Dear PAIRCARS user,\n\nBandpass self-calibration for : '+\
							os.path.basename(msname)+'\nMessage : '+error_msgs(13)+'\n\nBest regards,\nPAIRCARS developing team'
						msg_subject='Notification from PAIRCARS : Intensity Selfcal : OBSID = '+str(OBSID)
						if inputs.send_notification==True:
							send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
							os.system('rm -rf '+quickimage)
						os.system('touch '+touch_file)
						if inputs.keep_logger==False:
							os.system('rm -rf '+working_dir+'/*.log '+working_dir+'/TempLattice*')
						os.system('rm -rf '+working_dir+'/'+file_str+'*')
						end_time=time.time()
						run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
						logger.info('Total runtime : '+str(run_time))
					return 13
						
# Function to run the script stand alone from command line
if __name__=='__main__':
	start_time=time.time()
	usage= ' Perform bandpass self calibration'
	parser = OptionParser(usage=usage)
	parser.add_option('--msname',dest="chantime_msname",default=None,help="Name of measurement set of a single time and frequency slice",metavar="Measurement Set")
	parser.add_option('--metafits',dest="metafits",default=None,help="Name of metafits file of the observation",metavar="Metafits file")
	parser.add_option('--workdir',dest='workdir',default=None,help='Name of the working directory',metavar='Directory path')
	parser.add_option('--verbose',dest="verbose",default=False,help="Verbose mode",metavar="Boolean")
	parser.add_option('--interactive',dest="interactive",default=False,help="Interactive mode",metavar="Boolean")
	parser.add_option('--fresh',dest="fresh",default=True,help="Start fresh self calibration loop",metavar="Boolean")
	(options, args) = parser.parse_args()
	if (os.path.isfile(str(options.workdir)+'/Bandpass_Selfcal.log') and eval(str(options.fresh))==True) or \
		(os.path.isfile(str(options.workdir)+'/Bandpass_Selfcal.log') and os.path.isdir(str(options.workdir)+'/junk1.ms')==False and eval(str(options.fresh))==False):
		os.system('rm -rf '+str(options.workdir)+'/Bandpass_Selfcal.log')
	if os.path.isfile(str(options.workdir)+'/Bandpass_Selfcal_verbose.log') and eval(str(options.verbose))==True:
			os.system('rm -rf '+str(options.workdir)+'/Bandpass_Selfcal_verbose.log')
	formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
	logger = logging.getLogger('bandpass_selfcal_log')
	logger.setLevel(logging.DEBUG)
	if eval(str(options.verbose))==True:
		console=logging.StreamHandler(sys.stdout)
		console.setFormatter(formatter)
		logger.addHandler(console)
	filehandle=logging.FileHandler(str(options.workdir)+'/Bandpass_Selfcal.log')
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
		touch_file=inputs.basedir+'/.Finished_bcal_'+str(OBSID)+'_'+basemsdir+'_'+msbasename+'_'+str('noms')
		end_time=time.time()
		run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
		logger.info('#############################\n')
		logger.info('Bandpass selfcal failed for ms : '+options.chantime_msname+'\n')
		logger.info('Total runtime : '+str(run_time)+'\n')
		logger.info('##############################\n')
		msg_str='Dear PAIRCARS user,\n\nBandpass self-calibration for : '+msbasename+'\nMessage : No measurement set is present\nTotal runtime : '+\
					str(run_time)+'\n\nBest regards,\nPAIRCARS developing team'
		msg_subject='Notification from PAIRCARS : Bandpass Selfcal : OBSID = '+str(OBSID)
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
		touch_file=inputs.basedir+'/.Finished_bcal_'+str(OBSID)+'_'+basemsdir+'_'+msbasename+'_'+str('nometa')
		end_time=time.time()
		run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
		logger.info('#############################\n')
		logger.info('Bandpass selfcal failed for ms : '+options.chantime_msname+'\n')
		logger.info('Total runtime : '+str(run_time)+'\n')
		logger.info('##############################\n')
		msg_str='Dear PAIRCARS user,\n\nBandpass self-calibration for : '+msbasename+'\nMessage : No metafits file is present\nTotal runtime : '+\
					str(run_time)+'\n\nBest regards,\nPAIRCARS developing team'
		msg_subject='Notification from PAIRCARS : Bandpass Selfcal : OBSID = '+str(OBSID)
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
		print ('\n\t##########################\n\tStarting Bandpass self-calibration.....\n\t##########################\n')
		print ('run_bandpass_selfcal(\''+options.chantime_msname+'\',\''+options.metafits+'\',\''+options.workdir+'\',verbose=\''+str(options.verbose)\
				+'\',interactive=\''+str(options.interactive)+'\',start_fresh=\''+str(options.fresh)+'\')\n')
		msg=run_bandpass_selfcal(options.chantime_msname,options.metafits,options.workdir,verbose=eval(str(options.verbose)),\
				interactive=eval(str(options.interactive)),start_fresh=eval(str(options.fresh)))
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
		touch_file=inputs.basedir+'/.Finished_bcal_'+str(OBSID)+'_'+basemsdir+'_'+msbasename+'_'+str(msg)
		end_time=time.time()
		run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
		logger.info('#############################\n')
		logger.info('Bandpass selfcal finished for ms : '+options.chantime_msname+'\n')
		logger.info('Total runtime : '+str(run_time)+'\n')
		logger.info('##############################\n')
		msg_str='Dear PAIRCARS user,\n\nBandpass self-calibration for : '+msbasename+'\n'+msg_str+'\nTotal runtime : '+str(run_time)+'\n\nBest regards,\nPAIRCARS developing team'
		msg_subject='Notification from PAIRCARS : Bandpass Selfcal : OBSID = '+str(OBSID)
		if inputs.send_notification==True:
			attachments=glob.glob(options.workdir+'/quick_image_*.png')
			send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=attachments)
			os.system('rm -rf '+options.workdir+'/quick_image_*.png')
		os.system('touch '+touch_file)
		if inputs.keep_logger==False:
			os.system('rm -rf '+options.workdir+'/*.log '+options.workdir+'/TempLattice*')
		file_str=msbasename.split('.ms')[0]
		os.system('rm -rf '+options.workdir+'/'+file_str+'*')
	except Exception as e:
		touch_file=inputs.basedir+'/.Finished_bcal_'+str(OBSID)+'_'+basemsdir+'_'+msbasename+'_'+str('error')
		end_time=time.time()
		run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
		logger.info('#############################\n')
		logger.info('Bandpass selfcal failed for ms : '+options.chantime_msname+'\n')
		logger.info('Error occured : '+str(e)+'\n')
		logger.info('Total runtime : '+str(run_time)+'\n')
		logger.info('##############################\n')
		msg_str='Dear PAIRCARS user,\n\nBandpass self-calibration for : '+msbasename+'\nMessage : Error in runtime : '+str(e)+'\nTotal runtime : '+\
					str(run_time)+'\n\nBest regards,\nPAIRCARS developing team'
		msg_subject='Notification from PAIRCARS : Bandpass Selfcal : OBSID = '+str(OBSID)
		if inputs.send_notification==True:
			attachments=glob.glob(options.workdir+'/quick_image_*.png')
			send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=attachments)
			os.system('rm -rf '+options.workdir+'/quick_image_*.png')
		os.system('touch '+touch_file)
		if inputs.keep_logger==False:
			os.system('rm -rf '+options.workdir+'/*.log '+options.workdir+'/TempLattice*')
		file_str=msbasename.split('.ms')[0]
		os.system('rm -rf '+options.workdir+'/'+file_str+'*')
		pass
