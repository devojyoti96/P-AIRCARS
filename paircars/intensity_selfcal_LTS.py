import numpy as np,os,sys,glob,copy
import logging
from casatools import *
from casatasks import *
from . import access_ms as am
from . import basic_func as B
from . import flagger as fg
from paircars_casatasks.poltclean import *
from astropy.io import fits
from astropy.wcs import WCS

# PyBDSF installation has problem
bdsf_import=False

'''
Code is written by Devojyoti Kansabanik, 17 Jan, 2021
'''

class IntensitySelfcal:
	'''
	Generic class to perform intensity self-calibration
	Attributes:
	msname = Name of the measurement set
	maximum_emission_scale = Maximum scale of the emission present in the image in arcsec
	verbose = False,If True keep all the intermediate images, model, residuals, caltables and details of the log to detailed analysis
	interactive = False, If True user have interactive control on self-calibration
	'''
	def __init__(self,msname,maximum_emission_scale,verbose=False,interactive=False):
		self.cwd=os.getcwd()
		if msname[-1]=='/':
			self.msname=msname[:-1]
		else:
			self.msname=msname
		self.mspath=os.path.dirname(os.path.realpath(msname))
		AM=am.AccessMS(self.msname)
		IB=B.ImageBasic(self.msname)
		self.max_baseline=AM.get_max_baseline()	
		self.cellsize=IB.calc_cellsize(3) # Assuming 3 pixels in one PSF
		self.imsize=IB.num_pixels(3)
		self.max_size=maximum_emission_scale
		self.multiscale_scales=IB.choose_scales(3,self.max_size)
		self.uvtaper=IB.calc_uvtaper()
		self.calib_uvrange=IB.calc_calib_uvrange(4)[0] # Short baselines sensitive to larger than 4 deg are excluded 
		self.rms_box='50,50,'+str(self.imsize-50)+','+str(int(self.imsize/4)) # CASA box to calculate the rms
		self.verbose=verbose
		self.interactive=interactive
		formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
		self.log_verbose = logging.getLogger('selfcal_verbose_log')
		self.log_verbose.setLevel(logging.DEBUG)
		if self.verbose:
			self.console=logging.StreamHandler(sys.stdout)
			self.console.setFormatter(formatter)
			self.log_verbose.addHandler(self.console)
		self.filehandle=logging.FileHandler(self.cwd+'/Intensity_Selfcal_verbose.log')
		self.filehandle.setFormatter(formatter)
		self.log_verbose.addHandler(self.filehandle)
		self.log_verbose.propagate = False
		self.log_verbose.info('Initiating Intensity selfcal object.\n')
		os.system('touch '+self.msname+'/.usedby_paircars')

	def negative_box(self,max_pix,box_width=3):
		'''
		Create a box about the maximum pixel of image to search negative.
		Parameters:
		max_pix= Maximum pixel [xxmax,yymax]
		box_width = Box width in degree (default : 3 degree)
		Return:
		CASA box 'xblc,yblc,xrtc,yrtc'
		'''
		max_pix_xx=max_pix[0]
		max_pix_yy=max_pix[1]
		box_length=(float(box_width)*3600.0)/self.cellsize # Taking a box about the max pixel to serach to minimum
		xblc=max_pix_xx-(box_length/2.0)
		xrtc=max_pix_xx+(box_length/2.0)
		yblc=max_pix_yy-(box_length/2.0)
		yrtc=max_pix_yy+(box_length/2.0)
		os.system('rm -rf casa*log')
		return str(int(xblc))+','+str(int(yblc))+','+str(int(xrtc))+','+str(int(yrtc))

	def calc_dyn_range(self,num_iter,sigma,box_width=3,stokes_list=['I']):
		'''
		Calculate the dynamic range of the full Stokes image cube
		Parameters:
		num_iter = Number of selfcal iteration
		sigma = nsigma value to put a mask for calculating total flux
		box_width = Negative box width around the maximum pixel in degree (default : 3 degree)
		stokes_list = ['I','Q','U','V','XX','YY'] list of stokes planes in the image
		Return:	
		Python dictionary {'STOKES':[rms dynamic range,rms,total_flux(non-negative)]},negative dynamic range for Stokes I
		'''
		file_str=self.msname.split('.ms')[0]+'_'+str(num_iter) # File string prefix
		ia=image()
		imagename=file_str+'.image'
		if os.path.isdir(imagename):
			imageheader=imhead(imagename=imagename,mode='summary')
			if imageheader['ndim']>2 and imageheader['ndim']==4:
				out_dict={}
				neg_dyn=0
				for stokes in stokes_list:
					if stokes=='I' or stokes=='XX' or stokes =='YY':
						maxpos=imstat(imagename=imagename,stokes=stokes)['maxpos']
						negative_box=self.negative_box(maxpos,box_width=box_width)
						max_pix=imstat(imagename=imagename,stokes=stokes)['max'][0]
						rms=imstat(imagename=imagename,box=self.rms_box,stokes=stokes)['rms'][0]
						min_pix=imstat(imagename=imagename,box=negative_box,stokes=stokes)['min'][0]
						rms_dyn_range=max_pix/rms
						if min_pix!=0:
							neg_dyn+=max_pix/abs(min_pix)
						else:
							neg_dyn=rms_dyn_range
						ia.open(imagename)
						ia.calcmask('\"'+imagename+'\">'+str(sigma*rms),'mymask')
						ia.close()
						try:
							total_flux=imstat(imagename=imagename,stokes=stokes)['flux'][0]
							out_dict[stokes]=[rms_dyn_range,rms,total_flux]							
						except:
							out_dict[stokes]=[rms_dyn_range,rms,np.nan]
						makemask(mode='delete',inpmask=imagename+':mymask')
					else:
						max_pix=imstat(imagename=imagename,stokes=stokes)['max'][0]
						min_pix=imstat(imagename=imagename,stokes=stokes)['min'][0]
						rms=imstat(imagename=imagename,box=self.rms_box,stokes=stokes)['rms'][0]
						if abs(min_pix)>max_pix:
							max_pix=abs(min_pix)
						rms_dyn_range=max_pix/rms
						if os.path.isdir(stokes+'.image'):
							os.system('rm -rf '+stokes+'.image')
						immath(imagename=imagename,outfile=stokes+'.image',mode='evalexpr',expr='abs(IM0)',stokes=stokes)
						ia.open(stokes+'.image')
						ia.calcmask(stokes+'.image'+'>'+str(sigma*rms),'mymask')
						ia.close()
						try:
							total_flux=imstat(imagename=stokes+'.image',stokes=stokes)['flux'][0]
							out_dict[stokes]=[rms_dyn_range,rms,total_flux]
						except:
							out_dict[stokes]=[rms_dyn_range,rms,np.nan]	
						os.system('rm -rf '+stokes+'.image')			
		else:
			out_dict={}
			rms_dyn_range=np.nan
			neg_dyn=np.nan
			out_dict['NAN']=[np.nan,np.nan]
		negative_dyn_range=neg_dyn/len(stokes_list)
		os.system('rm -rf casa*log')
		return out_dict,negative_dyn_range

	def change_start_sigma(self,nsigma,sigma_step,min_sigma):
		'''
		Function to calculate start sigma
		Parameters:
		nsigma = Present sigma value
		sigma_step = Sigma reduction step
		min_sigma = Minimum value of sigma to reduce
		Return:
		Value of the start sigma
		'''
		if nsigma-sigma_step > min_sigma: # If initial CLEAN could not pick up any flux at start sigma, we lower down start sigma upto minsigma
			if sigma_step>1.0:
				print ('WARNING : Choosing sigma step greater than 1 is too risky. Selfcal may diverge\n')
				if self.interactive==True:
					want_to_continue=input('Do you still want to continue? Y/y/Yes/yes')
					self.log_verbose.warning('Choosing sigma step 1 is too risky. Selfcal may diverge\n')
					self.log_verbose.info('Interactive=True\n')
					self.log_verbose.info('Continue? Y/y/Yes/yes:'+str(want_to_continue)+'\n')
					if want_to_continue=='Y' or want_to_continue=='y' or want_to_continue=='Yes' or want_to_continue=='yes':		
						nsigma-=sigma_step
						self.log_verbose.info('Continuing with sigma step :'+str(sigma_step)+'\n')
				else:
					self.log_verbose.info('Interactive=False\n')
					self.log_verbose.info('Continuing with sigma step :'+str(sigma_step)+'\n')
					if self.verbose==False:
						print ('Continuing with sigma step :'+str(sigma_step)+'\n')
					nsigma-=sigma_step
			else:
				nsigma-=sigma_step
		else:
			if self.interactive==True:
				print ('WARNING : Start sigma is below '+str(min_sigma)+'. It is too risky. Initial model can pickup spurious emissions in the model and selfcal may not progress\n')
				want_to_continue=input('Do you still want to continue? Y/y/Yes/yes')
				self.log_verbose.warning('Start sigma is below '+str(min_sigma)+\
							'. It is too risky. Initial model can pickup spurious emissions in the model and selfcal may not progress\n')
				self.log_verbose.info('Interactive=True\n')
				self.log_verbose.info('Continue? Y/y/Yes/yes:'+str(want_to_continue)+'\n')
				if want_to_continue=='Y' or want_to_continue=='y' or want_to_continue=='Yes' or want_to_continue=='yes':
					nsigma-=sigma_step
					self.log_verbose.info('Continuing with sigma :'+str(nsigma)+'\n')
			else:
				self.log_verbose.info('Interactive=False\n')
				self.log_verbose.info('Stop here, since start sigma is lower than '+str(min_sigma)+'\n')
				nsigma=np.nan
		os.system('rm -rf casa*log')
		return nsigma
		
	def reduce_sigma(self,imagename,nsigma,sigma_step,minsigma,residual_frac=0.1,stokes_list=['I']):
		'''
		Function to determine whether reduce the CLEAN sigma or not.
		Parameters:
		imagename = Name of the image
		nsigma = Value of the present n-sigma
		sigma_step = Step to reduce sigma value
		minsigma = Minimum allowed sigma
		residual_frac = Residual flux fraction to reduce sigma
		stokes_list = ['I'], stokes plane list
		Return:
		Reduced value of n-sigma if residual flux is more than given percentage (default : 10%) of the total flux in Stokes I or in all Stokes Q,U,V.
		'''
		imagename=imagename
		residual=imagename.split('.image')[0]+'.residual'
		do_reduce_list=[]
		ia=image()
		imagename_path=os.path.dirname(os.path.realpath(imagename))
		cwd=os.getcwd()
		if imagename_path!='':
			os.chdir(imagename_path)
		os.system('rm -rf reduce_sigma_*')
		for stokes in stokes_list:
			self.log_verbose.info('imstat(imagename=\''+imagename+'\',box=\''+self.rms_box+'\',stokes=\''+stokes+'\')[\'rms\'][0]\n')
			rms=imstat(imagename=imagename,box=self.rms_box,stokes=stokes)['rms'][0]
			os.system('rm -rf reduce_sigma_*')
			immath(imagename=imagename,mode='evalexpr',expr='abs(IM0)',stokes=stokes,outfile='reduce_sigma_'+stokes+'.image')
			immath(imagename=residual,mode='evalexpr',expr='abs(IM0)',stokes=stokes,outfile='reduce_sigma_'+stokes+'.residual')
			ia.open('reduce_sigma_'+stokes+'.image')			
			ia.calcmask('\"reduce_sigma_'+stokes+'.image\">'+str(nsigma*rms),'mymask')
			ia.close()
			makemask(inpimage='reduce_sigma_'+stokes+'.image',inpmask='reduce_sigma_'+stokes+'.image:mymask',output='reduce_sigma_'+stokes+'.residual:mymask',mode='copy')
			try:
				image_pix_sum=imstat(imagename='reduce_sigma_'+stokes+'.image')['sum'][0]
				residual_pix_sum=imstat(imagename='reduce_sigma_'+stokes+'.residual')['sum'][0]
			except:
				image_pix_sum=0
				residual_pix_sum=1
			if residual_pix_sum/image_pix_sum>residual_frac:
				do_reduce_list.append(1)
		os.system('rm -rf reduce_sigma_*')
		os.chdir(cwd)
		if int(np.sum(np.array(do_reduce_list)))>=2:
			if sigma_step>1.0:
				self.log_verbose.info('WARNING : Choosing sigma step 1 is too risky. Selfcal may diverge\n')
				if self.verbose==False:
					print ('WARNING : Choosing sigma step greater than 1 is too risky. Selfcal may diverge\n')
				if self.interactive==True:
					want_to_continue=input('Do you still want to continue? Y/y/Yes/yes')
					self.log_verbose.info('Interactive=True\n')
					self.log_verbose.info('Do you still want to continue? Y/y/Yes/yes:'+str(want_to_continue)+'\n')
					if want_to_continue=='Y' or want_to_continue=='y' or want_to_continue=='Yes' or want_to_continue=='yes':	
						self.log_verbose.info('Continuing with sigma step :'+str(sigma_step)+'\n')
						os.system('rm -rf casa*log')
						return nsigma-sigma_step
					else:
						os.system('rm -rf casa*log')
						return nsigma
			elif nsigma-sigma_step<minsigma:
				self.log_verbose.info('WARNING : Choosing sigma less than '+str(minsigma)+'\n')
				if self.verbose==False:
					print ('WARNING : Choosing sigma less than '+str(minsigma)+'\n')
				if self.interactive==True:
					want_to_continue=input('Do you still want to continue? Y/y/Yes/yes')
					if self.verbose:
						self.log_verbose.info('Interactive=True\n')
						self.log_verbose.info('Do you still want to continue? Y/y/Yes/yes:'+str(want_to_continue)+'\n')
					if want_to_continue=='Y' or want_to_continue=='y' or want_to_continue=='Yes' or want_to_continue=='yes':		
						self.log_verbose.info('Continuing with sigma step :'+str(sigma_step)+'\n')
						os.system('rm -rf casa*log')
						return nsigma-sigma_step
					else:
						os.system('rm -rf casa*log')
						return nsigma
				else:
					self.log_verbose.info('Interactive=False\n')
					self.log_verbose.info('Continuing with sigma step :'+str(sigma_step)+'\n')
					if self.verbose==False:
						print ('Continuing with sigma step :'+str(sigma_step)+'\n')
					os.system('rm -rf casa*log')
					return nsigma-sigma_step
			else:
				self.log_verbose.info('Reducing sigma to:'+str(nsigma-sigma_step)+'\n')
				os.system('rm -rf casa*log')
				return nsigma-sigma_step
		else:
			self.log_verbose.info('Sigma value is not changed. Sigma is at :'+str(nsigma)+'\n')
			os.system('rm -rf casa*log')
			return nsigma
	
	def calc_iter_num(self,safety_factor,quality_factor,scratch=True,bandpass_selfcal=False):
		'''
		Function to calculate minimum number of selfcal iteration based on safety standard and quality factor
		Parameters:
		safety_factor = Factor to determine the robustness of the selfcal
		quality_factor = Factor to determine the quality of the images
		scratch : True, whether start the selfcal from scratch or not
		bandpass_selfcal = False, performinh bandpass selfcal or not
		Return:
		Minimum iteration at fixed sigma, Minimum iteration, Maximum iteration , Number of antenna bins
		'''
		if bandpass_selfcal==True:
			min_num_iter_fixed_sigma=1
			min_iteration=1
			max_iteration=100
			antenna_bin=1
		else:
			if quality_factor==0:     # Low quality (Quick look image making)
				if (safety_factor==0):
					min_num_iter_fixed_sigma=0
					if (scratch==True):
						min_iteration=3
						max_iteration=50
						antenna_bin=2
					else:
						min_iteration=3
						max_iteration=10
						antenna_bin=1
				elif (safety_factor==1):
					min_num_iter_fixed_sigma=0
					if (scratch==True):
						min_iteration=10
						max_iteration=70
						antenna_bin=3
					else:
						min_iteration=5
						max_iteration=30
						antenna_bin=1
				else:
					min_num_iter_fixed_sigma=0
					if (scratch==True):
						min_iteration=20
						max_iteration=100
						antenna_bin=5
					else:
						min_iteration=10
						max_iteration=60
						antenna_bin=1
			elif quality_factor==1:  # Medium quality imaging (Computing speed medium)
				if (safety_factor==0):
					min_num_iter_fixed_sigma=0
					if (scratch==True):
						min_iteration=5
						max_iteration=200
						antenna_bin=3
					else:
						min_iteration=3
						max_iteration=150
						antenna_bin=1
				elif (safety_factor==1):
					min_num_iter_fixed_sigma=5
					if (scratch==True):
						min_iteration=10
						max_iteration=400
						antenna_bin=5
					else:
						min_iteration=5
						max_iteration=350
						antenna_bin=1
				else:
					min_num_iter_fixed_sigma=10
					if (scratch==True):
						min_iteration=20
						max_iteration=600
						antenna_bin=7
					else:
						min_iteration=10
						max_iteration=550
						antenna_bin=1
			else:  # Best quality imaging (Computing slow)
				if (safety_factor==0):
					max_iteration=700
					min_num_iter_fixed_sigma=0
					if (scratch==True):
						min_iteration=10
						max_iteration=600
						antenna_bin=5
					else:
						max_iteration=550
						min_iteration=5
						antenna_bin=1
				elif (safety_factor==1):
					min_num_iter_fixed_sigma=5
					if (scratch==True):
						max_iteration=650
						min_iteration=20
						antenna_bin=7
					else:
						max_iteration=600
						min_iteration=10
						antenna_bin=1
				else:
					min_num_iter_fixed_sigma=10
					if (scratch==True):
						max_iteration=700
						min_iteration=30
						antenna_bin=9
					else:
						max_iteration=650
						min_iteration=20
						antenna_bin=1
		self.log_verbose.info('Quality factor : '+str(quality_factor)+', Safety standard : '+str(safety_factor)+', Scratch : '+str(scratch)+\
				', Minimum number of iteration at fixed sigma : '+str(min_num_iter_fixed_sigma)+', Minimum iteration : '+str(min_iteration)+', Antenna bins : '+str(antenna_bin)+'\n')
		os.system('rm -rf casa*log')
		return min_num_iter_fixed_sigma,min_iteration,max_iteration,antenna_bin

	def antenna_string(self,antenna_list,antenna_list_index):
		'''
		Function to return antenna string from antenna list
		Parameters:
		antenna_list = Antenna list or array
		antenna_list_index = Bin number of antenna list
		Return:
		Antenna string
		'''
		antenna_string=''
		for ant in antenna_list[antenna_list_index]:
			antenna_string+=str(ant)+','
		antenna_string=antenna_string[:-1]
		os.system('rm -rf casa*log')
		return antenna_string
		
	def DR_record(self,DR,record_filename,init=True):
		'''
		Function to keep the record of dynamic ranges at different self calibration steps
		Parameters:
		DR = Dynamic range of the current selfcal iteration to record
		record_filename = Name of the file to stro dynamic ranges
		init = True, initiating a new record from the current selfcal iteration
		Return:
		Dynamic range record array
		'''
		if init==True:
			os.system('rm -rf '+record_filename+'.npy')
			DR_array=np.empty([1])
			DR_array[0]=DR
			np.save(record_filename,DR_array)
		else:
			if os.path.isfile(record_filename+'.npy')==False:
				DR_array=np.empty([1])
				DR_array[0]=DR
			else:
				DR_array=np.load(record_filename+'.npy')
				DR_array=np.append(DR_array,DR)
			np.save(record_filename,DR_array)
		os.system('rm -rf casa*log')
		return DR_array
			
	def DR_delta(self,quality_factor,DR_array,frac_increase,safety_factor): # Not tested and not used
		'''
		Function to calculate the dynamic range increse step
		Parameters:
		quality_factor = Factor to determine the image quality
		DR_array = Dynamic range records array
		frac_increase = Fractional increase of the DR step from previous DR steps
		safety-factor = Factor to determine the robustness of the selfcal
		Return:
		Dynamic range increase step
		'''
		if quality_factor==0 and safety_factor==0 and len(DR_array)<3:
			DR_delta_step=20
		else:
			DR_delta_step=np.polyfit(DR_array,1)[0]*frac_increase 
		os.system('rm -rf casa*log')
		return DR_delta_step
			
	def calc_selfcal_snr(self,imagename,sigma,num_antenna_to_use):
		'''
		Function to calculate the SNR of the present iteration for selfcal
		Parameters:
		imagename = Name of the image
		sigma = Present sigma value used for making the image
		num_antenna_to_use = Number of antenna used for imaging
		Return:
		SNR of the image for selfcal
		'''
		if os.path.isdir(imagename):
			try:
				rms=imstat(imagename=imagename,box=self.rms_box)['rms'][0]
				ia=image()
				ia.open(imagename)
				ia.calcmask('\''+imagename+'\'>'+str(sigma*rms),'selfcal_snr_mask')
				ia.close()
				flux=imstat(imagename)['flux'][0]
				makemask(inpmask=imagename+':selfcal_snr_mask',mode='delete')
				selfcal_snr=flux/(rms*np.sqrt(num_antenna_to_use-3))
				os.system('rm -rf casa*log')
				return selfcal_snr
			except:
				os.system('rm -rf casa*log')
				return np.nan
		else:
			os.system('rm -rf casa*log')
			return np.nan

	def dirty_image(self,start_sigma,box_width=3,antenna_to_use=''):
		'''
		Function make dirty image and return dynamic range and selfcal SNR.
		Parameters:
		start_sigma = Start sigma
		box_width = Negative box width around the maximum pixel in degree (default : 3 degree)
		antenna_to_use = List of antennas
		Return:
		rms based dynamic range,negative based dynamic range,selfcal SNR,error message code
		'''
		imagename=self.msname.split('.ms')[0]+'_dirty'
		present_file=glob.glob(imagename+'*')
		if len(present_file)!=0:
			os.system('rm -rf '+imagename+'*')
			self.log_verbose.info('rm -rf '+imagename+'*\n')
		self.log_verbose.info('Making dirty image...........\n')
		self.log_verbose.info('tclean(vis=\''+self.msname+'\',imagename=\''+imagename+'\',imsize=['+str(self.imsize)+'],cell=\''+\
								str(self.cellsize)+'\',niter=0,antenna=\''+antenna_to_use+'\',uvtaper=\''+str(self.uvtaper)+'\',weighting=\'natural\')\n')
		tclean(vis=self.msname,imagename=imagename,imsize=[self.imsize],cell=self.cellsize,niter=0,antenna=antenna_to_use,uvtaper=self.uvtaper,weighting='natural')
		out_dict,negative_dyn_range=self.calc_dyn_range('dirty',start_sigma,box_width=box_width,stokes_list=['I']) # Calculating the dynamic range of the image
		selfcal_snr=self.calc_selfcal_snr(imagename+'.image',start_sigma,len(antenna_to_use.split(',')))
		if self.verbose==False:
			os.system('rm -rf '+imagename+'.image '+imagename+'.model '+imagename+'.residual '+imagename+'.sumwt '+imagename+'.pb '+imagename+'.psf '+imagename+'.mask ')
		out_dict_keys=out_dict.keys()
		if 'NAN' in out_dict_keys:
			os.system('rm -rf casa*log')
			os.chdir(self.cwd)
			os.system('rm -rf casa*log')
			return 2,out_dict,negative_dyn_range,selfcal_snr	
		else:
			os.system('rm -rf casa*log')
			return 0,out_dict,negative_dyn_range,selfcal_snr
		
	def file_remover_and_keeper(self,num_iter,msg_code,do_bandpass=False,ref_time_chan=True):
		'''
		This function keep and remove caltables, ms, imaging related files based on the need
		Parameters:
		num_iter = Number of self-calibration iteration
		msg_code = Selfcal message code
		do_bandpass = False, performing bandpass or not
		ref_timechan = True , reference time channel or not
		'''
		msname_str=am.splited_ms_rename(self.msname,ref_time_chan=ref_time_chan,change_msname=False)
		freqstr=os.path.basename(msname_str).split('.ms')[0].split('_freq_')[1].split('_')[0]  # Frequency string in MHz
		datestr_list=os.path.basename(msname_str).split('.ms')[0].split('_freq_')[0].split('time_')[1].split('_')
		datestr_for_file='_'.join(datestr_list[:3])+'_'+'_'.join(datestr_list[3:]) # Datetime string 
		cwd=os.getcwd()
		file_str=os.path.basename(self.msname).split('.ms')[0]+'_'+str(num_iter) # File string prefix
		caltable_name=self.msname.split('.ms')[0]+'.cal' # Caltable name
		if do_bandpass==True:
			file_str_prefix='freq_'+freqstr+'_datetime_'+datestr_for_file+'_bp'
		else:
			file_str_prefix='freq_'+freqstr+'_datetime_'+datestr_for_file
		if self.verbose==True and os.path.isdir(self.mspath+'/'+file_str_prefix)==False: 
						# If verbose=True, directory to keep all intermediate images, caltables, models, residuals 
			os.mkdir(self.mspath+'/'+file_str_prefix)
		if self.verbose:
			if os.path.isdir(self.mspath+'/'+file_str_prefix+'/backup_ms')==False:
				os.makedirs(self.mspath+'/'+file_str_prefix+'/backup_ms')
			if os.path.isdir(self.mspath+'/'+file_str_prefix+'/backup_imagemodel')==False:
				os.makedirs(self.mspath+'/'+file_str_prefix+'/backup_imagemodel')
			if os.path.isdir(self.mspath+'/'+file_str_prefix+'/backup_cal')==False:
				os.makedirs(self.mspath+'/'+file_str_prefix+'/backup_cal')
		if num_iter=='dirty':
			if self.verbose and (msg_code==0):
				os.system('cp -r '+file_str+'.model '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.model') # If Verbose=True, keep all the models
				os.system('cp -r '+file_str+'.image '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.image') # If Verbose=True, keep all the image
			os.system('rm -rf '+file_str+'.image '+file_str+'.model '+file_str+'.residual '+file_str+'.sumwt '+file_str+'.pb '+file_str+'.psf') 
			os.chdir(self.cwd)
			os.system('rm -rf casa*log')
			return
		if num_iter==0 and (msg_code==0 or msg_code==8 or msg_code==9):
			os.system('cp -r '+caltable_name+' junk0.cal') # Copying num_iter=0 caltable to junk0.cal
			os.system('cp -r '+self.msname+' junk0.ms') # Copying num_iter=0 ms to junk0.ms
			os.system('cp -r '+file_str+'.image junk0.image') # Copying num_iter=0 image to junk0.image
			os.system('cp -r '+file_str+'.model junk0.model') # Copying num_iter=0 model to junk0.model
			os.system('cp -r '+file_str+'.mask junk0.mask') # Copying num_iter=0 mask to junk0.mask
			os.system('cp -r '+file_str+'.residual junk0.residual') # Copying num_iter=0 residual to junk0.residual
			if self.verbose and (msg_code==0 or msg_code==8 or msg_code==9):	
				os.system('cp -r '+caltable_name+' '+self.mspath+'/'+file_str_prefix+'/backup_cal/'+file_str+'.cal') # Verbose=True, keep all the caltables
				os.system('cp -r '+self.msname+' '+self.mspath+'/'+file_str_prefix+'/backup_ms/'+file_str+'.ms') # If Verbose=True, keep all the ms
				os.system('cp -r '+file_str+'.model '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.model') # If Verbose=True, keep all the models
				os.system('cp -r '+file_str+'.mask '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.mask') # If Verbose=True, keep all the masks
				os.system('cp -r '+file_str+'.image '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.image') # If Verbose=True, keep all the image
				os.system('cp -r '+file_str+'.residual '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.residual') # If Verbose=True, keep all the residuals
			os.system('rm -rf '+caltable_name)
			os.system('rm -rf '+file_str+'.image '+file_str+'.model '+file_str+'.residual '+file_str+'.sumwt '+file_str+'.pb '+file_str+'.psf '+file_str+'.mask ') 
				# Removing all imaging related files
		elif num_iter==1 and (msg_code==0 or msg_code==8 or msg_code==9):
			os.system('cp -r '+caltable_name+' junk1.cal') # Copying num_iter=1 caltable to junk1.cal
			os.system('cp -r '+self.msname+' junk1.ms') # Copying num_iter=1 ms to junk1.ms
			os.system('cp -r '+file_str+'.image junk1.image') # Copying num_iter=1 image to junk1.image
			os.system('cp -r '+file_str+'.model junk1.model') # Copying num_iter=1 model to junk1.model
			os.system('cp -r '+file_str+'.mask junk1.mask') # Copying num_iter=1 model to junk1.mask
			os.system('cp -r '+file_str+'.residual junk1.residual') # Copying num_iter=1 residual to junk1.residual
			if self.verbose and (msg_code==0 or msg_code==8 or msg_code==9):	
				os.system('cp -r '+caltable_name+' '+self.mspath+'/'+file_str_prefix+'/backup_cal/'+file_str+'.cal') # Verbose=True, keep all the caltables
				os.system('cp -r '+self.msname+' '+self.mspath+'/'+file_str_prefix+'/backup_ms/'+file_str+'.ms') # If Verbose=True, keep all the ms
				os.system('cp -r '+file_str+'.model '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.model') # If Verbose=True, keep all the models
				os.system('cp -r '+file_str+'.mask '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.mask') # If Verbose=True, keep all the masks
				os.system('cp -r '+file_str+'.image '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.image') # If Verbose=True, keep all the image
				os.system('cp -r '+file_str+'.residual '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.residual') # If Verbose=True, keep all the residuals
			os.system('rm -rf '+caltable_name)
			os.system('rm -rf '+file_str+'.image '+file_str+'.model '+file_str+'.residual '+file_str+'.sumwt '+file_str+'.pb '+file_str+'.psf '+file_str+'.mask ') 
				# Removing all imaging related files
		elif num_iter>1 and (msg_code==0 or msg_code==8 or msg_code==9):
			if os.path.isdir('junk1.cal'):
				os.system('rm -rf junk0.cal')
				os.system('cp -r junk1.cal junk0.cal') # Move the previous round caltable to pre-previous round for num_iter>1
			else:
				os.system('cp -r '+caltable_name+' junk0.cal') # Copying caltable to junk0.cal
			if os.path.isdir('junk1.ms'):
				os.system('rm -rf junk0.ms')
				os.system('cp -r junk1.ms junk0.ms')	# Move the previous round ms to pre-previous round for num_iter>1
			else:
				os.system('cp -r '+self.msname+' junk0.ms') # Copying ms to junk0.ms
			if os.path.isdir('junk1.model'):
				os.system('rm -rf junk0.model')
				os.system('cp -r junk1.model junk0.model') # Move the previous round model to pre-previous round for num_iter>1
			else:
				os.system('cp -r '+file_str+'.model junk0.model') # Copying model to junk0.model
			if os.path.isdir('junk1.image'):
				os.system('rm -rf junk0.image')
				os.system('cp -r junk1.image junk0.image') # Move the previous round image to pre-previous round for num_iter>1
			else:
				os.system('cp -r '+file_str+'.image junk0.image') # Copying image to junk0.image
			if os.path.isdir('junk1.mask'):
				os.system('rm -rf junk0.mask')
				os.system('cp -r junk1.mask junk0.mask') # Move the previous round mask to pre-previous round for num_iter>1
			else:
				os.system('cp -r '+file_str+'.model junk0.mask') # Copying model to junk0.mask
			if os.path.isdir('junk1.residual'):
				os.system('rm -rf junk0.residual')
				os.system('cp -r junk1.residual junk0.residual') # Move the previous round residual to pre-previous round for num_iter>1
			else:
				os.system('cp -r '+file_str+'.residual junk0.residual') # Copying model to junk0.residual
			os.system('rm -rf junk1.ms junk1.cal junk1.model junk1.mask junk1.image junk1.residual')
			os.system('cp -r '+caltable_name+' junk1.cal') # Copying caltable to junk1.cal
			os.system('cp -r '+self.msname+' junk1.ms') # Copying ms to junk1.ms
			os.system('cp -r '+file_str+'.model junk1.model') # Copying model to junk1.model
			os.system('cp -r '+file_str+'.mask junk1.mask') # Copying mask to junk1.mask
			os.system('cp -r '+file_str+'.image junk1.image') # Copying image to junk1.image
			os.system('cp -r '+file_str+'.residual junk1.residual') # Copying residual to junk1.resuidual
			if self.verbose and (msg_code==0 or msg_code==8 or msg_code==9):	
				os.system('cp -r '+caltable_name+' '+self.mspath+'/'+file_str_prefix+'/backup_cal/'+file_str+'.cal') # Verbose=True, keep all the caltables
				os.system('cp -r '+self.msname+' '+self.mspath+'/'+file_str_prefix+'/backup_ms/'+file_str+'.ms') # If Verbose=True, keep all the ms
				os.system('cp -r '+file_str+'.model '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.model') # If Verbose=True, keep all the models
				os.system('cp -r '+file_str+'.mask '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.mask') # If Verbose=True, keep all the masks
				os.system('cp -r '+file_str+'.image '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.image') # If Verbose=True, keep all the image
				os.system('cp -r '+file_str+'.residual '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.residual') # If Verbose=True, keep all the residuals
			os.system('rm -rf '+file_str+'.image '+file_str+'.model '+file_str+'.residual '+file_str+'.sumwt '+file_str+'.pb '+file_str+'.psf '+file_str+'.mask ') 
				# Removing all imaging related files
		os.chdir(self.cwd)
		os.system('rm -rf casa*log')
		return

	def initial_bpass(self,modelname,num_iter,rms_thresh,ref_ant,minsnr,calmode,calibrator_caltable=[]):
		'''
		Function to perform and apply bandpass calibration
		Parameters:
		modelname = Name of the initial model
		num_iter = Selfcal iteration number
		rms_thresh = [], list of n-sigma value for rms based flagging (Do not go below 6)
		ref_ant = Reference antenna number
		calmode ='p' pr 'ap'
		calibrator_caltable=[], list of any previous caltables
		Return:
		Name of the initial bandpass table
		'''
		AM=am.AccessMS(self.msname)
		bandwidth=AM.calc_bandwidth()
		cent_freq=AM.calc_meanfreq()
		header=imhead(imagename=modelname,mode='list')
		caltable_name=self.msname.split('.ms')[0]+'.cal' # Caltable name
		if header['ctype1']=='Frequency':
			header_key='crval1'
			header_key_1='cdelt1'
		elif header['ctype2']=='Frequency':
			header_key='crval2'
			header_key_1='cdelt2'
		elif header['ctype3']=='Frequency':
			header_key='crval3'
			header_key_1='cdelt3'
		elif header['ctype4']=='Frequency':
			header_key='crval4'
			header_key_1='cdelt4'
		self.log_verbose.info('imhead(imagename=\''+modelname+'\',mode=\'put\',hdkey=\''+header_key+'\',hdvalue={\'value\':'+str(cent_freq)+',\'unit\':\'Hz\'},verbose=False)\n')
		imhead(imagename=modelname,mode='put',hdkey=header_key,hdvalue={'value':cent_freq,'unit':'Hz'},verbose=False)
		self.log_verbose.info('imhead(imagename=\''+modelname+'\',mode=\'put\',hdkey=\''+header_key_1+'\',hdvalue={\'value\':'+str(bandwidth)+',\'unit\':\'Hz\'},verbose=False)\n')
		imhead(imagename=modelname,mode='put',hdkey=header_key_1,hdvalue={'value':bandwidth,'unit':'Hz'},verbose=False)
		self.log_verbose.info('delmod(vis=\''+self.msname+'\',scr=True)\n')
		delmod(vis=self.msname,scr=True)
		self.log_verbose.info('ft(vis=\''+self.msname+'\',model=\''+modelname+'\',usescratch=True)\n')
		ft(vis=self.msname,model=modelname,usescratch=True)
		self.log_verbose.info('bpass_solver(\''+caltable_name+'\',spw=\'\',timerange=\'\',calmode=\''+calmode+'\',uvrange=\''+self.calib_uvrange+'\',solnorm=True,refant=\''+\
							str(ref_ant)+'\',minsnr='+str(minsnr)+',solmode=\'R\',rmsthresh=[10,8,6],gaintable='+str(calibrator_caltable)+')\n')
		caltable_name=self.bpass_solver(caltable_name,spw='',timerange='',calmode=calmode,uvrange=self.calib_uvrange,solnorm=True,refant=str(ref_ant),minsnr=minsnr,\
						solmode='R',rmsthresh=[10,8,6],gaintable=calibrator_caltable) # Perform bandpass calibration
		calibrator_caltable.append(caltable_name)
		self.log_verbose.info('applycal(vis=\''+self.msname+'\',gaintable='+str(calibrator_caltable)+',applymode=\'calflag\',flagbackup=True)\n')
		applycal(vis=self.msname,gaintable=calibrator_caltable,applymode='calflag',flagbackup=True)
		os.system('rm -rf casa*log')
		return caltable_name

	def bpass_solver(self,caltable,spw='',timerange='',calmode='ap',uvrange='',solnorm=True,refant='1',minsnr=3,solmode='',rmsthresh=[],gaintable=[]):
		'''
		Function to solve bandpass phase-only or amplitude-phase
		Parameters:
		caltable = Name of the caltable
		spw= '' , spectral window range
		timerange = '', timerange
		calmode = 'p' or 'ap'
		uvrange = '' , UV-range for calibration
		solnorm = True, normalise solution or not
		refant = Reference antenna
		minsnr = Minimum gain SNR
		solmode = 'R' for robust calibration
		rmsthresh = [], list of n-sigma values for rms based flagging (Do not put below 6)
		gaintable = [] , list of any previous caltables
		Return : Name of the caltable
		'''
		cb=calibrater()
		if solmode!='R' or len(rmsthresh)==0:
			cb.open(self.msname)
			cb.selectvis(spw=spw,uvrange=uvrange,time=timerange)
			if len(gaintable)!=0:
				for gt in gaintable:
					cb.setapply(table=gt)
			cb.setsolve(type='B',t='inf',refant=str(refant),apmode=calmode,table=caltable,append=False,minsnr=minsnr,solnorm=solnorm)
			cb.solve()
			cb.close()
		elif solmode=='R' and len(rmsthresh)!=0:
			os.system('cp -r '+self.msname+' '+self.msname.split('.ms')[0]+'_temp.ms')
			visname=self.msname.split('.ms')[0]+'_temp.ms'
			for rms in rmsthresh:
				c=0
				self.log_verbose.info('Calibrating and flagging on threshold :'+str(rms)+' sigma\n')
				while c==0:
					cb.open(visname)
					cb.selectvis(spw=spw,uvrange=uvrange,time=timerange)
					if len(gaintable)!=0:
						for gt in gaintable:
							cb.setapply(table=gt)
					cb.setsolve(type='B',t='inf',refant=str(refant),apmode=calmode,table=caltable,append=False,minsnr=minsnr,solnorm=solnorm)
					cb.solve()
					cb.close()
					applycal(vis=visname,gaintable=gaintable+[caltable],applymode='calflag',flagbackup=False)
					num_flag=fg.flagger(visname,float(rms))
					if num_flag==0:
						c=1
			os.system('rm -rf '+visname)
		os.system('rm -rf casa*log')
		return	caltable

	def cal_solar_phaseshift(self,imagename,sigma=10): # TODO :PyBDSF part is not tested
		'''
		Function to correct an average phase shift of the Sun if no mask is used. Only shifted when the shift between phasecenter and the Sun is greater than 25 arcmin.
		Parameters:
		imagename = Name of the image
		sigma = Threshold for model fitting (default =10)
		Return:
		New RA, DEC in degree, (True/False) for shift required or not
		'''
		AM=am.AccessMS(self.msname)
		radec_str,radeg,decdeg=AM.get_phasecenter()
		imsmooth(imagename=imagename,outfile='convolved.image',beam={'major':'1920arcsec','minor':'1920arcsec','pa':'0deg'},targetres=False)
		modelname=imagename.split('.image')[0]+'.model'
		if bdsf_import==True:
			exportfits(imagename='convolved.image',fitsimage='convolved.fits')
			bdsf.process_image('convolved.fits',thresh_isl=sigma,thresh_pix=sigma-2,output_opts=True,output_all=True)
			fitsfile=glob.glob('convolved_pybdsm/*/cata*/*srl.FITS')[0]
			data=fits.getdata(fitsfile)
			os.system('rm -rf convolved*')
			ra=data['RA']
			dec=data['DEC']
			major=data['Major']
			minor=data['Minor']
			if len(major)>1:
				maxpos=np.argmax(major)
				ra=ra[maxpos]
				dec=dec[maxpos]
			if np.sqrt((ra-radeg)**2+(dec-decdeg)**2)<(16/60.0):
				os.system('rm -rf casa*log')
				return radeg,decdeg,False
			else:
				os.system('rm -rf casa*log')
				return ra,dec,True
		elif bdsf_import==False: # Using CASA imfit
			rms=imstat(imagename='convolved.image',box=self.rms_box)['rms'][0]
			ia=image()			
			ia.open('convolved.image')
			ia.calcmask('convolved.image>'+str(sigma*rms),'mymask')
			ia.close()
			imfit(imagename='convolved.image',summary='convolved.image.summary')
			data=np.loadtxt('convolved.image.summary')
			if self.verbose:
				os.system('cp -r convolved* freq_*datetime_*')
			os.system('rm -rf convolved*')
			ra=data[5]
			dec=data[6]
			if np.sqrt((ra-radeg)**2+(dec-decdeg)**2)<(16/60.0):
				os.system('rm -rf casa*log')
				return radeg,decdeg,False
			else:
				os.system('rm -rf casa*log')
				return ra,dec,True

	def shift_phasecenter(self,imagename,ra,dec):
		'''
		Function to shift image reference RA DEC to a certain value
		Parameters:
		imagename = Name of the image
		ra = New RA in degree
		dec = New DEC in degree
		'''
		cwd=os.getcwd()
		AM=am.AccessMS(self.msname)
		radec_str,radeg,decdeg=AM.get_phasecenter()
		IB=B.ImageBasic(self.msname)
		psf=IB.calc_psf()/3600.0
		if np.sqrt((ra-radeg)**2+(dec-decdeg)**2)>psf:
			try:
				exportfits(imagename=imagename,fitsimage='wcs_model.fits',dropdeg=True,dropstokes=True)
				w=WCS('wcs_model.fits')
				pix=np.mean(w.all_world2pix(np.array([[ra,dec],[ra,dec]]),0),axis=0)
				ra_pix=int(pix[0])
				dec_pix=int(pix[1])
				os.system('rm -rf wcs_model.fits')
				try:
					imsubimage(imagename=imagename,outfile='shift_model.model',stokes='I',dropdeg=False)
					exportfits(imagename='shift_model.model',fitsimage='shift_model.fits',dropstokes=False,dropdeg=False)												 		
					hdul=fits.open('shift_model.fits')
					hdr=hdul[0].header
					data=(hdul[0].data)
					hdr['CRPIX1']=ra_pix
					hdr['CRPIX2']=dec_pix
					fits.writeto('shift_model.fits',data=data,header=hdr,overwrite=True)
					os.system('rm -rf '+imagename+' shift_model.model')
					importfits(fitsimage='shift_model.fits',imagename=imagename)
					self.log_verbose.info('Image phase center shifted to , RA : '+str(radec_str[0])+', DEC : '+str(radec_str[1])+'\n')
					os.system('rm -rf shift_model.fits')
					os.system('rm -rf casa*log')
					os.chdir(cwd)
					return 0
				except:
					self.log_verbose.info('Image phase center could not be shifted. Please provide only Stokes I image.\n')
			except:
				try:
					importfits(fitsimage=imagename,imagename=imagename+'.model')
					try:
						imsubimage(imagename=imagename+'.model',outfile='I.model',stokes='I',dropdeg=False)
						os.system('rm -rf '+imagename+'.model')
						exportfits(imagename='I.model',fitsimage='wcs_model.fits',dropdeg=True,dropstokes=True)
						w=WCS('wcs_model.fits')
						pix=np.mean(w.all_world2pix(np.array([[ra,dec],[ra,dec]]),0),axis=0)
						ra_pix=int(pix[0])
						dec_pix=int(pix[1])
						os.system('rm -rf wcs_model.fits')
						exportfits(imagename='I.model',fitsimage='shift_model.fits',dropstokes=False,dropdeg=False)												 		
						hdul=fits.open('shift_model.fits')
						hdr=hdul[0].header
						data=(hdul[0].data)
						hdr['CRPIX1']=ra_pix
						hdr['CRPIX2']=dec_pix
						fits.writeto('shift_model.fits',data=data,header=hdr,overwrite=True)
						os.system('rm -rf I.model shift_model.model '+imagename)
						os.system('mv shift_model.fits '+imagename)
						self.log_verbose.info('Image phase center shifted to , RA : '+str(radec_str[0])+', DEC : '+str(radec_str[1])+'\n')
						os.system('rm -rf shift_model.fits')
						os.system('rm -rf casa*log')
						os.chdir(cwd)
						return 0
					except:
						self.log_verbose.info('Image phase center could not be shifted. Please provide only Stokes I image.\n')	 
				except:
					self.log_verbose.info('Image is not either in CASA or fits format.\n')
					os.system('rm -rf casa*log')
					os.chdir(cwd)
					return 2
		else:
			os.system('rm -rf casa*log')
			os.chdir(cwd)
			return 1
		
	def image_source_true_loc(self,outdir,do_bandpass=False):
		'''
		Function to make quick image for finding source true location with respect to reference time and channel
		Parameters:
		Outdir = Output directory
		do_bandpass = False, Source true location for bandpass self-calibration or not 
		Return:
		Output imagename
		'''
		os.system('rm -rf test_loc*')
		poltclean(vis=self.msname,imagename='test_loc',selectdata=True,startmodel='',stokes='I',antenna='',imsize=[self.imsize],cell=self.cellsize,\
				niter=0,gain=0.5,threshold=['10Jy'],deconvolver='multiscale',scales=self.multiscale_scales,uvtaper=self.uvtaper,weighting='natural',\
				interactive=False,mask='')
		maxpix=imstat(imagename='test_loc.image')['max'][0]
		thresh=str(maxpix/5)+'Jy'
		os.system('rm -rf test_loc*')
		if do_bandpass==True:
			poltclean(vis=self.msname,imagename='test_loc',selectdata=True,startmodel='',stokes='I',antenna='',imsize=[self.imsize],cell=self.cellsize,\
					niter=100000,gain=0.5,threshold=[thresh],specmode='cube',deconvolver='multiscale',scales=self.multiscale_scales,uvtaper=self.uvtaper,weighting='natural',\
					interactive=False,mask='')
			self.log_verbose.info('poltclean(vis=\''+self.msname+'\',imagename=\'test_loc\',selectdata=True,startmodel=\'\',stokes=\'I\',antenna=\'\',imsize=['\
					+str(self.imsize)+'],cell=\''+str(self.cellsize)+'\',niter=100000,gain=0.5,threshold=[\''+str(thresh)+'\'],specmode=\'cube\',deconvolver=\'multiscale\',scales='+\
					str(self.multiscale_scales)+',uvtaper=\''+str(self.uvtaper)+'\',weighting=\'natural\',interactive=False,mask=\'\')\n')
		else:
			poltclean(vis=self.msname,imagename='test_loc',selectdata=True,startmodel='',stokes='I',antenna='',imsize=[self.imsize],cell=self.cellsize,\
					niter=100000,gain=0.5,threshold=[thresh],deconvolver='multiscale',scales=self.multiscale_scales,uvtaper=self.uvtaper,weighting='natural',\
					interactive=False,mask='')
			self.log_verbose.info('poltclean(vis=\''+self.msname+'\',imagename=\'test_loc\',selectdata=True,startmodel=\'\',stokes=\'I\',antenna=\'\',imsize=['\
					+str(self.imsize)+'],cell=\''+str(self.cellsize)+'\',niter=100000,gain=0.5,threshold=[\''+str(thresh)+'\'],deconvolver=\'multiscale\',scales='+\
					str(self.multiscale_scales)+',uvtaper=\''+str(self.uvtaper)+'\',weighting=\'natural\',interactive=False,mask=\'\')\n')
		if os.path.isdir(outdir)==False:
			os.makedirs(outdir)
		if outdir[-1]=='/':
			outdir=outdir[:-1]
		imagename=outdir+'/'+os.path.basename(self.msname).split('.ms')[0]+'_true_loc.fits' # Imagename prefix
		if os.path.isfile(imagename):
			os.system('rm -rf '+imagename)
		exportfits(imagename='test_loc.image',fitsimage=imagename)
		os.system('rm -rf test_loc.* casa*log')
		return imagename

	def selfcal_iteration(self,num_iter,rms_thresh,sigma,maskstr,antenna_to_use,startmodel,startmask,ref_ant,minsnr,calmode,maskfile='',want_auto_masking=False,\
							stokes='I',interactive=False,do_bandpass=False,solmode='R',correct_phasecenter=False,ra=0,dec=0,box_width=3,calibrator_caltable=[]):
		'''
		Function to perform a self-calibration loop, make an image, put the model in the measurement set, and perform the calibration
		Parameters:
		num_iter = Number of self-calibration iteration
		rms_thresh = RMS for threshold, list
		sigma = Threshold sigma
		maskstr = Mask string for CLEANing
		antenna_to_use = List of antennas for CLEANing
		startmodel = Model to start the CLEANing
		startmask = Mask to start 
		ref_ant = Reference antenna
		minsnr = Minimum gain SNR
		calmode= = 'p' for phase-only and 'ap' for amplitude-phase calibration
		maskfile = Name of the maskfile
		want_auto_masking = False, if True use CASA auto-multithresh for auto masking
		stokes = 'I', Stokes plane to image
		interactive= False, Perform interactive selfcal, change options on the fly
		do_bandpass = False, Perform bandpass calibration
		solmode ='R', 'L1R', 'L1' for gaincal and only 'R' for bandpass
		correct_phasecenter = False, correct the phase center
		ra = New RA to change phasecenter
		dec = New DEC to change phasecenter
		box_width = Negative box width around the maximum pixel in degree (default : 3 degree)
		calibrator_caltable = List of calilbrator caltables
		Return:
		rms based dynamic range,negative based dynamic range, error message code
		'''
		os.chdir(self.mspath)
		imagename=self.msname.split('.ms')[0]+'_'+str(num_iter) # Imagename prefix
		present_file=glob.glob(imagename+'*')
		if len(present_file)!=0:
			os.system('rm -rf '+imagename+'*')
			self.log_verbose.info('rm -rf '+imagename+'*\n')
		caltable_name=self.msname.split('.ms')[0]+'.cal' # Caltable name
		if os.path.isdir(caltable_name):
			os.system('rm -rf '+caltable_name)
		threshold=[str(rms*sigma)+'Jy' for rms in rms_thresh]
		# Making image
		self.log_verbose.info('============================\n')
		self.log_verbose.info('Iteration number : '+str(num_iter)+'\n')
		self.log_verbose.info('============================\n')
		if calmode=='p':
			clean_gain=0.05
		else:
			clean_gain=0.1
		if maskfile=='':
			maskfile=maskstr
		if maskfile!='': # Using the user provided mask file
			self.log_verbose.info('poltclean(vis=\''+self.msname+'\',imagename=\''+imagename+'\',selectdata=True,startmodel=\''+startmodel+'\',startmask=\''+\
					startmask+'\',stokes=\''+stokes+'\',antenna=\''+antenna_to_use+'\',imsize=['+str(self.imsize)+'],cell=\''+str(self.cellsize)\
					+'\',niter=100000000000,gain='+str(clean_gain)+',threshold='+str(threshold)+',deconvolver=\'multiscale\',scales='+str(self.multiscale_scales)
					+',uvtaper=\''+self.uvtaper+'\',weighting=\'natural\',interactive=False,mask=\''+str(maskfile)+'\')\n')
			poltclean(vis=self.msname,imagename=imagename,selectdata=True,startmodel=startmodel,startmask=startmask,stokes=stokes,antenna=antenna_to_use,imsize=[self.imsize],\
				cell=self.cellsize,niter=100000000000,gain=clean_gain,threshold=threshold,deconvolver='multiscale',scales=self.multiscale_scales,uvtaper=self.uvtaper,\
				weighting='natural',interactive=False,mask=maskfile)
		elif want_auto_masking==True and maskfile=='' and maskstr=='': # Use auto-masking
			try_count=0
			while True:
				if try_count==0:
					self.log_verbose.info('Normal auto-masking.\n')
					self.log_verbose.info('poltclean(vis=\''+self.msname+'\',imagename=\''+imagename+'\',selectdata=True,startmodel=\''+startmodel+'\',startmask=\''\
						+startmask+'\',stokes=\''+stokes+'\',antenna=\''+antenna_to_use+'\',imsize=['+str(self.imsize)+'],cell=\''+str(self.cellsize)+\
						'\',niter=100000000000,gain='+str(clean_gain)+',threshold='+str(threshold)+',deconvolver=\'multiscale\',scales='+str(self.multiscale_scales)+\
						',uvtaper=\''+str(self.uvtaper)+'\',weighting=\'natural\',interactive=False,usemask='+\
						'\'auto-multithresh\',mask=\'\',pbmask=0.0,sidelobethreshold=1.5,noisethreshold=3.0,lownoisethreshold=1.5,negativethreshold=0.0,smoothfactor=1.0,'+\
						'minbeamfrac=0.1,growiterations=75,minpercentchange=5.0)\n')
					poltclean(vis=self.msname,imagename=imagename,selectdata=True,startmodel=startmodel,startmask=startmask,stokes=stokes,antenna=antenna_to_use,imsize=[self.imsize],\
					cell=self.cellsize,niter=100000000000,gain=clean_gain,threshold=threshold,deconvolver='multiscale',scales=self.multiscale_scales,uvtaper=self.uvtaper,\
					weighting='natural',interactive=False,usemask='auto-multithresh',mask='',pbmask=0.0,sidelobethreshold=1.5,noisethreshold=3.0,lownoisethreshold=1.5,\
					negativethreshold=0.0,smoothfactor=1.0,minbeamfrac=0.5,growiterations=75,minpercentchange=5.0)
				elif try_count==1:
					self.log_verbose.info('Trying with auto-masking with no restriction of minimum beam fraction.\n')
					self.log_verbose.info('poltclean(vis=\''+self.msname+'\',imagename=\''+imagename+'\',selectdata=True,startmodel=\''+startmodel\
						+'\',startmask=\''+startmask+'\',stokes=\''+stokes+',antenna=\''+antenna_to_use+'\',imsize=['+str(self.imsize)+'],cell=\''+str(self.cellsize)+\
						'\',niter=100000000000,gain='+str(clean_gain)+',threshold='+str(threshold)+',deconvolver=\'multiscale\',scales='+str(self.multiscale_scales)+\
						',uvtaper=\''+str(self.uvtaper)+'\',weighting=\'natural\',interactive=False,usemask=\''+\
						'auto-multithresh\',mask=\'\',pbmask=0.0,sidelobethreshold=1.5,noisethreshold=3.0,lownoisethreshold=1.5,negativethreshold=0.0,smoothfactor=1.0,'+\
						'minbeamfrac=0.1,growiterations=75,minpercentchange=-1.0)\n')
					poltclean(vis=self.msname,imagename=imagename,selectdata=True,startmodel=startmodel,startmask=startmask,stokes=stokes,antenna=antenna_to_use,imsize=[self.imsize],\
					cell=self.cellsize,niter=100000000000,gain=clean_gain,threshold=threshold,deconvolver='multiscale',scales=self.multiscale_scales,uvtaper=self.uvtaper,\
					weighting='natural',interactive=False,usemask='auto-multithresh',mask='',pbmask=0.0,sidelobethreshold=1.5,noisethreshold=3.0,lownoisethreshold=1.5,\
					negativethreshold=0.0,smoothfactor=1.0,minbeamfrac=0.1,growiterations=75,minpercentchange=-1.0)
				elif try_count==2:
					self.log_verbose.info('Trying without masking.\n')
					self.log_verbose.info('poltclean(vis=\''+self.msname+'\',imagename=\''+imagename+'\',selectdata=True,startmodel=\''+startmodel+'\',startmask=\''+\
						startmask+'\'stokes=\''+stokes+'\',antenna=\''+antenna_to_use+'\',imsize=['+str(self.imsize)+'],cell=\''+str(self.cellsize)+\
						'\',niter=100000000000,gain='+str(clean_gain)+',threshold='+str(threshold)+',deconvolver=\'multiscale\',scales='+str(self.multiscale_scales)+\
						',uvtaper=\''+str(self.uvtaper)+'\',weighting=\'natural\',interactive=False,usemask=\'user\',mask=\'\')\n')
					poltclean(vis=self.msname,imagename=imagename,selectdata=True,startmodel=startmodel,startmask=startmask,stokes=stokes,antenna=antenna_to_use,imsize=[self.imsize],\
					cell=self.cellsize,niter=100000000000,gain=clean_gain,threshold=threshold,deconvolver='multiscale',scales=self.multiscale_scales,uvtaper=self.uvtaper,\
					weighting='natural',interactive=False,usemask='user',mask='')
				else:
					break
				modelflux=imstat(imagename=imagename+'.model')['sum'][0]
				if modelflux==0.0:
					try_count+=1
				else:
					break
		else: # If no masking
			self.log_verbose.info('poltclean(vis=\''+self.msname+'\',imagename=\''+imagename+'\',selectdata=True,startmodel=\''+startmodel+'\',dtartmask=\''+startmask\
					+'\',stokes=\''+stokes+'\',antenna=\''+antenna_to_use+'\',imsize=['+str(self.imsize)+'],cell=\''+str(self.cellsize)\
					+'\',niter=100000000000,gain='+str(clean_gain)+',threshold=\''+str(threshold)+'\',deconvolver=\'multiscale\',scales='+str(self.multiscale_scales)
					+',uvtaper=\''+self.uvtaper+'\',weighting=\'natural\',interactive=False,mask='+str([maskfile])+')\n')		
			poltclean(vis=self.msname,imagename=imagename,selectdata=True,startmodel=startmodel,startmask=startmask,stokes=stokes,antenna=antenna_to_use,imsize=[self.imsize],\
				cell=self.cellsize,niter=100000000000,gain=clean_gain,threshold=threshold,deconvolver='multiscale',scales=self.multiscale_scales,uvtaper=self.uvtaper,\
				weighting='natural',interactive=False,mask='')
		if stokes=='I':
			stokes_list=['I']
		elif stokes=='RR':
			stokes_list=['RR']
		elif stokes=='LL':
			stokes_list=['LL']
		elif stokes=='XX':
			stokes_list=['XX']
		elif stokes=='YY':
			stokes_list=['YY']
		elif stokes=='RRLL':
			stokes_list=['RR','LL']
		elif stokes=='XXYY':
			stokes_list=['XX','YY']
		out_dict,negative_dyn_range=self.calc_dyn_range(num_iter,sigma,box_width=box_width,stokes_list=stokes_list) # Calculating the dynamic range of the image
		out_dict_keys=out_dict.keys()
		if 'NAN' in out_dict_keys:
			self.log_verbose.error(B.error_msgs(3))
			os.system('rm -rf casa*log')
			os.chdir(self.cwd)
			os.system('rm -rf casa*log')
			return 3     # If image is not made, no point in continuing
		if os.path.isdir(imagename+'.model')==False:
			self.log_verbose.error(B.error_msgs(4)+'\n')
			os.system('rm -rf casa*log')
			os.chdir(self.cwd)
			os.system('rm -rf casa*log')
			return 4	   # If model is not present no point in continuing
		else:
			modelflux=imstat(imagename=imagename+'.model')['sum'][0]
			if modelflux==0.0:
				self.log_verbose.error(B.error_msgs(5)+'\n')
				os.system('rm -rf casa*log')
				os.chdir(self.cwd)
				os.system('rm -rf casa*log')	
				return 5 # If modelflux is 0, no point in continuing
			else:
				flaglist=flagmanager(vis=self.msname,mode='list')
				if len(flaglist)>1:
					flaglist_keys=list(flaglist.keys())
					flaglist_keys.remove('MS')
					last_flag_key=len(flaglist_keys)-1
					last_flagversion=flaglist[last_flag_key]['name']
					# Restore the flag and delete the present flag version
					if 'applycal' in last_flagversion:  # Restore flagversion if the last flagversion is from applycal
						self.log_verbose.info('flagmanager(vis=\''+self.msname+'\',mode=\'restore\',versionname=\''+str(last_flagversion)+'\',merge=\'replace\')\n')
						self.log_verbose.info('flagmanager(vis=\''+self.msname+'\',mode=\'delete\',versionname=\''+str(last_flagversion)+'\')\n')
						flagmanager(vis=self.msname,mode='restore',versionname=last_flagversion,merge='replace')
						flagmanager(vis=self.msname,mode='delete',versionname=last_flagversion)
				self.log_verbose.info('clearcal(vis=\''+self.msname+'\')\n')
				clearcal(vis=self.msname)
				self.log_verbose.info('delmod(vis=\''+self.msname+'\',scr=True)\n') 
				delmod(vis=self.msname,scr=True) # Clear the MODEL column
				if correct_phasecenter==True:
					if ra==0 or dec==0:
						AM=am.AccessMS(self.msname)
						radec_str,ra,dec=AM.get_phasecenter()
					shifted=self.shift_phasecenter(imagename+'.model',ra,dec)
					if shifted==0:
						self.log_verbose.info('Image phasecenter changed to : RA = '+str(ra)+' deg, DEC = '+str(dec)+' deg\n')
					elif shifted==1:
						self.log_verbose.info('Shifting is not required. Continuing with previous phasecenter.\n')
					else:
						self.log_verbose.info('Error in image phsecenter shifting. Continuing with previous phasecenter.\n')
				self.log_verbose.info('ft(vis=\''+self.msname+'\',model=\''+imagename+'.model\',nterms=1,usescratch=True)\n') 
				ft(vis=self.msname,model=imagename+'.model',nterms=1,usescratch=True) # Putting the model into MS
				if do_bandpass==True:
					if solmode=='L1' or solmode=='L1R':
						solmode=''
					self.log_verbose.info('Doing bandpass.....\n')
					self.log_verbose.info('bpass_solver(\''+caltable_name+'\',spw=\'\',timerange=\'\',calmode=\''+calmode+'\',uvrange=\''+self.calib_uvrange+\
				'\',solnorm=True,refant=\''+str(ref_ant)+'\',minsnr='+str(minsnr)+',solmode=\''+solmode+'\',rmsthresh=[15,10,8],gaintable='+str(calibrator_caltable)+')\n')
					self.bpass_solver(caltable_name,spw='',timerange='',calmode=calmode,uvrange=self.calib_uvrange,solnorm=True,refant=str(ref_ant),minsnr=minsnr,\
									solmode=solmode,rmsthresh=[15,10,8],gaintable=calibrator_caltable) # Perform bandpass calibration
				else:
					self.log_verbose.info('Doing gaincal.....\n')
					self.log_verbose.info('gaincal(vis=\''+self.msname+'\',caltable=\''+caltable_name+'\',refant=\''+str(ref_ant)+'\',minsnr='+str(minsnr)
							+',calmode=\''+calmode+'\',solnorm=True,uvrange=\''+self.calib_uvrange+'\',gaintype=\'G\',solmode=\''+solmode+'\',rmsthresh=[10,7,5,3.5],gaintable='+
							str(calibrator_caltable)+')\n') 
					gaincal(vis=self.msname,caltable=caltable_name,refant=str(ref_ant),minsnr=minsnr,calmode=calmode,solnorm=True,uvrange=self.calib_uvrange,\
						gaintype='G',solmode=solmode,rmsthresh=[10,7,5,3.5],gaintable=calibrator_caltable) # Performing gain calibration
				if os.path.isdir(caltable_name)==False:
					self.log_verbose.info('No good solution found. No caltable made.\n')
					os.system('rm -rf casa*log')
					os.chdir(self.cwd)
					os.system('rm -rf casa*log')
					return 7
				applycal_caltable=copy.deepcopy(calibrator_caltable)
				applycal_caltable.append(caltable_name)
				self.log_verbose.info('Applying solutions from :'+str(applycal_caltable)+'\n')
				self.log_verbose.info('applycal(vis=\''+self.msname+'\',gaintable='+str(applycal_caltable)+',applymode=\'calflag\',flagbackup=True,calwt=[False])\n')	
				applycal(vis=self.msname,gaintable=applycal_caltable,applymode='calflag',flagbackup=True,calwt=[False]) # Applying the solution
				self.log_verbose.info('Success.\n')
				os.system('rm -rf casa*log')
				os.chdir(self.cwd)
				os.system('rm -rf casa*log')
				return 0,out_dict,negative_dyn_range
			
#########################################
# Finished Intensity Selfcal Class
#########################################


