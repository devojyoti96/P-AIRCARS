import os
import numpy as np,copy,sys,glob
import logging,paircars
from casatools import msmetadata,table,measures,quanta,agentflagger,image,calibrater,ms
from casatasks import *
from . import access_ms as am
from . import basic_func as B
from . import flagger as fg
from CALIBRATE.access_calibrate import *
from numpy.linalg import inv,det
from scipy.linalg import polar
from mwa_pb.mwapb import *
from paircars_casatasks.poltclean import *
import matplotlib,matplotlib.pyplot as plt
import scipy.linalg
from mpl_toolkits.mplot3d import Axes3D
matplotlib.use('Agg')
os.system('rm -rf casa*log')
'''
Code is written by Devojyoti Kansabanik, 01 Mar, 2021
'''

datadir = os.path.dirname(__file__)
class PolSelfcal:
	'''
	Generic class to perform polarisation self-calibration (Using Andre Offringa's CALIBRATE code on based Mitchcal algorithm)

	Parameters
	----------
	msname : str
		Name of the measurement set
	metafits : str
		Name of the MWA metafits file
	num_pixel_in_psf : int
		Number of pixels side one point-spread-funtion
	maximum_emission_scale : float 
		Maximum scale of the emission present in the image
	largest_scale : float 
		Largest spatial scale in degree used for self calibration
	verbose : bool
 		If True keep all the intermediate images, model, residuals, caltables and details of the log to detailed analysis
	interactive : bool 
		If True user have interactive control on self-calibration
	savelog : bool
		Save log
	use_wsclean : bool
		Use WSClean or not
	'''
	
	def __init__(self,msname,metafits,maximum_emission_scale,num_pixel_in_psf=5,largest_scale=12,verbose=False,interactive=False,savelog=True,use_wsclean=True):
		self.cwd=os.getcwd()
		if msname[-1]=='/':
			self.msname=msname[:-1]
		else:
			self.msname=msname
		self.mspath=os.path.dirname(os.path.realpath(msname))
		AM=am.AccessMS(self.msname)
		IB=B.ImageBasic(self.msname)
		self.metafits=metafits
		self.max_baseline=AM.get_max_baseline()
		self.num_pixel_in_psf=int(num_pixel_in_psf)
		self.cellsize=IB.calc_cellsize(self.num_pixel_in_psf) 
		self.imsize=IB.num_pixels(self.num_pixel_in_psf)
		self.max_size=maximum_emission_scale
		self.multiscale_scales=IB.choose_scales(self.num_pixel_in_psf,self.max_size)
		self.uvtaper=IB.calc_uvtaper()
		calib_uvrange_min=IB.calc_calib_uvrange(largest_scale)
		self.calib_uvrange_min=calib_uvrange_min[1]
		self.calib_uvrange_max=calib_uvrange_min[2]
		self.imaging_uvrange=calib_uvrange_min[0]
		self.imaging_minuv=calib_uvrange_min[3]
		self.imaging_maxuv=calib_uvrange_min[4]
		self.rms_box='50,50,'+str(self.imsize-50)+','+str(int(self.imsize/4)) # CASA box to calculate the rms
		self.verbose=verbose
		self.interactive=interactive
		self.wsclean=use_wsclean
		self.savelog=savelog
		formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
		self.pollog_verbose = logging.getLogger('polselfcal_verbose_log')
		self.pollog_verbose.setLevel(logging.DEBUG)
		if self.verbose:
			self.console=logging.StreamHandler(sys.stdout)
			self.console.setFormatter(formatter)
			self.pollog_verbose.addHandler(self.console)
		if self.savelog==True:
			self.filehandle=logging.FileHandler(self.cwd+'/Pol_Selfcal_verbose.log')
			self.filehandle.setFormatter(formatter)
			self.pollog_verbose.addHandler(self.filehandle)
		self.pollog_verbose.propagate = False
		self.pollog_verbose.info('Initiating Polarisation selfcal object.\n')
		if self.wsclean==True:
			datadir=os.path.abspath(os.path.dirname(paircars.__file__))
			try:
				wsclean_path=str(np.load(datadir+'/wsclean_path.npy',allow_pickle=True))
				os.path.join(wsclean_path)
				a=os.system('wsclean > wsclean_test')
				if a!=0:
					self.log_verbose.info('WSClean is not installed. Using CASA for imaging.\n')
					self.wsclean=False
			except:
				self.log_verbose.info('WSClean is not installed. Using CASA for imaging.\n')
				self.wsclean=False
			os.system('rm -rf wsclean_test')
		
	def negative_box(self,max_pix,box_width=3):
		'''
		Create a 3 degree box about the maximum pixel of image to search negative

		Parameters
		----------
		max_pix : list 
			Maximum pixel [xxmax,yymax]
		box_width : int 
			Box width in degree (default : 3 degree)
		Returns
		---------
		str
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

	def calc_dyn_range(self,imagename,sigma,box_width=3,stokes_list=['I']):
		'''
		Calculate the dynamic range of the full Stokes image cube

		Parameters
		-----------
		imagename : str 
			Name of the CASA image
		sigma : float 
			nsigma value to put a mask for calculating total flux
		box_width : float 
			Negative box width around the maximum pixel in degree (default : 3 degree)
		stokes_list : list 
			List of stokes planes in the image
		Returns
		-----------	
		dict
			{'STOKES':[rms dynamic range,rms,total_flux(non-negative)]}
		float
			negative dynamic range for Stokes I
		'''
		ia=image()
		try:
			if os.path.isdir(imagename):
				imageheader=imhead(imagename=imagename,mode='summary')
				if imageheader['ndim']>2 and imageheader['ndim']==4:
					out_dict={}
					neg_dyn=0
					for stokes in stokes_list:
						if stokes=='I' or stokes=='XX' or stokes =='YY':
							if os.path.isdir('I.image')==True:
								os.system('rm -rf I.image')
							immath(imagename=imagename,outfile='I.image',mode='evalexpr',stokes=stokes)
							maxpos=imstat(imagename='I.image',stokes=stokes)['maxpos']
							negative_box=self.negative_box(maxpos,box_width=box_width)
							max_pix=imstat(imagename='I.image',stokes=stokes)['max'][0]
							rms=imstat(imagename='I.image',box=self.rms_box,stokes=stokes)['rms'][0]
							min_pix=imstat(imagename='I.image',box=negative_box,stokes=stokes)['min'][0]
							rms_dyn_range=max_pix/rms
							if min_pix!=0:
								neg_dyn+=max_pix/abs(min_pix)
							else:
								neg_dyn=rms_dyn_range
							ia.open('I.image')
							ia.calcmask('\"I.image\">'+str(sigma*rms),'mymask')
							ia.close()
							try:
								total_flux=imstat(imagename='I.image',stokes=stokes)['flux'][0]
								out_dict[stokes]=[rms_dyn_range,rms,total_flux]							
							except:
								out_dict[stokes]=[rms_dyn_range,rms,np.nan]
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
							makemask(inpimage='I.image',inpmask='I.image:mymask',output=stokes+'.image:mymask',mode='copy')
							try:
								total_flux=imstat(imagename=stokes+'.image',stokes=stokes)['flux'][0]
								out_dict[stokes]=[rms_dyn_range,rms,total_flux]
							except:
								out_dict[stokes]=[rms_dyn_range,rms,np.nan]	
							os.system('rm -rf '+stokes+'.image')		
					os.system('rm -rf I.image')		
			else:
				out_dict={}
				rms_dyn_range=np.nan
				neg_dyn=np.nan
				out_dict['NAN']=[np.nan,np.nan]
			negative_dyn_range=neg_dyn/len(stokes_list)
		except:
			negative_dyn_range=np.nan
			out_dict['NAN']=[np.nan,np.nan]
		os.system('rm -rf casa*log')
		return out_dict,negative_dyn_range

	def calc_iter_num(self,safety_factor,quality_factor,scratch=True):
		'''
		Function to calculate minimum number of selfcal iteration based on safety standard and quality factor

		Parameters
		----------
		safety_factor : int 
			Factor to determine the robustness of the selfcal
		quality_factor : int 
			Factor to determine the quality of the images
		scratch : bool
			Whether start the selfcal from scratch or not
		Returns
		----------
		int
			Minimum iteration at fixed sigma
		int
			Minimum iteration
		int
			Maximum iteration
		int
			Number of antenna bins
		float
			Fraction change in flux for convergence
		float
			Minimum value of allowed sigma
		'''
		if quality_factor==0:     # Low quality (Quick look image making)
			frac_flux_change=0.03
			pol_frac_change=0.1
			if (safety_factor==0):
				min_sigma=9.0
				min_num_iter_fixed_sigma=1
				if (scratch==True):
					min_iteration=2
					max_iteration=20
				else:
					min_iteration=1
					max_iteration=10
			elif (safety_factor==1):
				min_sigma=8.0
				min_num_iter_fixed_sigma=1
				if (scratch==True):
					min_iteration=4
					max_iteration=30
				else:
					min_iteration=1
					max_iteration=20
			else:
				min_sigma=7.0
				min_num_iter_fixed_sigma=1
				if (scratch==True):
					min_iteration=6
					max_iteration=40
				else:
					min_iteration=1
					max_iteration=30
		elif quality_factor==1:  # Medium quality imaging (Computing speed medium)
			frac_flux_change=0.015
			pol_frac_change=0.08
			if (safety_factor==0):
				min_sigma=8.0
				min_num_iter_fixed_sigma=2
				if (scratch==True):
					min_iteration=3
					max_iteration=40
				else:
					min_iteration=2
					max_iteration=30
			elif (safety_factor==1):
				min_sigma=7.0
				min_num_iter_fixed_sigma=2
				if (scratch==True):
					min_iteration=5
					max_iteration=50
				else:
					min_iteration=2
					max_iteration=40
			else:
				min_sigma=6.0
				min_num_iter_fixed_sigma=2
				if (scratch==True):
					min_iteration=7
					max_iteration=60
				else:
					min_iteration=2
					max_iteration=50
		else:  # Best quality imaging (Computing slow)
			frac_flux_change=0.01
			pol_frac_change=0.05
			if (safety_factor==0):
				min_sigma=7.0
				max_iteration=60
				min_num_iter_fixed_sigma=3
				if (scratch==True):
					min_iteration=4
				else:
					min_iteration=3
			elif (safety_factor==1):
				min_sigma=6.0
				min_num_iter_fixed_sigma=3
				if (scratch==True):
					min_iteration=6
				else:
					min_iteration=3
			else:
				min_sigma=5.0
				min_num_iter_fixed_sigma=3
				if (scratch==True):
					min_iteration=8
				else:
					min_iteration=3
		antenna_bin=1
		self.pollog_verbose.info('Quality factor :'+str(quality_factor)+', Safety standard :'+str(safety_factor)+\
				', Minimum number of iteration at fixed sigma : '+str(min_num_iter_fixed_sigma)+', Minimum iteration :'+str(min_iteration)+', Antenna bins :'+str(antenna_bin)+\
				', Fraction flux change for convergence : '+str(frac_flux_change)+'\n')
		os.system('rm -rf casa*log')
		return min_num_iter_fixed_sigma,min_iteration,max_iteration,antenna_bin,frac_flux_change,pol_frac_change,min_sigma

	def antenna_string(self,antenna_list,antenna_list_index):
		'''
		Function to return antenna string from antenna list

		Parameters
		-----------
		antenna_list : list
			Antenna list or array
		antenna_list_index : int 
			Bin number of antenna list
		Returns
		-----------
		str
			Antenna string
		'''
		antenna_string=''
		for ant in antenna_list[antenna_list_index]:
			antenna_string+=str(ant)+','
		antenna_string=antenna_string[:-1]
		os.system('rm -rf casa*log')
		return antenna_string

	def cal_poldistortion(self,gaintable,poldistortion_matrix='UH'): # TODO : Diagnostic plots
		'''
		Function to calculate the estmated Jones matrices for poldistortion after correcting for instrumental Jones matrix
		Note : Saved X matrix is for B'=XBX^\dagger, which is inverse of poldistortion_matrix of correct_poldistortion function

		Parameters
		-----------
		gaintable : str 
			Name of the gaintable (Assuming CALIBRATE gaintable format only right now)
		poldistortion _matrix : str
			'UH' or 'HU', where H is polconversion matrix and U is the polrotation matrix
		Returns
		------------
		numpy.matrix
			Poldistortion matrix
		numpy.matrix 
			Inverse of poldistortion matrix
		numpy.matrix
			Polconversion
		numpy.matrix
			Inversion of polconversion
		numpy.matrix 
			Polrotation
		numpy.matrix
			Inverse of polrotation
		str
			Filename to save poldistortion matrix 
		'''
		cal=CALIBRATE()
		if gaintable[-1]=='/':
			gaintable=gaintable[:-1]
		bin_caltable=cal.modify_caltable_for_ms(self.msname,gaintable,gaintable+'.calibrate_bin')
		npytable=cal.convert_gaintable_bin2npy(bin_caltable,gaintable+'.calibrate_bin.temp')
		jones_array=np.load(npytable,allow_pickle=True)[1]
		jones_array=jones_array.reshape(2,2,-1)
		jones_array_copy=copy.deepcopy(jones_array)
		nanpos=np.where(np.isnan(jones_array[0,0,:])==True)
		jones_array=np.delete(jones_array,nanpos,axis=2)
		jones_array=[inv(np.matrix(jones_array[:,:,i])) for i in range(jones_array.shape[-1])] # CAIBRATE caltables save inverses of Jones matrices
		d1=[np.matmul(x.H,x) for x in jones_array]
		d2=[x.H for x in jones_array]
		d1_sum=np.matrix(np.sum(d1,axis=0))
		d2_sum=np.matrix(np.sum(d2,axis=0))
		x_inv=np.matmul(inv(d1_sum),d2_sum)
		x=inv(x_inv)
		if poldistortion_matrix=='UH':
			U,H=polar(x,side='right')
		else:
			H,U=polar(x,side='left')
		os.system('rm -rf '+gaintable+'.calibrate_bin*')
		os.system('rm -rf casa*log')
		np.save(gaintable+'.poldist',np.array([inv(x)]))  
						# Saving X matrix for B'=XBX^\dagger, which is inverse of poldistortion_matrix of correct_poldistortion function
		return x,inv(x),np.matrix(H),np.matrix(inv(H)),np.matrix(U),np.matrix(inv(U)),gaintable+'.poldist'

	def correct_poldistortion(self,gaintable,outfile,poldistortion_matrix):
		'''
		Function to applycal poldistortion correction (either polconversion or polrotation or both) to the gaintable

		Parameters
		----------
		gaintable : str 
			Name of the gaintable
		outfile : str 
			Name of the output poldistortion corrected gaintable		
		poldistortion_matrix : numpy.array 
			Poldistortion matrix
		Returns
		-------
		str
			Poldistortion corrected gaintable name
		'''
		cal=CALIBRATE()
		if gaintable[-1]=='/':
			gaintable=gaintable[:-1]
		outfile_path=os.path.dirname(outfile)
		bin_caltable=cal.modify_caltable_for_ms(self.msname,gaintable,gaintable+'.calibrate_bin')
		npytable=cal.convert_gaintable_bin2npy(bin_caltable,gaintable+'.calibrate_bin.temp')
		numpy_table=np.load(npytable,allow_pickle=True)
		bin_header=numpy_table[2]
		data=numpy_table[1]
		header=numpy_table[0]
		nint=int(header[2])
		nant=int(header[3])
		nchan=int(header[4])
		for i in range(nint):
			for j in range(nant):
				for k in range(nchan):
					data[:,:,i,j,k]=np.matmul(poldistortion_matrix,np.matrix(data[:,:,i,j,k])) # Correcting poldistortion, CALIBRATE caltables has inverses of Jones matrices
		numpy_table[1]=data
		np.save(npytable,numpy_table)
		if outfile_path=='':
			outputfile=gaintable_path+'/'+outfile
		elif os.path.isdir(outfile_path)==False:
			outputfile=gaintable_path+'/'+os.path.basename(outfile)
		else:
			outputfile=outfile
		outputfile_bin,bad_flags=cal.convert_gaintable_npy2bin(npytable,outputfile+'.calibrate_bin',remove_nan=False)
		bin_data=np.fromfile(outputfile_bin,dtype=np.float64)
		data=np.load(gaintable,allow_pickle=True)
		data[0]=bin_data
		np.save(outputfile+'.temp',data)
		os.system('mv '+outputfile+'.temp.npy '+outputfile)
		os.system('rm -rf '+npytable+' *.calibrate_bin*')
		os.system('rm -rf casa*log')
		return outputfile

	def remove_model_negative(self,imagename,modelname,sigma=10,overwrite=False):
		'''
		Function to remove negatives from model image

		Parameters
		----------
		imagename : str 
			Name of the image
		modelname : str 
			Name of the model
		sigma : float 
			Sigma value for thresholding
		overwrite : bool 
			Overwrite the model image or not
		Returns
		-------
		str					
			Model imagename without negatives
		'''
		ia=image()
		if overwrite==False:
			if os.path.isdir(modelname+'.nonegative')==True:
				os.system('rm -rf '+modelname+'.nonegative')
			os.system('cp -r '+modelname+' '+modelname+'.nonegative')
			modelname=modelname+'.nonegative'
		ia.open(imagename)
		data=ia.getchunk()
		ia.close()
		rmsi=imstat(imagename=imagename,box=self.rms_box,stokes='I')['rms'][0]		
		sigma_pos=np.where(data[:,:,0,:]>=(sigma*rmsi))
		sigma_pos1=np.where(data[:,:,0,:]<(sigma*rmsi))
		ia.open(modelname)
		data=ia.getchunk()
		data_copy=copy.deepcopy(data)
		data[:,:,0,:][sigma_pos]=0
		pos=np.where(data[:,:,0,:]<0)
		data_copy[:,:,0,:][pos]=0
		data_copy[:,:,0,:][sigma_pos1]=0
		ia.putchunk(data_copy)
		ia.done()
		return modelname

	def file_remover_and_keeper(self,num_iter,msg_code,ref_time_chan=True):
		'''
		This function keep and remove caltables, ms, imaging related files based on the need

		Parameters
		----------
		num_iter : int
			Number of self-calibration iteration
		msg_code : str 
			Selfcal message code
		ref_timechan : bool
			Reference time channel or not
		'''
		msname_str=am.splited_ms_rename(self.msname,ref_time_chan=ref_time_chan,change_msname=False)
		freqstr=os.path.basename(msname_str).split('.ms')[0].split('_freq_')[1].split('_')[0]  # Frequency string in MHz
		datestr_list=os.path.basename(msname_str).split('.ms')[0].split('_freq_')[0].split('time_')[1].split('_')
		datestr_for_file='_'.join(datestr_list[:3])+'_'+'_'.join(datestr_list[3:]) # Datetime string 
		cwd=os.getcwd()
		file_str=os.path.basename(self.msname).split('.ms')[0]+'_'+str(num_iter) # File string prefix
		caltable_name=self.msname.split('.ms')[0]+'.bin' # Caltable name
		file_str_prefix='freq_'+freqstr+'_datetime_'+datestr_for_file+'_pol'
		if self.verbose==True and os.path.isdir(self.mspath+'/'+file_str_prefix)==False: # If verbose=True, directory to keep all intermediate images, caltables, models, residuals 
			os.mkdir(self.mspath+'/'+file_str_prefix)
		if self.verbose:
			if os.path.isdir(self.mspath+'/'+file_str_prefix+'/backup_ms')==False:
				os.makedirs(self.mspath+'/'+file_str_prefix+'/backup_ms')
			if os.path.isdir(self.mspath+'/'+file_str_prefix+'/backup_imagemodel')==False:
				os.makedirs(self.mspath+'/'+file_str_prefix+'/backup_imagemodel')
			if os.path.isdir(self.mspath+'/'+file_str_prefix+'/backup_cal')==False:
				os.makedirs(self.mspath+'/'+file_str_prefix+'/backup_cal')
		if num_iter==0 and (msg_code==0 or msg_code==8 or msg_code==9):
			os.system('cp -r '+caltable_name+' junk0.bin') # Copying num_iter=0 caltable to junk0.bin
			os.system('cp -r '+self.msname+' junk0.ms') # Copying num_iter=0 ms to junk0.ms
			os.system('cp -r '+file_str+'.image junk0.image') # Copying num_iter=0 image to junk0.image
			os.system('cp -r '+file_str+'.model junk0.model') # Copying num_iter=0 model to junk0.model
			if os.path.isdir(file_str+'.mask'):
				os.system('cp -r '+file_str+'.mask junk0.mask') # Copying num_iter=0 mask to junk0.mask
			os.system('cp -r '+file_str+'.residual junk0.residual') # Copying num_iter=0 residual to junk0.residual
			if self.verbose and (msg_code==0 or msg_code==8 or msg_code==9):	
				os.system('cp -r '+caltable_name+' '+self.mspath+'/'+file_str_prefix+'/backup_cal/'+file_str+'.bin') # Verbose=True, keep all the caltables
				os.system('cp -r '+self.msname+' '+self.mspath+'/'+file_str_prefix+'/backup_ms/'+file_str+'.ms') # If Verbose=True, keep all the ms
				os.system('cp -r '+file_str+'.model '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.model') # If Verbose=True, keep all the models
				if os.path.isdir(file_str+'.mask'):
					os.system('cp -r '+file_str+'.mask '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.mask') # If Verbose=True, keep all the masks
				os.system('cp -r '+file_str+'.image '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.image') # If Verbose=True, keep all the image
				os.system('cp -r '+file_str+'.residual '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.residual') # If Verbose=True, keep all the residuals
			os.system('rm -rf '+caltable_name)
			os.system('rm -rf '+file_str+'.image '+file_str+'.model '+file_str+'.residual '+file_str+'.sumwt '+file_str+'.pb '+file_str+'.psf '+file_str+'.mask *psf*') 
				# Removing all imaging related files
		elif num_iter==1 and (msg_code==0 or msg_code==8 or msg_code==9):
			os.system('cp -r '+caltable_name+' junk1.bin') # Copying num_iter=1 caltable to junk1.bin
			os.system('cp -r '+self.msname+' junk1.ms') # Copying num_iter=1 ms to junk1.ms
			os.system('cp -r '+file_str+'.image junk1.image') # Copying num_iter=1 image to junk1.image
			os.system('cp -r '+file_str+'.model junk1.model') # Copying num_iter=1 model to junk1.model
			if os.path.isdir(file_str+'.mask'):
				os.system('cp -r '+file_str+'.mask junk1.mask') # Copying num_iter=1 model to junk1.mask
			os.system('cp -r '+file_str+'.residual junk1.residual') # Copying num_iter=1 residual to junk1.residual
			if self.verbose and (msg_code==0 or msg_code==8 or msg_code==9):	
				os.system('cp -r '+caltable_name+' '+self.mspath+'/'+file_str_prefix+'/backup_cal/'+file_str+'.bin') # Verbose=True, keep all the caltables
				os.system('cp -r '+self.msname+' '+self.mspath+'/'+file_str_prefix+'/backup_ms/'+file_str+'.ms') # If Verbose=True, keep all the ms
				os.system('cp -r '+file_str+'.model '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.model') # If Verbose=True, keep all the models
				if os.path.isdir(file_str+'.mask'):
					os.system('cp -r '+file_str+'.mask '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.mask') # If Verbose=True, keep all the masks
				os.system('cp -r '+file_str+'.image '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.image') # If Verbose=True, keep all the image
				os.system('cp -r '+file_str+'.residual '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.residual') # If Verbose=True, keep all the residuals
			os.system('rm -rf '+caltable_name)
			os.system('rm -rf '+file_str+'.image '+file_str+'.model '+file_str+'.residual '+file_str+'.sumwt '+file_str+'.pb '+file_str+'.psf '+file_str+'.mask *psf*') 
				# Removing all imaging related files
		elif num_iter>1 and (msg_code==0 or msg_code==8 or msg_code==9):
			if os.path.isdir('junk1.bin'):
				os.system('rm -rf junk0.bin')
				os.system('cp -r junk1.bin junk0.bin') # Move the previous round caltable to pre-previous round for num_iter>1
			else:
				os.system('cp -r '+caltable_name+' junk0.bin') # Copying caltable to junk0.bin
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
				if os.path.isdir(file_str+'.mask'):	
					os.system('cp -r '+file_str+'.residual junk0.residual') # Copying model to junk0.residual
			os.system('rm -rf junk1.ms junk1.bin junk1.model junk1.mask junk1.image junk1.residual')
			os.system('cp -r '+caltable_name+' junk1.bin') # Copying caltable to junk1.bin
			os.system('cp -r '+self.msname+' junk1.ms') # Copying ms to junk1.ms
			os.system('cp -r '+file_str+'.model junk1.model') # Copying model to junk1.model
			if os.path.isdir(file_str+'.mask'):
				os.system('cp -r '+file_str+'.mask junk1.mask') # Copying mask to junk1.mask
			os.system('cp -r '+file_str+'.image junk1.image') # Copying image to junk1.image
			os.system('cp -r '+file_str+'.residual junk1.residual') # Copying residual to junk1.resuidual
			if self.verbose and (msg_code==0 or msg_code==8 or msg_code==9):	
				if os.path.exists(self.mspath+'/'+file_str_prefix+'/backup_cal/'+file_str+'.bin'):
					os.system('rm -rf '+self.mspath+'/'+file_str_prefix+'/backup_cal/'+file_str+'.bin')
				os.system('cp -r '+caltable_name+' '+self.mspath+'/'+file_str_prefix+'/backup_cal/'+file_str+'.bin') # Verbose=True, keep all the caltables
				if os.path.exists(self.mspath+'/'+file_str_prefix+'/backup_ms/'+file_str+'.ms'):
					os.system('rm -rf '+self.mspath+'/'+file_str_prefix+'/backup_ms/'+file_str+'.ms')
				os.system('cp -r '+self.msname+' '+self.mspath+'/'+file_str_prefix+'/backup_ms/'+file_str+'.ms') # If Verbose=True, keep all the ms
				if os.path.exists(self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.model'):	
					os.system('rm -rf '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.model')
				os.system('cp -r '+file_str+'.model '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.model') # If Verbose=True, keep all the models
				if os.path.exists(self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.mask'):
					os.system('rm -rf '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.mask')
				if os.path.isdir(file_str+'.mask'):
					os.system('cp -r '+file_str+'.mask '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.mask') # If Verbose=True, keep all the masks
				if os.path.exists(self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.image'):
					os.system('rm -rf '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.image')
				os.system('cp -r '+file_str+'.image '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.image') # If Verbose=True, keep all the image
				if os.path.exists(self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.residual'):
					os.system('rm -rf '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.residual')
				os.system('cp -r '+file_str+'.residual '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.residual') # If Verbose=True, keep all the residuals
			os.system('rm -rf '+file_str+'.image '+file_str+'.model '+file_str+'.residual '+file_str+'.sumwt '+file_str+'.pb '+file_str+'.psf '+file_str+'.mask *psf*') 
				# Removing all imaging related files
		os.chdir(self.cwd)
		os.system('rm -rf casa*log')
		return

	def estimateSkyBrightnessMatrix(self,beam_jones,Vij):
		''' 
		Return beam corrected brightness matrix

		Parameters
		----------
		beam_jones : numpy.array
			Beam Jones matrix
		Vij : numpy.array
			Instrumental brightness matrix
		Returns
		-------
		numpy.array
			Beam corrected brightness matrix in instrumental basis
		'''
		J=np.matrix(beam_jones)
		J1=np.array(inv(J))
		J2=np.array(inv(J.H))
		J1_11=J1[0,0]
		J1_12=J1[0,1]
		J1_21=J1[1,0]
		J1_22=J1[1,1]

		J2_11=J2[0,0]
		J2_12=J2[0,1]
		J2_21=J2[1,0]
		J2_22=J2[1,1]

		XX=Vij[0,0]
		XY=Vij[0,1]
		YX=Vij[1,0]
		YY=Vij[1,1]
		XX_out=J2_11*(J1_11*XX+J1_12*YX)+J2_21*(J1_11*XY+J1_12*YY)
		XY_out=J2_12*(J1_11*XX+J1_12*YX)+J2_22*(J1_11*XY+J1_12*YY)
		YX_out=J2_11*(J1_21*XX+J1_22*YX)+J2_21*(J1_21*XY+J1_22*YY)
		YY_out=J2_12*(J1_21*XX+J1_22*YX)+J2_22*(J1_21*XY+J1_22*YY)
		B_corr=np.array([[XX_out,XY_out],[YX_out,YY_out]])
		os.system('rm -rf casa*log')
		return B_corr
	
	def get_IQUV(self,imagename,imagetype='FITS'):
		'''
		Stokes I,Q,U,V from a Stokes IQUV image cube
	
		Parameters
		----------
		imagename : str 
			Name of the image
		imagetype : str 
			Type of the image, CASA or FITS
		Returns
		-------
		dict
			{'STOKES':imagedata}
		'''
		if imagetype=='CASA':
			if os.path.isfile(imagename.split('.image')[0]+'.fits'):
				os.system('rm -rf '+imagename.split('.image')[0]+'.fits')
			exportfits(imagename=imagename,fitsimage=imagename.split('.image')[0]+'.fits')
		fitsimage=imagename.split('.image')[0]+'.fits'
		data=fits.getdata(fitsimage)
		stokes = {}
		stokes['I'] = data[0, 0, :, :]
		stokes['Q'] = data[1, 0, :, :]
		stokes['U'] = data[2, 0, :, :]
		stokes['V'] = data[3, 0, :, :]
		os.system('rm -rf '+fitsimage)
		os.system('rm -rf casa*log')
		return stokes	

	def get_inst_pols(self,stokes_image,imagetype='FITS',pol_basis='Linear'): #TODO : Circular basis
		'''
		Return instrumental polarisation matrix (Vij)

		Parameters
		----------
		stokes_image : str 
			Name of the Stokes IQUV image cube
		imagetype : str
			Type of the image, CASA or FITS
		pol_basis : str 
			Polarisation basis of the instrument, Linear or Circular (Circular basis not implemented)
		Returns
		-------
		numpy.array
			Instrumental polarisation matrix
		'''
		stokes=self.get_IQUV(stokes_image,imagetype=imagetype)
		XX = stokes['I'] + stokes['Q']
		XY = stokes['U'] + stokes['V'] * 1j
		YX = stokes['U'] - stokes['V'] * 1j
		YY = stokes['I'] - stokes['Q']
		Vij = np.array([[XX, XY], [YX, YY]])
		os.system('rm -rf casa*log')
		return Vij

	def B_to_IQUV(self,B,pol_basis='Linear'): # TODO : Circular Basis
		'''
		Convert brightness matrix in instrumental basis to I, Q, U, V

		Parameters
		----------
		B : numpy.array 
			Brightness matrix in instrumental basis
		pol_basis : str 
			Polarisation basis of the instrument, Linear or Circular (Circular basis not implemented)
		Returns
		-------
		dict
			Stokes I, Q, U,V
		'''
		B11 = B[0, 0, :, :]
		B12 = B[0, 1, :, :]
		B21 = B[1, 0, :, :]
		B22 = B[1, 1, :, :]
		stokes = {}
		stokes['I'] = (B11 + B22) / 2.
		stokes['Q'] = (B11 - B22) / 2.
		stokes['U'] = (B12 + B21) / 2.
		stokes['V'] = 1j * (B21 - B12) / 2.
		os.system('rm -rf casa*log')
		return stokes

	def correct_for_single_beam_jones(self,imagename,outfile,beam_jones,imagetype='FITS',outtype='FITS',pol_basis='Linear'): #TODO : Circular basis
		'''
		Correct Stokes IQUV image cube for full Stokes Beam Jones at a single pointing

		Parameters
		----------
		imagename : str 
			Name of the image of model
		outfile : str 
			Name of the beam corrected image or model
		beam_jones : numpy.array 
			Beam jones matrix
		imagetype : str 
			Type of the image, CASA or FITS
		outtype : str 
			Output image type, CASA or FITS
		pol_basis : str 
			Polarisation basis of the instrument, Linear or Circular (Circular basis not implemented)
		Returns
		-------
		str
			Beam corrected image or model
		'''
		Vij=self.get_inst_pols(imagename,imagetype=imagetype,pol_basis=pol_basis)
		B=self.estimateSkyBrightnessMatrix(beam_jones,Vij)
		stokes=self.B_to_IQUV(B,pol_basis='Linear')
		if os.path.exists(outfile):
			os.system('rm -rf '+outfile)
			os.system('rm -rf '+'.'.join(imagename.split('.')[:-1])+'.fits')
		if imagetype=='CASA' and os.path.isfile('.'.join(imagename.split('.')[:-1])+'.fits')==False:
			exportfits(imagename=imagename,fitsimage='.'.join(imagename.split('.')[:-1])+'.fits',stokeslast=False)
			os.system('cp -r '+imagename+' temp_org.image')
		imagename='.'.join(imagename.split('.')[:-1])+'.fits'
		data=fits.getdata(imagename)
		header=fits.getheader(imagename)
		data[0,0,:,:]=np.real(stokes['I'])
		data[0,1,:,:]=np.real(stokes['Q'])
		data[0,2,:,:]=np.real(stokes['U'])
		data[0,3,:,:]=np.real(stokes['V'])
		fits.writeto(outfile,data=data,header=header,overwrite=True)
		if outtype=='CASA':
			if os.path.exists('temp.image'):
				os.system('rm -rf temp.image')
			importfits(fitsimage=outfile,imagename='temp.image')
			os.system('rm -rf '+outfile+' '+imagename)
			ia=image()
			ia.open('temp.image')
			pbcor_data=ia.getchunk()
			ia.close()
			ia.open('temp_org.image')
			ia.putchunk(pbcor_data)
			ia.done()
			ia.close()
			os.system('mv temp_org.image '+outfile)
		os.system('rm -rf casa*log temp*')
		return outfile

	def correct_image_for_cross_phase(self,imagename,modelname,outfile,cross_phase=15,imagetype='FITS',outtype='FITS',pol_basis='Linear',do_fluxcal=False):
		'''
		Correct Stokes IQUV image cube for full Stokes Beam Jones at a single pointing

		Parameters
		----------
		imagename : str 
			Name of the image
		modelname : str 
			Name of the model
		outfile : str 
			Prefix name of the beam corrected image and model
		cross_phase : str 
			Cross hand phase in degree
		imagetype : str 
			Type of the image, CASA or FITS
		outtype : str 
			Output image type, CASA or FITS
		pol_basis : str 
			Polarisation basis of the instrument, Linear or Circular
		do_fluxcal : str 
			Perform flux scaling or not
		Returns
		-------
		str
			Cross hand phase corrected image and model
		'''
		cross_phase=np.deg2rad(cross_phase)/2.0
		cross_jones=np.matrix([[np.cos(cross_phase)+1j*np.sin(cross_phase),0],[0,np.cos(cross_phase)-1j*np.sin(cross_phase)]])
		outfile_image=self.correct_for_single_beam_jones(imagename,outfile+'.image',cross_jones,imagetype=imagetype,outtype=outtype,pol_basis=pol_basis)
		outfile_model=self.correct_for_single_beam_jones(modelname,outfile+'.model',cross_jones,imagetype=imagetype,outtype=outtype,pol_basis=pol_basis)
		os.system('rm -rf casa*log')
		return outfile_image,outfile_model

	def uncorrect_for_single_beam_jones(self,imagename,outfile,inv_beam_jones,imagetype='FITS',outtype='FITS',pol_basis='Linear'): # TODO : Circular basis
		'''
		Undo the beam correction for Stokes IQUV image cube for full Stokes Beam Jones at a single pointing

		Parameters
		----------
		imagename : str 
			Name of the image of model
		outfile : str 
			Name of the beam corrected image or model
		inv_beam_jones : numpy.array 
			Inverse of Beam jones matrix
		imagetype : str 
			Type of the image, CASA or FITS
		outtype : str 
			Output image type, CASA or FITS
		pol_basis : str 
			Polarisation basis of the instrument, Linear or Circular
		Returns
		-------
		str
			Beam un-corrected image or model
		'''
		Vij=self.get_inst_pols(imagename,imagetype=imagetype,pol_basis=pol_basis)
		B=self.estimateSkyBrightnessMatrix(inv_beam_jones,Vij)
		stokes=self.B_to_IQUV(B,pol_basis='Linear')
		if os.path.exists(outfile):
			os.system('rm -rf '+outfile)
			os.system('rm -rf '+'.'.join(imagename.split('.')[:-1])+'.fits')
		if imagetype=='CASA' and os.path.isfile('.'.join(imagename.split('.')[:-1])+'.fits')==False:
			exportfits(imagename=imagename,fitsimage='.'.join(imagename.split('.')[:-1])+'.fits',stokeslast=False)
			os.system('cp -r '+imagename+' temp_org.image')
		imagename='.'.join(imagename.split('.')[:-1])+'.fits'
		data=fits.getdata(imagename)
		header=fits.getheader(imagename)
		data[0,0,:,:]=np.real(stokes['I'])
		data[0,1,:,:]=np.real(stokes['Q'])
		data[0,2,:,:]=np.real(stokes['U'])
		data[0,3,:,:]=np.real(stokes['V'])
		fits.writeto(outfile,data=data,header=header,overwrite=True)
		if outtype=='CASA':
			importfits(fitsimage=outfile,imagename='temp.image')
			os.system('rm -rf '+outfile+' '+imagename)
			ia=image()
			ia.open('temp.image')
			pbcor_data=ia.getchunk()
			ia.close()
			ia.open('temp_org.image')
			ia.putchunk(pbcor_data)
			ia.done()
			ia.close()
			os.system('mv temp_org.image '+outfile)
		os.system('rm -rf casa*log temp*')
		return outfile

	def uncorrect_visibility_model_single_beam_jones(self,force=False,skip_freq=1.28):	
		'''
		Undo Correct visibility data for a single pointing beam jones

		Parameters
		----------
		force : bool 
			Undo beam correct forcefully avoiding ms header info
		skip_freq : float 
			Frequency interval in MHz to make independent beams (default : 1.28 MHz). If anything greater than 1.28 MHz is given it will be overwritten to 1.28 MHz
		Returns
		-------
		str
			Name of the beam jones file
		'''
		mwapb=MWA_PrimaryBeam(self.msname,self.metafits,inverse_beam=True)
		cal=CALIBRATE()
		beamfile=self.msname+'.beam.bin'
		beamfile,beamjones=mwapb.MWA_phasecenter_beam_jones(outputfile=beamfile,skip_freq=float(skip_freq))
		tb=table()
		tb.open(self.msname,nomodify=False)
		data=tb.getcol('DATA')
		try:
			cor_data=tb.getcol('CORRECTED_DATA')
		except:
			cor_data=data
		model_data=tb.getcol('MODEL_DATA')
		tb.putcol('DATA',model_data)
		tb.flush()
		tb.close()
		cal.applycal(msname=self.msname,gaintable=beamfile,applymode='calonly') # Applying the inverse beam correction
		tb.open(self.msname,nomodify=False)
		model_data=tb.getcol('CORRECTED_DATA')
		tb.putcol('MODEL_DATA',model_data)
		try:
			tb.putcol('CORRECTED_DATA',cor_data)
		except:
			pass
		tb.putcol('DATA',data)
		tb.flush()
		tb.close()
		del data,cor_data,model_data
		os.system('rm -rf casa*log')
		return beamfile

	def correct_visibility_single_beam_jones(self,datacolumn='DATA',modify_datacolumn=True,force=False,skip_freq=1.28,save_beamfile=''):
		'''
		Correct visibility data for a single pointing beam jones

		Parameters
		----------
		datacolumn : str 
			'DATA', datacolumn to apply beam correction
		modify_datacolumn :bool
			Modify the DATA column, otherwise beam corrected visibilities will be saved on CORRECTED_DATA
		force : bool 
			Beam correct forcefully avoiding ms header info
		skip_freq : float 
			Frequency interval in MHz to make independent beams (default : 1.28 MHz). If anything greater than 1.28 MHz is given it will be overwritten to 1.28 MHz
		save_beamfile : str 
			Save beam file in this given name
		Returns
		-------
		str
			Name of the beam jones file, Beam Jones matrix.
		'''
		if os.path.exists(save_beamfile)==True:
			os.system('rm -rf '+save_beamfile)
		mwapb=MWA_PrimaryBeam(self.msname,self.metafits,inverse_beam=False) 
		cal=CALIBRATE()
		beamfile=self.msname+'.beam.bin'
		beamfile,beamjones=mwapb.MWA_phasecenter_beam_jones(outputfile=beamfile,skip_freq=float(skip_freq))
		if save_beamfile!='' and save_beamfile!=beamfile and os.path.exists(save_beamfile)==False:
			os.system('cp -r '+beamfile+' '+save_beamfile)
		code=vishead(vis=self.msname,mode='get',hdkey='fld_code')[0][0]
		code_list=code.split(',')
		if 'S_PBCOR' not in code_list or 'S_PBUNCOR' in code_list:
			cal.applycal(msname=self.msname,gaintable=beamfile,datacolumn=datacolumn,applymode='calonly') # Applying the beam correction
			if modify_datacolumn==True:
				tb=table()
				tb.open(self.msname,nomodify=False)
				cor_data=tb.getcol('CORRECTED_DATA')
				tb.putcol('DATA',cor_data)
				tb.flush()
				tb.close()
				self.pollog_verbose.info('Modified DATA column.\n')
			if save_beamfile!=beamfile:
				os.system('rm -rf '+self.msname+'.beam*')
			if len(code_list)==1 and code_list[0]=='':
				code+='S_PBCOR'
			else:
				code+=',S_PBCOR'
			vishead(vis=self.msname,mode='put',hdkey='fld_code',hdvalue=np.array([code]))
			os.system('rm -rf casa*log')
			self.pollog_verbose.info('Beam correction done. Beam file is at : '+beamfile+'\n')
			return beamfile,beamjones
		elif force==True:
			cal.applycal(msname=self.msname,gaintable=beamfile,datacolumn=datacolumn,applymode='calonly') # Applying the beam correction
			if modify_datacolumn==True:
				tb=table()
				tb.open(self.msname,nomodify=False)
				cor_data=tb.getcol('CORRECTED_DATA')
				tb.putcol('DATA',cor_data)
				tb.flush()
				tb.close()
				self.pollog_verbose.info('Modified DATA column.\n')
			if save_beamfile!=beamfile:
				os.system('rm -rf '+self.msname+'.beam*')
			os.system('rm -rf casa*log')
			self.pollog_verbose.info('Beam correction done. Beam file is at : '+beamfile+'\n')
			return beamfile,beamjones
		else:
			self.pollog_verbose.info('Beam correction has already been applied.\n')
			os.system('rm -rf casa*log')
			return beamfile,beamjones

	def uncorrect_visibility_single_beam_jones(self,modify_datacolumn=True,force=False,skip_freq=1.28):	
		'''
		Undo Correct visibility data for a single pointing beam jones

		Parameters
		----------
		modify : bool 
			Modify the DATA column, otherwise beam corrected visibilities will be saved on CORRECTED_DATA
		force : bool 
			Undo beam correct forcefully avoiding ms header info
		skip_freq : float 
			Frequency interval in MHz to make independent beams (default : 1.28 MHz). If anything greater than 1.28 MHz is given it will be overwritten to 1.28 MHz
		Returns
		-------
		str
			Name of the beam jones file
		'''
		mwapb=MWA_PrimaryBeam(self.msname,self.metafits,inverse_beam=True)
		cal=CALIBRATE()
		beamfile=self.msname+'.beam.bin'
		beamfile,beamjobes=mwapb.MWA_phasecenter_beam_jones(outputfile=beamfile,skip_freq=float(skip_freq))
		cal.applycal(msname=self.msname,gaintable=beamfile,applymode='calonly') # Applying the inverse beam correction
		code=vishead(vis=self.msname,mode='get',hdkey='fld_code')[0][0]
		code_list=code.split(',')
		if modify_datacolumn==True:
			if 'S_PBCOR' in code_list:
				if os.path.isdir(self.msname+'.beam.ms'):
					os.system('rm -rf '+self.msname+'.beam.ms')
				split(vis=self.msname,outputvis=self.msname+'.beam.ms',datacolumn='corrected')
				if self.msname[-1]=='/':
					msname=self.msname[:-1]
				else:
					msname=self.msname
				os.system('rm -rf '+msname)
				os.system('mv '+self.msname+'.beam.ms '+msname)
				os.system('rm -rf '+self.msname+'.beam*')
				if len(code_list)==1 and code_list[0]=='':
					code+='S_PBUNCOR'
				else:
					code+=',S_PBUNCOR'
				vishead(vis=self.msname,mode='put',hdkey='fld_code',hdvalue=np.array([code]))
				self.pollog_verbose.info('Undo beam correction\n')
			elif force==True:
				if os.path.isdir(self.msname+'.beam.ms'):
					os.system('rm -rf '+self.msname+'.beam.ms')
				split(vis=self.msname,outputvis=self.msname+'.beam.ms',datacolumn='corrected')
				if self.msname[-1]=='/':
					msname=self.msname[:-1]
				else:
					msname=self.msname
				os.system('rm -rf '+msname)
				os.system('mv '+self.msname+'.beam.ms '+msname)
				os.system('rm -rf '+self.msname+'.beam*')
				self.pollog_verbose.info('Undo beam correction\n')
			else:
				self.pollog_verbose.info('No beam correction was done on this measurement set. Thus not undoing any beam correction.\n')
				if self.verbose==False:
					print('No beam correction was done on this measurement set. Thus not undoing any beam correction.\n')
		os.system('rm -rf casa*log')
		return beamfile

	def IMSTAT_record(self,DRI,DR_neg,FXI,FXQ,FXU,FXV,FXT,FXP,record_filename,init=True):
		'''
		Function to keep the record of image statistics at different self calibration steps

		Parameters
		----------
		DRI : float 
			RMS based dynamic range of the Stokes I
		DR_neg : float 
			Negative based dynamic range of the Stokes I
		FXI : float 
			Total Stokes I flux
		FXQ : float 
			Total Stokes Q flux
		FXU : float 
			Total Stokes U flux
		FXV : float 
			Total Stokes V flux
		FXT : float 
			Total Stokes T flux
		FXP : float 
			Total Stokes P flux
		record_filename : str 
			Name of the file to stro dynamic ranges
		init : bool
			Initiating a new record from the current selfcal iteration
		Returns
		-------
		numpy.array
			Image statistic record array; shape [7,num_of_record]
		'''
		if init==True:
			os.system('rm -rf '+record_filename+'.npy')
			IMSTAT_array=np.empty([7,1])
			IMSTAT_array=np.array([DRI,FXI,FXQ,FXU,FXV,FXT,FXP]).reshape(7,1)
			np.save(record_filename,IMSTAT_array)
		else:
			if os.path.isfile(record_filename+'.npy')==False:
				IMSTAT_array=np.empty([7,1])
				IMSTAT_array=np.array([DRI,FXI,FXQ,FXU,FXV,FXT,FXP]).reshape(7,1)
			else:
				IMSTAT_array=np.load(record_filename+'.npy')
				DR_array=np.array([DRI,FXI,FXQ,FXU,FXV,FXT,FXP]).reshape(7,1)
				IMSTAT_array=np.append(IMSTAT_array,DR_array,axis=1)
			np.save(record_filename,IMSTAT_array)
		os.system('rm -rf casa*log')
		return IMSTAT_array

	def reduce_sigma(self,imagename,nsigma,sigma_step,minsigma,pre_residual=0.0,residual_frac=0.01,stokes_list=['I']):
		'''
		Function to determine whether reduce the CLEAN sigma or not

		Parameters
		----------
		imagename : str 
			Name of the image
		nsigma : float 
			Value of the present n-sigma
		sigma_step : float 
			Step to reduce sigma value
		minsigma : float 
			Minimum allowed sigma
		pre_residual : float 
			Previous residual fraction to compare (default : 0.0)
		residual_frac : float 
			Residual flux fraction to reduce sigma (default : 0.01)
		stokes_list : list 		
			Stokes plane list
		Returns
		-------
		float
			Reduced value of n-sigma and median residual fraction if residual flux is more than given percentage (default : 1%) of the total flux in Stokes I or in all Stokes Q,U,V.
		float
			Median residual flux fraction over all stokes planes
		'''
		imagename=imagename
		residual=imagename.split('.image')[0]+'.residual'
		do_reduce_list=[]
		residual_frac_list=[]
		ia=image()
		imagename_path=os.path.dirname(os.path.realpath(imagename))
		cwd=os.getcwd()
		if imagename_path!='':
			os.chdir(imagename_path)
		os.system('rm -rf reduce_sigma_*')
		for stokes in stokes_list:
			if stokes=='I' or stokes=='XX' or stokes=='YY':
				if os.path.isdir('reduce_sigma_I.image')==True:
					os.system('rm -rf reduce_sigma_I.image')
				if os.path.isdir('reduce_sigma_I.residual')==True:
					os.system('rm -rf reduce_sigma_I.residual')
				self.pollog_verbose.info('immath(imagename=\''+imagename+'\',mode=\'evalexpr\',stokes=\''+stokes+'\',outfile=\'reduce_sigma_I.image\')\n')
				immath(imagename=imagename,mode='evalexpr',stokes=stokes,outfile='reduce_sigma_I.image')
				self.pollog_verbose.info('immath(imagename=\''+residual+'\',mode=\'evalexpr\',stokes=\''+stokes+'\',outfile=\'reduce_sigma_I.residual\')\n')
				immath(imagename=residual,mode='evalexpr',stokes=stokes,outfile='reduce_sigma_I.residual')
				self.pollog_verbose.info('imstat(imagename=\''+imagename+'\',box=\''+self.rms_box+'\',stokes=\''+stokes+'\')[\'rms\'][0]\n')
				rms=imstat(imagename='reduce_sigma_I.image',box=self.rms_box,stokes=stokes)['rms'][0]
				ia.open('reduce_sigma_I.image')			
				ia.calcmask('\"reduce_sigma_I.image\">'+str(nsigma*rms),'mymask')
				ia.close()
				makemask(inpimage='reduce_sigma_I.image',inpmask='reduce_sigma_I.image:mymask',output='reduce_sigma_I.residual:mymask',mode='copy')
				try:
					image_pix_sum=imstat(imagename='reduce_sigma_I.image')['sum'][0]
					residual_pix_sum=imstat(imagename='reduce_sigma_I.residual')['sum'][0]
				except:
					image_pix_sum=1
					residual_pix_sum=0
				try:
					maxval=imstat(imagename='reduce_sigma_I.residual')['max'][0]
				except:
					maxval=0
				try:
					minval=imstat(imagename='reduce_sigma_I.residual')['min'][0]
				except:
					minval=0					
				if maxval>(nsigma-sigma_step)*rms:
					max_frac_diff=(maxval-(nsigma-sigma_step)*rms)/maxval
				else:
					max_frac_diff=0
				min_frac_diff=0
			else:
				self.pollog_verbose.info('immath(imagename=\''+imagename+'\',mode=\'evalexpr\',stokes=\''+stokes+'\',expr=\'abs(IM0)\',outfile=\'reduce_sigma_'+stokes+'.image\')\n')
				immath(imagename=imagename,mode='evalexpr',stokes=stokes,expr='abs(IM0)',outfile='reduce_sigma_'+stokes+'.image')
				self.pollog_verbose.info('immath(imagename=\''+residual+'\',mode=\'evalexpr\',stokes=\''+stokes+'\',outfile=\'reduce_sigma_'+stokes+'.residual\')\n')
				immath(imagename=residual,mode='evalexpr',stokes=stokes,outfile='reduce_sigma_'+stokes+'.residual')
				self.pollog_verbose.info('imstat(imagename=\''+imagename+'\',box=\''+self.rms_box+'\',stokes=\''+stokes+'\')[\'rms\'][0]\n')
				rms=imstat(imagename='reduce_sigma_'+stokes+'.image',box=self.rms_box,stokes=stokes)['rms'][0]
				ia.open('reduce_sigma_'+stokes+'.image')			
				ia.calcmask('\"reduce_sigma_'+stokes+'.image\">'+str(nsigma*rms),'mymask')
				ia.close()
				makemask(inpimage='reduce_sigma_'+stokes+'.image',inpmask='reduce_sigma_'+stokes+'.image:mymask',output='reduce_sigma_'+stokes+'.residual:mymask',mode='copy')
				try:
					image_pix_sum=imstat(imagename='reduce_sigma_'+stokes+'.image')['sum'][0]
					residual_pix_sum=imstat(imagename='reduce_sigma_'+stokes+'.residual')['sum'][0]
				except:
					image_pix_sum=1
					residual_pix_sum=0
				try:
					maxval=imstat(imagename='reduce_sigma_'+stokes+'.residual')['max'][0]
				except:
					maxval=0
				try:
					minval=imstat(imagename='reduce_sigma_'+stokes+'.residual')['min'][0]
				except:
					minval=0					
				if maxval>(nsigma-sigma_step)*rms:
					max_frac_diff=(maxval-(nsigma-sigma_step)*rms)/maxval
				else:
					max_frac_diff=0
				if abs(minval)>(nsigma-sigma_step)*rms:
					min_frac_diff=(abs(minval)-(nsigma-sigma_step)*rms)/abs(minval)
				else:
					min_frac_diff=0
			residual_frac_list.append(residual_pix_sum/image_pix_sum)
			if (residual_pix_sum/image_pix_sum>residual_frac) or ((max_frac_diff>0 and max_frac_diff>residual_frac) or ((min_frac_diff>=0 and min_frac_diff>residual_frac) \
					and stokes!='I' and stokes!='XX' and stokes!='YY')):
				if (pre_residual>0 and residual_pix_sum/image_pix_sum<pre_residual) or pre_residual==0:
					do_reduce_list.append(1)
		os.system('rm -rf reduce_sigma_*')
		os.chdir(cwd)
		residual_frac_median=np.median(np.array(residual_frac_list))
		if int(np.sum(np.array(do_reduce_list)))>=1:
			if sigma_step>1.0:
				self.pollog_verbose.info('WARNING : Choosing sigma step 1 is too risky. Selfcal may diverge\n')
				if self.verbose==False:
					print ('WARNING : Choosing sigma step greater than 1 is too risky. Selfcal may diverge\n')
				if self.interactive==True:
					want_to_continue=input('Do you still want to continue? Y/y/Yes/yes')
					self.pollog_verbose.info('Interactive=True\n')
					self.pollog_verbose.info('Do you still want to continue? Y/y/Yes/yes:'+str(want_to_continue)+'\n')
					if want_to_continue=='Y' or want_to_continue=='y' or want_to_continue=='Yes' or want_to_continue=='yes':	
						self.pollog_verbose.info('Continuing with sigma step :'+str(sigma_step)+'\n')
						os.system('rm -rf casa*log')
						return nsigma-sigma_step,residual_frac_median
					else:
						os.system('rm -rf casa*log')
						return nsigma,residual_frac_median
			elif nsigma-sigma_step<minsigma:
				self.pollog_verbose.info('WARNING : Choosing sigma less than '+str(minsigma)+'\n')
				if self.verbose==False:
					print ('WARNING : Choosing sigma less than '+str(minsigma)+'\n')
				if self.interactive==True:
					want_to_continue=input('Do you still want to continue? Y/y/Yes/yes')
					self.pollog_verbose.info('Interactive=True\n')
					self.pollog_verbose.info('Do you still want to continue? Y/y/Yes/yes:'+str(want_to_continue)+'\n')
					if want_to_continue=='Y' or want_to_continue=='y' or want_to_continue=='Yes' or want_to_continue=='yes':		
						self.pollog_verbose.info('Continuing with sigma step :'+str(sigma_step)+'\n')
						os.system('rm -rf casa*log')
						return nsigma-sigma_step,residual_frac_median
					else:
						os.system('rm -rf casa*log')
						return nsigma,residual_frac_median
				else:
					self.pollog_verbose.info('Interactive=False\n')
					self.pollog_verbose.info('Continuing with sigma step :'+str(sigma_step)+'\n')
					if self.verbose==False:
						print ('Continuing with sigma step :'+str(sigma_step)+'\n')
					os.system('rm -rf casa*log')
					return nsigma-sigma_step,residual_frac_median
			else:
				self.pollog_verbose.info('Reducing sigma to:'+str(nsigma-sigma_step)+', because residual flux is more than '+str(residual_frac*100)+' %\n')
				os.system('rm -rf casa*log')
				return nsigma-sigma_step,residual_frac_median
		else:
			self.pollog_verbose.info('Sigma value is not changed, because residual flux less more than '+str(residual_frac*100)+' %. Sigma is at :'+str(nsigma)+'\n')
			os.system('rm -rf casa*log')
			return nsigma,residual_frac_median

	def solarlin_pol_minimise(self,datai,datal,l,rmsl,i_flux):
		'''
		Polarisation minimisation function for Sun

		Parameters
		----------
		datai : numpy.array 
			Stokes I image data
		datal : numpy.array 
			Stokes Q or U image data	
		l : float 
			Trial leakage (-1 to 1)
		rmsl : float 
			RMS of the Stokes map
		i_flux : float 
			Mean brightness in Jy/beam
		Returns
		-------
		int
			The number of pixels having polarisation fraction greater than rmsl/i_flux
		'''
		x1=np.abs((datal-l*datai)/datai)
		pos=np.where(x1.flatten()>rmsl/i_flux)
		f_out=(len(pos[0]))
		del x1,datai,datal,l,rmsl
		return f_out

	def solarcir_pol_minimise(self,datai,datav,l,rmsv,mean_i_flux,sigma):
		'''
		Polarisation minimisation function for Sun

		Parameters
		----------
		datai : numpy.array 
			Stokes I image data
		datav : numpy.array 
			Stokes V image data	
		l : float 
			Trial leakage (-1 to 1)
		rmsv : float 
			RMS of the Stokes V map
		mean_i_flux : float 
			Mean brightness in Jy/beam
		sigma : float 
			Sigma value for thresholding
		Return:
		int
			The number of pixels having polarisation fraction greater than rmsl/i_flux
		'''
		x1=np.abs((datav-l*datai)/datai)
		if (sigma*rmsv)/mean_i_flux>0.1:
			threshold=(sigma*rmsv)/mean_i_flux
		else:
			threshold=0.1
		pos=np.where(x1.flatten()>threshold)
		f_out=(len(pos[0]))
		del x1,datai,datav,l
		return f_out

	def subtract_leakage_surface(self,imagename,modelname,sigma=10,do_fluxcal=False,overwrite=False):
		'''
		Function to subtract quadratic leakage surface

		Parameters
		----------
		imagename : str 
			Name of the image
		modelname : str 
			Name of the model
		sigma : float 
			N-sigma value above which any emission is considered to be real
		do_fluxcal : bool
			Perform polynomial based flux calibration (See details Kansabanik et al. 2021, submitted to ApJ)
		overwrite : bool
			Overwrite the input image and model
		Returns
		-------
		str
			Leakage surface subtracted image
		str 
			Leakage surface subtracted model
		float 
			Stokes Q fractional change
		float
			Stokes U fractional change
		'''
		self.pollog_verbose.info('Correcting Stokes I to Stokes Q, U leakage surface.\n')
		if os.path.exists(imagename.split('.image')[0]+'_quvcor_surface.image'):
			os.system('rm -rf '+imagename.split('.image')[0]+'_quvcor_surface.image')
		if os.path.exists(modelname.split('.model')[0]+'_quvcor_surface.model'):
			os.system('rm -rf '+modelname.split('.model')[0]+'_quvcor_surface.model')
		if overwrite==False:
			os.system('cp -r '+imagename+' '+imagename.split('.image')[0]+'_quvcor_surface.image')
			os.system('cp -r '+modelname+' '+modelname.split('.model')[0]+'_quvcor_surface.model')
			imagename=imagename.split('.image')[0]+'_quvcor_surface.image'
			modelname=modelname.split('.model')[0]+'_quvcor_surface.model'
		outfile_path=os.path.dirname(os.path.realpath(imagename))
		os.system('rm -rf '+outfile_path+'/I*')
		imsubimage(imagename=imagename,outfile=outfile_path+'/I.image',stokes='I',dropdeg=False)
		if do_fluxcal==True:
			fluxcal_image=self.mwa_solar_fluxcal(imagename=outfile_path+'/I.image',outfile=outfile_path+'/I_fluxcal.image')
		else:
			fluxcal_image=outfile_path+'/I.image'
		major=imhead(imagename=fluxcal_image)['restoringbeam']['major']['value'] # In arcsec
		minor=imhead(imagename=fluxcal_image)['restoringbeam']['minor']['value'] # In arcsec
		freq=imhead(imagename=fluxcal_image)['refval'][-1]/10**9 # In GHz
		ref_beam_axes_multi=500000
		beam_axes_multi=major*minor	
		scale_tb_limit=(beam_axes_multi/ref_beam_axes_multi)
		expr='1.222e6*IM0/'+str(freq)+'^2/('+str(major*minor)+')'
		immath(imagename=fluxcal_image,outfile=outfile_path+'/I_Tb.image',mode='evalexpr',expr=expr)
		imhead(imagename=outfile_path+'/I_Tb.image', mode='put', hdkey='bunit', hdvalue='K')
		ia=image()
		ia.open(imagename)
		data=ia.getchunk()
		ia.close()
		a=image()
		ia.open(modelname)
		modeldata=ia.getchunk()
		ia.close()
		model_datai=modeldata[:,:,0,0]
		model_dataq=modeldata[:,:,1,0]
		model_datau=modeldata[:,:,2,0]
		rmsq=imstat(imagename=imagename,box=self.rms_box,stokes='Q')['rms'][0]
		rmsu=imstat(imagename=imagename,box=self.rms_box,stokes='U')['rms'][0]
		rmsi=imstat(imagename=imagename,box=self.rms_box,stokes='I')['rms'][0]
		i_flux=imstat(imagename=imagename,stokes='I')['flux'][0]
		dataq=data[:,:,1,0]
		datau=data[:,:,2,0]
		datai=data[:,:,0,0]
		datai_copy=copy.deepcopy(datai)
		dataq_copy=copy.deepcopy(dataq)
		datau_copy=copy.deepcopy(datau)
		dataq1=copy.deepcopy(dataq)
		datau1=copy.deepcopy(datau)
		posi=np.where(datai<(sigma*rmsi))
		maxpos=imstat(imagename=imagename,stokes='I')['maxpos']
		box=self.negative_box(maxpos,box_width=2).split(',')
		box_coords=[int(i) for i in box]
		posq=np.where(np.abs(dataq)<(sigma*rmsq))
		posu=np.where(np.abs(datau)<(sigma*rmsu))
		ia.open(outfile_path+'/I_Tb.image')
		datatb=ia.getchunk()
		ia.close()
		datai_mask=copy.deepcopy(datai)
		for i in range(datai_mask.shape[0]):
			for j in range(datai_mask.shape[1]):
				if i<box_coords[0] or i>box_coords[2] or j<box_coords[1] or j>box_coords[3]:
					datai_mask[i,j]=np.nan
		tb_nonpos=np.where(datatb[:,:,0,0]<=0)
		datatb[:,:,0,0][tb_nonpos]=np.nan
		datatblog=np.log10(datatb[:,:,0,0])
		pos=np.where(np.isnan(datatblog)==True)
		datatblog[pos]=0
		postb=np.where(datatblog>=7)
		postb1=np.where(datatblog<6)
		postb2=np.where((datatblog>=6) & (datatblog<7))
		datatblog[postb]=np.nan
		datatblog[postb1]=np.nan
		datatblog[posi]=np.nan
		datai[posi]=np.nan
		dataq[posi]=np.nan
		datai[postb1]=np.nan
		dataq[postb1]=np.nan
		area=np.nansum(np.isnan(datai)==False)
		datai[posq]=np.nan
		dataq[posq]=np.nan
		datai[postb]=np.nan
		dataq[postb]=np.nan	
		unmasked_pixels=np.nansum(np.isnan(dataq)==False)
		x=[]
		y=[]
		z=[]
		q_by_i=dataq/datai
		for k in range(datai.shape[0]):
			for l in range(datai.shape[1]):
				if np.isnan(q_by_i[k,l])==False:
					x.append(k)
					y.append(l)
					z.append(q_by_i[k,l])
		q_stack=np.vstack((x,y,z)).T
		del x,y,z
		pre_q_image=np.nanmean(np.abs(dataq_copy))
		if (unmasked_pixels/area)>0.3:
			AQ = np.c_[np.ones(q_stack.shape[0]), q_stack[:,:2], np.prod(q_stack[:,:2], axis=1), q_stack[:,:2]**2]
			CQ,_,_,_ = scipy.linalg.lstsq(AQ, q_stack[:,2])	
			for k in range(datai_copy.shape[0]):
				for l in range(datai_copy.shape[1]):
					if np.isnan(datai_mask[k,l])==False:
						dataq_copy[k,l] -= (CQ[4]*k**2. + CQ[5]*l**2. + CQ[3]*k*l + CQ[1]*k + CQ[2]*l + CQ[0])*datai_copy[k,l]
		else:
			CQ=[0,0,0,0,0,0]
		new_q_image=np.nanmean(np.abs(dataq_copy))
		datai=copy.deepcopy(datai_copy)
		datau=copy.deepcopy(datau_copy)
		datai[posi]=np.nan
		datau[posi]=np.nan
		datai[postb1]=np.nan
		datau[postb1]=np.nan
		area=np.nansum(np.isnan(datai)==False)
		datai[posu]=np.nan
		datau[posu]=np.nan
		datai[postb]=np.nan
		datau[postb]=np.nan
		unmasked_pixels=np.nansum(np.isnan(datau)==False)
		u_by_i=datau/datai
		x=[]
		y=[]
		z=[]
		for k in range(datai.shape[0]):
			for l in range(datai.shape[1]):
				if np.isnan(u_by_i[k,l])==False:
					x.append(k)
					y.append(l)
					z.append(u_by_i[k,l])
		u_stack=np.vstack((x,y,z)).T
		del x,y,z
		pre_u_image=np.nanmean(np.abs(datau_copy))
		if (unmasked_pixels/area)>0.3:
			AU = np.c_[np.ones(u_stack.shape[0]), u_stack[:,:2], np.prod(u_stack[:,:2], axis=1), u_stack[:,:2]**2]
			CU,_,_,_ = scipy.linalg.lstsq(AU, u_stack[:,2])
			for k in range(datai_copy.shape[0]):
				for l in range(datai_copy.shape[1]):
					if np.isnan(datai_mask[k,l])==False:
						datau_copy[k,l] -= (CU[4]*k**2. + CU[5]*l**2. + CU[3]*k*l + CU[1]*k + CU[2]*l + CU[0])*datai_copy[k,l]
		else:
			CU=[0,0,0,0,0,0]
		new_u_image=np.nanmean(np.abs(datau_copy))
		'''
		X,Y = np.meshgrid(np.arange(630, 650, 1), np.arange(630, 650, 1))
		XX = X.flatten()
		YY = Y.flatten()
		fig = plt.figure()
		ax = plt.axes(projection='3d')
		Z=(CU[4]*X**2. + CU[5]*Y**2. + CU[3]*X*Y + CU[1]*X + CU[2]*Y + CU[0])
		ax.plot_surface(X, Y, Z*100, rstride=1, cstride=1, alpha=0.3)
		ax.scatter(u_stack[:,0], u_stack[:,1], u_stack[:,2]*100, c='r', s=10)
		plt.xlabel('RA',fontsize=10)
		plt.ylabel('DEC',fontsize=10)
		ax.set_zlabel('Stokes U (%)',fontsize=10)
		ax.axis('tight')
		plt.show()
		X,Y = np.meshgrid(np.arange(625, 660, 1), np.arange(625, 660, 1))
		XX = X.flatten()
		YY = Y.flatten()
		Z=(CQ[4]*X**2. + CQ[5]*Y**2. + CQ[3]*X*Y + CQ[1]*X + CQ[2]*Y + CQ[0])
		fig = plt.figure()
		ax = plt.axes(projection='3d')
		ax.plot_surface(X, Y, Z*100, rstride=1, cstride=1, alpha=0.3)
		ax.scatter(q_stack[:,0], q_stack[:,1], q_stack[:,2]*100, c='r', s=10)
		plt.xlabel('RA',fontsize=10)
		plt.ylabel('DEC',fontsize=10)
		ax.set_zlabel('Stokes Q (%)',fontsize=10)
		ax.axis('tight')
		plt.show()
		'''
		posq=np.where(np.abs(dataq_copy)<(sigma*rmsq))
		posu=np.where(np.abs(datau_copy)<(sigma*rmsu))
		data[:,:,0,0]=datai_copy
		data[:,:,1,0]=dataq_copy
		data[:,:,2,0]=datau_copy
		ia.open(imagename)
		ia.putchunk(data)
		ia.close()
		for k in range(model_datai.shape[0]):
			for l in range(model_datai.shape[1]):
				if np.isnan(datai_mask[k,l])==False:
					model_dataq[k,l] -= (CU[4]*k**2. + CU[5]*l**2. + CU[3]*k*l + CU[1]*k + CU[2]*l + CU[0])*model_datai[k,l]
		for k in range(model_datai.shape[0]):
			for l in range(model_datai.shape[1]):
				if np.isnan(datai_mask[k,l])==False:
					model_datau[k,l] -= (CU[4]*k**2. + CU[5]*l**2. + CU[3]*k*l + CU[1]*k + CU[2]*l + CU[0])*model_datai[k,l]
		modelq=copy.deepcopy(model_dataq)
		modelu=copy.deepcopy(model_datau)
		modeldata[:,:,0,0]=model_datai
		modeldata[:,:,1,0]=modelq
		modeldata[:,:,2,0]=modelu
		dataq_by_i=dataq_copy/datai_copy
		datau_by_i=datau_copy/datai_copy
		dataq_by_i[posi]=np.nan
		dataq_by_i[posq]=np.nan
		dataq_by_i[postb]=np.nan
		datau_by_i[posi]=np.nan
		datau_by_i[posu]=np.nan
		datau_by_i[postb]=np.nan
		posqbyi=np.where(np.abs(dataq_by_i)<(sigma*rmsq/np.nanmean(datai)))
		posubyi=np.where(np.abs(datau_by_i)<(sigma*rmsu/np.nanmean(datai)))
		model_dataq[postb1]=0
		model_dataq[posq]=0
		model_dataq[posqbyi]=0
		model_datau[postb1]=0
		model_datau[posu]=0
		model_datau[posubyi]=0
		model_dataq[postb2]=0
		model_datau[postb2]=0
		modeldata[:,:,1,0]=model_dataq
		modeldata[:,:,2,0]=model_datau
		q_change=abs(new_q_image-pre_q_image)/pre_q_image
		u_change=abs(new_u_image-pre_u_image)/pre_u_image
		ia.open(modelname)
		ia.putchunk(modeldata)
		ia.close()
		if overwrite==True:
			if os.path.isdir(imagename.split('.image')[0]+'_quvcor_surface.image')==True:
				os.system('rm -rf '+imagename.split('.image')[0]+'_quvcor_surface.image')
			if os.path.isdir(modelname.split('.model')[0]+'_quvcor_surface.model')==True:
				os.system('rm -rf '+modelname.split('.model')[0]+'_quvcor_surface.model')
			os.system('cp -r '+imagename+' '+imagename.split('.image')[0]+'_quvcor_surface.image')
			os.system('cp -r '+modelname+' '+modelname.split('.model')[0]+'_quvcor_surface.model')
		os.system('rm -rf casa*log I*.image')
		os.system('rm -rf '+outfile_path+'/I*')
		del modeldata,data,datai_mask,dataq,datai,datau,model_datai,model_dataq,model_datau,datai_copy,dataq_copy,datau_copy
		return imagename,modelname,q_change,u_change

	def subtract_stokesV_solar_leakage_surface(self,imagename,modelname,sigma=10,do_fluxcal=False,overwrite=False):
		'''
		Function to subtract quadratic leakage surface

		Parameters
		----------
		imagename : str 
			Name of the image
		modelname : str 
			Name of the model
		sigma : float 
			N-sigma value above which any emission is considered to be real
		do_fluxcal : bool
			Perform polynomial based flux calibration (See details Kansabanik et al. 2021, submitted to ApJ)
		overwrite : bool
			Overwrite the input image and model
		Returns
		-------
		str
			Leakage surface subtracted image 
		str
			Leakage surface subtracted model
		float
			Stokes V change fraction
		'''
		self.pollog_verbose.info('Correcting Stokes I to Stokes V leakage surface.\n')
		if os.path.exists(imagename.split('.image')[0]+'_vcor_surface.image'):
			os.system('rm -rf '+imagename.split('.image')[0]+'_vcor_surface.image')
		if os.path.exists(modelname.split('.model')[0]+'_vcor_surface.model'):
			os.system('rm -rf '+modelname.split('.model')[0]+'_vcor_surface.model')
		if overwrite==False:
			os.system('cp -r '+imagename+' '+imagename.split('.image')[0]+'_vcor_surface.image')
			os.system('cp -r '+modelname+' '+'vcor_surface.model')
			imagename=imagename.split('.image')[0]+'_vcor_surface.image'
			modelname=modelname.split('.model')[0]+'_vcor_surface.model'
		outfile_path=os.path.dirname(os.path.realpath(imagename))
		os.system('rm -rf '+outfile_path+'/I*')
		imsubimage(imagename=imagename,outfile=outfile_path+'/I.image',stokes='I',dropdeg=False)
		if do_fluxcal==True:
			fluxcal_image=self.mwa_solar_fluxcal(imagename=outfile_path+'/I.image',outfile=outfile_path+'/I_fluxcal.image')
		else:
			fluxcal_image=outfile_path+'/I.image'
		major=imhead(imagename=fluxcal_image)['restoringbeam']['major']['value'] # In arcsec
		minor=imhead(imagename=fluxcal_image)['restoringbeam']['minor']['value'] # In arcsec
		freq=imhead(imagename=fluxcal_image)['refval'][-1]/10**9 # In GHz
		ref_beam_axes_multi=500000
		beam_axes_multi=major*minor	
		scale_tb_limit=(beam_axes_multi/ref_beam_axes_multi)
		expr='1.222e6*IM0/'+str(freq)+'^2/('+str(major*minor)+')'
		immath(imagename=fluxcal_image,outfile=outfile_path+'/I_Tb.image',mode='evalexpr',expr=expr)
		imhead(imagename=outfile_path+'/I_Tb.image', mode='put', hdkey='bunit', hdvalue='K')
		ia=image()
		ia.open(imagename)
		data=ia.getchunk()
		ia.close()
		a=image()
		ia.open(modelname)
		modeldata=ia.getchunk()
		ia.close()
		model_datai=modeldata[:,:,0,0]
		model_datav=modeldata[:,:,3,0]
		rmsv=imstat(imagename=imagename,box=self.rms_box,stokes='V')['rms'][0]
		rmsi=imstat(imagename=imagename,box=self.rms_box,stokes='I')['rms'][0]
		i_flux=imstat(imagename=imagename,stokes='I')['flux'][0]
		datai=data[:,:,0,0]
		datav=data[:,:,3,0]
		datai_copy=copy.deepcopy(datai)
		datav_copy=copy.deepcopy(datav)
		datav1=copy.deepcopy(datav)
		posi=np.where(datai<(sigma*rmsi))
		maxpos=imstat(imagename=imagename,stokes='I')['maxpos']
		box=self.negative_box(maxpos,box_width=2).split(',')
		box_coords=[int(i) for i in box]
		posv=np.where(np.abs(datav)<(sigma*rmsv))
		ia.open(outfile_path+'/I_Tb.image')
		datatb=ia.getchunk()
		ia.close()
		datai_mask=copy.deepcopy(datai)
		for i in range(datai_mask.shape[0]):
			for j in range(datai_mask.shape[1]):
				if i<box_coords[0] or i>box_coords[2] or j<box_coords[1] or j>box_coords[3]:
					datai_mask[i,j]=np.nan
		tb_nonpos=np.where(datatb[:,:,0,0]<=0)
		datatb[:,:,0,0][tb_nonpos]=np.nan
		datatblog=np.log10(datatb[:,:,0,0])
		pos=np.where(np.isnan(datatblog)==True)
		datatblog[pos]=0
		postb=np.where(datatblog>=7)
		postb1=np.where(datatblog<6)
		datai[posi]=np.nan
		datav[posi]=np.nan
		datai[postb1]=np.nan
		datav[postb1]=np.nan
		area=np.nansum(np.isnan(datai)==False)
		datai[posv]=np.nan
		datav[posv]=np.nan
		datai[postb]=np.nan
		datav[postb]=np.nan	
		unmasked_pixels=np.nansum(np.isnan(datav)==False)
		v_by_i=datav/datai
		v_by_i[abs(v_by_i)<0.05]=np.nan
		x=[]
		y=[]
		z=[]
		for k in range(datai.shape[0]):
			for l in range(datai.shape[1]):
				if np.isnan(v_by_i[k,l])==False:
					x.append(k)
					y.append(l)
					z.append(v_by_i[k,l])
		v_stack=np.vstack((x,y,z)).T
		del x,y,z
		pre_v_image=np.nansum(np.abs(datav_copy))
		if (unmasked_pixels/area)>0.3:
			AV = np.c_[np.ones(v_stack.shape[0]), v_stack[:,:2], np.prod(v_stack[:,:2], axis=1), v_stack[:,:2]**2]
			CV,_,_,_ = scipy.linalg.lstsq(AV, v_stack[:,2])	
			for k in range(datai_copy.shape[0]):
				for l in range(datai_copy.shape[1]):
					if np.isnan(datai_mask[k,l])==False:
						datav_copy[k,l] -= (CV[4]*k**2. + CV[5]*l**2. + CV[3]*k*l + CV[1]*k + CV[2]*l + CV[0])*datai_copy[k,l]
		else:
			CV=[0,0,0,0,0,0]
		new_v_image=np.nansum(np.abs(datav_copy))
		'''
		X,Y = np.meshgrid(np.arange(630, 650, 1), np.arange(630, 650, 1))
		XX = X.flatten()
		YY = Y.flatten()
		fig = plt.figure()
		ax = plt.axes(projection='3d')
		Z=(CV[4]*X**2. + CV[5]*Y**2. + CV[3]*X*Y + CV[1]*X + CV[2]*Y + CV[0])
		ax.plot_surface(X, Y, Z*100, rstride=1, cstride=1, alpha=0.3)
		ax.scatter(v_stack[:,0], v_stack[:,1], v_stack[:,2]*100, c='r', s=10)
		plt.xlabel('RA',fontsize=10)
		plt.ylabel('DEC',fontsize=10)
		ax.set_zlabel('Stokes V (%)',fontsize=10)
		ax.axis('tight')
		plt.show()
		'''
		posv=np.where(np.abs(datav_copy)<(sigma*rmsv))
		data[:,:,0,0]=datai_copy
		data[:,:,3,0]=datav_copy
		ia.open(imagename)
		ia.putchunk(data)
		ia.close()
		for k in range(model_datai.shape[0]):
			for l in range(model_datai.shape[1]):
				if np.isnan(datai_mask[k,l])==False:
					model_datav[k,l] -= (CV[4]*k**2. + CV[5]*l**2. + CV[3]*k*l + CV[1]*k + CV[2]*l + CV[0])*model_datai[k,l]
		modelv=copy.deepcopy(model_datav)
		v_change=abs(pre_v_image-new_v_image)/pre_v_image
		modeldata[:,:,0,0]=model_datai
		modeldata[:,:,3,0]=modelv
		datav_by_i=datav_copy/datai_copy
		datav_by_i[posi]=np.nan
		datav_by_i[posv]=np.nan
		model_datav[posv]=0
		modeldata[:,:,3,0]=model_datav
		ia.open(modelname)
		ia.putchunk(modeldata)
		ia.close()
		if overwrite==True:
			if os.path.isdir(imagename.split('.image')[0]+'_vcor_surface.image')==True:
				os.system('rm -rf '+imagename.split('.image')[0]+'_vcor_surface.image')
			if os.path.isdir(modelname.split('.model')[0]+'_vcor_surface.model')==True:
				os.system('rm -rf '+modelname.split('.model')[0]+'_vcor_surface.model')
			os.system('cp -r '+imagename+' '+imagename.split('.image')[0]+'_vcor_surface.image')
			os.system('cp -r '+modelname+' '+modelname.split('.model')[0]+'_vcor_surface.model')
		os.system('rm -rf casa*log I*.image')
		os.system('rm -rf '+outfile_path+'/I*')
		del modeldata,data,datai_mask,datav,datai,model_datai,model_datav,datai_copy,datav_copy
		return imagename,modelname,v_change	

	def mwa_solar_fluxcal(self,imagename,outfile,atten=10): #TODO : use reference band shape
		'''
		Function to flux calibrate MWA solar observations using method described in Kansabanik et al. 2021

		Parameters
		----------
		imagename : str 
			Name of the image
		outfile : str
			Output file name
		atten : float 
			Attenuator gain value in dB
		Returns
		-------
		str
			Flux calibrated image
		''' 
		freq=imhead(imagename)['refval'][-1]/10**6 # MHz
		fluxscale_poly=np.poly1d(np.load(datadir+'/flux_scale_polyfit.npy',allow_pickle=True)[0])
		fluxscale=fluxscale_poly(freq)
		if atten!=10: # Valid for 10 and 14 dB, not tested for any other values and typically other values not used for solar observations
			fluxscale=fluxscale*10**((atten-10)/20.0)
		immath(imagename=imagename,outfile=outfile,mode='evalexpr',expr='IM0*'+str(fluxscale))
		return outfile

	def cal_solar_qu_leakage(self,imagename,sigma=10,do_fluxcal=False): 
		'''
		Function to calculate Stokes QU leakage for solar observation (Not vaild for any other astrophysical observation)

		Parameters
		----------
		imagename : str 
			Name of the image
		do_fluxcal : bool 
			Perform flux calibration or not
		outfile_path : str 
			Name of the directory to save leakage data numpy table (default : image directory)
		Returns
		-------
		float
			Stokes Q leakage
		float 
			Stokes U leakage
		(Two numpy table will also be saved) 
		'''
		outfile_path=os.path.dirname(os.path.realpath(imagename))
		os.system('rm -rf '+outfile_path+'/I*')
		imsubimage(imagename=imagename,outfile=outfile_path+'/I.image',stokes='I',dropdeg=False)
		if do_fluxcal==True:
			fluxcal_image=self.mwa_solar_fluxcal(imagename=outfile_path+'/I.image',outfile=outfile_path+'/I_fluxcal.image')
		else:
			fluxcal_image=outfile_path+'/I.image'
		major=imhead(imagename=fluxcal_image)['restoringbeam']['major']['value'] # In arcsec
		minor=imhead(imagename=fluxcal_image)['restoringbeam']['minor']['value'] # In arcsec
		freq=imhead(imagename=fluxcal_image)['refval'][-1]/10**9 # In GHz
		ref_beam_axes_multi=500000
		beam_axes_multi=major*minor	
		scale_tb_limit=(beam_axes_multi/ref_beam_axes_multi)
		expr='1.222e6*IM0/'+str(freq)+'^2/('+str(major*minor)+')'
		immath(imagename=fluxcal_image,outfile=outfile_path+'/I_Tb.image',mode='evalexpr',expr=expr)
		imhead(imagename=outfile_path+'/I_Tb.image', mode='put', hdkey='bunit', hdvalue='K')
		ia=image()
		ia.open(imagename)
		data=ia.getchunk()
		ia.close()
		rmsq=imstat(imagename=imagename,box=self.rms_box,stokes='Q')['rms'][0]
		rmsu=imstat(imagename=imagename,box=self.rms_box,stokes='U')['rms'][0]
		rmsi=imstat(imagename=imagename,box=self.rms_box,stokes='I')['rms'][0]
		i_flux=imstat(imagename=imagename,stokes='I')['flux'][0]
		dataq=data[:,:,1,0]
		datau=data[:,:,2,0]
		datai=data[:,:,0,0]
		posi=np.where(datai<(sigma*rmsi))
		posq=np.where(dataq<(5*rmsq))
		posu=np.where(datau<(5*rmsu))
		datai[posi]=np.nan
		dataq[posi]=np.nan
		datau[posi]=np.nan
		datai[posq]=np.nan
		dataq[posq]=np.nan
		datai[posu]=np.nan
		datau[posu]=np.nan
		ia.open(outfile_path+'/I_Tb.image')
		datatb=ia.getchunk()
		ia.close()
		tb_nonpos=np.where(datatb[:,:,0,0]<=0)
		datatb[:,:,0,0][tb_nonpos]=np.nan
		datatblog=np.log10(datatb[:,:,0,0])
		pos=np.where(np.isnan(datatblog)==True)
		datatblog[pos]=0
		postb=np.where(datatblog>=7)
		postb1=np.where(datatblog<6)
		datai[postb]=np.nan
		dataq[postb]=np.nan
		datau[postb]=np.nan
		datai[postb1]=np.nan
		dataq[postb1]=np.nan
		datau[postb1]=np.nan
		i_flux=np.nanmean(datai)
		leakage_list=[]		
		if np.sum(np.isnan(i_flux))==len(i_flux.flatten()):
			leakage_list=[0,0]
			return leakage_list		
		for stokes in ['Q','U']:
			if stokes=='Q':
				datal=dataq
				rmsl=rmsq
			elif stokes=='U':
				datal=datau
				rmsl=rmsu
			step_range=[0.1,0.01,0.001]
			start_range=-1
			end_range=1
			x=[]
			y=[]
			for step in step_range:
				l_range=np.arange(start_range,end_range,step)
				results=[]
				for l in l_range:
					x.append(l)
					r=self.solarlin_pol_minimise(datai,datal,l,rmsl,i_flux)
					y.append(r)
					results.append(r)
				results=np.array(results)
				minval=l_range[np.argmin(results)]
				start_range=(start_range+minval)/3.0
				end_range=(end_range+minval)/3.0
			y=np.array(y)
			leakage=x[np.argmin(y)]
			np.save(outfile_path+'/'+os.path.basename(imagename)+'_'+stokes+'_leakage',np.array([x,y,leakage],dtype=object))
			leakage_list.append(leakage)
		os.system('rm -rf casa*log I*.image')
		os.system('rm -rf '+outfile_path+'/I*')
		self.pollog_verbose.info('Calculation of Stokes leakages has been done.\n')
		return leakage_list[0],leakage_list[1]

	def create_circular_mask(self,h, w, center=None, radius=None):
		'''
		Function to create a circular mask

		Parameters
		----------
		h : int 
			Number of pixels along Y axis
		w : int 
			Number of pixels along X axis
		center : tuple 
			(x_cen,y_cen), center of the mask circle
		radius : float 
			Radius in number of pixels
		Returns
		-------
		numpy.array
			Cirular mask
		'''
		if center is None:
		    center=(int(w/2),int(h/2))
		if radius is None:
		    radius=min(center[0],center[1],w-center[0],h-center[1])
		Y,X=np.ogrid[:h,:w]
		dist_from_center=np.sqrt((X-center[0])**2+(Y-center[1])**2)
		mask=dist_from_center<=radius
		return mask

	def cal_solar_v_leakage(self,imagename,sigma=10,do_fluxcal=False): 
		'''
		Function to calculate Stokes V leakage for solar observation (Not vaild for any other astrophysical observation)

		Parameters
		----------
		imagename : str 
			Name of the image
		sigma : float
			Sigma value for thresholding
		do_fluxcal : bool 
			Perform flux calibration or not
		Returns
		-------
		float
			Stokes V leakage (A numpy table will also be saved) 
		'''
		outfile_path=os.path.dirname(os.path.realpath(imagename))
		os.system('rm -rf '+outfile_path+'/I*')
		self.pollog_verbose.info('Correcting Stokes V leakage.\n')
		imsubimage(imagename=imagename,outfile=outfile_path+'/I.image',stokes='I',dropdeg=False)
		if do_fluxcal==True:
			fluxcal_image=self.mwa_solar_fluxcal(imagename=outfile_path+'/I.image',outfile=outfile_path+'/I_fluxcal.image')
		else:
			fluxcal_image=outfile_path+'/I.image'
		major=imhead(imagename=fluxcal_image)['restoringbeam']['major']['value'] # In arcsec
		minor=imhead(imagename=fluxcal_image)['restoringbeam']['minor']['value'] # In arcsec
		freq=imhead(imagename=fluxcal_image)['refval'][-1]/10**9 # In GHz
		ref_beam_axes_multi=500000
		beam_axes_multi=major*minor	
		scale_tb_limit=(beam_axes_multi/ref_beam_axes_multi)
		expr='1.222e6*IM0/'+str(freq)+'^2/('+str(major*minor)+')'
		immath(imagename=fluxcal_image,outfile=outfile_path+'/I_Tb.image',mode='evalexpr',expr=expr)
		imhead(imagename=outfile_path+'/I_Tb.image', mode='put', hdkey='bunit', hdvalue='K')
		ia=image()
		ia.open(imagename)
		data=ia.getchunk()
		ia.close()
		rmsv=imstat(imagename=imagename,box=self.rms_box,stokes='V')['rms'][0]
		rmsi=imstat(imagename=imagename,box=self.rms_box,stokes='I')['rms'][0]
		datav=data[:,:,3,0]
		datai=data[:,:,0,0]
		datav_copy=copy.deepcopy(datav)
		posi=np.where(datai<(sigma*rmsi))
		posv=np.where(np.abs(datav)<(sigma*rmsv))
		datai[posi]=np.nan
		datav[posi]=np.nan
		datai[posv]=np.nan
		datav[posv]=np.nan	
		area=np.nansum(np.isnan(datav)==False)
		ia.open(outfile_path+'/I_Tb.image')
		datatb=ia.getchunk()
		ia.close()
		tb_nonpos=np.where(datatb[:,:,0,0]<=0)
		datatb[:,:,0,0][tb_nonpos]=np.nan
		datatblog=np.log10(datatb[:,:,0,0])
		pos=np.where(np.isnan(datatblog)==True)
		datatblog[posi]=np.nan
		datatblog[posv]=np.nan
		datatblog[pos]=np.nan
		postb=np.where(datatblog>=7)
		postb1=np.where(datatblog<6)
		datai[postb]=np.nan
		datav[postb]=np.nan
		datai[postb1]=np.nan
		datav[postb1]=np.nan
		v_by_i=datav/datai
		mean_i_flux=np.nanmean(datai)
		unmasked_pixels=np.nansum(np.isnan(datav)==False)
		if (unmasked_pixels/area)<0.1:
			np.save(outfile_path+'/'+os.path.basename(imagename)+'_V_leakage',np.array([0,0,0],dtype=object))
			os.system('rm -rf casa*log I*.image')
			os.system('rm -rf '+outfile_path+'/I_*')
			return 0,0
		step_range=[0.1,0.01,0.001]
		start_range=-1
		end_range=1
		x=[]
		y=[]
		for step in step_range:
			l_range=np.arange(start_range,end_range,step)
			results=[]
			for l in l_range:
				x.append(l)
				r=self.solarcir_pol_minimise(datai,datav,l,rmsv,mean_i_flux,sigma)
				y.append(r)
				results.append(r)
			results=np.array(results)
			minval=l_range[np.argmin(results)]
			minvals=l_range[np.where(l_range==minval)]
			minval=np.nanmedian(minvals)
			start_range=(start_range+minval)/3.0
			end_range=(end_range+minval)/3.0
		y=np.array(y)
		x=np.array(x)
		ymin=np.nanmin(y)
		pos=np.where(y==ymin)
		x_min=x[pos]
		leakage=np.nanmedian(x_min)
		np.save(outfile_path+'/'+os.path.basename(imagename)+'_V_leakage',np.array([x,y,leakage],dtype=object))
		os.system('rm -rf casa*log I*.image')
		os.system('rm -rf '+outfile_path+'/I*')
		return leakage

	def correct_solar_quv_leakage(self,imagename,modelname,sigma,overwrite=False,stokes='QUV'):
		'''
		Function to correction for solar Stokes Q, U and V leakage 
		(It is based on the fact that we do not expect any linear polarisation from the Quiet Sun emission)

		Parameters
		----------
		imagename : str 
			Name of image
		modelname : str 
			Name of the model
		sigma : float 
			N-sigma threshold to choose Stokes I emission region
		overwrite : bool 
			Overwrite the image and model or not
		stokes : str
			Stokes plane to correct
		Returns
		-------
		str
			Stokes Q ,U, V leakage corrected image name
		str 
			Stokes Q ,U, V leakage corrected model name 
		float
			Stokes Q fractional change
		float 
			Stokes U fractional change
		float
			Stokes V fractional change	 
		'''
		q_change_frac,u_change_frac,v_change_frac=0,0,0
		self.pollog_verbose.info('Correcting Stokes I to Stokes '+','.join(list(stokes))+'.....\n')
		self.pollog_verbose.info('Image name : '+imagename+'\n')
		self.pollog_verbose.info('Model name : '+modelname+'\n')
		if os.path.exists(imagename.split('.image')[0]+'_quvcor_surface.image'):
			os.system('rm -rf '+imagename.split('.image')[0]+'_quvcor_surface.image')
		if os.path.exists(modelname.split('.model')[0]+'_quvcor_surface.model'):
			os.system('rm -rf '+modelname.split('.model')[0]+'_quvcor_surface.model')
		if overwrite==False:
			os.system('cp -r '+imagename+' '+imagename.split('.image')[0]+'_quvcor_surface.image')
			os.system('cp -r '+modelname+' '+modelname.split('.model')[0]+'_quvcor_surface.model')
			imagename=imagename.split('.image')[0]+'_quvcor_surface.image'
			modelname=modelname.split('.model')[0]+'_quvcor_surface.model'
		outfile_path=os.path.dirname(os.path.realpath(imagename))
		imagename,modelname,q_change_frac,u_change_frac=self.subtract_leakage_surface(imagename,modelname,sigma=sigma,do_fluxcal=True,overwrite=True)
		if 'V'in stokes:
			imagename,modelname,v_change_frac=self.subtract_stokesV_solar_leakage_surface(imagename,modelname,sigma=sigma,do_fluxcal=True,overwrite=True)
		ia=image()
		ia.open(imagename) # Correcting image
		data=ia.getchunk()
		datai=data[:,:,0,:]
		dataq=data[:,:,1,:]
		datau=data[:,:,2,:]
		datav=data[:,:,3,:]
		rmsq=imstat(imagename=imagename,box=self.rms_box,stokes='Q')['rms'][0]
		rmsu=imstat(imagename=imagename,box=self.rms_box,stokes='U')['rms'][0]
		rmsv=imstat(imagename=imagename,box=self.rms_box,stokes='V')['rms'][0]
		rmsi=imstat(imagename=imagename,box=self.rms_box,stokes='I')['rms'][0]
		ia.putchunk(data)
		ia.close()
		imsubimage(imagename=imagename,outfile=outfile_path+'/I.image',stokes='I',dropdeg=False)
		fluxcal_image=self.mwa_solar_fluxcal(imagename=outfile_path+'/I.image',outfile=outfile_path+'/I_fluxcal.image')
		major=imhead(imagename=fluxcal_image)['restoringbeam']['major']['value'] # In arcsec
		minor=imhead(imagename=fluxcal_image)['restoringbeam']['minor']['value'] # In arcsec
		freq=imhead(imagename=fluxcal_image)['refval'][-1]/10**9 # In GHz
		ref_beam_axes_multi=500000
		beam_axes_multi=major*minor	
		scale_tb_limit=(beam_axes_multi/ref_beam_axes_multi)
		expr='1.222e6*IM0/'+str(freq)+'^2/('+str(major*minor)+')'
		immath(imagename=fluxcal_image,outfile=outfile_path+'/I_Tb.image',mode='evalexpr',expr=expr)
		imhead(imagename=outfile_path+'/I_Tb.image', mode='put', hdkey='bunit', hdvalue='K')
		ia.open(outfile_path+'/I_Tb.image')
		datatb=ia.getchunk()
		ia.close()
		tb_nonpos=np.where(datatb[:,:,0,0]<=0)
		datatb[:,:,0,0][tb_nonpos]=np.nan
		datatblog=np.log10(datatb[:,:,0,0])
		pos=np.where(np.isnan(datatblog)==True)
		datatblog[pos]=0
		postb=np.where(datatblog>=7)
		postb1=np.where(datatblog<6)
		postb2=np.where((datatblog>=6) & (datatblog<7))
		ia.open(modelname) # Correcting model
		datam=ia.getchunk()
		im=datam[:,:,0,:]
		qm=datam[:,:,1,:]
		um=datam[:,:,2,:]
		vm=datam[:,:,3,:]
		posqm=np.where(qm==0)
		posum=np.where(um==0)
		posvm=np.where(vm==0)
		posi=np.where(datai<(sigma*rmsi))
		posq=np.where(np.abs(dataq)<(1.5*sigma*rmsq))
		posu=np.where(np.abs(datau)<(1.5*sigma*rmsu))
		posv=np.where(np.abs(datav)<(1.5*sigma*rmsv))
		qm[posi]=0
		qm[posq]=0
		um[posi]=0
		um[posu]=0
		vm[posi]=0
		vm[posv]=0
		qm[postb2]=0
		um[postb2]=0
		qm[posqm]=0
		um[posum]=0
		vm[posvm]=0
		if 'Q' in stokes:
			datam[:,:,1,:]=qm
		if 'U' in stokes:
			datam[:,:,2,:]=um
		if 'V' in stokes:
			datam[:,:,3,:]=vm
		ia.putchunk(datam)
		ia.close()
		if overwrite==True:
			if os.path.isdir(imagename.split('.image')[0]+'_quvcor_surface.image')==True:
				os.system('cp -r '+imagename.split('.image')[0]+'_quvcor_surface.image '+imagename)
				os.system('rm -rf '+imagename.split('.image')[0]+'_quvcor_surface.image')
			if os.path.isdir(modelname.split('.model')[0]+'_quvcor_surface.model')==True:
				os.system('cp -r '+modelname.split('.model')[0]+'_quvcor_surface.model '+modelname)
				os.system('rm -rf '+modelname.split('.model')[0]+'_quvcor_surface.model')
		os.system('rm -rf casa*log '+modelname.split('.model')[0]+'_quvcor_surface.model '+imagename.split('.image')[0]+'_quvcor_surface.image'+\
				modelname.split('.model')[0]+'_vcor_surface.model '+imagename.split('.image')[0]+'_vcor_surface.image')
		os.system('rm -rf '+outfile_path+'/I*')
		return imagename,modelname,q_change_frac,u_change_frac,v_change_frac

	def pol_model_threshold(self,imagename,modelname,sigma,polmodel_thresh):
		'''
		Function to put rms based threshold on polarisation model

		Parameters
		----------
		imagename : str 
			Name of the image
		modelname : str 
			Name of the model
		sigma : float 
			Sigma value for threshold 
		polmodel_thresh : float 
			Multiplying factor to sigma for polarisation model threshold
		Returns
		-------
		str
			Imagename 
		str
			Modelname
		'''
		self.pollog_verbose.info('Threshold polarisation models with sigma = '+str(sigma)+'\n')
		rmsi=imstat(imagename=imagename,box=self.rms_box,stokes='I')['rms'][0]
		rmsq=imstat(imagename=imagename,box=self.rms_box,stokes='Q')['rms'][0]
		rmsu=imstat(imagename=imagename,box=self.rms_box,stokes='U')['rms'][0]
		rmsv=imstat(imagename=imagename,box=self.rms_box,stokes='V')['rms'][0]
		ia=image()
		ia.open(imagename)
		data=ia.getchunk()
		ia.close()
		ia.open(modelname)
		modeldata=ia.getchunk()
		ia.close()
		nonI_pos=np.where(modeldata[:,:,0,:]<=0)
		posi=np.where(data[:,:,0,:]<=sigma*rmsi)
		posq=np.where(np.abs(data[:,:,1,:])<=(polmodel_thresh*sigma*rmsq))
		posu=np.where(np.abs(data[:,:,2,:])<=(polmodel_thresh*sigma*rmsu))
		posv=np.where(np.abs(data[:,:,3,:])<=(polmodel_thresh*sigma*rmsv))
		modeldata[:,:,0,:][nonI_pos]=0
		modeldata[:,:,1,:][nonI_pos]=0
		modeldata[:,:,2,:][nonI_pos]=0
		modeldata[:,:,3,:][nonI_pos]=0
		modeldata[:,:,0,:][posi]=0
		modeldata[:,:,1,:][posq]=0
		modeldata[:,:,2,:][posu]=0
		modeldata[:,:,3,:][posv]=0
		ia.open(modelname)
		ia.putchunk(modeldata)
		ia.close()
		return imagename,modelname

	def subtract_background_sources(self,stokes,rms_thresh,sigma,modelthres,maskregion='',includeregion=False,overwrite=False,modify_datacolumn=False,cpus=5,absmem=2,\
									weight='briggs',robust=0.5):
		'''
		Function to subtract background sources outside the mask region

		Parameters
		----------
		stokes : str 
			Stokes parameters to image
		rms_thresh : list 
			RMS list for each Stokes plane
		sigma : float 
			Sigma value for rms based thresholding
		modelthresh : float
			Sigma value to put a threshold on model
		maskregion : str 
			Region outside which the sources are subtracted.
		includeregion : bool 
			Include the maskregion for subtraction (True) or exclude the mask region for subtraction (False).
		overwrite : bool 
			Overwrite the corrected datacolumn with the subtracted visibilities.
		modelify_datacolumn : bool 
			Modify the datacolumn
		weight : str
			Visibility weighting during imaging , 'uniform', 'natural' or 'briggs'. Default is 'briggs'
		robust : float
			Robust parameter for briggs weighting, default : 0.5
		Returns
		-------
		str
			Backgroud source subtracted ms
		bool 
			Background source subtracted or not
		'''
		imagename=self.msname.split('.ms')[0]+'_background_sources' # Imagename prefix
		threshold=[str(rms*sigma)+'Jy' for rms in rms_thresh]
		os.system('rm -rf '+imagename+'*')

		self.pollog_verbose.info('poltclean(vis=\''+self.msname+'\',imagename=\''+imagename+'.maskregion\',selectdata=True,stokes=\''+str(stokes)+'\',imsize=[\''+\
		str(self.imsize)+'\'],cell=\''+str(self.cellsize)+'arcsec\',niter=0,gain=0.1,threshold=\''+str(threshold)+'\',deconvolver=\'multiscale\',scales='+\
		str(self.multiscale_scales)+',uvtaper=\''+self.uvtaper+'\',weighting=\''+weight+'\',robust='+str(robust*2)+'interactive=False,usemask=\'user\',mask=\''+maskregion+'\')\n')	
		poltclean(vis=self.msname,imagename=imagename+'.maskregion',selectdata=True,stokes=stokes,imsize=[self.imsize],cell=str(self.cellsize)+'arcsec',niter=0,\
		gain=0.1,threshold=threshold,deconvolver='multiscale',scales=self.multiscale_scales,uvtaper=self.uvtaper,weighting=weight,robust=robust*2,interactive=False,usemask='user',\
		mask=maskregion)
		if self.wsclean==False:
			robust=robust*2
			self.pollog_verbose.info('poltclean(vis=\''+self.msname+'\',imagename=\''+imagename+'\',selectdata=True,stokes=\''+str(stokes)+'\',imsize=['+str(self.imsize)+'],cell=\''+\
			str(self.cellsize)+'arcsec\',niter=10000,gain=0.1,threshold='+str(threshold)+',deconvolver=\'multiscale\',scales='+\
			str(self.multiscale_scales)+',uvtaper=\''+self.uvtaper+'\',weighting=\''+weight+'\',robust='+str(robust)+',interactive=False,usemask=\'user\')\n')	
			poltclean(vis=self.msname,imagename=imagename,selectdata=True,stokes=stokes,imsize=[self.imsize],cell=str(self.cellsize)+'arcsec',niter=10000,\
			gain=0.1,threshold=threshold,deconvolver='multiscale',scales=self.multiscale_scales,uvtaper=self.uvtaper,weighting=weight,robust=robust,interactive=False,usemask='user')
			casa_imagename=imagename+'.image'
			casa_modelname=imagename+'.model'
			casa_residualname=imagename+'.residual'
		else:
			scales=[str(i) for i in self.multiscale_scales]
			if weight=='briggs':
				weight=weight+' '+str(robust)
			wsclean_args=['-scale '+str(self.cellsize)+'asec','-size '+str(self.imsize)+' '+str(self.imsize),'-no-dirty','-j '+str(cpus),\
			'-abs-mem '+str(absmem),'-weight '+weight,'-taper-tukey 10','-name '+imagename,'-pol iquv','-nwlayers '+str(10),'-maxuv-l '+str(self.imaging_maxuv),\
			'-minuv-l '+str(self.imaging_minuv),'-niter 10000','-mgain 0.8','-auto-threshold '+str(sigma),'-auto-mask '+str(sigma+0.5),'-gain 0.1','-multiscale',\
			'-multiscale-scales '+','.join(scales)]
			wsclean_args.append('-quiet')
			if self.verbose:
				self.pollog_verbose.info('wsclean '+' '.join(wsclean_args)+' '+self.msname)			
			os.system('wsclean '+' '.join(wsclean_args)+' '+self.msname)			
			wsclean_images=glob.glob(imagename+'*image.fits')
			wsclean_models=glob.glob(imagename+'*model.fits')
			wsclean_residuals=glob.glob(imagename+'*residual.fits')
			casa_imagename=self.wsclean_to_casaimage(wsclean_images=wsclean_images,casaimage_prefix=imagename,imagetype='image',keep_wsclean_images=False)
			casa_modelname=self.wsclean_to_casaimage(wsclean_images=wsclean_models,casaimage_prefix=imagename,imagetype='model',keep_wsclean_images=False)
			casa_residualname=self.wsclean_to_casaimage(wsclean_images=wsclean_residuals,casaimage_prefix=imagename,imagetype='residual',keep_wsclean_images=False)
			os.system('rm -rf *psf*')

		ia=image()
		ia.open(imagename+'.maskregion.mask')
		maskregion=ia.getchunk()
		ia.close()
		if includeregion==False:
			pos_zero=np.where(maskregion==0)
			pos_one=np.where(maskregion==1)
			maskregion[pos_zero]=1
			maskregion[pos_one]=0
		ia.open(casa_modelname)
		modeldata=ia.getchunk()
		modeldata=modeldata*maskregion
		if np.nansum(modeldata[:,:,0,:])!=0:
			subtracted=True
		else:
			subtracted=False
		ia.putchunk(modeldata)
		ia.close()
		self.pol_model_threshold(casa_imagename,casa_modelname,sigma,1)
		if overwrite==False:
			if os.path.exists(self.msname+'.bsub'):
				os.system('rm -rf '+self.msname+'.bsub')
			os.system('cp -r '+self.msname+' '+self.msname+'.bsub')
			msname=self.msname+'.bsub'
		else:
			msname=self.msname
		self.pollog_verbose.info('delmod(vis=\''+msname+'\',scr=True,otf=True)\n')
		delmod(vis=msname,scr=True,otf=True)
		self.pollog_verbose.info('ft(vis=\''+msname+'\',model=\''+casa_modelname+'\',usescratch=True)\n')
		ft(vis=msname,model=casa_modelname,usescratch=True)
		self.pollog_verbose.info('uvsub(vis=\''+msname+'\',reverse=False)\n')
		uvsub(vis=msname,reverse=False)
		if modify_datacolumn==True:
			tb=table()
			tb.open(msname)
			cor_data=tb.getcol('CORRECTED_DATA')
			tb.close()
			tb.open(msname,nomodify=False)
			tb.putcol('DATA',cor_data)
			tb.flush()
			tb.close()
		os.system('rm -rf casa*log *background_sources*')
		return msname,subtracted

	def compare_leakages(self,present_image='',previous_image='',present_model='',previous_model='',outputfile_prefix='',overwrite=False,TB_limit=-1):
		'''
		Function to compare Stokes I to Stokes Q,U,V leakages with the previous image

		Parameters
		----------
		present_image : str 
			Name of the present CASA image
		previous_image : str 
			Name of the previous CASA image
		present_model : str 
			Name of the present CASA model
		previous_model : str 
			Name of the previous CASA model
		outputfile_prefix : str 
			Final output file prefix name
		overwrite : bool 
			Overwrite the present image or not
		TB_limit : float 
			Brightness temperature limit to calculate polarised flux density
		Returns
		-------
		str
			Output file name
		'''
		if os.path.isdir(present_image)==False or os.path.isdir(previous_image)==False or os.path.isdir(present_model)==False or os.path.isdir(previous_model)==False:
			if self.verbose==False:
				print ('Either present or previous image or model is not present.\n')
			self.pollog_verbose.info('Either present or previous image or model is not present.\n')
			return 
		self.pollog_verbose.info('Compairing '+present_image+' with '+previous_image+'\n')
		if present_image[-1]=='/':
			present_image=present_image[:-1]
		if previous_image[-1]=='/':
			previous_image=previous_image[:-1]
		try:
			outfile_path=os.path.basedir(present_image)
		except:
			outfile_path=os.getcwd()
		os.system('rm -rf casa*log '+outfile_path+'/I* '+outfile_path+'/'+os.path.basename(present_image)+'.temp '+outfile_path+'/'+os.path.basename(previous_image)+'.temp')
		ia=image()
		ia.open(present_image)
		present_data=ia.getchunk()
		ia.close()
		ia.open(previous_image)
		previous_data=ia.getchunk()
		ia.close()
		ia=image()
		ia.open(present_model)
		present_modeldata=ia.getchunk()
		ia.close()
		ia.open(previous_model)
		previous_modeldata=ia.getchunk()
		ia.close()
		present_data_copy=copy.deepcopy(present_data)
		previous_data_copy=copy.deepcopy(previous_data)
		present_model_copy=copy.deepcopy(present_modeldata)
		previous_model_copy=copy.deepcopy(previous_modeldata)
		os.system('cp -r '+present_image+' '+outfile_path+'/'+os.path.basename(present_image)+'.temp')
		os.system('cp -r '+previous_image+' '+outfile_path+'/'+os.path.basename(previous_image)+'.temp')
		if TB_limit!=-1:
			imsubimage(imagename=present_image,outfile=outfile_path+'/I.image',stokes='I',dropdeg=False)
			fluxcal_image=self.mwa_solar_fluxcal(imagename=outfile_path+'/I.image',outfile=outfile_path+'/I_fluxcal.image')
			major=imhead(imagename=fluxcal_image)['restoringbeam']['major']['value'] # In arcsec
			minor=imhead(imagename=fluxcal_image)['restoringbeam']['minor']['value'] # In arcsec
			freq=imhead(imagename=fluxcal_image)['refval'][-1]/10**9 # In GHz
			ref_beam_axes_multi=500000
			beam_axes_multi=major*minor	
			expr='1.222e6*IM0/'+str(freq)+'^2/('+str(major*minor)+')'
			immath(imagename=fluxcal_image,outfile=outfile_path+'/I_Tb.image',mode='evalexpr',expr=expr)
			imhead(imagename=outfile_path+'/I_Tb.image', mode='put', hdkey='bunit', hdvalue='K')
			ia.open(outfile_path+'/I_Tb.image')
			datatb=ia.getchunk()
			ia.close()
			tb_nonpos=np.where(datatb[:,:,0,0]<=0)
			datatb[:,:,0,0][tb_nonpos]=np.nan
			datatblog=np.log10(datatb[:,:,0,0])
			pos=np.where(np.isnan(datatblog)==True)
			datatblog[pos]=0
			postb=np.where(datatblog>int(np.log10(TB_limit)))
			present_data[:,:,1,:][postb]=0
			present_data[:,:,2,:][postb]=0
			present_data[:,:,3,:][postb]=0
			previous_data[:,:,1,:][postb]=0
			previous_data[:,:,2,:][postb]=0
			previous_data[:,:,3,:][postb]=0
			ia.open(outfile_path+'/'+os.path.basename(previous_image)+'.temp')
			ia.putchunk(previous_data)
			ia.close()
			ia.open(outfile_path+'/'+os.path.basename(present_image)+'.temp')
			ia.putchunk(present_data)
			ia.close()
		dyn_range_present=self.calc_dyn_range(outfile_path+'/'+os.path.basename(present_image)+'.temp',10,stokes_list=['I','Q','U','V'])
		dyn_range_previous=self.calc_dyn_range(outfile_path+'/'+os.path.basename(previous_image)+'.temp',10,stokes_list=['I','Q','U','V'])
		if (dyn_range_present[0]['Q'][-1]-dyn_range_previous[0]['Q'][-1])/dyn_range_previous[0]['Q'][-1]>0.05:
			self.pollog_verbose.info('Q leakage increases.\n')
			present_data_copy[:,:,1,:]=previous_data_copy[:,:,1,:]
			present_model_copy[:,:,1,:]=previous_model_copy[:,:,1,:]
		if (dyn_range_present[0]['U'][-1]-dyn_range_previous[0]['U'][-1])/dyn_range_previous[0]['U'][-1]>0.05:
			self.pollog_verbose.info('U leakage increases.\n')
			present_data_copy[:,:,2,:]=previous_data_copy[:,:,2,:]
			present_model_copy[:,:,2,:]=previous_model_copy[:,:,2,:]
		if (dyn_range_present[0]['V'][-1]-dyn_range_previous[0]['V'][-1])/dyn_range_previous[0]['V'][-1]>0.05:
			self.pollog_verbose.info('V leakage increases.\n')
			present_data_copy[:,:,3,:]=previous_data_copy[:,:,3,:]
			present_model_copy[:,:,3,:]=previous_model_copy[:,:,3,:]
		if (pres_q_frac-prev_q_frac)/prev_q_frac>0.01 or (pres_u_frac-prev_u_frac)/prev_q_frac>0.01 or (pres_v_frac-prev_v_frac)/prev_v_frac>0.01:
			pol_increase=True
		else:
			pol_increase=False
		if overwrite==True or outputfile_prefix=='':
			ia.open(present_image)
			ia.putchunk(present_data_copy)
			ia.close()
			ia.open(present_model)
			ia.putchunk(present_model_copy)
			ia.close()
			os.system('rm -rf casa*log '+outfile_path+'/I* '+outfile_path+'/'+os.path.basename(present_image)+'.temp '+outfile_path+'/'+os.path.basename(previous_image)+'.temp')
			return present_image,present_model,pol_increase
		else:
			os.system('cp -r '+present_image+' '+outfile_path+'/'+outputfile_prefix+'.image')
			os.system('cp -r '+present_model+' '+outfile_path+'/'+outputfile_prefix+'.model')
			ia.open(outfile_path+'/'+outputfile_prefix+'.image')
			ia.putchunk(present_data_copy)
			ia.close()
			ia.open(outfile_path+'/'+outputfile_prefix+'.model')
			ia.putchunk(present_model_copy)
			ia.close()		
			os.system('rm -rf casa*log '+outfile_path+'/I* '+outfile_path+'/'+os.path.basename(present_image)+'.temp '+outfile_path+'/'+os.path.basename(previous_image)+'.temp')
			return outfile_path+'/'+outputfile_prefix+'.image',outfile_path+'/'+outputfile_prefix+'.model',pol_increase

	def compare_leakage_for_sun(self,present_image='',previous_image='',present_model='',previous_model='',outputfile_prefix='',sigma=10,overwrite=False,qucor_step=False):
		'''
		Function to compare Stokes I to Stokes Q,U,V leakage from quiet Sun part with the previous image

		Parameters
		----------
		present_image : str 
			Name of the present CASA image
		previous_image : str 
			Name of the previous CASA image
		present_model : str 
			Name of the present CASA model
		previous_model : str 
			Name of the previous CASA model
		outputfile_prefix : str 
			Final output file prefix name
		sigma : float 
			Sigma value for rms based thresholding
		overwrite : bool
			Overwrite the present image or not
		qucor_step : bool
 			Performed images based Stokes Q, U correction at last selfcal iteration
		Returns
		-------
		str
			Output file name
		'''
		if os.path.isdir(present_image)==False or os.path.isdir(previous_image)==False or os.path.isdir(present_model)==False or os.path.isdir(previous_model)==False:
			if self.verbose==False:
				print ('Either present or previous image or model is not present.\n')
			self.pollog_verbose.info('Either present or previous image or model is not present.\n')
			return 
		self.pollog_verbose.info('Compairing solar images : '+present_image+' with '+previous_image+'\n')
		self.pollog_verbose.info('QU correction step : '+str(qucor_step)+'\n')
		if present_image[-1]=='/':
			present_image=present_image[:-1]
		if previous_image[-1]=='/':
			previous_image=previous_image[:-1]
		try:
			outfile_path=os.path.basedir(present_image)
		except:
			outfile_path=os.getcwd()
		os.system('rm -rf casa*log '+outfile_path+'/I* '+outfile_path+'/'+os.path.basename(present_image)+'.temp '+outfile_path+'/'+os.path.basename(previous_image)+'.temp')
		ia=image()
		ia.open(present_image)
		present_data=ia.getchunk()
		ia.close()
		ia.open(previous_image)
		previous_data=ia.getchunk()
		ia.close()
		ia=image()
		ia.open(present_model)
		present_modeldata=ia.getchunk()
		ia.close()
		ia.open(previous_model)
		previous_modeldata=ia.getchunk()
		ia.close()
		present_data_copy=copy.deepcopy(present_data)
		previous_data_copy=copy.deepcopy(previous_data)
		present_model_copy=copy.deepcopy(present_modeldata)
		previous_model_copy=copy.deepcopy(previous_modeldata)
		os.system('cp -r '+present_image+' '+outfile_path+'/'+os.path.basename(present_image)+'.temp')
		os.system('cp -r '+previous_image+' '+outfile_path+'/'+os.path.basename(previous_image)+'.temp')
		if qucor_step==False:
			imsubimage(imagename=present_image,outfile=outfile_path+'/I.image',stokes='I',dropdeg=False)
			fluxcal_image=self.mwa_solar_fluxcal(imagename=outfile_path+'/I.image',outfile=outfile_path+'/I_fluxcal.image')
			major=imhead(imagename=fluxcal_image)['restoringbeam']['major']['value'] # In arcsec
			minor=imhead(imagename=fluxcal_image)['restoringbeam']['minor']['value'] # In arcsec
			freq=imhead(imagename=fluxcal_image)['refval'][-1]/10**9 # In GHz
			ref_beam_axes_multi=500000
			beam_axes_multi=major*minor	
			expr='1.222e6*IM0/'+str(freq)+'^2/('+str(major*minor)+')'
			immath(imagename=fluxcal_image,outfile=outfile_path+'/I_Tb.image',mode='evalexpr',expr=expr)
			imhead(imagename=outfile_path+'/I_Tb.image', mode='put', hdkey='bunit', hdvalue='K')	
			scale_tb_limit=(beam_axes_multi/ref_beam_axes_multi)
			ia.open(outfile_path+'/I_Tb.image')
			datatb=ia.getchunk()
			ia.close()
			tb_nonpos=np.where(datatb[:,:,0,0]<=0)
			datatb[:,:,0,0][tb_nonpos]=np.nan
			datatblog=np.log10(datatb[:,:,0,0])
			pos=np.where(np.isnan(datatblog)==True)
			datatblog[pos]=0
			postb=np.where(datatblog>=7)
			postb1=np.where(datatblog<6)
			postb2=np.where((datatblog>=6) & (datatblog<7))
			present_data[:,:,1,:][postb]=0
			present_data[:,:,2,:][postb]=0
			present_data[:,:,3,:][postb]=0
			previous_data[:,:,1,:][postb]=0
			previous_data[:,:,2,:][postb]=0
			previous_data[:,:,3,:][postb]=0
			ia.open(outfile_path+'/'+os.path.basename(previous_image)+'.temp')
			ia.putchunk(previous_data)
			ia.close()
			ia.open(outfile_path+'/'+os.path.basename(present_image)+'.temp')
			ia.putchunk(present_data)
			ia.close()
		dyn_range_present=self.calc_dyn_range(outfile_path+'/'+os.path.basename(present_image)+'.temp',sigma,stokes_list=['I','Q','U','V'])
		dyn_range_previous=self.calc_dyn_range(outfile_path+'/'+os.path.basename(previous_image)+'.temp',sigma,stokes_list=['I','Q','U','V'])
		prev_q_frac=dyn_range_previous[0]['Q'][-1]/dyn_range_previous[0]['I'][-1]
		pres_q_frac=dyn_range_present[0]['Q'][-1]/dyn_range_present[0]['I'][-1]
		prev_u_frac=dyn_range_previous[0]['U'][-1]/dyn_range_previous[0]['I'][-1]
		pres_u_frac=dyn_range_present[0]['U'][-1]/dyn_range_present[0]['I'][-1]
		prev_v_frac=dyn_range_previous[0]['V'][-1]/dyn_range_previous[0]['I'][-1]
		pres_v_frac=dyn_range_present[0]['V'][-1]/dyn_range_present[0]['I'][-1]
		if qucor_step==False:
			if ((dyn_range_present[0]['Q'][-1]-dyn_range_previous[0]['Q'][-1])/dyn_range_previous[0]['Q'][-1]>0.05) or (pres_q_frac-prev_q_frac)>0.05:
				self.pollog_verbose.info('Q leakage increases.\n')
				present_data_copy[:,:,1,:]=previous_data_copy[:,:,1,:]
				present_model_copy[:,:,1,:]=previous_model_copy[:,:,1,:]
			if ((dyn_range_present[0]['U'][-1]-dyn_range_previous[0]['U'][-1])/dyn_range_previous[0]['U'][-1]>0.05) or (pres_u_frac-prev_u_frac)>0.05:
				self.pollog_verbose.info('U leakage increases.\n')
				present_data_copy[:,:,2,:]=previous_data_copy[:,:,2,:]
				present_model_copy[:,:,2,:]=previous_model_copy[:,:,2,:]
		if ((dyn_range_present[0]['V'][-1]-dyn_range_previous[0]['V'][-1])/dyn_range_previous[0]['V'][-1]>0.05) or (pres_v_frac-prev_v_frac)>0.05:
			self.pollog_verbose.info('V leakage increases.\n')
			present_data_copy[:,:,3,:]=previous_data_copy[:,:,3,:]
			present_model_copy[:,:,3,:]=previous_model_copy[:,:,3,:]
		if (pres_q_frac-prev_q_frac)>0.05 or (pres_u_frac-prev_u_frac)>0.05 or (pres_v_frac-prev_v_frac)>0.05:
			pol_increase=True
		else:
			pol_increase=False
		if overwrite==True or outputfile_prefix=='':
			ia.open(present_image)
			ia.putchunk(present_data_copy)
			ia.close()
			ia.open(present_model)
			ia.putchunk(present_model_copy)
			ia.close()
			os.system('rm -rf casa*log '+outfile_path+'/I* '+outfile_path+'/'+os.path.basename(present_image)+'.temp '+outfile_path+'/'+os.path.basename(previous_image)+'.temp')
			return present_image,present_model,pol_increase
		else:
			os.system('cp -r '+present_image+' '+outfile_path+'/'+outputfile_prefix+'.image')
			os.system('cp -r '+present_model+' '+outfile_path+'/'+outputfile_prefix+'.model')
			ia.open(outfile_path+'/'+outputfile_prefix+'.image')
			ia.putchunk(present_data_copy)
			ia.close()
			ia.open(outfile_path+'/'+outputfile_prefix+'.model')
			ia.putchunk(present_model_copy)
			ia.close()		
			os.system('rm -rf casa*log '+outfile_path+'/I* '+outfile_path+'/'+os.path.basename(present_image)+'.temp '+outfile_path+'/'+os.path.basename(previous_image)+'.temp')
			return outfile_path+'/'+outputfile_prefix+'.image',outfile_path+'/'+outputfile_prefix+'.model',pol_increase

	def make_crosshand_phase_caltable(self,cross_phase=15,caltable='',polbasis='Linear'): # TODO : Implement circular basis as well
		'''
		Function to make cross hand phase Jones matrices

		Parameters
		----------
		cross_phase : float 
			Cross hand phase different in degree
		caltable : str 
			Name of the caltable to store cross hand Jones matrices
		polbasis : str 
			'Linear' or 'Circular' (Circular basis not implemented)
		Returns
		-------
		str
			Caltable name
		'''
		cross_phase=np.deg2rad(cross_phase)/2.0
		cross_jones=np.matrix([[np.cos(cross_phase)+1j*np.sin(cross_phase),0],[0,np.cos(cross_phase)-1j*np.sin(cross_phase)]])
		AM=AccessMS(self.msname)
		nant=AM.get_num_antenna()
		nchan=AM.get_num_channels()
		nint=AM.get_num_timestamps()
		freqs=AM.get_freqs()/10**6
		start_freq=freqs[0]
		end_freq=freqs[-1]
		mjdsecs=AM.get_timestamps_in_mjdsecs()[0]
		startmjd=mjdsecs[0]
		endmjd=mjdsecs[-1]
		cwd=os.getcwd()
		if caltable!='':
			if os.path.dirname(caltable)=='':
				caltable=cwd+'/'+caltable
			header = struct.pack("8s",b"MWAOCAL")+struct.pack("i",0)+struct.pack("i",0)+struct.pack("i",int(nint))+struct.pack("i",int(nant))+struct.pack("i",int(nchan))+\
					struct.pack("i",4)+struct.pack("d",0.0)+struct.pack("d",0.0)
			jones_array=np.array([np.array([np.real(inv(cross_jones)).flatten(),np.imag(inv(cross_jones)).flatten()]).flatten(order='F') for a in range(nchan)])
			jones_array=jones_array.reshape((1,1,jones_array.shape[0],jones_array.shape[-1]))
			numpy_data=np.repeat(np.repeat(jones_array,nint,axis=0),nant,axis=1)
			numpy_data=numpy_data.flatten(order='C')
			fil=open(caltable,'wb')
			fil.write(header)
			fil.close()
			with open(caltable,mode='ba+') as f:
				numpy_data.tofile(f,format='np.float64')
			bin_data=np.fromfile(caltable,dtype=np.float64)
			os.system('rm -rf '+caltable)
			np.save(caltable,np.array([bin_data,start_freq,end_freq,startmjd,endmjd,nchan,nint],dtype='object'))
			os.system('mv '+caltable+'.npy '+caltable)
			return caltable,cross_jones
		else:
			return None,cross_jones

	def apply_cross_hand_phase(self,cross_phase=15,caltable='',polbasis='Linear',datacolumn='DATA',modify_datacolumn=False):
		'''
		Function to applycal cross hand phase solution

		Parameters
		----------
		cross_phase : float 
			Cross hand phase different in degree
		caltable : str 
			Name of the caltable to store cross hand Jones matrices
		polbasis : str 
			'Linear' or 'Circular' (Circular basis not implemented)
		datacolumn : str
			Datacolumn to apply the solution ('DATA', or 'CORRECTED_DATA')
		modify_datacolumn : bool 
			Modify the DATA column or not
		Returns
		-------
		str
			Caltable name
		'''	
		if caltable=='':
			caltable=os.getcwd()+'/'+os.path.basename(self.msname).split('.ms')[0]+'_cross_phase.temp'
		self.pollog_verbose.info('make_crosshand_phase_caltable(cross_phase='+str(cross_phase)+',caltable=\''+caltable+'\',polbasis=\''+polbasis+'\')\n')
		caltable_name,cross_jones=self.make_crosshand_phase_caltable(cross_phase=cross_phase,caltable=caltable,polbasis=polbasis)
		cal=CALIBRATE()
		self.pollog_verbose.info('cal.applycal(msname=\''+self.msname+'\',gaintable=\''+caltable_name+'\',datacolumn=\''+datacolumn+'\',applymode=\'calonly\',verbose='+\
			str(self.verbose)+')\n')
		cal.applycal(msname=self.msname,gaintable=caltable_name,datacolumn=datacolumn,applymode='calonly',verbose=self.verbose) # Applying the solution
		if modify_datacolumn==True:
			tb=table()
			tb.open(self.msname,nomodify=False)
			cor_data=tb.getcol('CORRECTED_DATA')
			tb.putcol('DATA',cor_data)
			tb.flush()
			tb.close()
		os.system('rm -rf '+caltable)
		return caltable_name

	def wsclean_to_casaimage(self,wsclean_images=[],casaimage_prefix='CASA',imagetype='image',keep_wsclean_images=True): #TODI : Include other stokes mode
		'''
		Function to convert WSClean image in CASA image (Stokes modes : 'IQUV', 'XXYY', 'I', 'QU', 'IV')

		Parameters
		----------
		wsclean_images : list 
			List of WSClean images
		casaimage_prefix : str 
			Output CASA image name prefix (default : 'CASA\_')
		imagetype : str 
			'image', 'model', 'residual' or 'dirty' (default : 'image')
		keep_wsclean_images : bool 
			Keep the WSClean images or not
		Returns
		-------
		str
			Output CASA imagename
		'''
		stokes=[]
		wsclean_images=sorted(wsclean_images)
		for i in wsclean_images:
			name_split=i.split('.fits')[0].split('-')
			if len(name_split)>=3:
				if name_split[-2] not in stokes:
					stokes.append(name_split[-2])
			else:
				if 'I' not in stokes:
					stokes.append('I')
		stokes=sorted(stokes)
		if stokes!=['I','Q','U','V'] and stokes!=['XX','YY'] and stokes!=['I','V'] and stokes!=['Q','U'] and stokes!=['I']:
			if self.verbose:
				self.pollog_verbose.info('Stokes axes are not in \'IQUV\',\'I\',\'QU\',\'IV\' or \'XX,YY\'\n')
			else:
				print('Stokes axes are not in \'IQUV\',\'I\',\'QU\',\'IV\' or \'XX,YY\'\n')
		elif stokes==['I']:
			imagename=casaimage_prefix+'.'+imagetype
			if os.path.isdir(imagename):
				os.system('rm -rf '+imagename)
			importfits(fitsimage=wsclean_images[0],imagename=imagename,defaultaxes=True,defaultaxesvalues=['ra','dec','stokes','freq'])
			if keep_wsclean_images==False:
				for i in wsclean_images:
					os.system('rm -rf '+i)
			return imagename
		elif stokes==['I','V']:
			for i in wsclean_images:
				if i.split('-')[1]=='I':
					data=fits.getdata(i)
					header=fits.getheader(i)
				else:
					data=np.append(data,fits.getdata(i),axis=0)
			header['NAXIS4']=2.
			header['CRVAL4']=1.
			header['CDELT4']=3.
			fits.writeto(casaimage_prefix+'_IV.fits',data=data,header=header,overwrite=True)
			imagename=casaimage_prefix+'.'+imagetype
			if os.path.isdir(imagename):
				os.system('rm -rf '+imagename)
			importfits(fitsimage=casaimage_prefix+'_IV.fits',imagename=imagename,defaultaxes=True,defaultaxesvalues=['ra','dec','stokes','freq'])
			os.system('rm -rf '+casaimage_prefix+'_IV.fits')
			if keep_wsclean_images==False:
				for i in wsclean_images:
					os.system('rm -rf '+i)
			return imagename
		elif stokes==['I','Q','U','V']:
			for i in wsclean_images:
				if i.split('-')[1]=='I':
					data=fits.getdata(i)
					header=fits.getheader(i)
				else:
					data=np.append(data,fits.getdata(i),axis=0)
			header['NAXIS4']=4.
			header['CRVAL4']=1.
			header['CDELT4']=1.
			fits.writeto(casaimage_prefix+'_IQUV.fits',data=data,header=header,overwrite=True)
			imagename=casaimage_prefix+'.'+imagetype
			if os.path.isdir(imagename):
				os.system('rm -rf '+imagename)
			importfits(fitsimage=casaimage_prefix+'_IQUV.fits',imagename=imagename,defaultaxes=True,defaultaxesvalues=['ra','dec','stokes','freq'])
			os.system('rm -rf '+casaimage_prefix+'_IQUV.fits')
			if keep_wsclean_images==False:
				for i in wsclean_images:
					os.system('rm -rf '+i)
			return imagename
		elif stokes==['XX','YY']:
			for i in wsclean_images:
				if i.split('-')[1]=='XX':
					data=fits.getdata(i)
					header=fits.getheader(i)
				else:
					data=np.append(data,fits.getdata(i),axis=0)
			header['NAXIS4']=2.
			header['CRVAL4']=-5.
			header['CDELT4']=-1.
			fits.writeto(casaimage_prefix+'_XXYY.fits',data=data,header=header,overwrite=True)
			imagename=casaimage_prefix+'.'+imagetype
			if os.path.isdir(imagename):
				os.system('rm -rf '+imagename)
			importfits(fitsimage=casaimage_prefix+'_XXYY.fits',imagename=imagename,defaultaxes=True,defaultaxesvalues=['ra','dec','stokes','freq'])
			os.system('rm -rf '+casaimage_prefix+'_XXYY.fits')
			if keep_wsclean_images==False:
				for i in wsclean_images:
					os.system('rm -rf '+i)
			return imagename
		elif stokes==['Q','U']:
			for i in wsclean_images:
				if i.split('-')[1]=='Q':
					data=fits.getdata(i)
					header=fits.getheader(i)
				else:
					data=np.append(data,fits.getdata(i),axis=0)
			header['NAXIS4']=2.
			header['CRVAL4']=2.
			header['CDELT4']=1.
			fits.writeto(casaimage_prefix+'_QU.fits',data=data,header=header,overwrite=True)
			imagename=casaimage_prefix+'.'+imagetype
			if os.path.isdir(imagename):
				os.system('rm -rf '+imagename)
			importfits(fitsimage=casaimage_prefix+'_QU.fits',imagename=imagename,defaultaxes=True,defaultaxesvalues=['ra','dec','stokes','freq'])
			os.system('rm -rf '+casaimage_prefix+'_QU.fits')
			if keep_wsclean_images==False:
				for i in wsclean_images:
					os.system('rm -rf '+i)
			return imagename

	def polselfcal_iteration(self,num_iter,rms_thresh,mask_str,sigma,maskfile,antenna_to_use,startmodel,startmask,want_auto_masking=False,TB_limit=-1,solar_imaging=True,\
							stokes='',interactive=False,use_ankflagger=False,do_flag=False,poldistortion_correction=True,poldistortion_type='poldistortion',crossphase=-1,\
							poldistortion_matrix='UH',do_solarqu_cor=False,box_width=3,previous_image='',previous_model='',calibrator_caltable=[],maskregion='',\
							quvcor_stokes='QUV',polmodel_threshold=-1,cpus=5,absmem=2,wlayers=1,weight='briggs',robust=0.5):
		'''
		Function to perform a polarisation self-calibration loop, make an image, put the model in the measurement set, and perform the calibration

		Parameters
		----------
		num_iter : int 
			Number of self-calibration iteration
		rms_thresh : float 
			RMS for threshold
		maskstr : str 
			Mask string for CLEANing (Only for CASA)
		sigma : float 
			Threshold sigma
		maskfile : str 
			Maskfile for CLEANing (Only for CASA)
		antenna_to_use : list 
			List of antennas for CLEANing
		startmodel : str 
			Model to start the CLEANing (Only for CASA)
		startmask : str 
			Mask to start (Only for CASA)
		want_auto_masking : bool 
			If True use CASA auto-multithresh for auto masking
		TB_limit : float 
			Brightness temperature limit to calculate polarised flux to compare between two polarisation self-calibration rounds (Only for non MWA solar observations)
		solar_imaging : bool 
			Performing Solar image calibration
		stokes : str
			Stokes plane to image
		interactive : bool 
			Perform interactive CLEAN (Only for CASA)
		use_ankflagger : bool 
			Use aNKflagger for flagging after each selfcal round
		do_flag: bool
			Flag after selfcal round
		poldistortion_correction : bool 
			Correct poldistortion using the known ideal Jones matrix of the instrument
		poldistortion_type : str 
			'polconversion ; Stokes I to STOKES Q,U,V leakages' or 'polrotation; changes between Stokes Q,U,V' or 'poldistortion' (default : poldistortion)		
		crossphase : float 
			Cross hand phase in degree for image based correction
		poldistortion_matrix : str 
			'UH or HU ' , where H is polconversion and U is polrotation
		do_solarqu_cor : bool
			Correct solar Stokes I to Q,U imaged based leakage correction
		box_width : float 
			Length of negative box width in degree (default : 3 degree)
		previous_image : str 
			Name of the previous round image to compare leakage
		previous_model : str 
			Name of the previouos round model
		calibrator_caltable : list 
			List of calilbrator caltables
		maskregion : str 
			Mask region in case of auto-masking (Only for CASA)
		quvcor_stokes : str
			Stokes plane for images based leakage correction
		polmodel_threshold : float 
			Sigma value for thresholding on the polarisation model
		cpus : int
			Number of cpu threads to use
		absmem : float
			Memory in GB to use
		wlayers : int
			Number of w-stacking layers (For wsclean only)
		weight : str
			Visibility weighting during imaging , 'uniform', 'natural' or 'briggs'. Default is 'briggs'
		robust : float
			Robust parameter for briggs weighting, default : 0.5 
		Returns
		-------
		int
			Message code 
		dict
			Dynamic range information [DR dictionary : {'STOKES':[rms dynamic range,rms,total_flux]}]
		float
			Negative based dynamic range 
		'''
		os.chdir(self.mspath)
		cal=CALIBRATE()	
		available_cpus=psutil.cpu_count()
		available_memory=psutil.virtual_memory()[1]/10**9
		if cpus>=available_cpus:
			cpus=int(available_cpus*0.7)
		if absmem>=available_memory:
			absmem=int(available_memory*0.7)
		imagename=self.msname.split('.ms')[0]+'_'+str(num_iter) # Imagename prefix
		present_file=glob.glob(imagename+'*')
		if len(present_file)!=0:
			os.system('rm -rf '+imagename+'*')
			self.pollog_verbose.info('rm -rf '+imagename+'*\n')
		caltable_name=self.msname.split('.ms')[0]+'.bin' # Caltable name
		if os.path.isfile(caltable_name):
			os.system('rm -rf '+caltable_name)
		threshold=[str(rms*sigma)+'Jy' for rms in rms_thresh]
		self.pollog_verbose.info('==============================\n')
		self.pollog_verbose.info('Iteration number : '+str(num_iter)+'\n')
		self.pollog_verbose.info('==============================\n')
		if self.wsclean==False:
			# Making image
			if maskfile=='':
				maskfile=mask_str
			robust=robust*2
			if maskfile!='':
				self.pollog_verbose.info('poltclean(vis=\''+self.msname+'\',imagename=\''+imagename+'\',selectdata=True,startmodel=\''+startmodel+'\',startmask=\''\
						+startmask+'\',stokes=\''+stokes+'\',antenna=\''+antenna_to_use+'\',imsize=['+str(self.imsize)+'],cell=\''+str(self.cellsize)\
						+'arcsec\',niter=10000,gain=0.08,threshold='+str(threshold)+',deconvolver=\'multiscale\',scales='+str(self.multiscale_scales)
						+',uvtaper=\''+self.uvtaper+'\',weighting=\''+weight+'\',robust='+str(robust)+',interactive=False,mask=\''+str([maskfile])+'\')\n')
				poltclean(vis=self.msname,imagename=imagename,selectdata=True,startmodel=startmodel,startmask=startmask,stokes=stokes,antenna=antenna_to_use,imsize=[self.imsize],\
				cell=str(self.cellsize)+'arcsec',niter=10000,gain=0.08,threshold=threshold,deconvolver='multiscale',scales=self.multiscale_scales,\
				uvtaper=self.uvtaper,weighting=weight,robust=robust,interactive=False,mask=[maskfile])
			elif want_auto_masking==True and maskfile=='': # Use auto-masking
				try_count=0
				if startmodel!='': # Add auto masking safety
					automask_trials=2
				else:
					automask_trials=10
				while True:
					if try_count==0:
						self.pollog_verbose.info('Normal auto-masking.\n')
						if startmodel!='': # Add auto masking safety
							automask_trials=2
						else:
							automask_trials=10
						self.pollog_verbose.info('poltclean(vis=\''+self.msname+'\',imagename=\''+imagename+'\',selectdata=True,startmodel=\''+startmodel+\
							'\',startmask=\''+startmask+'\',stokes=\''+stokes+'\',antenna=\''+antenna_to_use+'\',imsize=['+str(self.imsize)+'],cell=\''+str(self.cellsize)+\
							'arcsec\',niter=10000,gain=0.08,threshold='+str(threshold)+',deconvolver=\'multiscale\',scales='+str(self.multiscale_scales)+\
							',uvtaper=\''+str(self.uvtaper)+'\',weighting=\''+weight+'\',robust='+str(robust)+',interactive=False,usemask='+\
							'\'auto-multithresh\',mask=\'\',pbmask=0.0,sidelobethreshold='+str(float(3.0))+',noisethreshold='+str(float(sigma))+\
							',lownoisethreshold='+str(float(sigma/3.0))+',negativethreshold='+str(float(sigma))+',smoothfactor=1.0,minbeamfrac=0.1,growiterations=75,'+\
							'minpercentchange=5.0,automask_trials='+str(automask_trials)+',maskregion=\''+maskregion+'\')\n')
						poltclean(vis=self.msname,imagename=imagename,selectdata=True,startmodel=startmodel,startmask=startmask,stokes=stokes,antenna=antenna_to_use,\
						imsize=[self.imsize],cell=str(self.cellsize)+'arcsec',niter=10000,gain=0.08,threshold=threshold,deconvolver='multiscale',scales=self.multiscale_scales,\
						uvtaper=self.uvtaper,weighting=weight,robust=robust,interactive=False,usemask='auto-multithresh',mask='',pbmask=0.0,sidelobethreshold=float(3.0),\
						noisethreshold=float(sigma),lownoisethreshold=float(sigma/3.0),negativethreshold=float(sigma),smoothfactor=1.0,minbeamfrac=0.1,growiterations=75,\
						minpercentchange=5.0,automask_trials=automask_trials,maskregion=maskregion)
					elif try_count==1:
						self.pollog_verbose.info('Trying with auto-masking with no restriction of minimum beam fraction.\n')
						self.pollog_verbose.info('poltclean(vis=\''+self.msname+'\',imagename=\''+imagename+'\',selectdata=True,startmodel=\''+startmodel+\
							'\'startmask=\''+startmask+'\',,stokes=\''+stokes+'\',antenna=\''+antenna_to_use+'\',imsize=['+str(self.imsize)+'],cell=\''+str(self.cellsize)+\
							'arcsec\',niter=10000,gain=0.08,threshold='+str(threshold)+',deconvolver=\'multiscale\',scales='+str(self.multiscale_scales)+\
							',uvtaper=\''+str(self.uvtaper)+'\',weighting=\''+weight+'\',robust='+str(robust)+',interactive=False,usemask=\''+\
							'auto-multithresh\',mask=\'\',pbmask=0.0,sidelobethreshold='+str(float(3.0))+',noisethreshold='+str(float(sigma))+\
							',lownoisethreshold='+str(float(sigma/3.0))+',negativethreshold='+str(float(sigma))+\
							',smoothfactor=1.0,minbeamfrac=0.1,growiterations=75,minpercentchange=-1.0,automask_trials=\''+str(automask_trials)+'\')\n')
						poltclean(vis=self.msname,imagename=imagename,selectdata=True,startmodel=startmodel,startmask=startmask,stokes=stokes,antenna=antenna_to_use,\
						imsize=[self.imsize],cell=str(self.cellsize)+'arcsec',niter=10000,gain=0.08,threshold=threshold,deconvolver='multiscale',scales=self.multiscale_scales,\
						uvtaper=self.uvtaper,weighting=weight,robust=robust,interactive=False,usemask='auto-multithresh',mask='',pbmask=0.0,sidelobethreshold=float(3.0),\
						noisethreshold=float(sigma),lownoisethreshold=float(sigma/3.0),negativethreshold=float(sigma),smoothfactor=1.0,minbeamfrac=0.1,growiterations=75,\
						minpercentchange=-1.0,automask_trials=automask_trials,maskregion=maskregion)
					elif try_count==2:
						self.pollog_verbose.info('Trying without masking.\n')
						self.pollog_verbose.info('poltclean(vis=\''+self.msname+'\',imagename=\''+imagename+'\',selectdata=True,startmodel=\''+startmodel+\
							'\',startmask=\''+startmask+'\',stokes=\''+stokes+'\',antenna=\''+antenna_to_use+'\',imsize=['+str(self.imsize)+'],cell=\''+str(self.cellsize)+\
							'arcsec\',niter=10000,gain=0.08,threshold='+str(threshold)+',deconvolver=\'multiscale\',scales='+str(self.multiscale_scales)+\
							',uvtaper=\''+str(self.uvtaper)+'\',weighting=\''+weight+'\',robust='+str(robust)+',interactive=False,usemask=\'user\',mask=\'\')\n')
						poltclean(vis=self.msname,imagename=imagename,selectdata=True,startmodel=startmodel,startmask=startmask,stokes=stokes,antenna=antenna_to_use,\
						imsize=[self.imsize],cell=str(self.cellsize)+'arcsec',niter=10000,gain=0.08,threshold=threshold,deconvolver='multiscale',scales=self.multiscale_scales,\
						uvtaper=self.uvtaper,weighting=weight,robust=robust,interactive=False,usemask='user',mask='')
					else:
						break
					modelflux=imstat(imagename=imagename+'.model')['sum'][0]
					if modelflux==0.0:
						try_count+=1
					else:
						break
			else: # If no masking	
				self.pollog_verbose.info('poltclean(vis=\''+self.msname+'\',imagename=\''+imagename+'\',selectdata=True,startmodel=\''+startmodel+\
						'\',startmask=\''+startmask+'\',stokes=\''+stokes+'\',antenna=\''+antenna_to_use+'\',imsize=['+str(self.imsize)+'],cell=\''+str(self.cellsize)\
						+'arcsec\',niter=10000,gain=0.08,threshold=\''+str(threshold)+'\',deconvolver=\'multiscale\',scales='+str(self.multiscale_scales)
						+',uvtaper=\''+self.uvtaper+'\',weighting=\''+weight+'\',robust='+str(robust)+',interactive=False,mask='+str([maskfile])+')\n')	
				poltclean(vis=self.msname,imagename=imagename,selectdata=True,startmodel=startmodel,startmask=startmask,stokes=stokes,antenna=antenna_to_use,imsize=[self.imsize],\
				cell=str(self.cellsize)+'arcsec',niter=10000,gain=0.08,threshold=threshold,deconvolver='multiscale',scales=self.multiscale_scales,\
				uvtaper=self.uvtaper,weighting=weight,robust=robust,interactive=False,mask='')
			casa_imagename=imagename+'.image'
			casa_modelname=imagename+'.model'
			casa_residualname=imagename+'.residual'
		else:
			scales=[str(i) for i in self.multiscale_scales]
			if weight=='briggs':
				weight=weight+' '+str(robust)
			wsclean_args=['-scale '+str(self.cellsize)+'asec','-size '+str(self.imsize)+' '+str(self.imsize),'-no-dirty','-j '+str(cpus),\
			'-abs-mem '+str(absmem),'-weight '+weight,'-taper-tukey 10','-name '+imagename,'-nwlayers '+str(wlayers),'-pol iquv','-maxuv-l '+str(self.imaging_maxuv),\
			'-minuv-l '+str(self.imaging_minuv),'-niter 100000','-mgain 0.8','-auto-threshold '+str(sigma),'-auto-mask '+str(sigma+0.5),'-gain 0.08','-multiscale',\
			'-multiscale-scales '+','.join(scales)]
			wsclean_args.append('-quiet')
			if self.verbose:
				self.pollog_verbose.info('wsclean '+' '.join(wsclean_args)+' '+self.msname)			
			os.system('wsclean '+' '.join(wsclean_args)+' '+self.msname)			
			wsclean_images=glob.glob(imagename+'*image*')
			wsclean_models=glob.glob(imagename+'*model*')
			wsclean_residuals=glob.glob(imagename+'*residual*')
			casa_imagename=self.wsclean_to_casaimage(wsclean_images=wsclean_images,casaimage_prefix=imagename,imagetype='image',keep_wsclean_images=False)
			casa_modelname=self.wsclean_to_casaimage(wsclean_images=wsclean_models,casaimage_prefix=imagename,imagetype='model',keep_wsclean_images=False)
			casa_residualname=self.wsclean_to_casaimage(wsclean_images=wsclean_residuals,casaimage_prefix=imagename,imagetype='residual',keep_wsclean_images=False)
			os.system('rm -rf *psf*')
		if do_solarqu_cor==True and solar_imaging==True:
			quvcor_imagename,quvcor_modelname,q_change_frac,u_change_frac,v_change_frac=self.correct_solar_quv_leakage(casa_imagename,casa_modelname,sigma,overwrite=True,\
																														stokes=quvcor_stokes)
			poldistortion_correction=False
			if os.path.isdir('quvcor.image')==False:
				os.system('cp -r '+quvcor_imagename+' quvcor.image')
			if os.path.isdir('quvcor.model')==False:
				os.system('cp -r '+quvcor_modelname+' quvcor.model')
		if polmodel_threshold!=-1:
			self.pol_model_threshold(casa_imagename,casa_modelname,sigma,polmodel_threshold/sigma)
			if do_solarqu_cor==True and solar_imaging==True:
				self.pol_model_threshold('quvcor.image','quvcor.model',sigma,polmodel_threshold/sigma)
		if previous_image!='' and previous_model!='' and solar_imaging==True:
			self.compare_leakage_for_sun(present_image=casa_imagename,previous_image=previous_image,present_model=casa_modelname,previous_model=previous_model,overwrite=True,\
										qucor_step=do_solarqu_cor)
			if do_solarqu_cor==True:
				self.compare_leakage_for_sun(present_image='quvcor.image',previous_image=previous_image,present_model='quvcor.model',previous_model=previous_model,\
						overwrite=True,qucor_step=True)
		elif previous_image!='' and previous_model!='' and solar_imaging==False:
			self.compare_leakages(present_image=casa_imagename,previous_image=previous_image,present_model=casa_modelname,previous_model=previous_model,overwrite=True,\
				TB_limit=TB_limit)
		if crossphase!=-1:
			outimage,outmodel=self.correct_image_for_cross_phase(casa_imagename,casa_modelname,casa_imagename+'.image.temp',cross_phase=crossphase,\
							imagetype='CASA',outtype='CASA',pol_basis='Linear',do_fluxcal=True)
			os.system('mv '+outimage+' '+casa_imagename)
			os.system('mv '+outmodel+' '+casa_modelname)	
			if do_solarqu_cor==True:
				outimage,outmodel=self.correct_image_for_cross_phase('quvcor.image','quvcor.model','quvcor.image.temp',cross_phase=crossphase,\
						imagetype='CASA',outtype='CASA',pol_basis='Linear',do_fluxcal=True)
				os.system('mv '+outimage+' quvcor.image')
				os.system('mv '+outmodel+' quvcor.model')
		out_dict,negative_dyn_range=self.calc_dyn_range(casa_imagename,sigma,box_width=box_width,stokes_list=['I','Q','U','V']) # Calculating the dynamic range of the image
		out_dict_keys=out_dict.keys()
		if 'NAN' in out_dict_keys:
			self.pollog_verbose.info(B.error_msgs(3))
			os.system('rm -rf casa*log')
			os.chdir(self.cwd)
			os.system('rm -rf casa*log')
			return 3     # If image is not made, no point in continuing
		if os.path.isdir(casa_modelname)==False:
			self.pollog_verbose.info(B.error_msgs(4))
			os.system('rm -rf casa*log')
			os.chdir(self.cwd)
			os.system('rm -rf casa*log')
			return 4	   # If model is not present no point in continuing
		else:
			modelflux=imstat(imagename=casa_modelname)['sum'][0]
			if modelflux==0.0:
				self.pollog_verbose.info(B.error_msgs(5))
				os.system('rm -rf casa*log')
				os.chdir(self.cwd)
				os.system('rm -rf casa*log')
				return 5 # If modelflux is 0, no point in continuing
			else:
				flaglist=flagmanager(vis=self.msname,mode='list') # Restore flagversion if thr last flag version is from CALIBRATE
				if len(flaglist)>1:
					flaglist_keys=list(flaglist.keys())
					flaglist_keys.remove('MS')
					if len(flaglist_keys)>=2:
						last_flag_key=len(flaglist_keys)-2  # Last flag is either from uvsub_Flagger or aNKflagger
					else:
						last_flag_key=len(flaglist_keys)-1
					last_flagversion=flaglist[last_flag_key]['name']
					# Restore the flag and delete the present flag version
					if 'CALIBRATE_applycal' in last_flagversion:
						self.pollog_verbose.info('flagmanager(vis=\''+self.msname+'\',mode=\'restore\',versionname=\''+str(last_flagversion)+'\',merge=\'replace\')\n')
						flagmanager(vis=self.msname,mode='restore',versionname=last_flagversion,merge='replace')
					for key in flaglist_keys:
						flagversion=flaglist[key]['name']
						self.pollog_verbose.info('flagmanager(vis=\''+self.msname+'\',mode=\'delete\',versionname=\''+str(flagversion)+'\')\n')
						flagmanager(vis=self.msname,mode='delete',versionname=flagversion)
				if num_iter==0:
					tb=table()
					tb.open(self.msname,nomodify=False)
					flag=tb.getcol('FLAG')*False
					tb.putcol('FLAG',flag)
					tb.flush()
					tb.close()	
				self.pollog_verbose.info('clearcal(vis=\''+self.msname+'\')\n')
				clearcal(vis=self.msname)
				self.pollog_verbose.info('delmod(vis=\''+self.msname+'\',scr=True,otf=True)\n') 
				delmod(vis=self.msname,scr=True,otf=True) # Clear the MODEL column
				self.remove_model_negative(casa_imagename,casa_modelname,sigma=sigma,overwrite=True) # Removing negatives from model
				self.pollog_verbose.info('ft(vis=\''+self.msname+'\',model=\''+casa_modelname+'\',nterms=1,usescratch=True)\n') 
				ft(vis=self.msname,model=casa_modelname,nterms=1,usescratch=True) # Putting the model into MS
				if self.verbose==False: # Performing Full Jones calibration
					self.pollog_verbose.info('cal.calibrate(msname=\''+self.msname+'\',caltable=\''+caltable_name+'\',minuv='+str(self.calib_uvrange_min)+',quiet=True,maxuv='\
							+str(self.calib_uvrange_max)+',j=1,absmem=1,solmode=\'\')\n')
					cal.calibrate(msname=self.msname,caltable=caltable_name,minuv=self.calib_uvrange_min,quiet=True,maxuv=self.calib_uvrange_max,\
									j=1,absmem=1,solmode='')
				else:
					self.pollog_verbose.info('cal.calibrate(msname=\''+self.msname+'\',caltable=\''+caltable_name+'\',minuv='+str(self.calib_uvrange_min)+',quiet=False,maxuv='\
							+str(self.calib_uvrange_max)+',j=1,absmem=1,solmode=\'\')\n') 
					cal.calibrate(msname=self.msname,caltable=caltable_name,minuv=self.calib_uvrange_min,quiet=False,maxuv=self.calib_uvrange_max,\
									j=1,absmem=1,solmode='')					
				if os.path.isfile(caltable_name)==False:
					self.pollog_verbose.info('DR_I:'+str(out_dict['I'][0])+', DR_Q:'+str(out_dict['Q'][0])+', DR_U:'+str(out_dict['U'][0])\
												+', DR_V:'+str(out_dict['V'][0])+', DR_neg:'+str(negative_dyn_range)+'\n')
					self.pollog_verbose.info('No caltable made.\n')
					os.system('rm -rf casa*log')
					os.chdir(self.cwd)
					os.system('rm -rf casa*log')
					return 7
				self.pollog_verbose.info('poldistortion_correction='+str(poldistortion_correction)+'\n')
				self.pollog_verbose.info('cal.applycal(msname=\''+self.msname+'\',gaintable=\''+caltable_name+'\',applymode=\'calflag\',flagbackup=True)\n')
				cal.applycal(msname=self.msname,gaintable=caltable_name,applymode='calflag',flagbackup=True) # Applying the solution
				if use_ankflagger==True and do_flag==True:
					try:
						self.pollog_verbose.info('do_uvsub_ankflag(\''+self.msname+'\',nthread=1,verbose='+str(False)+')\n')
						fg.do_uvsub_ankflag(self.msname,nthread=1,verbose=False)
					except Exception as e:
						self.pollog_verbose.info('Error in aNKflagger : '+str(e)+'\n')
						self.pollog_verbose.info('Error in running aNKflagger. Using rms threshold flagging.\n')
						self.pollog_verbose.info('do_uvsub_flagger(\''+self.msname+'\',mode=\'uvsub_flag\',rmsthresh=[10,7,5,3.5])\n')
						fg.do_uvsub_flagger(self.msname,mode='uvsub_flag',rmsthresh=[10,7,5,3.5])
				elif do_flag==True:
					self.pollog_verbose.info('do_uvsub_flagger(\''+self.msname+'\',mode=\'uvsub_flag\',rmsthresh=[10,7,5,3.5])\n')
					fg.do_uvsub_flagger(self.msname,mode='uvsub_flag',rmsthresh=[10,7,5,3.5])
				if poldistortion_correction==True:
					self.pollog_verbose.info('cal_poldistortion(\''+caltable_name+'\',poldistortion_matrix='+poldistortion_matrix+')\n')
					X,inv_X,H,inv_H,U,inv_U,poldist_file=self.cal_poldistortion(caltable_name,poldistortion_matrix=poldistortion_matrix)
					if poldistortion_type=='polconversion':
						self.pollog_verbose.info('poldistortion_type=\'polconversion\'\n')
						self.pollog_verbose.info('correct_poldistortion(\''+caltable_name+'\',\''+caltable_name+'\',H)\n')
						corrected_gaintable=self.correct_poldistortion(caltable_name,caltable_name,H) # Correct for polconversion
					elif poldistortion_type=='polrotation':
						self.pollog_verbose.info('poldistortion_type=\'polrotation\'\n')
						self.pollog_verbose.info('correct_poldistortion(\''+caltable_name+'\',\''+caltable_name+'\',U)\n')
						corrected_gaintable=self.correct_poldistortion(caltable_name,caltable_name,U) # Correct for polrotation
					else:
						self.pollog_verbose.info('poldistortion_type=\'poldistortion\'\n')					
						self.pollog_verbose.info('correct_poldistortion(\''+caltable_name+'\',\''+caltable_name+'\',X)\n')
						corrected_gaintable=self.correct_poldistortion(caltable_name,caltable_name,X) # Correct for full poldistortion
					caltable_name=corrected_gaintable
					self.pollog_verbose.info('Applying poldistortion corrected caltable ........\n')
					self.pollog_verbose.info('cal.applycal(msname=\''+self.msname+'\',gaintable=\''+caltable_name+\
									'\',applymode=\'calflag\',flagbackup=True)\n')
					cal.applycal(msname=self.msname,gaintable=caltable_name,applymode='calflag',flagbackup=True) # Applying the solution
				self.pollog_verbose.info('DR_I:'+str(out_dict['I'][0])+', DR_Q:'+str(out_dict['Q'][0])+', DR_U:'+str(out_dict['U'][0])\
												+', DR_V:'+str(out_dict['V'][0])+', DR_neg:'+str(negative_dyn_range)+'\n')
				if do_solarqu_cor==True:
					tb=table()
					tb.open(self.msname)
					cor_data=tb.getcol('CORRECTED_DATA')
					tb.close()
					tb.open(self.msname,nomodify=False)
					tb.putcol('DATA',cor_data)
					tb.flush()
					tb.close()
				self.pollog_verbose.info('Success.\n')
				os.system('rm -rf casa*log')
				os.chdir(self.cwd)
				os.system('rm -rf casa*log')
				return 0,out_dict,negative_dyn_range

#########################################
# Finished PolSelfcal Class
#########################################
