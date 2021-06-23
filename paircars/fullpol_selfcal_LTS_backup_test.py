import numpy as np,os,copy,sys,glob
import logging
from casatools import *
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
#matplotlib.use('Agg')
'''
Code is written by Devojyoti Kansabanik, 01 Mar, 2021
'''

class PolSelfcal:
	'''
	Generic class to perform polarisation self-calibration (Using Andre Offringa's CALIBRATE code on based Mitchcal algorithm)
	Attributes:
	msname = Name of the measurement set
	maximum_emission_scale = Maximum scale of the emission present in the image
	verbose = False,If True keep all the intermediate images, model, residuals, caltables and details of the log to detailed analysis
	interactive = False, If True user have interactive control on self-calibration
	'''
	
	def __init__(self,msname,metafits,maximum_emission_scale,verbose=False,interactive=False):
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
		self.cellsize=IB.calc_cellsize(3) # Assuming 3 pixels in one PSF
		self.imsize=IB.num_pixels(3)
		self.max_size=maximum_emission_scale
		self.multiscale_scales=IB.choose_scales(3,self.max_size)
		self.uvtaper=IB.calc_uvtaper()
		self.calib_uvrange_min=IB.calc_calib_uvrange(2)[1]
		self.calib_uvrange_max=IB.calc_calib_uvrange(2)[2]
		self.rms_box='50,50,'+str(self.imsize-50)+','+str(int(self.imsize/4)) # CASA box to calculate the rms
		self.verbose=verbose
		self.interactive=interactive
		formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
		self.pollog_verbose = logging.getLogger('polselfcal_verbose_log')
		self.pollog_verbose.setLevel(logging.DEBUG)
		if self.verbose:
			self.console=logging.StreamHandler(sys.stdout)
			self.console.setFormatter(formatter)
			self.pollog_verbose.addHandler(self.console)
		self.filehandle=logging.FileHandler(self.cwd+'/Pol_Selfcal_verbose.log')
		self.filehandle.setFormatter(formatter)
		self.pollog_verbose.addHandler(self.filehandle)
		self.pollog_verbose.propagate = False
		self.pollog_verbose.info('Initiating Polarisation selfcal object.\n')
		
	def negative_box(self,max_pix,box_width=3):
		'''
		Create a 3 degree box about the maximum pixel of image to search negative.
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
		os.system('rm -rf casa*log')
		return out_dict,negative_dyn_range

	def calc_iter_num(self,safety_factor,quality_factor,scratch=True):
		'''
		Function to calculate minimum number of selfcal iteration based on safety standard and quality factor
		Parameters:
		safety_factor = Factor to determine the robustness of the selfcal
		quality_factor = Factor to determine the quality of the images
		scratch : True, whether start the selfcal from scratch or not
		Return:
		Minimum iteration at fixed sigma, Minimum iteration, Maximum iteration , Number of antenna bins, Fraction change in flux for convergence
		'''
		if quality_factor==0:     # Low quality (Quick look image making)
			frac_flux_change=0.2
			if (safety_factor==0):
				min_num_iter_fixed_sigma=1
				if (scratch==True):
					min_iteration=1
					max_iteration=20
				else:
					min_iteration=1
					max_iteration=10
			elif (safety_factor==1):
				min_num_iter_fixed_sigma=1
				if (scratch==True):
					min_iteration=3
					max_iteration=30
				else:
					min_iteration=1
					max_iteration=20
			else:
				min_num_iter_fixed_sigma=1
				if (scratch==True):
					min_iteration=5
					max_iteration=40
				else:
					min_iteration=1
					max_iteration=30
		elif quality_factor==1:  # Medium quality imaging (Computing speed medium)
			frac_flux_change=0.15
			if (safety_factor==0):
				min_num_iter_fixed_sigma=2
				if (scratch==True):
					min_iteration=2
					max_iteration=40
				else:
					min_iteration=2
					max_iteration=30
			elif (safety_factor==1):
				min_num_iter_fixed_sigma=2
				if (scratch==True):
					min_iteration=4
					max_iteration=50
				else:
					min_iteration=2
					max_iteration=40
			else:
				min_num_iter_fixed_sigma=2
				if (scratch==True):
					min_iteration=6
					max_iteration=60
				else:
					min_iteration=2
					max_iteration=50
		else:  # Best quality imaging (Computing slow)
			frac_flux_change=0.1
			if (safety_factor==0):
				max_iteration=60
				min_num_iter_fixed_sigma=3
				if (scratch==True):
					min_iteration=3
				else:
					min_iteration=3
			elif (safety_factor==1):
				min_num_iter_fixed_sigma=3
				if (scratch==True):
					min_iteration=5
				else:
					min_iteration=3
			else:
				min_num_iter_fixed_sigma=3
				if (scratch==True):
					min_iteration=7
				else:
					min_iteration=3
		antenna_bin=1
		self.pollog_verbose.info('Quality factor :'+str(quality_factor)+', Safety standard :'+str(safety_factor)+\
				', Minimum number of iteration at fixed sigma : '+str(min_num_iter_fixed_sigma)+', Minimum iteration :'+str(min_iteration)+', Antenna bins :'+str(antenna_bin)+\
				', Fraction flux change for convergence : '+str(frac_flux_change)+'\n')
		os.system('rm -rf casa*log')
		return min_num_iter_fixed_sigma,min_iteration,max_iteration,antenna_bin,frac_flux_change

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

	def cal_poldistortion(self,gaintable,poldistortion_matrix='UH'): # TODO : Diagnostic plots
		'''
		Function to calculate the estmated Jones matrices for poldistortion after correcting for instrumental Jones matrix
		Parameters:
		gaintable = Name of the gaintable (Assuming CALIBRATE gaintable format only right now)
		poldistortion _matrix ='UH' or 'HU', where H is polconversion matrix and U is the polrotation matrix
		Return:
		Poldistortion matrix, Inverse of poldistortion matrix, Polconversion, Inversion of polconversion, Polrotation, Inverse of polrotation, Filename to save poldistortion matrix 
		NOte : Saved X matrix is for B'=XBX^\dagger, which is inverse of poldistortion_matrix of correct_poldistortion function
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
		s=np.sqrt(det(x))
		x=x/s
		if poldistortion_matrix=='UH':
			U,H=polar(x,side='right')
		else:
			H,U=polar(x,side='left')
		os.system('rm -rf '+gaintable+'.calibrate_bin*')
		os.system('rm -rf casa*log')
		np.save(gaintable+'.poldist',np.array([inv(x)]))   # Saving X matrix for B'=XBX^\dagger, which is inverse of poldistortion_matrix of correct_poldistortion function
		return x,inv(x),np.matrix(H),np.matrix(inv(H)),np.matrix(U),np.matrix(inv(U)),gaintable+'.poldist'

	def correct_poldistortion(self,gaintable,outfile,poldistortion_matrix):
		'''
		Function to applycal poldistortion correction (either polconversion or polrotation or both) to the gaintable
		Parameters:
		gaintable = Name of the gaintable
		outfile = Name of the output poldistortion corrected gaintable		
		poldistortion_matrix = Poldistortion matrix, either in numpy array or numpy matrix form
		Return:
		Poldistortion corrected gaintable
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

	def file_remover_and_keeper(self,num_iter,msg_code,ref_time_chan=True):
		'''
		This function keep and remove caltables, ms, imaging related files based on the need
		Parameters:
		num_iter = Number of self-calibration iteration
		msg_code = Selfcal message code
		ref_timechan = True , reference time channel or not
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
			os.system('cp -r '+file_str+'.mask junk0.mask') # Copying num_iter=0 mask to junk0.mask
			os.system('cp -r '+file_str+'.residual junk0.residual') # Copying num_iter=0 residual to junk0.residual
			if self.verbose and (msg_code==0 or msg_code==8 or msg_code==9):	
				os.system('cp -r '+caltable_name+' '+self.mspath+'/'+file_str_prefix+'/backup_cal/'+file_str+'.bin') # Verbose=True, keep all the caltables
				os.system('cp -r '+self.msname+' '+self.mspath+'/'+file_str_prefix+'/backup_ms/'+file_str+'.ms') # If Verbose=True, keep all the ms
				os.system('cp -r '+file_str+'.model '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.model') # If Verbose=True, keep all the models
				os.system('cp -r '+file_str+'.mask '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.mask') # If Verbose=True, keep all the masks
				os.system('cp -r '+file_str+'.image '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.image') # If Verbose=True, keep all the image
				os.system('cp -r '+file_str+'.residual '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.residual') # If Verbose=True, keep all the residuals
			os.system('rm -rf '+caltable_name)
			os.system('rm -rf '+file_str+'.image '+file_str+'.model '+file_str+'.residual '+file_str+'.sumwt '+file_str+'.pb '+file_str+'.psf '+file_str+'.mask ') 
				# Removing all imaging related files
		elif num_iter==1 and (msg_code==0 or msg_code==8 or msg_code==9):
			os.system('cp -r '+caltable_name+' junk1.bin') # Copying num_iter=1 caltable to junk1.bin
			os.system('cp -r '+self.msname+' junk1.ms') # Copying num_iter=1 ms to junk1.ms
			os.system('cp -r '+file_str+'.image junk1.image') # Copying num_iter=1 image to junk1.image
			os.system('cp -r '+file_str+'.model junk1.model') # Copying num_iter=1 model to junk1.model
			os.system('cp -r '+file_str+'.mask junk1.mask') # Copying num_iter=1 model to junk1.mask
			os.system('cp -r '+file_str+'.residual junk1.residual') # Copying num_iter=1 residual to junk1.residual
			if self.verbose and (msg_code==0 or msg_code==8 or msg_code==9):	
				os.system('cp -r '+caltable_name+' '+self.mspath+'/'+file_str_prefix+'/backup_cal/'+file_str+'.bin') # Verbose=True, keep all the caltables
				os.system('cp -r '+self.msname+' '+self.mspath+'/'+file_str_prefix+'/backup_ms/'+file_str+'.ms') # If Verbose=True, keep all the ms
				os.system('cp -r '+file_str+'.model '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.model') # If Verbose=True, keep all the models
				os.system('cp -r '+file_str+'.mask '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.mask') # If Verbose=True, keep all the masks
				os.system('cp -r '+file_str+'.image '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.image') # If Verbose=True, keep all the image
				os.system('cp -r '+file_str+'.residual '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.residual') # If Verbose=True, keep all the residuals
			os.system('rm -rf '+caltable_name)
			os.system('rm -rf '+file_str+'.image '+file_str+'.model '+file_str+'.residual '+file_str+'.sumwt '+file_str+'.pb '+file_str+'.psf '+file_str+'.mask ') 
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
				os.system('cp -r '+file_str+'.residual junk0.residual') # Copying model to junk0.residual
			os.system('rm -rf junk1.ms junk1.bin junk1.model junk1.mask junk1.image junk1.residual')
			os.system('cp -r '+caltable_name+' junk1.bin') # Copying caltable to junk1.bin
			os.system('cp -r '+self.msname+' junk1.ms') # Copying ms to junk1.ms
			os.system('cp -r '+file_str+'.model junk1.model') # Copying model to junk1.model
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
				os.system('cp -r '+file_str+'.mask '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.mask') # If Verbose=True, keep all the masks
				if os.path.exists(self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.image'):
					os.system('rm -rf '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.image')
				os.system('cp -r '+file_str+'.image '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.image') # If Verbose=True, keep all the image
				if os.path.exists(self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.residual'):
					os.system('rm -rf '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.residual')
				os.system('cp -r '+file_str+'.residual '+self.mspath+'/'+file_str_prefix+'/backup_imagemodel/'+file_str+'.residual') # If Verbose=True, keep all the residuals
			os.system('rm -rf '+file_str+'.image '+file_str+'.model '+file_str+'.residual '+file_str+'.sumwt '+file_str+'.pb '+file_str+'.psf '+file_str+'.mask ') 
				# Removing all imaging related files
		os.chdir(self.cwd)
		os.system('rm -rf casa*log')
		return

	def estimateSkyBrightnessMatrix(self,beam_jones,Vij):
		''' 
		Return beam corrected brightness matrix
		Parameters:
		beam_jones = Beam Jones matrix
		Vij = Instrumental brightness matrix
		Return:
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
	
	def get_IQUV(self,imagename,imagetype='FITS',stokes='IQUV'):
		'''
		Stokes I,Q,U,V from a Stokes IQUV image cube.
		Parameters:
		imagename = Name of the image
		imagetype= Type of the image, CASA or FITS
		stokes = 'IQUV' , to options 'IQUV' for full Stokes or 'I' for only intensity for unpolarised source
		Return:
		Python dictionary {'STOKES':imagedata}
		'''
		if imagetype=='CASA':
			if os.path.isfile(imagename.split('.image')[0]+'.fits'):
				os.system('rm -rf '+imagename.split('.image')[0]+'.fits')
			exportfits(imagename=imagename,fitsimage=imagename.split('.image')[0]+'.fits')
			fitsimage=imagename.split('.image')[0]+'.fits'
		else:
			fitsimage=imagename
		data=fits.getdata(fitsimage)
		stokes_data = {}
		if stokes=='IQUV' and data.shape[0]==4:
			stokes_data['I'] = data[0, 0, :, :]
			stokes_data['Q'] = data[1, 0, :, :]
			stokes_data['U'] = data[2, 0, :, :]
			stokes_data['V'] = data[3, 0, :, :]
			os.system('rm -rf casa*log')
			return stokes_data
		elif stokes=='I' and data.shape[0]==1:
			stokes_data['I'] = data[0, 0, :, :]
			stokes_data['Q'] = data[0, 0, :, :]*0
			stokes_data['U'] = data[0, 0, :, :]*0
			stokes_data['V'] = data[0, 0, :, :]*0
			os.system('rm -rf casa*log')
			return stokes_data
		else:
			print ('Stokes parameter does not match with the Stokes parameters of your image.\n')
			os.system('rm -rf casa*log')
			return stokes_data
		
	def get_inst_pols(self,stokes_image,imagetype='FITS',pol_basis='Linear',stokes='IQUV'): #TODO : Circular basis
		'''
		Return instrumental polarisation matrix (Vij)
		Parameters:
		stokes_image = Name of the Stokes IQUV image cube
		imagetype= Type of the image, CASA or FITS
		pol_basis = Polarisation basis of the instrument, Linear or Circular
		stokes = 'IQUV' , to options 'IQUV' for full Stokes or 'I' for only intensity for unpolarised source
		Return:
		Instrumental polarisation matrix
		'''
		stokes_data=self.get_IQUV(stokes_image,imagetype=imagetype,stokes=stokes)
		XX = stokes_data['I'] + stokes_data['Q']
		XY = stokes_data['U'] + stokes_data['V'] * 1j
		YX = stokes_data['U'] - stokes_data['V'] * 1j
		YY = stokes_data['I'] - stokes_data['Q']
		Vij = np.array([[XX, XY], [YX, YY]])
		os.system('rm -rf casa*log')
		return Vij

	def B_to_IQUV(self,B,pol_basis='Linear'): # TODO : Circular Basis
		'''
		Convert brightness matrix in instrumental basis to I, Q, U, V
		Parameters:
		B = Brightness matrix in instrumental basis
		pol_basis = Polarisation basis of the instrument, Linear or Circular
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

	def correct_for_single_beam_jones(self,imagename,outfile,beam_jones,imagetype='FITS',outtype='FITS',pol_basis='Linear',stokes='IQUV'): #TODO : Circular basis
		'''
		Correct Stokes IQUV image cube for full Stokes Beam Jones at a single pointing
		Parameters:
		imagename = Name of the image of model
		outfile = Name of the beam corrected image or model
		beam_jones = Beam jones matrix
		imagetype= Type of the image, CASA or FITS
		outtype = Output image type, CASA or FITS
		pol_basis = Polarisation basis of the instrument, Linear or Circular
		stokes = 'IQUV' , to options 'IQUV' for full Stokes or 'I' for only intensity for unpolarised source
		Return:
		Beam corrected image or model
		'''
		Vij=self.get_inst_pols(imagename,imagetype=imagetype,pol_basis=pol_basis,stokes=stokes)
		B=self.estimateSkyBrightnessMatrix(beam_jones,Vij)
		stokes=self.B_to_IQUV(B,pol_basis='Linear')
		if os.path.exists(outfile):
			os.system('rm -rf '+outfile)
			os.system('rm -rf '+'.'.join(outfile.split('.')[:-1])+'.fits')
		if imagetype=='CASA' and os.path.isfile('.'.join(imagename.split('.')[:-1])+'.fits')==False:
			exportfits(imagename=imagename,fitsimage='.'.join(imagename.split('.')[:-1])+'.fits',stokeslast=False)
			os.system('cp -r '+imagename+' temp_org.image')
		imagename='.'.join(imagename.split('.')[:-1])+'.fits'
		fitsdata=fits.getdata(imagename)
		header=fits.getheader(imagename)
		data=np.empty((1,4,fitsdata.shape[-2],fitsdata.shape[-1]))
		data[0,0,:,:]=np.real(stokes['I'])
		data[0,1,:,:]=np.real(stokes['Q'])
		data[0,2,:,:]=np.real(stokes['U'])
		data[0,3,:,:]=np.real(stokes['V'])
		fits.writeto(outfile,data=data,header=header,overwrite=True)
		if outtype=='CASA':
			#if os.path.isdir('temp_org.image')==False:
			#	importfits(fitsimage=imagename,imagename='temp_org.image')
			importfits(fitsimage=outfile,imagename='temp.image')
			os.system('rm -rf '+outfile+' '+imagename)
			ia=image()
			ia.open('temp.image')
			pbcor_data=ia.getchunk()
			ia.close()
			#ia.open('temp_org.image')
			#ia.putchunk(pbcor_data)
			#ia.done()
			#ia.close()
			os.system('mv temp.image '+outfile)
		os.system('rm -rf casa*log temp*')
		return outfile

	def uncorrect_for_single_beam_jones(self,imagename,outfile,inv_beam_jones,imagetype='FITS',outtype='FITS',pol_basis='Linear',stokes='IQUV'): # TODO : Circular basis
		'''
		Undo the beam correction for Stokes IQUV image cube for full Stokes Beam Jones at a single pointing
		Parameters:
		imagename = Name of the image of model
		outfile = Name of the beam corrected image or model
		inv_beam_jones = Inverse of Beam jones matrix
		imagetype = Type of the image, CASA or FITS
		outtype = Output image type, CASA or FITS
		pol_basis = Polarisation basis of the instrument, Linear or Circular
		stokes = 'IQUV' , to options 'IQUV' for full Stokes or 'I' for only intensity for unpolarised source
		Return:
		Beam un-corrected image or model
		'''
		Vij=self.get_inst_pols(imagename,imagetype=imagetype,pol_basis=pol_basis,stokes=stokes)
		B=self.estimateSkyBrightnessMatrix(inv_beam_jones,Vij)
		stokes=self.B_to_IQUV(B,pol_basis='Linear')
		if os.path.exists(outfile):
			os.system('rm -rf '+outfile)
			os.system('rm -rf '+'.'.join(outfile.split('.')[:-1])+'.fits')
		if imagetype=='CASA' and os.path.isfile('.'.join(imagename.split('.')[:-1])+'.fits')==False:
			exportfits(imagename=imagename,fitsimage='.'.join(imagename.split('.')[:-1])+'.fits')
			os.system('cp -r '+imagename+' temp_org.image')
		imagename='.'.join(imagename.split('.')[:-1])+'.fits'
		fitsdata=fits.getdata(imagename)
		header=fits.getheader(imagename)
		print (header)
		data=np.empty((4,1,fitsdata.shape[-2],fitsdata.shape[-1]))
		data[0,0,:,:]=np.real(stokes['I'])
		data[1,0,:,:]=np.real(stokes['Q'])
		data[2,0,:,:]=np.real(stokes['U'])
		data[3,0,:,:]=np.real(stokes['V'])
		fits.writeto(outfile,data=data,header=header,overwrite=True)
		if outtype=='CASA':
			#if os.path.isdir('temp_org.image')==False:
			#	importfits(fitsimage=imagename,imagename='temp_org.image')
			importfits(fitsimage=outfile,imagename='temp.image')
			os.system('rm -rf '+outfile+' '+imagename)
			ia=image()
			ia.open('temp.image')
			pbcor_data=ia.getchunk()
			ia.close()
			#ia.open('temp_org.image')
			#ia.putchunk(pbcor_data)
			#ia.done()
			#ia.close()
			os.system('mv temp.image '+outfile)
		os.system('rm -rf casa*log temp*')
		return outfile

	def correct_visibility_single_beam_jones(self,modify_datacolumn=True,force=False,skip_freq=1.28,save_beamfile=''):
		'''
		Correct visibility data for a single pointing beam jones
		Parameters:
		modify_datacolumn = True, modify the DATA column, otherwise beam corrected visibilities will be saved on CORRECTED_DATA
		force = False, beam correct forcefully avoiding ms header info
		skip_freq = Frequency interval in MHz to make independent beams (default : 1.28 MHz). If anything greater than 1.28 MHz is given it will be overwritten to 1.28 MHz
		save_beamfile= = Save beam file in this given name
		Return:
		Name of the beam jones file, Beam Jones matrix.
		'''
		if os.path.exists(save_beamfile)==True:
			os.system('rm -rf '+save_beamfile)
		mwapb=MWA_PrimaryBeam(self.msname,self.metafits,inverse_beam=False)  #TODO : Beam per coarse channel for multi coarse chan ms
		cal=CALIBRATE()
		beamfile=self.msname+'.beam.bin'
		beamfile,beamjones=mwapb.MWA_phasecenter_beam_jones(outputfile=beamfile,skip_freq=float(skip_freq))
		if save_beamfile!='' and save_beamfile!=beamfile and os.path.exists(save_beamfile)==False:
			os.system('cp -r '+beamfile+' '+save_beamfile)
		code=vishead(vis=self.msname,mode='get',hdkey='fld_code')[0][0]
		code_list=code.split(',')
		if 'S_PBCOR' not in code_list or 'S_PBUNCOR' in code_list:
			cal.applycal(msname=self.msname,gaintable=beamfile,applymode='calonly') # Applying the beam correction
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
			cal.applycal(msname=self.msname,gaintable=beamfile,applymode='calonly') # Applying the beam correction
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
		Parameters:
		modify = True, modify the DATA column, otherwise beam corrected visibilities will be saved on CORRECTED_DATA
		force = False, undo beam correct forcefully avoiding ms header info
		skip_freq = Frequency interval in MHz to make independent beams (default : 1.28 MHz). If anything greater than 1.28 MHz is given it will be overwritten to 1.28 MHz
		Return:
		Name of the beam jones file
		'''
		mwapb=MWA_PrimaryBeam(self.msname,self.metafits,inverse_beam=True)
		cal=CALIBRATE()
		beamfile=self.msname+'.beam.bin'
		beamfile,beamjones=mwapb.MWA_phasecenter_beam_jones(outputfile=beamfile,skip_freq=float(skip_freq))
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

	def uncorrect_visibility_model_single_beam_jones(self,force=False,skip_freq=1.28):	
		'''
		Undo Correct visibility data for a single pointing beam jones
		Parameters:
		force = False, undo beam correct forcefully avoiding ms header info
		skip_freq = Frequency interval in MHz to make independent beams (default : 1.28 MHz). If anything greater than 1.28 MHz is given it will be overwritten to 1.28 MHz
		Return:
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

	def IMSTAT_record(self,DRI,DR_neg,FXI,FXQ,FXU,FXV,FXT,FXP,record_filename,init=True):
		'''
		Function to keep the record of image statistics at different self calibration steps
		Parameters:
		DRI = RMS based dynamic range of the Stokes I
		DR_neg = Negative based dynamic range of the Stokes I
		FXI = Total Stokes I flux
		FXQ = Total Stokes Q flux
		FXU = Total Stokes U flux
		FXV = Total Stokes V flux
		FXT = Total Stokes T flux
		FXP = Total Stokes P flux
		record_filename = Name of the file to store dynamic ranges
		init = True, initiating a new record from the current selfcal iteration
		Return:
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

	def reduce_sigma(self,imagename,nsigma,sigma_step,minsigma,pre_residual=0.0,residual_frac=0.1,stokes_list=['I']):
		'''
		Function to determine whether reduce the CLEAN sigma or not.
		Parameters:
		imagename = Name of the image
		nsigma = Value of the present n-sigma
		sigma_step = Step to reduce sigma value
		minsigma = Minimum allowed sigma
		pre_residual = Previous residual fraction to compare (default : 0.0)
		residual_frac = Residual flux fraction to reduce sigma (default : 0.1)
		stokes_list = ['I'], stokes plane list
		Return:
		Reduced value of n-sigma and median residual fraction if residual flux is more than given percentage (default : 10%) of the total flux in Stokes I or in all Stokes Q,U,V.
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
					image_pix_sum=0
					residual_pix_sum=1
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
			if (residual_pix_sum/image_pix_sum>residual_frac and residual_pix_sum/image_pix_sum<pre_residual) or ((max_frac_diff>0 and max_frac_diff>residual_frac) or \
				((min_frac_diff>0 and min_frac_diff>residual_frac) and stokes!='I' and stokes!='XX' and stokes!='YY')) :
				do_reduce_list.append(1)
		os.system('rm -rf reduce_sigma_*')
		os.chdir(cwd)
		residual_frac_median=np.median(np.array(residual_frac_list))
		if int(np.sum(np.array(do_reduce_list)))>1:
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
		Polarisation minimisation function fir Sun
		Parameters:
		datai = Stokes I image data
		datal = Stokes Q or U image data	
		l = Trial leakage (-1,1)
		rmsl = RMS of the Stokes map
		i_flux = Mean brightness in Jy/beam
		Return:
		The number of pixels having polarisation fraction greater than rmsl/i_flux
		'''
		x1=np.abs((datal-l*datai)/datai)
		pos=np.where(x1.flatten()>rmsl/i_flux)
		f_out=(len(pos[0]))
		del x1,datai,datal,l,rmsl
		return f_out

	def solarcir_pol_minimise(self,datai,datav,l,rmsv,mean_i_flux):
		'''
		Polarisation minimisation function fir Sun
		Parameters:
		datai = Stokes I image data
		datal = Stokes Q or U image data	
		l = Trial leakage (-1,1)
		rmsv = RMS of the Stokes V map
		mean_i_flux = Mean brightness in Jy/beam
		Return:
		The number of pixels having polarisation fraction greater than rmsl/i_flux
		'''
		x1=np.abs((datav-l*datai)/datai)
		if (3*rmsv)/mean_i_flux>0.001:
			threshold=(3*rmsv)/mean_i_flux
		else:
			threshold=0.001
		pos=np.where(x1.flatten()>threshold)
		f_out=(len(pos[0]))
		del x1,datai,datav,l
		return f_out

	def subtract_leakage_surface(self,imagename,modelname,sigma=10,do_fluxcal=False,overwrite=False):
		'''
		Function to subtract quadratic leakage surface
		Parameters:
		imagename = Name of the image
		modelname = Name of the model
		sigma = N-sigma value above which any emission is considered to be real
		Return:
		Leakage surface subtracted image and model
		'''
		if os.path.exists('qucor_surface.image'):
			os.system('rm -rf qucor_surface.image')
		if os.path.exists('qucor_surface.model'):
			os.system('rm -rf qucor_surface.model')
		if overwrite==False:
			os.system('cp -r '+imagename+' '+'qucor_surface.image')
			os.system('cp -r '+modelname+' '+'qucor_surface.model')
			imagename='qucor_surface.image'
			modelname='qucor_surface.model'
		self.pollog_verbose.info('Correcting Stokes Q,U leakage surface.\n')
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
		if scale_tb_limit!=0:
			s='%.2E' % (10**6/scale_tb_limit)
			limit=int(float(s.split('E')[0]))*10**(int(float(s.split('E')[-1])))
			postb=np.where(((datatb[:,:,0,0]/(10**6)).astype('int'))*10**6>=limit)
			postb1=np.where(((datatb[:,:,0,0]/(10**6)).astype('int'))*10**6<limit)
		else:
			limit=10**6
			postb=np.where(((datatb[:,:,0,0]/(10**6)).astype('int'))*10**6>=limit)
			postb1=np.where(((datatb[:,:,0,0]/(10**6)).astype('int'))*10**6<limit)
		datai[posi]=np.nan
		dataq[posi]=np.nan
		datai[posq]=np.nan
		dataq[posq]=np.nan
		datai[postb]=np.nan
		dataq[postb]=np.nan
		q_by_i=dataq/datai
		pos_meanq=np.where(np.abs(q_by_i)>(np.nanmedian(np.abs(q_by_i))+5*np.nanstd(q_by_i)))
		datai[pos_meanq]=np.nan
		dataq[pos_meanq]=np.nan
		q_by_i=dataq/datai
		q_by_i=np.delete(q_by_i,np.where(np.isnan(q_by_i)))
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
		AQ = np.c_[np.ones(q_stack.shape[0]), q_stack[:,:2], np.prod(q_stack[:,:2], axis=1), q_stack[:,:2]**2]
		CQ,_,_,_ = scipy.linalg.lstsq(AQ, q_stack[:,2])	
		for k in range(datai_copy.shape[0]):
			for l in range(datai_copy.shape[1]):
				if np.isnan(datai_mask[k,l])==False:
					dataq_copy[k,l] -= (CQ[4]*k**2. + CQ[5]*l**2. + CQ[3]*k*l + CQ[1]*k + CQ[2]*l + CQ[0])*datai_copy[k,l]
		datai=copy.deepcopy(datai_copy)
		datau=copy.deepcopy(datau_copy)
		datai[posi]=np.nan
		dataq[posi]=np.nan
		datai[posu]=np.nan
		datau[posu]=np.nan
		datai[postb]=np.nan
		datau[postb]=np.nan
		u_by_i=datau/datai
		pos_meanu=np.where(np.abs(u_by_i)>(np.nanmedian(np.abs(u_by_i))+5*np.nanstd(u_by_i)))
		datai[pos_meanu]=np.nan
		datau[pos_meanu]=np.nan
		u_by_i=datau/datai
		u_by_i=np.delete(u_by_i,np.where(np.isnan(u_by_i)))
		x=[]
		y=[]
		z=[]
		u_by_i=datau/datai
		for k in range(datai.shape[0]):
			for l in range(datai.shape[1]):
				if np.isnan(u_by_i[k,l])==False:
					x.append(k)
					y.append(l)
					z.append(u_by_i[k,l])
		u_stack=np.vstack((x,y,z)).T
		del x,y,z
		AU = np.c_[np.ones(u_stack.shape[0]), u_stack[:,:2], np.prod(u_stack[:,:2], axis=1), u_stack[:,:2]**2]
		CU,_,_,_ = scipy.linalg.lstsq(AU, u_stack[:,2])
		for k in range(datai_copy.shape[0]):
			for l in range(datai_copy.shape[1]):
				if np.isnan(datai_mask[k,l])==False:
					datau_copy[k,l] -= (CU[4]*k**2. + CU[5]*l**2. + CU[3]*k*l + CU[1]*k + CU[2]*l + CU[0])*datai_copy[k,l]
		'''
		X,Y = np.meshgrid(np.arange(625, 665, 1), np.arange(625, 665, 1))
		XX = X.flatten()
		YY = Y.flatten()
		fig = plt.figure()
		ax = plt.axes(projection='3d')
		Z=(CU[4]*X**2. + CU[5]*Y**2. + CU[3]*X*Y + CU[1]*X + CU[2]*Y + CU[0])
		ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
		ax.scatter(u_stack[:,0], u_stack[:,1], u_stack[:,2], c='r', s=10)
		plt.xlabel('X')
		plt.ylabel('Y')
		ax.set_zlabel('Z')
		ax.axis('tight')
		plt.show()
		plt.imshow(datau_copy)
		plt.colorbar()
		plt.show()'''
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
		modelq[postb1]=0
		modelu[postb1]=0
		modeldata[:,:,0,0]=model_datai
		modeldata[:,:,1,0]=modelq
		modeldata[:,:,2,0]=modelu
		ia.open(modelname)
		ia.putchunk(modeldata)
		ia.close()
		model_datau[postb]=np.nan
		if overwrite==True:
			if os.path.isdir('qucor_surface.image')==True:
				os.system('rm -rf qucor_surface.image')
			if os.path.isdir('qucor_surface.model')==True:
				os.system('rm -rf qucor_surface.model')
			os.system('cp -r '+imagename+' '+'qucor_surface.image')
			os.system('cp -r '+modelname+' '+'qucor_surface.model')
		os.system('rm -rf casa*log I*.image')
		os.system('rm -rf '+outfile_path+'/I*')
		del modeldata,data,datai_mask,dataq,datai,datau,model_datai,model_dataq,model_datau,datai_copy,dataq_copy,datau_copy
		return imagename,modelname		

	def mwa_solar_fluxcal(self,imagename,outfile): # TODO :Flux scale needs to include
		immath(imagename=imagename,outfile=outfile,mode='evalexpr',expr='IM0*50')
		return outfile

	def cal_solar_qu_leakage(self,imagename,sigma=10,do_fluxcal=False): 
		'''
		Function to calculate Stokes QU leakage for solar observation (Not vaild for any other astrophysical observation)
		Parameters:
		imagename = Name of the image
		do_fluxcal = Do flux calibration or not
		outfile_path = Name of the directory to save leakage data numpy table (default : image directory)
		Return:
		Stokes Q leakage, Stokes U leakage (Two numpy table will also be saved) 
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
		posq=np.where(dataq<(sigma*rmsq))
		posu=np.where(datau<(sigma*rmsu))
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
		if scale_tb_limit!=0:
			s='%.2E' % (10**6/scale_tb_limit)
			limit=int(float(s.split('E')[0]))*10**(int(float(s.split('E')[-1])))
			postb=np.where(((datatb[:,:,0,0]/(10**6)).astype('int'))*10**6>=limit)
			postb1=np.where(((datatb[:,:,0,0]/(10**6)).astype('int'))*10**6<limit)
		else:
			limit=10**6
			postb=np.where(((datatb[:,:,0,0]/(10**6)).astype('int'))*10**6>=limit)
			postb1=np.where(((datatb[:,:,0,0]/(10**6)).astype('int'))*10**6<limit)
		datai[postb]=np.nan
		dataq[postb]=np.nan
		datau[postb]=np.nan
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
		os.system('rm -rf '+outfile_path+'/I_*')
		self.pollog_verbose.info('Calculation of Stokes leakages has been done.\n')
		return leakage_list[0],leakage_list[1]

	def create_circular_mask(self,h, w, center=None, radius=None):
		'''
		Function to create a circular mask
		Parameters:
		h = Number of pixels in Y
		w = Number of pixels in X
		center = (x_cen,y_cen), center of the mask circle
		radius = Radius in number of pixels
		Return:
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
		Parameters:
		imagename = Name of the image
		do_fluxcal = Do flux calibration or not
		outfile_path = Name of the directory to save leakage data numpy table (default : image directory)
		Return:
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
		posv=np.where(np.abs(datav)<(5*rmsv))
		datai[posi]=np.nan
		datav[posi]=np.nan
		datai[posv]=np.nan
		datav[posv]=np.nan
		ia.open(outfile_path+'/I_Tb.image')
		datatb=ia.getchunk()
		ia.close()
		if scale_tb_limit!=0:
			s='%.2E' % (10**6/scale_tb_limit)
			limit=int(float(s.split('E')[0]))*10**(int(float(s.split('E')[-1])))
			postb=np.where(((datatb[:,:,0,0]/(10**6)).astype('int'))*10**6>=limit)
			postb1=np.where(((datatb[:,:,0,0]/(10**6)).astype('int'))*10**6<limit)
		else:
			limit=10**6
			postb=np.where(((datatb[:,:,0,0]/(10**6)).astype('int'))*10**6>=limit)
			postb1=np.where(((datatb[:,:,0,0]/(10**6)).astype('int'))*10**6<limit)
		datai[postb]=np.nan
		datav[postb]=np.nan
		mean_i_flux=np.nanmean(datai)
		h,w=datav.shape[:2]
		center=(int(h/2),int(w/2))
		radius=int((16*60)/float(self.cellsize))
		mask=self.create_circular_mask(h,w,center=center,radius=radius)
		datav[~mask]=np.nan
		datai[~mask]=np.nan
		unmasked_pixels=np.sum(np.isnan(datav)==False)		
		area=3.14*radius**2	
		leakage_list=[]					
		if (unmasked_pixels/area)<=0.5:
			np.save(outfile_path+'/'+os.path.basename(imagename)+'_V_leakage',np.array([0,0,0],dtype=object))
			leakage_list.append(0)
			os.system('rm -rf casa*log I*.image')
			os.system('rm -rf '+outfile_path+'/I_*')
			return 0	
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
				r=self.solarcir_pol_minimise(datai,datav,l,rmsv,mean_i_flux)
				y.append(r)
				results.append(r)
			results=np.array(results)
			minval=l_range[np.argmin(results)]
			start_range=(start_range+minval)/3.0
			end_range=(end_range+minval)/3.0
		y=np.array(y)
		leakage=x[np.argmin(y)]
		np.save(outfile_path+'/'+os.path.basename(imagename)+'_V_leakage',np.array([x,y,leakage],dtype=object))
		leakage_list.append(leakage)
		os.system('rm -rf casa*log I*.image')
		os.system('rm -rf '+outfile_path+'/I_*')
		return leakage

	def correct_solar_quv_leakage(self,imagename,modelname,sigma,overwrite=False):
		'''
		Function to correction for solar Stokes Q, U and V leakage 
		(It is based on the fact that we do not expect any linear polarisation from the Quiet Sun emission)
		(We also do expect less than 0.5% polarisation from Quiet Sun about 1 solar radio from disc center)
		Parameters:
		imagename = Name of image
		modelname = Name of the model
		sigma = N-sigma threshold to choose Stokes I emission region
		overwrite = False, overwrite the image and model or not
		Return:
		Stokes Q ,U and V leakage corrected image and model name
		'''
		if os.path.exists('quvcor.image'):
			os.system('rm -rf quvcor.image')
		if os.path.exists('quvcor.model'):
			os.system('rm -rf quvcor.model')
		if overwrite==False:
			os.system('cp -r '+imagename+' '+'quvcor.image')
			os.system('cp -r '+modelname+' '+'quvcor.model')
			imagename='quvcor.image'
			modelname='quvcor.model'
		imagename,modelname=self.subtract_leakage_surface(imagename,modelname,sigma=sigma,do_fluxcal=True,overwrite=True)
		v_leakage=self.cal_solar_v_leakage(imagename,sigma=sigma,do_fluxcal=True) #TODO: fluxcal or not
		ia=image()
		ia.open(imagename) # Correcting image
		data=ia.getchunk()
		datai=data[:,:,0,:]
		dataq=data[:,:,1,:]
		datau=data[:,:,2,:]
		datav=data[:,:,3,:]-(v_leakage*datai)
		data[:,:,3,:]-=(v_leakage*datai)
		ia.putchunk(data)
		ia.close()
		ia.open(modelname) # Correcting model
		datam=ia.getchunk()
		im=datam[:,:,0,:]
		qm=datam[:,:,1,:]
		um=datam[:,:,2,:]
		vm=datam[:,:,3,:]-(v_leakage*im)
		posqm=np.where(qm==0)
		posum=np.where(um==0)
		posvm=np.where(vm==0)
		rmsq=imstat(imagename=imagename,box=self.rms_box,stokes='Q')['rms'][0]
		rmsu=imstat(imagename=imagename,box=self.rms_box,stokes='U')['rms'][0]
		rmsv=imstat(imagename=imagename,box=self.rms_box,stokes='V')['rms'][0]
		rmsi=imstat(imagename=imagename,box=self.rms_box,stokes='I')['rms'][0]
		posi=np.where(datai<(sigma*rmsi))
		posq=np.where(np.abs(dataq)<(sigma*rmsq))
		posu=np.where(np.abs(datau)<(sigma*rmsu))
		posv=np.where(np.abs(datau)<(sigma*rmsv))
		qm[posi]=0
		qm[posq]=0
		um[posi]=0
		um[posu]=0
		vm[posi]=0
		vm[posv]=0
		datam[:,:,1,:]=qm
		datam[:,:,2,:]=um
		datam[:,:,3,:]=vm
		ia.putchunk(datam)
		ia.close()
		if overwrite==True:
			if os.path.isdir('quvcor.image')==True:
				os.system('rm -rf quvcor.image')
			if os.path.isdir('quvcor.model')==True:
				os.system('rm -rf quvcor.model')
			os.system('cp -r '+imagename+' quvcor.image')
			os.system('cp -r '+modelname+' quvcor.model')
		os.system('rm -rf casa*log qucor_surface*')
		return imagename,modelname

	def polselfcal_iteration(self,num_iter,rms_thresh,mask_str,sigma,maskfile,antenna_to_use,startmodel,startmask,want_auto_masking=False,\
							stokes='',interactive=False,use_ankflagger=False,poldistortion_correction=True,poldistortion_type='poldistortion',\
							poldistortion_matrix='UH',do_solarquv_cor=False,box_width=3,calibrator_caltable=[]):
		'''
		Function to perform a polarisation self-calibration loop, make an image, put the model in the measurement set, and perform the calibration
		Parameters:
		num_iter = Number of self-calibration iteration
		rms_thresh = RMS for threshold
		maskstr = Mask string for CLEANing
		sigma = Threshold sigma
		maskfile = Maskfile for CLEANing
		antenna_to_use = List of antennas for CLEANing
		startmodel = Model to start the CLEANing
		startmask = Mask to start
		want_auto_masking = False, if True use CASA auto-multithresh for auto masking
		stokes = '', Stokes plane to image
		interactive= False, Perform interactive CLEAN
		use_ankflagger = False, use aNKflagger for flagging after each selfcal round
		poldistortion_correction = True, Correct poldistortion using the known ideal Jones matrix of the instrument
		poldistortion_type = 'polconversion ; Stokes I to STOKES Q,U,V leakages' or 'polrotation; changes between Stokes Q,U,V' or 'poldistortion' (default : poldistortion)
		poldistortion_matrix = 'UH or HU ' , where H is polconversion and U is polrotation
		do_solarquv_cor = False, correct solar Stokes I to Q,U imaged based leakage correction
		box_width = Length of negative box width in degree (default : 3 degree)
		calibrator_caltable = List of calilbrator caltables
		Return:
		Message code, DR dictionary, negative based dynamic range [DR dictionary : {'STOKES':[rms dynamic range,rms,total_flux]}]
		'''
		os.chdir(self.mspath)
		cal=CALIBRATE()	
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
		# Making image
		if maskfile=='':
			maskfile=mask_str
		if maskfile!='':
			self.pollog_verbose.info('poltclean(vis=\''+self.msname+'\',imagename=\''+imagename+'\',selectdata=True,startmodel=\''+startmodel+'\',startmask=\''\
					+startmask+'\',stokes=\''+stokes+'\',antenna=\''+antenna_to_use+'\',imsize=['+str(self.imsize)+'],cell=\''+str(self.cellsize)\
					+'arcsec\',niter=100000000000,gain=0.1,threshold='+str(threshold)+',deconvolver=\'multiscale\',scales='+str(self.multiscale_scales)
					+',uvtaper=\''+self.uvtaper+'\',weighting=\'natural\',interactive=False,mask=\''+str([maskfile])+'\')\n')
			poltclean(vis=self.msname,imagename=imagename,selectdata=True,startmodel=startmodel,startmask=startmask,stokes=stokes,antenna=antenna_to_use,imsize=[self.imsize],\
			cell=str(self.cellsize)+'arcsec',niter=100000000000,gain=0.1,threshold=threshold,deconvolver='multiscale',scales=self.multiscale_scales,\
			uvtaper=self.uvtaper,weighting='natural',interactive=False,mask=[maskfile])
		elif want_auto_masking==True and maskfile=='': # Use auto-masking
			try_count=0
			while True:
				if try_count==0:
					self.pollog_verbose.info('Normal auto-masking.\n')
					self.pollog_verbose.info('poltclean(vis=\''+self.msname+'\',imagename=\''+imagename+'\',selectdata=True,startmodel=\''+startmodel+\
						'\',startmask=\''+startmask+'\',stokes=\''+stokes+'\',antenna=\''+antenna_to_use+'\',imsize=['+str(self.imsize)+'],cell=\''+str(self.cellsize)+\
						'arcsec\',niter=100000000000,gain=0.1,threshold='+str(threshold)+',deconvolver=\'multiscale\',scales='+str(self.multiscale_scales)+\
						',uvtaper=\''+str(self.uvtaper)+'\',weighting=\'natural\',interactive=False,usemask='+\
						'\'auto-multithresh\',mask=\'\',pbmask=0.0,sidelobethreshold=1.5,noisethreshold=3.0,lownoisethreshold=1.5,negativethreshold=3.0,smoothfactor=1.0,'+\
						'minbeamfrac=0.1,growiterations=75,minpercentchange=5.0)\n')
					poltclean(vis=self.msname,imagename=imagename,selectdata=True,startmodel=startmodel,startmask=startmask,stokes=stokes,antenna=antenna_to_use,imsize=[self.imsize],\
					cell=str(self.cellsize)+'arcsec',niter=100000000000,gain=0.1,threshold=threshold,deconvolver='multiscale',scales=self.multiscale_scales,\
					uvtaper=self.uvtaper,weighting='natural',interactive=False,usemask='auto-multithresh',mask='',pbmask=0.0,sidelobethreshold=1.5,noisethreshold=3.0,\
					lownoisethreshold=1.5,negativethreshold=3.0,smoothfactor=1.0,minbeamfrac=0.1,growiterations=75,minpercentchange=5.0)
				elif try_count==1:
					self.pollog_verbose.info('Trying with auto-masking with no restriction of minimum beam fraction.\n')
					self.pollog_verbose.info('poltclean(vis=\''+self.msname+'\',imagename=\''+imagename+'\',selectdata=True,startmodel=\''+startmodel+\
						'\'startmask=\''+startmask+'\',,stokes=\''+stokes+'\',antenna=\''+antenna_to_use+'\',imsize=['+str(self.imsize)+'],cell=\''+str(self.cellsize)+\
						'arcsec\',niter=100000000000,gain=0.05,threshold='+str(threshold)+',deconvolver=\'multiscale\',scales='+str(self.multiscale_scales)+\
						',uvtaper=\''+str(self.uvtaper)+'\',weighting=\'natural\',interactive=False,usemask=\''+\
						'auto-multithresh\',mask=\'\',pbmask=0.0,sidelobethreshold=1.5,noisethreshold=3.0,lownoisethreshold=1.5,negativethreshold=3.0,smoothfactor=1.0,'+\
						'minbeamfrac=0.1,growiterations=75,minpercentchange=-1.0)\n')
					poltclean(vis=self.msname,imagename=imagename,selectdata=True,startmodel=startmodel,startmask=startmask,stokes=stokes,antenna=antenna_to_use,imsize=[self.imsize],\
					cell=str(self.cellsize)+'arcsec',niter=100000000000,gain=0.1,threshold=threshold,deconvolver='multiscale',scales=self.multiscale_scales,uvtaper=self.uvtaper,\
					weighting='natural',interactive=False,usemask='auto-multithresh',mask='',pbmask=0.0,sidelobethreshold=1.5,noisethreshold=3.0,lownoisethreshold=1.5,\
					negativethreshold=3.0,smoothfactor=1.0,minbeamfrac=0.1,growiterations=75,minpercentchange=-1.0)
				elif try_count==2:
					self.pollog_verbose.info('Trying without masking.\n')
					self.pollog_verbose.info('poltclean(vis=\''+self.msname+'\',imagename=\''+imagename+'\',selectdata=True,startmodel=\''+startmodel+\
						'\',startmask=\''+startmask+'\',stokes=\''+stokes+'\',antenna=\''+antenna_to_use+'\',imsize=['+str(self.imsize)+'],cell=\''+str(self.cellsize)+\
						'arcsec\',niter=100000000000,gain=0.1,threshold='+str(threshold)+',deconvolver=\'multiscale\',scales='+str(self.multiscale_scales)+\
						',uvtaper=\''+str(self.uvtaper)+'\',weighting=\'natural\',interactive=False,usemask=\'user\',mask=\'\')\n')
					poltclean(vis=self.msname,imagename=imagename,selectdata=True,startmodel=startmodel,startmask=startmask,stokes=stokes,antenna=antenna_to_use,imsize=[self.imsize],\
					cell=str(self.cellsize)+'arcsec',niter=100000000000,gain=0.1,threshold=threshold,deconvolver='multiscale',scales=self.multiscale_scales,\
					uvtaper=self.uvtaper,weighting='natural',interactive=False,usemask='user',mask='')
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
					+'arcsec\',niter=100000000000,gain=0.05,threshold=\''+str(threshold)+'\',deconvolver=\'multiscale\',scales='+str(self.multiscale_scales)
					+',uvtaper=\''+self.uvtaper+'\',weighting=\'natural\',interactive=False,mask='+str([maskfile])+')\n')	
			poltclean(vis=self.msname,imagename=imagename,selectdata=True,startmodel=startmodel,startmask=startmask,stokes=stokes,antenna=antenna_to_use,imsize=[self.imsize],\
			cell=str(self.cellsize)+'arcsec',niter=100000000000,gain=0.05,threshold=threshold,deconvolver='multiscale',scales=self.multiscale_scales,\
			uvtaper=self.uvtaper,weighting='natural',interactive=False,mask='')
		if do_solarquv_cor==True:
			self.pollog_verbose.info('Correcting solar Stokes I to Q,U leakage based on image.\n')	
			self.correct_solar_quv_leakage(imagename+'.image',imagename+'.model',sigma,overwrite=True)
		out_dict,negative_dyn_range=self.calc_dyn_range(num_iter,sigma,box_width=box_width,stokes_list=['I','Q','U','V']) # Calculating the dynamic range of the image
		out_dict_keys=out_dict.keys()
		if 'NAN' in out_dict_keys:
			self.pollog_verbose.info(B.error_msgs(3))
			os.system('rm -rf casa*log')
			os.chdir(self.cwd)
			os.system('rm -rf casa*log')
			return 3     # If image is not made, no point in continuing
		if os.path.isdir(imagename+'.model')==False:
			self.pollog_verbose.info(B.error_msgs(4))
			os.system('rm -rf casa*log')
			os.chdir(self.cwd)
			os.system('rm -rf casa*log')
			return 4	   # If model is not present no point in continuing
		else:
			modelflux=imstat(imagename=imagename+'.model')['sum'][0]
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
				self.pollog_verbose.info('delmod(vis=\''+self.msname+'\',scr=True)\n') 
				delmod(vis=self.msname,scr=True) # Clear the MODEL column
				self.pollog_verbose.info('ft(vis=\''+self.msname+'\',model=\''+imagename+'.model\',nterms=1,usescratch=True)\n') 
				ft(vis=self.msname,model=imagename+'.model',nterms=1,usescratch=True) # Putting the model into MS
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
				if poldistortion_correction==False:
					self.pollog_verbose.info('cal.applycal(msname=\''+self.msname+'\',gaintable=\''+caltable_name+'\',applymode=\'calflag\',flagbackup=True)\n')
					cal.applycal(msname=self.msname,gaintable=caltable_name,applymode='calflag',flagbackup=True) # Applying the solution
					if use_ankflagger==True:
						try:
							self.pollog_verbose.info('do_uvsub_ankflag(\''+self.msname+'\',nthread=1,verbose='+str(False)+')\n')
							fg.do_uvsub_ankflag(self.msname,nthread=1,verbose=False)
						except Exception as e:
							self.pollog_verbose.info('Error in aNKflagger : '+str(e)+'\n')
							self.pollog_verbose.info('Error in running aNKflagger. Using rms threshold flagging.\n')
							self.pollog_verbose.info('do_uvsub_flagger(\''+self.msname+'\',mode=\'uvsub_flag\',rmsthresh=[10,7,5,3.5])\n')
							fg.do_uvsub_flagger(self.msname,mode='uvsub_flag',rmsthresh=[10,7,5,3.5])
					else:
						self.pollog_verbose.info('do_uvsub_flagger(\''+self.msname+'\',mode=\'uvsub_flag\',rmsthresh=[10,7,5,3.5])\n')
						fg.do_uvsub_flagger(self.msname,mode='uvsub_flag',rmsthresh=[10,7,5,3.5])
					if do_solarquv_cor==True:
						tb=table()
						tb.open(self.msname)
						cor_data=tb.getcol('CORRECTED_DATA')
						tb.close()
						tb.open(self.msname,nomodify=False)
						tb.putcol('DATA',cor_data)
						tb.flush()
						tb.close()
				elif poldistortion_correction==True:
					self.pollog_verbose.info('self.cal_poldistortion(\''+caltable_name+'\',poldistortion_matrix='+poldistortion_matrix+')\n')
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
					self.pollog_verbose.info('Applying poldistoryion corrected caltable ........\n')
					self.pollog_verbose.info('cal.applycal(msname=\''+self.msname+'\',gaintable=\''+caltable_name+'\',applymode=\'calflag\',flagbackup=True)\n')
					cal.applycal(msname=self.msname,gaintable=caltable_name,applymode='calflag',flagbackup=True) # Applying the solution
					self.pollog_verbose.info('do_uvsub_flagger(\''+self.msname+'\',mode=\'uvsub_flag\',rmsthresh=[10,7,5,3.5])\n')
					fg.do_uvsub_flagger(self.msname,mode='uvsub_flag',rmsthresh=[10,7,5,3.5])
					if do_solarquv_cor==True:
						tb=table()
						tb.open(self.msname)
						cor_data=tb.getcol('CORRECTED_DATA')
						tb.close()
						tb.open(self.msname,nomodify=False)
						tb.putcol('DATA',cor_data)
						tb.flush()
						tb.close()
				self.pollog_verbose.info('DR_I:'+str(out_dict['I'][0])+', DR_Q:'+str(out_dict['Q'][0])+', DR_U:'+str(out_dict['U'][0])\
												+', DR_V:'+str(out_dict['V'][0])+', DR_neg:'+str(negative_dyn_range)+'\n')
				self.pollog_verbose.info('Success.\n')
				os.system('rm -rf casa*log')
				os.chdir(self.cwd)
				os.system('rm -rf casa*log')
				return 0,out_dict,negative_dyn_range

#########################################
# Finished PolSelfcal Class
#########################################
