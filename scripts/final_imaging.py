import os,time,psutil
from optparse import OptionParser
os.system('rm -rf casa*log')
cwd=os.getcwd()
start_time=time.time()
usage= ' Perform final imaging\n'
parser = OptionParser(usage=usage)
parser.add_option('--msname',dest="msname",default=None,help="Name of measurement set of a single time anf frequency slice",metavar="Measurement Set")
parser.add_option('--metafits',dest='metafits',default=None,help='Name of the metafits file',metavar='Metafits file')
parser.add_option('--basedir',dest='basedir',default=None,help='Name of the base directory',metavar='Directory path')
parser.add_option('--workdir',dest='workdir',default=None,help='Name of the working directory',metavar='Directory path')
parser.add_option('--savedir',dest='savedir',default=None,help='Directory name to save final images',metavar="Directory path")
parser.add_option('--savemodel',dest='savemodel',default=False,help='Want to save final models',metavar="Boolean")
parser.add_option('--saveres',dest="saveresidual",default=False,help="Want to save residual images",metavar="Boolean")
parser.add_option('--stokes',dest='stokes',default='pseudoI',help='Stokes planes to image',metavar="String")
parser.add_option('--cutoutbox',dest='cutoutbox',default='',help='Cutout box \'X_width,Y_width\' in degree',metavar="Comma separated string")
parser.add_option('--threshold',dest='threshold',default=0.1,help='RMS threshold for cleaning for each Stokes plane',metavar="Comma separated string")
parser.add_option('--sigma',dest='sigma',default=10,help='Sigma value for thresholding',metavar="Float")
parser.add_option('--want_automask',dest='want_automask',default=False,help='Want auto masking or not',metavar="Boolean")
parser.add_option('--maskfile',dest='maskfile',default=None,help='Mask for imaging when auto masking is off',metavar="Maskfile or CASA mask string")
parser.add_option('--quality_factor',dest='quality_factor',default=1,help='Quality factor of imaging',metavar="Integer")
parser.add_option('--use_ankflag',dest='use_ankflag',default=False,help='Use aNKflag for flagging or not',metavar="Boolean")
parser.add_option('--residual_frac',dest='resfrac',default=0.1,help='Residual flux fraction',metavar="Float")
parser.add_option('--casa_caltables',dest='casacals',default='',help='CASA caltables',metavar="Comma separated string")
parser.add_option('--calibrate_caltables',dest='calibratecals',default='',help='CALIBRATE caltables',metavar="Comma separated string")
parser.add_option('--wsclean',dest="use_wsclean",default=True,help="Use WSClean for imaging or not",metavar="Boolean")
parser.add_option('--do_diffcal',dest="do_diffcal",default=False,help="Use WSClean for imaging or not",metavar="Boolean")
parser.add_option('--cpu_frac',dest='cpu_frac',default=0.5,help='Fraction of cpu to use',metavar="Float")
parser.add_option('--imaging_mode',dest='mode',default='final',help='Imaging for database or final imaging',metavar="String")
parser.add_option('--inputfile',dest='inputfile',default=None,help='Path of the P-AIRCARS input file',metavar="File path")
parser.add_option('--major_axis',dest='maj',default=None,help='Final image restoring beam major axis (FWHM) in arcsec',metavar="Float")
parser.add_option('--minor_axis',dest='minor',default=None,help='Final image restoring beam minor axis (FWHM) in arcsec',metavar="Float")
parser.add_option('--pa',dest='pa',default=None,help='Final image restoring beam position angle in degree',metavar="Float")
(options, args) = parser.parse_args()

while True:
	available_cpu=int(psutil.cpu_count()*(1-(psutil.cpu_percent()/100.0))*float(options.cpu_frac))
	if available_cpu>0:
		break
os.environ['OPENBLAS_NUM_THREADS'] = str(available_cpu)
os.environ['OMP_NUM_THREADS'] = str(available_cpu)
from casatasks import *
from casatools import msmetadata,table,measures,quanta,agentflagger,image,calibrater,ms,imager
from paircars.basic_func import *
from paircars.access_ms import *
from paircars_casatasks.poltclean import *
from paircars.flagger import *
from astropy.io import fits
from CALIBRATE.access_calibrate import *
from paircars.fullpol_selfcal_LTS import *
import time
inputfile=str(options.inputfile)
if inputfile[-1]=='/':
	inputfile=inputfile[:-1]
sys.path.append(os.path.dirname(os.path.abspath(inputfile)))
import selfcal_inputs as inputs

def calc_residual_flux(imagename,residual,nsigma,rms_box,stokes_list=['I']):
	'''
	Function to calculate residual flux.
	Parameters:
	imagename = Name of the image
	nsigma = Sigma of present CLEAN cycle
	rms_box = Box to calculate rms
	stokes_list = ['I'], stokes plane list
	Return:
	Median residual flux over the stokes plane
	'''
	imagename_prefix=imagename.split('.image')[0]
	do_reduce_list=[]
	residual_frac_list=[]
	rms_list=[]
	ia=image()
	imagename_path=os.path.dirname(os.path.realpath(imagename))
	cwd=os.getcwd()
	if imagename_path!='':
		os.chdir(imagename_path)
	os.system('rm -rf '+imagename_prefix+'.reduce_sigma_*')
	for stokes in stokes_list:
		if stokes=='I' or stokes=='XX' or stokes=='YY':
			if os.path.isdir(imagename_prefix+'.reduce_sigma_I.image')==True:
				os.system('rm -rf '+imagename_prefix+'.reduce_sigma_I.image')
			if os.path.isdir(imagename_prefix+'.reduce_sigma_I.residual')==True:
				os.system('rm -rf '+imagename_prefix+'.reduce_sigma_I.residual')
			immath(imagename=imagename,mode='evalexpr',stokes=stokes,outfile=imagename_prefix+'.reduce_sigma_I.image')
			immath(imagename=residual,mode='evalexpr',stokes=stokes,outfile=imagename_prefix+'.reduce_sigma_I.residual')
			rms=imstat(imagename=imagename_prefix+'.reduce_sigma_I.image',box=rms_box,stokes=stokes)['rms'][0]
			rms_list.append(rms)
			ia.open(imagename_prefix+'.reduce_sigma_I.image')			
			ia.calcmask('\"'+imagename_prefix+'.reduce_sigma_I.image\">'+str(nsigma*rms),'mymask')
			ia.close()
			makemask(inpimage=imagename_prefix+'.reduce_sigma_I.image',inpmask=imagename_prefix+'.reduce_sigma_I.image:mymask',\
						output=imagename_prefix+'.reduce_sigma_I.residual:mymask',mode='copy')
			try:
				image_pix_sum=imstat(imagename=imagename_prefix+'.reduce_sigma_I.image')['sum'][0]
				residual_pix_sum=imstat(imagename=imagename_prefix+'.reduce_sigma_I.residual')['sum'][0]
			except:
				image_pix_sum=1
				residual_pix_sum=0
		else:
			immath(imagename=imagename,mode='evalexpr',stokes=stokes,expr='abs(IM0)',outfile=imagename_prefix+'.reduce_sigma_'+stokes+'.image')
			immath(imagename=residual,mode='evalexpr',stokes=stokes,outfile=imagename_prefix+'.reduce_sigma_'+stokes+'.residual')
			rms=imstat(imagename=imagename_prefix+'.reduce_sigma_'+stokes+'.image',box=rms_box,stokes=stokes)['rms'][0]
			rms_list.append(rms)
			ia.open(imagename_prefix+'.reduce_sigma_'+stokes+'.image')			
			ia.calcmask('\"'+imagename_prefix+'.reduce_sigma_'+stokes+'.image\">'+str(nsigma*rms),'mymask')
			ia.close()
			makemask(inpimage=imagename_prefix+'.reduce_sigma_'+stokes+'.image',inpmask=imagename_prefix+'.reduce_sigma_'+stokes+'.image:mymask',\
					output=imagename_prefix+'.reduce_sigma_'+stokes+'.residual:mymask',mode='copy')
			try:
				image_pix_sum=imstat(imagename=imagename_prefix+'.reduce_sigma_'+stokes+'.image')['sum'][0]
				residual_pix_sum=imstat(imagename=imagename_prefix+'.reduce_sigma_'+stokes+'.residual')['sum'][0]
			except:
				image_pix_sum=1
				residual_pix_sum=0
		if image_pix_sum==0:
			residual_frac_list.append(0)
		else:
			residual_frac_list.append(residual_pix_sum/image_pix_sum)
	os.system('rm -rf '+imagename_prefix+'.reduce_sigma_*')
	os.chdir(cwd)
	residual_frac_median=np.median(np.array(residual_frac_list))
	return residual_frac_median,rms_list


def modify_header(imagename,astrometry=False,imaging_cmd=''):
	header=fits.getheader(imagename)
	header['PIPELINE']='P-AIRCARS'
	header['WRITTER']='Devojyoti Kansabanik, Surajit Mondal'
	if astrometry==True:
		astrometry=1
	else:
		astrometry=0
	header['ASTRO']=astrometry
	if imaging_cmd!='':
		header['HISTORY']=imaging_cmd
	fits.writeto(imagename,data=fits.getdata(imagename),header=header,overwrite=True)
	return 

def export_images(imagename,OBSID,cell,imsize,imaging_cmd='',savedir='',savemodel=False,saveresidual=False,cutoutbox=[],astrometry=False): 
	'''
	Function to save final images
	imagename = Name of the image to export
	OBSID = OBSID of the observation	
	cell = Cell size of the image in arcsecond
	imsize = Number of pixels in image
	savedir = Directory to save the final images
	savemodel = Whether save model images or not
	saveresidual = Whether save residual images or not
	cutoutbox = [] , cutout box for final image [x_width,y_width] in degree
	astrometry = Astrometry corrected or not
	Return :
	List of image,model,residual
	'''
	cwd=os.getcwd()
	if os.path.isdir(savedir+'/All_final_images/'+str(OBSID))==False:	
		os.makedirs(savedir+'/All_final_images/'+str(OBSID))
	imagedir=savedir+'/All_final_images/'+str(OBSID)
	if savemodel==True:
		if os.path.isdir(savedir+'/All_final_models/'+str(OBSID))==False:
			os.makedirs(savedir+'/All_final_models/'+str(OBSID))
	modeldir=savedir+'/All_final_models/'+str(OBSID)
	if saveresidual==True:
		if os.path.isdir(savedir+'/All_final_residuals/'+str(OBSID))==False:
			os.makedirs(savedir+'/All_final_residuals/'+str(OBSID))
	resdir=savedir+'/All_final_residuals/'+str(OBSID)
	output=[]
	if len(cutoutbox)!=0:
		x_pix=int((float(cutoutbox[0])*3600)/float(cell))
		y_pix=int((float(cutoutbox[1])*3600)/float(cell))
		x_cen=int(imsize/2)
		y_cen=x_cen
		box=str(int(x_cen-x_pix/2))+','+str(int(y_cen-y_pix/2))+','+str(int(x_cen+x_pix/2))+','+str(int(y_cen+y_pix/2))
		os.system('rm -rf '+imagename+'.cutout*')
		if os.path.isdir(imagedir):		
			os.system('rm -rf '+imagedir+'/'+os.path.basename(imagename)+'_*.fits')
		if os.path.isdir(modeldir):
			os.system('rm -rf '+modeldir+'/'+os.path.basename(imagename)+'_*.fits')
		if os.path.isdir(resdir):
			os.system('rm -rf '+resdir+'/'+os.path.basename(imagename)+'_*.fits')
		imsubimage(imagename=imagename+'.image',outfile=imagename+'.cutout.image',box=box)
		exportfits(imagename=imagename+'.cutout.image',fitsimage=imagedir+'/'+os.path.basename(imagename)+'_image.fits',history=False)
		output.append(imagedir+'/'+os.path.basename(imagename)+'_image.fits')
		modify_header(imagedir+'/'+os.path.basename(imagename)+'_image.fits',astrometry=astrometry,imaging_cmd=imaging_cmd)
		if savemodel==True:
			imsubimage(imagename=imagename+'.model',outfile=imagename+'.cutout.model',box=box)
			exportfits(imagename=imagename+'.cutout.model',fitsimage=modeldir+'/'+os.path.basename(imagename)+'_model.fits',history=False)
			output.append(modeldir+'/'+os.path.basename(imagename)+'_model.fits')
			modify_header(modeldir+'/'+os.path.basename(imagename)+'_model.fits',astrometry=astrometry,imaging_cmd=imaging_cmd)
		if saveresidual==True:
			imsubimage(imagename=imagename+'.residual',outfile=imagename+'.cutout.residual',box=box)
			exportfits(imagename=imagename+'.cutout.residual',fitsimage=resdir+'/'+os.path.basename(imagename)+'_res.fits',history=False)
			output.append(resdir+'/'+os.path.basename(imagename)+'_res.fits')
			modify_header(resdir+'/'+os.path.basename(imagename)+'_res.fits',astrometry=astrometry,imaging_cmd=imaging_cmd)
		os.system('rm -rf '+imagename+'.cutout*')
	else:
		if os.path.isdir(imagedir):
			os.system('rm -rf '+imagedir+'/'+os.path.basename(imagename)+'_*.fits')
		if os.path.isdir(modeldir):		
			os.system('rm -rf '+modeldir+'/'+os.path.basename(imagename)+'_*.fits')
		if os.path.isdir(resdir):
			os.system('rm -rf '+resdir+'/'+os.path.basename(imagename)+'_*.fits')
		exportfits(imagename=imagename+'.image',fitsimage=imagedir+'/'+os.path.basename(imagename)+'_image.fits',history=False)
		output.append(imagedir+'/'+os.path.basename(imagename)+'_image.fits')
		modify_header(imagedir+'/'+os.path.basename(imagename)+'_image.fits',astrometry=astrometry,imaging_cmd=imaging_cmd)
		if savemodel==True:
			exportfits(imagename=imagename+'.model',fitsimage=modeldir+'/'+os.path.basename(imagename)+'_model.fits',history=False)
			output.append(modeldir+'/'+os.path.basename(imagename)+'_model.fits')
			modify_header(modeldir+'/'+os.path.basename(imagename)+'_model.fits',astrometry=astrometry,imaging_cmd=imaging_cmd)
		if saveresidual==True:
			exportfits(imagename=imagename+'.residual',fitsimage=resdir+'/'+os.path.basename(imagename)+'_res.fits',history=False)
			output.append(resdir+'/'+os.path.basename(imagename)+'_res.fits')
			modify_header(resdir+'/'+os.path.basename(imagename)+'_res.fits',astrometry=astrometry,imaging_cmd=imaging_cmd)
	os.chdir(imagedir)
	os.system('rm -rf *.pb *.mask *.model *.image *.flux *.sumwt *.residual *psf*')
	os.chdir(cwd)
	return output

def get_stokes(stokes):
	if stokes=='I':
		return ['I']
	elif stokes=='Q':
		return ['Q']
	elif stokes=='U':
		return ['U']
	elif stokes=='V':
		return ['V']
	elif stokes=='IV':
		return ['I','V']
	elif stokes=='QU':
		return ['Q','U']
	elif stokes=='IQ':
		return ['I','Q']
	elif stokes=='UV':
		rerun ['U','V']
	elif stokes=='IQUV':
		return ['I','Q','U','V']
	elif stokes=='RR':
		return ['RR']
	elif stokes=='LL':
		return ['LL']
	elif stokes=='XX':
		return ['XX']
	elif stokes=='YY':
		return ['YY']
	elif stokes=='RRLL':
		return ['RR','LL']
	elif stokes=='XXYY':
		return ['XX','YY']
	else:
		return


def make_image(msname,metafits,workdir,sigma=10,stokes='I',savedir='',want_automask=False,maskfile='',quality_factor=1,threshold=[0.1],do_diffcal=False,\
				savemodel=False,saveresidual=False,cutoutbox='',use_ankflagger=False,residual_frac=0.1,use_wsclean=True,cpus=3,absmem=5,\
				clean_beam=[]): #TODO : Wide FOV Differential beam correction
	'''
	Function to make final images
	Parameters:		
	msname = Name of the measurement set
	metafits = Name of the metafits file
	workdir = Name of the working directory
	sigma = Sigma value for thresholding
	stokes = Stokes planes to image
	savedir = Directory to save final directory
	threshold = [], rms list for all Stokes plane for thresholding
	do_diffcal = Perform differential calibration
	want_automask = Want to use auto masking
	maskfile = Name of any previous mask file or CASA mask string
	uvtaper = UV-taper for imaging
	quality_factor = Quality factor for imaging
	savemodel = Save model images or not
	saveresidual = Save residual images or not
	cutoutbox = Cutout box of the final image [x_width,y_width] in degree (default : [], no cutout, save full image)
	use_ankflagger = Whether use aNKflagger for flagging or not
	residual_frac = Residual flux fraction to stop CLEANing
	use_wsclean = True, use wsclean or not
	cpus = Number of cpus to use in wsclean
	absmem = Absolute memory in GB for wsclean
	clean_beam = [], clean beam [maj,min,pa]
	Result:
	Name of final image,model,residual
	'''
	print ('Making image.....\n')
	if msname[-1]=='/':
		msname=msname[:-1]
	os.system('cp -r '+msname+' '+msname+'.backup')
	backup_ms=msname+'.backup'
	wsclean_check_count=0
	if use_wsclean==True:
		a=os.system('wsclean > wsclean_test')
		while a!=0:
			if wsclean_check_count>5:
				print('WSClean is not installed. Using CASA for imaging.\n')
				use_wsclean=False
			else:
				time.sleep(1.0)
		os.system('rm -rf wsclean_test')
		tempdir=workdir+'/tempdir'
		if os.path.isdir(tempdir)==False:
			os.makedirs(tempdir)
	casalog=False
	cwd=os.getcwd()
	os.chdir(workdir)
	if cutoutbox!='':
		cutoutbox=cutoutbox.split(',')
	else:
		cutoutbox=[]
	if stokes!='I' or stokes!='XXYY':
		cpus=5
	else:
		cpus=3
	AM=AccessMS(msname)
	try_count=0
	done_imaging=False
	while done_imaging==False:
		try:
			freqs=AM.calc_meanfreq()/10**6	
			OBSID=str(fits.getheader(metafits)['GPSTIME'])	
			file_str=workdir+'/'+os.path.basename(splited_ms_rename(msname,ref_time_chan=False,change_msname=False)).split('.ms')[0]
			if quality_factor==0:
				gain=0.2
				mgain=0.8
				multiscale_gain=0.12
			elif quality_factor==1:
				gain=0.1
				mgain=0.75
				multiscale_gain=0.1
			else:
				gain=0.08	
				mgain=0.7
				multiscale_gain=0.07			
			stokes_list=get_stokes(stokes)
			mask=file_str+'.mask'
			PSC=PolSelfcal(options.msname,options.metafits,32*60,verbose=False,interactive=False,use_wsclean=use_wsclean,savelog=False) 
			cal=CALIBRATE()
			IB=ImageBasic(msname)
			if inputs.quality_factor==0 and inputs.safety_factor==0:
				num_pixel_in_psf=3
			else:
				num_pixel_in_psf=5
			if inputs.calc_image_parameters==True:
				calc_calib_uvrange=IB.calc_calib_uvrange(12,includeflag=False)
				uvrange=calc_calib_uvrange[0]	
				imaging_minuv=calc_calib_uvrange[3]
				imaging_maxuv=calc_calib_uvrange[4]	
				cell=IB.calc_cellsize(num_pixel_in_psf) 
				imsize=IB.num_pixels(num_pixel_in_psf)
				scales=IB.choose_scales(num_pixel_in_psf,32*60)
				uvtaper=IB.calc_uvtaper()
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
				uvrange_to_cal=inputs.uvrange_to_cal
				uvrange=uvrange_to_cal
				if uvrange_to_cal!='':
					if '~' in uvrange_to_cal:
						uvrange=uvrange_to_cal.split('~')
						imaging_minuv=float(uvrange[0])
						imaging_maxuv=float(uvrange[1])
					elif '>' in uvrange_to_cal:
						uvrange=uvrange_to_cal.split('>')
						imaging_minuv=float(uvrange[1])
						imaging_maxuv=AM.get_max_baseline()/AM.calc_meanwavelength()
					elif '<' in uvrange_to_cal:
						uvrange=uvrange_to_cal.split('<')
						imaging_minuv=0
						imaging_maxuv=float(uvrange[1])
				else:
					calc_calib_uvrange=IB.calc_calib_uvrange(12,includeflag=False)
					uvrange=calc_calib_uvrange[0]	
					imaging_minuv=calc_calib_uvrange[3]
					imaging_maxuv=calc_calib_uvrange[4]	
				cell=inputs.cellsize
				imsize=inputs.imsize[0]
				scales=inputs.multiscale_scales		
				if inputs.uvtaper!='':
					uvtaper=inputs.uvtaper
				else:
					uvtaper=IB.calc_uvtaper()
				weight='briggs'
				robust=1.0
			rms_box='50,50,'+str(imsize-50)+','+str(int(imsize/4))
			os.system('rm -rf '+file_str+'.image '+file_str+'.model '+file_str+'.mask '+file_str+'.psf '+file_str+'.pb '+file_str+'.sumwt '+file_str+'.flux '+file_str+'.residual '+\
							file_str+'*.fits')
			count=0
			clean_continue=False
			while True:
				AM=AccessMS(msname)
				flagfrac=AM.calc_flagfrac()
				if flagfrac<=0.95:
					try:
						if use_wsclean==False:
							casa_imagename=file_str+'.image'
							casa_residualname=file_str+'.residual'
							casa_modelname=file_str+'.model'
							if len(threshold)!=len(stokes_list):	
								threshold=[threshold[0]]*len(stokes_list)
							if maskfile=='':
								mask_rad=int((60*60)/float(cell)) # Creating a mask with 60 arcmin radius centered on the image
								mask_str='circle[['+str(imsize/2)+'pix,'+str(imsize/2)+'pix],'+str(mask_rad)+'pix]'
							else:
								mask_str=maskfile
							threshold_list=[str(rms*5)+'Jy' for rms in threshold]
							if len(clean_beam)==3:
								restoring_beam=[str(clean_beam[0])+'arcsec',str(clean_beam[1])+'arcsec',str(clean_beam[2])+'deg']
							else:
								restoring_beam=[]
							if os.path.isdir(mask)==True:
								imaging_cmd='poltclean(vis=\''+msname+'\',selectdata=True,datacolumn="corrected",uvrange=\''+str(uvrange)+'\',imagename=\''+file_str+\
								'\',imsize=['+str(imsize)+'],cell=\''+str(cell)+'\',stokes=\''+str(stokes)+'\',gridder=\'standard\',pblimit=-1,deconvolver="multiscale",scales='+\
								str(scales)+',nterms=1,weighting=\''+weight+'\',robust='+str(robust)+',uvtaper=\''+str(uvtaper)+'\',casalogger='+\
								str(casalog)+',niter=100000000000,gain='+str(gain)+',threshold='+str(threshold_list)+',interactive=False,usemask="user",startmask='+\
								str(mask)+',savemodel=\'modelcolumn\',restoringbeam='+str(restoring_beam)+')'
								print (imaging_cmd+'\n')
								poltclean(vis=msname,selectdata=True,datacolumn="corrected",uvrange=uvrange,imagename=file_str,imsize=[imsize],cell=cell,stokes=stokes,gridder='standard',\
								pblimit=-1,deconvolver="multiscale",scales=scales,nterms=1,weighting=weight,robust=robust,uvtaper=uvtaper,casalogger=casalog,\
								niter=100000000000,gain=gain,threshold=threshold_list,interactive=False,usemask="user",startmask=mask,savemodel='modelcolumn',restoringbeam=restoring_beam)
							elif want_automask==True and maskfile=='':
								imaging_cmd='poltclean(vis=\''+msname+'\',selectdata=True,datacolumn="corrected",imagename=\''+file_str+'\',imsize=['+str(imsize)+'],cell=\''+str(cell)+\
								'\',stokes=\''+stokes+'\',gridder=\'standard\',pblimit=-1,deconvolver="multiscale",scales='+str(scales)+',nterms=1,weighting=\''+weight+'\',robust='+\
								str(robust)+'uvtaper='+str(uvtaper)+',niter=100000000000,gain='+str(gain)+',casalogger='+str(casalog)+',threshold='+str(threshold_list)+\
								',interactive=False,usemask=\'auto-multithresh\',negativethreshold=3.0,savemodel=\'modelcolumn\',restoringbeam='+str(restoring_beam)+')'
								print (imaging_cmd+'\n')
								poltclean(vis=msname,selectdata=True,datacolumn="corrected",uvrange=uvrange,imagename=file_str,imsize=[imsize],cell=cell,stokes=stokes,gridder='standard',\
								pblimit=-1,deconvolver="multiscale",scales=scales,nterms=1,weighting=weight,robust=robust,uvtaper=uvtaper,niter=100000000000,gain=gain,casalogger=casalog,\
								threshold=threshold_list,interactive=False,usemask='auto-multithresh',negativethreshold=3.0,savemodel='modelcolumn',restoringbeam=restoring_beam)
							else:
								imaging_cmd='poltclean(vis=\''+msname+'\',selectdata=True,datacolumn="corrected",imagename=\''+str(file_str)+'\',imsize=['+str(imsize)+'],cell=\''+\
								str(cell)+'\',stokes=\''+str(stokes)+'\',gridder=\'standard\',pblimit=-1,deconvolver="multiscale",scales='+str(scales)+\
								',nterms=1,weighting=\''+weight+'\',robust='+str(robust)+',uvtaper=\''+str(uvtaper)+'\',casalogger='+str(casalog)+\
								',niter=100000000000,threshold='+str(threshold_list)+',interactive=False,usemask="user",mask=\''+str(mask_str)+\
								',savemodel=\'modelcolumn\',restoringbeam='+str(restoring_beam)+')'
								print (imaging_cmd+'\n')
								poltclean(vis=msname,selectdata=True,datacolumn="corrected",uvrange=uvrange,imagename=file_str,imsize=[imsize],cell=cell,stokes=stokes,gridder='standard',\
								pblimit=-1,deconvolver="multiscale",scales=scales,nterms=1,weighting=weight,robust=robust,uvtaper=uvtaper,casalogger=casalog,\
								niter=100000000000,threshold=threshold_list,interactive=False,usemask="user",mask=mask_str,savemodel='modelcolumn',restoringbeam=restoring_beam)
						else:
							if stokes=='I':
								pol='i'
							elif stokes=='RR':
								pol='rr'
							elif stokes=='LL':
								pol='ll'
							elif stokes=='XX':
								pol='xx'
							elif stokes=='YY':
								pol='yy'
							elif stokes=='RRLL':
								pol='rrll'
							elif stokes=='XXYY':
								pol='xxyy'
							elif stokes=='IQUV':
								pol='iquv'
							scales=[str(i) for i in scales]
							if weight=='briggs':
								weight=weight+' '+str(robust)
							wsclean_args=['-scale '+str(cell)+'asec','-size '+str(imsize)+' '+str(imsize),'-no-dirty','-j '+str(cpus),'-abs-mem '+str(absmem),'-weight '+weight,\
							'-taper-tukey 10','-name '+str(file_str),'-nwlayers 1','-pol '+str(pol),'-maxuv-l '+str(imaging_maxuv),'-minuv-l '+str(imaging_minuv),'-niter 1000000',\
							'-mgain '+str(mgain),'-quiet','-auto-threshold '+str(sigma),'-auto-mask '+str(sigma+0.1),'-multiscale-gain '+str(multiscale_gain),\
							'-multiscale-scale-bias 0.7','-gain '+str(gain),'-multiscale','-multiscale-scales '+','.join(scales),'-temp-dir '+tempdir]
							if clean_continue==True:
								wsclean_args.append('-continue')
							if len(clean_beam)==3:
								wsclean_args.append('-beam-shape '+str(clean_beam[0])+' '+str(clean_beam[1])+' '+str(clean_beam[2]))
							imaging_cmd='wsclean '+' '.join(wsclean_args)+' '+msname
							print (imaging_cmd+'\n')			
							os.system('wsclean '+' '.join(wsclean_args)+' '+msname)
							time.sleep(2)
							wsclean_images=sorted(glob.glob(file_str+'*image.fits'))
							wsclean_models=sorted(glob.glob(file_str+'*model.fits'))
							wsclean_residuals=sorted(glob.glob(file_str+'*residual.fits'))
							print ('Converting WSClean images : '+','.join([os.path.basename(i) for i in wsclean_images])+' to CASA image : '+file_str+'.image\n')
							casa_imagename=PSC.wsclean_to_casaimage(wsclean_images=wsclean_images,casaimage_prefix=file_str,imagetype='image',keep_wsclean_images=True)
							print ('Converting WSClean models : '+','.join([os.path.basename(i) for i in wsclean_models])+' to CASA model : '+file_str+'.model\n')
							casa_modelname=PSC.wsclean_to_casaimage(wsclean_images=wsclean_models,casaimage_prefix=file_str,imagetype='model',keep_wsclean_images=True)
							print ('Converting WSClean residuals : '+','.join([os.path.basename(i) for i in wsclean_residuals])+' to CASA residual : '+file_str+'.residual\n')
							casa_residualname=PSC.wsclean_to_casaimage(wsclean_images=wsclean_residuals,casaimage_prefix=file_str,imagetype='residual',keep_wsclean_images=True)
						if count==0:
							qucor_image,qucor_model,qchange,uchange,vchange=PSC.correct_solar_quv_leakage(casa_imagename,casa_modelname,sigma,overwrite=False)
							qucor_image,qucor_model=PSC.pol_model_threshold(qucor_image,qucor_model,sigma,1)
							if qchange>=0.01 or uchange>=0.01 or vchange>=0.01 or do_diffcal==True:
								casa_modelname=qucor_model
								if do_diffcal==True:
									cal_cause='Differential calibration'
								elif qchange>=0.01 or uchange>=0.01 or vchange>=0.01:
									cal_cause='Stokes leakage decreases.'
								print ('Performing differential calibration because : '+cal_cause+'\n')
								print ('delmod(vis=\''+msname+'\',scr=True,otf=True)\n')
								delmod(vis=msname,scr=True,otf=True)
								k=0
								while True:
									if k>3:
										print ('Model import tried three times. Model import failed.\n')
										return 2
									print ('ft(vis=\''+msname+'\',model=\''+casa_modelname+'\',usescratch=True)\n')
									ft(vis=msname,model=casa_modelname,usescratch=True)
									AMmodel=AccessMS(msname)
									model_imported=AM.model_imported()
									if model_imported==True:
										return 2
									else:
										print ('Try to import model. Trial number : '+str(k)+'\n')
										time.sleep(1.0)
									k+=1
								if stokes!='I':
									print ('cal.calibrate(msname=\''+msname+'\',caltable=\''+file_str+'.cal\',verbose=False,j=3,'+\
									',datacolumn=\'data\')\n')
									cal.calibrate(msname=msname,caltable=file_str+'.cal',verbose=False,j=3,datacolumn='data')
									print ('cal.applycal(msname=\''+msname+'\',gaintable=\''+file_str+'.cal\',datacolumn=\'data\')\n')
									cal.applycal(msname=msname,gaintable=file_str+'.cal',datacolumn='data')
								else:
									print ('gaincal(vis=\''+msname+'\',caltable=\''+file_str+'.cal\',solmode=\'R\',rmsthresh=[10,7,5,3.5])\n')
									gaincal(vis=msname,caltable=file_str+'.cal',solmode='R',rmsthresh=[10,7,5,3.5])
									print ('applycal(vis=\''+msname+'\',gaintable=[\''+file_str+'.cal\'],applymode=\'calflag\')\n')
									applycal(vis=msname,gaintable=[file_str+'.cal'],applymode='calflag')
							else:
								print ('delmod(vis=\''+msname+'\',scr=True,otf=True)\n')
								delmod(vis=msname,scr=True,otf=True)
								print ('ft(vis=\''+msname+'\',model=\''+casa_modelname+'\',usescratch=True)\n')
								ft(vis=msname,model=casa_modelname,usescratch=True)
							if inputs.use_ankflagger:
								do_uvsub_ankflag(msname,model=casa_modelname,verbose=False,nthread=3,extendpols=False,chantime_minfrac=0.8,casaflag='tfcrop')
							else:
								do_uvsub_flagger(msname,model=casa_modelname,rmsthresh=[10,7,5])
							count+=1
							continue
						median_res_frac,threshold=calc_residual_flux(casa_imagename,casa_residualname,sigma,rms_box,stokes_list=stokes_list)
						if median_res_frac>=residual_frac and sigma>=inputs.min_sigma and count>=1:
							print ('Continuing CLEANing, since residual fraction is more than '+str(residual_frac*100)+'%\n')
							clean_continue=True
							sigma-=inputs.sigma_step
							count+=1
							if quality_factor==0 and count>=1:
								done_imaging=True
								break
							elif quality_factor==0 and count>=5:
								done_imaging=True
								break
							elif quality_factor==2 and count>=9:
								done_imaging=True
								break
							else:
								continue
						else:
							done_imaging=True
							break
					except Exception as e:
						print ('Error occured in final imaging : '+str(e)+'\n')
						if try_count>=2:
							print ('Maximum 2 tries failed. Quit.\n')
							done_imaging=True
							return 2
						else:
							antenna=AM.get_antenna_string()
							flagdata(vis=msname,mode='unflag',antenna=antenna)
							flag_MWA_coarse(msname,edgewidth=280,do_flag=True,force=True,flagbackup=False)
							try_count+=1
							print ('Trying count " '+str(try_count)+'\n')
							continue
			else:
				done_imaging=True
				return 1
		except Exception as e:
			print ('Error occured in final imaging : '+str(e)+'\n')
			if try_count>=2:
				print ('Maximum 2 tries failed. Quit.\n')
				done_imaging=True
				return 2
			else:
				os.system('rm -rf '+msname+' '+msname+'.flagversions')
				os.system('cp -r '+backup_ms+' '+msname)
				try_count+=1
				print ('Trying count : '+str(try_count)+'\n')
				continue
			
	# Exporting images
	##################		
	output=export_images(file_str,OBSID,cell,imsize,imaging_cmd=imaging_cmd,savedir=savedir,savemodel=savemodel,saveresidual=saveresidual,cutoutbox=cutoutbox,astrometry=False)
	#os.system('rm -rf *.image *.model *pb *psf* *.sumwt *.flux *.fits')
	#os.system('rm -rf '+file_str+'*')	
	os.system('cd ../')	
	#if savedir!=workdir:
	#	os.system('rm -rf '+workdir)
	os.chdir(cwd)		
	return 0

try:
	if str(options.msname)[-1]=='/':
		msname=str(options.msname)[:-1]
	else:
		msname=str(options.msname)

	if os.path.isdir(str(options.msname))==False or options.msname==None:
		print ('Measurement set is not present.\n')
		os.system('touch '+cwd+'/.Finished_final_imaging_'+str(options.mode)+'_'+os.path.basename(str(msname))+'_noms')
		os._exit(1)
	elif os.path.isfile(str(options.metafits))==False or options.metafits==None:
		print ('Metafits file is not present.\n')
		os.system('touch '+cwd+'/.Finished_final_imaging_'+str(options.mode)+'_'+os.path.basename(str(msname))+'_nometa')
		os._exit(1)
	else:
		if os.path.isdir(options.workdir+'/imaging_data')==False:
			os.makedirs(options.workdir+'/imaging_data')
		if str(options.casacals)=='':
			casacaltables=[]
		else:
			casacaltables=str(options.casacals).split(',')
			if str(options.calibratecals)=='':
				calibratecaltables=[]
			else:
				calibratecaltables=str(options.calibratecals).split(',')
			print ('Applying gain calibration.....\n')
			print ('clearcal(vis=\''+options.msname+'\')\n')
			clearcal(vis=options.msname)
			if len(casacaltables)!=0 or len(calibratecaltables)!=0:
				cal=CALIBRATE()
				casacals=[]
				calibratecals=[]
				for i in casacaltables:
					os.system('cp -r '+i+' '+options.workdir+'/imaging_data/')
					casacals.append(options.workdir+'/imaging_data/'+os.path.basename(i))
				for i in calibratecaltables:
					os.system('cp -r '+i+' '+options.workdir+'/imaging_data/')
					calibratecals.append(options.workdir+'/imaging_data/'+os.path.basename(i))
				print ('applycal(vis=\''+options.msname+'\',gaintable='+str(casacals)+',applymode=\'calflag\',calwt=[False],flagbackup=False)\n')
				a=applycal(vis=options.msname,gaintable=casacals,applymode='calflag',calwt=[False],flagbackup=False)
				
				# Applying cross-hand phase correction
				######################################
				mwa_config=get_MWA_phase(options.metafits) # TODO : Include from cross phase cal solutions
				PSC=PolSelfcal(options.msname,options.metafits,32*60,verbose=False,interactive=False,savelog=False) 
				if mwa_config=='MWAPhaseI':
					crossphase=15
				elif mwa_config=='MWAPhaseIILB' or mwa_config=='MWAPhaseIICOMPACT':
					crossphase=135
				print('Applying cross hand phase solution. Cross hand phase : '+str(crossphase)+' deg.\n')
				PSC.apply_cross_hand_phase(cross_phase=crossphase,caltable='',polbasis='Linear',modify_datacolumn=False,datacolumn='CORRECTED')
				print ('Applying ideal beam correction..\n')
				PSC.correct_visibility_single_beam_jones(datacolumn='CORRECTED_DATA',modify_datacolumn=False,force=True,skip_freq=1.28)
				for i in calibratecals:
					if '.beam' not in i: 
						print ('cal.applycal(msname=\''+options.msname+'\',gaintable=\''+str(i)+'\',applymode=\'calflag\',flagbackup=False,verbose=False)\n')
						cal.applycal(msname=options.msname,gaintable=i,applymode='calflag',flagbackup=False,verbose=False,datacolumn='corrected')
				casa_autoflag(options.msname,mode='rflag',datacolumn='corrected',sigma_thresh=10.0,flagbackup=False,verbose=False)
			os.system('rm -rf '+options.workdir+'/imaging_data')

		print ('#############################\nStart imaging for ms : '+str(msname)+'\n')

		if options.basedir==None:
			basedir=os.path.dirname(os.path.abspath(str(options.msname)))+'/final_imaging_basedir'
		else:
			basedir=str(options.basedir)
		print ('Base directory : '+basedir+'\n')

		if options.workdir==None:
			workdir=os.path.dirname(os.path.abspath(str(options.msname)))+'/final_imaging_workdir'
		else:
			workdir=str(options.workdir)
			
		print ('Working directory : '+workdir+'\n')	
		if options.savedir==None:
			savedir=basedir
			print ('Save directory is not given. Saving final images in base directory : '+basedir+'\n')
		else:
			savedir=str(options.savedir)

		if os.path.isdir(basedir)==False:
			os.makedirs(basedir)
		if os.path.isdir(workdir)==False:
			os.makedirs(workdir)
		if os.path.isdir(savedir)==False:
			os.makedirs(savedir)
		
		error_files=glob.glob(basedir+'/.Finished_final_imaging_'+str(options.mode)+'_*error*')
		for i in error_files:
			os.system('rm -rf '+i)
		touch_file=basedir+'/.Finished_final_imaging_'+str(options.mode)+'_'+os.path.basename(msname)+'_success'
		touch_file_moreflag=basedir+'/.Finished_final_imaging_'+str(options.mode)+'_'+os.path.basename(msname)+'_moreflag'
		if os.path.exists(touch_file) or os.path.exists(touch_file_moreflag):
			print ('Imaging is already done or attempted.\n#############################\n')
			os._exit(0)

		rmsthresh=str(options.threshold).split(',')
		threshold=[float(i) for i in rmsthresh]
		IB=ImageBasic(str(options.msname))
		AM=AccessMS(str(options.msname))
		cent_freq=AM.calc_meanfreq()/10**6
		coarse_chan_freq=freq_to_MWA_coarse(cent_freq)*1.28

		if options.maskfile==None:
			maskfile=''
		else:
			maskfile=str(options.maskfile)

		if options.maj!=None and options.minor!=None and options.pa!=None:
			clean_beam=[float(options.maj),float(options.minor),float(options.pa)]
		else:
			clean_beam=[]
		output=make_image(str(options.msname),str(options.metafits),str(workdir),sigma=float(options.sigma),stokes=str(options.stokes),savedir=str(savedir),threshold=threshold,\
				want_automask=eval(str(options.want_automask)),maskfile=maskfile,quality_factor=int(options.quality_factor),savemodel=eval(str(options.savemodel)),\
				saveresidual=eval(str(options.saveresidual)),cutoutbox=str(options.cutoutbox),use_ankflagger=eval(str(options.use_ankflag)),residual_frac=float(options.resfrac),\
				use_wsclean=eval(str(options.use_wsclean)),do_diffcal=eval(str(options.do_diffcal)),clean_beam=clean_beam)
		if output==0:
			print ('\nImaging finished.\n#############################\n')
			os.system('touch '+touch_file)
		elif output==1:
			print ('More than 95% data are flagged. Image is not made.\n############################\n')
			touch_file=basedir+'/.Finished_final_imaging_'+str(options.mode)+'_'+os.path.basename(msname)+'_moreflag'
			os.system('touch '+touch_file)
		elif output==2:
			print ('Error occured during final imaging.\n############################\n')
			touch_file=basedir+'/.Finished_final_imaging_'+str(options.mode)+'_'+os.path.basename(msname)+'_error'
			os.system('touch '+touch_file)
	#	os.system('rm -rf '+str(msname)+' '+str(msname)+'*')
		print ('Total run time : '+str(time.time()-start_time)+' s\n#############################\n')
except Exception as e:
	print ('Error occured during final imaging. Error : '+str(e)+'\n############################\n')
	touch_file=basedir+'/.Finished_final_imaging_'+str(options.mode)+'_'+os.path.basename(msname)+'_error'
	os.system('touch '+touch_file)
	#os.system('rm -rf '+str(msname)+' '+str(msname)+'*')
	print ('Total run time : '+str(time.time()-start_time)+' s\n#############################\n')







