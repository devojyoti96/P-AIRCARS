from casatasks import *
from casatools import *
from paircars.basic_func import *
from paircars.access_ms import *
from paircars_casatasks.poltclean import *
from optparse import OptionParser
from paircars.flagger import *
from astropy.io import fits
from CALIBRATE.access_calibrate import *
from paircars.fullpol_selfcal_LTS import *
import time


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
				image_pix_sum=0
				residual_pix_sum=1
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
				savemodel=False,saveresidual=False,cutoutbox='',use_ankflagger=False,residual_frac=0.1,use_wsclean=True,cpus=3,absmem=5): #TODO : Wide FOV Differential beam correction
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
	Result:
	Name of final image,model,residual
	'''
	print ('Making image.....\n')
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
	freqs=AM.calc_meanfreq()/10**6	
	OBSID=str(fits.getheader(metafits)['GPSTIME'])	
	file_str=workdir+'/'+os.path.basename(splited_ms_rename(msname,ref_time_chan=False,change_msname=False)).split('.ms')[0]
	if quality_factor==0:
		gain=0.3
	elif quality_factor==1:
		gain=0.15
	else:
		gain=0.1				
	stokes_list=get_stokes(stokes)
	mask=file_str+'.mask'
	PSC=PolSelfcal(options.msname,options.metafits,32*60,verbose=False,interactive=False,use_wsclean=use_wsclean,savelog=False) 
	cal=CALIBRATE()
	IB=ImageBasic(msname)
	uvrange=IB.calc_calib_uvrange(12)[0]	
	imaging_minuv=IB.calc_calib_uvrange(12)[3]
	imaging_maxuv=IB.calc_calib_uvrange(12)[4]
	cell=IB.calc_cellsize(3) # Assuming 3 pixels in one PSF
	imsize=IB.num_pixels(3)
	scales=IB.choose_scales(3,32*60)
	uvtaper=IB.calc_uvtaper()
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
					if os.path.isdir(mask)==True:
						imaging_cmd='poltclean(vis=\''+msname+'\',selectdata=True,datacolumn="corrected",uvrange=\''+str(uvrange)+'\',imagename=\''+file_str+'\',imsize=['+str(imsize)+\
						'],cell=\''+str(cell)+'\',stokes=\''+str(stokes)+'\',gridder=\'standard\',pblimit=-1,deconvolver="multiscale",scales='+str(scales)\
						+',nterms=1,weighting="briggs",robust=1.0,uvtaper=\''+str(uvtaper)+'\',casalogger='+str(casalog)+',niter=100000000000,gain='+str(gain)+',threshold='+\
						str(threshold_list)+',interactive=False,usemask="user",startmask='+str(mask)+',savemodel=\'modelcolumn\')'
						print (imaging_cmd+'\n')
						poltclean(vis=msname,selectdata=True,datacolumn="corrected",uvrange=uvrange,imagename=file_str,imsize=[imsize],cell=cell,stokes=stokes,gridder='standard',\
						pblimit=-1,deconvolver="multiscale",scales=scales,nterms=1,weighting="briggs",robust=1.0,uvtaper=uvtaper,casalogger=casalog,\
						niter=100000000000,gain=gain,threshold=threshold_list,interactive=False,usemask="user",startmask=mask,savemodel='modelcolumn')
					elif want_automask==True and maskfile=='':
						imaging_cmd='poltclean(vis=\''+msname+'\',selectdata=True,datacolumn="corrected",imagename=\''+file_str+'\',imsize=['+str(imsize)+'],cell=\''+str(cell)+\
						'\',stokes=\''+stokes+'\',gridder=\'standard\',pblimit=-1,deconvolver="multiscale",scales='+str(scales)+',nterms=1,weighting="natural",uvtaper='+\
						str(uvtaper)+',niter=100000000000,gain='+str(gain)+',casalogger='+str(casalog)+',threshold='+str(threshold_list)+\
						',interactive=False,usemask=\'auto-multithresh\',negativethreshold=3.0,savemodel=\'modelcolumn\')'
						print (imaging_cmd+'\n')
						poltclean(vis=msname,selectdata=True,datacolumn="corrected",uvrange=uvrange,imagename=file_str,imsize=[imsize],cell=cell,stokes=stokes,gridder='standard',\
						pblimit=-1,deconvolver="multiscale",scales=scales,nterms=1,weighting="briggs",robust=1.0,uvtaper=uvtaper,niter=100000000000,gain=gain,casalogger=casalog,\
						threshold=threshold_list,interactive=False,usemask='auto-multithresh',negativethreshold=3.0,savemodel='modelcolumn')
					else:
						imaging_cmd='poltclean(vis=\''+msname+'\',selectdata=True,datacolumn="corrected",imagename=\''+str(file_str)+'\',imsize=['+str(imsize)+'],cell=\''+str(cell)+\
						'\',stokes=\''+str(stokes)+'\',gridder=\'standard\',pblimit=-1,deconvolver="multiscale",scales='+str(scales)+',nterms=1,weighting="briggs",robust=1.0,uvtaper=\''+\
						str(uvtaper)+'\',casalogger='+str(casalog)+',niter=100000000000,threshold='+str(threshold_list)+',interactive=False,usemask="user",mask=\''+str(mask_str)+\
						',savemodel=\'modelcolumn\')'
						print (imaging_cmd+'\n')
						poltclean(vis=msname,selectdata=True,datacolumn="corrected",uvrange=uvrange,imagename=file_str,imsize=[imsize],cell=cell,stokes=stokes,gridder='standard',\
						pblimit=-1,deconvolver="multiscale",scales=scales,nterms=1,weighting="briggs",robust=1.0,uvtaper=uvtaper,casalogger=casalog,\
						niter=100000000000,threshold=threshold_list,interactive=False,usemask="user",mask=mask_str,savemodel='modelcolumn')
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
					wsclean_args=['-scale '+str(cell)+'asec','-size '+str(imsize)+' '+str(imsize),'-no-dirty','-j '+str(cpus),'-abs-mem '+str(absmem),'-weight briggs 0.8',\
					'-taper-tukey 10','-name '+str(file_str),'-nwlayers 1','-pol '+str(pol),'-maxuv-l '+str(imaging_maxuv),\
					'-minuv-l '+str(imaging_minuv),'-niter 1000000','-mgain 0.9','-quiet','-auto-threshold '+str(sigma),'-auto-mask '+str(sigma+0.5),\
					'-gain 0.1','-multiscale','-multiscale-scales '+','.join(scales)]
					if clean_continue==True:
						wsclean_args.append('-continue')
					imaging_cmd='wsclean '+' '.join(wsclean_args)+' '+msname
					print (imaging_cmd+'\n')			
					os.system('wsclean '+' '.join(wsclean_args)+' '+msname)
					wsclean_images=glob.glob(file_str+'*image*.fits')
					wsclean_models=glob.glob(file_str+'*model*.fits')
					wsclean_residuals=glob.glob(file_str+'*residual*.fits')
					print ('Converting WSClean images : '+','.join([os.path.basename(i) for i in wsclean_images])+' to CASA image : '+file_str+'.image\n')
					casa_imagename=PSC.wsclean_to_casaimage(wsclean_images=wsclean_images,casaimage_prefix=file_str,imagetype='image',keep_wsclean_images=True)
					print ('Converting WSClean models : '+','.join([os.path.basename(i) for i in wsclean_models])+' to CASA model : '+file_str+'.model\n')
					casa_modelname=PSC.wsclean_to_casaimage(wsclean_images=wsclean_models,casaimage_prefix=file_str,imagetype='model',keep_wsclean_images=True)
					print ('Converting WSClean residuals : '+','.join([os.path.basename(i) for i in wsclean_residuals])+' to CASA residual : '+file_str+'.residual\n')
					casa_residualname=PSC.wsclean_to_casaimage(wsclean_images=wsclean_residuals,casaimage_prefix=file_str,imagetype='residual',keep_wsclean_images=True)
				if count==0:
					qucor_image,qucor_model,qchange,uchange,vchange=PSC.correct_solar_quv_leakage(casa_imagename,casa_modelname,sigma,overwrite=False)
					qucor_image,qucor_model=PSC.pol_model_threshold(qucor_image,qucor_model,sigma,1)
					casa_modelname=qucor_model
					if qchange>=0.05 or uchange>=0.05 or vchange>=0.05 or do_diffcal==True:
						if do_diffcal==True:
							cal_cause='Differential calibration'
						elif qchange>=0.05 or uchange>=0.05 or vchange>=0.05:
							cal_cause='Stokes leakage decreases.'
						print ('Performing differential calibration because : '+cal_cause+'\n')
						delmod(vis=msname,scr=True,otf=True)
						ft(vis=msname,model=casa_modelname,usescratch=True)
						AM=AccessMS(msname)
						timeres=AM.calc_timeres()
						tstamps=AM.get_num_timestamps()
						tintg=int((timeres*tstamps)/10)
						cal.calibrate(msname=msname,caltable=file_str+'.cal',verbose=False,j=3,t=tintg,a='0.001,0.0001')
						cal.applycal(msname=msname,gaintable=file_str+'.cal')
					do_uvsub_ankflag(msname,model=casa_modelname,verbose=False,nthread=3,extendpols=True,chantime_minfrac=0.8,casaflag='rflag')
					count+=1
					if qchange>=0.05 or uchange>=0.05 or vchange>=0.05 or do_diffcal==True:
						clean_continue=False
					else:
						clean_continue=True
					continue
				median_res_frac,threshold=calc_residual_flux(casa_imagename,casa_residualname,5,rms_box,stokes_list=stokes_list)
				if median_res_frac>=residual_frac and sigma>=5 and count>=1:
					print ('Continuing CLEANing, since residual fraction is more than '+str(residual_frac*100)+'%\n')
					clean_continue=True
					sigma-=1.0
					count+=1
					continue
				else:
					break
			except Exception as e:
				print ('Error occured in final imaging : '+str(e)+'\n')
				return 2
		else:
			return 1
	# Exporting images
	##################		
	output=export_images(file_str,OBSID,cell,imsize,imaging_cmd=imaging_cmd,savedir=savedir,savemodel=savemodel,saveresidual=saveresidual,cutoutbox=cutoutbox,astrometry=False)
	os.system('rm -rf *.image *.model *pb *psf* *.sumwt *.flux *.fits')
	os.system('rm -rf '+file_str+'*')	
	os.system('cd ../')	
	if savedir!=workdir:
		os.system('rm -rf '+workdir)
	os.chdir(cwd)		
	return 0
	
# Function to run the script stand alone from command line
if __name__=='__main__':
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
	parser.add_option('--residual_frac',dest='resfrac',default=0.15,help='Residual flux fraction',metavar="Float")
	parser.add_option('--casa_caltables',dest='casacals',default='',help='CASA caltables',metavar="Comma separated string")
	parser.add_option('--calibrate_caltables',dest='calibratecals',default='',help='CALIBRATE caltables',metavar="Comma separated string")
	parser.add_option('--wsclean',dest="use_wsclean",default=True,help="Use WSClean for imaging or not",metavar="Boolean")
	parser.add_option('--do_diffcal',dest="do_diffcal",default=False,help="Use WSClean for imaging or not",metavar="Boolean")
	(options, args) = parser.parse_args()

	if str(options.msname)[-1]=='/':
		msname=str(options.msname)[:-1]
	else:
		msname=str(options.msname)

	if os.path.isdir(str(options.msname))==False or options.msname==None:
		print ('Measurement set is not present.\n')
		os.system('touch '+cwd+'/.Finished_final_imaging_'+os.path.basename(str(msname))+'_noms')
		os._exit(1)
	elif os.path.isfile(str(options.metafits))==False or options.metafits==None:
		print ('Metafits file is not present.\n')
		os.system('touch '+cwd+'/.Finished_final_imaging_'+os.path.basename(str(msname))+'_nometa')
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
				print ('applycal(vis=\''+options.msname+'\',gaintable='+str(casacals)+',applymode=\'calflag\',flagbackup=False)\n')
				a=applycal(vis=options.msname,gaintable=casacals,applymode='calflag',flagbackup=False)
				
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
				for i in calibratecals:
					print ('cal.applycal(msname=\''+options.msname+'\',gaintable=\''+str(i)+'\',applymode=\'calflag\',flagbackup=False,verbose=False)\n')
					cal.applycal(msname=options.msname,gaintable=i,applymode='calflag',flagbackup=False,verbose=False,datacolumn='corrected')
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
		
		error_files=glob.glob(basedir+'/.Finished_final_imaging_*error*')
		for i in error_files:
			os.system('rm -rf '+i)
		touch_file=basedir+'/.Finished_final_imaging_'+os.path.basename(msname)+'_success'
		touch_file_moreflag=basedir+'/.Finished_final_imaging_'+os.path.basename(msname)+'_moreflag'
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

		output=make_image(str(options.msname),str(options.metafits),str(workdir),sigma=float(options.sigma),stokes=str(options.stokes),savedir=str(savedir),threshold=threshold,\
				want_automask=eval(str(options.want_automask)),maskfile=maskfile,quality_factor=int(options.quality_factor),savemodel=eval(str(options.savemodel)),\
				saveresidual=eval(str(options.saveresidual)),cutoutbox=str(options.cutoutbox),use_ankflagger=eval(str(options.use_ankflag)),residual_frac=float(options.resfrac),\
				use_wsclean=eval(str(options.use_wsclean)),do_diffcal=eval(str(options.do_diffcal)))
		if output==0:
			print ('\nImaging finished.\n#############################\n')
			os.system('touch '+touch_file)
		elif output==1:
			print ('More than 95% data are flagged. Image is not made.\n############################\n')
			touch_file=basedir+'/.Finished_final_imaging_'+os.path.basename(msname)+'_moreflag'
			os.system('touch '+touch_file)
		elif output==2:
			print ('Error occured during final imaging.\n############################\n')
			touch_file=basedir+'/.Finished_final_imaging_'+os.path.basename(msname)+'_error'
			os.system('touch '+touch_file)
		print ('Total run time : '+str(time.time()-start_time)+' s\n#############################\n')





