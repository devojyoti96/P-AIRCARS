from casatasks import *
from paircars.basic_func import *
from paircars.access_ms import *
from paircars_casatasks.poltclean import *
from paircars.flagger import *


def export_images(imagename,savedir='',savemodel=False,saveresidual=False,savecutout=False,cutoutbox=[5,5],imageformat='FITS',poltclean_dict={},astrometry=False): #TODO :Other formats in heliocentric coordinate
	'''
	Function to save final images
	'''
	inplines=[]
	with open('selfcal_input.py','r') as fil:
		lines=fil.readlines()
	for line in lines:
		if 'safety_standard' in line:
			skip_header=int(lines.index(line))
	input_params,input_vals=np.genfromtxt('selfcal_inputs.py',delimiter='=',dtype=str,autostrip=True,skip_header=skip_header,skip_footer=True,usecols=(0,1),unpack=True) 
	for key,value in zip(inp_params,inp_vals):
		imhead(imagename=imagename+'.image',mode='add',hdkey=key,hdvalue=value)
	tclean_keys=poltclean_dict.keys()
	tclean_values=poltclean_dict.values()
	for key,value in zip(tclean_keys,tclean_values):
		imhead(imagename=imagename+'.image',mode='add',hdkey=key,hdvalue=value)
	if savedir=='':
		os.makedirs(basedir+'/All_final_images')
		imagedir=basedir+'/All_final_images'
		if savemodel==True:
			os.makedirs(basedir+'/All_final_models')
			modeldir=basedir+'/All_final_models'
		if saveresidual==True:
			os.makedirs(basedir+'/All_final_residuals')
			resdir=basedir+'/All_final_residuals'
	else:
		os.makedirs(savedir+'/All_final_images')
		imagedir=savedir+'/All_final_images'
		if savemodel==True:
			os.makedirs(savedir+'/All_final_models')
			modeldir=savedir+'/All_final_models'
		if saveresidual==True:
			os.makedirs(savedir+'/All_final_residuals')
			resdir=savedir+'/All_final_residuals'

	if imageformat=='CASA' and savecutout==False:
		os.system('mv '+imagename+'.image '+imagedir)
		if savemodel==True:
			os.system('mv '+imagename+'.model '+modeldir)
		if saveresidual==True:
			os.system('mv '+imagename+'.residual '+resdir)
	elif imageformat=='CASA' and savecutout==True:
		x_pix=int((cutoutbox[0]*3600)/self.cellsize)
		y_pix=int((cutoutbox[1]*3600)/self.cellsize)
		x_cen=int(self.imsize[0]/2)
		y_cen=x_cen
		box=str(int(x_cen-x_pix/2))+','+str(int(x_cen+x_pix/2))+','+str(int(y_cen-y_pix/2))+','+str(int(y_cen+y_pix/2))
		imsubimage(imagename=imagename+'.image',outfile=imagedir+'/'+imagename+'.image',box=box)
		if savemodel==True:
			imsubimage(imagename=imagename+'.model',outfile=modeldir+'/'+imagename+'.model',box=box)
		if saveresidual==True:
			imsubimage(imagename=imagename+'.residual',outfile=resdir+'/'+imagename+'.residual',box=box)
	elif savecutout==False:
		exportfits(imagename=imagename+'.image ',fitsimage=imagedir+'/'+imagename+'_image.fits',history=False)
		if savemodel==True:
			exportfits(imagename=imagename+'.model ',fitsimage=modeldir+'/'+imagename+'_model.fits',history=False)
		if saveresidual==True:
			exportfits(imagename=imagename+'.residual ',fitsimage=resdir+'/'+imagename+'_res.fits',history=False)
	else:
		x_pix=int((cutoutbox[0]*3600)/self.cellsize)
		y_pix=int((cutoutbox[1]*3600)/self.cellsize)
		x_cen=int(self.imsize[0]/2)
		y_cen=x_cen
		box=str(int(x_cen-x_pix/2))+','+str(int(x_cen+x_pix/2))+','+str(int(y_cen-y_pix/2))+','+str(int(y_cen+y_pix/2))
		imsubimage(imagename=imagename+'.image',outfile='cutout.image',box=box)
		exportfits(imagename='cutout.image',fitsimage=imagedir+'/'+imagename+'_image.fits',history=False)
		if savemodel==True:
			imsubimage(imagename='cutout.model',outfile='cutout.model',box=box)
			exportfits(imagename='cutout.model',fitsimage=modeldir+'/'+imagename+'_model.fits',history=False)
		if saveresidual==True:
			imsubimage(imagename=imagename+'.residual',outfile='cutout.residual',box=box)
			exportfits(imagename='cutout.residual',fitsimage=resdir+'/'+imagename+'_res.fits',history=False)
		os.system('rm -rf cutout*')

	#if other_format!='':
	os.system('rm -rf *.pb *.mask *.model *.image *.flux *.sumwt *.residual *.psf')
	return 

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


def make_image(msname,workdir,stokes='I',savedir='',threshold=[0.1],imsize=[],cell='',scales=[],want_automask=False,maskfile='',uvtaper='',quality_factor=1,\
				savemodel=False,saveresidual=False,cutoutbox='5,5',imageformat='FITS'): #TODO : Wide FOV beam correction
	'''
	Function to make final images
	Parameters:		
	savedir = Name of the directory to save the image
	savecutout = False, save cut out images or not
	imageformat ='CASA' or 'FITS'
	'''
	cwd=os.getcwd()
	os.chdir(workdir)
	cutoutbox=cutoutbox.split(',')
	AM=AccessMS(msname)
	freqs=AM.calc_meanfreq()/10**6				
	file_str=os.path.basename(splited_ms_rename(msname,ref_time_chan=False,change_msname=False)).split('.ms')[0]
	if quality_factor==0:
		gain=0.3
		sigma=10
		cycleniter=6000
	elif quality_factor==1:
		gain=0.15
		sigma=7
		cycleniter=3000
	else:
		gain=0.1				
		sigma=5
		cycleniter=-1

	rmsthresh=[str(thresh*sigma)+'Jy' for i in threshold]
	mask_rad=int((32*60)/ISC.cellsize) # Creating a mask with 32 arcmin radius centered on the image
	mask_str='circle[['+str(ISC.imsize/2)+'pix,'+str(ISC.imsize/2)+'pix],'+str(mask_rad)+'pix]'
	poltclean(vis=msname,selectdata=True,datacolumn="corrected",imagename='temp_image',imsize=imsize,cell=cellsize,stokes=stokes,gridder='standard',\
	pblimit=-1,deconvolver="multiscale",scales=scales,nterms=1,weighting="natural",uvtaper=uvtaper,\
	niter=100000000000,gain=0.5,threshold=rmsthresh,cycleniter=-1,cyclefactor=1.0,interactive=False,usemask="user",mask=mask_str,savemodel='modelcolumn')

	rms_box='50,50,'+str(imsize-50)+','+str(int(imsize/4)) # CASA box to calculate the rms
	stokes_list=get_stokes(stokes)
	rms_list=[]
	for s in stokes_list: 
		rms_list.append(imstat(imagename='temp_image.residual',box=rms_box,stokes=s)['rms'][0])

	threshold=[str(rms*sigma)+'Jy' for rms in rms_list]

	if use_ankflagger==True:
		do_uvsub_ankflag(msname,model='',nthread=1,verbose=False,flagbackup=False)
	else:
		do_uvsub_flagger(msname,model='',mode='uvsub',rmsthresh=[10,8,6,4],flagbackup=False)

	os.system('rm -rf temp_image*')

	imagename=file_str+'.image'
	residual=file_str+'.residual'
	mask=file_str+'.mask'
	if mask!='':
		poltclean(vis=msname,selectdata=True,datacolumn="corrected",imagename=file_str,imsize=imsize,cell=cellsize,stokes=stokes,gridder='standard',\
		pblimit=-1,deconvolver="multiscale",scales=scales,nterms=1,weighting="natural",uvtaper=uvtaper,\
		niter=100000000000,gain=gain,threshold=threshold,cycleniter=cycleniter,cyclefactor=1.0,interactive=False,usemask="user",mask=[maskfile],savemodel='modelcolumn')
		poltclean_dict={'vis':msname,'datacolumn':"corrected",'imsize':imsize,'cell':cellsize,\
		'stokes':stokes,'gridder':gridder,'pblimit':-1,'deconvolver':"multiscale",'scales':scales,\
		'nterms':1,'weighting':"natural",'uvtaper':uvtaper,\
		'niter':100000000000,'gain':gain,'threshold':threshold,'cycleniter':cycleniter,'cyclefactor':1.0,'usemask':"user",'mask':maskfile}
	elif want_auto_masking==True and maskfile=='':
		poltclean(vis=msname,selectdata=True,datacolumn="corrected",imagename=file_str,imsize=self.imsize,cell=self.cellsize,stokes=self.stokes,gridder='standard',\
		pblimit=-1,deconvolver="multiscale",scales=scales,nterms=1,weighting="natural",uvtaper=uvtaper,niter=100000000000,gain=gain,\
		threshold=threshold,cycleniter=cycleniter,cyclefactor=1.0,interactive=False,usemask='auto-multithresh',negativethreshold=3.0,savemodel='modelcolumn')
		poltclean_dict={'vis':msname,'datacolumn':"corrected",'imsize':imsize,'cell':cellsize,'stokes':stokes,'gridder':gridder,\
		'pblimit':-1,'deconvolver':"multiscale",'scales':scales,'nterms':1,'smallscalebias':0.0,'weighting':"natural",'uvtaper':self.uvtaper,\
		'niter':100000000000,'gain':gain,'threshold':threshold,'cycleniter':cycleniter,'cyclefactor':1.0,'usemask':'auto-multithresh','negativethreshold':0.0}
	else:
		mask_rad=int((32*60)/ISC.cellsize) # Creating a mask with 32 arcmin radius centered on the image
		mask_str='circle[['+str(ISC.imsize/2)+'pix,'+str(ISC.imsize/2)+'pix],'+str(mask_rad)+'pix]'
		poltclean(vis=msname,selectdata=True,datacolumn="corrected",imagename=file_str,imsize=imsize,cell=cellsize,stokes=stokes,gridder='standard',\
		pblimit=-1,deconvolver="multiscale",scales=scales,nterms=1,weighting="natural",uvtaper=uvtaper,\
		niter=100000000000,gain=0.5,threshold=threshold,cycleniter=cycleniter,cyclefactor=1.0,interactive=False,usemask="user",mask=mask_str,savemodel='modelcolumn')
		poltclean_dict={'vis':msname,'datacolumn':"corrected",'imsize':imsize,'cell':cellsize,'stokes':stokes,'gridder':gridder,'pblimit':-1,'deconvolver':"multiscale",\
		'scales':scales,'nterms':1,'weighting':"natural",'uvtaper':uvtaper,'niter':100000000000,'gain':gain,'threshold':threshold,'cycleniter':cycleniter,\
		'cyclefactor':1.0,'usemask':"user",'mask':mask_str}

	# Exporting images
	##################		
	export_images(file_str,savedir=savedir,savemodel=savemodel,saveresidual=saveresidual,cutoutbox=cutoutbox,imageformat=imageformat,poltclean_dict=poltclean_dict)
	os.system('rm -rf '+file_str+'*')	
	os.system('cd ../')	
	os.system('rm -rf '+workdir)
	os.chdir(cwd)		
	return 0
	

# Function to run the script stand alone from command line
if __name__='__main__':
	usage= ' Perform final imaging.....\n'
	parser = OptionParser(usage=usage)
	parser.add_option('--msname',dest="chantime_msname",default=None,help="Name of measurement set of a single time anf frequency slice",metavar="Measurement Set")
	parser.add_option('--workdir',dest='workdir',default=None,help='Name of the working directory',metavar='Directory path')
	parser.add_option('--savedir',dest='savedir',default=None,help='Directory name to save final images',metavar="Boolean")
	parser.add_option('--savemodel',dest='savemodel',default=False,help='Want to save final models',metavar="Boolean")
	parser.add_option('--saveres',dest="saveresidual",default=False,help="Want to save residual images",metavar="Boolean")
	parser.add_option('--savecutout',dest='savecutout',default=False,help='Want to save curout images',metavar="Boolean")
	parser.add_option('--cutoutbox',dest='cutoutbox',default='5,5',help='Cutout box \'X_width,Y_width\' in degree',metavar="Comma separated string")
	parser.add_option('--imageformat',dest='imageformat',default='FITS',help='Output image format',metavar="String")
	parser.add_option('--other_format',dest='other_format',default='png',help='Other image formats to save',metavar="Comma separated string")
	parser.add_option('--plotcontour',dest='plotcontour',default=False,help='Plot image contours',metavar="Boolean")
	parser.add_option('--contour_levels',dest='contour_levels',default='0,0.2,0.4,0.6,0.8',help='Contour levels',metavar="Comma separated string")
	parser.add_option('--image_delta_freq',dest='delta_f',default=40.0,help='Image frequency resolution in kHz',metavar="Float")
	parser.add_option('--image_delta_time',dest='delta_t',default=0.5,help='Image time resolution in second',metavar="Float")
	(options, args) = parser.parse_args()
	make_image(workdir,delta_t,delta_f,savedir='',savecutout=False,imageformat='FITS',cutoutformat='png',plotcontour=False,contour_levels=[])
	touch_file=basedir+'/.Finished_final_'+workdir
	os.system('touch '+touch_file)
	





