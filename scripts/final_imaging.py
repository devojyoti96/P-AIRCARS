from casatasks import *
from paircars.basic_func import *
from paircars.access_ms import *
from paircars_casatasks.poltclean import *
from paircars.flagger import *
from astropy.io import fits

def modify_header(imagename,inputfile,tclean_dic={},astrometry=False):
	header=fits.getheader(imagename)
	if inputfile!='' or os.path.isfile(inputfile)==True:
		with open(inputfile,'r') as fil:
			lines=fil.readlines()
		for line in lines:
			if 'timerange' in line:
				skip_header=int(lines.index(line))
		input_params,input_vals=np.genfromtxt(inputfile,delimiter='=',dtype=str,autostrip=True,skip_header=skip_header,skip_footer=True,usecols=(0,1),unpack=True) 
		for key,value in zip(inp_params,inp_vals):
			header[key]=str(value)
	if len(tclean_dic)!=0:
		tclean_keys=tclean_dic.keys()
		tclean_values=tclean_dic.values()
		for key,value in zip(tclean_keys,tclean_values):
			header[key]=str(value)
	header['PIPELINE']='P-AIRCARS'
	header['Devoloper']='Devojyoti Kansabanik, Surajit Mondal'
	header=['astrometry_corrected']=str(astrometry)
	fits.writeto(imagename,data=fits.getdata(imagename),header=header,overwrite=True)
	return 

def export_images(imagename,OBSID,savedir='',savemodel=False,saveresidual=False,cutoutbox=[5,5],poltclean_dict={},inputfile='',astrometry=False): #TODO :Other formats in heliocentric coordinate
	'''
	Function to save final images
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

	if len(cutoutbox)!=0:
		x_pix=int((cutoutbox[0]*3600)/self.cellsize)
		y_pix=int((cutoutbox[1]*3600)/self.cellsize)
		x_cen=int(self.imsize[0]/2)
		y_cen=x_cen
		box=str(int(x_cen-x_pix/2))+','+str(int(x_cen+x_pix/2))+','+str(int(y_cen-y_pix/2))+','+str(int(y_cen+y_pix/2))
		imsubimage(imagename=imagename+'.image',outfile='cutout.image',box=box)
		exportfits(imagename='cutout.image',fitsimage=imagedir+'/'+imagename+'_image.fits',history=False)
		modify_header(imagedir+'/'+imagename+'_image.fits',inputfile,tclean_dic=poltclean_dict,astrometry=astrometry)
		if savemodel==True:
			imsubimage(imagename='cutout.model',outfile='cutout.model',box=box)
			exportfits(imagename='cutout.model',fitsimage=modeldir+'/'+imagename+'_model.fits',history=False)
			modify_header(modeldir+'/'+imagename+'_model.fits',inputfile,tclean_dic=poltclean_dict,astrometry=astrometry)
		if saveresidual==True:
			imsubimage(imagename=imagename+'.residual',outfile='cutout.residual',box=box)
			exportfits(imagename='cutout.residual',fitsimage=resdir+'/'+imagename+'_res.fits',history=False)
			modify_header(resdir+'/'+imagename+'_res.fits',inputfile,tclean_dic=poltclean_dict,astrometry=astrometry)
		os.system('rm -rf cutout*')
	else:
		exportfits(imagename=imagename+'.image',fitsimage=imagedir+'/'+imagename+'_image.fits',history=False)
		modify_header(imagedir+'/'+imagename+'_image.fits',inputfile,tclean_dic=poltclean_dict,astrometry=astrometry)
		if savemodel==True:
			exportfits(imagename=imagename+'.model',fitsimage=modeldir+'/'+imagename+'_model.fits',history=False)
			modify_header(modeldir+'/'+imagename+'_model.fits',inputfile,tclean_dic=poltclean_dict,astrometry=astrometry)
		if saveresidual==True:
			exportfits(imagename=imagename+'.residual',fitsimage=resdir+'/'+imagename+'_res.fits',history=False)
			modify_header(resdir+'/'+imagename+'_res.fits',inputfile,tclean_dic=poltclean_dict,astrometry=astrometry)
	os.chdir(imagedir)
	os.system('rm -rf *.pb *.mask *.model *.image *.flux *.sumwt *.residual *.psf')
	os.chdir(cwd)
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


def make_image(msname,workdir,stokes='I',savedir='',threshold=0.1,imsize=[],cell='',scales=[],want_automask=False,maskfile='',uvtaper='',quality_factor=1,\
				savemodel=False,saveresidual=False,cutoutbox='',inputfile=''): #TODO : Wide FOV beam correction
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
	stokes_list=get_stokes(stokes)

	rmsthresh=[str(thresh*sigma)+'Jy']*len(stokes_list)
	mask_rad=int((32*60)/ISC.cellsize) # Creating a mask with 32 arcmin radius centered on the image
	mask_str='circle[['+str(ISC.imsize/2)+'pix,'+str(ISC.imsize/2)+'pix],'+str(mask_rad)+'pix]'
	poltclean(vis=msname,selectdata=True,datacolumn="corrected",imagename='temp_image',imsize=imsize,cell=cellsize,stokes=stokes,gridder='standard',\
	pblimit=-1,deconvolver="multiscale",scales=scales,nterms=1,weighting="natural",uvtaper=uvtaper,\
	niter=100000000000,gain=0.5,threshold=rmsthresh,cycleniter=-1,cyclefactor=1.0,interactive=False,usemask="user",mask=mask_str,savemodel='modelcolumn')

	rms_box='50,50,'+str(imsize-50)+','+str(int(imsize/4)) # CASA box to calculate the rms
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
	export_images(imagename,OBSID,savedir=savedir,savemodel=False,saveresidual=False,cutoutbox=cutoutbox,poltclean_dict=poltclean_dict,inputfile=inputfile,astrometry=False)
	os.system('rm -rf '+file_str+'*')	
	os.system('cd ../')	
	os.system('rm -rf '+workdir)
	os.chdir(cwd)		
	return 0
	

# Function to run the script stand alone from command line
if __name__='__main__':
	usage= ' Perform final imaging\n'
	parser = OptionParser(usage=usage)
	parser.add_option('--msname',dest="chantime_msname",default=None,help="Name of measurement set of a single time anf frequency slice",metavar="Measurement Set")
	parser.add_option('--basedir',dest='basedir',default=None,help='Name of the base directory',metavar='Directory path')
	parser.add_option('--workdir',dest='workdir',default=None,help='Name of the working directory',metavar='Directory path')
	parser.add_option('--savedir',dest='savedir',default=None,help='Directory name to save final images',metavar="Boolean")
	parser.add_option('--savemodel',dest='savemodel',default=False,help='Want to save final models',metavar="Boolean")
	parser.add_option('--saveres',dest="saveresidual",default=False,help="Want to save residual images",metavar="Boolean")
	parser.add_option('--stokes',dest='stokes',default='pseudoI',help='Stokes planes to image',metavar="String")
	parser.add_option('--cutoutbox',dest='cutoutbox',default='5,5',help='Cutout box \'X_width,Y_width\' in degree',metavar="Comma separated string")
	parser.add_option('--threshold',dest='thresh',default=0.1,help='RMS threshold for cleaning',metavar="Float")
	parser.add_option('--imsize',dest='imsize',default=1024,help='Number of pixels in the image',metavar="Integer")
	parser.add_option('--cell',dest='cell',default=1.0,help='Pixel size in arcsec',metavar="Float")
	parser.add_option('--scales',dest='scales',default='0,3,6,9',help='Multiscale scales in number of pixels',metavar="Comma separated string")
	parser.add_option('--want_automask',dest='want_automask',default=False,help='Want auto masking or not',metavar="Boolean")
	parser.add_option('--maskfile',dest='maskfile',default='',help='Mask for imaging when auto masking is off',metavar="Maskfile or CASA mask dtring")
	parser.add_option('--uvtaper',dest='uvtaper',default='',help='UV taper string',metavar="String")
	parser.add_option('--quality_factor',dest='quality_factor',default=1,help='Quality factor of imaging',metavar="Integer")
	parser.add_option('--inputfile',dest='inputfile',default='',help='Path of the P-AIRCARS input file',metavar="File path")
	(options, args) = parser.parse_args()

	if os.path.isdir(str(options.msname))==False or options.msname==None:
		



	multiscales=str(options.scales).split(',')
	make_image(str(options.msname),str(options.workdir),stokes=str(options.stokes),savedir=str(options.savedir),threshold=float(options.threshold),imsize=[int(options.imsize)],\
			cell=float(options.cell),scales=multiscales,want_automask=eval(str(options.want_automask)),maskfile=str(options.maskfile),uvtaper=str(options.uvtaper),\
			quality_factor=int(options.quality_factpr),savemodel=eval(str(options.savemodel)),saveresidual=eval(str(options.saveresidual)),\
			cutoutbox=str(options.cutoutbox),inputfile=str(options.inputfile))
	touch_file=basedir+'/.Finished_final_imaging_'+os.path.basename(msname).split('.ms')[0]
	os.system('touch '+touch_file)
	





