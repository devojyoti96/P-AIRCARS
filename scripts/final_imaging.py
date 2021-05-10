from casatasks import *
from paircars.basic_func import *
from paircars.access_ms import *
from paircars_casatasks.poltclean import *
from optparse import OptionParser
from paircars.flagger import *
from astropy.io import fits
from CALIBRATE.access_calibrate import *
import time


def calc_residual_flux(imagename,nsigma,rms_box,stokes_list=['I']):
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
	imagename=imagename
	residual=imagename.split('.image')[0]+'.residual'
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


def modify_header(imagename,inputfile,tclean_dic={},astrometry=False):
	header=fits.getheader(imagename)
	if inputfile!='' or os.path.isfile(inputfile)==True:
		with open(inputfile,'r') as fil:
			lines=fil.readlines()
		for line in lines:
			if 'final_image_dir' in line:
				skip_header=int(lines.index(line))
		input_params,input_vals=np.genfromtxt(inputfile,delimiter='=',dtype=str,autostrip=True,skip_header=skip_header,skip_footer=True,usecols=(0,1),unpack=True) 
		for key,value in zip(input_params,input_vals):
			if key!='email':
				header[key]=str(value)
	if len(tclean_dic)!=0:
		tclean_keys=tclean_dic.keys()
		tclean_values=tclean_dic.values()
		for key,value in zip(tclean_keys,tclean_values):
			header[key]=str(value)
	header['PIPELINE']='P-AIRCARS'
	header['Devoloper']='Devojyoti Kansabanik, Surajit Mondal'
	header['astrometry_corrected']=str(astrometry)
	fits.writeto(imagename,data=fits.getdata(imagename),header=header,overwrite=True)
	return 

def export_images(imagename,OBSID,cell,imsize,savedir='',savemodel=False,saveresidual=False,cutoutbox=[],poltclean_dict={},inputfile='',astrometry=False): 
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
	poltclean_dic = poltclean parameters dictionary
	inputfile = P-AIRCARS input file
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
		x_cen=int(imsize[0]/2)
		y_cen=x_cen
		box=str(int(x_cen-x_pix/2))+','+str(int(y_cen-y_pix/2))+','+str(int(x_cen+x_pix/2))+','+str(int(y_cen+y_pix/2))
		os.system('rm -rf '+imagename+'.cutout*')
		os.system('rm -rf '+imagedir+'/'+os.path.basename(imagename)+'_*.fits')
		os.system('rm -rf '+modeldir+'/'+os.path.basename(imagename)+'_*.fits')
		os.system('rm -rf '+resdir+'/'+os.path.basename(imagename)+'_*.fits')
		imsubimage(imagename=imagename+'.image',outfile=imagename+'.cutout.image',box=box)
		exportfits(imagename=imagename+'.cutout.image',fitsimage=imagedir+'/'+os.path.basename(imagename)+'_image.fits',history=False)
		output.append(imagedir+'/'+os.path.basename(imagename)+'_image.fits')
		modify_header(imagedir+'/'+os.path.basename(imagename)+'_image.fits',inputfile,tclean_dic=poltclean_dict,astrometry=astrometry)
		if savemodel==True:
			imsubimage(imagename=imagename+'.model',outfile=imagename+'.cutout.model',box=box)
			exportfits(imagename=imagename+'.cutout.model',fitsimage=modeldir+'/'+os.path.basename(imagename)+'_model.fits',history=False)
			output.append(modeldir+'/'+os.path.basename(imagename)+'_model.fits')
			modify_header(modeldir+'/'+os.path.basename(imagename)+'_model.fits',inputfile,tclean_dic=poltclean_dict,astrometry=astrometry)
		if saveresidual==True:
			imsubimage(imagename=imagename+'.residual',outfile=imagename+'.cutout.residual',box=box)
			exportfits(imagename=imagename+'.cutout.residual',fitsimage=resdir+'/'+os.path.basename(imagename)+'_res.fits',history=False)
			output.append(resdir+'/'+os.path.basename(imagename)+'_res.fits')
			modify_header(resdir+'/'+os.path.basename(imagename)+'_res.fits',inputfile,tclean_dic=poltclean_dict,astrometry=astrometry)
		os.system('rm -rf '+imagename+'.cutout*')
	else:
		os.system('rm -rf '+imagedir+'/'+os.path.basename(imagename)+'_*.fits')
		os.system('rm -rf '+modeldir+'/'+os.path.basename(imagename)+'_*.fits')
		os.system('rm -rf '+resdir+'/'+os.path.basename(imagename)+'_*.fits')
		exportfits(imagename=imagename+'.image',fitsimage=imagedir+'/'+os.path.basename(imagename)+'_image.fits',history=False)
		output.append(imagedir+'/'+os.path.basename(imagename)+'_image.fits')
		modify_header(imagedir+'/'+os.path.basename(imagename)+'_image.fits',inputfile,tclean_dic=poltclean_dict,astrometry=astrometry)
		if savemodel==True:
			exportfits(imagename=imagename+'.model',fitsimage=modeldir+'/'+os.path.basename(imagename)+'_model.fits',history=False)
			output.append(modeldir+'/'+os.path.basename(imagename)+'_model.fits')
			modify_header(modeldir+'/'+os.path.basename(imagename)+'_model.fits',inputfile,tclean_dic=poltclean_dict,astrometry=astrometry)
		if saveresidual==True:
			exportfits(imagename=imagename+'.residual',fitsimage=resdir+'/'+os.path.basename(imagename)+'_res.fits',history=False)
			output.append(resdir+'/'+os.path.basename(imagename)+'_res.fits')
			modify_header(resdir+'/'+os.path.basename(imagename)+'_res.fits',inputfile,tclean_dic=poltclean_dict,astrometry=astrometry)
	os.chdir(imagedir)
	os.system('rm -rf *.pb *.mask *.model *.image *.flux *.sumwt *.residual *.psf')
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


def make_image(msname,metafits,workdir,sigma=10,stokes='I',savedir='',threshold=[0.1],imsize=[],cell=1.0,scales=[],want_automask=False,maskfile='',uvtaper='',quality_factor=1,\
				savemodel=False,saveresidual=False,cutoutbox='',inputfile='',use_ankflagger=False,residual_frac=0.1): #TODO : Wide FOV beam correction
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
	imsize = [] ,Number of pixel of the image
	cell = Cellsize in arcsecond
	scales = [] , multiscale scales
	want_automask = Want to use auto masking
	maskfile = Name of any previous mask file or CASA mask string
	uvtaper = UV-taper for imaging
	quality_factor = Quality factor for imaging
	savemodel = Save model images or not
	saveresidual = Save residual images or not
	cutoutbox = Cutout box of the final image [x_width,y_width] in degree (default : [], no cutout, save full image)
	inputfile = P-AIRCARS input file
	use_ankflagger = Whether use aNKflagger for flagging or not
	residual_frac = Residual flux fraction to stop CLEANing
	Result:
	Name of final image,model,residual
	'''
	print ('Making image.....\n')
	casalog=False
	cwd=os.getcwd()
	os.chdir(workdir)
	if cutoutbox!='':
		cutoutbox=cutoutbox.split(',')
	else:
		cutoutbox=[]
	AM=AccessMS(msname)
	freqs=AM.calc_meanfreq()/10**6	
	multiscales=str(scales).split(',')		
	scales=[int(i) for i in multiscales]
	OBSID=str(fits.getheader(metafits)['GPSTIME'])	
	file_str=workdir+'/'+os.path.basename(splited_ms_rename(msname,ref_time_chan=False,change_msname=False)).split('.ms')[0]
	if quality_factor==0:
		gain=0.3
		cyclefactor=0.5
	elif quality_factor==1:
		gain=0.15
		cyclefactor=0.7
	else:
		gain=0.1				
		cyclefactor=1.0
	sigma+=1
	stokes_list=get_stokes(stokes)
	imagename=file_str+'.image'
	residual=file_str+'.residual'
	mask=file_str+'.mask'

	if len(threshold)!=len(stokes_list):	
		threshold=[threshold[0]]*len(stokes_list)

	threshold_list=[str(rms*sigma)+'Jy' for rms in threshold]

	os.system('rm -rf '+file_str+'.image '+file_str+'.model '+file_str+'.mask '+file_str+'.psf '+file_str+'.pb '+file_str+'.sumwt '+file_str+'.flux '+file_str+'.residual')

	if maskfile=='':
		mask_rad=int((32*60)/float(cell)) # Creating a mask with 32 arcmin radius centered on the image
		mask_str='circle[['+str(imsize[0]/2)+'pix,'+str(imsize[0]/2)+'pix],'+str(mask_rad)+'pix]'
	else:
		mask_str=inputs.maskstr
	while True:
		if os.path.isdir(mask)==True:
			poltclean(vis=msname,selectdata=True,datacolumn="corrected",imagename=file_str,imsize=imsize,cell=cell,stokes=stokes,gridder='standard',\
			pblimit=-1,deconvolver="multiscale",scales=scales,nterms=1,weighting="natural",uvtaper=uvtaper,casalogger=casalog,\
			niter=100000000000,gain=gain,threshold=threshold_list,cycleniter=-1,cyclefactor=cyclefactor,interactive=False,usemask="user",mask=mask,savemodel='modelcolumn')
			poltclean_dict={'vis':msname,'datacolumn':"corrected",'imsize':imsize,'cell':cell,\
			'stokes':stokes,'pblimit':-1,'deconvolver':"multiscale",'scales':scales,'nterms':1,'weighting':"natural",'uvtaper':uvtaper,\
			'niter':100000000000,'gain':gain,'threshold':threshold,'cycleniter':-1,'cyclefactor':cyclefactor,'usemask':"user",'mask':maskfile}
		elif want_automask==True and maskfile=='':
			poltclean(vis=msname,selectdata=True,datacolumn="corrected",imagename=file_str,imsize=imsize,cell=cell,stokes=stokes,gridder='standard',\
			pblimit=-1,deconvolver="multiscale",scales=scales,nterms=1,weighting="natural",uvtaper=uvtaper,niter=100000000000,gain=gain,casalogger=casalog,\
			threshold=threshold_list,cycleniter=-1,cyclefactor=cyclefactor,interactive=False,usemask='auto-multithresh',negativethreshold=3.0,savemodel='modelcolumn')
			poltclean_dict={'vis':msname,'datacolumn':"corrected",'imsize':imsize,'cell':cell,'stokes':stokes,\
			'pblimit':-1,'deconvolver':"multiscale",'scales':scales,'nterms':1,'smallscalebias':0.0,'weighting':"natural",'uvtaper':uvtaper,\
			'niter':100000000000,'gain':gain,'threshold':threshold,'cycleniter':-1,'cyclefactor':cyclefactor,'usemask':'auto-multithresh','negativethreshold':0.0}
		else:
			poltclean(vis=msname,selectdata=True,datacolumn="corrected",imagename=file_str,imsize=imsize,cell=cell,stokes=stokes,gridder='standard',\
			pblimit=-1,deconvolver="multiscale",scales=scales,nterms=1,weighting="natural",uvtaper=uvtaper,casalogger=casalog,\
			niter=100000000000,gain=0.5,threshold=threshold_list,cycleniter=-1,cyclefactor=cyclefactor,interactive=False,usemask="user",mask=mask_str,savemodel='modelcolumn')
			poltclean_dict={'vis':msname,'datacolumn':"corrected",'imsize':imsize,'cell':cell,'stokes':stokes,'pblimit':-1,'deconvolver':"multiscale",\
			'scales':scales,'nterms':1,'weighting':"natural",'uvtaper':uvtaper,'niter':100000000000,'gain':gain,'threshold':threshold,'cycleniter':-1,\
			'cyclefactor':cyclefactor,'usemask':"user",'mask':mask_str}
		rms_box='50,50,'+str(imsize[0]-50)+','+str(int(imsize[0]/4))
		median_res_frac,threshold=calc_residual_flux(file_str+'.image',sigma,rms_box,stokes_list=stokes_list)
		if median_res_frac>=residual_frac:
			if use_ankflagger==True:
				do_uvsub_ankflag(msname,model='',nthread=1,verbose=False,flagbackup=False)
			else:
				do_uvsub_flagger(msname,model='',mode='uvsub',rmsthresh=[10,8,6,4],flagbackup=False)
			print ('Continuing CLEANing, since residual fraction is more than '+str(residual_frac*100)+'%\n')
			sigma-=1.0
			continue
		else:
			break
	# Exporting images
	##################		
	print (cell)
	output=export_images(file_str,OBSID,cell,imsize,savedir=savedir,savemodel=savemodel,saveresidual=saveresidual,cutoutbox=cutoutbox,\
					poltclean_dict=poltclean_dict,inputfile=inputfile,astrometry=False)
	os.system('rm -rf '+file_str+'*')	
	os.system('cd ../')	
	os.system('rm -rf '+workdir)
	os.chdir(cwd)		
	return output
	
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
	parser.add_option('--imsize',dest='imsize',default=1024,help='Number of pixels in the image',metavar="Integer")
	parser.add_option('--cell',dest='cell',default=1.0,help='Pixel size in arcsec',metavar="Float")
	parser.add_option('--scales',dest='scales',default='0,3,6,9',help='Multiscale scales in number of pixels',metavar="Comma separated string")
	parser.add_option('--want_automask',dest='want_automask',default=False,help='Want auto masking or not',metavar="Boolean")
	parser.add_option('--maskfile',dest='maskfile',default='',help='Mask for imaging when auto masking is off',metavar="Maskfile or CASA mask string")
	parser.add_option('--uvtaper',dest='uvtaper',default='',help='UV taper string',metavar="String")
	parser.add_option('--quality_factor',dest='quality_factor',default=1,help='Quality factor of imaging',metavar="Integer")
	parser.add_option('--inputfile',dest='inputfile',default='',help='Path of the P-AIRCARS input file',metavar="File path")
	parser.add_option('--use_ankflag',dest='use_ankflag',default=False,help='Use aNKflag for flagging or not',metavar="Boolean")
	parser.add_option('--residual_frac',dest='resfrac',default=0.15,help='Residual flux fraction',metavar="Float")
	(options, args) = parser.parse_args()

	if str(options.msname)[-1]=='/':
		msname=str(options.msname)[:-1]
	else:
		msname=str(options.msname)

	if os.path.isdir(str(options.msname))==False or options.msname==None:
		print ('Measurement set is not present.\n')
		os.system('touch '+cwd+'/.Finished_final_imaging_'+os.path.basename(str(msname))+'_noms')
		os._exit(0)
	elif os.path.isfile(str(options.metafits))==False or options.metafits==None:
		print ('Metafits file is not present.\n')
		os.system('touch '+cwd+'/.Finished_final_imaging_'+os.path.basename(str(msname))+'_nometa')
		os._exit(0)
	else:
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
		
		touch_file=basedir+'/.Finished_final_imaging_'+os.path.basename(msname)+'_success'
		if os.path.exists(touch_file):
			print ('Image is already made.\n#############################\n')
			os._exit(0)

		rmsthresh=str(options.threshold).split(',')
		threshold=[float(i) for i in rmsthresh]

		make_image(str(options.msname),str(options.metafits),str(workdir),sigma=float(options.sigma),stokes=str(options.stokes),\
				savedir=str(savedir),threshold=threshold,imsize=[int(options.imsize)],cell=float(options.cell),\
				scales=str(options.scales),want_automask=eval(str(options.want_automask)),maskfile=str(options.maskfile),uvtaper=str(options.uvtaper),\
				quality_factor=int(options.quality_factor),savemodel=eval(str(options.savemodel)),saveresidual=eval(str(options.saveresidual)),\
				cutoutbox=str(options.cutoutbox),inputfile=str(options.inputfile),use_ankflagger=eval(str(options.use_ankflag)),residual_frac=float(options.resfrac))
		print ('\nImaging finished.\n#############################\n')
		os.system('touch '+touch_file)
		print ('Total run time : '+str(time.time()-start_time)+' s\n#############################\n')
	





