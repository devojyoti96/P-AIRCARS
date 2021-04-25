from casatasks import *
from selfcal_inputs import *
from basic_func import *
from access_ms import *
from poltclean import *
from flagger import *


class FinalImage:
	'''
	Class to make final images (Either Stokes I or Full Stokes images)
	'''
	def __init__(self,msname,workdir):
		print ('Initiating final imaging\n')
		AM=AccessMS(self.msname)
		IB=ImageBasic(self.msname)
		self.msname=msname
		self.max_baseline=AM.get_max_baseline()	
		self.cellsize=IB.calc_cellsize(3) # Assuming 3 pixels in one PSF
		self.imsize=IB.num_pixels(3)
		self.max_size=maximum_emission_scale
		self.multiscale_scales=IB.choose_scales(3,self.max_size)
		self.uvtaper=IB.calc_uvtaper()
		self.calib_uvrange=IB.calc_calib_uvrange()[0]
		self.rms_box='50,50,'+str(self.imsize[0]-50)+','+str(int(self.imsize[0]/4)) # CASA box to calculate the rms
		self.sigma_step=sigma_step
		self.residual_frac=residual_frac
		self.safety_standard=safety_standard
		self.quality_factor=quality_factor
		self.verbose=verbose
		self.interactive=interactive
		self.min_sigma=min_sigma
		self.do_polcal=do_polcal
		self.wprojection=wprojection
		self.wprojplanes=wprojplanes
		self.start_sigma=int(np.load(basedir+'/Ref_time_chan_sigma.npy',allow_pickle=True)[0])
		self.rmsthresh=list(np.load(basedir+'/Ref_time_chan_sigma.npy',allow_pickle=True)[1])
		if self.verbose:
			log_verbose.info('Initiating final imaging')

	def export_images(self,imagename,savedir='',savemodel=False,saveresidual=False,savecutout=False,cutoutbox=[5,5],imageformat='FITS',\
					other_format=[],plotcontour=False,contour_levels=[],poltclean_dict={}): #TODO :Other formats in heliocentric coordinate
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

	def make_image(self,workdir,delta_t,delta_f,savedir='',savemodel=False,saveresidual=False,savecutout=False,cutoutbox='5,5',imageformat='FITS',\
					other_format='',plotcontour=False,contour_levels=''): #TODO : Wide FOV beam correction
		'''
		Function to make final images
		Parameters:
		delta_t = Image temporal resolution
		delta_f = Image frequency resolution
		savedir = Name of the directory to save the image
		savecutout = False, save cut out images or not
		imageformat ='CASA' or 'FITS'
		cutoutformat = 'fits','png','pdf','jpg','eps'
		'''
		cwd=os.getcwd()
		os.chdir(workdir)
		cutoutbox=cutoutbox.split(',')
		other_format=otherformat.split(',')
		contour_levels=contour_levels.split(',')
		AM=AccessMS(self.msname)
		tot_nchan=AM.get_num_channels()
		tot_ntime=AM.get_num_timestamps()
		timestamps=AM.get_timestamps()
		freqs=AM.get_freqs()
		nchan=int((delta_f*10**3)/AM.calc_freqres)
		ntime=int(delta_t/AM.calc_timeres())
		if self.do_polcal==True:
			self.stokes='IQUV'
			self.start_sigma_list=[self.start_sigma]*4
			corr=''
		else:
			self.stokes='I'
			self.start_sigma_list=[self.start_sigma]
			corr='XX,YY'
		for i in range(0,tot_ntime,ntime):
			for j in range(0,tot_chan,nchan):
				timestamp=timestamps[int(i+(ntime/2))]
				freq=freqs[int(j+(nchan/2))]
				file_str='time_'+timestamp+'_freq_'+str(freq/10**6)+'_'+stokes
				print ('Spliting : '+file_str+'.ms')
				split(vis=self.msname,outputvis=file_str+'.ms',spw='0:'+str(j)+'~'+str(j+nchan),timerange=timestamps[i]+'~'+timestamps[i+ntime],correlation=corr)
				msname=file_str+'.ms'
				self.threshold=[str(self.rmsthresh[i]*self.start_sigma_list[i])+'Jy' for i in range(len(self.start_sigma_list))]
				if wprojection==True:
					gridder='wproject'
					wprojplanes=wprojplanes
				else:
					gridder='wproject'
					wprojplanes=wprojplanes
				if maskfile!='':
					poltclean(vis=msname,selectdata=True,datacolumn="corrected",imagename='temp_image',imsize=self.imsize,cell=self.cellsize,stokes=self.stokes,gridder='standard',\
					wprojplanes=1,pblimit=-1,deconvolver="multiscale",scales=self.multiscale_scales,nterms=1,smallscalebias=0.0,weighting="natural",uvtaper=self.uvtaper,\
					niter=100000000000,gain=0.1,threshold=self.threshold,cycleniter=-1,cyclefactor=1.0,interactive=False,usemask="user",mask=[maskfile],savemodel='modelcolumn')
				elif want_auto_masking==True and maskfile=='':
					poltclean(vis=msname,selectdata=True,datacolumn="corrected",imagename='temp_image',imsize=self.imsize,cell=self.cellsize,stokes=self.stokes,gridder='standard',\
					wprojplanes=1,pblimit=-1,deconvolver="multiscale",scales=self.multiscale_scales,nterms=1,smallscalebias=0.0,weighting="natural",uvtaper=self.uvtaper,\
					niter=100000000000,gain=0.1,threshold=self.threshold,cycleniter=-1,cyclefactor=1.0,interactive=False,usemask='auto-multithresh',\
					pbmask=0.0,sidelobethreshold=5.0,noisethreshold=sigma,lownoisethreshold=3.0,negativethreshold=0.0,smoothfactor=1.0,\
					minbeamfrac=0.5,growiterations=75,minpercentchange=10.0,savemodel='modelcolumn')
				else:
					mask_rad=int((20*60)/self.cellsize) # Creating a mask with 20 arcmin radius centered on the image
					maskstr='circle[['+str(self.imsize[0]/2)+','+str(self.imsize[0]/2)+'],'+str(mask_rad)+'pix]'
					poltclean(vis=msname,selectdata=True,datacolumn="corrected",imagename='temp_image',imsize=self.imsize,cell=self.cellsize,stokes=self.stokes,gridder='standard',\
					wprojplanes=1,pblimit=-1,deconvolver="multiscale",scales=self.multiscale_scales,nterms=1,smallscalebias=0.0,weighting="natural",uvtaper=self.uvtaper,\
					niter=100000000000,gain=0.1,threshold=self.threshold,cycleniter=-1,cyclefactor=1.0,interactive=False,usemask="user",mask=mask_str,savemodel='modelcolumn')
				if use_ankflagger==True:
					do_uvsub_ankflag(msname,model='',nthread=1,verbose=False)
				else:
					do_uvsub_flagger(msname,model='',mode='uvsub_flag',rmsthresh=[10,8,6,4])
				os.system('rm -rf temp_image*')
				while True:
					if self.do_polcal==True:
						self.threshold=[str(self.rmsthresh[i]*self.start_sigma_list[i])+'Jy' for i in range(len(self.rmsthresh))]
					else:
						self.threshold=[str(self.rmsthresh[i]*self.start_sigma_list[i])+'Jy' for i in range(len(self.rmsthresh))]
					imagename=file_str+'.image'
					residual=file_str+'.residual'
					mask=file_str+'.mask'
					residual_frac_list=[]
					if mask!='':
						poltclean(vis=msname,selectdata=True,datacolumn="corrected",imagename=file_str,imsize=self.imsize,cell=self.cellsize,\
						stokes=self.stokes,gridder=gridder,wprojplanes=wprojplanes,pblimit=-1,deconvolver="multiscale",scales=self.multiscale_scales,nterms=1,smallscalebias=0.0,\
						weighting="natural",uvtaper=self.uvtaper,niter=100000000000,gain=0.05,threshold=self.threshold,cycleniter=2000,cyclefactor=1.0,\
						interactive=False,usemask="user",mask=[maskfile])
						poltclean_dict={'vis':msname,'datacolumn':"corrected",'imsize':self.imsize,'cell':self.cellsize,\
						'stokes':self.stokes,'gridder':gridder,'wprojplanes':wprojplanes,'pblimit':-1,'deconvolver':"multiscale",'scales':self.multiscale_scales,\
						'nterms':1,'smallscalebias':0.0,'weighting':"natural",'uvtaper':self.uvtaper,\
						'niter':100000000000,'gain':0.05,'threshold':self.threshold,'cycleniter':2000,'cyclefactor':1.0,'usemask':"user",'mask':maskfile}
					elif want_auto_masking==True and maskfile=='':
						poltclean(vis=msname,selectdata=True,datacolumn="corrected",imagename=file_str,imsize=self.imsize,cell=self.cellsize,\
						stokes=self.stokes,gridder=gridder,wprojplanes=wprojplanes,pblimit=-1,deconvolver="multiscale",scales=self.multiscale_scales,nterms=1,smallscalebias=0.0,\
						weighting="natural",uvtaper=self.uvtaper,niter=100000000000,gain=0.05,threshold=self.threshold,cycleniter=2000,cyclefactor=1.0,\
						interactive=False,usemask='auto-multithresh',pbmask=0.0,sidelobethreshold=5.0,\
						noisethreshold=self.start_sigma,lownoisethreshold=3.0,negativethreshold=0.0,smoothfactor=1.0,\
						minbeamfrac=0.5,growiterations=75,minpercentchange=10.0)
						poltclean_dict={'vis':msname,'datacolumn':"corrected",'imsize':self.imsize,'cell':self.cellsize,\
						'stokes':self.stokes,'gridder':gridder,'wprojplanes':wprojplanes,'pblimit':-1,'deconvolver':"multiscale",'scales':self.multiscale_scales,\
						'nterms':1,'smallscalebias':0.0,'weighting':"natural",'uvtaper':self.uvtaper,\
						'niter':100000000000,'gain':0.05,'threshold':self.threshold,'cycleniter':2000,'cyclefactor':1.0,\
						'usemask':'auto-multithresh','pbmask':0.0,'sidelobethreshold':5.0,'noisethreshold':self.start_sigma,\
						'lownoisethreshold':3.0,'negativethreshold':0.0,'smoothfactor':1.0,\
						'minbeamfrac':0.5,'growiterations':75,'minpercentchange':10.0}
					else:
						poltclean(vis=msname,selectdata=True,datacolumn="corrected",imagename=file_str,imsize=self.imsize,cell=self.cellsize,\
						stokes=self.stokes,gridder=gridder,wprojplanes=wprojplanes,pblimit=-1,deconvolver="multiscale",scales=self.multiscale_scales,nterms=1,smallscalebias=0.0,\
						weighting="natural",uvtaper=self.uvtaper,niter=100000000000,gain=0.05,threshold=self.threshold,cycleniter=2000,cyclefactor=1.0,\
						interactive=False,usemask="user",mask=[mask_str])
						poltclean_dict={'vis':msname,'datacolumn':"corrected",'imsize':self.imsize,'cell':self.cellsize,\
						'stokes':self.stokes,'gridder':gridder,'wprojplanes':wprojplanes,'pblimit':-1,'deconvolver':"multiscale",'scales':self.multiscale_scales,\
						'nterms':1,'smallscalebias':0.0,'weighting':"natural",'uvtaper':self.uvtaper,\
						'niter':100000000000,'gain':0.05,'threshold':self.threshold,'cycleniter':2000,'cyclefactor':1.0,'usemask':"user",'mask':mask_str}
					if do_polcal==True:
						stokes_list=['I','Q','U','V']
						c=0
						for i in range(len(stokes_list)):
							s=self.stokes[i]
							image_pix_sum=imstat(imagename=imagename,mask=mask,stokes=s)['sum'][0]
							residual_pix_sum=imstat(imagename=residual,mask=mask,stokes=s)['sum'][0]
							residual_frac_list.append(residual_pix_sum/image_pix_sum)	
							if residual_pix_sum/image_pix_sum>self.residual_frac:
								self.threshold[i]=str((float(self.threshold[i])/self.start_sigma_list[i])*(self.start_sigma_list[i]-self.sigma_step))+'Jy'
								start_sigma_list[i]=start_sigma_list[i]-self.sigma_step
							else:
								c+=1
						if c==len(stokes_list):
							break
					else:
						image_pix_sum=imstat(imagename=imagename,mask=mask,stokes='I')['sum'][0]
						residual_pix_sum=imstat(imagename=residual,mask=mask,stokes='I')['sum'][0]
						residual_frac_list.append(residual_pix_sum/image_pix_sum)						
						if (residual_pix_sum/image_pix_sum)>self.residual_frac:
							self.threshold[0]=str((float(self.threshold[0])/self.start_sigma_list[0])*(self.start_sigma_list[0]-self.sigma_step))+'Jy'
							start_sigma_list[0]=start_sigma_list[0]-self.sigma_step
						else:
							break
				
				# Exporting images
				##################		
				self.export_images(self,file_str,savedir=savedir,savemodel=savemodel,saveresidual=saveresidual,savecutout=savecutout,cutoutbox=cutoutbox,imageformat=imageformat,\
					other_format=other_format,plotcontour=plotcontour,contour_levels=contour_levels,poltclean_dict=poltclean_dict)
				os.system('rm -rf '+file_str+'*')	
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
	





