import numpy as np
import os,logging
import time as tm
from astropy.io import fits
from casatasks import importuvfits,exportuvfits
from casatools import *
from datetime import datetime 
from . import convertfits as cf
from . import inputs

'''
Code is written by Apurba Bera (NCRA-TIFR)
Wrapper for PAIRCARS is written by Devojyoti Kansabanik, 23 Jan, 2021
'''
os.system('rm -rf casa*log')

class ANKFLAG():
	def __init__(self):
		pathname=os.path.dirname(os.path.realpath(inputs.__file__))
		self.path=pathname
		os.system('rm -rf casa*log')

	def flag_params(self,flagmode):
		'''
		Function to generate flag parameters based on number of channels and timeslices in the dataset
		Parameters :
		flagmode = int
			0 : Single channel and single time slice data
			1 : Multi channel and single time slice data
			2 : Single channel and multi timeslice data
			3 : Multi channel and time slice data
		Return:
		Flag parameters in list
		''' 
		os.system('rm -rf casa*log')
		if flagmode==0:
			return [['vis_ind','mean','median','re',1.8,0.0,1,'',0,0,0.0,0.0],['vis_ind','mean','median','im',1.8,0.0,1,'',0,0,0.0,0.0]]
		elif flagmode==1:
			return [['chan_ind','mean','median','re',1.8,0.0,1,'',0,0,0.0,0.0],['chan_ind','mean','median','im',1.8,0.0,1,'',0,0,0.0,0.0],\
					['vis_ind','mean','median','re',1.8,0.0,1,'',0,0,0.0,0.0],['vis_ind','mean','median','im',1.8,0.0,1,'',0,0,0.0,0.0]]
		elif flagmode==2:
			return [['rec_ind','mean','median','re',1.8,0.0,1,'',0,0,0.0,0.0],['rec_ind','mean','median','im',1.8,0.0,1,'',0,0,0.0,0.0],\
					['vis_ind','mean','median','re',1.8,0.0,1,'',0,0,0.0,0.0],['vis_ind','mean','median','im',1.8,0.0,1,'',0,0,0.0,0.0]]
		else:
			return [['chan_ind','mean','median','re',1.8,0.0,1,'',0,0,0.0,0.0],['chan_ind','mean_rms','median','re',1.8,0.0,1,'',0,0,0.0,0.0],\
					['chan_ind','mean','median','im',1.8,0.0,1,'',0,0,0.0,0.0],['chan_ind','mean_rms','median','im',1.8,0.0,1,'',0,0,0.0,0.0],\
					['rec_ind','mean','median','re',1.8,0.0,1,'',0,0,0.0,0.0],['rec_ind','mean_rms','median','re',1.8,0.0,1,'',0,0,0.0,0.0],\
					['rec_ind','mean','median','im',1.8,0.0,1,'',0,0,0.0,0.0],['rec_ind','mean_rms','median','im',1.8,0.0,1,'',0,0,0.0,0.0],\
					['vis_ind','mean','median','re',1.8,0.0,1,'',0,0,0.0,0.0],['vis_ind','mean_rms','median','re',1.8,0.0,1,'',0,0,0.0,0.0],\
					['vis_ind','mean','median','im',1.8,0.0,1,'',0,0,0.0,0.0],['vis_ind','mean_rms','median','im',1.8,0.0,1,'',0,0,0.0,0.0]]


	def runankflag(self,inpfilename,ANTS,THREADS,nchan,ntime,npols,inp_fileformat='ms',out_fileformat='ms',overwrite=False,automode=True,datacolumn='corrected',\
					flagpars=[],verbose=False,**kwargs):
		'''
		Function to run the aNKflag for full polarization data. Advanced options can be changed from aNKflag installed directory inputs.py file.
		Parameters:
		inpfilename = Input fits file name
		inp_fileformat = 'ms' or 'uvfits' (default = 'ms')
		overwrite = False, overwrite the present data with flagged data or not (default : False)
		ANTS = Total number of antennas (Calculate using CASA msmd tool, do not use the number of antennas shown in listobs)
		THREADS = Number of CPU threads to be used for flagging
		nchan = Number of frequency channels of the dataset
		ntime = Number of time slices of the dataset
		npols = Number of polarisations in dataset
		automode = True, use default flagging settings
		datacolumn = 'CORRECTED_DATA', datacolumn to perform flagging
		flagpars = [], use aNKflag specific flag parameters if automode=False
		verbose = False, verbose output
		###############################
		Advanced options:
		CLEARSCRATCH : True ,Clear the scratch directory (default : True)
		ugrids : Number of grids in u axis (default : 10)
		vgrids : Number of grids in v axis (default : 10)
		plotuv : Plot uv gridding (default : False)
		CONVERTFITS	: Convert FITS to binary (1 for True or 0 for False) (default : 1)
		DOFLAG : Do flagging (1 for True or 0 for False) (default : 1)
		FLAGMODE : 'baseline' / 'uvbin' (default : 'uvbin') 
		BASEFLAGMEAN : [FLAGON,	tolerance_mean,	tolearnce_rms,	min fraction]]	Only for 'baseline' (default : ['mean_rms',	1.5,	1.4,	0.01])
		READBACK : Read back baselines (1 for True or 0 for False) (default : 1)
		SHOWBASE : SHOW baseline stats (1 for True or 0 for False) (default : 1)
		SHOWTF : Show time-frequency plots (1 for True or 0 for False) (default : 1)
		WRITEOUT : Write output (1 for True or 0 for False) (default : 1)
		BLOCKPOW : Power low for Block non-Gaussianity (DON'T CHANGE UNLESS YOU KNOW WHAT IT IS !) (default : 0.8)
		###############################
		Return:
		Flagged output file 
		if overwrite = True, output will be in the same format (ignoring the out_fileformat). 
		Otherwise, for out_fileformat='ms', output will be inpfilename+'.aNKoutms' and inpfilename+'.aNKoutfits' for out_fileformat ='uvfits' if inp_fileformat='uvfits'
		For inp_fileformat='ms' output will be saved in the same ms
		'''
		os.system('rm -rf casa*log')
		pwd=os.getcwd()
		CLEARSCRATCH	=	inputs.CLEARSCRATCH											#	Clear the scratch directory
		ukey			=	inputs.ukey
		vkey			=	inputs.vkey
		wkey			=	inputs.wkey
		ugrids			=	inputs.ugrids
		vgrids			=	inputs.vgrids
		plotuv			=	inputs.plotuv
		CONVERTFITS		=	inputs.CONVERTFITS											#	Convert FITS to binary ?	
		DOFLAG			=	inputs.DOFLAG											#	Do flagging ?
		FLAGMODE		=	inputs.FLAGMODE 									#	'baseline' / 'uvbin' 
		BASEFLAGMEAN	=	inputs.BASEFLAGMEAN		#	[FLAGON,	tolerance_mean,	tolearnce_rms,	min fraction]]			ONLY for 'baseline'
		READBACK		=	inputs.READBACK											#	Read back baselines ?
		SHOWBASE		=	inputs.SHOWBASE											#	SHOW baseline stats ?
		SHOWTF			=	inputs.SHOWTF											#	Show time-frequency plots ?
		WRITEOUT		=	inputs.WRITEOUT											#	Write output ?
		BLOCKPOW		=	inputs.BLOCKPOW											#	Power low for Block non-Gaussianity (DON'T CHANGE UNLESS YOU KNOW WHAT IT IS !)
		os.chdir(self.path)
		LDPATH=str(np.load('LDPATH.npy',allow_pickle=True))
		os.environ['LD_LIBRARY_PATH']=LDPATH
		cwd=os.getcwd()
		os.system('rm -rf casa*log')
		kwords=list(kwargs.keys())
		if len(kwords)!=0:
			inpfil=open('inputs.py','r+')
			lines=inpfil.readlines()
			if 'CLEARSCRATCH' in kwords:
				CLEARSCRATCH=kwargs['CLEARSCRATCH']
				for i in range(len(lines)):
					if 'CLEARSCRATCH' in lines[i]:
						lines[i]='CLEARSCRATCH\t=\t\''+str(kwargs['CLEARSCRATCH'])+'\'\n'	

			if 'ugrids' in kwords:
				ugrids=kwargs['ugrids']
				for i in range(len(lines)):
					if 'ugrids' in lines[i]:
						lines[i]='ugrids\t\t\t=\t'+str(kwargs['ugrids'])+'\n'	 

			if 'vgrids' in kwords:
				vgrids=kwargs['vgrids']
				for i in range(len(lines)):
					if 'vgrids' in lines[i]:
						lines[i]='vgrids\t\t\t=\t'+str(kwargs['vgrids'])+'\n'	 
		
			if 'plotuv' in kwords:
				plotuv=kwargs['plotuv']
				for i in range(len(lines)):
					if 'plotuv' in lines[i]:
						lines[i]='plotuv\t\t\t=\t'+str(kwargs['plotuv'])+'\n'
		
			if 'CONVERTFITS' in kwords:
				CONVERTFITS=kwargs['CONVERTFITS']
				for i in range(len(lines)):
					if 'CONVERTFITS' in lines[i]:
						lines[i]='CONVERTFITS\t\t=\t'+str(kwargs['CONVERTFITS'])+'\n'	 

			if 'DOFLAG' in kwords:
				DOFLAG=kwargs['DOFLAG']
				for i in range(len(lines)):
					if 'DOFLAG' in lines[i]:
						lines[i]='DOFLAG\t\t\t=\t'+str(kwargs['DOFLAG'])+'\n'	 
		
			if 'FLAGMODE' in kwords:
				FLAGMODE=kwargs['FLAGMODE']
				for i in range(len(lines)):
					if 'FLAGMODE' in lines[i]:
						lines[i]='FLAGMODE\t\t=\t\''+str(kwargs['FLAGMODE'])+'\'\n'	 
			
			if 'BASEFLAGMEAN' in kwords:
				BASEFLAGMEAN=kwargs['BASEFLAGMEAN']
				for i in range(len(lines)):
					if 'BASEFLAGMEAN' in lines[i]:
						lines[i]='BASEFLAGMEAN\t=\t'+str(kwargs['BASEFLAGMEAN'])+'\n'	

			if 'READBACK' in kwords:
				READBACK=kwargs['READBACK']
				for i in range(len(lines)):
					if 'READBACK' in lines[i]:
						lines[i]='READBACK\t\t=\t'+str(kwargs['READBACK'])+'\n' 

			if 'SHOWBASE' in kwords:
				SHOWBASE=kwargs['SHOWBASE']
				for i in range(len(lines)):
					if 'SHOWBASE' in lines[i]:
						lines[i]='SHOWBASE\t\t=\t'+str(kwargs['READBACK'])+'\n' 
		
			if 'SHOWTF' in kwords:
				SHOWTF=kwargs['SHOWTF']
				for i in range(len(lines)):
					if 'SHOWTF' in lines[i]:
						lines[i]='SHOWTF\t\t\t=\t'+str(kwargs['SHOWTF'])+'\n' 

			if 'WRITEOUT' in kwords:
				WRITEOUT=kwargs['WRITEOUT']
				for i in range(len(lines)):
					if 'WRITEOUT' in lines[i]:
						lines[i]='WRITEOUT\t\t=\t'+str(kwargs['WRITEOUT'])+'\n' 

			if 'BLOCKPOW' in kwords:
				BLACKPOW=kwargs['BLACKPOW']
				for i in range(len(lines)):
					if 'BLOCKPOW' in lines[i]:
						lines[i]='BLOCKPOW\t\t=\t'+str(kwargs['BLOCKPOW'])+'\n' 

			inpfil.seek(0)
			inpfil.writelines(lines)
			inpfil.close()

		inp_filepath=os.path.dirname(os.path.realpath(inpfilename))
		inpfile_basename=os.path.basename(inpfilename)
		if verbose==True:
			logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
			logger = logging.getLogger('ankflag_logger')
			logger.setLevel(logging.DEBUG)
			if os.path.isfile(inp_filepath+'/aNKflagger.log')==True:
				os.system('rm -rf '+inp_filepath+'/aNKflagger.log')
			fh = logging.FileHandler(inp_filepath+'/aNKflagger.log')
			fh.setLevel(logging.DEBUG)
			logger.addHandler(fh)
			logger.info('Input file name : '+inpfilename+'\n')
			logger.info('Starting aNKflagger..........\n')
			logger.info('Flagging datacolumn : '+datacolumn+'\n')
			if datacolumn!='corrected' and datacolumn!='data':
				logger.error('Datacolumn is not either data or corrected.\n')

		if inp_fileformat=='uvfits':
			try:
				header=fits.getheader(inpfilename)
				inpfits=inpfilename
			except:
				if os.path.isdir(inpfilename):
					try:
						if os.path.isfile(inp_filepath+'/'+inpfile_basename.split('.ms')[0]+'.fits'):
							os.system('rm -rf '+inp_filepath+'/'+inpfile_basename.split('.ms')[0]+'.fits')
						exportuvfits(vis=inpfilename,fitsfile=inp_filepath+'/'+inpfile_basename.split('.ms')[0]+'.fits',datacolumn=datacolumn)	
						if verbose:
							logger.info('exportuvfits(vis='+inpfilename+',fitsfile='+inp_filepath+'/'+inpfile_basename.split('.ms')[0]+'.fits)\n')
						inpfits=inp_filepath+'/'+inpfile_basename.split('.ms')[0]+'.fits'
						inp_fileformat='ms'
					except:
						if verbose:
							logger.error('Input file format is not either uvfits or ms.\n')
						os._exit(1)
				else:
					if verbose:
						logger.error('Input file format is not either uvfits or ms.\n')
					os._exit(1)

		elif inp_fileformat=='ms':
			if os.path.isdir(inpfilename):
				try:
					if os.path.isfile(inp_filepath+'/'+inpfile_basename.split('.ms')[0]+'.fits'):
						os.system('rm -rf '+inp_filepath+'/'+inpfile_basename.split('.ms')[0]+'.fits')
					exportuvfits(vis=inpfilename,fitsfile=inp_filepath+'/'+inpfile_basename.split('.ms')[0]+'.fits',datacolumn=datacolumn)
					if verbose:
						logger.info('exportuvfits(vis='+inpfilename+',fitsfile='+inp_filepath+'/'+inpfile_basename.split('.ms')[0]+'.fits)\n')
					inpfits=inp_filepath+'/'+inpfile_basename.split('.ms')[0]+'.fits'
				except:
					if verbose:
						logger.error('Input file format is not either uvfits or ms.\n')
					os._exit(1)
			else:
				try:
					header=fits.getheader(inpfilename)
					inpfits=inpfilename
					inp_fileformat='uvfits'
				except:
					logger.error('Input file format is not either uvfits or ms.\n')
					os._exit(1)

		outfits=inp_filepath+'/temp_ankflag.uvfits'
		if os.path.isfile(outfits):
			os.system('rm -rf '+outfits)

		scratchdir=cwd+'/scratch/'
		if (CLEARSCRATCH):
			if os.path.isdir(scratchdir):
				if verbose:
					print('\nClearing scratch directory....\n')
				os.system('rm -rf '+scratchdir+'*')
			else:
				os.system('mkdir '+scratchdir)
		# Choosing flag parameters based on number of channels and time slices in the dataset
		if automode==True or len(flagpars)==0:
			if nchan==1 and ntime==1:
				FLAGPARS=self.flag_params(0)
				flagmode=0
			elif nchan!=1 and ntime==1:
				FLAGPARS=self.flag_params(1)
				flagmode=1
			elif nchan==1 and ntime!=1:
				FLAGPARS=self.flag_params(2)
				flagmode=2
			else:
				FLAGPARS=self.flag_params(3)
				flagmode=3
		else:
			FLAGPARS=flagpars
			flagmode=4

		if verbose==True:
			logger.info('Flag mode : '+str(flagmode)+'\n')
			logger.info(str(FLAGPARS)+'\n')

		exmode		=	['baseline', 'uvbin']
		flagwhat	=	['vis_ind', 'chan_ind', 'rec_ind', 'vis_block', 'chan_block', 'rec_block']
		flagon		=	['mean', 'rms', 'mean_rms']
		statused	=	['median', 'mean']
		datype		=	['re', 'im', 'am', 'ph']
		blkorder	=	['ascending', '', 'descending']
		if verbose:
			logger.info('Total flagging rounds	=	%d'%len(FLAGPARS))
		#	--------------		Convert inputs to numbers
		flagparfile	=	open(scratchdir+'flagpars.pars','w')
		flagparfile.write('%d	%d	%d	%d	%d	'%(ANTS, exmode.index(FLAGMODE), ugrids, flagon.index(BASEFLAGMEAN[0])+1, vgrids))
		flagparfile.write('%f	%f	%d	%d	%d	%f	%f\n'%(BASEFLAGMEAN[1], BASEFLAGMEAN[2], len(FLAGPARS), THREADS, WRITEOUT, BLOCKPOW, BASEFLAGMEAN[3]))

		for flpar in FLAGPARS:
			flagparfile.write('%d	%d	%d	%d	%d	'%(flagwhat.index(flpar[0]), flagon.index(flpar[1])+1, statused.index(flpar[2]), \
								datype.index(flpar[3]), -(blkorder.index(flpar[7])-1)))
			flagparfile.write('%f	%f	%d	%d	%d	%f	%f\n'%(flpar[4], flpar[5], flpar[6]+1, flpar[8], flpar[9], flpar[10], flpar[11]))

		flagparfile.close()

		#	----------------------------	Convert FITS to binary files
		start0	=	tm.time()	

		if (CONVERTFITS):

			infile		=	fits.open(inpfits)
			data		=	infile[0].data
			
			if (FLAGMODE==exmode[1]):		
				cf.uvfitstobinary(data,scratchdir,ugrids,vgrids,plotuv,npols,verbose)

			elif (FLAGMODE==exmode[0]):		
				cf.baselinestobinary(ANTS,data,scratcOUThdir,npols,verbose)
					
			else:
				if verbose:
					logger.info("Unknown flagging mode !!!!			Please tell me how to execute it ........\n")
			
			infile.close()	
			if verbose:
				logger.info("Convertion done in 		%d seconds\n"%(tm.time()-start0))
			
		#	------------------------------		Flag data
		start1	=	tm.time()

		if (DOFLAG):
			if verbose==False:
				status	=	os.system('./ankflag %d'%npols+' > ankflag.out')	
				#print("\nFlagging done in 		%d seconds\n"%(tm.time()-start1))
			else:
				status	=	os.system('./ankflag %d'%npols)	
				logger.info("Flagging done in 		%d seconds\n"%(tm.time()-start1))
				
			
		#	------------------------------		Convert back binary files to FITS	

		if (READBACK):

			infile2		=	fits.open(inpfits)		
			data2		=	infile2[0].data	

			if (FLAGMODE==exmode[1]):
				bintofits	=	cf.uvfitsfrombinary(data2,scratchdir,ugrids,vgrids,npols,verbose)

			elif (FLAGMODE==exmode[0]):
				cf.baselinesfrombinary(ANTS,data2,scratchdir,npols,verbose)

			#	---------------------------------------------------------------------
			#					Plot Baselines 
			#	---------------------------------------------------------------------

			if (SHOWBASE):

				infile	=	fits.open(inpfits)
				data	=	infile[0].data

				blid	=	[]
				for a in range (1,ANTS):
					for b in range (a+1,ANTS+1):
						blid.append([a,b,256*a+b])
				blid	=	np.array(blid)
				nbase	=	len(blid)
				if verbose:
					logger.info('Ideally total baselines =	%d'%nbase)

				flaggingstatus	=	[]
				for i in range (0,nbase):	
					flaggingstatus.append(cf.showbasecomparison(blid[i],data,data2,SHOWTF,npols,verbose))

				flaggingstatus	=	np.array(flaggingstatus)
				avgflag			=	np.mean(flaggingstatus,axis=0)
				
				if verbose:
					logger.info('Average flagging fraction	'),
					for p in range (0,npols):
						logger.info('%.3f %.3f	'%(avgflag[p],avgflag[p+npols])),
					logger.info('\n')
					
					logger.info('Flagged data			'),
					for p in range (0,npols):
						logger.info('%.3f	'%((avgflag[p+npols]-avgflag[p])/(1.0-avgflag[p]))),
					logger.info('\n')

				infile.close()

			if (WRITEOUT):
				infile2.writeto(outfits,output_verify='warn',overwrite=True)

			infile2.close()
			if verbose:
				logger.info("Everything done in 		%d seconds\n"%(tm.time()-start0))
			#	----------------------------------------------------------------------------

		if (WRITEOUT):
			if overwrite==True and inp_fileformat=='uvfits' and out_fileformat=='uvfits':
				os.system('mv '+outfits+' '+inpfits)
				finalout=inpfits
			else:
				if (inp_fileformat=='uvfits' and out_fileformat=='uvfits') or (inp_fileformat=='ms' and out_fileformat=='uvfits'):
					finalout=inpfits+'.aNKoutfits'
					os.system('mv '+outfits+' '+inpfits+'.aNKoutfits')
					if verbose:
						logger.info('Final outputfile : '+inpfits+'.aNKoutfits\n')
				elif inp_fileformat=='uvfits' and out_fileformat=='ms':
					if os.path.isdir(inpfits+'.aNKoutms'):
						os.system('rm -rf '+inpfits+'.aNKoutms')
					importuvfits(fitsfile=outfits,vis=inpfits+'.aNKoutms')
					finalout=inpfits+'.aNKoutms'
					if verbose:
						logger.info('Final outputfile : '+inpfits+'.aNKoutms\n')
				elif inp_fileformat=='ms' and out_fileformat=='ms':
					afank=agentflagger()
					afank.open(inpfilename)
					versionlist=afank.getflagversionlist()
					if len(versionlist)!=0:
						for version_name in versionlist:
							if 'aNKflag' in version_name:
								version_num=int(version_name.split(':')[0].split(' ')[0].split('_')[-1])+1
							else:
								version_num=1
					else:
						version_num=1
					now = datetime.now()
					dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
					afank.saveflagversion('aNKflag_'+str(version_num),'Flags autosave on '+dt_string)
					if os.path.isdir(inp_filepath+'/temp_aNKflag.ms'):
						os.system('rm -rf '+inp_filepath+'/temp_aNKflag.ms')
					importuvfits(fitsfile=outfits,vis=inp_filepath+'/temp_aNKflag.ms')
					tbank=table()
					tbank.open(inp_filepath+'/temp_aNKflag.ms')
					data=tbank.getcol('DATA')
					tbank.close()
					os.system('rm -rf '+inp_filepath+'/temp_aNKflag.ms')
					tbank.open(inpfilename,nomodify=False)
					if datacolumn=='corrected':
						try:
							tbank.putcol('CORRECTED_DATA',data)
						except:
							tbank.putcol('DATA',data)
					elif datacolumn=='data':
						tbank.putcol('DATA',data)
					tbank.flush()
					tbank.close()
					if verbose:
						logger.info('Final outputfile : '+inpfilename+'\n')
					os.system('rm -rf '+outfits)
					finalout=inpfilename

			if os.path.isfile(outfits):
				os.system('rm -rf '+outfits)
	
		else:
			finalout=''
		if (CLEARSCRATCH):
			if verbose:
				print('Clearing scratch directory....\n')
			os.system('rm '+scratchdir+'*')
		os.system('rm -rf casa*log')
		os.chdir(pwd)
		os.system('rm -rf casa*log')
		os.unsetenv('LD_LIBRARY_PATH')
		return finalout






























