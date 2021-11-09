import os
import numpy as np,psutil,copy,logging
from . import access_ms as am
from datetime import datetime 
try:
	from aNKflag import runank
except:
	pass
from casatools import msmetadata,table,measures,quanta,agentflagger,image,calibrater,ms
from casatasks import *
from paircars.access_ms import *
from astropy.io import fits
os.system('rm -rf casa*log')
'''
Code is written by Devojyoti Kansabanik, 26 Jan ,2021
'''

def flag_MWA_coarse(msname,edgewidth=80,do_flag=True,force=False,flagbackup=True):
	'''
	A function to generate the list of coarse-channels edges and the 
	central channel in each coarse channel to be flagged.

	Parameters
	----------
	msname : str 
		Name of the masurement set
	edgewidth : float 
		Flag edge channels width in kHz
	do_flag : bool 
		If true flag edge channels, otherwise return only the good channels list
	force : bool 
		If True flag again if even if the flag keyword is in header
	flagbackup : bool
		Keep flag backup or not
	Returns
	-------
	str
		Unflagged channels
	dict
		One central channel per coarse channel
	'''
	# Function is written by Divya Oberoi, 07 Apr, 2016
	# Modified by Devojyoti Kansabanik, 23 Jan, 2021
	AM=am.AccessMS(msname)
	ncoarse_chan=AM.calc_ncoarse()
	freqres=AM.calc_freqres()
	nchan_per_coarse=int(1.28*10**3/freqres)
	M = int(edgewidth/freqres) # No. of channels to be flagged at the start of the coarse channel
	N = int(edgewidth/freqres) # No. of channels to be flagged at the tail of the coarse channel
	a = int(nchan_per_coarse/2)	# The central channel which occassionally shows the DC spike
	CHAN_FLAG_STR='0:'
	CHAN_UNFLAG_STR='0:'
	i = 0
	ch0 = 0
	ch1 = nchan_per_coarse
	channels_per_coarse={}
	if ncoarse_chan!=0:
		while i < ncoarse_chan:
			# The 0th coarse channel requires special treatment (one less ';')
			if i == 0:
				CHAN_FLAG_STR=CHAN_FLAG_STR+str(ch0)+'~'+str(ch0+M)+';'+str(a)+';'+str(ch1-N)+'~'+str(ch1-1)	
				CHAN_UNFLAG_STR+=str(ch0+M+1)+'~'+str(a-1)+';'+str(a+1)+'~'+str(ch1-(N+1)-1)+';'+str(ch1-(N+1)+1)+'~'+str(ch1-1)		
			else:
				CHAN_FLAG_STR=CHAN_FLAG_STR+';'+str(ch0)+'~'+str(ch0+M)+';'+str(a)+';'+str(ch1-N)+'~'+str(ch1-1)
				CHAN_UNFLAG_STR+=';'+str(ch0+M+1)+'~'+str(a-1)+';'+str(a+1)+'~'+str(ch1-(N+1)-1)+';'+str(ch1-(N+1)+1)+'~'+str(ch1-1)
			channels_per_coarse[i]=a-1
			a = a+nchan_per_coarse
			ch0 = ch0 + nchan_per_coarse
			ch1 = ch1 + nchan_per_coarse
			i = i + 1
		code=vishead(vis=msname,mode='get',hdkey='fld_code')[0][0]
		code_list=code.split(',')
		if ('C_FLAG_'+str(edgewidth) not in code_list and do_flag==True) or (force==True and do_flag==True):
			print ('Flagging coarse channel edges and central DC-spike channels:'+CHAN_FLAG_STR+'\n')
			flagdata(vis=msname,spw=CHAN_FLAG_STR,mode='manual',flagbackup=flagbackup)
			flagdata(vis=msname,autocorr=True,flagbackup=flagbackup)
			if (force==True and do_flag==True and 'C_FLAG_'+str(edgewidth) not in code_list) or (force==False and do_flag==True):
				if len(code_list)==1 and code_list[0]=='':
					code+='C_FLAG_'+str(edgewidth)
				else:
					code+=',C_FLAG_'+str(edgewidth)
			vishead(vis=msname,mode='put',hdkey='fld_code',hdvalue=np.array([code]))
	else:
		if ncoarse_chan==0:
			print ('Number of coarse channel is less than 1. No coarse channel flagging is required.\n')
		else:
			print ('Coarse channel edges are already flagged.\n')
	return CHAN_UNFLAG_STR,channels_per_coarse

def flag_MWA_quack(msname,metafits,quacktime=0.0,flagbackup=False,force=False):
	'''
	Function to flag MWA Quack times

	Parameters
	----------
	msname : str 
		Name of the measurment set
	metafits : str 
		Name of the metafits file
	quacktime : float
		Quack time
	flagbackup : bool
		Backup flags or not
	force : bool
		If True flag again if even if the flag keyword is in header
	Returns
	-------
	float
		Quack time
	'''
	AM=AccessMS(msname)
	timeres=AM.calc_timeres()
	if timeres<1:
		timeres=1
	if quacktime==0:
		quacktime=float(fits.getheader(metafits)['QUACKTIM'])+timeres
	code=vishead(vis=msname,mode='get',hdkey='fld_code')[0][0]
	code_list=code.split(',')
	if 'QUACK_'+str(quacktime) not in code_list or force==True:
		flagdata(vis=msname,mode='quack',quackinterval=quacktime,quackmode='beg',flagbackup=flagbackup)
		flagdata(vis=msname,mode='quack',quackinterval=quacktime,quackmode='endb',flagbackup=flagbackup)
		if len(code_list)==1 and code_list[0]=='':
			code+='QUACK_'+str(quacktime)
		else:
			code+=',QUACK_'+str(quacktime)
		vishead(vis=msname,mode='put',hdkey='fld_code',hdvalue=np.array([code]))
	else:
		print ('Quack times are already flagged.\n')
	return quacktime

def do_uvsub_ankflag(msname,model='',nthread=0,verbose=False,flagbackup=True,extendpols=False,chantime_minfrac=0.5,casaflag=''): 
	'''
	Perform flagging on uv sub data using aNKflagger

	Parameters
	----------
	msname : str 
		Name of the measurement set
	model : str 
		Model image name, Keep blank if model is already in modelcolumn
	nthread : int 
		Number of CPU threads to be used by aNKflag. If 0, it will use 25% of the total available CPU threads.
	verbose : bool 
		If True keep all records
	flagbackup : bool
		Keep flagbackup
	extendpols : bool 
		Extend flag if one polarisation is flagged
	chantime_minfrac : float 
		Minimum fraction of data flagged for a single channel and time to extend the flag
	casaflag : str 
		Perform CASA flags ('tfcrop' or 'rflag')
	Returns
	-------
	str
		Flagged measurement set name
	'''
	print ('Using aNKflagger....\n')
	mdflag=msmetadata()
	tbflag=table()
	mdflag.open(msname)
	nants=mdflag.nantennas()
	nchan=mdflag.nchan(0)
	nbaseline=mdflag.nbaselines(ac=True)
	npols=mdflag.ncorrforpol()[0]
	ntimes=mdflag.timesforfield(0).size
	total_cpus=psutil.cpu_count()
	if nthread==0:
		available_cpus=int((total_cpus-(psutil.cpu_percent()*total_cpus)))
		if available_cpus<=0:
			available_cpus=1
	else:
		available_cpus=nthread
	mdflag.close()
	if model!='':
		delmod(vis=msname,scr=True)
		ft(vis=msname,model=model,usescratch=True)
	uvsub(vis=msname,reverse=False)
	ankflagger=runank.ANKFLAG()
	if verbose==True:
		print ('ankflagger.runankflag('+msname+','+str(nants)+','+str(available_cpus)+','+str(nchan)+','+str(nbaseline)+','+str(ntimes)+','+str(npols)+',inp_fileformat=\'ms\','+\
		'out_fileformat=\'ms\',automode=True,datacolumn=\'corrected\',verbose=verbose,flagbackup='+str(flagbackup)+',extendpols='+str(extendpols)+\
		',extendchan=0.0,extendtime=0.0,chantime_flagfrac='+str(chantime_minfrac)+')\n')
	outfile=ankflagger.runankflag(msname,nants,available_cpus,nchan,nbaseline,ntimes,npols,inp_fileformat='ms',out_fileformat='ms',automode=True,datacolumn='corrected',\
			verbose=verbose,flagbackup=flagbackup,extendpols=extendpols,extendchan=0.0,extendtime=0.0,chantime_flagfrac=chantime_minfrac)
	if casaflag!='' and (casaflag=='rflag' or casaflag=='tfcrop'):
		if casaflag=='rflag' and (nchan>=30 or ntimes>=30):
			if nchan<=30:
				freqdevscale=100000
				timedevscale=5.0
			if ntimes<=30:
				timedevscale=100000
				freqdevscale=5.0
			else:
				freqdevscale=5.0
				timedevscale=5.0
			if verbose:
				print ('flagdata(vis=\''+msname+'\',mode=\'rflag\',datacolumn=\'corrected\',extendflags=False,freqdevscale='+str(freqdevscale)+',timedevscale='+str(timedevscale)+\
						',flagbackup='+str(flagbackup)+')\n')
			try:
				flagdata(vis=msname,mode=casaflag,datacolumn='corrected',extendflags=False,freqdevscale=freqdevscale,timedevscale=timedevscale,flagbackup=flagbackup)
			except Exception as e:
				print ('Error occured during CASA flag : '+str(e)+'\n')
				pass
		elif casaflag=='tfcrop' and (nchan>=30 or ntimes>=30):
			if nchan<=30:
				freqcutoff=100000
				timecutoff=5.0
				flagdimension='time'
			if ntimes<=30:
				flagdimension='freq'
				timecutoff=100000
				freqcutoff=5.0
			else:
				flagdimension='freqtime'
				timecutoff=5.0
				freqcutoff=5.0
			if verbose:
				print ('flagdata(vis=\''+msname+'\',mode=\'tfcrop\',datacolumn=\'corrected\',extendflags=False,freqcutoff='+str(freqcutoff)+',timecutoff='+str(timecutoff)+\
						',flagdimension=\''+flagdimension+'\',flagbackup='+str(flagbackup)+')\n')
			try:
				flagdata(vis=msname,mode=casaflag,datacolumn='corrected',extendflags=False,freqcutoff=freqcutoff,timecutoff=timecutoff,flagdimension=flagdimension,flagbackup=flagbackup)
			except Exception as e:
				print ('Error occured during CASA flag : '+str(e)+'\n')
				pass
		else:
			print ('Number of time and frequency slices is less than 30. Do not perform CASA flag.\n')
	else:
		print ('CASA flag should be either tfcrop or rflag.\n')
	uvsub(vis=msname,reverse=True)
	os.system('rm -rf aNKflagger.log casa*log ankflag.out '+msname.split('.ms')[0]+'.fits')
	return outfile

def do_uvsub_flagger(msname,model='',mode='',rmsthresh=[],flagbackup=True):
	'''
	Flagger on residual data

	Parameters
	----------
	msname : str 
		Name of the measurement set
	model : str 
		Name of the model image, keep blank if model is alrealy in modelcolumn
	rmsthresh : list
		List of rms threshold 
	flagbackup : bool
		Keep flagbackup
	Returns
	-------
	int
		New number of flaggs
	'''
	tb=table()
	md=msmetadata()
	print ('Using uvsub_flagged.....\n')
	AM=am.AccessMS(msname)
	unflag_chan=AM.get_unflag_chan()
	antenna=''
	unflag_spws='0:'
	nant=AM.get_num_antenna()
	if len(unflag_chan)>0:
		for i in unflag_chan:
			unflag_spws+=str(i)+';'
		unflag_spws=unflag_spws[:-1]
		for i in range(nant):
			antenna+=str(i)+','
		antenna=antenna[:-1]
	# Keeping flag backup
	if flagbackup==True:
		af=agentflagger()
		af.open(msname)
		versionlist=af.getflagversionlist()
		if len(versionlist)!=0:
			for version_name in versionlist:
				if mode=='applycal':
					if 'applycal' in version_name:
						version_num=int(version_name.split(':')[0].split(' ')[0].split('_')[-1])+1
					else:
						version_num=1
				elif mode=='flagdata':
					if 'flagdata' in version_name:
						version_num=int(version_name.split(':')[0].split(' ')[0].split('_')[-1])+1
					else:
						version_num=1
				else:
					if mode=='':
						mode='uvsub_flag'
					if mode in version_name:
						version_num=int(version_name.split(':')[0].split(' ')[0].split('_')[-1])+1
					else:
						version_num=1
		else:
			version_num=1
		now = datetime.now()
		dt_string = now.strftime("%Y-%m-%d %H:%M:%S")
		if mode=='applycal':
			af.saveflagversion('applycal_'+str(version_num),'Flags autosave on '+dt_string)
		elif mode=='flagdata':
			af.saveflagversion('flagdata_'+str(version_num),'Flags autosave on '+dt_string)
		else:
			af.saveflagversion(mode+'_'+str(version_num),'Flags autosave on '+dt_string)
		af.done()
	if model!='':
		delmod(vis=msname,scr=True)
		ft(vis=msname,model=model,usescratch=True)
	uvsub(vis=msname,reverse=False)
	tb.open(msname,nomodify=False)
	cor_col=tb.getcol('CORRECTED_DATA')
	cor_col_copy=copy.deepcopy(cor_col)
	flag_col=tb.getcol('FLAG')
	md.open(msname)
	nchan=md.nchan(0)
	npols=md.ncorrforpol()[0]
	md.close()
	cor_col[flag_col]=np.nan
	old_flags=np.sum(flag_col)
	for nsigma in rmsthresh:
		print ('Flagging in '+str(nsigma)+' sigma threshold.\n')
		flag_count=True
		flag_round=0
		while flag_count==True:
			flag_size=0
			for chan in range(nchan):			
				for pol in range(npols):
					data=cor_col[pol,chan,:]
					pos_flag=np.append(np.where(np.real(data)>np.nanstd(np.real(data))*nsigma),np.where(np.imag(data)>np.nanstd(np.imag(data))*nsigma))
					flag_size+=pos_flag.size
					if pos_flag.size!=0:
						flag_col[pol,chan,pos_flag]=True
						cor_col[flag_col]=np.nan
			for pol in range(npols):
				data=cor_col[pol,:,:]
				pos_flag=np.append(np.where(np.real(data)>np.nanstd(np.real(data))*nsigma),np.where(np.imag(data)>np.nanstd(np.imag(data))*nsigma))
				flag_size+=pos_flag.size
				if pos_flag.size!=0:
					flag_col[pol,pos_flag[0],pos_flag[1]]=True
					cor_col[flag_col]=np.nan
			flag_round+=1
			if flag_size==0 or flag_round>500:
				flag_count=False
	tb.putcol('FLAG',flag_col)
	tb.flush()
	tb.close()
	uvsub(vis=msname,reverse=True)
	new_flags=np.sum(flag_col)
	os.system('rm -rf casa*log')
	return (new_flags-old_flags)

def get_bad_ants(msname,flagfrac=1.0):
	'''
	Function to get the antennas for which a certain amount of data are flagged

	Parameters
	----------
	msname : str 
		Name of the measurement set
	flagfrac : float 
		Fraction of data flagged to consider the antennas as bad (default : 1.0)
	Returns
	--------
	str
		Bad antenna strings
	'''
	mstool=ms()
	msmd=msmetadata()
	msmd.open(msname)
	nant=msmd.nantennas()
	msmd.close()
	ant_flag_percent={}
	for i in range(nant):
		mstool.open(msname)
		mstool.select({'antenna1':i})
		flag=mstool.getdata('FLAG')['flag']
		mstool.close()
		mstool.open(msname)
		mstool.select({'antenna2':i})
		flag=np.append(flag,mstool.getdata('FLAG')['flag'])
		mstool.close()
		ant_flag_percent[i]=(np.sum(flag.flatten())/float(len(flag.flatten())))
	bad_ants=''
	for i in range(nant):
		if ant_flag_percent[i]>flagfrac:
			bad_ants+=str(i)+','
	bad_ants=bad_ants[:-1]
	return bad_ants

def flag_zeros(msname,flagbackup=False,force=False):
	'''
	Function to flag zero data

	Parameters
	----------
	msname : str
		Name of the measurement set
	flagbackup : bool
		Backup flags or not
	force : bool
		Force to flag even if the header saying flagging is already done
	'''
	code=vishead(vis=msname,mode='get',hdkey='fld_code')[0][0]
	code_list=code.split(',')
	if 'ZEROFLAG' not in code_list or force==True:
		flagdata(vis=msname,mode='clip',clipzeros=True,flagbackup=flagbackup)
		if len(code_list)==1 and code_list[0]=='':
			code+='ZEROFLAG'
		else:
			code+=',ZEROFLAG'
		vishead(vis=msname,mode='put',hdkey='fld_code',hdvalue=np.array([code]))
	else:
		print ('Zero data flags are already flagged.\n')
	return

def casa_autoflag(msname,mode='rflag',datacolumn='corrected',sigma_thresh=5.0,flagbackup=False,verbose=False):
	'''
	Function to perform CASA flag 
	Parameters
	----------
	msname : str
		Name of the measurement set
	mode : str
		CASA flag mode, 'rflag' or 'tfcrop'
	datacolumn : str
		Datacolumn to flag, 'data', 'corrected' or 'residual'
	sigma_thresh : float
		Flagging threshold
	flagbackup : bool
		Backup flags or not
	verbose : bool
		Verbose output
	Returns
	-------
	int
		Success code, 0 or 1
	'''
	mdflag=msmetadata()
	mdflag.open(msname)
	nants=mdflag.nantennas()
	nchan=mdflag.nchan(0)
	nbaseline=mdflag.nbaselines(ac=True)
	npols=mdflag.ncorrforpol()[0]
	ntimes=mdflag.timesforfield(0).size
	mdflag.close()
	if mode!='' and (mode=='rflag' or mode=='tfcrop'):
		if mode=='rflag' and (nchan>=10 or ntimes>=10):
			if nchan<=10:
				freqdevscale=100000
				timedevscale=sigma_thresh
			if ntimes<=10:
				timedevscale=100000
				freqdevscale=sigma_thresh
			else:
				freqdevscale=sigma_thresh
				timedevscale=sigma_thresh
			if verbose:
				print ('flagdata(vis=\''+msname+'\',mode=\''+str(mode)+'\',datacolumn=\''+datacolumn+'\',extendflags=False,freqdevscale='+\
				str(freqdevscale)+',timedevscale='+str(timedevscale)+',flagbackup='+str(flagbackup)+')\n')
			try:
				flagdata(vis=msname,mode=mode,datacolumn=datacolumn,extendflags=False,freqdevscale=freqdevscale,timedevscale=timedevscale,flagbackup=flagbackup)
				return 0
			except Exception as e:
				print ('Error occured during CASA flag : '+str(e)+'\n')
				return 1
				pass
		elif mode=='tfcrop' and (nchan>=10 or ntimes>=10):
			if nchan<=10:
				freqcutoff=100000
				timecutoff=sigma_thresh
				flagdimension='time'
			if ntimes<=10:
				flagdimension='freq'
				timecutoff=100000
				freqcutoff=sigma_thresh
			else:
				flagdimension='freqtime'
				timecutoff=sigma_thresh
				freqcutoff=sigma_thresh
			if verbose:
				print ('flagdata(vis=\''+msname+'\',mode=\''+str(mode)+'\',datacolumn=\''+datacolumn+'\',extendflags=False,freqcutoff='+str(freqcutoff)+',timecutoff='+str(timecutoff)+\
						',flagdimension=\''+flagdimension+'\',flagbackup='+str(flagbackup)+')\n')
			try:
				flagdata(vis=msname,mode=mode,datacolumn=datacolumn,extendflags=False,freqcutoff=freqcutoff,timecutoff=timecutoff,flagdimension=flagdimension,flagbackup=flagbackup)
				return 0
			except Exception as e:
				print ('Error occured during CASA flag : '+str(e)+'\n')
				return 1
		else:
			print ('Number of time and frequency slices is less than 10. Do not perform CASA flag.\n')
			return 1
	else:
		print ('CASA flag should be either tfcrop or rflag.\n')
		return 1





















