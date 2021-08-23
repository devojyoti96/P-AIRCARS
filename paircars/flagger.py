import numpy as np,os,psutil,copy,logging
from . import access_ms as am
from datetime import datetime 
try:
	from aNKflag import runank
except:
	pass
from casatools import *
from casatasks import *
'''
Code is written by Devojyoti Kansabanik, 26 Jan ,2021
'''

def flag_MWA_coarse(msname,edgewidth=80,do_flag=True,force=False,flagbackup=True):
	'''
	A function to generate the list of coarse-channels edges + the 
	central channel in each coarse channel to be flagged.
	Parameters:
	msname= Name of the masurement set
	edgewidth = Flag edge channels width in kHz
	do_flag = True, If true flag edge channels, otherwise return only the good channels list
	force = False, If True flag again if even if the flag keyword is in header
	flagbackup = True, keep flag backup or not
	Return:
	Unflagged channels, One central channel per coarse channel
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


def do_uvsub_ankflag(msname,model='',nthread=0,verbose=False,flagbackup=True): 
	'''
	Perform flagging on uv sub data using aNKflagger
	Parameters:
	msname = Name of the measurement set
	model = Model image name, Keep blank if model is already in modelcolumn
	nthread = Number of CPU threads to be used by aNKflag. If 0, it will use 25% of the total available CPU threads.
	verbose = False, If True keep all records
	flagbackup = True, Keep flagbackup
	Return:
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
	print ('ankflagger.runankflag('+msname+','+str(nants)+','+str(available_cpus)+','+str(nchan)+','+str(nbaseline)+','+str(ntimes)+','+str(npols)+',inp_fileformat=\'ms\','+\
		'out_fileformat=\'ms\',automode=True,datacolumn=\'corrected\',verbose=verbose,flagbackup=True,extendpols=True,extendchan=0.0,extendtime=0.0,chantime_flagfrac=0.5)\n')
	outfile=ankflagger.runankflag(msname,nants,available_cpus,nchan,nbaseline,ntimes,npols,inp_fileformat='ms',out_fileformat='ms',automode=True,datacolumn='corrected',\
			verbose=verbose,flagbackup=True,extendpols=False,extendchan=0.0,extendtime=0.0,chantime_flagfrac=0.5)
	uvsub(vis=msname,reverse=True)
	os.system('rm -rf aNKflagger.log casa*log ankflag.out '+msname.split('.ms')[0]+'.fits')
	return outfile

def do_uvsub_flagger(msname,model='',mode='',rmsthresh=[],flagbackup=True):
	'''
	Flagger on residual data
	Parameters:
	msname = Name of the measurement set
	model = Name of the model image, keep blank if model is alrealy in modelcolumn
	rmsthresh = [], rms threshold lists
	flagbackup = True, keep flagbackup
	Return:
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
		print ('flagdata(vis=\''+msname+'\',antenna=\''+antenna+'\',spw=\''+unflag_spws+'\',mode=\'unflag\')\n')
		flagdata(vis=msname,antenna=antenna,spw=unflag_spws,mode='unflag')
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

def get_num_flag_baselines(msname,flagfrac=0.5):
	'''
	Function to get the antennas for which a certain amount of data are flagged
	Parameters:
	msname = Name of the measurement set
	flagfrac = Fraction of data flagged to consider the antennas as bad (default : 0.5)
	Return:
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

def flagger(msname,rms):
	'''
	Function to flag real and imaginary part of visibility based on rms threshold
	Parameters:
	msname = Name of the measurement set
	rmsthresh = RMS threshold for n-sigma flagging
	Return:
	Total flag points, new flag fraction
	'''
	mstool=ms()
	mstool.open(msname,nomodify=False)
	res_data=mstool.getdata('residual_data')['residual_data']
	flag=mstool.getdata('flag')
	flag_data=flag['flag']
	pos_flag=np.where(flag_data==True)
	num_flag=np.nansum(flag_data)
	res_data[pos_flag]=np.nan
	xx=res_data[0,:,:]
	xx_flag=flag_data[0,:,:]
	xy=res_data[1,:,:]
	xy_flag=flag_data[1,:,:]
	yx=res_data[2,:,:]
	yx_flag=flag_data[2,:,:]
	yy=res_data[3,:,:]
	yy_flag=flag_data[3,:,:]

	sigma_xxre=np.nanstd(np.real(xx))
	sigma_xxim=np.nanstd(np.imag(xx))
	sigma_xyre=np.nanstd(np.real(xy))
	sigma_xyim=np.nanstd(np.imag(xy))
	sigma_yxre=np.nanstd(np.real(yx))
	sigma_yxim=np.nanstd(np.imag(yx))
	sigma_yyre=np.nanstd(np.real(yy))
	sigma_yyim=np.nanstd(np.imag(yy))

	pos_xxre=np.where(np.real(xx)>rms*sigma_xxre)
	pos_xxim=np.where(np.imag(xx)>rms*sigma_xxim)
	pos_xyre=np.where(np.real(xy)>rms*sigma_xyre)
	pos_xyim=np.where(np.imag(xy)>rms*sigma_xyim)
	pos_yxre=np.where(np.real(yx)>rms*sigma_yxre)
	pos_yxim=np.where(np.imag(yx)>rms*sigma_yxim)
	pos_yyre=np.where(np.real(yy)>rms*sigma_yyre)
	pos_yyim=np.where(np.imag(yy)>rms*sigma_yyim)

	pos_xx=np.append(pos_xxre,pos_xxim)
	pos_xy=np.append(pos_xyre,pos_xyim)
	pos_yx=np.append(pos_yxre,pos_yxim)
	pos_yy=np.append(pos_yyre,pos_yyim)
	pos=np.unique(np.append(np.append(np.append(pos_xx,pos_xy),pos_yx),pos_yy))
	for i in pos:
			flag_data[:,:,i]=True
	flag['flag']=flag_data
	final_num_flag=np.nansum(flag_data)-num_flag
	mstool.putdata(flag)
	mstool.close()
	return final_num_flag


























