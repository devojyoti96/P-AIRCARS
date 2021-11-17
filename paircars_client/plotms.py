import os,glob,subprocess,numpy as np
from casatools import msmetadata

def plotms(msname='',savedir='',plotflag=False,savelog=True,xaxis='TIME',yaxis='amp',xdatacolumn='DATA',ydatacolumn='DATA',\
			iteration='',ants='',baseline='',spw='',chan='',field=0,scan=1,corr=''):
	'''
	Function to MS data using shadems
	NOte : error checks for wrong input parameters has not been implemented, so please give correct inputs.
	Parameters
	----------
	msname : str
		Measurement set name
	savedir : str
		Directory to save final png files (default : current directory)
	plotflag : bool
		Plot flagged data or not (default : False)
	savelog : bool
		Save log file or not (default : False)
	xaxis : str
		X axis (TIME, CHAN, FREQ, CORR, U, V, W, UV, 'amp', 'phase', 'real', 'imag')
	yaxis : str
		Y axis (TIME, CHAN, FREQ, CORR, U, V, W, UV, 'amp', 'phase', 'real', 'imag')
	xdatacolumn : str
		X data column ('DATA','CORRECTED_DATA','MODEL_DATA','DATA-MODEL_DATA','CORRECTED_DATA-MODEL_DATA')	
	iteration : str
		Iteration axis ('field','scan','spw','corr','ant','baseline')
	ants : str
		Antenna numbers; comma separated
	baseline : str
		Baselines to plot, 'ant1-ant2', comma separated
	spw : str
		Spectral windows (comma separated)
	chan : str
		Channels start:stop:step format
	field : int
		Field ID 
	scan : int
		Scan number
	corr : str
		Correction or Stokes, 'XX,XY,YX,YY,RR,RL,LR,LL,I,Q,U,V'	
	Returns
	-------
	list
		Output image name list
	int
		Success message (0 or 1)
	'''	
	cwd=os.getcwd()
	if msname=='':
		print ('No measurement set is given.\n')
		return [],1
	else:
		os.chdir(msname) 
	if xaxis=='amp' or xaxis=='phase' or  xaxis=='real' or xaxis=='imag':
		shadems_args=['shadems','--xaxis='+xdatacolumn+':'+xaxis]
	else:
		shadems_args=['shadems','--xaxis='+xaxis]
	if yaxis=='amp' or yaxis=='phase' or  yaxis=='real' or yaxis=='imag':
		shadems_args.append('--yaxis='+ydatacolumn+':'+yaxis)
	else:
		shadems_args.append('--yaxis='+yaxis)
	if iteration!='' and iteration!=None:
		if iteration!='None':
			shadems_args.append('--iter-'+iteration)
	if ants!='':
		shadems_args.append('--ant-num='+ants)
	if baseline!='':
		shadems_args.append('--baseline='+baseline)
	if spw!='':
		shadems_args.append('--spw='+spw)
	if chan!='':
		shadems_args.append('--chan='+chan)
	if field!=0:
		shadems_args.append('--field='+str(field))
	if scan!=1:
		shadems_args.append('--scan='+str(scan))
	if corr!='':
		shadems_args.append('--corr='+corr)
	md=msmetadata()
	md.open(msname)
	nrows=int(md.nrows())
	nchan=md.nchan(0)
	chunk_size=int(nrows/3)
	shadems_args.append('-j0')
	shadems_args.append('--cmap=kb')
	shadems_args.append('--no-lim-save')
	shadems_args.append('-z'+str(chunk_size))
	if xaxis=='CHAN' or xaxis=='FREQ' or yaxis=='CHAN' or yaxis=='FREQ':
		if nchan<=1 and xaxis=='CHAN':
			shadems_args.append('--xmin=-50')
			shadems_args.append('--xmax=50')
		elif nchan<=1 and xaxis=='CHAN':
			shadems_args.append('--ymin=-50')
			shadems_args.append('--ymax=50')
		elif nchan<=1 and xaxis=='FREQ':
			freq=md.chanfreqs(0)[0]
			freqres=md.chanres(0)[0]
			shadems_args.append('--xmin='+str(freq-5*freqres))
			shadems_args.append('--xmax='+str(freq+5*freqres))
		elif nchan<=1 and yaxis=='FREQ':
			freq=md.chanfreqs(0)[0]
			freqres=md.chanres(0)[0]
			shadems_args.append('--ymin='+str(freq-5*freqres))
			shadems_args.append('--ymax='+str(freq+5*freqres))
		if xaxis=='CHAN' and nchan>1:
			shadems_args.append('--colour-by=CHAN')
			shadems_args.append('--cnum='+str(nchan))
		elif xaxis=='FREQ' and nchan>1:
			minfreq=np.min(md.chanfreqs(0))
			maxfreq=np.max(md.chanfreqs(0))
			freqres=md.chanres(0)[0]
			shadems_args.append('--xmin='+str(minfreq-5*freqres))
			shadems_args.append('--xmax='+str(maxfreq+5*freqres))
	md.close()
	if savedir!='':
		if os.path.isdir(savedir)==False:
			try:
				os.makedirs(savedir)
				savedir=savedir
			except:
				savedir=cwd
				pass
		else:
			savedir=savedir
	else:
		savedir=cwd
	shadems_args.append('--dir='+savedir)
	if plotflag==True:
		shadems_args.append('--noflags')
	shadems_args.append('--fontsize=30')
	shadems_args.append(msname)
	shadems_cmd=' '.join(shadems_args)
	if os.path.exists('log-shadems.txt'):
		os.system('rm -rf log-shadems.txt')
	a=os.system(shadems_cmd)
	if a!=0:
		return [],1
	outputname=[]
	with open('log-shadems.txt','r') as fil:
		for line in fil:
			if 'wrote ' in line:
				outputname.append(line.split('wrote ')[-1].split('\n')[0])
	if savelog==False:
		os.system('rm -rf '+savedir+'/log-shadems.txt '+savedir+'/casa*log')
	os.chdir(cwd)
	return outputname,0









