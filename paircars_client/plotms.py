import os,glob,subprocess


def plotms(msname='',xaxis='TIME',yaxis='amp',xdatacolumn='DATA',ydatacolumn='DATA',iteration='',ants='',baseline='',spw='',chan='',field=0,scan=1,corr=''):
	'''
	Function to MS data using shadems
	Parameters
	----------
	msname : str
		Measurement set name
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
	'''	
	shadems_args=['shadems','--xaxis='+xdatacolumn+':'+xaxis,'--yaxis='+ydatacolumn+':'+yaxis]
	if iteration!='':
		shadems_args.append('--iter-'+iteration)
	if ants!='':
		shadems_args.append('--ant-num='+ants)
	if baseline!='':
		shadems_args.append('--baseline='+baseline)
	if spw='':
		shadems_args.append('--spw='+spw)
	if chan!='':
		shadems_args.append('--chan='+chan)
	if field!=0:
		shadems_args.append('--field='+str(field))
	if scan!=1:
		shadems_args.append('--scan='+str(scan))
	











