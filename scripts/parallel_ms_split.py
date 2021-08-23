from paircars.access_ms import *
from paircars.basic_func import *
from paircars.flagger import *
from casatasks import *
from casatools import *
import multiprocessing as mp,os,time,psutil,datetime as dttime
from optparse import OptionParser
from CALIBRATE.access_calibrate import *

def split_func(params):
	'''
	Function to split an measurment set
	Parameters:
	params = List of the split parameters in the following order
				[vis,datacolumn,spw,timerange,chanwidth,timewidth,casacals,calibratecals], timewidth is in 't s' format
	Return:
	LIst of splited measurement sets
	'''
	casalog.showconsole(False)
	vis,datacolumn,chan,timestamp,chanwidth,timewidth,savedir,casacals,calibratecals=params
	outputvis=savedir+'/time_'+str(timestamp_to_mjdsec(timestamp.split('+')[0],format=0))+'_chan_'+str(chan)+'.ms'
	if os.path.exists(outputvis):
		os.system('rm -rf '+outputvis)
	split(vis=vis,outputvis=outputvis,spw=chan,timerange=timestamp,datacolumn=datacolumn,width=int(chanwidth),timebin=str(timewidth))
	outputvis=splited_ms_rename(outputvis,ref_time_chan=False,change_msname=True)
	applycal(vis=outputvis,gaintable=casacals,applymode='calflag',flagbackup=False)
	cal=CALIBRATE()
	for i in calibratecals:
		print ('cal.applycal(msname=\''+outputvis+'\',gaintable='+str(i)+',applymode=\'calflag\',flagbackup=False)\n')
		cal.applycal(msname=outputvis,gaintable=i,applymode='calflag',flagbackup=False)
	tb=table()
	tb.open(outputvis)
	cor_data=tb.getcol('CORRECTED_DATA')
	tb.close()
	tb.open(outputvis,nomodify=False)
	tb.putcol('DATA',cor_data)
	tb.flush()
	tb.close()
	os.system('rm -rf casa*log')
	return outputvis

def split_ms_parallel(paramList):
	'''
	Function to run multiple split jobs in parallel from same or different measurement sets
	Parameters:
	paramList = Lists of list of parameters for each split, format [[vis1,datacolumn1,spw1,timerange1,chanwidth1,timewidth1,casacals1,calibratecals1],\
	[vis2,datacolumn2,spw2,timerange2,chanwidth2,timewidth2,casacals2,calibratecals2],....]
	Return:
	List of splited ms
	'''
	pool=mp.Pool(processes=int(mp.cpu_count()*0.5))	
	results = pool.map(split_func,paramList)
	os.system('rm -rf casa*log')
	return results

# Function to run the script stand alone from command line
if __name__=='__main__':
	start_time=time.time()
	cwd=os.getcwd()
	usage= ' Split measurement set in parallel\n'
	parser = OptionParser(usage=usage)
	parser.add_option('--msname',dest="msname",default=None,help="Name of measurement set to split from",metavar="Measurement Set")
	parser.add_option('--savedir',dest='savedir',default=None,help='Name of directory to save splited measurement sets',metavar="Directory path")
	parser.add_option('--freq_interval',dest='skip_freq',default='MS frequency resolution',help='Spliting frequency interval in kHz',metavar="Float")
	parser.add_option('--time_interval',dest='skip_time',default='MS time resolution',help='Spliting temporal interval in second',metavar="Float")
	parser.add_option('--freq_list',dest='freq_list',default=None,help='Spliting frequency list in MHz (If given frequency interval will not be used)',metavar="Comma separated string")
	parser.add_option('--time_list',dest='time_list',default=None,help='Spliting timestamp list in \'yyyy/mm/dd/hh:mm:ss\' format (If given time interval will not be used)',\
						metavar="Comma separated string")
	parser.add_option('--freq_width',dest='freq_width',default='MS frequency resolution',help='Bandwidth of each splited ms in kHz start from frequency stamp',metavar="Float")
	parser.add_option('--time_width',dest='time_width',default='MS time resolution',help='Temporal width of each splited ms in second start from time stamp',metavar="Float")
	parser.add_option('--freq_avg',dest='freq_avg',default='MS frequency resolution',help='Frequency averaging in kHz',metavar="Float")
	parser.add_option('--time_avg',dest='time_avg',default='MS time resolution',help='Temporal averaging in second',metavar="Float")
	parser.add_option('--datacolumn',dest='datacolumn',default='data',help='Datacolumn to split',metavar="String")
	parser.add_option('--casa_caltables',dest='casacals',default='',help='CASA caltables',metavar="Comma separated string")
	parser.add_option('--calibrate_caltables',dest='calibratecals',default='',help='CALIBRATE caltables',metavar="Comma separated string")
	(options, args) = parser.parse_args()

	if str(options.msname)[-1]=='/':
		msname=str(options.msname)[:-1]
	else:
		msname=str(options.msname)
	if os.path.isdir(str(options.msname))==False or options.msname==None:
		print ('Measurement set is not present.\n')
		os.system('touch '+cwd+'/.Finished_spliting_'+os.path.basename(str(msname))+'_noms')
		os._exit(1)
	try:
		if options.savedir==None:
			savedir=os.path.dirname(os.path.abspath(str(options.msname)))+'/splted_ms'
		elif os.path.isdir(str(options.savedir)+'/splited_ms')==False:
			os.makedirs(str(options.savedir)+'/splited_ms')
			savedir=str(options.savedir)+'/splited_ms'
		else:
			savedir=str(options.savedir)+'/splited_ms'

		AM=AccessMS(str(options.msname))
		freqres=AM.calc_freqres()
		timeres=AM.calc_timeres()
		param_list=[]
		split_spws=[]
		split_timestamps=[]
		if str(options.freq_avg)=='MS frequency resolution' and str(options.time_avg)=='MS time resolution':
			chanwidth=1
			timebin=str(timeres)+'s'
		elif str(options.freq_avg)=='MS frequency resolution' and str(options.time_avg)!='MS time resolution':
			chanwidth=1
			timebin=str(options.time_avg)+'s'
		elif str(options.freq_avg)!='MS frequency resolution' and str(options.time_avg)=='MS time resolution':
			chanwidth=int(float(options.freq_avg)/freqres)
			timebin=str(timeres)+'s'
		elif str(options.freq_avg)!='MS frequency resolution' and str(options.time_avg)!='MS time resolution':
			chanwidth=int(float(options.freq_avg)/freqres)
			timebin=str(options.time_avg)+'s'

		if str(options.freq_width)=='MS frequency resolution' and str(options.time_width)=='MS time resolution':
			if str(options.skip_freq)!='MS frequency interval':
				freq_width=float(options.skip_freq)/10**3
			else:
				freq_width=float(freqres/10**3)
			if str(options.skip_time)!='MS frequency interval':
				time_width=float(options.skip_time)
			else:
				time_width=float(0)
		elif str(options.freq_width)=='MS frequency resolution' and str(options.time_width)!='MS time resolution':
			if str(options.skip_freq)!='MS frequency interval':
				freq_width=float(options.skip_freq)/10**3
			else:
				freq_width=float(freqres/10**3)
			time_width=float(options.time_width)
		elif str(options.freq_width)!='MS frequency resolution' and str(options.time_width)=='MS time resolution':
			freq_width=float(options.freq_width)/10**3
			if str(options.skip_time)!='MS frequency interval':
				time_width=float(options.skip_time)
			else:
				time_width=float(0)
		elif str(options.freq_width)!='MS frequency resolution' and str(options.time_width)!='MS time resolution':
			freq_width=float(options.freq_width)/10**3
			time_width=float(options.time_width)
		
		if options.freq_list!=None and options.time_list!=None:
			freq_list=str(options.freq_list).split(',')
			time_list=str(options.time_list).split(',')
		elif options.freq_list!=None and options.time_list==None:
			freq_list=str(options.freq_list).split(',')
			timestamps=AM.get_timestamps()
			skip_timestamps=int(float(options.skip_time)/timeres)
			time_list=[]
			for i in range(0,len(timestamps),skip_timestamps):
				time_list.append(timestamps[i])
		elif options.freq_list==None and options.time_list!=None:
			freqs=AM.get_freqs()/10**6
			skip_chan=int(float(options.skip_freq)/freqres)
			freq_list=[]
			for i in range(0,len(freqs),skip_chan):
				freq_list.append(freqs[i])
			time_list=str(options.time_list).split(',')
		elif options.freq_list==None and options.time_list==None:
			freqs=AM.get_freqs()/10**6
			skip_chan=int(float(options.skip_freq)/freqres)
			freq_list=[]
			for i in range(0,len(freqs),skip_chan):
				freq_list.append(freqs[i])
			timestamps=AM.get_timestamps()
			skip_timestamps=int(float(options.skip_time)/timeres)
			time_list=[]
			for i in range(0,len(timestamps),skip_timestamps):
				time_list.append(timestamps[i])

		for i in freq_list:
			split_spws.append('0:'+str(i)+'~'+str(float(i)+float(freq_width))+'MHz')
		for j in time_list:
			if time_width!=0:
				split_timestamps.append(j+'+'+str(dttime.timedelta(seconds=time_width)))
			else:
				split_timestamps.append(j)
		
		if str(options.casacals)=='':
			casacaltables=[]
		else:
			casacaltables=str(options.casacals).split(',')

		if str(options.calibratecals)=='':
			calibratecaltables=[]
		else:
			calibratecaltables=str(options.calibratecals).split(',')

		for i in split_spws:
			for j in split_timestamps:
				param_list.append([str(options.msname),str(options.datacolumn),i,j,chanwidth,timebin,savedir,casacaltables,calibratecaltables])

		print ('#########################################\n')
		print ('Spliting ms : '+str(options.msname)+' into total '+str(len(freq_list)*len(time_list))+' time frequency chuncks ..............\n')
		outputvis=split_ms_parallel(param_list)	
		print ('Total spliting time : '+str(time.time()-start_time)+'s\n')
		os.system('rm -rf casa*log')
		os.system('touch '+cwd+'/.Finished_spliting_'+os.path.basename(str(msname))+'_success')
		print ('#########################################\n')
	except Exception as e:
		os.system('touch '+cwd+'/.Finished_spliting_'+os.path.basename(str(msname))+'_'+str(e))
		os._exit(1)




