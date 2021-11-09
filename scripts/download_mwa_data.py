import os,time
from optparse import OptionParser
from paircars.basic_func import get_OBSID_for_timerange

if __name__=='__main__':
	start_time=time.time()
	usage= 'Submit and download MWA data'
	parser = OptionParser(usage=usage)
	parser.add_option('--API_key',dest="apikey",default=None,help="MWA ASVO API key",metavar="String")
	parser.add_option('--dest_dir',dest="destdir",default=os.getcwd(),help="Directory to download",metavar="Directory path")
	parser.add_option('--obsids',dest='obsids',default=None,help='Observation IDs',metavar='String, comma separated')
	parser.add_option('--cal_obsids',dest='cal_obsids',default=None,help='Observation IDs for calibration data',metavar="String, comma separated")
	parser.add_option('--timerange',dest='timerange',default=None,\
		help='Time range to download MWA data (format : yy0/mm0/dd0/hh0:mm0:ss0.ff0~yy1/mm1/dd1/hh1:mm1:ss1.ff1,yy2/mm2/dd2/hh2:mm2:ss2.ff2~....)',metavar="String, comma separated")
	parser.add_option('--cal_download',dest='caldownload',default=True,help='Download calibration data or not',metavar="Boolean")
	parser.add_option('--project_ID',dest='project_ID',default='G0002',help='MWA project ID',metavar="String")
	(options, args) = parser.parse_args()

if options.apikey==None:
	print ('Please provide your MWA ASVO API key.\n')
	os._exit(0)
os.environ['MWA_ASVO_API_KEY']=str(options.apikey)

if options.obsids==None:
	obsids=[]
else:
	obsids=str(options.obsids).split(',')

if options.cal_obsids==None:
	cal_obsids=[]
else:
	cal_obsids=str(options.cal_obsids).split(',')
if options.timerange!=None:
	obsid_list,cal_obsid_list=get_OBSID_for_timerange(timerange=str(options.timerange),caldata=eval(str(options.caldownload)),project_ID=str(options.project_ID))
	if len(obsid_list)!=0:
		obsids+=obsid_list
	if len(cal_obsid_list)!=0:
		cal_obsids+=cal_obsid_list
del obsid_list
del cal_obsid_list

if len(obsids)==0 and len(cal_obsids)==0:
	print ('No observation to download.\n')
	os._exit(0)

cwd=os.getcwd()
if os.path.isdir(str(options.destdir))==False:
	try:
		os.makedirs(str(options.destdir))
	except:
		print ('Please provide correct download directory.\n')
		os._exit(0)
os.chdir(str(options.destdir))
if os.path.isfile('mwa_download.csv'):
	os.system('rm -rf mwa_download.csv')
fil=open('mwa_download.csv','a+')

if len(obsids)!=0:
	obsids=[int(i) for i in obsids]
	for obsid in obsids:
		fil.write('obs_id='+str(obsid)+',job_type=c,timeres=0.5,freqres=40,edgewidth=80,conversion=ms,allowmissing=true,flagdcchannels=true,norfi=true,'+\
					'noprecomputedflags=true,noflagmissings=true,noantennapruning=true\n')

if len(cal_obsids)!=0:
	cal_obsids=[int(i) for i in cal_obsids]
	for obsid in cal_obsids:
		fil.write('obs_id='+str(obsid)+',job_type=c,timeres=0.5,freqres=40,edgewidth=80,conversion=ms,allowmissing=true,flagdcchannels=true,norfi=true,'+\
					'noprecomputedflags=true,noflagmissings=true,noantennapruning=true\n')

fil.seek(0)
fil.close()

print ('mwa_client -c '+str(options.destdir)+'/mwa_download.csv -d '+str(options.destdir)+'\n')
os.system('mwa_client -c '+str(options.destdir)+'/mwa_download.csv -d '+str(options.destdir))

if len(cal_obsids)!=0:
	if os.path.isdir(destdir+'/caldata')==False:
		os.makedirs(destdir+'/caldata')
	for obsid in cal_obsids:
		os.system('mv '+str(obsid)+'* '+destdir+'/caldata')

print ('\n##########\nData download completed.\n##########\n')
end_time=time.time()
run_time=time.strftime('%Hh %Mm %Ss',time.gmtime(end_time-start_time))
os.chdir(cwd)
print ('Total time : '+str(run_time)+'\n##########\n')





