from flask import Flask, render_template, request
import os,glob,jprq,numpy as np,psutil,socket,copy,paircars_client
from paircars_client.run_carta import *
from optparse import OptionParser
from paircars_client.plotms import *

requested_ms_copy=''
x_axis_copy=''
template_path=os.path.abspath(os.path.dirname(paircars_client.__file__))
app = Flask(__name__,template_folder=template_path+'/templates',static_folder=template_path+'/static')
def main(job_id=0,basedir='',template_path=''):
	carta_url=get_carta_url(basedir=basedir)
	logdir=basedir+'/Logs_and_Errors'
	@app.route('/')
	@app.route('/home')
	def home_page(job_id=job_id, methods=['GET', 'POST']):
		request_html=request.args.get('type')
		return render_template('home.html',job_id=job_id)

	@app.route('/logs', methods=['GET', 'POST'])
	def log_page(path=logdir,job_id=job_id):
		log_list=sorted(glob.glob(path+'/Logs/*'))
		log_list_nameonly=[os.path.basename(log).split('.log')[0] for log in log_list]
		error_list=sorted(glob.glob(path+'/Errors/*'))
		error_list_nameonly=[os.path.basename(error).split('.error')[0] for error in error_list]
		logerror='Logs and Errors'
		try:
			requested_file=request.args.get('type')
			if os.path.exists(requested_file):
				with open(requested_file, 'r') as f: 
					text=f.read()
				logerror_type=os.path.basename(requested_file).split('.')[-1]
				if logerror_type=='log':
					logerror='Log'
				elif logerror_type=='error':
					logerror='Error'
			else:
				text='No log files'
		except:
			text=''
		return render_template('log.html',job_id=job_id,logs=log_list_nameonly,logs_paths=log_list,errors=error_list_nameonly,errors_paths=error_list,text=text,logerror=logerror)

	@app.route('/carta')
	def open_carta(carta_url=carta_url):	
		return render_template('carta.html',carta_public_url=carta_url)

	@app.route('/plotms', methods=['GET', 'POST'])
	def open_plotms(path=basedir,job_id=job_id,template_path=template_path):
		message=0
		images=[]
		image_list=[]
		image=''
		requested_ms,x_axis,y_axis,xcolumn,ycolumn,iteration,plotflag,ants,baseline,spw,chan,field,scan,corr=np.load(path+'/plotms_input.npy',allow_pickle=True)
		all_dirs=[x[0] for x in os.walk(path)]
		x_list=['Time','Channel','Frequency','U','V','W','UVwave','Amplitude','Phase','Real','Imaginary']
		y_list=['Time','Channel','Frequency','U','V','W','UVwave','Amplitude','Phase','Real','Imaginary']	
		iteration_list = [None,'field','scan','spw','corr','ant','baseline']
		xdata_column=['DATA','CORRECTED_DATA','MODEL_DATA','DATA-MODEL_DATA','CORRECTED_DATA-MODEL_DATA']
		ydata_column=['DATA','CORRECTED_DATA','MODEL_DATA','DATA-MODEL_DATA','CORRECTED_DATA-MODEL_DATA']
		plotflag_list=[True,False]
		ms_list=[]
		ms_fullpath_list=[]
		noms=True
		if os.path.isdir(template_path+'/static/temp')==False:
			os.makedirs(template_path+'/static/temp')
		else:
			os.system('rm -rf '+template_path+'/static/temp/*')
		for i in all_dirs:
			if i.split('.')[-1]=='ms':
				ms_list.append(i.replace(path,''))
				ms_fullpath_list.append(i)
		if request.args.get('type')!=None:
			requested_ms=str(request.args.get('type'))
			noms=False
			np.save(path+'/plotms_input',[requested_ms,x_axis,y_axis,xcolumn,ycolumn,iteration,plotflag,ants,baseline,spw,chan,field,scan,corr])
		if request.args.get('xtype')!=None:
			x_axis=str(request.args.get('xtype'))
			np.save(path+'/plotms_input',[requested_ms,x_axis,y_axis,xcolumn,ycolumn,iteration,plotflag,ants,baseline,spw,chan,field,scan,corr])
		if request.args.get('ytype')!=None:
			y_axis=str(request.args.get('ytype'))
			np.save(path+'/plotms_input',[requested_ms,x_axis,y_axis,xcolumn,ycolumn,iteration,plotflag,ants,baseline,spw,chan,field,scan,corr])
		if request.args.get('xcolumntype')!=None:
			xcolumn=str(request.args.get('xcolumntype'))
			np.save(path+'/plotms_input',[requested_ms,x_axis,y_axis,xcolumn,ycolumn,iteration,plotflag,ants,baseline,spw,chan,field,scan,corr])
		if request.args.get('ycolumntype')!=None:
			ycolumn=str(request.args.get('ycolumntype'))
			np.save(path+'/plotms_input',[requested_ms,x_axis,y_axis,xcolumn,ycolumn,iteration,plotflag,ants,baseline,spw,chan,field,scan,corr])
		if request.args.get('flagtype')!=None:
			plotflag=str(request.args.get('flagtype'))
			np.save(path+'/plotms_input',[requested_ms,x_axis,y_axis,xcolumn,ycolumn,iteration,plotflag,ants,baseline,spw,chan,field,scan,corr])
		iteration=str(request.args.get('itertype'))
		np.save(path+'/plotms_input',[requested_ms,x_axis,y_axis,xcolumn,ycolumn,iteration,plotflag,ants,baseline,spw,chan,field,scan,corr])
		if requested_ms!=None:
			ms=os.path.basename(requested_ms)
		else:
			ms=None
		if request.method == 'POST':
			if len(images)!=0:
				for i in images:
					os.system('rm -rf '+i)
			if len(image_list)!=0:
				for i in image_list:
					if os.path.islink(i):
						os.unlink(i)
					else:
						os.system('rm -rf '+i)
			images=[]
			image_list=[]
			image=''
			os.system('rm -rf '+path+'/plot-*')
			if requested_ms!=None:
				form_data = request.form
				ants=form_data['antennas']
				if ants=='ant1,ant2....':
					ants=''
				baseline=form_data['baseline']
				if baseline=='ant1-ant2,ant2-ant3,....':
					baseline=''
				spw=form_data['spw']
				if spw=='spw0,spw1,...':
					spw=''
				chan=form_data['chans']
				if chan=='start:stop:step':
					chan=''
				field=form_data['field']
				scan=form_data['scan']
				corr=form_data['corr']
				if corr=='XX, XY, YX, YY, I, Q, U, V':
					corr=''
				np.save(path+'/plotms_input',[requested_ms,x_axis,y_axis,xcolumn,ycolumn,iteration,plotflag,ants,baseline,spw,chan,field,scan,corr])
				if x_axis=='Time':
					x_axis='TIME'
				elif x_axis=='Frequency':
					x_axis='FREQ'
				elif x_axis=='Channel':
					x_axis='CHAN'
				elif x_axis=='UVwave':
					x_axis='UV'
				elif x_axis=='Amplitude':
					x_axis='amp'
				elif x_axis=='Phase':
					x_axis='phase'
				elif x_axis=='Real':
					x_axis='real'
				elif x_axis=='Imaginary':
					x_axis='imag'
				if y_axis=='Time':
					y_axis='TIME'
				elif y_axis=='Frequency':
					y_axis='FREQ'
				elif y_axis=='Channel':
					y_axis='CHAN'
				elif y_axis=='UVwave':
					y_axis='UV'
				elif y_axis=='Amplitude':
					y_axis='amp'
				elif y_axis=='Phase':
					y_axis='phase'
				elif y_axis=='Real':
					y_axis='real'
				elif y_axis=='Imaginary':
					y_axis='imag'
				msname=ms_fullpath_list[ms_list.index(requested_ms)]
				if os.path.isdir(template_path+'/static/temp')==False:
					os.makedirs(template_path+'/static/temp')
				else:
					os.system('rm -rf '+template_path+'/static/temp/*')
				image_output,message=plotms(msname=msname,savedir=path,plotflag=plotflag,savelog=False,xaxis=x_axis,yaxis=y_axis,xdatacolumn=xcolumn,ydatacolumn=ycolumn,\
							iteration=iteration,ants=ants,baseline=baseline,spw=spw,chan=chan,field=field,scan=scan,corr=corr)
				for i in image_output:
					os.system('ln -s '+i+' '+template_path+'/static/temp/'+os.path.basename(i))
				noms=False
				image_list=['temp/'+os.path.basename(i) for i in image_output]
				if len(image_list)!=0:
					image=image_list[0]
			else:
				noms=True
		print (image_list)
		return render_template('plotms.html',ms_list=sorted(ms_list),xtypes=x_list,ytypes=y_list,ms=ms,xaxis=x_axis,yaxis=y_axis,job_id=job_id,xdata_column=xdata_column,\
					ydata_column=ydata_column,xcolumn=xcolumn,ycolumn=ycolumn,iteration_list=iteration_list,iteration=iteration,plotflag=plotflag,plotflag_list=plotflag_list,\
					ants=ants,baseline=baseline,spw=spw,chan=chan,field=field,scan=scan,corr=corr,noms=noms,image_list=image_list,message=message)

def getfreeport():
	port = np.random.randint(49152,65535)
	portsinuse=[]
	while True:
	    conns = psutil.net_connections()
	    for conn in conns:
	        portsinuse.append(conn.laddr[1])
	    if port in portsinuse:
	        port = np.random.randint(49152,65535)
	    else:
	        break
	return port

if __name__ == "__main__":
	template_path=os.path.abspath(os.path.dirname(paircars_client.__file__))
	usage= 'Start P-AIRCARS server for remote access'
	parser = OptionParser(usage=usage)
	parser.add_option('--basedir',dest="basedir",default=None,help="Name of P-AIRCARS base directory",metavar="Directory path")
	parser.add_option('--job_id',dest="job_id",default=0,help="P-AIRCARS Job ID",metavar="Integer")
	(options, args) = parser.parse_args()
	os.environ['FLASK_ENV']='development'
	main(job_id=int(options.job_id),basedir=str(options.basedir),template_path=template_path)
	port=getfreeport()
	hostname='-'.join(socket.gethostname().split('.'))
	public_paircars_url=hostname+'-'+str(int(options.job_id))+'.paircars.jprq.io'
	os.chdir(str(options.basedir))
	screen_cmd='jprq http '+str(port)+' -s='+str(hostname)+'-'+str(int(options.job_id))+'.paircars'
	screen_name=public_paircars_url+'_paircars_server_screen'
	if os.path.exists(str(options.basedir)+'/http.output'):
		os.system('rm -rf '+str(options.basedir)+'/http.output')
	os.system('screen -S '+screen_name+' -X quit')
	os.system('screen -mdS '+screen_name)
	os.system('screen -S '+screen_name+' -X stuff "'+screen_cmd+'\n"')
	os.chdir(os.getcwd())
	time.sleep(5)	
	with open(str(options.basedir)+'/http.output','r') as fil2:
		for line in fil2:
			if 'https://' in line:
				paircars_public_url=line
				break
	print ('Access P-AIRCARS remotely at : '+paircars_public_url+'\n')
	if os.path.exists(str(options.basedir)+'/plotms_input.npy'):
		os.system('rm -rf '+str(options.basedir)+'/plotms_input.npy')
	np.save(str(options.basedir)+'/plotms_input',[None,'Time','Amplitude','DATA','DATA',None,False,'','','','',0,1,''])
	app.run(host="localhost",use_reloader=False,port=port)
	
