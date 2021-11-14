from flask import Flask, render_template, request
import os,glob,jprq,numpy as np,psutil,socket,copy
from run_carta import *
from optparse import OptionParser
requested_ms_copy=''
x_axis_copy=''
app = Flask(__name__)
def main(job_id=0,basedir=''):
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
	def open_plotms(path=basedir):
		requested_ms,x_axis,y_axis=np.load('plotms_input.npy',allow_pickle=True)
		all_dirs=[x[0] for x in os.walk(path)]
		x_list=['Time','Channel','Frequency','U','V','W','UVwave','Amplitude','Phase','Real','Imaginary']
		y_list=['Time','Channel','Frequency','U','V','W','UVwave','Amplitude','Phase','Real','Imaginary']
		ms_list=[]
		for i in all_dirs:
			if i.split('.')[-1]=='ms':
				ms_list.append(i.replace(path,''))
		if request.args.get('type')!=None:
			requested_ms=str(request.args.get('type'))
			np.save('plotms_input',[requested_ms,x_axis,y_axis])
		if request.args.get('xtype')!=None:
			x_axis=str(request.args.get('xtype'))
			np.save('plotms_input',[requested_ms,x_axis,y_axis])
		if request.args.get('ytype')!=None:
			x_axis=str(request.args.get('ytype'))
			np.save('plotms_input',[requested_ms,x_axis,y_axis])
		print ('Requested MS : '+requested_ms+'\n')
		print ('X axis : '+x_axis+'\n')
		print ('Y axis : '+y_axis+'\n')
		return render_template('plotms.html',ms_list=ms_list,xtypes=x_list,ytypes=y_list)

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
	usage= 'Start P-AIRCARS server for remote access'
	parser = OptionParser(usage=usage)
	parser.add_option('--basedir',dest="basedir",default=None,help="Name of P-AIRCARS base directory",metavar="Directory path")
	parser.add_option('--job_id',dest="job_id",default=0,help="P-AIRCARS Job ID",metavar="Integer")
	(options, args) = parser.parse_args()
	os.environ['FLASK_ENV']='development'
	main(job_id=int(options.job_id),basedir=str(options.basedir))
	port=getfreeport()
	hostname='-'.join(socket.gethostname().split('.'))
	public_paircars_url=hostname+'-'+str(int(options.job_id))+'.paircars.jprq.io'
	screen_cmd='jprq http '+str(port)+' -s='+str(hostname)+'-'+str(int(options.job_id))+'.paircars'
	screen_name=os.path.basename(str(options.basedir))+'_paircars_server_screen'
	os.chdir(os.getcwd())
	if os.path.exists('http.output'):
		os.system('rm -rf http.output')
	os.system('screen -S '+screen_name+' -X quit')
	os.system('screen -mdS '+screen_name)
	os.system('screen -S '+screen_name+' -X stuff "'+screen_cmd+'\n"')
	time.sleep(5)	
	with open('http.output','r') as fil2:
		for line in fil2:
			if 'https://' in line:
				paircars_public_url=line
				break
	print ('Access P-AIRCARS remotely at : '+paircars_public_url+'\n')
	np.save('plotms_input',['','Time','Amplitude'])
	app.run(host="localhost",use_reloader=False,port=port)
	
