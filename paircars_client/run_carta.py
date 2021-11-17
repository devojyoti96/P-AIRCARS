import os,glob,jprq,time,numpy as np
import paircars_client

datadir=os.path.abspath(os.path.dirname(paircars_client.__file__))
try:
	screen_path=str(np.load(datadir+'/screen_path.npy',allow_pickle=True))
except:
	screen_path='/usr/bin'
os.path.join(screen_path)
def get_carta_url(basedir=''):
	os.environ['APPIMAGE_EXTRACT_AND_RUN']='1'
	if basedir!='':
		if basedir[-1]=='/':
			basedir=basedir[:-1]
		carta_public_url='#'
		if os.path.exists(basedir+'/cartalog'):
			with open(basedir+'/cartalog','r') as fil1:
				for line in fil1:
					if 'CARTA is accessible at ' in line:
						carta_url=line.split('CARTA is accessible at ')[-1].split('\n')[0]
			try:
				carta_port=carta_url.split('/?')[0].split(':')[-1]
				os.system('fuser -k -n tcp '+str(carta_port))
				del carta_port
			except:
				pass
			os.system('rm -rf '+basedir+'/cartalog')
		if os.path.exists(basedir+'/cartaerror'):
			os.system('rm -rf '+basedir+'/cartaerror')
		os.system('nohup '+datadir+'/CARTA.AppImage --no_browser --top_level_folder '+basedir+' '+basedir+' > '+basedir+'/cartalog 2> '+basedir+'/cartaerror < /dev/null &')
		time.sleep(5)
		with open(basedir+'/cartalog','r') as fil1:
			for line in fil1:
				if 'CARTA is accessible at ' in line:
					carta_url=line.split('CARTA is accessible at ')[-1].split('\n')[0]
		carta_port=carta_url.split('/?')[0].split(':')[-1]
		carta_token='/?'+carta_url.split('/?')[-1]
		if os.path.exists(basedir+'/tcp.output'):
			os.system('rm -rf '+basedir+'/tcp.output')
		cwd=os.getcwd()
		os.chdir(basedir)
		screen_cmd='jprq tcp '+carta_port
		screen_name=os.path.basename(basedir)+'_carta_screen'
		os.system('screen -S '+screen_name+' -X quit')
		os.system('screen -mdS '+screen_name)
		os.system('screen -S '+screen_name+' -X stuff "'+screen_cmd+'\n"')
		os.chdir(cwd)
		time.sleep(5)
		with open(basedir+'/tcp.output','r') as fil2:
			for line in fil2:
				if 'tcp.jprq.io' in line:
					carta_public_url=line.split(', ')[-1].split('\n')[0]+carta_token
		return 'http://'+carta_public_url
	else:
		return '#'

