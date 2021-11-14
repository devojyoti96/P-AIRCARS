import os,glob,jprq,time,numpy as np

screen_path='/usr/bin'
os.path.join(screen_path)

def get_carta_url(basedir=''):
	if basedir!='':
		if basedir[-1]=='/':
			basedir=basedir[:-1]
		carta_public_url='#'
		if os.path.exists(basedir+'/cartalog'):
			with open(basedir+'/cartalog','r') as fil1:
				for line in fil1:
					if 'CARTA is accessible at ' in line:
						carta_url=line.split('CARTA is accessible at ')[-1].split('\n')[0]
			carta_port=carta_url.split('/?')[0].split(':')[-1]
			os.system('fuser -k -n tcp '+str(carta_port))
			del carta_port
			os.system('rm -rf '+basedir+'/cartalog')
		if os.path.exists(basedir+'/cartaerror'):
			os.system('rm -rf '+basedir+'/cartaerror')
		os.system('nohup ./CARTA-v2.0-redhat.AppImage --no_browser --top_level_folder '+basedir+' '+basedir+' > '+basedir+'/cartalog 2> '+basedir+'/cartaerror < /dev/null &')
		time.sleep(5)
		with open(basedir+'/cartalog','r') as fil1:
			for line in fil1:
				if 'CARTA is accessible at ' in line:
					carta_url=line.split('CARTA is accessible at ')[-1].split('\n')[0]
		carta_port=carta_url.split('/?')[0].split(':')[-1]
		carta_token='/?'+carta_url.split('/?')[-1]
		if os.path.exists('tcp.output'):
			os.system('rm -rf tcp.output')
		screen_cmd='jprq tcp '+carta_port
		screen_name=os.path.basename(basedir)+'_carta_screen'
		os.system('screen -S '+screen_name+' -X quit')
		os.system('screen -mdS '+screen_name)
		os.system('screen -S '+screen_name+' -X stuff "'+screen_cmd+'\n"')
		time.sleep(5)
		with open('tcp.output','r') as fil2:
			for line in fil2:
				if 'tcp.jprq.io' in line:
					carta_public_url=line.split(', ')[-1].split('\n')[0]+carta_token
		return 'http://'+carta_public_url
	else:
		return '#'

