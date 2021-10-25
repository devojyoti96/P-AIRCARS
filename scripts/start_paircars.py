import os,time,sys
os.chdir(os.getcwd())
sys.path.append(os.getcwd())
from selfcal_inputs import *
import selfcal_inputs as inputs
from paircars.access_ms import *
from paircars.basic_func import *

from optparse import OptionParser
if __name__=='__main__':
	usage= ' P-AIRCARS database manager'
	parser = OptionParser(usage=usage)
	parser.add_option('--fresh',dest="fresh",default=False,help="Want to start fresh calibration with P-AIRCARS",metavar="Boolean")
	parser.add_option('--restart',dest="restart",default=False,help="Want to re-start calibration with P-AIRCARS",metavar="Boolean")
	(options, args) = parser.parse_args()
	
print ('############################\n Starting P-AIRCARS..............\n############################\n')

# Validating measurement set dir 
################################
if os.path.isdir(msdir)==False:
	print('Measurement set directory does not exist. Check the measurement set path and re run. Exiting P-AIRCARS....\n')
	os.system('rm -rf casa*log')
	os._exit(1)

# Organising measurement sets
#############################
if basedir[-1]=='/':
	basedir=basedir[:-1]
print('Searching for measurment sets......\n')
file_list=sorted(glob.glob(msdir+'/*.ms'))
measurement_set_list=[]
basedir_list=[]
remove_count=0
prebasedir=''
for f in file_list:
	if os.path.isdir(f)==True:
		try:	
			msname=splited_ms_rename(f,ref_time_chan=False,change_msname=False)
			datestamp='_'.join(msname.split('time_')[-1].split('_')[:3])
			if eval(str(options.fresh))==True:
				if os.path.exists(basedir+'/basedir_for_'+datestamp)==True and basedir+'/basedir_for_'+datestamp!=prebasedir:
					print ('Removing existing base directory : '+basedir+'/basedir_for_'+datestamp)
					os.system('rm -rf '+basedir+'/basedir_for_'+datestamp)		
				prebasedir=basedir+'/basedir_for_'+datestamp
			if eval(str(options.restart))==True:
				os.system('rm -rf '+basedir+'/basedir_for_'+datestamp+'/.*paircars* ')	
			if os.path.isdir(basedir+'/basedir_for_'+datestamp)==False:
				os.makedirs(basedir+'/basedir_for_'+datestamp+'/data/')
			if basedir+'/basedir_for_'+datestamp not in basedir_list:
				basedir_list.append(basedir+'/basedir_for_'+datestamp)
			measurement_set_list.append(basedir+'/basedir_for_'+datestamp+'/data/'+os.path.basename(msname))
			if os.path.islink(basedir+'/basedir_for_'+datestamp+'/data/'+os.path.basename(msname))==False:
				if os.path.isdir(basedir+'/basedir_for_'+datestamp+'/data/')==False:
					os.makedirs(basedir+'/basedir_for_'+datestamp+'/data/')
				print('Linking '+f+' to '+basedir+'/basedir_for_'+datestamp+'/data/'+os.path.basename(msname)+'\n')
				os.system('ln -s '+f+' '+basedir+'/basedir_for_'+datestamp+'/data/'+os.path.basename(msname))
			elif os.path.islink(basedir+'/basedir_for_'+datestamp+'/data/'+os.path.basename(msname))==True:
				if os.path.isdir(os.path.realpath(basedir+'/basedir_for_'+datestamp+'/data/'+os.path.basename(msname)))==False:
					os.system('unlink '+basedir+'/basedir_for_'+datestamp+'/data/'+os.path.basename(msname))
					measurement_set_list.remove(basedir+'/basedir_for_'+datestamp+'/data/'+os.path.basename(msname))
		except Exception as e: 
			print ('Error occured : '+str(e)+'\n')
if len(measurement_set_list)==0:
	print('No valid measurement set is present. Put the correct data. Exiting P-AIRCARS.....\n')
else:
	print(str(len(measurement_set_list))+' measurement set has been found.\n')

finished_file_list=[]
for base_dir in basedir_list:
	os.system('cp -r selfcal_inputs.py '+base_dir+'/selfcal_inputs.py')
	inpfil=open(base_dir+'/selfcal_inputs.py','r+')
	lines=inpfil.readlines()
	for i in range(len(lines)):
		if 'basedir' in lines[i]:
			lines[i]='basedir\t\t\t\t=\t\''+base_dir+'\'\n'
	inpfil.seek(0)
	inpfil.writelines(lines)
	inpfil.close()
	os.system('rm -rf casa*log')
	screen_name='Control_'+base_dir.split('/')[-1]+'_'+str(inputs.job_id)
	cmd='control_paircars --basedir '+base_dir
	if os.path.isdir(base_dir)==False:
		os.makedirs(base_dir)
	finished_file='/'.join(base_dir.split('/')[:-1])+'/.Finished_'+base_dir.split('/')[-1]
	datestamp='/'.join(base_dir.split('/')[-1].split('_')[2:])
	finished_file_list.append(finished_file)
	print ('Spawning jobs for : '+datestamp+'\n')
	batch_file=paircars_instance_runner(cmd,base_dir,inputs.paircars_dir,screen_name,finished_file,inputs.job_id,prefix_cmds=[])
	os.system('sh '+batch_file)
os.system('touch '+inputs.basedir+'/.'+str(inputs.job_id)+'_paircars_running')







