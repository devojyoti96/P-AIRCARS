import os,sys
sys.path.append(os.getcwd())
from casatools import *
from casatasks import *
import selfcal_inputs as inputs
from selfcal_inputs import *
from datetime import datetime as dt, timedelta
import logging,numpy as np,sys,copy,glob,psutil,json,urllib.request
from paircars.basic_func import *
from paircars.access_ms import *
from astropy.io import fits
'''
Code is written by Devojyoti Kansabanik, 01 Feb, 2021

Code to validate the PAIRCARS input parameters
#######################
# Here we are validating the user given parameters. 
# If some parameters are found to be unsuitable, code will take default suitable value or show error and stop.
# If interactive mode is on, code will ask the user to give correct value
#######################
'''
if __name__=='__main__':
	print ('#######################################\n')
	print ('Starting PAIRCARS......................\n')
	print ('#######################################\n')

	print ('Validating inputs........\n')
	# Validating basedir path
	#########################
	if os.path.isfile('selfcal_inputs.py')==False:
		print ('Input file does not exist. Exiting PAIRCARS......\n')
		os.system('rm -rf casa*log')
		os._exit(1)
	else:
		inpfil=open('selfcal_inputs.py','r+')
		lines=inpfil.readlines()
	if basedir=='':
		print ('Path of the base directory is empty\n')
		if interactive==True:
			basedir=input('Give the base directory path:\n')
			if basedir=='':
				print ('Base directory path is still empty. Please become stable first and then run the code.\n')
				os.system('rm -rf casa*log')
				os._exit(1)
		else:
			print ('Base directory path is empty. Exititng PAIRCARS.....\n')
			os.system('rm -rf casa*log')
			os._exit(1)
		
	if os.path.isdir(basedir)==False: # If basedir is not present making it
		print ('Base directory is not present. Making base directory at :'+basedir)
		try:
			os.makedirs(basedir)
		except:
			print('Give correct path. Base directory could not be made')

	if basedir[-1]=='/':
		basedir=basedir[:-1]

	if os.path.isdir(basedir+'/data')==False:
		os.makedirs(basedir+'/data')

	# Logger initiating
	###################
	formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
	mainlog = logging.getLogger('paircars_main_log')
	mainlog.setLevel(logging.DEBUG)
	console=logging.StreamHandler(sys.stdout)
	console.setFormatter(formatter)
	mainlog.addHandler(console)
	mainlog.propagate = False

	# Starting PAIRCARS
	###################

	if os.path.isfile(basedir+'/.paircars_running'):
		mainlog.error('PAIRCARS is already running in this base directory. Choose a different directory. Exiting PAIRCARS......\n')
		os.system('rm -rf casa*log')
		os._exit(1)
	elif os.path.isfile(basedir+'/.paircars_finished'):
		mainlog.error('PAIRCARS is already have final results in this base directory.\n')
		want_to_continue=input('Do you want to run it again? Y/y/N/n')
		if want_to_continue=='Y' or want_to_continue=='y':
			os.system('mv '+basedir+'/PAIRCARS_mainlog.log temp.log')
			os.system('rm -rf '+basedir+'/*')
			os.system('rm -rf '+basedir+'.paircars_failed')
			loglist=len(glob.glob(basedir+'/PAIRCARS_mainlog*.log'))
			os.system('mv temp.log '+basedir+'/PAIRCARS_mainlog_'+str(loglist)+'.log')
			filehandle=logging.FileHandler(basedir+'/PAIRCARS_mainlog.log')
			filehandle.setFormatter(formatter)
			mainlog.addHandler(filehandle)
		else:
			mainlog.info('Exiting the code.\n')
			os.system('rm -rf casa*log')
			os._exit(1)
	elif os.path.isfile(basedir+'/.paircars_failed'):
		mainlog.error('PAIRCARS have failed in this base directory. Cleaning the base directory for fresh start up.\n')
		os.system('mv '+basedir+'/PAIRCARS_mainlog.log temp.log')
		file_list=glob.glob(basedir+'/*')
		loglist=len(glob.glob(basedir+'/PAIRCARS_mainlog*.log'))
		for f in file_list:
			if f!=msdir or f!=nasedir+'/data':
				os.system('rm -rf '+f)
		os.system('rm -rf '+basedir+'.paircars_failed')
		os.system('mv temp.log '+basedir+'/PAIRCARS_mainlog_'+str(loglist)+'.log')
		filehandle=logging.FileHandler(basedir+'/PAIRCARS_mainlog.log')
		filehandle.setFormatter(formatter)
		mainlog.addHandler(filehandle)
	else:
		filehandle=logging.FileHandler(basedir+'/PAIRCARS_mainlog.log')
		filehandle.setFormatter(formatter)
		mainlog.addHandler(filehandle)
		mainlog.info('Starting PAIRCARS in fresh base directory.\n')


	# Validating calibrator caltables
	#################################
	if len(calibrator_caltable)==0:
		mainlog.info('No calibrator caltable is provided\n')
	else:
		for i in range(len(calibrator_caltable)): # If wrong calibrator caltable path is given, make the calibrator caltable list empty
			if os.path.isdir(calibrator_caltable[i])==False:
				mainlog.info('Calibrator caltable '+str(calibrator_caltable[i])+' is not present\n')
				calibrator_caltable.remove(calibrator_caltable[i])
			else:
				mainlog.info('Calibrator caltable '+str(calibrator_caltable[i])+' is present\n')
		for i in range(len(lines)):
			if 'calibrator_caltable' in lines[i]:
				lines[i]='calibrator_caltable\t=\t'+str(calibrator_caltable)+'\n'

	
	# Validating Safety standard and quality factor
	###############################################
	if type(safety_factor)==int and safety_factor>2:
		mainlog.info('Safety factor is outside possible range. Setting it to default value : 1\n')
		safety_factor=1
		for i in range(len(lines)):
			if 'safety_factor' in lines[i]:
				lines[i]='safety_factor\t\t=\t'+str(safety_factor)+'\n'
	elif type(safety_factor)!=int:
		mainlog.info('Safety factor was not integer. Setting it to default value : 1\n')
		safety_factor=1
		for i in range(len(lines)):
			if 'safety_factor' in lines[i]:
				lines[i]='safety_factor\t\t=\t'+str(safety_factor)+'\n'
	###
	if type(quality_factor)==int and quality_factor>2:
		mainlog.info('Quality factor is outside possible range. Setting it to default value : 1\n')
		quality_factor=1
		for i in range(len(lines)):
			if 'quality_factor' in lines[i]:
				lines[i]='quality_factor\t\t=\t'+str(quality_factor)+'\n'
	elif type(quality_factor)!=int:
		mainlog.info('Quality factor was not integer. Setting it to default value : 1\n')
		quality_factor=1
		for i in range(len(lines)):
			if 'quality_factor' in lines[i]:
				lines[i]='quality_factor\t\t=\t'+str(quality_factor)+'\n'

	# Validating verbose and interactive and do_decor_correction
	############################################################
	if type(interactive)!=bool:
		mainlog.info('Interactive parameter was not boolean. Setting it to default value False.\n')
		interactive=False
		for i in range(len(lines)):
			if 'interactive' in lines[i]:
				lines[i]='interactive\t\t\t=\tFalse\n'
	if type(verbose)!=bool:
		mainlog.info('Verbose parameter was not boolean. Setting it to default value False.\n')
		verbose=False
		for i in range(len(lines)):
			if 'verbose' in lines[i]:
				lines[i]='verbose\t\t\t\t=\tFalse\n'
	if type(do_decor_correction)!=bool:
		mainlog.info('do_decor_correction parameter was not boolean. Setting it to default value True.\n')
		do_decor_correction=True
		for i in range(len(lines)):
			if 'do_decor_correction' in lines[i]:
				lines[i]='do_decor_correction\t=\tTrue\n'

	# Validating ref ant
	####################
	if type(ref_ant)==int and ref_ant>20:
		mainlog.info('Reference antenna was chosen beyond core 20 antennas. Setting it to default value 1.\n')
		ref_ant=1
		for i in range(len(lines)):
			if 'ref_ant' in lines[i]:
				lines[i]='ref_ant\t\t\t\t=\t'+str(ref_ant)+'\n'
	elif type(ref_ant)!=int:
		mainlog.info('Given reference antenna parameter is not an integer.\n')
		try:
			ref_ant=int(ref_ant)
		except:
			ref_ant=1
		mainlog.info('Change reference antenna to :'+str(ref_ant)+'\n')
		for i in range(len(lines)):
			if 'ref_ant' in lines[i]:
				lines[i]='ref_ant\t\t\t\t=\t'+str(ref_ant)+'\n'

	# Validating calc_image_parameters and calc_selfcalib_params
	############################################################
	if type(calc_image_parameters)!=bool:
		mainlog.info('calc_image_parameters was not boolean. Setting it to default value True.\n')
		calc_image_parameters=True
		for i in range(len(lines)):
			if 'calc_image_parameters' in lines[i]:
				lines[i]='calc_image_parameters\t=\tTrue\n'

	if type(calc_selfcalib_params)!=bool:
		mainlog.info('calc_selfcalib_params was not boolean. Setting it to default value True.\n')
		calc_selfcalib_params=True
		for i in range(len(lines)):
			if 'calc_selfcalib_params' in lines[i]:
				lines[i]='calc_selfcalib_params\t=\tTrue\n'

	# Validating the exsistence of maskfile
	#######################################
	if os.path.isdir(maskfile)==False:
		mainlog.info('Given maskfile is not present. Setting it to default value \' \'\n')
		maskfile=''
		for i in range(len(lines)):
			if 'maskfile' in lines[i]:
				lines[i]='maskfile\t\t\t\t=\t\''+maskfile+'\'\n'

	# Validating want_auto_masking
	##############################
	if type(want_auto_masking)!=bool:
		mainlog.info('want_auto_masking was not boolean. Seeting it to default value False.\n')
		want_auto_masking=False
		for i in range(len(lines)):
			if 'want_auto_masking' in lines[i]:
				lines[i]='want_auto_masking\t\t=\tFalse\n'


	# Validating cpu_frac
	#####################
	if cpu_frac>0.8:
		cpu_frac=0.8
		mainlog.info('CPU fraction was chosen more than 80%. Setting it to default value : 0.8.\n')
		for i in range(len(lines)):
			if 'cpu_frac' in lines[i]:
				lines[i]='cpu_frac\t\t\t\t=\t'+str(cpu_frac)+'\n'

	if cpu_frac==0:
		cpu_frac=0.5
		mainlog.info('CPU fraction was chosen to be 0%. Setting it to default value : 0.5.\n')
		for i in range(len(lines)):
			if 'cpu_frac' in lines[i]:
				lines[i]='cpu_frac\t\t\t\t=\t'+str(cpu_frac)+'\n'

	# Validating do_bandpass and do_polcal
	######################################
	if type(do_bandpass)!=bool and quality_factor==0:
		mainlog.info('do_bandpass was not boolean. Quality factor = '+str(quality_factor)+'. Setting it to default value False.\n')
		do_bandpass=False
	elif type(do_bandpass)!=bool:
		mainlog.info('do_bandpass was not boolean. Quality factor = '+str(quality_factor)+'. Setting it to default value True.\n')
		do_bandpass=True
	for i in range(len(lines)):
		if 'do_bandpass' in lines[i]:
			lines[i]='do_bandpass\t\t\t\t=\t'+str(do_bandpass)+'\n'

	if type(do_polcal)!=bool:
		mainlog.info('do_polcal was not boolean. Setting it to default value True.\n')
		do_polcal=True
	for i in range(len(lines)):
		if 'do_polcal' in lines[i]:
			lines[i]='do_polcal\t\t\t\t=\tTrue\n'

	# Validating save_true_loc_image
	################################
	if type(save_true_loc_image)!=bool:
		mainlog.info('save_true_loc_image was not boolean. Setting it to default value False.\n')
		save_true_loc_image=False
	for i in range(len(lines)):
		if 'save_true_loc_image' in lines[i]:
			lines[i]='save_true_loc_image\t\t=\tFalse\n'

	# Validating send_notification and email address
	################################################
	if type(send_notification)!=bool or email=='':
		mainlog.info('send_notification was not boolean or no valid email is given. Setting it to default value False.\n')
		send_notification=False
		for i in range(len(lines)):
			if 'send_notification' in lines[i]:
				lines[i]='send_notification\t\t=\t'+str(send_notification)+'\n'

	if email!='' and send_notification==False:
		mainlog.info('send_notification = False. Setting email address to to default value None.\n')	
		email=''
		for i in range(len(lines)):
			if 'email' in lines[i]:
				lines[i]='email\t\t\t\t\t=\t\'\'\n'

	# Validating basic imaging parameters
	######################################
	if calc_image_parameters==False:
		if float(cellsize)==0 or imsize[0]==0:
			mainlog.info('calc_image_parameters=False, Have you drunk? Pixel size or image size can not 0. Setting it to default value.\n')	
			cellsize='nan'
			imsize='nan'
			for i in range(len(lines)):
				if 'cellsize' in lines[i]:
					lines[i]='cellsize\t\t\t\t=\t'+str(cellsize)+'\n'
				if 'imsize' in lines[i]:
					lines[i]='imsize\t\t\t\t\t=\t['+str(imsize)+']\n'

		if len(multiscale_scales)!=0:
			types=[type(i) for i in multiscale_scales]
			c=0
			for i in types:
				if i!=int:
					c+=1
			if c!=0:
				mainlog.info('Multiscale scale list should have anything other than integer. Setting it to default.')
				multiscale_scales='nan'
				for i in range(len(lines)):
					if 'multiscale_scales' in lines[i]:
						lines[i]='multiscale_scales\t\t=\t'+str(multiscale_scales)+'\n'

		uvtaper_format_err=False
		if uvtaper!='':
			if 'lambda' not in uvtaper:
				uvtaper_format_err=True
			else:
				uvtaper_list=uvtaper.split('lambda')
				if uvtaper_list[-1]!='':
					uvtaper_format_err=True
				else:
					uvtaper_str=uvtaper_list[0]
					try:
						uvtaper=float(uvtaper_str)
					except:
						try:
							uvtaper=float(uvtaper_str[:-1])
						except:
							uvtaper_format_err=True
			if uvtaper_format_err==True:
				mainlog.info('uvtaper format is wrong. Setting it to default value.\n')
				for i in range(len(lines)):
					if 'uvtaper' in lines[i]:
						lines[i]='uvtaper\t\t\t\t\t=\t\'\'\n'

	# Validating selfcal parameters
	###############################
	if calc_selfcalib_params==False:
		#
		try:
			start_sigma=float(start_sigma)
			if start_sigma<8.0 and interactive==False:
				mainlog.info('Start sigma is very low at :'+str(start_sigma)+'. Setting it to default value 8.0.\n')
				start_sigma=8.0
				for i in range(len(lines)):
					if 'start_sigma' in lines[i]:
						lines[i]='start_sigma\t\t=\t'+str(8.0)+'\n'
		except:
			mainlog.info('Start sigma was not float or integer. Setting it to default value 10.0.\n')
			start_sigma=10.0
			for i in range(len(lines)):
				if 'start_sigma' in lines[i]:
					lines[i]='start_sigma\t\t=\t'+str(10.0)+'\n'
		#
		try:
			sigma_step=float(sigma_step)
			if sigma_step>2:
				mainlog.info('Sigma step was greater than 2. Setting it to default value 1.5.\n')
				sigma_step=1.5
				for i in range(len(lines)):
					if 'sigma_step' in lines[i]:
						lines[i]='sigma_step\t\t=\t'+str(1.5)+'\n'
			elif sigma_step==0:
				mainlog.info('Sigma step was 0.0. Setting it to default value 0.5.\n')
				sigma_step=0.5
				for i in range(len(lines)):
					if 'sigma_step' in lines[i]:
						lines[i]='sigma_step\t\t=\t'+str(0.5)+'\n'
		except:
			mainlog.info('Sigma step was not float or integer. Setting it to default value 0.5.\n')
			sigma_step=0.5
			for i in range(len(lines)):
				if 'sigma_step' in lines[i]:
					lines[i]='sigma_step\t\t=\t'+str(0.5)+'\n'
		#
		try:
			residual_frac=float(residual_frac)
			if residual_frac>=0.5:
				mainlog.info('residual_frac was greater than or equal to 0.5. Setting it to default value 0.1.\n')
				residual_frac=0.1
				for i in range(len(lines)):
					if 'residual_frac' in lines[i]:
						lines[i]='residual_frac\t=\t'+str(0.1)+'\n'
			elif residual_frac==0:
				mainlog.info('residual_frac was 0. Setting it to default value 0.1.\n')
				residual_frac=0.1
				for i in range(len(lines)):
					if 'residual_frac' in lines[i]:
						lines[i]='residual_frac\t=\t'+str(0.1)+'\n'
		except:
			mainlog.info('residual_frac was not float or integer. Setting it to default value 0.1.\n')
			residual_frac=0.1
			for i in range(len(lines)):
				if 'residual_frac' in lines[i]:
					lines[i]='residual_frac\t=\t'+str(0.1)+'\n'
		#
		try:
			min_sigma=float(min_sigma)
			if min_sigma<5:
				mainlog.info('min_sigma is less than 5. Setting it to default value 5.0.\n')
				min_sigma=5.0
				for i in range(len(lines)):
					if 'min_sigma' in lines[i]:
						lines[i]='min_sigma\t\t=\t'+str(5.0)+'\n'
		except:
			mainlog.info('min_sigma was not float or integer. Setting it to default value 8.0.\n')
			min_sigma=8.0
			for i in range(len(lines)):
				if 'min_sigma' in lines[i]:
					lines[i]='min_sigma\t\t=\t'+str(8.0)+'\n'

		#
		uv_format_err=False
		if uvrange_to_cal!='':
			if 'lambda' not in uvrange_to_cal:
				uv_format_err=True
			else:
				uv_list=uvrange_to_cal.split('lambda')
				if uv_list[-1]!='':
					uv_format_err=True
				else:
					uv_str=uv_list[0]
					if uv_str[-1]=='k':
						uv_str=uv_str[:-1]
					if '~' not in uv_str and '>' not in uv_str and '<' not in uv_str:
						uv_format_err=True
					elif '~' in uv_str:
						print (uv_format_err)
						uv_str_list=uv_str.split('~')
						if len(uv_str_list)!=3:
							uv_format_err=True
						else:
							try:
								a=float(uv_str_list[0])
								b=float(uv_str_list[1])
							except:
								uv_format_err=True
					elif '>' in uv_str:
						uv_str_list=uv_str.split('>')
						if len(uv_str_list)!=2:
							uv_format_err=True
						else:
							try:
								a=float(uv_str_list[1])
							except:
								uv_format_err=True
					elif '<' in uv_str:
						uv_str_list=uv_str.split('<')
						if len(uv_str_list)!=2:
							uv_format_err=True
						else:
							try:
								a=float(uv_str_list[1])
							except:
								uv_format_err=True

			if uv_format_err==True:
				uvrange_to_cal=''
				mainlog.info('uvrange format is wrong. Setting it to default value.\n')
				for i in range(len(lines)):
					if 'uvrange_to_cal' in lines[i]:
						lines[i]='uvrange_to_cal\t=\t\''+str(uvrange_to_cal)+'\'\n'
		
		###
		try:
			gain_minsnr=float(gain_minsnr)
			if gain_minsnr<3.0:
				mainlog.info('Gain SNR is less than 3.0. Setting it to default 3.0.\n')
				gain_minsnr=3.0
				for i in range(len(lines)):
					if 'gain_minsnr' in lines[i]:
						lines[i]='gain_minsnr\t\t=\t'+str(3.0)+'\n'
		except:
			mainlog.info('gain_minsnr was not float or integer. Setting it to default value 3.0.\n')
			gain_minsnr=3.0
			for i in range(len(lines)):
				if 'gain_minsnr' in lines[i]:
					lines[i]='gain_minsnr\t\t=\t'+str(3.0)+'\n'
		###

		try:
			if type(DR_delta_rms)!=int or type(DR_delta_rms)!=float:
				DR_delta_rms=int(DR_delta_rms)
		except:
			mainlog.info('DR_delta_rms could not be convereted to integer. Setting it to default value 20.\n')
			DR_delta_rms=20
			for i in range(len(lines)):
				if 'DR_delta_rms' in lines[i]:
					lines[i]='DR_delta_rms\t=\t'+str(20)+'\n'
		try:
			if type(DR_delta_neg)!=int or type(DR_delta_neg)!=float:
				DR_delta_neg=int(DR_delta_neg)
		except:
			DR_delta_neg=10
			if DR_delta_neg>DR_delta_rms:
				mainlog.info('DR_delta_neg is greater than DR_delta_rms. Setting it to DR_delta_rms.\n')
				DR_delta_neg=DR_delta_rms
				for i in range(len(lines)):
					if 'DR_delta_neg' in lines[i]:
						lines[i]='DR_delta_neg\t=\t'+str(DR_delta_rms)+'\n'
			else:
				mainlog.info('DR_delta_neg was not an integer. Setting it to default value 10.\n')
				for i in range(len(lines)):
					if 'DR_delta_neg' in lines[i]:
						lines[i]='DR_delta_neg\t=\t'+str(10)+'\n'
		###
		
		if type(max_DR)!=int or type(max_DR)!=float:
			try:
				max_DR=int(max_DR)
				if max_DR==0:
					mainlog.info('max_DR = 0. Setting it to default value 10000.\n')
					max_DR=10000
					for i in range(len(lines)):
						if 'max_DR' in lines[i]:
							lines[i]='max_DR\t\t\t=\t'+str(10000)+'\n'
			except:
				mainlog.info('max_DR was not an integer or float or could not be converted to integer/float. Setting it to default value 10000.\n')
				max_DR=10000
				for i in range(len(lines)):
					if 'max_DR' in lines[i]:
						lines[i]='max_DR\t\t\t=\t'+str(10000)+'\n'

		if type(min_DR)!=int or type(min_DR)!=float:
			try:
				min_DR=int(min_DR)
				if max_DR==0:
					mainlog.info('min_DR = 0. Setting it to default value 20.\n')
					min_DR=20
					for i in range(len(lines)):
						if 'min_DR' in lines[i]:
							lines[i]='min_DR\t\t\t=\t'+str(20)+'\n'
			except:
				min_DR=20
				if min_DR>max_DR:
					mainlog.info('min_DR is greater than max_DR. Setting it to default 20.\n')
					min_DR=20
					for i in range(len(lines)):
						if 'min_DR' in lines[i]:
							lines[i]='min_DR\t\t\t=\t'+str(min_DR)+'\n'
				else:
					mainlog.info('min_DR was not an integer or float or could not be converted to integer/float. Setting it to default value 20.\n')
					for i in range(len(lines)):
						if 'min_DR' in lines[i]:
							lines[i]='min_DR\t\t\t=\t'+str(20)+'\n'
		
		if type(min_selfcal_snr)!=int or type(min_selfcal_snr)!=float:
			try:
				min_selfcal_snr=int(min_selfcal_snr)
				if min_selfcal_snr<3:
					mainlog.info('min_selfcal_snr is less than 3. Setting it to default value 3.\n')
					min_selfcal_snr=3
					for i in range(len(lines)):
						if 'min_selfcal_snr' in lines[i]:
							lines[i]='min_selfcal_snr\t=\t'+str(min_selfcal_snr)+'\n'
			except:
				min_selfcal_snr=3
				mainlog.info('min_selfcal_snr was not an integer or float or could not be converted to integer/float. Setting it to default value 3.\n')
				for i in range(len(lines)):
					if 'min_selfcal_snr' in lines[i]:
						lines[i]='min_selfcal_snr\t=\t'+str(min_selfcal_snr)+'\n'

	# Validating flagging parameters
	################################

	if type(want_uvsub_flag)!=bool:
		mainlog.info('want_uvsub_flag was not an boolean. Setting it to default value False.\n')
		want_uvsub_flag=False
		for i in range(len(lines)):
			if 'want_uvsub_flag' in lines[i]:
				lines[i]='want_uvsub_flag\t=\t'+str(False)+'\n'
	if type(use_ankflagger)!=bool:
		mainlog.info('use_ankflagger was not an boolean. Setting it to default value False.\n')
		use_ankflagger=False
		for i in range(len(lines)):
			if 'use_ankflagger' in lines[i]:
				lines[i]='use_ankflagger\t=\t'+str(False)+'\n'
	elif want_uvsub_flag==True:
		try:
			from aNKflag import runank
		except:
			mainlog.info('aNKflag can not be imported. Not using aNKflagger.\n')
			use_ankflagger=False
			for i in range(len(lines)):
				if 'use_ankflagger' in lines[i]:
					lines[i]='use_ankflagger\t=\t'+str(False)+'\n'

	# Validating image export options
	#################################
	if savedir!='':
		if os.path.isdir(savedir)==False:
			mainlog.info('savedir is not present. Setting it to basedir.\n')
			savedir=basedir
		for i in range(len(lines)):
			if 'savedir' in lines[i]:
				lines[i]='savedir\t\t\t=\t'+str(savedir)+'\n'

	if type(savemodel)!=bool:
		mainlog.info('savemodel was not boolean. Set it default value False.\n')
		for i in range(len(lines)):
			if 'savemodel' in lines[i]:
				lines[i]='savemodel\t\t=\t'+str(False)+'\n'
		
	if type(saveresidual)!=bool:
		mainlog.info('saveresidual was not boolean. Set it default value False.\n')
		for i in range(len(lines)):
			if 'saveresidual' in lines[i]:
				lines[i]='saveresidual\t=\t'+str(False)+'\n'

	# Writing corrected parameters
	##############################

	inpfil.seek(0)
	inpfil.writelines(lines)
	inpfil.close()
	os.system('rm -rf casa*log')
	mainlog.info('Input file validation complete.\n')	
	







