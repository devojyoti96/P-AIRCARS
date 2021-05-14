import os,psutil,sys,struct,copy,numpy as np
from datetime import datetime 
from casatools import *
from casatasks import *
from paircars.access_ms import *
from paircars.basic_func import *

'''
Code is written by Devojyoti Kansabanik, 16 Feb, 2021; 22 Feb, 2021
'''

class CALIBRATE():
	def __init__(self):
		pathname=os.path.dirname(os.path.realpath(__file__))
		self.path=pathname
		os.system('rm -rf casa*log')
	
	def calibrate(self,**kwargs):
		'''
		Function to perform calibration with the full-Jones Mitchcal algorithm
		Parameters:
		Keyword arguments:
		msname = Name of the measurement set
		caltable = Name of the Jones calibration table
		calmode = scalar/diag/rotation. Type of Jones matrix. If nothing is given code will assume generalised Jones matrix
		quiet = True, Do not print verbose output (default : True)
		absmem = Amount of memory to use (in GB)
		minuv = Minimum uv distance to use for calibration in meter
		maxuv = Maximum uv distance to use for calibration in meter
		gaintable = Previous gaintables to apply (Either CASA gaintable or CALIBRATE gaintable)
		a = minimum_accuracy,stopping_accuracy
		i = Number of calibration iteration (integer)
		j = Number of cpu threads to use (integer)
		solmode = 'R' for robust calibration
		rmsthresh = [] List of rms threshold based flagging 
		'''
		cwd=os.getcwd()
		os.chdir(self.path+'/calibrate_tools')
		kwords=kwargs.keys()
		arg_str=''
		if ('msname' not in kwords or 'caltable' not in kwords) or (kwargs['msname']=='' or kwargs['caltable']==''):
			print ('Argument msname or caltable is missing.\n')
			return 
		else:
			msname=kwargs['msname'] # MS name
			if os.path.isdir(msname)==False:
				print ('Measurement set is not present.\n')
				return
			else:
				ms_dirname=os.path.dirname(os.path.realpath(msname))
				if os.path.isdir(msname.split('.ms')[0]+'.temp.ms'):
					os.system('rm -rf '+msname.split('.ms')[0]+'.temp.ms')
				os.system('cp -r '+msname+' '+msname.split('.ms')[0]+'.temp.ms')
				msname=msname.split('.ms')[0]+'.temp.ms'
			caltable=kwargs['caltable'] # Caltable name
			
			AM=AccessMS(msname)
			freqs=AM.get_freqs()/10**6
			start_freq=freqs[0]
			end_freq=freqs[-1]
			mjdsecs=AM.get_timestamps_in_mjdsecs()
			startmjd=mjdsecs[0]
			endmjd=mjdsecs[-1]
			nchan=AM.get_num_channels()
			ntime=AM.get_num_timestamps()

			if 'calmode' in kwords:
				if kwargs['calmode']=='scalar':
					arg_str+=' -scalar '
				elif kwargs['calmode']=='diag':
					arg_str+=' -diag '
				elif kwargs['calmode']=='rotation':
					arg_str+=' -rotation '

			if 'quiet' in kwords:
				if kwargs['quiet']!=False:
					arg_str+=' -quiet '
			else:
				arg_str+=' -quiet'

			if 'absmem' in kwords: # Absolute memory usuage
				try:
					absmem=float(kwargs.get('absmem'))
					avail_memory=int(psutil.virtual_memory().available/10**9) # In GB o nearest integer
					if absmem>0.5*avail_memory:
						absmem=int(0.5*avail_memory)
				except:
					avail_memory=int(psutil.virtual_memory().available/10**9) # In GB o nearest integer
					absmem=int(0.5*avail_memory)
				arg_str+=' -absmem '+str(absmem)

			if 'minuv' in kwords: # Minimum uv distance to use in meter
				try:
					minuv=float(kwargs.get('minuv'))
					arg_str+=' -minuv '+str(minuv)
				except:
					print ('Wrong minuv format.')
			if 'maxuv' in kwords: # Maximum uv distance in meter
				try:
					maxuv=float(kwargs.get('maxuv'))
					arg_str+=' -maxuv '+str(maxuv)
				except:
					print ('Wrong maxuv format.')	
			if 'a' in kwords: # Minimum and stopping accuracy
				try:
					a=kwargs.get('a').split(',')
					if len(a)>2:
						print ('Wrong format of -a')
						a0=''
						a1=''
					else:
						a0=a[0]
						a1=a[1]
				except:
					print ('Wrong format of -a')
					a0=''
					a1=''
				if a0!='' and a1!='':
					arg_str+=' -a '+str(a0)+' '+str(a1)
			if 'i' in kwords: # Number of iterations
				try:
					i=int(kwargs.get('i'))
					arg_str+=' -i '+str(i)
				except:
					print ('Wrong format of -i')		
			if 'j' in kwords: # Number of threads
				try:
					i=int(kwargs.get('j'))
					total_threads=psutil.cpu_count()
					if i>int(0.5*total_threads):
						i=int(0.5*total_threads)
					arg_str+=' -j '+str(i)
				except:
					print ('Wrong format')	
			if 'm' in kwords: # Model
				m=kwargs.get('m')
				if os.path.isfile(m)==False:
					print ('Model file is not present.')
				else:
					arg_str+=' -m '+str(m)
			if 't' in kwords: # Timesteps
				try:
					t=int(kwargs.get('t'))
					arg_str+=' -t '+str(t)
					ntime=int(ntime/t)
				except:
					print ('Wrong format of -t')
			if 'interval' in kwords: # Time interval
				try:
					interval=kwargs.get('interval').split(',')
					if len(interval)>2:
						print ('Wrong format of -interval')
						int0=''
						int1=''
					else:
						int0=interval[0]
						int1=interval[1]
				except:
					print ('Wrong format of -interval')
					int0=''
					int1=''
				if int0!='' and int1!='':
					arg_str+=' -interval '+str(int0)+' '+str(int1)
			if 'ch' in kwords: # Channel resolution
				try:
					ch=int(kwargs.get('ch'))
					arg_str+=' -ch '+str(ch)
					nchan=int(nchan/ch)
				except:
					print ('Wrong format of -ch')
			if 'gaintable' in kwords: # Previous gaintable
				gaintable=kwargs.get('gaintable')
				if len(gaintable)>0:
					casa_gaintable=[]
					CALIBRATE_gaintable=[]
					for table in gaintable:
						if os.path.isfile(table)==False:
							if os.path.isdir(table)==True:
								casa_gaintable.append(table)
							else:
								print ('Gaintable '+table+' not found\n')
								continue 
						else:
							CALIBRATE_gaintable.append(table)
					if len(casa_gaintable)>0:
						print ('Applying CASA gaintables :',casa_gaintable)
						applycal(vis=msname,gaintable=casa_gaintable,flagbackup=True)	
					if len(CALIBRATE_gaintable)>0:
						for table in CALIBRATE_gaintable:
							print ('Applying CALIBRATE gaintable :'+table+'..........')
							applycal_result=self.applycal(msname=msname,gaintable=table,datacolumn='CORRECTED_DATA')
							if applycal_result!=0:
								print ('Error in applying CALIBRATE gaintable '+table)
								return
			else:
				gaintable=[]
			if 'datacolumn' in kwords: # Datacolumn
				datacolumn=kwargs.get('datacolumn')
				if len(gaintable)>0:
					arg_str+=' -datacolumn CORRECTED_DATA'
					datacolumn='CORRECTED_DATA'
				elif datacolumn=='DATA' or datacolumn=='CORRECTED_DATA':
					arg_str+=' -datacolumn '+str(datacolumn)
				else:
					print ('Wrong datacolumn.')
			elif len(gaintable)>0:
				arg_str+=' -datacolumn CORRECTED_DATA'
				datacolumn='CORRECTED_DATA'
			else:
				arg_str+=' -datacolumn DATA'
				datacolumn='DATA'
			if ('solmode' not in kwords or 'rmsthresh' not in kwords) or (kwargs['solmode']!='R' or len(kwargs['rmsthresh'])==0): 
				print ('./calibrate '+arg_str+' '+msname+' '+caltable)
				os.system('./calibrate '+arg_str+' '+msname+' '+caltable)
				bin_data=np.fromfile(caltable,dtype=np.float64)
				np.save(caltable+'.temp',np.array([bin_data,start_freq,end_freq,startmjd,endmjd,nchan,ntime],dtype='object'))
				os.system('rm -rf '+caltable)
				os.system('mv '+caltable+'.temp.npy '+caltable)	
			elif kwargs['solmode']=='R' and len(kwargs['rmsthresh'])!=0:
				solmode=kwargs['solmode']
				rmsthresh=kwargs['rmsthresh']
				for rms in rmsthresh:
					c=0
					print ('Calibrating and flagging on threshold :'+str(rms)+' sigma\n')
					bad_ants=self.get_num_flag_baselines(msname,flagfrac=0.8)
					flagdata(vis=msname,antenna=bad_ants)
					print ('flagdata(vis=\''+msname+',antenna=\''+bad_ants+'\')')
					while c==0:
						print ('./calibrate '+arg_str+' '+msname+' '+caltable)
						os.system('./calibrate '+arg_str+' '+msname+' '+caltable)
						bin_data=np.fromfile(caltable,dtype=np.float64)
						np.save(caltable+'.temp',np.array([bin_data,start_freq,end_freq,startmjd,endmjd,nchan,ntime],dtype='object'))
						os.system('rm -rf '+caltable)
						os.system('mv '+caltable+'.temp.npy '+caltable)
						self.applycal(msname=msname,gaintable=caltable,datacolumn=datacolumn,applymode='calflag',flagbackup=False)
						num_flag,flag_fraction=self.flagger(msname,float(rms))
						if int(num_flag)==0:
							c=1				
			if os.path.isdir(msname.split('.ms')[0]+'.temp.ms'):
				os.system('rm -rf '+msname.split('.ms')[0]+'.temp.ms*')
			os.chdir(cwd)
			os.system('rm -rf casa*log')
		return caltable

	def applycal(self,**kwargs): # TODO: Add feature of applying multiple CASA and CALIBRATE gaintables
		'''
		Function to apply CALIBRATE solution (Right now only one CALIBRATE caltable accepts)
		Keyword arguments:
		datacolumn = Datacolumn to apply solution
		msname = Name of the measurement set
		gaintable = Name of the CALIBRATE gaintable
		applymode = 'calflag' or 'calonly'
		flagbackup = True, keep flag backup or not (default : True)
		'''
		cwd=os.getcwd()
		os.chdir(self.path+'/calibrate_tools')
		kwords=kwargs.keys()
		if ('msname' not in kwords or 'gaintable' not in kwords) or (kwargs['msname']=='' or kwargs['gaintable']==''):
			print ('Argument msname or gaintable is missing.\n')
			return 1
		else:
			arg_str=''
			msname=kwargs['msname'] # MS name
			if os.path.isdir(msname)==False:
				print ('Measurement set is not present.\n')
				return 1
			md=msmetadata()
			md.open(msname)
			nchan=md.nchan(0) # TODO :for multiple spw
			ntime=len(md.timesforfield(0))
			md.close()
			gaintable=kwargs['gaintable'] # Gaintable name
			if (os.path.isfile(gaintable)==False) and (os.path.isdir(gaintable)==False):
				print ('Gaintable is not found.')
				return 1
			if 'datacolumn' in kwords: # Datacolumn 
				datacolumn=kwargs.get('datacolumn')
				if datacolumn=='DATA' or datacolumn=='CORRECTED_DATA':
					arg_str+=' -datacolumn '+str(datacolumn)
				else:
					print ('Wrong datacolumn.')
			else:
				datacolumn='DATA'
				arg_str+=' -datacolumn DATA'
			if 'applymode' not in kwords:
				applymode='calflag'
			else:
				applymode=kwargs['applymode']

			if 'flagbackup' in kwords:
				flagbackup=kwargs['flagbackup']
			else:
				flagbackup=True
			gaintable_path=os.path.dirname(os.path.realpath(gaintable))
			result=self.modify_caltable_for_ms(msname,gaintable,gaintable+'.temp_nchan_ntime.bin')
			if result=='Nosol':
				os.system('rm -rf casa*log '+gaintable+'.temp_nchan_ntime.bin')
				os.chdir(cwd)
				os.system('rm -rf casa*log')
				return 'Nosol'
			gaintable=gaintable+'.temp_nchan_ntime.bin'
			if applymode=='calflag':
				if flagbackup==True:
					af=agentflagger()
					af.open(msname)
					versionlist=af.getflagversionlist()
					if len(versionlist)!=0:
						for version_name in versionlist:
							if 'CALIBRATE_applycal' in version_name:
								version_num=int(version_name.split(':')[0].split(' ')[0].split('_')[-1])+1
							else:
								version_num=1
					else:
						version_num=1
					now = datetime.now()
					dt_string = now.strftime("%d-%m-%Y %H:%M:%S")
					af.saveflagversion('CALIBRATE_applycal_'+str(version_num),'Flags autosave on '+dt_string)
					af.done()
				print ('./applysolutions '+arg_str+' '+msname+' '+gaintable)
				os.system('./applysolutions '+arg_str+' '+msname+' '+gaintable)
				tb=table()	
				tb.open(msname,nomodify=False)
				data=tb.getcol('DATA')
				cordata=tb.getcol('CORRECTED_DATA')
				flag=tb.getcol('FLAG')
				pos=np.isnan(cordata)
				cordata[pos]=data[pos]
				flag[pos]=True
				tb.putcol('FLAG',flag)
				tb.putcol('CORRECTED_DATA',cordata)
				tb.flush()
				tb.close()
			else:		
				print ('./applysolutions '+arg_str+' '+msname+' '+gaintable)
				os.system('./applysolutions '+arg_str+' '+msname+' '+gaintable)
			os.system('rm -rf casa*log '+gaintable+'.temp_nchan_ntime.bin')
		os.chdir(cwd)
		os.system('rm -rf casa*log')
		return 0
			
	def convert_gaintable_bin2npy(self,gaintable,outputfile):
		'''
		Function to convert CALIBRATE binary gaintable to numpy table
		Parameter:
		gaintable = Name of the CALIBRATE binary gaintable
		outputfile - Name of the numpy gaintable
		Return:
		A numpy table. Format [Gaintable header, Jones matrices, Binary header, Jones data]
		Gaintable header = ['MWAOCAL',0,number of intervals,number of antennas,number of channels,number of polarisation]
		Jones matrices array shape = [2,2,number_of_intervals,number_of_antennas,number_of_channels,8(4 polarisation,real and imaginary part)]
		Jones data shape = [number_of_intervals,number_of_antennas,number_of_channels,8(4 polarisation,real and imaginary part)]
		'''
		cwd=os.getcwd()
		gaintable_path=os.path.dirname(os.path.realpath(gaintable))
		outfile_path=os.path.dirname(outputfile)
		data=open(gaintable,'rb').read()
		data_header=data[:48]
		struct_unpack=struct.Struct('sssssss').unpack_from	# Header string	
		try:
			header_intro=''.join(struct_unpack(data_header[:7])) # Header intro 'MWAOCAL'
		except TypeError:
			header_intro=str(data_header[:7],'utf-8')
		struct_unpack=struct.Struct('i').unpack_from      
		filetype=struct_unpack(data_header[8:12])[0] # Filetype, always 0
		num_intervals=struct_unpack(data_header[16:20])[0] # Number of time intervals 
		num_antenna=struct_unpack(data_header[20:24])[0] # Number of antenna
		num_channels=struct_unpack(data_header[24:28])[0] # Number of channels
		num_polarisation=struct_unpack(data_header[28:32])[0] # Number of polarisation
		jones_header=np.array([header_intro,filetype,num_intervals,num_antenna,num_channels,num_polarisation])
		jones_data=np.fromfile(gaintable,dtype=np.float64)
		jones_data_copy=copy.deepcopy(jones_data)
		jones_data=jones_data[6:]
		# Structuring the data		
		jones_data=jones_data.reshape(num_intervals,num_antenna,num_channels,num_polarisation*2) # Shape = ninterval,nant,nchan,pol*2
		jones_data_array=np.array([[jones_data[:,:,:,0]+1j*jones_data[:,:,:,1],jones_data[:,:,:,2]+1j*jones_data[:,:,:,3]],\
				[jones_data[:,:,:,4]+1j*jones_data[:,:,:,5],jones_data[:,:,:,6]+1j*jones_data[:,:,:,7]]]) # Shape = [2,2,ninterval,nant,nchan]
		if outfile_path=='':
			np.save(gaintable_path+'/'+outputfile,np.array([jones_header,jones_data_array,jones_data_copy[:6],jones_data],dtype='object'))
			os.chdir(cwd)
			os.system('rm -rf casa*log')
			return gaintable_path+'/'+outputfile+'.npy'
		elif os.path.isdir(outfile_path)==False:
			print ('Output file directory does not exists. Saving into gaintbale directory.\n')
			np.save(gaintable_path+'/'+os.path.basename(outputfile),np.array([jones_header,jones_data_array,jones_data_copy[:6],jones_data],dtype='object'))
			os.chdir(cwd)
			os.system('rm -rf casa*log')
			return gaintable_path+'/'+os.path.basename(outputfile)+'.npy'
		else:
			np.save(outputfile,np.array([jones_header,jones_data_array,jones_data_copy[:6],jones_data],dtype='object'))
			os.chdir(cwd)
			os.system('rm -rf casa*log')
			return outputfile+'.npy'	

	def convert_gaintable_npy2bin(self,gaintable,outputfile,remove_nan=True):
		'''
		Function to convert CALIBRATE numpy gaintable to binary file
		Parameter:
		gaintable = Name of the CALIBRATE numpy gaintable
		outputfile = Name of the binary gaintable
		remove_nan = True, remove nan solutions and replace with identity
		Return:
		A CALIBRATE binary gaintable.
		'''	
		cwd=os.getcwd()
		gaintable_path=os.path.dirname(os.path.realpath(gaintable))
		outfile_path=os.path.dirname(outputfile)
		numpy_table=np.load(gaintable,allow_pickle=True)
		bin_header=numpy_table[2]
		data=numpy_table[1]
		header=numpy_table[0]
		nint=int(header[2])
		nant=int(header[3])
		nchan=int(header[4])
		numpy_data=np.empty([nint,nant,nchan,8])
		bad_flags=[]
		for i in range(nint):
			for j in range(nant):
				for k in range(nchan):
					if remove_nan==True:
						if np.sum(np.isnan(data[:,:,i,j,k]))==4:
							data[:,:,i,j,k]=np.array([[1,0],[0,1]])
							bad_flags.append([i,j,k])
					data_re=np.real(data[:,:,i,j,k].flatten())
					data_im=np.imag(data[:,:,i,j,k].flatten())
					numpy_data[i,j,k,:]=np.insert(data_im, np.arange(len(data_re)),data_re)
		numpy_data=numpy_data.flatten()
		final_data_header=np.append(bin_header,numpy_data)
		if outfile_path=='':
			final_data_header.tofile(gaintable_path+'/'+outputfile,format='np.float64')
			os.chdir(cwd)
			os.system('rm -rf casa*log')
			return gaintable_path+'/'+outputfile,bad_flags
		elif os.path.isdir(outfile_path)==False:
			print ('Output file directory does not exists. Saving into gaintbale directory.\n')
			final_data_header.tofile(gaintable_path+'/'+os.path.basename(outputfile),format='np.float64')
			os.chdir(cwd)
			os.system('rm -rf casa*log')
			return gaintable_path+'/'+os.path.basename(outputfile),bad_flags
		else:
			final_data_header.tofile(outputfile,format='np.float64')
			os.chdir(cwd)
			os.system('rm -rf casa*log')
			return outputfile,bad_flags

	def modify_caltable_for_ms(self,msname,caltable,outputname):
		'''
		Function to modify the CALIBRATE caltable to apply on a ms
		Parameters:
		msname = Name of the measurement set
		caltable = Name of the numpy caltable
		outputname= Name of the output caltable in CALIBRATE binary format
		Return:
		Modified caltable		
		'''
		print ('Arranging solutions to apply on the ms.........\n')
		cwd=os.getcwd()
		gaintable_path=os.path.dirname(os.path.realpath(caltable))
		outfile_path=os.path.dirname(outputname)
		bin_data,cal_start_freq,cal_end_freq,cal_startmjd,cal_endmjd,cal_nchan,cal_ntime=np.load(caltable,allow_pickle=True)
		AM=AccessMS(msname)
		ntime=AM.get_num_timestamps()
		nchan=AM.get_num_channels()
		nant=AM.get_num_antenna()
		times=AM.get_timestamps_in_mjdsecs()
		startmjd=times[0]
		endmjd=times[-1]
		freqs=AM.get_freqs()/10**6
		startfreq=float(freqs[0])
		endfreq=float(freqs[-1])
		bin_data=bin_data.astype('float64')
		bin_data.tofile(caltable+'.temp.bin',format='np.float64')
		npyfile=self.convert_gaintable_bin2npy(caltable+'.temp.bin',caltable+'.CALIBRATE_temp')
		numpy_table=np.load(npyfile,allow_pickle=True)
		data=numpy_table[3]
		new_data=np.empty((ntime,nant,nchan,8))
		bad_calchantime=[]
		bin_header = struct.pack("8s",b"MWAOCAL\n")+struct.pack("i",0)+struct.pack("i",0)+struct.pack("i",int(ntime))+struct.pack("i",int(nant))+struct.pack("i",int(nchan))+\
				struct.pack("i",4)+struct.pack("d",0.0)+struct.pack("d",0.0)
		cal_nchan=data.shape[2]
		cal_ntime=data.shape[0]
		cal_time_res=abs(cal_endmjd-cal_startmjd)/cal_ntime
		cal_freq_res=(cal_end_freq-cal_start_freq)/cal_nchan
		if cal_time_res>0:
			cal_times=np.arange(cal_startmjd,cal_endmjd,cal_time_res)
		else:
			cal_times=np.array([cal_startmjd])
		if cal_freq_res>0:
			cal_freqs=np.arange(cal_start_freq,cal_end_freq,cal_freq_res)
		elif cal_freq_res==0:
			cal_freqs=np.array([cal_start_freq])
		else:
			cal_freqs=np.arange(cal_end_freq,cal_start_freq,cal_freq_res)
		for i in range(data.shape[0]):
			for j in range(data.shape[2]):
				caltime=cal_times[i]
				calfreq=cal_freqs[j]
				if np.sum(np.isnan(data[i,:,j,:]))==4:
					bad_calchantime.append([caltime,calfreq])
		bad_calchantime=np.array(bad_calchantime)
		if len(bad_calchantime)==cal_ntime*cal_nchan:
			print ('No unflagged solutions in the caltable.\n')
			os.system('rm -rf casa*log '+caltable+'.temp.bin '+caltable+'.CALIBRATE_temp*')
			os.chdir(cwd)
			os.system('rm rf casa*log')
			return 'Nosol'
		for i in range(ntime):
			for j in range(nchan):
				mstime=times[i]
				msfreq=freqs[j]
				while True:
					min_cal_chan=min(getnearpos(np.array(cal_freqs),msfreq))
					max_cal_chan=max(getnearpos(np.array(cal_freqs),msfreq))
					if j>=max_cal_chan:
						nearest_cal_chan=max_cal_chan
					else:
						nearest_cal_chan=min_cal_chan
					min_cal_time=min(getnearpos(np.array(cal_times),mstime))
					max_cal_time=max(getnearpos(np.array(cal_times),mstime))
					if i>=max_cal_time:
						nearest_cal_time=max_cal_time
					else:
						nearest_cal_time=min_cal_time
					if [nearest_cal_time,nearest_cal_chan] in bad_calchantime:
						print ('Bad cal chantime encountered.\n')
						continue
					else:
						new_data[i,:,j,:]=data[nearest_cal_time,:,nearest_cal_chan,:]
						break	
		np.save(caltable+'.test',new_data)	
		new_data_flattened=new_data.flatten(order='C')
		if outfile_path=='':
			outfile=gaintable_path+'/'+outputname
		elif os.path.isdir(outfile_path)==False:
			outfile=gaintable_path+'/'+os.path.basename(outputname)
		else:
			outfile=outputname
		fil=open(outfile,'wb')
		fil.write(bin_header)
		fil.close()
		with open(outfile,mode='ba+') as f:
			new_data_flattened.tofile(f,format='np.float64')
		os.system('rm -rf '+npyfile)
		os.system('rm -rf casa*log '+caltable+'.temp.bin '+caltable+'.CALIBRATE_temp*')
		os.chdir(cwd)
		os.system('rm -rf casa*log')
		return outfile

	def flagger(self,msname,rms):
		'''
		Function to flag real and imaginary part of visibility based on rms threshold
		Parameters:
		msname = Name of the measurement set
		rmsthresh = RMS threshold for n-sigma flagging
		Return:
		Total flag points
		'''
		cwd=os.getcwd()
		mstool=ms()
		mstool.open(msname,nomodify=False)
		res_data=mstool.getdata('residual_data')['residual_data']
		flag=mstool.getdata('flag')
		flag_data=flag['flag']
		pos_flag=np.where(flag_data==True)
		num_flag=np.nansum(flag_data)
		res_data[pos_flag]=np.nan
		xx=res_data[0,:,:]
		xx_flag=flag_data[0,:,:]
		xy=res_data[1,:,:]
		xy_flag=flag_data[1,:,:]
		yx=res_data[2,:,:]
		yx_flag=flag_data[2,:,:]
		yy=res_data[3,:,:]
		yy_flag=flag_data[3,:,:]

		sigma_xxre=np.nanstd(np.real(xx))
		sigma_xxim=np.nanstd(np.imag(xx))
		sigma_xyre=np.nanstd(np.real(xy))
		sigma_xyim=np.nanstd(np.imag(xy))
		sigma_yxre=np.nanstd(np.real(yx))
		sigma_yxim=np.nanstd(np.imag(yx))
		sigma_yyre=np.nanstd(np.real(yy))
		sigma_yyim=np.nanstd(np.imag(yy))

		pos_xxre=np.where(np.real(xx)>rms*sigma_xxre)
		pos_xxim=np.where(np.imag(xx)>rms*sigma_xxim)
		pos_xyre=np.where(np.real(xy)>rms*sigma_xyre)
		pos_xyim=np.where(np.imag(xy)>rms*sigma_xyim)
		pos_yxre=np.where(np.real(yx)>rms*sigma_yxre)
		pos_yxim=np.where(np.imag(yx)>rms*sigma_yxim)
		pos_yyre=np.where(np.real(yy)>rms*sigma_yyre)
		pos_yyim=np.where(np.imag(yy)>rms*sigma_yyim)

		pos_xx=np.append(pos_xxre,pos_xxim)
		pos_xy=np.append(pos_xyre,pos_xyim)
		pos_yx=np.append(pos_yxre,pos_yxim)
		pos_yy=np.append(pos_yyre,pos_yyim)
		pos=np.unique(np.append(np.append(np.append(pos_xx,pos_xy),pos_yx),pos_yy))
		for i in pos:
				flag_data[:,:,i]=True
		flag['flag']=flag_data
		final_num_flag=np.nansum(flag_data)-num_flag
		new_flag_fraction=final_num_flag/num_flag
		mstool.putdata(flag)
		mstool.close()
		os.chdir(cwd)
		os.system('rm -rf casa*log')
		return final_num_flag,new_flag_fraction

	def get_num_flag_baselines(self,msname,flagfrac=0.5):
		'''
		Function to get the antennas for which a certain amount of data are flagged
		Parameters:
		msname = Name of the measurement set
		flagfrac = Fraction of data flagged to consider the antennas as bad (default : 0.5)
		Return:
		Bad antenna strings
		'''
		cwd=os.getcwd()
		mstool=ms()
		msmd=msmetadata()
		msmd.open(msname)
		nant=msmd.nantennas()
		msmd.close()
		ant_flag_percent={}
		for i in range(nant):
			mstool.open(msname)
			mstool.select({'antenna1':i})
			flag=mstool.getdata('FLAG')['flag']
			mstool.close()
			mstool.open(msname)
			mstool.select({'antenna2':i})
			flag=np.append(flag,mstool.getdata('FLAG')['flag'])
			mstool.close()
			ant_flag_percent[i]=(np.sum(flag.flatten())/float(len(flag.flatten())))
		bad_ants=''
		for i in range(nant):
			if ant_flag_percent[i]>flagfrac:
				bad_ants+=str(i)+','
		bad_ants=bad_ants[:-1]
		os.chdir(cwd)
		os.system('rm -rf casa*log')
		return bad_ants








		

