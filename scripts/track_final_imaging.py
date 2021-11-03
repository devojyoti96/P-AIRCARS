import os,glob,copy,time
from astropy.io import fits
from casatasks import *
from casatools import msmetadata,table,measures,quanta,agentflagger,image,calibrater,ms
from paircars.libpaircars import send_paircars_notification
from optparse import OptionParser


usage= 'Track and notifier final imaging'
parser = OptionParser(usage=usage)
parser.add_option('--OBSID',dest="obsid",default=0,help="Observation ID",metavar="Integer")
parser.add_option('--basedir',dest="basedir",default=None,help="Name of the base directory",metavar="Directory name")
parser.add_option('--savedir',dest="savedir",default=None,help="Name of the directory where images are saved",metavar="Directory name")
parser.add_option('--email',dest="email",default=None,help="E-mail ID to send notification",metavar="E-mail ID")
parser.add_option('--num_ms',dest="num_ms",default=0,help="Total number of measurement sets being imaged",metavar="Integer")
(options, args) = parser.parse_args()

def get_quicklook_image(imagename,outfile,freq,timestamp,field_of_view=2): 
	'''
	Function to get a quick look image
	Parameters:
	imagename = Name of the CASA image
	outfile = Output file name
	freq = Frequency in MHz
	timestamp = Timestamp string
	field_of_view = Field of view to cut the image in degree (default : 2)
	Return:
	Outfile name
	'''
	os.system('cp -r '+imagename+' '+'quick_look_'+os.path.basename(imagename))
	imagename='quick_look_'+os.path.basename(imagename)
	org_image=copy.deepcopy(imagename)
	header=imhead(imagename=imagename,mode='list')
	xcent=int(header['shape'][0]/2)
	ycent=int(header['shape'][1]/2)
	cell=np.rad2deg(abs(header['cdelt2'])) # In degree
	freq="{:.2f}".format(float(freq))
	xwidth=ywidth=int((field_of_view)/cell)
	box=str(xcent-int(xwidth/2))+','+str(ycent-int(ywidth/2))+','+str(xcent+int(xwidth/2))+','+str(ycent+int(ywidth/2))
	try:
		header=fits.getheader(imagename)
		if header['NAXIS']==4:
			if header['CTYPE3']=='STOKES':
				stokes_length=header['NAXIS3']
			elif header['CTYPE4']=='STOKES':
				stokes_length=header['NAXIS4']
		else:
			stokes_length=1
		importfits(fitsimage=imagename,imagename=imagename.split('.fits')[0]+'.image')
		imagename=imagename.split('.fits')[0]+'.image'
	except:
		header=imhead(imagename=imagename)
		stokes_axis=np.where(header['axisnames']=='Stokes')[0][0]
		stokes_length=header['shape'][stokes_axis]
	if stokes_length==1:
		stokes_list=['I']
	elif stokes_length==4:
		stokes_list=['I','Q','U','V']
	else:
		print ('Stokes axes are not I or IQUV.\n')
	fig = plt.figure(figsize=(8,8))
	plt.subplots_adjust(wspace=0.45, hspace=0.1)
	for i in range(len(stokes_list)):
		stokes=stokes_list[i]
		try:
			imsubimage(imagename=imagename,outfile='temp_'+stokes+'_'+os.path.basename(imagename)+'.image',box=box,stokes=stokes)
		except:
			return 
		exportfits(imagename='temp_'+stokes+'_'+os.path.basename(imagename)+'.image',fitsimage='temp_'+stokes+'_'+os.path.basename(imagename)+'.fits',dropdeg=True,dropstokes=True)
		data=fits.getdata('temp_'+stokes+'_'+os.path.basename(imagename)+'.fits')
		wlist=fits.getheader('temp_'+stokes+'_'+os.path.basename(imagename)+'.fits')
		w = wcs.WCS(wlist)
		if i==0:
			ax1 = fig.add_subplot(221, projection = w)
			im=ax1.imshow(data,cmap='seismic',origin='lower')
			fig.colorbar(im, ax=ax1)
			ax1.set_title('Stokes : '+stokes)
			ax1.set_xlabel('RA')
			ax1.set_ylabel('DEC')
		elif i==1:		
			ax2 = fig.add_subplot(222, projection = w)
			im=ax2.imshow(data,cmap='seismic',origin='lower')
			fig.colorbar(im, ax=ax2)
			ax2.set_title('Stokes : '+stokes)
			ax2.set_xlabel('RA')
			ax2.set_ylabel('DEC')
		elif i==2:
			ax3 = fig.add_subplot(223, projection = w)
			im=ax3.imshow(data,cmap='seismic',origin='lower')
			fig.colorbar(im, ax=ax3)
			ax3.set_title('Stokes : '+stokes)
			ax3.set_xlabel('RA')
			ax3.set_ylabel('DEC')
		elif i==3:
			ax4 = fig.add_subplot(224, projection = w)
			im=ax4.imshow(data,cmap='seismic',origin='lower')
			fig.colorbar(im, ax=ax4)
			ax4.set_title('Stokes : '+stokes)
			ax4.set_xlabel('RA')
			ax4.set_ylabel('DEC')
	title='Frequency : '+str(freq)+' MHz, Timestamp : '+str(timestamp)+' UTC'
	plt.suptitle(title,fontsize=12)	
	cwd=os.getcwd()
	outfile_dir=os.path.dirname(outfile)
	if outfile_dir=='':
		outfile=cwd+'/'+outfile
	plt.savefig(outfile)
	os.system('rm -rf casa*log temp_'+stokes+'_'+os.path.basename(imagename)+'* '+imagename+' '+org_image+' '+imagename.split('.fits')[0]+'.image')
	return outfile

while True:
	time.sleep(6000)
	if os.path.isdir(options.savedir+'/'+str(options.obsid)):
		os.chdir(options.savedir+'/'+str(options.obsid))
		break
while True:
	time.sleep(6000)
	imagelist=glob.glob('*.fits')
	if len(imagelist)>0:
		finished_list=glob.glob(str(options.basedir)+'/.Imaging_done_*')
		if len(finished_list)==int(options.num_ms):
			imagename=imagelist[-1]
			freqstr=imagename.split('freq_')[-1].split('_image')[0]
			datestrfile=imagename.split('time_')[-1].split('_freq')[0].split('_')
			datetimestr='/'.join(datestrfile[:3])+'/'+':'.join(datestrfile[3:])
			quickimage=get_quicklook_image(imagename,'sample_image_freq_'+freqstr+'_time_'+str('_'.join(datestrfile))+'.png',freqstr,datetimestr,field_of_view=2)
			msg_str='Dear PAIRCARS User,\n\nFinal imaging for : '+str(OBSID)+' is finished.\n Total number of images made : '+str(len(imagelist))+\
					'\n\nBest regards,\nPAIRCARS developing team'
			msg_subject='Notification from PAIRCARS : Final imaging : OBSID = '+str(OBSID)
			send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
			os.system('rm -rf '+quickimage)
			break
		else:
			imagename=imagelist[-1]
			freqstr=imagename.split('freq_')[-1].split('_image')[0]
			datestrfile=imagename.split('time_')[-1].split('_freq')[0].split('_')
			datetimestr='/'.join(datestrfile[:3])+'/'+':'.join(datestrfile[3:])
			quickimage=get_quicklook_image(imagename,'sample_image_freq_'+freqstr+'_time_'+str('_'.join(datestrfile))+'.png',freqstr,datetimestr,field_of_view=2)
			msg_str=''
			if len(finished_list)>0:
				for i in finished_list:
					msname=os.path.basename(i).split('.Imaging_done_')[-1].split('.ms')[0]+'.ms'
					msg_str+='Imaging finished for ms : '+msname+'.\n'
			if msg_str!='': 
				msg_str='Dear PAIRCARS User,\n\nFinal imaging for : '+str(OSBID)+' is running.\n Number of images made : '+str(len(imagelist))+\
					'.\n'+msg_str+'\n\nBest regards,\nPAIRCARS developing team'
			msg_subject='Notification from PAIRCARS : Final imaging : OBSID = '+str(OBSID)
			send_paircars_notification(inputs.email,msg_subject,msg_str,attachments=[quickimage])
			os.system('rm -rf '+quickimage)
	
