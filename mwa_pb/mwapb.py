import numpy as np,sys,struct,time
from . import mwa_sweet_spots
from . import config
from . import primary_beam
from . import beam_tools
from paircars.access_ms import *
from paircars.basic_func import *
from paircars.flagger import *
from casatools import msmetadata,table,measures,quanta,agentflagger,image,calibrater,ms
from casatasks import *
from astropy.io import fits
from CALIBRATE.access_calibrate import *
from numpy.linalg import inv
from astropy.coordinates import SkyCoord, FK5
import astropy.wcs as pywcs
from astropy.time import Time
os.system('rm -rf casa*log')
'''
Code is written by Devojyoti Kansabanik, 23 Feb, 2021
'''

class MWA_PrimaryBeam:
	'''
	Generic class for correcting ideal beam response of MWA
	Parameters:
	msname = Name of the measurement set
	metafits = Name of the observation metafits
	inverse_beam = False, whether want beam jones or inverse of the beam jones
	'''
	def __init__(self,msname,metafits,inverse_beam=False):
		self.msname=msname
		self.metafits=metafits
		self.invbeam=inverse_beam

	def calc_beamjones_phasecenter(self,outputfile='',skip_freq=1.28):
		'''
		Function to calculate the Beam Jones matrix for the phasecenter.
		Parameters:
		outputfile = Name of the file to save the Jones matrix. If blank the Beam jones will not be writen to disk.
		skip_freq = Frequency interval in MHz to make independent beams (default : 1.28 MHz). If anything greater than 1.28 MHz is given it will be overwritten to 1.28 MHz
		Return :
		Return Jones matrix in IAU coordinate and parallactic angle corrected as a function of frequencies in an array
		'''
		AM=AccessMS(self.msname)
		freqs=AM.get_freqs()/10**6
		if skip_freq>1.28:
			print ('Frequency interval is chosen more than a coarse channel. Beam can change significantly, thus setting frequency interval to 1.28 MHz.\n')
			skip_freq=1.28
		coarse_freqs=[]
		for i in freqs:
			try:
				coarse_freqs.append(freq_to_MWA_coarse(i)*float(skip_freq))
			except:
				pass
		coarse_freqs=np.array(coarse_freqs)
		freqlist=[str(i) for i in list(np.unique(coarse_freqs))]
		print ('Calculating beam jones for frequencies : '+','.join(freqlist)+' MHz\n')
		radecstr,ra,dec=AM.get_phasecenter()
		LAT,LON,ALT=AM.get_observatory_loc()
		parang=np.deg2rad(AM.get_phasecenter_parang(source_field=0,combine='field')-360.0)
		alt,az=AM.get_altaz(source_field=0,source_scan=1)
		metaheader=fits.getheader(self.metafits)
		gridpoint=metaheader['GRIDNUM']
		delays=[int(h) for h in metaheader['DELAYS'].split(',')]
		delays = np.vstack((delays, delays))
		beam_jones=[]
		prefreq=0
		for freq in freqs:
			nearest_coarse_freq=coarse_freqs[np.argmin(abs(coarse_freqs-freq))]
			if nearest_coarse_freq!=prefreq:
				Jones_FullEE=primary_beam.MWA_Tile_full_EE(np.deg2rad(90-np.rad2deg(alt)),az,nearest_coarse_freq*10**6,delays=delays,zenithnorm=True,jones=True,interp=True)
				coord_transform=np.matrix([[0.0,-1.0],[-1.0,0.0]])
				Jones_FEE_2D_IAU=np.matmul(coord_transform,np.matrix(Jones_FullEE))
				pa_matrix=np.matrix([[np.cos(parang),np.sin(parang)],[-np.sin(parang),np.cos(parang)]])
				Jones_FEE_2D_IAU_parang=np.matrix(np.matmul(pa_matrix,Jones_FEE_2D_IAU))
				if self.invbeam==True:
					Jones_FEE_2D_IAU_parang=inv(Jones_FEE_2D_IAU_parang)
				prefreq=nearest_coarse_freq
			else:
				Jones_FEE_2D_IAU_parang=beam_jones[-1]
			beam_jones.append(Jones_FEE_2D_IAU_parang)
		if outputfile!='':
			np.save(outputfile,np.array(beam_jones))
		return np.array(beam_jones)


	def MWA_phasecenter_beam_jones(self,outputfile='',skip_freq=1.28):
		'''
		Function to save MWA beam Jones for phasecenter as CALIBRATE caltable format
		Parameters:
		outputfile= = Name of the output file
		skip_freq = Frequency interval in MHz to make independent beams (default : 1.28 MHz). If anything greater than 1.28 MHz is given it will be overwritten to 1.28 MHz
		Return:
		Beam jones of the phasecenter as CALIBRATE caltable format, Beam Jones matrix
		'''
		if outputfile=='':
			return
		beam_jones=self.calc_beamjones_phasecenter(outputfile='',skip_freq=float(skip_freq))
		AM=AccessMS(self.msname)
		nant=AM.get_num_antenna()
		nchan=AM.get_num_channels()
		nint=AM.get_num_timestamps()
		freqs=AM.get_freqs()/10**6
		start_freq=freqs[0]
		end_freq=freqs[-1]
		mjdsecs=AM.get_timestamps_in_mjdsecs()[0]
		startmjd=mjdsecs[0]
		endmjd=mjdsecs[-1]
		header = struct.pack("8s",b"MWAOCAL")+struct.pack("i",0)+struct.pack("i",0)+struct.pack("i",int(nint))+struct.pack("i",int(nant))+struct.pack("i",int(nchan))+\
				struct.pack("i",4)+struct.pack("d",0.0)+struct.pack("d",0.0)
		beam_array=np.array([np.array([np.real(inv(a)).flatten(),np.imag(inv(a)).flatten()]).flatten(order='F') for a in beam_jones])
		beam_array=beam_array.reshape((1,1,beam_array.shape[0],beam_array.shape[-1]))
		numpy_data=np.repeat(np.repeat(beam_array,nint,axis=0),nant,axis=1)
		numpy_data=numpy_data.flatten(order='C')
		fil=open(outputfile,'wb')
		fil.write(header)
		fil.close()
		with open(outputfile,mode='ba+') as f:
			numpy_data.tofile(f,format='np.float64')
		bin_data=np.fromfile(outputfile,dtype=np.float64)
		os.system('rm -rf '+outputfile)
		np.save(outputfile,np.array([bin_data,start_freq,end_freq,startmjd,endmjd,nchan,nint],dtype='object'))
		os.system('mv '+outputfile+'.npy '+outputfile)
		return outputfile,beam_jones

