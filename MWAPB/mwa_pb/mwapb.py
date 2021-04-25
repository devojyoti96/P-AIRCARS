import numpy as np,sys,struct,time
from . import mwa_sweet_spots
from . import config
from . import primary_beam
from . import beam_tools
from paircars.access_ms import *
from paircars.basic_func import *
from casatools import *
from casatasks import *
from astropy.io import fits
from CALIBRATE.access_calibrate import *
from numpy.linalg import inv
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

	def calc_beamjones_phasecenter(self,outputfile=''):
		'''
		Function to calculate the Beam Jones matrix for the phasecenter.
		Parameters:
		outputfile = Name of the file to save the Jones matrix. If blank the Beam jones will not be writen to disk.
		Return :
		Return Jones matrix in IAU coordinate and parallactic angle corrected
		'''
		AM=AccessMS(self.msname)
		radecstr,ra,dec=AM.get_phasecenter()
		LAT,LON,ALT=AM.get_observatory_loc()
		parang=AM.get_phasecenter_parang(source_field=0,combine='field')
		alt,az=AM.get_altaz(source_field=0,source_scan=1)
		metaheader=fits.getheader(self.metafits)
		gridpoint=metaheader['GRIDNUM']
		delays=[int(h) for h in metaheader['DELAYS'].split(',')]
		delays = np.vstack((delays, delays))
		freq=AM.calc_meanfreq()
		Jones_FullEE=primary_beam.MWA_Tile_full_EE(np.deg2rad(90-np.rad2deg(alt)),az,freq,delays=delays,zenithnorm=True,jones=True,interp=True)
		coord_transform=np.matrix([[0.0,-1.0],[-1.0,0.0]])
		Jones_FEE_2D_IAU=np.matmul(coord_transform,np.matrix(Jones_FullEE))
		pa_matrix=np.matrix([[np.cos(parang),np.sin(parang)],[-np.sin(parang),np.cos(parang)]])
		Jones_FEE_2D_IAU=np.matrix(np.matmul(pa_matrix,Jones_FEE_2D_IAU))
		if self.invbeam==True:
			Jones_FEE_2D_IAU=inv(Jones_FEE_2D_IAU)
		if outputfile!='':
			np.save(outputfile,Jones_FEE_2D_IAU)
		return Jones_FEE_2D_IAU


	def MWA_phasecenter_beam_jones(self,outputfile='',nant=128,nchan=1024,nint=1000):
		'''
		Function to save MWA beam Jones for phasecenter as CALIBRATE caltable format
		Parameters:
		outputfile= = Name of the output file
		nant= Number of total antennas (Including the flag once)
		nchan = Number of frequency channels
		nint = Number of time slices
		Return:
		Beam jones of the phasecenter as CALIBRATE caltable format, Beam JOnes matrix
		'''
		if outputfile=='':
			return
		beam_jones=self.calc_beamjones_phasecenter(outputfile='')
		beam_jones=inv(np.matrix(beam_jones))
		header = struct.pack("8s",b"MWAOCAL")+struct.pack("i",0)+struct.pack("i",0)+struct.pack("i",int(nint))+struct.pack("i",int(nant))+struct.pack("i",int(nchan))+\
				struct.pack("i",4)+struct.pack("d",0.0)+struct.pack("d",0.0)
		numpy_data=np.empty([nint,nant,nchan,8])
		data=np.empty([2,2,nint,nant,nchan],dtype=np.complex128)
		for i in range(nint):
			for j in range(nant):
				for k in range(nchan):
					data[:,:,i,j,k]=np.array(beam_jones)
					data_re=np.real(data[:,:,i,j,k].flatten())
					data_im=np.imag(data[:,:,i,j,k].flatten())
					numpy_data[i,j,k,:]=np.insert(data_im, np.arange(len(data_re)),data_re)
		numpy_data=numpy_data.flatten()
		fil=open(outputfile,'wb')
		fil.write(header)
		fil.close()
		with open(outputfile,mode='ba+') as f:
			numpy_data.tofile(f,format='np.float64')
		return outputfile,beam_jones





