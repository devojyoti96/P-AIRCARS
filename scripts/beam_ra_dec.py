import numpy as np,sys,struct,time,os
from mwa_pb import mwa_sweet_spots
from mwa_pb import config
from mwa_pb import primary_beam
from mwa_pb import beam_tools
from mwa_pb.calc_beam import calc_beamjones_info
from casatools import msmetadata,table,measures,quanta,agentflagger,image,calibrater,ms
from casatasks import *
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import EarthLocation,SkyCoord,AltAz
from astropy import units as u
from numpy.linalg import inv
from optparse import OptionParser
os.system('rm -rf casa*log')
'''
Code is written by Devojyoti Kansabanik, 24 April, 2021
'''

if __name__=='__main__':
	usage= 'Calculate beam Jones and Stokes beam value at a given ra, dec and epoch for a list of frequencies. Stokes I beam is defined as I_PB=(XX^2+YY^2)/2'
	parser = OptionParser(usage=usage)
	parser.add_option('--ra',dest="ra",default=0,help="RA in degree",metavar="Float")
	parser.add_option('--dec',dest="dec",default=0,help="DEC in degree",metavar="Float")
	parser.add_option('--epoch',dest="epoch",default='2014-05-10 00:00:00',help="Epoch in yyyy-mm-dd hh:mm:ss format",metavar="String")
	parser.add_option('--freq',dest="freq",default='150.0',help="Frequency in MHz",metavar="Comma separated string")
	parser.add_option('--metafits',dest="metafits",default=None,help="Name of the metafits file",metavar="Metafits file")
	parser.add_option('--gridpoint',dest="gridpoint",default=0,help="Gridpoint number. If metafits is given gridpoint parameter will not be used",metavar="Integer")
	parser.add_option('--convention',dest="convention",default='MWA',help="Coordinate convention (IAU or MWA)",metavar="String")
	(options, args) = parser.parse_args()

	if os.path.isfile(str(options.metafits))==False or options.metafits==None:
		print ('Metafits file is not given.\n')
		metafits=''
	else:
		metafits=str(options.metafits)
	if str(options.convention)!='MWA' and str(options.convention)!='IAU':
		convention='MWA'
	else:
		convention=str(options.convention)

	freqlist=freq.aplit(',')

	beam_jones,I_beam=calc_beamjones_info(ra=float(options.ra),dec=float(options.dec),obstime=str(options.epoch),gridpoint=int(options.gridpoint),freq=freqlist,\
							metafits=options.metafits,coord=str(options.convention))
	print ('Frequency list : ',freqlist)
	print ('Beam Jones : ',beam_jones)
	print ('Stokes I beam value : ',I_beam)


