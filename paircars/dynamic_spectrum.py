import os
from casatools import msmetadata,table,measures,quanta,agentflagger,image,calibrater,ms
from casatasks import *
os.system('rm -rf casa*log')

class DynamicSpectrum:
	def __init__(self,msname):
		self.msname=msname
		self.ms=ms()

	def cal_norm_crosscorr(self):
		'''
		Function to obtain normalised cross correlation amplitude
		'''
		self.ms.open(self.msname)
		

