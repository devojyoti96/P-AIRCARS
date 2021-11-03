import numpy as np,os,sys
from paircars.basic_func import *
from optparse import OptionParser


if __name__=='__main__':
	usage= ' Compress lists of caltables in P-AIRCARS format'
	parser = OptionParser(usage=usage)
	parser.add_option('--caltables',dest="caltables",default=None,help="List of caltables",metavar="Comma separated string")
	parser.add_option('--compressed_file',dest="final",default=None,help="Compressed final name",metavar="String")
	(options, args) = parser.parse_args()

caltable_list=str(options.caltables).split(',')
compressed_file=compress_files(caltable_list,str(options.final))
print ('Compressed caltable saved as : '+compressed_file+'\n')
