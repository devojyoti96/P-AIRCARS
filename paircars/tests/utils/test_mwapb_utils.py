import traceback
import os
from casatasks import casalog
from paircars.utils.mwapb_utils import *

try:
    casalogfile = casalog.logfile()
    os.system("rm -rf " + casalogfile)
except BaseException:
    traceback.print_exc()
    pass
