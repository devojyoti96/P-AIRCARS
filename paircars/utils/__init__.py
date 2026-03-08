import os
os.environ["PYTHONWARNINGS"] = "ignore"
import logging
from astropy.utils import iers
from casatasks import casalog
from .udocker_utils import set_udocker_env
logging.getLogger("distributed").setLevel(logging.ERROR)
logging.getLogger("tornado.application").setLevel(logging.CRITICAL)
logging.getLogger("tzlocal").setLevel(logging.ERROR)
set_udocker_env()
try:
    logfile = casalog.logfile()
    os.remove(logfile)
except BaseException:
    pass
iers.conf.auto_download = False
iers.conf.auto_max_age = None
