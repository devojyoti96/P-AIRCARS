import os
os.environ["PYTHONWARNINGS"] = "ignore"
import logging
logging.getLogger("distributed").setLevel(logging.ERROR)
logging.getLogger("tornado.application").setLevel(logging.CRITICAL)
logging.getLogger("tzlocal").setLevel(logging.ERROR)
from .udocker_utils import set_udocker_env
set_udocker_env()
from casatasks import casalog
try:
    logfile = casalog.logfile()
    os.remove(logfile)
except BaseException:
    pass
from astropy.utils import iers
iers.conf.auto_download = False
iers.conf.auto_max_age = None
