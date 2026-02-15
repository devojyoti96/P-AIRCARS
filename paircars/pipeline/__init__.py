import os

os.environ["PYTHONWARNINGS"] = "ignore"
import logging

logging.getLogger("distributed").setLevel(logging.ERROR)
logging.getLogger("tornado.application").setLevel(logging.CRITICAL)
logging.getLogger("tzlocal").setLevel(logging.ERROR)
from casatasks import casalog
from astropy.utils import iers

try:
    logfile = casalog.logfile()
    os.remove(logfile)
except BaseException:
    pass

iers.conf.auto_download = False
iers.conf.auto_max_age = None
