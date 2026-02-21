# ---------------------------------------------------------
# Global warning & logging policy for PairCARS
# ---------------------------------------------------------

import warnings
import logging
import os

# -----------------------------
# Suppress Prefect / pydantic noise
# -----------------------------
warnings.filterwarnings(
    "ignore",
    message=".*pyproject_toml_table_header.*",
)

warnings.filterwarnings(
    "ignore",
    message=".*toml_file.*",
)

# -----------------------------
# Suppress tzlocal deprecated timezone warning
# -----------------------------
warnings.filterwarnings(
    "ignore",
    message=".*timezone is deprecated.*",
)

# -----------------------------
# Optional: suppress generic UserWarnings from libraries
# (comment out if you want stricter behavior)
# -----------------------------
warnings.filterwarnings("ignore", category=UserWarning)

# ---------------------------------------------------------
# Reduce noisy third-party logging
# ---------------------------------------------------------
logging.getLogger("distributed").setLevel(logging.ERROR)
logging.getLogger("distributed.scheduler").setLevel(logging.ERROR)
logging.getLogger("distributed.worker").setLevel(logging.ERROR)
logging.getLogger("tornado.application").setLevel(logging.CRITICAL)
logging.getLogger("tzlocal").setLevel(logging.ERROR)
logging.getLogger("prefect").setLevel(logging.ERROR)

# ---------------------------------------------------------
# PairCARS environment setup
# ---------------------------------------------------------
from .udocker_utils import set_udocker_env

set_udocker_env()

# ---------------------------------------------------------
# CASA log cleanup
# ---------------------------------------------------------
from casatasks import casalog

try:
    logfile = casalog.logfile()
    if os.path.exists(logfile):
        os.remove(logfile)
except Exception:
    pass

# ---------------------------------------------------------
# Astropy IERS configuration (offline safe)
# ---------------------------------------------------------
from astropy.utils import iers

iers.conf.auto_download = False
iers.conf.auto_max_age = None
