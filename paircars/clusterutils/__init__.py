import os

os.environ["PYTHONWARNINGS"] = "ignore"
import logging

logging.getLogger("distributed").setLevel(logging.ERROR)
logging.getLogger("tornado.application").setLevel(logging.CRITICAL)
logging.getLogger("tzlocal").setLevel(logging.ERROR)

