# whenever someone purely imports dppa, core main functions are also imported.
from .core import run, export, set_debug_mode
# basic config for logging
import logging
logging.basicConfig(level=logging.INFO, format='%(message)s')