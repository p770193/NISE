from . import trive
from . import tsf
from . import TA

import os
import glob

# include all written experiments in the from import
modules = glob.glob(os.path.dirname(__file__)+"/*.py")
__all__ = [ os.path.basename(f)[:-3] for f in modules]

