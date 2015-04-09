# misc module contains effective globals
import numpy as np
from os.path import realpath
folders = realpath(__file__).split('\\')
NISE_path = '\\'.join(folders[:-3])

"""
from os.path import realpath
import sys
folders = realpath(__file__).split('\\')
NISE_path = '\\'.join(folders[:-2])
# make this path the priority when searching for files
if path not in sys.path:
    sys.path.insert(0, path)
    # for just adding to the python path var
    print 'added {0} to python path'.format(path)
"""
timestep = 4.0
# buffers for integration
early_buffer = 100.0
late_buffer  = 400.0

wn_to_omega = 2*np.pi*3*10**-5  #omega is radians / fs

# set the current working directory relative to this

def rotor(theta):
    """
    returns 1 rotated CW about the complex plane by angle theta i.e. e^(-i*theta)
    can work on scalars as well as higher rank objects
    note:  can be replaced by np.exp(-j*theta)
    """
    #return complex(np.cos(theta), -np.sin(theta))
    return np.cos(theta)*complex(1,0) - np.sin(theta)*complex(0,1)
    
import os, errno
def mkdir_p(path):
    try:
        os.mkdir(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise
