# -*- coding: utf-8 -*-
"""
Created on Sun Jun 17 14:35:52 2012
@author: Dan

module for the IR part of TSF:  functions that assume:
    two pulses
    w_rw = w1 default
    
so:  free vars to consider:  
    d1, d2, (d3)
    w1, w2
    
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import importlib

from .. import lib
from ..lib.misc import *
from ..lib import scan as S
#import matplotlib.pyplot as plt
#from NISE.lib.axis import Axis

importlib.reload(S)

w1 = S.Axis(0, 'w', name=r'$\mathsf{\omega_1}$')
w2 = S.Axis(1, 'w', name=r'$\mathsf{\omega_2}$')
d1 = S.Axis(1, 'd', name=r'$\mathsf{\tau_{21}}$')
ws = S.Axis(0, 'w', also=[1], name=r'$\mathsf{\omega_1 = \omega_2}$')
A1 = S.Axis(0, 'A', name=r'$k_1$ Fluence (a.u.)')
A2 = S.Axis(1, 'A', name=r'$k_2$ Fluence (a.u.)')
As = S.Axis(0, 'A', also=[1], name=r'$k_1, k_2$ Fluence (a.u.)')
# cycle phases
p1 = S.Axis(0, 'p')
p2 = S.Axis(1, 'p')
s1 = S.Axis(0, 's', name=r'$\mathsf{FWHM (\vec{k}_1, fs)}$')
s2 = S.Axis(1, 's', name=r'$\mathsf{FWHM (\vec{k}_2, fs)}$')
ss = S.Axis(0, 's', also=[1], name=r'$\mathsf{FWHM (\vec{k}_1, \vec{k}_2, fs)}$')

# communicates the appropriate rotating wave
exp = S.Experiment(pm=[1,1])
exp.early_buffer = 1500.
exp.late_buffer = 3000.
