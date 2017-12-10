# -*- coding: utf-8 -*-
"""
Created on Sun Jun 17 14:35:52 2012
@author: Dan

module for TRSF:  
    three pulses
    
free vars to consider:  
    d1 = tau_13
    d2 = tau_23
    w1, w2, w3 (w4 is an output variable, not input)
    
"""

from __future__ import absolute_import, division, print_function, unicode_literals

from .. import lib
from ..lib.misc import *
from ..lib import scan as S
#import matplotlib.pyplot as plt
#from NISE.lib.axis import Axis

w1 = S.Axis(0, 'w', name=r'$\mathsf{\omega_1}$')
w2 = S.Axis(1, 'w', name=r'$\mathsf{\omega_2}$')
w3 = S.Axis(2, 'w', name=r'$\mathsf{\omega_3}$')
wIR = S.Axis(0, 'w', also=[1], name=r'$\mathsf{\omega_1 = \omega_2}$')

d1 = S.Axis(0, 'd', name=r'$\mathsf{\tau_{13}}$')
d2 = S.Axis(1, 'd', name=r'$\mathsf{\tau_{23}}$')

#A1 = S.Axis(0, 'A', name=r'$E_1$ Fluence (a.u.)')
#A2 = S.Axis(1, 'A', name=r'$E_2$ Fluence (a.u.)')
#A3 = S.Axis(2, 'A', name=r'$E_2$ Fluence (a.u.)')
#As = S.Axis(0, 'A', also=[1,2], name=r'$E_1, E_2, E_3$ Fluence (a.u.)')

# cycle phases
#p1 = S.Axis(0, 'p')
#p2 = S.Axis(1, 'p')

#s1 = S.Axis(0, 's', name=r'$\mathsf{FWHM (\vec{k}_1, fs)}$')
#s2 = S.Axis(1, 's', name=r'$\mathsf{FWHM (\vec{k}_2, fs)}$')
#s3 = S.Axis(2, 's', name=r'$\mathsf{FWHM (\vec{k}_2, fs)}$')
ss = S.Axis(0, 's', also=[1,2], name=r'$\mathsf{FWHM (E_1, E_2, E_3, fs)}$')

# communicates the appropriate rotating wave
exp = S.Experiment(pm=[1,1,1])
#exp.early_buffer = 1500.
#exp.late_buffer = 3000.
