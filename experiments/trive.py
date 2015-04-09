# -*- coding: utf-8 -*-
"""
Created on Sun Jun 17 14:35:52 2012
@author: Dan

module for TRIVE:  functions that assume:
    three pulses
    w2=w3
    w_rw = w1 default
    
so:  free vars to consider:  
    d1, d2, (d3)
    w1, w2
    
"""

from NISE.lib.misc import *
import NISE.lib.scan as S
#import matplotlib.pyplot as plt
#from NISE.lib.axis import Axis

reload(S)

w1 = S.Axis(0, 'w', name=r'$\mathsf{\omega_1}$')
w2 = S.Axis(1, 'w', also=[2], name=r'$\mathsf{\omega_2 = \omega_{2^\prime}}$')
# write d2 in new delay coordinate system--all relative to w2
d1 = S.Axis(2, 'd', name=r'$\mathsf{\tau_{2^\prime 2}}$')
d2 = S.Axis(0, 'd', name=r'$\mathsf{\tau_{12}}$')
ws = S.Axis(0, 'w', also=[1,2], name=r'$\mathsf{\omega_1 = \omega_2 = \omega_{2^\prime}}$')
A1 = S.Axis(0, 'A', name=r'$k_1$ Fluence (a.u.)')
A2 = S.Axis(1, 'A', also=[2], name=r'$k_2, k_{2^\prime}$ Fluence (a.u.)')
As = S.Axis(0, 'A', also=[1,2], name=r'$k_1, k_2, k_{2^\prime}$ Fluence (a.u.)')
# cycle phases
p1 = S.Axis(0, 'p')
p2 = S.Axis(1, 'p')
p3 = S.Axis(2, 'p')
s1 = S.Axis(0, 's', name=r'$\mathsf{FWHM (\vec{k}_1, fs)}$')
s2 = S.Axis(1, 's', also=[2], name=r'$\mathsf{FWHM (\vec{k}_2, \vec{k}_{2^\prime}, fs)}$')
ss = S.Axis(0, 's', also=[1,2], name=r'$\mathsf{FWHM (\vec{k}_1, \vec{k}_2, \vec{k}_{2^\prime}, fs)}$')

# defaults generate TRIVE conditions
exp = S.Experiment(pm=[1,-1,1])
