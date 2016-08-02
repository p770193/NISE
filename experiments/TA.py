# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 12:59:14 2015

transient absorption common degrees of freedom

one delay

@author: Dan
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import importlib

from .. import lib
from ..lib.misc import *
from ..lib import scan as S
#import matplotlib.pyplot as plt
#from NISE.lib.axis import Axis

importlib.reload(S)

w1 = S.Axis(0, 'w', name=r'$\mathsf{\omega_1 (cm^{-1})}$')
w2 = S.Axis(1, 'w', also=[2], name=r'$\mathsf{\omega_p (cm^{-1})}$')
# write d2 in new delay coordinate system--all relative to w2
T = S.Axis(1, 'd', name=r'$\mathsf{\tau_{p} (fs)}$')
ws = S.Axis(0, 'w', also=[1,2], name=r'$\mathsf{\omega_1 = \omega_p (cm^{-1})}$')
A1 = S.Axis(0, 'A', name=r'$k_1$ Fluence (a.u.)')
A2 = S.Axis(1, 'A', also=[2], name=r'$k_p$ Fluence (a.u.)')
# cycle phases
p1 = S.Axis(0, 'p')
p2 = S.Axis(1, 'p')
p3 = S.Axis(2, 'p')
s1 = S.Axis(0, 's', name=r'$\mathsf{FWHM (\vec{k}_1, fs)}$')
s2 = S.Axis(1, 's', also=[2], name=r'$\mathsf{FWHM (\vec{k}_p, fs)}$')
ss = S.Axis(0, 's', also=[1,2], name=r'$\mathsf{FWHM (\vec{k}_1, \vec{k}_p, fs)}$')

# defaults generate TRIVE conditions
exp = S.Experiment(pm=[1,0])
