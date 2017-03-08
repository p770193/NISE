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

from __future__ import absolute_import, division, print_function, unicode_literals

from .. import lib
from ..lib.misc import *
from ..lib import scan as S

pulse_class_name = 'Gauss_rwa'

# introduce these all to the namespace
# k_1 - k_2 + k_3
def assert_pulse_class(pulse_class_name):
    global w1 
    w1 = S.Axis(0, 'w', name=r'$\mathsf{\bar\nu_1}$', 
                units = 'wn',
                pulse_class_name = pulse_class_name)
    global w2 
    w2 = S.Axis(1, 'w', also=[2], 
                name=r'$\mathsf{\bar\nu_2 = \bar\nu_{2^\prime}}$',
                units = 'wn',
                pulse_class_name = pulse_class_name)
    global d1 # tau_2'2
    d1 = S.Axis(2, 'd', name=r'$\mathsf{\tau_{2^\prime 2}}$',
                units = 'fs',
                pulse_class_name = pulse_class_name)
    # write d2 in new delay coordinate system--all relative to E_2
    global d2 # tau_12
    d2 = S.Axis(0, 'd', name=r'$\mathsf{\tau_{12}}$',
                units = 'fs',
                pulse_class_name = pulse_class_name)
    global ws
    ws = S.Axis(0, 'w', also=[1,2], 
                name=r'$\mathsf{\omega_1 = \omega_2 = \omega_{2^\prime}}$',
                units = 'wn',
                pulse_class_name = pulse_class_name)
    global A1 
    A1 = S.Axis(0, 'A', name=r'$k_1$ Fluence (a.u.)',
                pulse_class_name = pulse_class_name)
    global A2 
    A2 = S.Axis(1, 'A', also=[2], name=r'$k_2, k_{2^\prime}$ Fluence (a.u.)',
                pulse_class_name = pulse_class_name)
    global As 
    As = S.Axis(0, 'A', also=[1,2], 
                name=r'$k_1, k_2, k_{2^\prime}$ Fluence (a.u.)',
                pulse_class_name = pulse_class_name)
    # cycle phases
    global p1 
    p1 = S.Axis(0, 'p', pulse_class_name = pulse_class_name)
    global p2 
    p2 = S.Axis(1, 'p', pulse_class_name = pulse_class_name)
    global p3 
    p3 = S.Axis(2, 'p', pulse_class_name = pulse_class_name)
    
    # pulse width variations
    global s1 
    s1 = S.Axis(0, 's', name=r'$\mathsf{FWHM (\vec{k}_1, fs)}$',
                pulse_class_name = pulse_class_name)
    global s2
    s2 = S.Axis(1, 's', also=[2], 
                name=r'$\mathsf{FWHM (\vec{k}_2, \vec{k}_{2^\prime}, fs)}$',
                pulse_class_name = pulse_class_name)
    global ss
    ss = S.Axis(0, 's', also=[1,2], 
                name=r'$\mathsf{FWHM (\vec{k}_1, \vec{k}_2, \vec{k}_{2^\prime}, fs)}$',
                pulse_class_name = pulse_class_name)
    
    # defaults generate TRIVE conditions
    global exp 
    exp = S.Experiment(pm=[1,-1,1], pulse_class_name=pulse_class_name)

assert_pulse_class(pulse_class_name)

