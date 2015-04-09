# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 15:22:27 2015

test script for running a tsf demo

@author: Dan
"""

from NISE.lib.misc.__init__ import NISE_path
#import matplotlib.pyplot as plt
import numpy as np
# important:  we can't call on scan directly for some reason; use trive to do it
import NISE.lib.measure as m
import NISE.experiments.tsf as tsf
import NISE.hamiltonians.H_tsf1 as H_1
import NISE.hamiltonians.params.inhom as inhom
reload(tsf)
reload(H_1)
reload(m)

d1 = 1000.
autosave  = True
imported  = True
add_nrb   = False
use_inhom = False
nrb_level = 0.3 # amplitude, relative to sig max
slitwidth = 12.0 # width in wavenumbers; set to none if no slit
# takes ~ 50 pts per second (with no inhomogeneity)

rel_path = r'\data\2015.02.11 16-37-05 - TSF1 test'
if imported:
    filepath = r''.join([NISE_path, rel_path])
    out1 = tsf.S.Scan._import(filepath)
    w1 = out1.axis_objs[0]
    w2 = out1.axis_objs[1]
else:
    tsf.exp.set_coord(tsf.d1, d1)
    # set pulse widths to be ~1000 fs FWHM
    tsf.exp.set_coord(tsf.ss, 1000.)
    tsf.w1.points = np.linspace(1700, 2100, num=31)
    tsf.w2.points = np.linspace(1700, 2100, num=30)
    w1 = tsf.w1
    w2 = tsf.w2
    tsf.exp.timestep = 5.0 # fs
    if use_inhom:
        inhom_object = inhom.Inhom(inhom_sampling='linear', num=5, sigma=100.0)
    else:
        inhom_object = inhom.Inhom()
    H = H_1.Omega()
    out1 = tsf.exp.scan(w1, w2, H=H, inhom_object=inhom_object)
    out1.run(autosave=autosave)

if add_nrb:
    """
    need to rewrite this part?  maybe not rigorously correct...
    """
    nrb1 = out1.nrb()
    print nrb1.shape
    print out1.sig.shape
    nrb1_max = np.abs(nrb1).sum(axis=-1).max()
    out1_smax = np.abs(out1.sig).sum(axis=-1).max()
    #print nrb1_max, out1_smax
    nrb1 *= out1_smax / nrb1_max * nrb_level
    #print np.abs(nrb1).max(), np.abs(out1.sig).max()
    out1.sig = out1.sig.sum(axis=-2)[...,None,:]
    out1.sig += nrb1[...,None,:]

sig1 = m.Measure(out1, m.Mono, m.SLD)
# should be equivalent to the others if we set mono to pass
m.Mono.slitwidth=slitwidth
sig1.run()

print sig1.pol.shape
# run out1.axes() if you forgot which axis is which
sig1.plot(0,yaxis=1, zoom=3)
#sig1.plot(0,yaxis=2)
#sig1.plot(1,yaxis=2)
