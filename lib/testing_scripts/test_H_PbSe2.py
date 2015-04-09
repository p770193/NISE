# -*- coding: utf-8 -*-
"""
Created on Mon Feb 09 12:40:52 2015

Updated simulation with different parameters

@author: Dan
"""
from NISE.lib.misc.__init__ import NISE_path
#import matplotlib.pyplot as plt
import numpy as np
# important:  we can't call on scan directly for some reason; use trive to do it
import NISE.lib.measure as m
import NISE.experiments.trive as trive
import NISE.hamiltonians.H_PbSe2 as H_1
import NISE.hamiltonians.params.inhom as inhom
reload(trive)
reload(H_1)
reload(m)
1/0.
autosave  = True
imported  = True
add_nrb   = False
use_inhom = True
plot      = False
nrb_level = 0.1 # amplitude, relative to sig max
slitwidth = 120.0 # width in wavenumbers; set to none if no slit
# takes ~ 50 pts per second (with no inhomogeneity)

rel_path = r'\data\2015.03.18 00-53-44'
if imported:
    filepath = r''.join([NISE_path, rel_path])
    out1 = trive.S.Scan._import(filepath)
    # because I didn't set up the class correctly, coords_set was not written
    # on older files; i need to manually re-enter it
    try: out1.coords_set
    except AttributeError:
        print 'did not find coords_set'
        out1.coords_set = [[0, 3], [2, 3], [1, 3], [0, 2]]
        out1.save()
    w1 = out1.axis_objs[0]
    w2 = out1.axis_objs[1]
else:
    trive.exp.set_coord(trive.d1, 25.)
    trive.exp.set_coord(trive.d2, 0.)
    trive.exp.set_coord(trive.ss, 45.)
    trive.w1.points = np.linspace(6000, 9000, num=21)
    trive.w2.points = np.linspace(6000, 9000, num=20)
    trive.ws.points = np.linspace(6000, 8800, num=21)
    trive.d2.points = np.linspace(-125, 125, num=11)
    #trive.d2.points = np.array([-75, -25, 0, 25, 75])#, 150, 250, 500])
    #np.linspace(-75, 225, num=13)
    w1 = trive.w1
    w2 = trive.w2
    ws = trive.ws
    d2 = trive.d2
    trive.exp.timestep = 2.0
    if use_inhom:
        inhom_object = inhom.Inhom(inhom_sampling='linear', num=5, sigma=75.0)
    else:
        inhom_object = inhom.Inhom()
    H = H_1.Omega()
    H.out_group = [[13,14,15,16],[12,17]]
    #H.Gamma_bb  = 0. # eliminate coupling to examine cross-peak
    out1 = trive.exp.scan(w1, w2, d2, H=H, inhom_object=inhom_object)
    out1.run(autosave=autosave)
"""
outsig = out1.sig.copy()
print out1.sig.shape
out1.sig = outsig[...,1:2,:]
print out1.sig.shape
#"""
if add_nrb:
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
if plot:
    # run out1.axes() if you forgot which axis is which
    sig1.plot(0,yaxis=1, zoom=2)
    #sig1.plot(0,yaxis=2, zoom=2)
    #sig1.plot(1,yaxis=2)
