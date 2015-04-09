# -*- coding: utf-8 -*-
"""
Created on Fri Nov 07 00:59:03 2014

@author: Dan
"""

from NISE.lib.misc.__init__ import NISE_path
#import matplotlib.pyplot as plt
import numpy as np
# important:  we can't call on scan directly for some reason; use trive to do it
import NISE.lib.measure as m
import NISE.experiments.trive as trive
import NISE.hamiltonians.H_dk1 as H_dk
import NISE.hamiltonians.params.inhom as inhom
reload(trive)
reload(H_dk)
reload(m)

autosave  = True
imported  = False
add_nrb   = True
use_inhom = True
Da = 8
nrb_level = 0.3 # amplitude, relative to sig max
slitwidth = 120.0 # width in wavenumbers; set to none if no slit
# takes ~ 50 pts per second (with no inhomogeneity)

rel_path = r'\data\2015.02.05 02-40-47'
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
    trive.exp.set_coord(trive.d1, 0.)
    trive.exp.set_coord(trive.d2, 0.)
    trive.exp.set_coord(trive.ss, 40.)
    trive.w1.points = np.linspace(6000, 8800, num=21)
    trive.w2.points = np.linspace(6000, 8800, num=20)
    #trive.d2.points = np.linspace(-75, 225, num=13)
    trive.d2.points = np.array([-75, -25, 0, 25, 75, 150, 250, 500])
    #np.linspace(-75, 225, num=13)
    w1 = trive.w1
    w2 = trive.w2
    d2 = trive.d2
    trive.exp.timestep = 2.0
    if use_inhom:
        inhom_object = inhom.Inhom(inhom_sampling='linear', num=5, sigma=100.0)
    else:
        inhom_object = inhom.Inhom()
    H = H_dk.Omega(Da=Da)
    out1 = trive.exp.scan(w1, w2, d2, H=H, inhom_object=inhom_object)
    out1.run(autosave=autosave)

if add_nrb:
    nrb1 = out1.nrb()
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
sig1.plot(0,yaxis=1)
#sig1.plot(0,yaxis=2)
#sig1.plot(1,yaxis=2)
