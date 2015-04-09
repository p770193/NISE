# -*- coding: utf-8 -*-
"""
Created on Thu Aug 7 2014

testing 2d frequency processing
@author: Dan
"""

from NISE.lib.misc.__init__ import NISE_path
#import matplotlib.pyplot as plt
import numpy as np
# important:  we can't call on scan directly for some reason; use trive to do it
import NISE.lib.measure as m
import NISE.experiments.trive as trive
import NISE.hamiltonians.H_ab as H_ab
import NISE.hamiltonians.params.inhom as inhom
reload(trive)
reload(H_ab)
reload(m)

autosave = True
imported = False
add_nrb  = False
nrb_level = 0.3 # amplitude, relative to sig max
slitwidth = 120.0 # width in wavenumbers; set to none if no slit

rel_path = r'data\2014.08.12 16-44-46'
if imported:
    filepath = r'//'.join([NISE_path, rel_path])
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
    trive.w1.points = np.linspace(6000, 8000, num=21)
    trive.w2.points = np.linspace(6000, 8000, num=21)
    trive.d2.points = np.linspace(-50, 200, num=6)
    w1 = trive.w1
    w2 = trive.w2
    d2 = trive.d2
    trive.exp.timestep = 4.0
    inhom_object = inhom.Inhom()#inhom_sampling='linear', num=5, sigma=100.0)
    out1 = trive.exp.scan(w1, w2, d2, H=H_ab.Omega(), inhom_object=inhom_object)
    out1.run(autosave=autosave)

if add_nrb:
    nrb1 = out1.nrb()
    nrb1_max = np.abs(nrb1).sum(axis=-1).max()
    out1_smax = np.abs(out1.sig).sum(axis=-1).max()
    print nrb1_max, out1_smax
    nrb1 *= out1_smax / nrb1_max * nrb_level
    print np.abs(nrb1).max(), np.abs(out1.sig).max()
    out1.sig = out1.sig.sum(axis=-2)[...,None,:]
    out1.sig += nrb1[...,None,:]


sig1 = m.Measure(out1, m.Mono, m.SLD)
# should be equivalent to the others if we set mono to pass
m.Mono.slitwidth=slitwidth
sig1.run()

#if len(sig1.pol.shape) > 2.:
print sig1.pol.shape
# run out1.axes() to figure out which axis is which
sig1.plot(0,yaxis=1)
