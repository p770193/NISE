# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 15:27:51 2015

used for testing TOs of H0 hamiltonian

@author: Dan
"""
#if __name__ == '__main__':

from NISE.lib.misc.__init__ import NISE_path
#import matplotlib.pyplot as plt
import numpy as np
# important:  we can't call on scan directly for some reason; use trive to do it
import NISE.lib.measure as m
import NISE.experiments.trive as trive
import NISE.hamiltonians.H0 as H_1
import NISE.hamiltonians.params.inhom as inhom
reload(trive)
reload(H_1)
reload(m)

autosave  = True
imported  = False
add_nrb   = False
use_inhom = False
plot      = True
mp        = True
gh_quad   = False
TOs       = np.arange(6) + 1
nrb_level = 1. # amplitude, relative to sig max
slitwidth = 120.0 # width in wavenumbers; set to none if no slit
# takes ~ 50 pts per second (with no inhomogeneity)
if __name__ == '__main__':
    rel_path = r'\data\2015.03.28 18-08-31'
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
        trive.exp.set_coord(trive.ws, 7000.)
        trive.exp.set_coord(trive.ss, 45.)
        #trive.d2.points = np.array([-75, -25, 0, 25, 75])#, 150, 250, 500])
        #np.linspace(-75, 225, num=13)
        w1 = trive.w1
        w2 = trive.w2
        ws = trive.ws
        d2 = trive.d2
        d1 = trive.d1
        w1.points = np.linspace(6500, 7500, num=11)
        w2.points = np.linspace(6500, 7500, num=11)
        #trive.ws.points = np.linspace(6000, 8800, num=21)
        d2.points = np.linspace(-200, 200, num=21)
        d1.points = np.linspace(-200, 200, num=21)
        trive.exp.timestep = 2.0
        if use_inhom:
            #inhom_object = inhom.Inhom(inhom_sampling='linear', num=5, sigma=50.0)
            if gh_quad:
                inhom_object = inhom.Inhom(inhom_sampling='gh', 
                                           num=7, sigma=100.0)
            else:
                inhom_object = inhom.Inhom(inhom_sampling='linear', 
                                           num=100, sigma=100.0, 
                                           zeta_bound=2)
        else:
            inhom_object = inhom.Inhom()
        H = H_1.Omega()
        H.TOs = TOs
        H.tau_ag  = 10.
        H.tau_2aa = 10.
        H.tau_2ag = 10.
        trive.late_buffer = 250.
        #out1 = trive.exp.scan(w1, w2, d2, H=H, inhom_object=inhom_object)
        out1 = trive.exp.scan(d1, d2, H=H, inhom_object=inhom_object)
        if mp:
            out1.run(autosave=autosave, mp=True, chunk=False)
        else:
            out1.run(autosave=autosave, mp=False)
    
    if add_nrb:
        nrb1 = out1.nrb()
        print nrb1.shape
        print out1.sig.shape
        nrb1_max = np.abs(nrb1).sum(axis=-1).max()
        out1_smax = np.abs(out1.sig.sum(axis=-1).sum(axis=-2)).max()
        #print nrb1_max, out1_smax
        nrb1 *= out1_smax / nrb1_max * nrb_level
        #print np.abs(nrb1).max(), np.abs(out1.sig).max()
        # careful...this colapses out_groups
        out1.sig = out1.sig.sum(axis=-2)[...,None,:]
        out1.sig += nrb1[...,None,:]
    
    #outsig = out1.sig.copy()
    #for i in range(2):
    #out1.sig = outsig[...,i:i+1,:]
    sig1 = m.Measure(out1, m.Mono, m.SLD)
    # should be equivalent to the others if we set mono to pass
    m.Mono.slitwidth=slitwidth
    sig1.run()
    if plot:
        print sig1.pol.shape
        # run out1.axes() if you forgot which axis is which
        sig1.plot(0,yaxis=1, zoom=2)
        #sig1.plot(0,yaxis=2, zoom=2)
        #sig1.plot(1,yaxis=2)
