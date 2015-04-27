# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 19:30:47 2015

@author: Dan
"""

from NISE.lib.misc.__init__ import NISE_path
#import matplotlib.pyplot as plt
import numpy as np
# important:  we can't call on scan directly for some reason; use trive to do it
import NISE.experiments.trive as trive
import NISE.hamiltonians.H0 as H_1
import NISE.hamiltonians.params.inhom as inhom

folder_path = NISE_path

trive.exp.set_coord(trive.d1, 0.)
trive.exp.set_coord(trive.d2, 0.)
trive.exp.set_coord(trive.ss, 45.)
#trive.d2.points = np.array([-75, -25, 0, 25, 75])#, 150, 250, 500])
#np.linspace(-75, 225, num=13)

def get_pixel(
                TOs       = [1,2,3,4,5,6],
                out_group = None,
                mp        = True,
                use_inhom = False,
                w1s = None,
                w2s = None,
                d2s = None,
                d1s = None
            ):
    w1 = trive.w1
    w2 = trive.w2
    d2 = trive.d2
    d1 = trive.d1
    if w1s is None: w1.points=np.array([H_1.Omega.wa_central])
    else: w1.points = w1s
    if w2s is None: w2.points=np.array([H_1.Omega.wa_central])
    else: w2.points = w2s
    if d2s is None: d2.points=np.array([0])
    else: d2.points = d2s
    if d1s is None: d1.points=np.array([0])
    else: d1s.points = d1s
    
    #w1.points = np.linspace(5000, 9000, num=61)
    #w2.points = np.linspace(5000, 9000, num=61)
    #trive.ws.points = np.linspace(6000, 8800, num=21)
    #d2.points = np.linspace(-100, 100, num=11)
    #d1.points = np.linspace(-125, 125, num=15)
    trive.exp.timestep = 2.0
    if use_inhom:
        inhom_object = inhom.Inhom(inhom_sampling='linear', 
                                   num=100, sigma=100.0, 
                                   zeta_bound=2)
    else:
        inhom_object = inhom.Inhom()
    H = H_1.Omega()
    if out_group is None:  pass
    else: H.out_group = out_group
    H.TOs = TOs
    trive.exp.late_buffer = 250.
    out1 = trive.exp.scan(w1, w2, d1, d2, H=H, inhom_object=inhom_object)
    if out1.array.size < 4:
        mp = False
    out1.run(autosave=False, mp=True, chunk=False)
    return out1
