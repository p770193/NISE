# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 00:29:23 2015

@author: Dan
"""

import numpy as np
import NISE.hamiltonians.params.inhom as inhom
import NISE.experiments.trive as trive
import NISE.hamiltonians.H_PbSe2 as H_1
reload(H_1)

autosave  = False
imported  = False
add_nrb   = False
use_inhom = False
nrb_level = 1. # amplitude, relative to sig max
slitwidth = 120.0 # width in wavenumbers; set to none if no slit
# takes ~ 50 pts per second (with no inhomogeneity)

trive.exp.set_coord(trive.d1, 0.)
trive.exp.set_coord(trive.d2, 0.)
trive.exp.set_coord(trive.ss, 45.)
trive.w1.points = np.array([8300.])
trive.w2.points = np.array([6700.])
#trive.d2.points = np.array([-75, -25, 0, 25, 75])#, 150, 250, 500])
#np.linspace(-75, 225, num=13)
w1 = trive.w1
w2 = trive.w2
ws = trive.ws
d2 = trive.d2
trive.exp.timestep = 2.0
if use_inhom:
    inhom_object = inhom.Inhom(inhom_sampling='linear', num=10, sigma=100.0)
else:
    inhom_object = inhom.Inhom()
H = H_1.Omega()
H.out_group = [[13,14,15,16],[12,17]]
H.Gamma_bb  = 0. # eliminate coupling to examine cross-peak
out1 = trive.exp.scan(w1, w2, H=H, inhom_object=inhom_object)
out1.run(autosave=autosave)
