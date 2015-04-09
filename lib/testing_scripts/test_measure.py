# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 01:00:56 2014

testing measure functionality

@author: Dan
"""
from NISE.lib.misc.__init__ import NISE_path
import matplotlib.pyplot as plt
import numpy as np
#import NISE.lib.scan as S
import NISE.lib.measure as m
#import NISE.hamiltonians.H_sf as H_sf
import NISE.experiments.trive as trive
reload(trive)
reload(m)

imported = False
slitwidth = None # width in wavenumbers; set to none if no slit

trive.w1.points = np.linspace(6700, 7300, num=21)
trive.w2.points = np.linspace(6700, 7300, num=21)

trive.d1.points = np.linspace(-100, 200, num=21)
trive.d2.points = np.linspace(-100, 200, num=21)

if imported:
    filepath = ''.join([NISE_path, r'\data\2014.08.05 00-46-34'])
    out2 = trive.S.Scan._import(filepath)
else:
    trive.exp.set_coord(trive.ws, 7000.)
    #H = H_sf.Omega(D=4)
    out2 = trive.exp.scan(trive.d1, trive.d2)
    #out2.run()
    #trive.set_coord(d1, 0.)
    #trive.set_coord(d2, 0.)
    #out1 = trive.scan(w1, w2)
"""
# manually plotting the output of out2
sig2o = (np.abs(out2.sig.sum(axis=-2))**2).sum(axis=-1) * out2.timestep
plt.figure()
plt.subplot(131)
plt.contourf(trive.d1.points, trive.d2.points, sig2o, 200)
plt.title('sig2 manual')
plt.grid()
plt.colorbar()

# using measure processes - shoudl be equivalent to sig2o
sig2 = m.Measure(out2, m.SLD)
sig2.run()
plt.subplot(132)
plt.contourf(trive.d1.points, trive.d2.points, sig2.pol, 200)
plt.title('sig2 measure')
plt.grid()
plt.colorbar()

ratio1 = sig2.pol / sig2o
print ratio1.min(), ratio1.max()

# using measure processes, but with a mono
sig2m = m.Measure(out2, m.Mono, m.SLD)
# should be equivalent to the others if we set mono to pass
m.Mono.slitwidth=slitwidth
sig2m.run()
plt.subplot(133)
plt.contourf(trive.d1.points, trive.d2.points, sig2m.pol, 200)
plt.title('sig2 measure mono with slitwidth {0}'.format(slitwidth))
plt.grid()
plt.colorbar()

ratio2 = sig2m.pol / sig2o
print ratio2.min(), ratio2.max()
"""