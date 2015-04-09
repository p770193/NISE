# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 11:59:06 2014

testing 2d frequency processing

@author: Dan
"""

from NISE.lib.misc.__init__ import NISE_path
import matplotlib.pyplot as plt
import numpy as np
# important:  we can't call on scan directly for some reason; use trive to do it
#import NISE.lib.scan as S
import NISE.lib.measure as m
import NISE.experiments.trive as trive
reload(trive)
#reload(S)
reload(m)

imported = True
slitwidth = 120.0 # width in wavenumbers; set to none if no slit
rel_path = r'data\2014.08.05 00-51-19'

if imported:
    filepath = r'//'.join([NISE_path, rel_path])
    out1 = trive.S.Scan._import(filepath)
    w1 = out1.axis_objs[0]
    w2 = out1.axis_objs[1]
else:
    trive.exp.set_coord(trive.d1, 0.)
    trive.exp.set_coord(trive.d2, 0.)
    trive.w1.points = np.linspace(6700, 7300, num=6)
    trive.w2.points = np.linspace(6700, 7300, num=6)
    w1 = trive.w1
    w2 = trive.w2
    out1 = trive.exp.scan(w1, w2)
    out1.run()

sig1 = m.Measure(out1, m.Mono, m.SLD)
# should be equivalent to the others if we set mono to pass
m.Mono.slitwidth=slitwidth
sig1.run()
plt.figure()
plt.contourf(w1.points, w2.points, sig1.pol, 200)
plt.title('sig2 measure mono with slitwidth {0}'.format(slitwidth))
plt.grid()
plt.colorbar()
    

