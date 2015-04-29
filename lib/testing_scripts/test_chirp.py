# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 14:09:50 2015

testing the Gauss_chirp_rwa pulse class

@author: Dan
"""

if __name__ == '__main__':

    from NISE.lib.misc.__init__ import NISE_path
    #import matplotlib.pyplot as plt
    import numpy as np
    # important:  we can't call on scan directly for some reason; use trive to do it
    import NISE.lib.measure as m
    import NISE.experiments.trive as t
    # to change the pulse class from the default and use these definitions, 
    # we have to assert the new name
    t.assert_pulse_class('Gauss_chirp_rwa')
    import NISE.hamiltonians.H0 as H_1
    import NISE.hamiltonians.params.inhom as inhom
    
    autosave    = True
    imported    = False
    add_nrb     = False
    use_inhom   = False
    plot1       = True
    plot2       = True
    gh_quad     = False
    view_kernel = True
    TOs         = np.arange(6) + 1
    slitwidth = 120.
    
    folder = None
    rel_name = None
    rel_name2 = None

    # define the chirp axis for positioning
    dzs = t.S.Axis(0, 'dz', also=[1,2], name=r'$dz_1=dz_2=dz_{2^\prime}$',
                                      pulse_class_name='Gauss_chirp_rwa')

    if imported:
        filepath = r''.join([folder, rel_name])
        out1 = t.S.Scan._import(filepath)
        filepath = r''.join([folder, rel_name2])
        out2 = t.S.Scan._import(filepath)
        d1 = out1.axis_objs[0]
        d2 = out1.axis_objs[1]
    else:
        # going to do three scans:  without chirp, with chirp, with more chirp
        t.exp.set_coord(t.d1, 0.)
        t.exp.set_coord(t.d2, 0.)
        t.exp.set_coord(t.ss, 35.)
        #t.d2.points = np.array([-75, -25, 0, 25, 75])#, 150, 250, 500])
        #np.linspace(-75, 225, num=13)
        w1 = t.w1
        w2 = t.w2
        ws = t.ws
        d2 = t.d2
        d1 = t.d1
        w1.points = np.linspace(5000, 9000, num=121)
        w2.points = np.linspace(5000, 9000, num=121)
        #trive.ws.points = np.linspace(6000, 8800, num=21)
        d2.points = np.linspace(-150, 150, num=21)
        d1.points = np.linspace(-150, 150, num=21)
        t.exp.timestep = 2.0

        inhom_object = inhom.Inhom()

        H = H_1.Omega()
        H.tau_ag = 30
        H.TOs = TOs
        t.exp.late_buffer = H.tau_ag * 5
        #out1 = t.exp.scan(w1, w2, d2, H=H, inhom_object=inhom_object)
        t.exp.set_coord(dzs, 0.)
        out1 = t.exp.scan(d1, d2, H=H, inhom_object=inhom_object)
        out1.run(autosave=autosave, mp=True, chunk=False)
        t.exp.set_coord(dzs, 1.5e-5)
        out2 = t.exp.scan(d1,d2, H=H, inhom_object=inhom_object)
        out2.run(autosave=autosave, mp=True, chunk=False)

    m1 = m.Measure(out1, m.Mono, m.SLD)
    m.Mono.slitwidth=slitwidth
    m2 = m.Measure(out2, m.Mono, m.SLD)
    m.Mono.slitwidth=slitwidth

    if plot1:
        m1.run()
        m1.plot(0,yaxis=1, zoom=2)
    if plot2:
        m2.run()
        m2.plot(0,yaxis=1, zoom=2)

