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
    imported    = True
    add_nrb     = False
    use_inhom   = False
    plot1       = True
    plot2       = True
    gh_quad     = False

    TOs         = np.arange(6) + 1
    nrb_level   = 5
    tau = 25. #set all coherence times the same
    slitwidth = 120.
    
    folder = NISE_path + r'\data'
    # imports of 2d delay
    #rel_name = r'\2015.04.29 19-44-25'
    #rel_name2 = r'\2015.04.29 19-45-34'
    # imports of w1w2d1 scans
    rel_name = r'\2015.04.29 22-48-47 - w1w2d1 no chirp - driven'
    rel_name2 = r'\2015.04.29 22-55-36 - w1w2d1 chirp - driven'

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
        # going to do two scans:  without chirp, with chirp
        t.exp.set_coord(t.d1, 0.)
        t.exp.set_coord(t.d2, 0.)
        t.exp.set_coord(t.ss, 45.)
        #t.d2.points = np.array([-75, -25, 0, 25, 75])#, 150, 250, 500])
        #np.linspace(-75, 225, num=13)
        w1 = t.w1
        w2 = t.w2
        ws = t.ws
        d2 = t.d2
        d1 = t.d1
        w1.points = np.linspace(5000, 9000, num=41)
        w2.points = np.linspace(5000, 9000, num=41)
        #trive.ws.points = np.linspace(6000, 8800, num=21)
        d2.points = np.linspace(-200, 200, num=41)
        d1.points = np.linspace(-200, 200, num=11)
        t.exp.timestep = 2.0

        inhom_object = inhom.Inhom()
        H_1.Omega.tau_ag = tau
        H_1.Omega.tau_2aa = tau
        H_1.Omega.tau_2ag = tau
        H = H_1.Omega()
        H.TOs = TOs
        # make early buffer a little longer due to chirp
        t.early_buffer = 120.
        t.exp.late_buffer = H.tau_ag * 5
        #out1 = t.exp.scan(w1, w2, d2, H=H, inhom_object=inhom_object)
        t.exp.set_coord(dzs, 0.)
        #out1 = t.exp.scan(d1, d2, H=H, inhom_object=inhom_object)
        out1 = t.exp.scan(w1, w2, d1, H=H, inhom_object=inhom_object)
        out1.run(autosave=autosave, mp=True, chunk=False)
        t.exp.set_coord(dzs, 1.5e-5)
        out2 = t.exp.scan(w1, w2, d1, H=H, inhom_object=inhom_object)
        out2.run(autosave=autosave, mp=True, chunk=False)

    if add_nrb:
        for obj in [out1, out2]:
            nrb1 = obj.nrb()
            print nrb1.shape
            print obj.sig.shape
            nrb1_max = np.abs(nrb1).max()
            obj_smax = np.abs(obj.sig.sum(axis=-2)).max()
            #print nrb1_max, out1_smax
            nrb1 *= obj_smax / nrb1_max * nrb_level
            #print np.abs(nrb1).max(), np.abs(out1.sig).max()
            # careful...don't collapse the out_groups axis
            obj.sig = obj.sig.sum(axis=-2)[...,None,:]
            obj.sig += nrb1[...,None,:]

    m1 = m.Measure(out1, m.Mono, m.SLD)
    m.Mono.slitwidth=slitwidth
    m2 = m.Measure(out2, m.Mono, m.SLD)
    m.Mono.slitwidth=slitwidth

    if plot1:
        m1.run()
        m1.plot(0,yaxis=1, zoom=1)
    if plot2:
        m2.run()
        m2.plot(0,yaxis=1, zoom=1)

