# -*- coding: utf-8 -*-
"""
Created on Wed Apr 08 16:02:45 2015

@author: Dan
"""

if __name__ == '__main__':

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
    imported = True
    add_nrb   = False
    use_inhom = False
    plot      = False
    mp        = True
    gh_quad   = False
    view_kernel = False
    TOs       = [5]
    nrb_level = 1. # amplitude, relative to sig max
    smear_extent = 1000. # inohomogeneous broadening fwhm, in wn
    slitwidth = 120.#None#120.0 # width in wavenumbers; set to none if no slit
    # takes ~ 50 pts per second (with no inhomogeneity)
    angle = np.pi / 4.
    rel_path1 = r'\data\2015.04.08 16-18-04'
    #rel_path2 = r'\data\2015.04.07 20-55-42 - correlated'
    rel_path2 = r'\data\2015.04.08 16-18-37'
    if imported:
        filepath1 = r''.join([NISE_path, rel_path1])
        filepath2 = r''.join([NISE_path, rel_path2])
        out1 = trive.S.Scan._import(filepath1)
        out2 = trive.S.Scan._import(filepath2)
        w1 = out1.axis_objs[0]
        w2 = out1.axis_objs[1]
    else:
        trive.exp.set_coord(trive.d1, 0.)
        trive.exp.set_coord(trive.d2, 0.)
        trive.exp.set_coord(trive.ss, 45.)
        #trive.d2.points = np.array([-75, -25, 0, 25, 75])#, 150, 250, 500])
        #np.linspace(-75, 225, num=13)
        w1 = trive.w1
        w2 = trive.w2
        ws = trive.ws
        d2 = trive.d2
        d1 = trive.d1
        w1.points = np.linspace(5000, 9000, num=41)
        w2.points = np.linspace(5000, 9000, num=41)
        #trive.ws.points = np.linspace(6000, 8800, num=21)
        #d2.points = np.linspace(-50, 100, num=4)
        #d1.points = np.linspace(-125, 125, num=15)
        trive.exp.timestep = 2.0
        trive.exp.late_buffer = 250.
        inhom_object = inhom.Inhom()
        H = H_1.Omega()
        H.TOs = TOs
        H.wa_central = 7000
        out1 = trive.exp.scan(w1, w2, H=H, inhom_object=inhom_object)

        H2 = H_1.Omega()
        H2.TOs=TOs
        H2.wa_central = 6500
        out2 = trive.exp.scan(w1, w2, H=H2, inhom_object=inhom_object)
        if mp:
            out1.run(autosave=autosave, mp=True, chunk=False)
            out2.run(autosave=autosave, mp=True, chunk=False)
        else:
            out1.run(autosave=autosave, mp=False)
            out2.run(autosave=autosave, mp=False)

    if add_nrb:
        for obj in [out1, out2]:
            nrb1 = obj.nrb()
            print nrb1.shape
            print obj.sig.shape
            nrb1_max = np.abs(nrb1).sum(axis=-1).max()
            obj_smax = np.abs(obj.sig.sum(axis=-1).sum(axis=-2)).max()
            #print nrb1_max, out1_smax
            nrb1 *= obj_smax / nrb1_max * nrb_level
            #print np.abs(nrb1).max(), np.abs(out1.sig).max()
            # careful...this colapses out_groups
            obj.sig = obj.sig.sum(axis=-2)[...,None,:]
            obj.sig += nrb1[...,None,:]

    tprime = np.arange(out1.sig.shape[-1])*out1.timestep
    """
    import matplotlib.pyplot as plt
    from NISE.lib.misc import wn_to_omega
    obj = out1
    zmax=np.abs(obj.sig).max()
    for i in range(out1.sig.shape[-1]):
        plt.figure()
        plt.title('t={0}fs'.format(tprime-out1.early_buffer))
        plt.subplot(121, aspect='equal')
        #zz1 = np.imag(np.exp(1j*w1.points[None,:]
        #             * wn_to_omega * tprime[i]) * obj.sig[:,:,1,0,i].T)
        zz1 = np.imag(np.exp(1j*7000.* wn_to_omega * tprime[i]) * obj.sig[:,:,-1,0,i].T)
        #zz1 = np.fft.fft(zz1)
        plt.contourf(zz1/zmax, levels=np.linspace(-1,1, num=20))
        plt.grid(b=True)
        plt.title('imaginary')
        plt.subplot(122, aspect='equal')
        #zz2 = np.real(np.exp(1j*w1.points[None,:]
        #             * wn_to_omega * tprime[i]) * obj.sig[:,:,1,0,i].T)
        zz2 = np.real(np.exp(1j*7000.* wn_to_omega * tprime[i]) * obj.sig[:,:,-1,0,i].T)
        #zz2 = np.fft.fft(zz2)
        plt.contourf(zz2/zmax, levels=np.linspace(-1,1, num=20))
        plt.grid(b=True)
        plt.title('real')
        plt.savefig(obj.output_folder+r'\t{0}'.format(i))
        plt.close()
        print i
    1/0.
    #"""
    sig1 = m.Measure(out1, m.Mono, m.SLD)
    sig2 = m.Measure(out2, m.Mono, m.SLD)
    # should be equivalent to the others if we set mono to pass
    m.Mono.slitwidth=slitwidth
    if plot:
        sig1.run()
        sig2.run()
        print sig2.pol.shape
        # run out1.axes() if you forgot which axis is which
        sig1.plot(0,yaxis=1, zoom=2)
        sig2.plot(0,yaxis=1, zoom=2)
        #sig1.plot(0,yaxis=2, zoom=2)
        #sig1.plot(1,yaxis=2)
