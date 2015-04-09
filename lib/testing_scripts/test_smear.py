# -*- coding: utf-8 -*-
"""
Created on Tue Apr 07 17:30:13 2015

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
    imported_before = True
    imported_after  = True
    add_nrb   = False
    use_inhom = False
    plot      = True
    mp        = True
    gh_quad   = False
    view_kernel = True
    make_movie = False
    TOs       = [5]
    nrb_level = 0.5 # amplitude, relative to sig max
    smear_extent = 500. # inohomogeneous broadening fwhm, in wn
    slitwidth = None#120.0 # width in wavenumbers; set to none if no slit
    # takes ~ 50 pts per second (with no inhomogeneity)
    angle = np.pi / 4.
    rel_path1 = r'\data\2015.04.09 14-59-08'
    #rel_path2 = r'\data\2015.04.07 20-55-42 - correlated'
    rel_path2 = r'\data\2015.04.09 16-03-37'
    if imported_before:
        filepath = r''.join([NISE_path, rel_path1])
        out1 = trive.S.Scan._import(filepath)
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
        w1.points = np.linspace(5000, 9000, num=31)
        w2.points = np.linspace(5000, 9000, num=30)
        #trive.ws.points = np.linspace(6000, 8800, num=21)
        d2.points = np.linspace(-50, 100, num=4)
        #d1.points = np.linspace(-125, 125, num=15)
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
        trive.exp.late_buffer = 250.
        out1 = trive.exp.scan(w1, w2, d2, H=H, inhom_object=inhom_object)
        #out1 = trive.exp.scan(d1, d2, H=H, inhom_object=inhom_object)
        if mp:
            out1.run(autosave=autosave, mp=True, chunk=False)
        else:
            out1.run(autosave=autosave, mp=False)
    if imported_after:
        filepath = r''.join([NISE_path, rel_path2])
        out2 = trive.S.Scan._import(filepath)
    else:
        sig2 = out1.smear(0,1,smear_extent,50,theta=angle)
        from copy import deepcopy
        #print 'before deepcopy'
        #print out1.efp[15,15,0]
        out2 = deepcopy(out1)
        #print 'after deepcopy'
        #print out1.efp[15,15,0]
        out2.sig = sig2
        out2.save()
    if view_kernel:
        import matplotlib.pyplot as plt
        plt.figure()
        s_maj = smear_extent / 2*np.sqrt(2*np.log(2)) 
        s_min = 50 / 2*np.sqrt(2*np.log(2))
        x = out1.axis_objs[0].points.copy()
        y = out1.axis_objs[1].points.copy()
        x -= x.sum() / x.size
        y -= y.sum() / y.size
        uu = np.cos(angle)*x[None,:] + np.sin(angle) * y[:,None]
        vv = np.cos(angle)*y[:,None] - np.sin(angle) * x[None,:]
        uu /= np.sqrt(2) * s_maj
        vv /= np.sqrt(2) * s_min
        kk =  np.exp(-uu**2 -vv**2)
        # normalize the distribution
        kk /= kk.sum()
        import NISE.lib.fscolors_beta as f
        plt.contourf(x,y, kk, 200, cmap=f.plot_artist.mycm)
        plt.colorbar()
        plt.savefig(out2.output_folder + r'\smear kernel.png')

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
    
    # little bitty for looking at spectral evolution
    # ensuring that phase is working the correct way
    tprime = np.arange(out1.sig.shape[-1])*out1.timestep - out1.early_buffer

    import matplotlib.pyplot as plt
    from NISE.lib.misc import wn_to_omega
    if make_movie:
        for obj in [out1,out2]:
            zmax=np.abs(obj.sig).max()
            for i in range(out1.sig.shape[-1]):
                plt.figure()
                plt.title('t={0}fs'.format(tprime-out1.early_buffer))
                plt.subplot(121, aspect='equal')
                #zz1 = np.imag(np.exp(1j*w1.points[None,:]
                #             * wn_to_omega * tprime[i]) * obj.sig[:,:,1,0,i].T)
                zz1 = np.imag(np.exp(1j*w1.points[None,:]* wn_to_omega * tprime[i]) * obj.sig[:,:,1,0,i].T)
                #zz1 = np.fft.fft(zz1)
                plt.contourf(zz1/zmax, levels=np.linspace(-1,1, num=20))
                plt.grid(b=True)
                plt.title('imaginary')
                plt.subplot(122, aspect='equal')
                #zz2 = np.real(np.exp(1j*w1.points[None,:]
                #             * wn_to_omega * tprime[i]) * obj.sig[:,:,1,0,i].T)
                zz2 = np.real(np.exp(1j*w1.points[None,:]* wn_to_omega * tprime[i]) * obj.sig[:,:,1,0,i].T)
                #zz2 = np.fft.fft(zz2)
                plt.contourf(zz2/zmax, levels=np.linspace(-1,1, num=20))
                plt.grid(b=True)
                plt.title('real')
                plt.savefig(obj.output_folder+r'\t{0}'.format(i))
                plt.close()
                print i

    sig1 = m.Measure(out1, m.Mono, m.SLD)
    sig2 = m.Measure(out2, m.Mono, m.SLD)
    # should be equivalent to the others if we set mono to pass
    m.Mono.slitwidth=slitwidth
    if plot:
        sig1.run()
        sig2.run()
        #print sig2.pol.shape
        # run out1.axes() if you forgot which axis is which
        sig1.plot(0,yaxis=1, zoom=2)
        sig2.plot(0,yaxis=1, zoom=2)
        #sig1.plot(0,yaxis=2, zoom=2)
        #sig1.plot(1,yaxis=2)
