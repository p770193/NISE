"""
Created on Wed Jun 13 17:06:55 2012
@author: Dan
Using Domcke approach to numerically integrate TRIVE signal of diagonal and near-diagonal features
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
import matplotlib.pyplot as plt
wn_to_omega = 2*np.pi*3*10**-5  #omega is radians / fs
#simulation time properties, in fs
timestep = 4.0  
early_buffer = 100.0
late_buffer  = 400.0
#--------------------------QD Properties--------------------------
#inhomogeneous linewidth--stdev in wavenumbers
sigma = 610.0 / (2*np.sqrt(2*np.log(2)))
#exciton-exciton coupling
a_coupling = 150.0
#1S exciton central position
wa_central = 7725.
#dephasing times, 1/fs
Gamma_gg  = 0.0
Gamma_ag  = 1./80
Gamma_ga  = Gamma_ag
Gamma_aa  = 1./4000
Gamma_2aa = 1./30
Gamma_2ag = 1./30
#transition dipoles (a.u.)
mu_ag =  1.0
mu_2aa = 1.5 * mu_ag
#--------------------------Solvent (ISR) Properties--------------------------
#define frequencies of states (wavenumbers), lifetimes(rad/fs), and dipoles
#frequencies from NIST
#              ag     bg     cg     dg     eg
ctet_w      = [217.,  313.,  459.,  762.,  790.]
ctet_gamma  = [0.0002,0.0002,0.0002,0.0002,0.0002]
s_mu        = [0.45,  0.45,  1,     0.22, 0.22]
#--------------------------Pulse Properties--------------------------
#w1, w2, and w2p pulse parameters in time domain (and intensity scale)
E1amp       = 1.0
E2amp       = 1.0
E2pamp      = 1.0
#pulse stdev amplitude level
sigma_1     = 60 / (2*np.sqrt(np.log(2)))
sigma_2     = 60 / (2*np.sqrt(np.log(2)))
sigma_2p    = 60 / (2*np.sqrt(np.log(2)))
#--------------------------Function Definitions--------------------------

def rotor(theta):
    """
    returns 1 rotated CW about the complex plane by angle theta i.e. e^(-i*theta)
    can work on scalars as well as higher rank objects
    """
    #return complex(np.cos(theta), -np.sin(theta))
    return np.cos(theta)*complex(1,0) - np.sin(theta)*complex(0,1)

def pulse(vec, amp, sigma, mu, freq):
    """
    returns an array of gaussian values with arguments of 
        input array vec, 
        mean mu, 
        standard deviation sigma
        amplitude amp
        frequency freq
    """
    x = np.zeros((len(vec)), dtype=np.complex64)
    y = amp / (sigma*np.sqrt(2*np.pi))
    #incorporate pulse envelope function
    for i in range(len(vec)):
        x[i] = y * np.exp(-(vec[i] - mu)**2 / (2*sigma**2) ) * rotor(freq*(vec[i]-mu))
    return x

def sOmega(E1, E2, E2p, t, wvg, mu_vg, w1first = True):
    """
    creates the solvent coupling array given the input e-fields values for a specific time, t
    w1first selects whether w1 or w2p is the first interacting positive field
    accounts for CARS interactions
    """
    if w1first==True:
        first  = E1
    else:
        first  = E2p
    O = np.zeros((len(t),3,3),dtype=np.complex64)
    #from gg
    O[:,1,0] = 0.5*mu_vg * E2 * first * rotor(-wvg*t)
    O[:,2,0] = 0.5*mu_vg * E2 * first * rotor( wvg*t)
    #emitted coherences derived from appendix c of notes
    #these elements are not derivatives--they are the actual evolution of the coherence
    return O

def Omega(E1, E2, E2p, t, ws, w1first = True):
    """
    creates the coupling array given the input e-fields values for a specific time, t
    w1first selects whether w1 or w2p is the first interacting positive field
    
    Currently neglecting pathways where w2 and w2' require different frequencies
    """
    wag  = ws[1]
    w2aa = ws[7]
    if w1first==True:
        first  = E1
        second = E2p
    else:
        first  = E2p
        second = E1
    O = np.zeros((len(t),8,8), dtype=np.complex64)
    #from gg1
    O[:,1,0] =  mu_ag  * first  * rotor(-wag*t)
    O[:,2,0] = -mu_ag  * E2     * rotor(wag*t)
    #from ag1
    O[:,4,1] = -mu_ag  * E2     * rotor(wag*t)
    O[:,5,1] =  mu_2aa * second * rotor(-w2aa*t)
    O[:,3,1] =  mu_ag  * E2     * rotor(wag*t)
    #from ga
    O[:,4,2] =  mu_ag  * first  * rotor(-wag*t)
    O[:,3,2] = -mu_ag  * first  * rotor(-wag*t)
    #from gg2
    O[:,6,3] =  mu_ag  * second * rotor(-wag*t)      * mu_ag
    #from aa
    O[:,7,4] =  mu_2aa * second * rotor(-w2aa*t)     * mu_2aa
    O[:,6,4] = -mu_ag  * second * rotor(-wag*t)      * mu_ag
    #from 2ag
    O[:,7,5] = -mu_ag  * E2     * rotor(wag*t)       * mu_2aa
    O[:,6,5] =  mu_2aa * E2     * rotor(w2aa*t)      * mu_ag
    return O

def inhom_dist(inhom):
    """
    generates the list of sampling points in the distribution and their weights
    """
    if inhom == 'hermite':
        #uses 64th degree Hermite poly zeroes and Christoffel numbers
        inhom_dat = np.genfromtxt('C:\Users\Dan\_ipython\hermite 32v2.txt',dtype=np.float32)
        inhom_dat = inhom_dat.transpose()
        zeta = inhom_dat[0] * np.sqrt(2.)
        zweight = inhom_dat[1]
        dzeta = 1.0
    elif inhom == 'legendre':
        #uses 80th degree legendre poly zeroes and Christoffel numbers
        inhom_dat = np.genfromtxt('C:\Users\Dan\_ipython\legendre 40.txt', dtype=np.float32)
        inhom_dat = inhom_dat.transpose()
        zeta = inhom_dat[0]*3
        zweight = 1 / (np.sqrt(2*np.pi)) * np.exp(-(zeta**2)/2)*20
        zweight = zweight * inhom_dat[1]
        dzeta = 1.0
    elif inhom == 'linear':
        #currently the only inhomogeneity parameter that can normalize well with no inhomogeneity
        zeta = np.linspace(-2.0,2.0, num=15)
        #dzeta = zeta[1] - zeta[0]
        zweight = 1 / (np.sqrt(2*np.pi)) * np.exp(-(zeta**2)/2)
        dzeta = np.abs(zeta[1] - zeta[0])
    else:
        dzeta = 1.0
        zeta = np.array([0])
        zweight = [1.0]
    return zeta, zweight, dzeta
    
def energy_dist(zeta, wa_central, a_coupling, l):
    """
    creates the arrays of oscillator energies
    """
    zl = zeta[l]
    wg = 0.0
    wa = zl + wa_central
    w2a = 2*wa - a_coupling
    w_ag = wa - wg
    w_aa = wa - wa
    w_gg = wg - wg
    w_2ag = w2a - wg
    w_2aa = w2a - wa
    #array aggregates all frequencies to match state vectors
    w     = np.array( [w_gg, w_ag, -w_ag, w_gg, w_aa, w_2ag, w_ag, w_2aa] )
    #w     = np.array( [w_gg, w_ag, -w_ag, w_aa, w_2ag, w_2aa, w_gg, w_ag] )
    return w

def measure(data, timestep, slitwidth, offset=0.0, 
            plot_dimension=None, sig_type=None):
    """
        enacts the act of measurement of the emitted coherence, decides with plot_dimension:
            'time': a time-resolved photodetector (complex-valued array)
            'mono': frequency-windowed intensity (complex-valued array)
            'int' : integrated intensity (no mono)
            None  : a photodetector with mono window applied (single real number out) 
    """
    if plot_dimension == 'time':
        #all calculations for coherence as func of time is done automatically
        #so no need to call this function!
        return data
    elif plot_dimension == 'mono':
        #data freq is relative to wm
        out = np.fft.fft(data)
        out = np.fft.fftshift(out)
        fft_tprime = np.fft.fftfreq(len(data), d=timestep)
        var = np.fft.fftshift(fft_tprime) *2*np.pi / wn_to_omega
        #enact monochromator
        keep_lis = []
        for i in range(len(var)):
            if np.abs(var[i]-offset) <= slitwidth / 2.0:
                keep_lis.append(i)
        out = out[keep_lis]
        if sig_type == 'real':
            out = np.real(out)
        elif sig_type == 'imag':
            out = np.imag(out)
        elif sig_type == 'phase':
            outi = np.imag(out)
            outr = np.real(out)
            out = np.arctan(outi/outr)
        elif sig_type == 'amp':
            out = np.abs(out)
        else:
            out = np.abs(out)**2
            #not sure about why I multiply by timestep here...
        increment = timestep / len(data) #* 2*np.pi / wn_to_omega
        return out.sum()*increment
    elif plot_dimension == 'int':
        #integrate in time the electric field fluence
        if sig_type == 'real':
            out = np.real(data)
        elif sig_type == 'imag':
            out = np.imag(data)
        elif sig_type == 'phase':
            outi = np.imag(data)
            outr = np.real(data)
            out = np.arctan(outi/outr)
        elif sig_type == 'amp':
            out = np.abs(data)
        else:
            out = np.abs(data)**2
        return out.sum() * timestep
    else:
        print 'no plot_dimension specified!'

def get_lengths(timestep):
    """
    report back the length of t_prime array given timestep (buffer times are a fixed property of this file) 
    """
    t   = np.arange(-early_buffer, late_buffer, timestep)
    return len(t)

def get_t(t21,t2p1, timestep):
    """
        returns the t array used in pixel
        also returns tprime_i, the cutoff for when to record emitted coherences
    """    
    if t21 < 0 or t2p1 < 0:
        #w2 or w2p will come at t<0
        t_min = min(t21, t2p1) - early_buffer
    else:
        t_min = -early_buffer
    #last pulse dictates when emitted coherences should start to be watched for integration
    if t21 > 0 or t2p1 > 0:
        #w2 or w2p will come at t>0
        t_max = max(t21, t2p1) + late_buffer
    else:
        t_max = late_buffer
    #note:  dtype float16 can screw up t working properly!
    t   = np.arange(t_min, t_max, timestep)
    #tprime_i = t_max - late_buffer - early_buffer
    return t

def pixel(
    t21, t2p1, w1, w2, w2p, wm=None,
    wa_central=wa_central, a_coupling=a_coupling, sigma=sigma,
    timestep=timestep, 
    inhom=None, 
    solvent=False,
    scatter=None,
    plot_dimension=None
    ):
    """
    calculates the TRIVE signal of the 1S exciton given pulse and detection parameters
    output is a variety of different coherences rho_i(t) all in the same rotating frame
        t21 (d2) and t2p1 (d1) in fs
        w1,w2,w2p, in wavenumbers
        timestep in fs
        
        inhom dictates how inhomogeneous broadening is accounted for
            'linear' uses Reimann sums to integrate over the distribution
            'hermite' uses hermite mechanical quadrature (assume gaussian weight)
            'legendre' uses legendre mechanical quadrature (add weighting to function manually) 
            default (None) incorporates no inhomogeneous broadening
        solvent is a boolean--True if simulation should include solvent CARS contributions
    """
    print 'w1',w1
    print 'w2',w2
    print 't21',t21
    print 't2p1',t2p1 
    w1 = w1 *wn_to_omega
    w2 = w2 *wn_to_omega
    w2p= w2p*wn_to_omega
    #rounding to keep from errors in the array length
    #t21 = np.round(t21,3)
    #t2p1= np.round(t2p1,3)
    wa_central = wa_central * wn_to_omega
    a_coupling = a_coupling * wn_to_omega
    #use wm for w_RW for default
    if wm:
        w_RW = wm
    else:
        w_RW = w1
    #time interval to integrate over can be determined by the positions of the pulses
    #add in phase shift to account for different points at which coherence is emitted
    t = get_t(t21,t2p1,timestep)
    len_tprime = get_lengths(timestep)
    tprime = t[-len_tprime:]
    tprime_i = tprime[0]
    #tprime   = filter(lambda x: x >= tprime_i, t)
    E1  = 1./sigma_1* pulse(t,  E1amp,  sigma_1,  0.0,  w1 )
    E2  = 1./sigma_2* pulse(t,  E2amp,  sigma_2,  t21, -w2 )
    E2p = 1./sigma_2p*pulse(t,  E2pamp, sigma_2p, t2p1, w2p)
    #construct energy systems from which to sample
    zeta, zweight, dzeta = inhom_dist(inhom)
    sigma = sigma * wn_to_omega
    zeta = zeta * sigma
    dt=timestep
    out = np.array([0], dtype=np.complex64)
    #                 [gg1,      ag1,      ga,       aa,       2ag,       2aa,       gg2,      ag2]
    Gamma = np.array( [Gamma_gg, Gamma_ag, Gamma_ag, 
                       Gamma_aa, Gamma_aa, Gamma_2ag, 
                       Gamma_ag, Gamma_2aa ] )
    #compute NR signal as completely driven
    #first compute the solvent (CARS) and NR signal (E1*E2*E2p)--easy
    NR_sig = np.zeros((len(tprime)), dtype=np.complex64)
    index_NR = 0
    s_out = None
    if solvent:
        s_out = np.array([0], dtype=np.complex64)
        for z in range(len(ctet_w)):
            index_s = 0
            #calculate the output for each vibrational feature seperately
            s_rho_0 = np.zeros((3), dtype=np.complex)
            s_rho_0[0] = 1.0
            s_rho_i1 = s_rho_0
            s_rho_i2 = s_rho_0
            s_rho_emitted = np.zeros((len(tprime)), dtype=np.complex64)
            wvg = ctet_w[z]*wn_to_omega
            G=ctet_gamma[z]
            mu_vg=s_mu[z]
            sGamma= np.array([0,G,G])
            O1 = sOmega(E1, E2, E2p, t, wvg, mu_vg, w1first=True)
            O2 = sOmega(E1, E2, E2p, t, wvg, mu_vg, w1first=False)
            #need to implement Heun method here as well!
            for x in range(len(t)):
                tx = t[x]
                #calculate delta rho based on previous rho values
                delta_rho1 = -sGamma*s_rho_i1 + complex(0,0.5)*np.dot(O1[x], s_rho_i1)
                delta_rho2 = -sGamma*s_rho_i2 + complex(0,0.5)*np.dot(O2[x], s_rho_i2)
                s_rho_i1 = s_rho_i1 + delta_rho1*dt
                s_rho_i2 = s_rho_i2 + delta_rho2*dt
                eg1 = 0.5 * mu_vg * E2p[x] * rotor( wvg*tx)
                ev1 = 0.5 * mu_vg * E2p[x] * rotor(-wvg*tx)
                eg2 = 0.5 * mu_vg * E1[x]  * rotor( wvg*tx)
                ev2 = 0.5 * mu_vg * E1[x]  * rotor(-wvg*tx)
                if tprime_i <= tx:
                    #start storing ev and eg
                    #all components can be added together at same time since same rotating frame (namely, a non-rotating frame)
                    new_rho = (eg1*s_rho_i1[1] + ev1*s_rho_i1[2]
                        + eg2*s_rho_i2[1] + ev2*s_rho_i2[2])*rotor(-w_RW*tx)
                    s_rho_emitted[index_s] = s_rho_emitted[index_s] + new_rho
                    #print index_NR, s_rho_emitted[index_NR]
                    index_s = index_s + 1
            s_out = s_out + s_rho_emitted
        #s = np.abs(s_out).max()
        #return tprime, s_out / s
    else:
        pass
    for x in range(len(t)):
        tx = t[x]
        if tprime_i <= tx:
            NR_sig[index_NR] = E1[x]*E2[x]*E2p[x]*rotor(-w_RW*tx)
            #if signal gets too small, correct it to 0 instead of nan
            if np.isnan(NR_sig[index_NR]):
                NR_sig[index_NR] = 0.0
            index_NR = index_NR + 1
    NR_out = NR_sig #measure(NR_sig, tprime, w_RW, slitwidth, timestep, plot_dimension=plot_dimension)

    if inhom == 'linear':
        #initialize second distribution
        rho_emitted_prev = np.zeros((2,len(tprime)))
    for l in range(len(zeta)):
        #set initial conditions
        #generate state energies to use in simulation
        w = energy_dist(zeta, wa_central, a_coupling, l)
        #rho_0 = np.array( [1.0,   0,   0,   0,   0,   0,   0,   0  ], dtype=np.complex128)
        rho_0 = np.zeros((len(w)),dtype=np.complex128) #length should be 8
        #gg is populated
        rho_0[0] = 1.0
        #seperate state vectors keep track of both cases:  w1 first and w2p first
        rho_i1 = rho_0
        rho_i2 = rho_0
        rho_emitted = np.zeros((2,len(tprime)), dtype=np.complex64)
        #index to keep track of elements of rho_xx_emitted
        emitted_index = 0
        #integrate the signal at this pulse separation
        O1 = Omega(E1, E2, E2p, t, w, w1first=True)
        O2 = Omega(E1, E2, E2p, t, w, w1first=False)
        for k in range(len(t)):
            tk = t[k]
            #calculate delta rho based on previous rho values
            temp_delta_rho1 = -Gamma*rho_i1 + complex(0,0.5)*np.dot(O1[k], rho_i1)
            temp_delta_rho2 = -Gamma*rho_i2 + complex(0,0.5)*np.dot(O2[k], rho_i2)
            temp_rho_i1 = rho_i1 + temp_delta_rho1*dt
            temp_rho_i2 = rho_i2 + temp_delta_rho2*dt
            #make sure index is not out of range--iterative method fails on the last timestep
            try:
                delta_rho1 = -Gamma*temp_rho_i1 + complex(0,0.5)*np.dot(O1[k+1], temp_rho_i1)
                delta_rho2 = -Gamma*temp_rho_i2 + complex(0,0.5)*np.dot(O2[k+1], temp_rho_i2)
                #O1[k+1]
            except IndexError:
                rho_i1 = temp_rho_i1
                rho_i2 = temp_rho_i2
            else:
                rho_i1 = rho_i1 + dt/2 * (temp_delta_rho1 + delta_rho1)
                rho_i2 = rho_i2 + dt/2 * (temp_delta_rho2 + delta_rho2)
            #if we are close enough to final coherence emission, start storing these values
            if tprime_i <= tk:
                #transform all instances to the same rotating frame--namely, wm
                for s in range(len(rho_emitted)):
                    #each state vector has two elements that will emit (2aa and ag2)
                    #transform all instances to the same rotating frame--namely, wm
                    new_emitted = (rho_i1[s+6] + rho_i2[s+6])*rotor((w[s+6]-w_RW)*tk)
                    rho_emitted[s,emitted_index] = rho_emitted[s,emitted_index]  + new_emitted
                emitted_index = emitted_index+1
            else:
                pass
        if inhom == 'linear':
            #linearly spaced points, so we can use trapezoid rule
            #first make sure we have the two distributions needed to use the rule
            if l == 0:
                #first iteration does not work--store it
                rho_emitted_prev = rho_emitted
            else:
                out = out + dzeta/2 * (rho_emitted * zweight[l] + rho_emitted_prev * zweight[l-1])
                rho_emitted_prev = rho_emitted
        else:
            out = out + rho_emitted * zweight[l] * dzeta
    #now all coherences have been calculated, so convert to signal based on plot_dim
    if plot_dimension == 'time':
        var = tprime
        sig  = np.abs(out.sum(axis=0))
        sigmax = np.abs(sig).max()
        sig = sig #/ sigmax
        plt.plot(var, sig)
        #plt.plot(t,np.abs(E2))
        #plt.plot(t,np.abs(E1))
        #plt.plot(t,np.abs(E2p))
        plt.xlabel('t')
        plt.grid(b=True)
        return var, sig
    elif plot_dimension == 'freq':
        fft_tprime = np.fft.fftfreq(len(tprime), d=timestep)
        var        = np.fft.fftshift(fft_tprime) *2*np.pi / wn_to_omega
        #data freq is relative to w_RW (=wm)
        sig = np.array([var, out.sum(axis=0)])
        sig = np.abs(sig[1])**2
        sigmax  = np.abs(sig[1]).max()
        sig  = sig/sigmax
        plt.plot(var, sig)
        #sig_real = np.real(sig[1])
        #sig_imag = np.imag(sig[1])
        #plt.plot(var, sig_real)
        #plt.plot(var, sig_imag)
        plt.xlabel(r'$\nu\bar_m$')
        plt.grid(b=True)
        return var, sig
    elif plot_dimension == 'integrate':
        #sum up squared terms
        out = np.abs(out.sum(axis=0))**2
        return out.sum()*timestep
    else:
        #a quick hack to spit out the electric fields
        #gets rid of the ability to output ag and 2aa coherences seperately
        #because of concatenation, this doesn't work perfectly
        if scatter != None:
            l_tp = len(tprime) 
            if scatter == 'w1':
                scat_out = E1
            elif scatter == 'w2':
                scat_out = E2
            elif scatter == 'w2p':
                scat_out = E2p
            elif scatter == 'all':
                scat_out = E1+ E2 + E2p
            for i in range(len(scat_out)):
                scat_out[i] = scat_out[i]*rotor(-w_RW*t[i])
            out = out.sum(axis=0)
            return out, scat_out[0:l_tp], NR_out, s_out
        else:
            return out[0], out[1], NR_out, s_out
