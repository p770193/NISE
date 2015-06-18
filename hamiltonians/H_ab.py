# -*- coding: utf-8 -*-
"""
Created on Tue Aug 05 23:43:16 2014

@author: Dan
"""

from NISE.lib.misc import *

def gen_w_0(wa_central, a_coupling,
            wb_central, b_coupling, 
            ab_coupling):
    # convert nice system parameters into system vector indeces
    w_ag = wa_central
    w_bg = wb_central
    
    w_ba = w_bg - w_ag
    w_ab = -w_ba

    w_gg = 0.
    w_aa = w_gg
    w_bb = w_gg

    w_2ag = 2*w_ag - a_coupling
    w_2bg = 2*w_bg - b_coupling
    w_cg = w_ag + w_bg - ab_coupling
    
    w_2aa = w_ag - a_coupling
    w_2bb = w_bg - b_coupling

    w_ca = w_cg - w_ag
    w_cb = w_cg - w_bg

    return np.array( [w_gg, 
                      w_ag, w_bg, -w_ag, -w_bg,
                      w_aa, w_bb, w_ab, w_ba, 
                      w_2ag, w_2bg, w_cg,
                      w_ca, w_cb, 
                      w_2aa, w_2bb, w_ag, w_bg
                      ] )

def gen_Gamma_0(Gamma_ag, Gamma_bg,
                Gamma_aa, Gamma_bb, Gamma_ab,
                Gamma_dqc,
                Gamma_2aa, Gamma_2bb,
                Gamma_ca, Gamma_cb):
    # same as gen_w_0, but for dephasing/relaxation times
    Gamma = np.array( [0.0, Gamma_ag, Gamma_bg, Gamma_ag, Gamma_bg,
                       Gamma_aa, Gamma_bb, Gamma_ab, Gamma_ab, 
                       Gamma_dqc, Gamma_dqc, Gamma_dqc,
                       Gamma_ca, Gamma_cb, Gamma_2aa, Gamma_2bb,
                       Gamma_ag, Gamma_bg] )
    return Gamma

class Omega:
    # record the propagator module used to evolve this hamiltonian
    propagator = 'rk'
    pc = False
    kT = 210. # thermal energy parameter, in cm-1
    # 18 elements
    dm_vector = ['gg1',
                 'ag', 'bg', 'ga', 'gb', 
                 'aa', 'bb', 'ab', 'ba', 
                 '2ag', '2bg', 'a+b,g',
                 'a+b,a', 'a+b,b', '2a,a', '2b,b', 
                 'ag2','bg2']
    out_group = [[12,13,14,15],[16,17]]
    #--------------------------Oscillator Properties--------------------------
    rho_0 = np.zeros((len(dm_vector)), dtype=np.complex64)
    rho_0[0] = 1.
    # state central position
    wa_central = 6500.
    wb_central = 7500.
    # exciton-exciton coupling
    a_coupling = 200.0
    b_coupling = 200.0
    ab_coupling = 250.0
    # dephasing times, 1/fs
    Gamma_ag  = 1./25
    Gamma_bg  = 1./25
    Gamma_aa  = 0.0 #1./2000.
    Gamma_bb  = 1./100
    Gamma_ab  = 1./15
    Gamma_dqc = 1./15
    Gamma_2aa = 1./25
    Gamma_2bb = Gamma_2aa
    Gamma_ca  = Gamma_2aa
    Gamma_cb  = Gamma_2aa

    #transition dipoles (a.u.)
    mu_ag =  1.0
    mu_bg =  1.0

    mu_2bb = 1.0 * mu_bg
    mu_2aa = 1.0 * mu_ag 
    
    mu_ca = 1.0 * mu_bg
    mu_cb = 1.0 * mu_ag
    #--------------------------Recorded attributes--------------------------
    out_vars = ['dm_vector', 'out_group', 'rho_0', 
                'mu_ag', 'mu_2aa', 'mu_bg', 'mu_2bb',
                'wa_central', 'wb_central', 'a_coupling', 'b_coupling', 
                'ab_coupling', 'Gamma', 'w_0',
                'pc', 'propagator']
    #--------------------------Methods--------------------------

    def __init__(self, **kwargs):
        # inherit all class attributes unless kwargs has them; then use those 
        # values.  if kwargs is not an Omega attribute, it gets ignored
        # careful: don't redefine instance methods as class methods!
        for key, value in kwargs.items():
            if key in Omega.__dict__.keys(): 
                setattr(self, key, value)
            else:
                print 'did not recognize attribute {0}.  No assignment made'.format(key)
        # with this set, initialize parameter vectors
        # w_0 is never actually used for computations; only for reporting back...
        self.w_0 = gen_w_0(self.wa_central, self.a_coupling,
                           self.wb_central, self.b_coupling,
                           self.ab_coupling)
        self.Gamma = gen_Gamma_0(self.Gamma_ag, self.Gamma_bg, 
                                 self.Gamma_aa, self.Gamma_bb, self.Gamma_ab,
                                 self.Gamma_dqc, 
                                 self.Gamma_2aa, self.Gamma_2bb, 
                                 self.Gamma_ca, self.Gamma_cb)

    def o(self, efields, t, wl):
        # combine the two pulse permutations to produce one output array
        E1, E2, E3 = efields[0:3]
        
        # NOTICE:: w1first is false for both terms
        out1 = self._gen_matrix(E1, E2, E3, t, wl, w1first = False)
        out2 = self._gen_matrix(E1, E2, E3, t, wl, w1first = True)
        return np.array([out1, out2], dtype=np.complex64)
    
    # to get a dummy matrix that shows connectivities, run
    # _gen_matrix(1,1,1,0,np.zeros((len(dm_vector))))
    def _gen_matrix(self, E1, E2, E3, t, wl, w1first = True):
        """
        creates the coupling array given the input e-fields values for a specific time, t
        w1first selects whether w1 or w2p is the first interacting positive field
        
        Currently neglecting pathways where w2 and w3 require different frequencies
        (all TRIVE space, or DOVE on diagonal)
        
        Matrix formulated such that dephasing/relaxation is accounted for 
        outside of the matrix
        """
        wag  = wl[1]
        wbg  = wl[2]
        wca  = wl[12]
        wcb  = wl[13]
        w2aa = wl[14]
        w2bb = wl[15]
        
        mu_ag = self.mu_ag
        mu_2aa = self.mu_2aa
        mu_bg = self.mu_bg
        mu_2bb = self.mu_2bb
        mu_ca = self.mu_ca
        mu_cb = self.mu_cb
    
        if w1first==True:
            first  = E1
            second = E3
        else:
            first  = E3
            second = E1
        O = np.zeros((len(t), len(wl), len(wl)), dtype=np.complex64)

        # from gg1
        O[:,1,0] =  mu_ag  * first  * rotor(-wag*t)     # ag
        O[:,2,0] =  mu_bg  * first  * rotor(-wbg*t)     # bg
        O[:,3,0] = -mu_ag  * E2     * rotor(wag*t)      # ga
        O[:,4,0] = -mu_bg  * E2     * rotor(wbg*t)      # gb
        # from ag1
        O[:,5,1]  = -mu_ag  * E2     * rotor(wag*t)     # aa
        O[:,7,1]  = -mu_bg  * E2     * rotor(wbg*t)     # ab
        O[:,9,1]  =  mu_ag  * second * rotor(-w2aa*t)   # 2ag
        O[:,11,1] =  mu_bg  * second * rotor(-wca*t)    # cg
        # from bg1
        O[:,6,2]  = -mu_bg  * E2     * rotor(wbg*t)     # bb
        O[:,8,2]  = -mu_ag  * E2     * rotor(wag*t)     # ba
        O[:,10,2] =  mu_bg  * second * rotor(-w2bb*t)   # 2bg
        O[:,11,2] =  mu_ag  * second * rotor(-wcb*t)    # cg
        # from ga
        O[:,5,3] =  mu_ag  * first  * rotor(-wag*t)     # aa
        O[:,8,3] =  mu_bg  * first  * rotor(-wbg*t)     # ba
        # from gb
        O[:,6,4] =  mu_bg  * first  * rotor(-wbg*t)     # bb
        O[:,7,4] =  mu_ag  * first  * rotor(-wag*t)     # ab
        # from aa (-gg)
        O[:,14,5] =  mu_2aa * second * rotor(-w2aa*t) * mu_2aa      # 2aa
        O[:,12,5] =  mu_ca  * second * rotor(-wca*t)  * mu_ca       # ca
        O[:,17,5] = -mu_bg  * second * rotor(-wbg*t)  * mu_bg       # gg -> bg
        #   twice the weight because of gamma and alpha pathways
        O[:,16,5] = -mu_ag  * second * rotor(-wag*t)  * mu_ag * 2   # ag
        # from bb (-gg)
        O[:,15,6] =  mu_2bb * second * rotor(-w2bb*t) * mu_2bb      # 2bb
        O[:,13,6] =  mu_cb  * second * rotor(-wcb*t)  * mu_cb       # cb
        O[:,16,6] = -mu_ag  * second * rotor(-wag*t)  * mu_ag       # gg -> ag
        #   twice the weight because of gamma and alpha pathways
        O[:,17,6] = -mu_bg  * second * rotor(-wbg*t)  * mu_bg * 2   # bg
        # from ab
        O[:,13,7] =  mu_ca  * second * rotor(-wca*t)  * mu_cb       # cb
        O[:,16,7] = -mu_bg  * second * rotor(-wbg*t)  * mu_ag       # ag
        # from ba
        O[:,12,8] =  mu_cb  * second * rotor(-wcb*t)  * mu_ca       # ca
        O[:,17,8] = -mu_ag  * second * rotor(-wag*t)  * mu_bg       # bg
        # from 2ag
        O[:,14,9] = -mu_ag  * E2 * rotor(wag*t)  * mu_2aa       # 2aa
        O[:,16,9] =  mu_2aa * E2 * rotor(w2aa*t) * mu_ag        # ag
        # from 2bg
        O[:,15,10] = -mu_bg  * E2 * rotor(wbg*t)  * mu_2bb      # 2bb
        O[:,17,10] =  mu_2bb * E2 * rotor(w2bb*t) * mu_bg       # bg
        # from cg
        O[:,12,11] = -mu_ag  * E2 * rotor(wag*t)  * mu_ca       # ca
        O[:,13,11] = -mu_bg  * E2 * rotor(wbg*t)  * mu_cb       # cb
        O[:,16,11] =  mu_ca  * E2 * rotor(wca*t)  * mu_ag       # bg
        O[:,17,11] =  mu_cb  * E2 * rotor(wcb*t)  * mu_bg       # ag

        # make complex according to Liouville Equation
        O *= complex(0,0.5)
        for i in range(O.shape[-1]):
            O[:,i,i] = -self.Gamma[i]
        # include coherence/population transfers here:
        # bb --> aa
        ratio = np.exp((self.wa_central - self.wb_central) / self.kT)
        O[:,5,5] += -self.Gamma[6] * ratio
        O[:,6,5] = self.Gamma[6] * ratio
        #O[:,6,6] = -self.Gamma[6] 
        O[:,5,6] = self.Gamma[6]

        return O

    def ws(self, inhom_object):
        """
        creates the correspondence of oscillator energies to the state vector
        contains instructions for how energies change as subsets are changed
        """
        z = inhom_object.zeta
        # define bohr frequencies
        wg = 0.0 + 0*z
        wa = z + self.wa_central
        wb = z + self.wb_central
        w2a = 2*wa - self.a_coupling
        w2b = 2*wb - self.b_coupling
        wc = wa + wb - self.ab_coupling
        # define coherence frequencies
        w_gg = wg - wg
        w_ag = wa - wg
        w_bg = wb - wg
        w_aa = wa - wa
        w_bb = w_aa
        w_ab = wa - wb
        w_2ag = w2a - wg
        w_2bg = w2b - wg
        w_cg = wc - wg
        w_2aa = w2a - wa
        w_2bb = w2b - wb
        w_ca = wc - wa
        w_cb = wc - wb
        #array aggregates all frequencies to match state vectors
        w     = np.array( [w_gg, 
                           w_ag, w_bg, -w_ag, -w_bg,
                           w_aa, w_bb, w_ab, -w_ab,
                           w_2ag, w_2bg, w_cg,
                           w_ca, w_cb, w_2aa, w_2bb,
                           w_ag, w_bg] )
        return w

"""
Notes for generating a transparent colorbar; use later. From:
http://matplotlib.1069221.n5.nabble.com/colormap-from-blue-to-transparent-to-red-td26440.html

>>> theCM = cm.get_cmap('somemap')
>>> theCM._init()

Then, you need to fill in the alpha values for yourself.  For this, I am going to use a triangle function:

>>> alphas = np.abs(np.linspace(-1.0, 1.0, theCM.N))
>>> theCM._lut[:,-1] = alphas

"""