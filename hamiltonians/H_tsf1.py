# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 14:41:13 2015

tsf demo hamiltonian
only calculates the DQC amplitude:  leave the electronic excitation to 
a CW routine

@author: Dan
"""

from NISE.lib.misc import *

def gen_w_0(wa_central, a_coupling,
            wb_central, b_coupling, 
            ab_coupling):
    # convert nice system parameters into system vector indeces
    w_gg = 0.
    w_ag = wa_central
    w_bg = wb_central
    w_2ag = 2*w_ag - a_coupling
    w_2bg = 2*w_bg - b_coupling
    w_cg = w_ag + w_bg - ab_coupling
    return np.array( [w_gg, 
                      w_ag, w_bg,
                      w_2ag, w_2bg, w_cg
                      ] )

def gen_Gamma_0(Gamma_ag, Gamma_bg,
                Gamma_2ag, Gamma_2bg,
                Gamma_cg):
    # same as gen_w_0, but for dephasing/relaxation times
    Gamma = np.array( [0.0, 
                       Gamma_ag, Gamma_bg,
                       Gamma_2ag, Gamma_2bg, Gamma_cg] )
    return Gamma

class Omega:
    # record the propagator module used to evolve this hamiltonian
    propagator = 'rk'  # dictates the function we send this matrix to
    pc = False # hamiltonian uses phase-matched liouville pathways, not phase cycling
    # 18 elements
    dm_vector = ['g,g',
                 'a,g', 'b,g', 
                 '2a,g', '2b,g', 'ab,g']
    out_group = [[3],[4],[5]]
    #--------------------------Oscillator Properties--------------------------
    rho_0 = np.zeros((len(dm_vector)), dtype=np.complex64)
    rho_0[0] = 1.
    # state central position
    wa_central = 1800.  # cm-1
    wb_central = 2000.  # cm-1
    # exciton-exciton coupling
    a_coupling = 20.0
    b_coupling  = 20.0
    ab_coupling = 20.0  # cm-1
    # dephasing times, 1/fs
    Gamma_ag  = 1./1000
    Gamma_bg  = 1./1000
    Gamma_2ag = Gamma_ag
    Gamma_2bg = Gamma_bg
    Gamma_cg = 1./1000
    
    #transition dipoles (a.u.)
    mu_ag =  1.0
    mu_bg =  0.7

    mu_2bb = np.sqrt(2) * mu_bg
    mu_2aa = np.sqrt(2) * mu_ag 
    
    mu_ca = mu_bg
    mu_cb = mu_ag
    #--------------------------Recorded attributes--------------------------
    out_vars = ['dm_vector', 'out_group', 'rho_0',
                'mu_ag', 'mu_bg', 'mu_2aa', 'mu_2bb',
                'mu_ca', 'mu_cb',
                'wa_central', 'wb_central',
                'a_coupling', 'b_coupling', 
                'ab_coupling',
                'pc', 'propagator', 
                'Gamma_ag', 'Gamma_bg', 
                'Gamma_2ag', 'Gamma_2bg', 'Gamma_cg'
                ]

    #--------------------------Methods--------------------------

    def __init__(self, **kwargs):
        # inherit all class attributes unless kwargs has them; then use those 
        # values.  if kwargs is not an Omega attribute, it gets ignored
        # careful: don't redefine instance methods as class methods!
        for key, value in kwargs.items():
            if key in Omega.__dict__.keys(): 
                setattr(self, key, value)
            else:
                print('did not recognize attribute {0}.  No assignment made'.format(key))
        # with this set, initialize parameter vectors
        # w_0 is never actually used for computations; only for reporting back...
        self.w_0 = gen_w_0(self.wa_central, self.a_coupling,
                           self.wb_central, self.b_coupling,
                           self.ab_coupling
                           )
        self.Gamma = gen_Gamma_0(self.Gamma_ag, self.Gamma_bg, 
                                 self.Gamma_2ag, self.Gamma_2bg, 
                                 self.Gamma_cg
                                 )

    def o(self, efields, t, wl):
        # combine the two pulse permutations to produce one output array
        E1, E2 = efields[0:3]
        
        out1 = self._gen_matrix(E1, E2, t, wl, E1first = False)
        out2 = self._gen_matrix(E1, E2, t, wl, E1first = True)

        return np.array([out1, out2], dtype=np.complex64)
    
    # to get a dummy matrix that shows connectivities, run
    # _gen_matrix(1,1,1,0,np.zeros((len(dm_vector))))
    def _gen_matrix(self, E1, E2, t, wl, E1first = True):
        """
        creates the coupling array given the input e-fields values for a specific time, t
        w1first selects whether w1 or w2p is the first interacting positive field
        """
        wag, wbg, w2ag, w2bg, wcg = wl[1:]
        w2aa = w2ag - wag
        w2bb = w2bg - wbg
        wca = wcg - wag
        wcb = wcg - wbg
        
        mu_ag = self.mu_ag
        mu_2aa = self.mu_2aa
        mu_bg = self.mu_bg
        mu_2bb = self.mu_2bb
        mu_ca = self.mu_ca
        mu_cb = self.mu_cb
        
        if E1first==True:
            first  = E1
            second = E2
        else:
            first  = E2
            second = E1
        O = np.zeros((len(t), len(wl), len(wl)), dtype=np.complex64)
        # from gg
        O[:,1,0] =  mu_ag  * first  * rotor(-wag*t)
        O[:,2,0] =  mu_bg  * first  * rotor(-wbg*t)
        # from ag
        O[:,3,1] =  mu_2aa  * second * rotor(-w2aa*t)
        O[:,5,1] =  mu_ca   * second * rotor(-wca*t)
        # from bg
        O[:,4,2] =  mu_2bb  * second * rotor(-w2bb*t)
        O[:,5,2] =  mu_cb   * second * rotor(-wcb*t)

        # make complex according to Liouville Equation
        O *= complex(0,0.5)
        for i in range(O.shape[-1]):
            O[:,i,i] = -self.Gamma[i]

        return O

    def ws(self, inhom_object):
        """
        creates the correspondence of oscillator energies to the state vector
        contains instructions for how energies change as subsets are changed
        """
        z = inhom_object.zeta
        # define bohr frequencies of eigenstates--currently a and b are 
        # correlated 
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
        w_2ag = w2a - wg
        w_2bg = w2b - wg
        w_cg = wc - wg
        # array aggregates all frequencies to match state vectors
        w     = np.array( [w_gg, 
                           w_ag, w_bg,
                           w_2ag, w_2bg, w_cg
                           ])
        return w
