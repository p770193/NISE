# -*- coding: utf-8 -*-
"""
Created on Mon Jul 28 17:25:28 2014

@author: Dan

each instance of running this depends on a few initial conditions that have to 
be specified:
    out_group
    rho_0
    wa_central
    a_coupling
    gamma
    dipoles

so create a class where all these can describe the specific instance
"""

from NISE.lib.misc import *

def gen_w_0(wa_central, a_coupling):
    # convert nice system parameters into system vector indeces
    w_ag = wa_central
    w_2aa = w_ag - a_coupling
    w_2ag = 2*w_ag - a_coupling
    w_gg = 0.
    w_aa = w_gg
    return np.array( [w_gg, w_ag, -w_ag, w_gg, w_aa, w_2ag, w_ag, w_2aa] )

def gen_Gamma_0(Gamma_ag, Gamma_aa, Gamma_2ag, Gamma_2aa):
    # same as gen_w_0, but for dephasing/relaxation times
    Gamma = np.array( [0.0, Gamma_ag, Gamma_ag, 
                       Gamma_aa, Gamma_aa, Gamma_2ag, 
                       Gamma_ag, Gamma_2aa ] )
    return Gamma

class Omega:
    # record the propagator module used to evolve this hamiltonian
    propagator = 'rk'
    # phase cycling is not valuable in this hamiltonian
    pc = False
    # all attributes should have good initial guesses for parameters
    dm_vector = [
        '00><00',
        '11><00', '1m1><00', '00><11', '00><1m1', 
        '11><11', '00><00', '1m1><11', '11><1m1', '1m1><1m1',
        '22><00', '20><00', '2m2><00',
        '22><11', '2m2><1m1', '20><11', '20><1m1', '11><00', '1m1><00'
    ]
    #out_group = [[6,7]]#,[7]]
    out_group = [[6,7]]
    #--------------------------Oscillator Properties--------------------------
    rho_0 = np.zeros((8), dtype=np.complex64)
    rho_0[0] = 1.
    #1S exciton central position
    wa_central = 7000.
    #exciton-exciton coupling
    a_coupling = 0.0
    #dephasing times, 1/fs
    Gamma_1x0x  = 1./60
    Gamma_1x1x  = 0.
    Gamma_2x0x  = 1./60
    Gamma_2x1x  = 1./60
    Gamma_axbx  = 1./500
    #transition dipoles (a.u.)
    mu_ag =  1.0
    mu_2aa = 1.0 * mu_ag #HO approx (1.414) vs. uncorr. electron approx. (1.)
    #--------------------------Recorded attributes--------------------------
    out_vars = ['dm_vector', 'out_group', 'rho_0', 'mu_ag', 'mu_2aa', 
                'wa_central', 'a_coupling', 'pc', 'propagator', 'Gamma',
                'rlphase']
    #--------------------------Methods--------------------------
    def __init__(self, rlphase=[0,0,0],**kwargs):
        # inherit all class attributes unless kwargs has them; then use those 
        # values.  if kwargs is not an Omega attribute, it gets ignored
        # careful: don't redefine instance methods as class methods!
        for key, value in kwargs.items():
            if key in Omega.__dict__.keys(): 
                setattr(self, key, value)
            else:
                print 'did not recognize attribute {0}.  No assignment made'.format(key)
        # with this set, initialize parameter vectors
        self.rlphase=rlphase
        self.w_0 = gen_w_0(self.wa_central, self.a_coupling)
        self.Gamma = gen_Gamma_0(self.Gamma_ag, self.Gamma_aa, self.Gamma_2ag, self.Gamma_2aa)

    def o(self, efields, t, wl):
        # combine the two pulse permutations to produce one output array
        E1, E2, E3 = efields[0:3]
    
        out1 = self._gen_matrix(E1, E2, E3, t, wl, w1first = True)
        out2 = self._gen_matrix(E1, E2, E3, t, wl, w1first = False)

        return np.array([out1, out2], dtype=np.complex64)
    
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
        w2aa = wl[7]
        
        mu_ag = self.mu_ag
        mu_2aa = self.mu_2aa
    
        E1r = E1
        E2r = E2
        E3r = E3
        # decompose e-field into r and l components according to rlphase    
        rl_phase_factor = np.exp(1j*self.rlphase)
        E1l = E1 * rl_phase_factor[0]
        E2l = E2 * rl_phase_factor[1]
        E3l = E3 * rl_phase_factor[2]
    
        if w1first==True:
            firstr  = E1r
            firstl  = E1l
            secondr = E3r
            secondl = E3l
        else:
            firstr  = E3r
            firstl  = E3l
            secondr = E1r
            secondl = E1l
        O = np.zeros((len(t), len(wl), len(wl)), dtype=np.complex64)
        # from gg1
        O[:,1,0] =  mu_ag  * first  * rotor(-wag*t)
        O[:,2,0] = -mu_ag  * E2     * rotor(wag*t)
        # from ag1
        O[:,4,1] = -mu_ag  * E2     * rotor(wag*t)
        O[:,5,1] =  mu_2aa * second * rotor(-w2aa*t)
        O[:,3,1] =  mu_ag  * E2     * rotor(wag*t)
        # from ga
        O[:,4,2] =  mu_ag  * first  * rotor(-wag*t)
        O[:,3,2] = -mu_ag  * first  * rotor(-wag*t)
        # from gg2
        O[:,6,3] =  mu_ag  * second * rotor(-wag*t)      * mu_ag
        # from aa
        O[:,7,4] =  mu_2aa * second * rotor(-w2aa*t)     * mu_2aa
        O[:,6,4] = -mu_ag  * second * rotor(-wag*t)      * mu_ag
        # from 2ag
        O[:,7,5] = -mu_ag  * E2     * rotor(wag*t)       * mu_2aa
        O[:,6,5] =  mu_2aa * E2     * rotor(w2aa*t)      * mu_ag
        
        return O
    
    def ws(self, inhom_object):
        """
        creates the correspondence of oscillator energies to the state vector
        contains instructions for how energies change as subsets are changed
        """
        z = inhom_object.zeta
        
        wg = 0.0 + 0*z
        wa = z + self.wa_central
        w2a = 2*wa - self.a_coupling
    
        w_ag = wa - wg
        w_aa = wa - wa
        w_gg = wg - wg
        w_2ag = w2a - wg
        w_2aa = w2a - wa
        #array aggregates all frequencies to match state vectors
        
        w     = np.array( [w_gg, w_ag, -w_ag, w_gg, w_aa, w_2ag, w_ag, w_2aa] )
        
        return w
