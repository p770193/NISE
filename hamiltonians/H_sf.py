# -*- coding: utf-8 -*-
"""
Created on Wed Nov 05 03:18:27 2014

base hamiltonian for degenerate states (D parameter)

@author: Dan
each instance of running this depends on a few initial conditions that have to 
be specified:
    D
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
    dm_vector = ['gg1','ag','ga','gg2','aa','2ag','ag2','2aa']
    #out_group = [[6,7]]#,[7]]
    out_group = [[6,7]]
    #--------------------------Oscillator Properties--------------------------
    rho_0 = np.zeros((8), dtype=np.complex64)
    rho_0[0] = 1.
    D = 4
    #1S exciton central position
    wa_central = 7000.
    #exciton-exciton coupling
    a_coupling = 300.0
    #dephasing times, 1/fs
    Gamma_ag  = 1./25
    Gamma_aa  = 0.0 #1./2000.
    Gamma_2aa = 1./25
    Gamma_2ag = 1./15
    #transition dipoles (a.u.)
    mu_ag =  1.0
    mu_2aa = 1.0 * mu_ag #HO approx (1.414) vs. uncorr. electron approx. (1.)
    #--------------------------Recorded attributes--------------------------
    out_vars = ['dm_vector', 'out_group', 'rho_0', 'mu_ag', 'mu_2aa', 
                'wa_central', 'a_coupling', 'pc', 'propagator', 'D']
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
        self.w_0 = gen_w_0(self.wa_central, self.a_coupling)
        self.Gamma = gen_Gamma_0(self.Gamma_ag, self.Gamma_aa, self.Gamma_2ag, 
                                 self.Gamma_2aa)

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
        D = self.D
        D01 = D
        D12 = max(D-1, 0)
        D10 = 1.
        D21 = 2.
    
        if w1first==True:
            first  = E1
            second = E3
        else:
            first  = E3
            second = E1
        O = np.zeros((len(t), len(wl), len(wl)), dtype=np.complex64)
        # from gg1
        O[:,1,0] =  D01 * mu_ag  * first  * rotor(-wag*t)
        O[:,2,0] = -D01 * mu_ag  * E2     * rotor(wag*t)
        # from ag1
        O[:,4,1] = -1    * mu_ag  * E2     * rotor(wag*t)
        O[:,5,1] =  D12  * mu_2aa * second * rotor(-w2aa*t)
        O[:,3,1] =  1    * mu_ag  * E2     * rotor(wag*t)
        # from ga
        O[:,4,2] =  1 * mu_ag  * first  * rotor(-wag*t)
        O[:,3,2] = -1 * mu_ag  * first  * rotor(-wag*t)
        # from gg2
        O[:,6,3] =  D01 * mu_ag  * second * rotor(-wag*t)      * mu_ag
        # from aa
        O[:,7,4] =  D12 * mu_2aa * second * rotor(-w2aa*t)     * mu_2aa
        O[:,6,4] = -D10 * mu_ag  * second * rotor(-wag*t)      * mu_ag
        # from 2ag
        O[:,7,5] = -D21 * mu_ag  * E2     * rotor(wag*t)       * mu_2aa
        O[:,6,5] =  D21 * mu_2aa * E2     * rotor(w2aa*t)      * mu_ag
        
        # make complex according to Liouville Equation
        O *= complex(0,0.5)
        
        # include coherence decay rates:
        for i in range(O.shape[-1]):
            O[:,i,i] = -self.Gamma[i]

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
