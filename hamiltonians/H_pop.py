# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 17:29:11 2016

Derived from H0:  to use for refractive index incorporation.
Stop propagation after the populations are formed--use CW modifications after.
Need to distinguish between a 12,22' pop gratings, so use two different elements

---------------------
each instance of running this depends on a few initial conditions that have to 
be specified:
    out_group
    rho_0
    wa_central
    gamma
    dipole
"""

from NISE.lib.misc import *

def gen_w_0(wa_central):
    # convert nice system parameters into system vector indeces
    w_ag = wa_central
    w_gg = 0.
    w_aa = w_gg
    return np.array( [w_gg, w_ag, -w_ag, w_aa, w_aa] )

def gen_Gamma_0(tau_ag, tau_aa):
    # same as gen_w_0, but for dephasing/relaxation times
    tau = np.array( [np.inf, tau_ag, tau_ag, 
                     tau_aa, tau_aa] )
    Gamma = 1/tau
    return Gamma

class Omega:
    # record the propagator module used to evolve this hamiltonian
    propagator = 'rk'
    # phase cycling is not valuable in this hamiltonian
    pc = False
    # all attributes should have good initial guesses for parameters
    dm_vector = ['gg1','ag','ga','aa22', 'aa12']
    out_group = [[3],[4]] # use this to separate 12 from 22 population
    #--------------------------Oscillator Properties--------------------------
    rho_0 = np.zeros((len(dm_vector)), dtype=np.complex64)
    rho_0[0] = 1.
    # 1S exciton central position
    wa_central = 7000.
    # dephasing times, fs
    tau_ag  = 50.
    tau_aa  = np.inf #1./2000.
    mu_ag =  1.0
    # TOs sets which time-ordered pathways to include (1-6 for TrEE)
    # defaults to include all time-orderings included
    TOs = np.array([1,3,5,6]) 
    #--------------------------Recorded attributes--------------------------
    out_vars = ['dm_vector', 'out_group', 'rho_0', 'mu_ag', 
                'tau_ag', 'tau_aa',
                'wa_central', 'pc', 'propagator', 
                'TOs']
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
        self.w_0 = gen_w_0(self.wa_central)
        self.Gamma = gen_Gamma_0(self.tau_ag, self.tau_aa)

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
        
        mu_ag = self.mu_ag
    
        if w1first==True:
            first  = E1
            #second = E3
        else:
            first  = E3
            #second = E1

        O = np.zeros((len(t), len(wl), len(wl)), dtype=np.complex64)
        # from gg1
        O[:,1,0] =  mu_ag  * first  * rotor(-wag*t)
        O[:,2,0] = -mu_ag  * E2     * rotor(wag*t)
        # from ag1
        #   to pop
        if w1first and 1 in self.TOs:
            O[:,4,1] =  -mu_ag  * E2     * rotor(wag*t)
        if not w1first and 6 in self.TOs:
            O[:,3,1] =  -mu_ag  * E2     * rotor(wag*t)
        # from ga
        if w1first and 3 in self.TOs:
            O[:,4,2] =  mu_ag  * first  * rotor(-wag*t)
        if not w1first and 5 in self.TOs:
            O[:,3,2] =  mu_ag  * first  * rotor(-wag*t)
        
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
    
        w_ag = wa - wg
        w_aa = wa - wa
        w_gg = wg - wg
        #array aggregates all frequencies to match state vectors
        w     = np.array( [w_gg, w_ag, -w_ag, w_aa, w_aa] )
        return w
