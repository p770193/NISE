# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 20:24:20 2014

ctet hamiltonian for trive beams

@author: Dan
"""
from EOMPMA.misc import *

#--------------------------Solvent (ISR) Properties--------------------------
#define frequencies of states (wavenumbers), lifetimes(rad/fs), and dipoles
#frequencies from NIST

#              ag     bg     cg     dg     eg
ctet_w      = [217.,  313.,  459.,  762.,  790.]
ctet_gamma  = [0.0012,0.0012,0.0009,0.0012,0.0012]
s_mu        = [0.6, 0.5, 1, 0.22, 0.22]#[0.45,  0.45,  1,     0.22, 0.22]

def sOmega(E1, E2, E3, t, wvg, mu_vg, w1first = True):
    """
    creates the solvent coupling array given the input e-fields values for a specific time, t
    w1first selects whether w1 or w2p is the first interacting positive field
    accounts for CARS interactions
    """
    if w1first==True:
        first  = E1
    else:
        first  = E3
    O = np.zeros((len(t),3,3),dtype=np.complex64)
    # from gg
    O[:,1,0] = 0.5*mu_vg * E2 * first * rotor(-wvg*t)
    O[:,2,0] = 0.5*mu_vg * E2 * first * rotor( wvg*t)
    return O

def gen_H_ISR():
    # generate the Hamiltonian for the Raman modes
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
        O1 = sOmega(E1, E2, E3, t, wvg, mu_vg, w1first=True)
        O2 = sOmega(E1, E2, E3, t, wvg, mu_vg, w1first=False)
        # need to implement Heun method here as well!
        for x in range(len(t)):
            tx = t[x]
            # calculate delta rho based on previous rho values
            delta_rho1 = -sGamma*s_rho_i1 + complex(0,0.5)*np.dot(O1[x], s_rho_i1)
            delta_rho2 = -sGamma*s_rho_i2 + complex(0,0.5)*np.dot(O2[x], s_rho_i2)
            s_rho_i1 = s_rho_i1 + delta_rho1*dt
            s_rho_i2 = s_rho_i2 + delta_rho2*dt
            # eg and ev are NOT in their rotating frames
            # emitted coherences derived from appendix c of notes
            eg1 = 0.5 * mu_vg * E3[x] * rotor( wvg*tx)
            ev1 = 0.5 * mu_vg * E3[x] * rotor(-wvg*tx)
            eg2 = 0.5 * mu_vg * E1[x]  * rotor( wvg*tx)
            ev2 = 0.5 * mu_vg * E1[x]  * rotor(-wvg*tx)
            if tprime_i <= tx:
                # start storing ev and eg
                # all components can be added together at same time since 
                # same rotating frame (namely, a non-rotating frame)
                new_rho = (eg1*s_rho_i1[1] + ev1*s_rho_i1[2]
                    + eg2*s_rho_i2[1] + ev2*s_rho_i2[2])*rotor(-w_RW*tx)
                s_rho_emitted[index_s] = s_rho_emitted[index_s] + new_rho
                index_s = index_s + 1
        s_out = s_out + s_rho_emitted
    #s = np.abs(s_out).max()
    #return tprime, s_out / s
