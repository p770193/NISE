"""
Created on Wed Jun 13 17:06:55 2012
@author: Dan
Using Domcke approach to numerically integrate TRIVE signal of diagonal and 
near-diagonal features
"""

from .misc import *

def rk_1():
    # a function to use for working with raman signals, it does not integrate 
    # the whole matrix
    return

def rk(t, efields, iprime, inhom_object, omega_object):
    """
    Evolves the hamiltonian in time, given the three 1d array of E-field points 
    for each E-field and the hamiltonian

    H[pulse_permutations, t, ij_to, ij_from] is a rank 4 array of the hamiltonian 2D 
    matrix for each time point and necessary pulse permutations

    Uses Runge-Kutta method to integrate 
        --> assumes eqaully spaced t-points (t_int is spacing)

    Unlike earlier implementations, delays, phase, etc., must be incorporated
    externally to the function.  Rotating wave may still be specified.
    """
    # can only call on iprime and t after efield_object.E is called
    dt = np.abs(t[1]-t[0])
    # extract inhomogeneity parameters
    zeta = inhom_object.zeta
    zweight = inhom_object.zweight
    dzeta = inhom_object.dzeta
    # extract attributes of the system
    rho_0 = omega_object.rho_0
    #Gamma = omega_object.Gamma
    ws = omega_object.ws(inhom_object)
    # convert to rad / fs
    ws = ws * wn_to_omega
    # define out_groups
    out_group = omega_object.out_group
    rho_emitted = np.zeros((len(out_group),iprime), dtype=np.complex64)
    rho_emitted_prev = np.zeros((len(out_group),iprime), dtype=np.complex64)
    for l in range(len(zeta)):
        # identify the frequencies of the states?
        if len(ws.shape) == 1:
            w_l = ws
        else:
            w_l = ws[:,l]
        # H has 4 dimensions: permutations, time, and the 2 matrix dimensions
        H = omega_object.o(efields, t, w_l)
        # tile initial conditions for each permutation
        rho_i = np.tile(rho_0, (H.shape[0], 1))
        temp_rho_i = rho_i.copy()
        temp_delta_rho = rho_i.copy()
        delta_rho = rho_i.copy()
        # index to keep track of elements of rho_emitted
        emitted_index = 0
        out = rho_emitted.copy()
        for k in range(len(t)):
            tk = t[k]
            # now sum over p equivalent pulse permutations (e.g. TRIVE near-
            # diagonal has 2 permutations)
            for p in range(H.shape[0]):
                # calculate delta rho based on previous rho values
                temp_delta_rho[p] = np.dot(H[p,k], rho_i[p])
                temp_rho_i[p] = rho_i[p] + temp_delta_rho[p]*dt
                # test if index is out of range
                try:
                    delta_rho[p] = np.dot(H[p,k+1], temp_rho_i[p])
                # iterative method fails on the last timestep, so handle it
                except IndexError:
                    rho_i[p] = temp_rho_i[p]
                else:
                    rho_i[p] += dt/2 * (temp_delta_rho[p] + delta_rho[p])
                """
                if k % 30 == 0.0 and p==0.:
                    print k, tk
                    print rho_i[0]
                """
            # if we are close enough to final coherence emission, start 
            # storing these values
            if len(t) - iprime <= k:
                # each state vector has two elements that will emit (2aa and ag2)
                # transform all instances to the same rotating frame--namely, 0
                # conversions to a relevant frame are done outside of this program
                for oi in range(len(out_group)):
                    freqs = w_l[np.ix_(out_group[oi])]
                    perm_list = np.arange(H.shape[0])
                    # sum over permutations and put back in lab-time frame
                    new_emitted = rho_i[np.ix_(perm_list,out_group[oi])].sum(axis=0) * rotor(freqs*tk)
                    # f(t+dt) = f(x) + df/dt * dt
                    # add all elements within the out_group element oi together
                    rho_emitted[oi,emitted_index] = rho_emitted[oi,emitted_index] + new_emitted.sum()
                emitted_index += 1
            else:
                pass
        # if summing over ensemble, we can use the trapezoid formula
        # only use if inhom is not using quadrature rules
        if inhom_object.inhom_sampling not in ['gh']:
            if len(zeta) == 1:
                out = rho_emitted
            else:
                if l == 0:
                    #first iteration does not work--store it
                    rho_emitted_prev = rho_emitted
                else:
                    out += dzeta/2 * (rho_emitted * zweight[l] + rho_emitted_prev 
                            * zweight[l-1])
                    rho_emitted_prev = rho_emitted
                """
                # if spacing is linear, we can make this easier on ourselves
                if inhom_object.inhom_sampling in ['linear'] and len(zeta) %2==0:
                    #print l, zweight[l], w_l[1] / wn_to_omega
                    # can use simpsons rule if even; currently only works for 
                    # linearly spaced points
                    # first make sure we have the two distributions needed to 
                    # use the rule
                    if l == 0 or l == len(zeta)-1:
                        # first iteration gets a weight of 1
                        out += dzeta / 3. * rho_emitted * zweight[l]
                    elif l % 2 == 1:
                        out += 4 * dzeta/3. * rho_emitted * zweight[l]
                    else: 
                        #assert l%2 == 0
                        out += 2 * dzeta/3. * rho_emitted * zweight[l]
                else:  # use trapezoid rule
                    if l == 0:
                        #first iteration does not work--store it
                        rho_emitted_prev = rho_emitted
                    else:
                        out += dzeta/2 * (rho_emitted * zweight[l] + rho_emitted_prev 
                                * zweight[l-1])
                        rho_emitted_prev = rho_emitted
                """
        else:
            # example of gh quadrature; just take the points as they are
            if len(zeta) >  1:
                out += dzeta * rho_emitted * zweight[l]
            else: 
                out = rho_emitted
    # rho_emitted[s,t], s is out_group index, t is time index
    return out
