# -*- coding: utf-8 -*-
"""
Created on Sat Jun 21 14:07:53 2014

@author: Dan
"""

from __future__ import absolute_import, division, print_function, unicode_literals

from NISE.lib.misc import *

class Inhom():
    # class contains the list of weights and sampling values to use
    #--------------------------Recorded attributes--------------------------
    out_vars = ['inhom_sampling', 'dist_params']
    #--------------------------Methods--------------------------
    def __init__(self, inhom_sampling=None, **dist_params):
        """
        generates the list of sampling points in the distribution and their weights
        inhom dists should be normalized (int(f, dzeta) = 1.)
        """
        # inherit all class attributes unless kwargs has them; then use those 
        # values.  if kwargs is not an Omega attribute, it gets ignored
        for key, value in dist_params.items():
            setattr(self, key, value)
        #print self.__dict__.items()
        # eliminating other quadrature methods; linear works best anyways
        if inhom_sampling == 'linear':
            # currently the only inhomogeneity parameter that can normalize well 
            # in relation to the case of no inhomogeneity
            if isinstance(dist_params.get('num'), int):
                num = dist_params.get('num')
            else:
                try:
                    num = int(num)
                except TypeError:
                    print('no distribution sampling number specified; using 10 points as default')
                    num = 10
            if 'zeta_bound' in dist_params.keys():
                zeta_bound = dist_params.get('zeta_bound')
            else:
                zeta_bound = 3
            zeta = np.linspace(-zeta_bound, zeta_bound, num=num)
            # need parameter 'sigma'
            sigma = dist_params.get('sigma')
            # scale our sampling intervals according to sigma
            zeta = zeta * sigma
            self.zweight = 1 / (np.sqrt(2*np.pi)*sigma) * np.exp(- 0.5 * ((zeta / sigma)**2))
            self.dzeta = np.abs(zeta[1] - zeta[0])
            self.zeta = zeta
        elif inhom_sampling == 'rect':
            w = dist_params.get('w')
            if isinstance(dist_params.get('num'), int):
                num = dist_params['num']
            else:
                try:
                    num = int(num)
                except TypeError:
                    print('no distribution sampling number specified; using 10 points as default')
                    num = 10
            self.zeta = np.linspace(-w,w,num=num)
            self.dzeta = np.abs(self.zeta[1] - self.zeta[0])
            self.zweight = np.ones(self.zeta.shape)
        elif inhom_sampling == 'gh':
            import NISE.hamiltonians.params.gauss_hermite as gh
            # gaussian-hermite quadrature
            # see http://en.wikipedia.org/wiki/Gauss%E2%80%93Hermite_quadrature
            # for details
            n = dist_params.get('n')
            try:
                gh.quad[n]
            except KeyError:
                print('no table for quadrature of number {0} is available'.format(n))
                print('available quadrature numbers:  {0}'.format(str(gh.quad.keys())))
            sigma = dist_params.get('sigma')
            self.zeta = np.array(gh.quad[n])[0]
            self.zweight = np.array(gh.quad[n])[1]
            self.dzeta = 1.
            # substitution to inhom variables yields the following scaling:
            self.zeta*= np.sqrt(2) * sigma
            self.zweight*= np.pi**-0.5
        else:
            self.zeta = np.array([0])
            self.zweight = [1.0]
            self.dzeta = 1.0
        self.inhom_sampling = inhom_sampling
        self.dist_params = dist_params
