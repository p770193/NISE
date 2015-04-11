# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 12:17:46 2014

@author: Dan
"""

from .misc import *
from .pulse import pulse as pulse

class Axis:
    """
    a variable can have a specified scan points (self.points), 
    and associated pulses to follow with it (lock nums)
    """
    points = np.array([])
    default_units = None
    def __init__(self, pulse_ind, pulse_var, also=[], pulse_class_name='Gauss_rwa', 
                 name=None, units=None):
        self.pulse_ind = pulse_ind
        self.pulse_class_name = pulse_class_name
        self.pulse_var = pulse_var
        self.also = also
        if name:
            self.name = name
        else: self.name = None
        if units:
            self.default_units = units
        self._set_coords()
        # find the pulse class name and import the cols
        # no need to explicitly store the class--we can point to it whenever
        try:
            pulse_class = pulse.__dict__[pulse_class_name]
        except KeyError:
            print 'We could not find the specified pulse class {0}'.format(self.pulse_class_name)
        self.cols = pulse_class.cols
        
    def _set_coords(self):
        # report the coords of efield params that will be influenced by this this scan_var
        ind0 = self.also
        ind0.append(self.pulse_ind)
        print ind0
        ind1 = self.cols[self.pulse_var]
        self.coords = np.zeros((len(ind0),2), dtype=np.int)
        for ind in range(len(ind0)):
            # specify pulse and then efield param of pulse
            self.coords[ind] = [ind0[ind], ind1]
    """
    def export(self):
        self.pulse_ind
        self.pulse_var
        self.also
        self.points
        self.name
        self.pulse_class_name

    @classmethod
    def import_(cls):
        pass
    """
