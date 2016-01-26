# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 17:54:35 2013

@author: Dan
"""
#import csv
import os, sys
import csv
import cPickle as pickle
from time import clock, strftime
#from scipy.interpolate import griddata, interp2d

from .misc import *
from . import evolve as evolve
from . import pulse as pulse
#import NISE.hamiltonians.params.inhom as inhom
#import NISE.hamiltonians.H0 as H0

#reload(evolve)
#reload(pulse)
#reload(inhom)
#reload(H0)

# a default directory for dumping saves is the data folder
default_path = os.path.realpath(__file__).split('\\')[:-2]
default_path.append('data')
default_path = '\\'.join(default_path)

class Axis:
    """
    a variable can have a specified scan points (self.points), 
    and associated pulses to follow with it (lock nums)
    """
    points = np.array([])
    default_units = None
    def __init__(self, pulse_ind, pulse_var, also=None, 
                 pulse_class_name='Gauss_rwa', 
                 name=None, units=None):
        self.pulse_ind = pulse_ind
        self.pulse_class_name = pulse_class_name
        self.pulse_var = pulse_var
        if also is None:
            self.also = []
        else:
            self.also = also
        if name:
            self.name = name
        else: self.name = None
        if units:
            self.default_units = units
        # find the pulse class name and import the cols
        # no need to explicitly store the class--we can point to it whenever
        try:
            pulse_class = pulse.__dict__[pulse_class_name]
        except KeyError:
            print 'We could not find the specified pulse class {0}'.format(self.pulse_class_name)
        self.cols = pulse_class.cols
        # this initial call is just to initialize the coordinates used
        self._set_coords()
        
    def _set_coords(self):
        # report the coords of efield params that will be influenced by this 
        # scan_var
        ind0 = self.also
        ind0.append(self.pulse_ind)
        ind1 = self.cols[self.pulse_var]
        self.coords = np.zeros((len(ind0),2), dtype=np.int)
        for ind in range(len(ind0)):
            # specify pulse and then efield param of pulse
            self.coords[ind] = [ind0[ind], ind1]

# experiments are the actual running of the data
class Experiment:
    """
    Things that could potentially be set as default (class attributes):
        timestep, buffers
    Things to be specified by the experiment (instance attributes):
        the pulses module used
        number of pulses
        conjugate pulses
        pulse "positions"
    """
    timestep = timestep #4.0
    # buffers for integration
    early_buffer = early_buffer #100.0
    late_buffer  = late_buffer #400.0
    def __init__(self, pulse_class_name='Gauss_rwa', pm=[1,-1,1]):
        self.pm = pm
        self.npulses = len(pm)
        self.pulse_class_name = pulse_class_name
        # assign timestep values to the package
        # careful:  changing Experiment attributes will change the pulse defaults, too
        # if we want to change these bounds, what is the best way to do it?
        pulse_class = pulse.__dict__[pulse_class_name]
        # write time properties to the pulse class
        # do this now, or at the time of writing?
        pulse_class.timestep = self.__class__.timestep
        pulse_class.early_buffer = self.__class__.early_buffer
        pulse_class.late_buffer = self.__class__.late_buffer
        # extract properties from the pulse class--default positions are stored 
        # in pulse class
        defaults = pulse_class.defaults
        # initiate the pulse coordinates for all beams
        self.positions = np.tile(defaults, (self.npulses,1))
        # extract the lookup table
        self.cols = pulse_class.cols
        
    def set_coord(self, axis_object, pos):
        # set a default position for a laser beam
        for coord in axis_object.coords:
            self.positions[coord[0], coord[1]] = pos
            # figure out what pulse attribute this is
            # second coord index is cols lookup
            for key, value in self.cols.items():
                if value == coord[1]:
                    pulse_property = key
            pulse_number = coord[0]
            print '{0}{1} moved to {2}'.format(pulse_property, pulse_number, pos)

    def get_coords(self):
        # return the coordinates of all pulses
        for key, value in self.cols.items():
            print '{0}:  {1}'.format(key, self.positions[:,self.cols[key]])
    
    def acquire(self, H=None, inhom=None):
        # for a specified hamiltonian and using current efield positions, 
        # run the hammiltonian and time it for an estimate of runtime for scans
        pass
    
    def scan(self, *axis_objs, **sample_args):
        # "inherit" the properties needed
        Scan.timestep = self.timestep
        Scan.early_buffer = self.early_buffer
        Scan.late_buffer = self.late_buffer
        try:
            H = sample_args['H']
        except KeyError:
            print 'using default hamiltonian'
            import NISE.hamiltonians.H0 as H0
            H = H0.Omega()
        try:
            inhom_object = sample_args['inhom_object']
        except KeyError:
            # resort to the 'null' inhomogeneity
            print 'using no inhomogeneity profile as default'
            import NISE.hamiltonians.params.inhom as inhom
            inhom_object = inhom.Inhom()
        out = Scan(self.positions, self.pulse_class_name, self.pm, 
                   *axis_objs, H=H, inhom_object=inhom_object)
        # return the scan object, which can be run
        out.timestep = self.timestep
        out.early_buffer = self.early_buffer
        out.late_buffer = self.late_buffer
        return out

class Scan:
    """
    Things to be specified with each scan:
        -hamiltonian
        -inhomogeneity/correlations
        -scan axes
    things to inherit from the experiment
        -timestep, buffers
        -the pulses module used
        -number of pulses
        -conjugate pulses/phase matching
        -pulse "positions"
    """
    # a boolean that indicates whether the scan is run or not
    is_run = False
    # all the coords that are used for scanning
    npy_name = r'signal.npy'

    def __init__(self, positions, pulse_class_name, pm, *axis_objs, **kwargs):
        # create an array with the dimensions of the scan space
        # array will be as many dimensions as there are scan_objects
        # first determine creation method:  from file or from parent experiment
        if 'filepath' in kwargs.keys():
            # import the file and reinitialize all values
            return Scan._import(kwargs['filepath'])
        self.axis_objs = axis_objs
        # make positions independent of the experiment movements after set
        self.positions = positions.copy()
        self.pulse_class_name = pulse_class_name
        pulse_class = pulse.__dict__[pulse_class_name]
        self.cols = pulse_class.cols
        self.inv_cols = {v:k for k, v in pulse_class.cols.items()}
        self.npulses = self.positions.shape[0]
        self.pm = np.array(pm)
        array_shape = np.zeros((len(axis_objs)))
        # make sure we have some axes to scan
        # if not, make the array the default points
        self.coords_set = []
        if len(axis_objs) == 0:
            self.array = np.zeros((1))
        else:
            for i, obj in enumerate(self.axis_objs):
                # make sure axes refer to the same pulse class
                if obj.pulse_class_name != self.pulse_class_name:
                    raise TypeError(
                        'Scan pulse_class {0} is not the same as axis {1} pulse_class {2}'.format(
                        self.pulse_class_name, i, obj.pulse_class_name))
                # make sure all axes have points
                axis_size = obj.points.size
                if axis_size == 0:
                    # populate with the set positions
                    array_shape[i] = 1
                    point = self.positions[tuple(obj.coords[0])]
                    obj.points = np.array([point])
                else:
                    array_shape[i] = axis_size
                for coordi in obj.coords:
                    self.coords_set.append(list(coordi))
            self.array = np.zeros(array_shape)
        # Warning:  if no H/inhom objects specified, a new instance will be 
        # created from the parents
        # the above array can now be used to iterate through all objects
        # using np.ndindex
        self.H = kwargs['H']
        self.inhom_object = kwargs['inhom_object']
        self.iprime = np.arange(-self.early_buffer, 
                                self.late_buffer, 
                                self.timestep).size
        # if the hamiltonian is non-perturbative, then include phase cycling
        if self.H.pc:
            self.pc = True
            # phase cycling ranges are determined by the pm argument
            # reduced dimensionality by 1 by converting to relative phase
            pc_shape = np.zeros((self.npulses))
            # of which beam do we ignore the phase?  
            # we want to use the most perturbative pulse, because that is one 
            # that we have to account for linear terms in, so find the (first) 
            # minimum pm magnitude and make it the reference phase
            self.reference = np.where(np.abs(self.pm) == min(np.abs(self.pm)))[0][0]
            print 'pulse index {0} will be the referenced phase'.format(self.reference)
            # the problem:  we will be computing the linear polarization of 
            # this beam, so don't include it in the phase space search
            for i in range(self.npulses):
                # when phase cycling, pm should be used to indicate the highest
                # order signals to watch out for with each beam
                if i != self.reference:
                    # 2n + 1 should be safe to resolve all components, 
                    # provided no aliasing
                    # frequency components of +/- n, n-1, ..., 0
                    npts = np.abs(self.pm[i]) * 2 + 1
                    pc_shape[i] = npts
                else:
                    pc_shape[i] = 1.
            # now create the array shape for phase matching loops
            self.pc_shape = pc_shape
            self.pc_array = np.zeros((pc_shape))
            print 'encorporating phase cycling into dimensionality of scan'
            print 'array shape of {0} will be inserted into scan dimensions'.format(self.pc_shape)
        else:
            self.pc = False

    def axes(self):
        # report the axes names and their index positions
        for i, obj in enumerate(self.axis_objs):
            print i, obj.name

    def get_efield_params(self, indices=None):
        # another one-run function, like run
        # generate the array that has all paramters of efields explicitly 
        # listed
        # if this is already run, we don't have to do it again
        # optional to window certain indices
        try: 
            # prevent from re-running, as this is subject to errors if run 
            # again.  For some reason running this twice alters the params
            self.efp
            print 'stopped a recall of get_efield_params'
            return self.efp.copy()
        except AttributeError:
            if indices is None:
                indices = np.ndindex(self.array.shape)
            shape = list(self.array.shape)
            # add to each dims element a nxm array to store m params for all n pulses
            shape.append(self.npulses)
            num_args = len(self.cols)
            shape.append(num_args)
            efp = np.zeros(shape)
            for index in indices:
                this_point = np.zeros(shape[-2:])
                # loop through each axis to set points
                for i in range(len(index)):
                    # for each axis, loop through all vars changed by it
                    for j in range(len(self.axis_objs[i].coords)):
                        coords = self.axis_objs[i].coords[j]
                        this_point[coords[0], coords[1]] = self.axis_objs[i].points[index[i]]
                efp[index] = this_point
            # now, if we didn't fill in a value by going through the scan ranges, 
            # we now have to fill in with specified constants
            #filled_vals = np.zeros((self.npulses, num_args))
            for pi in range(self.npulses):
                for arg in range(num_args):
                    indices = [pi, arg]
                    if indices in self.coords_set:
                        pass
                    else:
                        # fill with default values, but where do i get the default 
                        # values from?
                        default_val = self.positions[pi, arg]
                        efp[...,pi,arg] = default_val
            #if efp.shape[:-2] == self.array.shape:
            self.efp = efp
            return efp

    def run(self, autosave=True, mp=True, chunk=False):
        if self.is_run:
            print 'scan has already run.  to execute again, make another scan object'
            return
        # initialize e-fields, using these parameters
        try:
            efp = self.efp
        except AttributeError:
            self.get_efield_params()
            efp = self.efp
        shape = list(self.array.shape)
        shape.append(len(self.H.out_group))
        shape.append(self.iprime)

        pulse_class = pulse.__dict__[self.pulse_class_name]
        pulse_class.timestep = self.timestep
        pulse_class.early_buffer = self.early_buffer
        pulse_class.late_buffer = self.late_buffer
        pulse_class.pm = self.pm
        # import the propagating function based on what the scan object tells us
        evolve_func = evolve.__dict__[self.H.propagator]
        # introduce a counter for reporting progress
        i=0
        nexti = 1
        tot = self.array.size
        if self.pc: # phase cycling
            outshape = list(shape)
            for i in xrange(len(self.pc_shape)):
                outshape.insert(-2,self.pc_shape[i])
            sig1 = np.zeros(outshape, dtype=np.complex64)
            # only loop through the scan axis coords still; output size has 
            # scaled up, however
            tot *= self.pc_array.size
            update_progress(0)
            with Timer():
                for indices in np.ndindex(self.array.shape):
                    t, efields = pulse_class.pulse(efp[indices], pm=self.pm)
                    # last two indices are outgroups and coherence time profiles
                    for jndex in np.ndindex(self.pc_array.shape):
                        phases = np.array(jndex, dtype=np.float)
                        phases*= 2*np.pi / np.array(self.pc_shape)
                        # now multiply efields by the desired phase factor
                        # don't rely on pulses to have the phase property
                        jfields = efields * rotor(phases)[:,None]
                        # signal indices are:
                        # [axes..., phases..., outgroups..., t]
                        sig1[indices][jndex] = evolve_func(t, jfields, 
                                                           self.iprime, 
                                                           self.inhom_object, 
                                                           self.H)
                        progress = i / float(tot) * 100 
                        if progress >= nexti:
                            update_progress(progress)
                            nexti = min(round(progress,0) + 1,100)
                        i += 1
                update_progress(100)
            # normalize to the number of samples taken
            sig1 = sig1 / self.pc_array.size
        else: # not phase cycling
            sig1 = np.zeros(shape, dtype=np.complex64)
            if mp:
                """
                def make_worker():
                    iprime = self.iprime
                    #eb, lb = self.early_buffer, self.late_buffer
                    pc = pulse_class
                    pm = self.pm
                    #timestep=self.timestep
                    def worker(index, efpi, 
                               H=self.H, 
                               inhom_object=self.inhom_object, 
                               pulse_class=pulse_class):
                        #pc.early_buffer = eb
                        #pc.late_buffer = lb
                        #pc.timestep = timestep
                        t, efields = pulse_class.pulse(efpi, pm=pm)
                        out = evolve_func(t, efields, iprime, inhom_object, H)
                        if indices[-1] == 0:
                            print str(indices) + str(pulse_class.timestep) + '\r',
                        return indices, out
                    return worker
                do_work = make_worker()
                arglist = [[ind, efp[ind]] for ind in np.ndindex(self.array.shape)]
                #"""    
                from multiprocessing import Pool, cpu_count
                #print self.H.__dict__
                arglist = [[ind, self.iprime, self.inhom_object, self.H, 
                            pulse_class, efp[ind], self.pm, self.timestep, 
                            self.early_buffer, self.late_buffer, evolve_func] 
                            for ind in np.ndindex(self.array.shape)]
                pool = Pool(processes=cpu_count())
                if chunk:
                    chunksize = int(self.array.size / cpu_count())
                else:
                    chunksize = 1
                print 'chunksize:', chunksize
                with Timer():
                    results = pool.map(do_work, arglist, chunksize=chunksize)
                pool.close()
                pool.join()
                # now write to the np array
                for i in range(len(results)):
                    #print results[i][0], np.abs(results[i][1]).sum()
                    sig1[results[i][0]] = results[i][1]
                del results
            else:
                update_progress(0)
                with Timer():
                    for indices in np.ndindex(self.array.shape):
                        t, efields = pulse_class.pulse(efp[indices], pm=self.pm)
                        # signal indices are:
                        # [axes..., outgroups..., t]
                        sig1[indices] = evolve_func(t, efields, self.iprime, 
                                                    self.inhom_object, self.H)
                        progress = i / float(tot) * 100 
                        if progress >= nexti:
                            update_progress(progress)
                            nexti = min(round(progress,0) + 1,100)
                        i += 1
                    update_progress(100)
        self.sig = sig1
        # reguster the scan as run
        self.is_run = True
        if autosave:
            self.save()

    def linear_signals(self):
        """
        for generating the linear signals from each beam to subtract off of 
        already acquired signals
        utilize the fact that linear signal should only have frequency of 1 =>
        only sample the beam 2 times to extract amplitude and phase
        """
        return

    def get_color(self):
        # returns array corresponding to driven signal frequency for each array 
        # point
        # supposed to be in units of *wavenumbers*
        w_axis = self.cols['w']
        wtemp = self.efp[...,w_axis].copy()
        wtemp *= self.pm
        wm = wtemp.sum(axis=-1) 
        return wm

    def kernel(self, ax, ay, fwhm_maj, fwhm_min, theta=np.pi/4., 
              center=None):
        """
        make a kernel (ellipsoidal 2D gaussian) to be used for smearing
        """
        x = self.axis_objs[ax].points.copy()
        y = self.axis_objs[ay].points.copy()
        if center is None:
            center = [x.sum() / x.size, y.sum()/y.size]
        x -= center[0]
        y -= center[1]
        s_maj = fwhm_maj / 2*np.sqrt(2*np.log(2))
        s_min = fwhm_min / 2*np.sqrt(2*np.log(2))
        # u is the major axis, y is the minor axis
        uu = np.cos(theta)*x[None,:] + np.sin(theta) * y[:,None]
        vv = np.cos(theta)*y[:,None] - np.sin(theta) * x[None,:]
        uu /= np.sqrt(2) * s_maj
        vv /= np.sqrt(2) * s_min
        kk =  np.exp(-uu**2 -vv**2)
        # normalize the distribution
        kk /= kk.sum()
        # analytical would be 2 * np.pi * s_maj * s_min, but factor 
        # in grid spacings as well (du, dv), but we'll forego this
        return kk

    def smear(self, ax, ay, fwhm_maj, fwhm_min, theta=np.pi/4., 
              center=None, save=True):
        """
        Perform a convolution of the data along scan axes using a 2-dimensional 
        elliptical gaussian.  Return a new scan object that is a copy of the 
        original but with the smeared sig array.

        ax,ay:  
            indices that specifies the two axes involved in the 
            convolution
        major:  
            number; defines the width of the gaussian along its major axis 
            (v argument)
        minor:  
            number; the width of the gaussian along its minor axis (the 
            axis perpendicular to major), in wavenumbers
        theta:  
            defines the angle of the major axis with respect to the first 
            axis.  This will be the angle required for transforming ax (ay) 
            into the major (minor) axis.
        center:  float
            specify the center of this function--not especially useful; 
            will default to the average of the axis values
        save_obj: boolean
            if True, will save the new object and record the smearing kernel
            applied
        phase: boolean
            if true, locks the phase for smearing.  Otherwise smears data 
            without any phase changes
        """
        if not self.is_run:
            print 'no data to perform convolution on!'
            return
        kk = self.kernel(ax, ay, fwhm_maj, fwhm_min, theta=theta,
                         center=center)
        # we've created the kernel; let's convolve now
        out = self.sig.copy()
        wm = self.get_color() * wn_to_omega
        tprime = np.arange(self.sig.shape[-1]) * self.timestep -self.early_buffer
        pulse_class= pulse.__dict__[self.pulse_class_name]
        d_ind = pulse_class.cols['d']
        w_ind = pulse_class.cols['w']
        for ind in np.ndindex(self.sig.shape[:-2]):
            # bring frames to their own local frequencies (all at wm)
            last_pulse = self.efp[ind][:,d_ind].max()
            ws = self.efp[ind][:,w_ind] * wn_to_omega
            p0 = np.dot(self.pm, ws * self.efp[ind][:, d_ind])
            out[ind] *= np.exp(1j * (wm[ind] * (tprime + last_pulse) - p0))
        transpose_order = list(range(len(self.sig.shape)))
        transpose_order.pop(max(ax,ay))
        transpose_order.pop(min(ax,ay))
        transpose_order.append(ay)        
        transpose_order.append(ax)
        #print transpose_order
        out = out.transpose(*transpose_order)
        
        # loop through and convolve
        #####
        from multiprocessing import Pool, cpu_count
        arglist = [[ind, out[ind], kk] for ind in np.ndindex(out.shape[:-2])]
        pool = Pool(processes=cpu_count())
        ax_size = np.product([self.axis_objs[i].points.size for i in [ax, ay]])       
        chunksize = int(self.sig.size / (cpu_count()*ax_size))
        print 'chunksize:', chunksize
        with Timer():
            results = pool.map(do_convolution, arglist, chunksize=chunksize)
        pool.close()
        pool.join()
        for i in range(len(results)):
            out[results[i][0]] = results[i][1]
        del results
        # return array to it's original order
        invert_transpose = [j for i in range(len(self.sig.shape)) 
                              for j in range(len(self.sig.shape)) 
                                  if transpose_order[j] == i]
        #print invert_transpose
        out = out.transpose(*invert_transpose)
        # undo phase adjustment
        # to minimize memory, loop through this explicitly
        for ind in np.ndindex(self.sig.shape[:-2]):
            last_pulse = self.efp[ind][:,d_ind].max()
            ws = self.efp[ind][:,w_ind] * wn_to_omega
            p0 = np.dot(self.pm, ws * self.efp[ind][:, d_ind])
            """
            # bring frames to their own local frequencies (all at wm)
            if phase_index is None:
                d0 = 0.
            else:
                last_pulse = self.efp[ind][:,d_ind].max()                    
                d0 = self.efp[ind][phase_index,d_ind] - last_pulse
            """
            #out[ind] *= np.exp(-1j*wm[ind] * (tprime-d0))
            out[ind] *= np.exp(-1j * (wm[ind] * (tprime + last_pulse) - p0))
        self.sig = out
        if save==True:
            self.save(name='smeared')
            # also record the kernel
            kernel = {
                'angle (deg)' : theta *180. / np.pi,
                'fwhm_maj'    : fwhm_maj,
                'fwhm_min'    : fwhm_min,
                'center'      : center,
                #'phase_index' : phase_index
            }
            p_name = r'kernel.csv'
            p_full_name = '\\'.join([self.output_folder, p_name])
            with open(p_full_name,'wb') as params:
                writer = csv.writer(params)
                writer.writerow(['----- kernel params -----'])
                for key, value in kernel.iteritems():
                    writer.writerow([key, value])
            from . import fscolors_beta as f
            import matplotlib.pyplot as plt
            plt.figure()
            plt.contourf(x, y, kk, 200, cmap=f.plot_artist.mycm)
            plt.grid()
            plt.colorbar()
            plt.savefig(self.output_folder + r'\smear kernel.png')
            plt.close()
        print 'smearing complete!'
        return 

    def efields(self, windowed=True):
        """
        return the e-fields used in the simulation
        if windowed == True, only returns values that are within the early 
        and late buffer
        """
        # [axes..., numpulses, nparams]
        efp = self.get_efield_params()
        # [axes..., numpulses, pulse field values]
        efields_shape = list(efp.shape)
        pulse_class= pulse.__dict__[self.pulse_class_name]
        if windowed:
            efields_shape[-1] = self.iprime
            efields = np.zeros((efields_shape), dtype=np.complex)
            with Timer():
                for ind in np.ndindex(tuple(efields_shape[:-2])):
                    ti, efi = pulse_class.pulse(efp[ind], pm=self.pm)
                    efields[ind] = efi[:,-self.iprime:]
        else: 
            # figure out the biggest array size we will get
            d_ind = pulse_class.cols['d']
            t = pulse_class.get_t(efp[...,d_ind])
            # now that we know t vals, we can set fixed bounds
            pulse_class.fixed_bounds_min =  t.min()
            pulse_class.fixed_bounds_max =  t.max()
            pulse_class.fixed_bounds = True
            efields_shape[-1] = t.size
            efields = np.zeros((efields_shape), dtype=np.complex)
            try:
                with Timer():
                    for ind in np.ndindex(tuple(efields_shape[:-2])):
                        ti, efi = pulse_class.pulse(efp[ind], pm=self.pm)
                        efields[ind] = efi
            finally:
                # set the class back to what it was before exiting
                pulse_class.fixed_bounds = False
        return efields
                    
    def nrb(self):
        # compute non-res 'forced' signal component 
        if self.is_run:
            efields = self.efields()
            efields = efields.prod(axis=-2)
            ratio = np.abs(self.sig.sum(axis=-2)).sum(axis=-1).max() / np.abs(efields).sum(axis=-1).max()
            print 'sig:nrb ratio of {0}'.format(ratio)
            return efields
        else:
            print 'cannot plot data, scan object has not been run yet'

    def plot(self):
        if self.is_run:
            # plot the obtained data
            print 'functionality has not been implemented yet'
        else:
            print 'cannot plot data, scan obejct {0} has not been run yet'

    def save(self, name=None, full_name=None, method='pickle'):
        """
        saves the scan instance
        name:
            name will appended to the timestamp
        full_name:
            saved folder will be named full_name with no timestamp
        """
        #-----------step 1:  create output_folder-----------
        # create a directory for all the data files
        foldertime = strftime("%Y.%m.%d %H-%M-%S")
        if full_name is None:
            # even if already saved, rewrite foldername property to the most recent
            # save
            if name is None:
                self.foldername = foldertime
            else:
                self.foldername = ' - '.join([foldertime, name])
            # create folder for output files 
            output_folder = r'\\'.join([default_path, self.foldername])
        else:
            output_folder = full_name
        # save the output folder for pointing in other exports
        self.output_folder = output_folder
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        #-----------step 2:  write signal to file (if exists)-----------
        if self.is_run:
            while True:
                npy_full_name = r'\\'.join([output_folder, self.__class__.npy_name])
                np.save(npy_full_name, self.sig)
                break
            # if we make it out of the loop, we can safely delete the data
            # 'hide' the sig data from pickle, but store a local copy for 
            # restoring after the save
            sig = self.sig
            del self.sig
        #-----------step 3:  copy subclass definition files-----------
        src_folder = output_folder + r'\src'
        os.makedirs(src_folder)
        # copy the files and save the module and name info
        # gather information from all relevant objects and write them to files
        H_file = _obj_classfile_copy(self.H, src_folder)
        inhom_file = _obj_classfile_copy(self.inhom_object, src_folder)
        old_dict = {
            'H':       [self.H.__class__.__name__, 
                        self.H.__class__.__module__,
                        H_file.split('\\')[-1]],
            'Inhom':   [self.inhom_object.__class__.__name__, 
                        self.inhom_object.__class__.__name__,
                        inhom_file.split('\\')[-1]]
        }
        # write dictionary to file for lookup when importing
        pickle_name = r'class_maps.p'
        pickle_full_name = '\\'.join([src_folder, pickle_name])
        try:
            with open(pickle_full_name, 'wb') as out_s:
                pickle.dump(old_dict, out_s, protocol=0)
        except IOError:
            print 'could not write class_maps properly'
        #-----------step 4:  serialization of scan object-----------
        # write the explicit objects to file
        pickle_name = r'scan_object.p'
        pickle_full_name = '\\'.join([output_folder, pickle_name])
        out_s = open(pickle_full_name, 'wb')
        try:
            pickle.dump(self, out_s)
        finally:
            out_s.close()
        if self.is_run:
            self.sig = sig
        #-----------step 5:  write useful properties to readable file-----------
        p_name = r'report.csv'
        p_full_name = '\\'.join([output_folder, p_name])
        with open(p_full_name,'wb') as params:
            writer = csv.writer(params)
            writer.writerow(['----- {0} -----'.format([self.__class__])])
            for index, axis in enumerate(self.axis_objs):
                writer.writerow(['----- Axis {0} -----'.format([index])])
                lis1 = []
                for i in range(len(axis.coords)):
                    stri = '{0}{1}'.format(self.inv_cols[axis.coords[i][1]], axis.coords[i][0]+1)
                    lis1.append(stri)
                writer.writerow(['vars', lis1])
                writer.writerow(['npts', axis.points.size])
                writer.writerow(['range', axis.points.min(), axis.points.max()])
            writer.writerow(['----- {0} -----'.format([self.__class__])])
            writer.writerow(['timestep', self.timestep])
            writer.writerow(['npulses', self.npulses])
            writer.writerow(['pulse_class', self.pulse_class_name])
            for i in range(self.npulses):
                positions_formatted = ['{:+.2f}'.format(member) for member in self.positions[i]]
                positions_formatted.insert(0, 'default positions')
                writer.writerow(positions_formatted)
            writer.writerow(['pm', self.pm])
            writer.writerow(['array shape', self.array.shape])
            writer.writerow(['pc', self.pc])
            if self.pc:
                writer.writerow(['pc shape', self.pc_shape])
            writer.writerow(['is_run', self.is_run])
            
            for obji in [self.H, self.inhom_object]:
                writer.writerow(['----- {0} -----'.format([obji.__class__])])
                for key in obji.__class__.out_vars:
                    # currently having issues with separating class defaults 
                    # and searching the instance dictionary.  This is the best
                    # I can do right now...
                    try:
                        value = getattr(obji, key)
                    except AttributeError:
                        value = getattr(obji.__class__, key)
                    writer.writerow([key, value])
        print 'save complete:  output directory is {0}'.format(output_folder)
        return output_folder
    
    @classmethod
    def _import(cls, foldername):
        """
        imports using pickle--uses NISE library if available, otherwise 
        attempts to import using saved modules
        """
        #-----------step 1:  import saved class definition files-----------
        #-----------and define class method/name reassignments-----------
        # find the tables (if they are there)
        try:
            scan_obj = pickle.load(open(foldername + r'\scan_object.p','rb'))
        except pickle.PicklingError:
            print 'using copied py files to reconstruct subobjects of scan instance'
            src_folder = foldername + r'\src'
            class_maps_file = src_folder + r'\class_maps.p' 
            with open(class_maps_file) as f:
                old_map = pickle.load(f)
            module_map = dict()
            cls_names = []
            for key, value in old_map.items():
                # names will remain unchanged, but modules will change
                key_name, key_module, key_file = value
                cls_names.append(key_name)
                key_filepath = '\\'.join([src_folder, key_file])
                if os.path.isfile(key_filepath):
                    # name the module uniquely so that overwrites are unlikely
                    new_module = '.'.join([foldername.split('\\')[-1],
                                           key_module.split('.')[-1]])
                    while new_module in sys.modules:
                        # just lengthen the string
                        new_module += 'f'
                    # import the file and aquire the module name
                    from imp import load_source
                    load_source(new_module, key_filepath)
                    # now modules can be mapped
                    module_map[key_module] = new_module
                    #sys.modules[new_module] accesses the module
            #-----------step 2:  import pickled scan-----------
            # adapted from 
            # http://stackoverflow.com/questions/13398462/
            # unpickling-python-objects-with-a-changed-module-path
            from importlib import import_module
            def map_cls(mod_name, kls_name):
                # determine whether or not name is one of our mapped names
                if kls_name in cls_names: 
                    print mod_name, kls_name
                    for key, value in old_map.items():
                        key_name, key_module, key_file = value
                        if kls_name == key_name:
                            new_module = module_map[key_module]
                            mod = sys.modules[new_module]
                            return getattr(mod, kls_name)
                else:
                    #print '{0} not mapped'.format(mod_name)
                    mod = import_module(mod_name)
                    return getattr(mod, kls_name)
    
            # need cpickle for find_global attribute
            def loads(filename):
                with open(filename,'rb') as sfile:
                    unpickler = pickle.Unpickler(sfile)
                    unpickler.find_global = map_cls
                    return unpickler.load() # object will now contain the new class path reference
    
            scan_obj = loads(foldername + r'\scan_object.p')
        #-----------step 3:  import sig data (if applicable)-----------
        # import the npy file if it was saved
        if scan_obj.is_run:
            npy_full_name = '\\'.join([foldername, cls.npy_name])
            try:
                scan_obj.sig = np.load(npy_full_name)
            except:
                print 'The scan data array could not be imported'  
                print 'You will have to rerun the simulation'
                scan_obj.is_run = False
        #----------step 4a:  classes that are not imported----------#
        pulse_class = pulse.__dict__[scan_obj.pulse_class_name]
        # write time properties to the pulse class
        pulse_class.timestep = scan_obj.timestep
        pulse_class.early_buffer = scan_obj.early_buffer
        pulse_class.late_buffer = scan_obj.late_buffer
        #----------step 4b:  update output_folder-------------#
        new_folder = os.path.normcase(
                        os.path.normpath(
                            os.path.abspath(foldername)))
        try: 
            old_folder = os.path.normcase(
                            os.path.normpath(
                                os.path.abspath(scan_obj.output_folder)))
            # check if output folder matches current directory
            if new_folder != old_folder:
                # rewrite the place where the new folder is
                scan_obj.output_folder = new_folder
        except AttributeError:
            # output_folder was never written, assign it the current folder
            scan_obj.output_folder = new_folder
        return scan_obj

    def k_filter(self, vector):
        # select one phased component of the polarization 
        # only works on phase-cycled scans
        pass
    def compare(self):
        # compare results to some other gridded data
        pass

    @classmethod
    def old_import_(cls, filepath):
        # method imports a pickle file and instantiates the object
        in_pickle = open(filepath, 'rb')
        try:
            # Write to the stream
            scan_obj = pickle.load(in_pickle)
        finally:
            in_pickle.close()
        # import the npy file if it was saved
        if scan_obj.is_run:
            import_folder = os.path.dirname(filepath)
            npy_full_name = r'\\'.join([import_folder, cls.npy_name])
            try:
                scan_obj.sig = np.load(npy_full_name)
            except:
                print 'could not import the scan data array.  You will have to rerun the simulation'
                scan_obj.is_run = False
        return scan_obj
        """
        # half-worked-out solution adapted from 
        # https://wiki.python.org/moin/UsingPickle/RenamingModules
        def mapname(name):
            print name
            # a simple checker to run through the dictionary
            if name in module_map:
                print name
                return module_map[name]
            return name

        def mapped_load_global(self):
            print 'in here!'
            # compiles the module and name to point to a given class
            module = mapname(self.readline()[:-1])
            name = mapname(self.readline()[:-1])
            klass = self.find_class(module, name)
            self.append(klass)

        def loads(filename):
            #sfile = StringIO(filename)
            with open(filename, 'rb') as sfile:
                unpickler = pickle.Unpickler(sfile)
                #print unpickler.load_global
                unpickler.load_global = mapped_load_global
                #unpickler.dispatch[pickle.GLOBAL] = mapped_load_global
                return unpickler.load()
        """

    def old_save(self, foldername=None, filename=None, method='pickle'):
        # gather information from all relevant objects and write them to files
        # create a directory for all the data files
        if foldername is None:
            foldername = strftime("%Y.%m.%d %H-%M-%S")
        # create folder for output files 
        output_folder = r'\\'.join([default_path, foldername])
        os.makedirs(output_folder)
        self.output_folder = output_folder
        # first, write the sig array to a npz file (if it exists)
        if self.is_run:
            while True:
                npy_full_name = r'\\'.join([output_folder, self.__class__.npy_name])
                np.save(npy_full_name, self.sig)
                break
            # if we make it out of the loop, we can safely delete the data
            # 'hide' the sig data from pickle, but store a local copy for 
            # after the saving
            # look into pytables for future storage options
            sig = self.sig
            del self.sig
        # write the explicit objects to file
        pickle_name = r'scan_object.pickle'
        pickle_full_name = r'\\'.join([output_folder, pickle_name])
        out_s = open(pickle_full_name, 'wb')
        try:
            # Write to the stream
            # coming:  don't save the instance, but save variables that can 
            # recreate the instance
            # all attributes that are not functions or classes/class instances
            # can be dumped directly (perhaps as a dictionary)
            # how to work with, say, axes, pulse, and such?
            pickle.dump(self, out_s)
        finally:
            out_s.close()
        self.sig = sig
        # create dictionaries of the information that should be well-formated
        # and readable; to be exported as crv or something similar
        p_name = r'report.csv'
        p_full_name = r'\\'.join([self.output_folder, p_name])
        params = open(p_full_name,'wb')
        writer = csv.writer(params)
        try:
            writer.writerow(['----- {0} -----'.format([self.__class__])])
            for index, axis in enumerate(self.axis_objs):
                writer.writerow(['----- Axis {0} -----'.format([index])])
                lis1 = []
                for i in range(len(axis.coords)):
                    stri = '{0}{1}'.format(self.inv_cols[axis.coords[i][1]], axis.coords[i][0]+1)
                    lis1.append(stri)
                writer.writerow(['vars', lis1])
                writer.writerow(['npts', axis.points.size])
                writer.writerow(['range', axis.points.min(), axis.points.max()])
            writer.writerow(['----- {0} -----'.format([self.__class__])])
            writer.writerow(['timestep', self.timestep])
            writer.writerow(['npulses', self.npulses])
            writer.writerow(['pulse_class', self.pulse_class_name])
            for i in range(self.npulses):
                positions_formatted = ['{:+.2f}'.format(member) for member in self.positions[i]]
                positions_formatted.insert(0, 'default positions')
                writer.writerow(positions_formatted)
            writer.writerow(['pm', self.pm])
            writer.writerow(['array shape', self.array.shape])
            writer.writerow(['pc', self.pc])
            if self.pc:
                writer.writerow(['pc shape', self.pc_shape])
            writer.writerow(['is_run', self.is_run])
            
            for obji in [self.H, self.inhom_object]:
                writer.writerow(['----- {0} -----'.format([obji.__class__])])
                for key in obji.__class__.out_vars:
                    # currently having issues with separating class defaults 
                    # and searching the instance dictionary.  This is the best
                    # I can do right now...
                    try:
                        value = getattr(obji, key)
                    except AttributeError:
                        value = getattr(obji.__class__, key)
                    writer.writerow([key, value])
        finally:
            params.close()
        print 'output sent to {0}'.format(output_folder)
    
    """
    @classmethod
    def importbeta_(cls, filepath):
        # method imports a pickle file and instantiates the object
        in_pickle = open(filepath, 'rb')
        try:
            # Write to the stream
            scan_obj = pickle.load(in_pickle)
        finally:
            in_pickle.close()
        # import the npy file if it was saved
        if scan_obj.is_run:
            import_folder = os.path.dirname(filepath)
            npy_full_name = r'\\'.join([import_folder, cls.npy_name])
            try:
                scan_obj.sig = np.load(npy_full_name)
            except:
                print 'could not import the scan data array.  You will have to rerun the simulation'
                scan_obj.is_run = False
        return scan_obj

    def important_export_vars(self):
        # reassign these class variables before input:
        self.timestep
        self.early_buffer
        self.late_buffer
        
        # initialization arguments
        # def __init__(self, positions, pulse_class_name, pm, *axis_objs, **kwargs):
        self.positions 
        self.pulse_class_name # will gather cols, inv_cols
        self.pm
        axes = []
        for axes_obj in self.axis_objs:
            axes.append(axes_obj.export())
        # I would really like to just pickle these, so that's what I'm going to do?
        self.H = kwargs['H']
        self.inhom_object = kwargs['inhom_object']
        # output folder will be determined by the inport command
        # just look for .npy files, or .pickle files--don't worry about internalizing the input
        self.output_folder, self.npy_name
        # to be gotten from the npy file
        self.is_run
        self.efp
        self.sig

        
        # these parameters will be derived and need not be imported
        self.pm
        self.npulses
        self.coords_set
        self.pc, self.pc_array, self.pc_shape        

        H_out = self.H.export()
        inhom_out = self.inhom_object.export()


        self.iprime
        self.cols
        self.inv_cols
    """
    

class Timer:
    def __enter__(self, progress=None, verbose=True):
        self.verbose = verbose        
        self.start = clock()
    def __exit__(self, type, value, traceback):
        self.end = clock()
        self.interval = self.end - self.start
        if self.verbose:
            print 'elapsed time: {0} sec'.format(self.interval)

def update_progress(progress, carriage_return = True, length = 50):
    '''
    prints a pretty progress bar to the console
    accepts 'progress' as a percentage
    carriage_return toggles overwrite behavior
    '''
    #make progress bar string
    progress_bar = ''
    num_oct = int(progress * (length/100.))
    progress_bar = progress_bar + '[{0}{1}]'.format('#'*num_oct, ' '*(length-num_oct))
    progress_bar = progress_bar + ' {}%'.format(np.round(progress, decimals = 2))
    if carriage_return:
        progress_bar = progress_bar + '\r'
        print progress_bar,
        return
    if progress == 100:
        progress_bar[-2:] = '\n'
    print progress_bar

def _obj_classfile_copy(obj, destination):
    # properly handle an object for export using pickle
    class_file = sys.modules[obj.__module__].__file__
    # copy the uncompiled version
    class_file = class_file.split('.')
    class_file[-1] = 'py'
    class_file = '.'.join(class_file)
    from shutil import copy2 as _copy
    _copy(class_file, destination)
    del _copy
    return class_file

from scipy.signal import convolve2d

def do_convolution(arglist):
    indices, signal, kk = arglist
    out = convolve2d(signal, kk, mode='same', fillvalue=0.0)
    return indices, out

def do_work(arglist):
    indices, iprime, inhom_object, H, pulse_class = arglist[:5]
    efpi, pm, timestep, eb, lb, evolve_func = arglist[5:]
    # need to declare these for each function
    pulse_class.early_buffer = eb
    pulse_class.late_buffer = lb
    pulse_class.timestep = timestep
    t, efields = pulse_class.pulse(efpi, pm=pm)
    out = evolve_func(t, efields, iprime, inhom_object, H)
    if indices[-1] == 0:
        #print str(inhom_object.zeta_bound) + str(len(inhom_object.zeta)) + '\r',    
        #print str(H.tau_ag), str(H.tau_2ag), '\r'
        print str(indices),  str(pulse_class.timestep), str(iprime) + '              \r',
    return indices, out
