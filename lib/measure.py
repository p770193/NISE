"""
    to do:  make an object that can store multiple scan objects, so all the 
    hamiltonians can collect here
"""

from .misc import *
from . import fscolors_beta as fscolors_beta
import matplotlib.pyplot as plt
from time import strftime

#wrightcm = ['#FFFFFF','#0000FF','#00FFFF','#00FF00','#FFFF00','#FF0000','#881111']
#mycm=mplcolors.LinearSegmentedColormap.from_list('wright',wrightcm)

class Measure:
    """
    a class that processes the scanned polarization in a specified way
    you can control the measurement methods, and make a string of commands to 
    issue resulting in an actual signal
    
    should phase cycling be external?
    """
    
    # placeholder:  a list of objects similar to that in Scan, that organizes 
    # all the axes encountered in the Measure().pol array
    axis_objs = None
    
    plot_lim = 100.
    # initialization should specify the measurement instructions
    def __init__(self, scan_obj, *instruction_classes, **kwargs):
        if not scan_obj.is_run:
            print 'cannot measure scan until signal is computed with \'run\' method'
            return
        # if all is well, internalize the scan objs
        # pointing will save space...
        self.scan_obj = scan_obj
        self.space = 't'
        self.instruction_classes = instruction_classes
        # the following commented out stuff is partial work for getting multiple
        # scans measured at once--not ready yet
        """
        # loop through the scan objects
        try:
            len(scan_objs)
        except:
            # only one scan object, so make it a list to accomodate loops
            scan_objs = list(scan_objs)
        shape = [None,None]
        for i, obj in scan_objs.enumerate():
            # for the moment, all scans must have same dimension in order to be combined
            if i==0: 
                shape[0] = obj.array.shape
                shape[1] = obj.sig.shape[-1]
            # though not sufficient, only check for correct shapes at this point
            if obj.array.shape != shape[0] or obj.sig.shape[-1] != shape[1]:
                raise TypeError('scans must be of the same dimension.')
            # if all checks out, add together the polarizations for processing?
            self.pol += 
        """
    def run(self, save=True):
        # store polarization of each step as a variable
        self.pol = self.scan_obj.sig.copy()
        # convert the polarization into emitted field via a phase shift
        self.pol *= 1j
        # loop through the manipulations
        for step in self.instruction_classes:
            # make sure operation is valid the space the signal is expressed in
            if step.check_space(self.space):
                self.pol = step.run(self.scan_obj, self.pol)
            # reset space to new space
            self.space = step.out_space
        # finally, save the data
        if save:
            output_folder = None
            try:
                output_folder = self.scan_obj.output_folder
            except AttributeError:
                # output folder is not defined
                output_folder = NISE_path + r'\data'
                print 'Could not find output folder.  Plots written to {0}'.format(
                    output_folder)
            pol_folder = output_folder + r'\measured'
            mkdir_p(pol_folder)
            pol_args = {obj.name:obj.points for obj in self.scan_obj.axis_objs}
            if 'pol' not in pol_args.keys():
                pol_args['pol'] = self.pol
            else:
                raise TypeError('pol name is reserved; cannot be used for axis_objs name')
            pol_name = pol_folder + '\\' + strftime("%Y.%m.%d %H-%M-%S")
            np.savez(pol_name, **pol_args)
            print 'pol saved as {0}'.format(pol_name)

    def plot(self, xaxis, yaxis=None, artist=None, zoom=None, folder_name=None):
        """
        plot the data (or a subset of it, if specified)
        x and y axis specify the axes to plot against
        """
        out = []
        # do i have an x and y axis, or just x?
        xaxis = int(xaxis)
        if yaxis is not None:
            yaxis = int(yaxis)
            plot_dim = 2
            var_list = [xaxis, yaxis]
        else:
            plot_dim = 1
            var_list = [xaxis]
        # gather the shape of pol:  neglect trivial array indices
        try:
            dim = len([i for i in range(len(self.pol.shape)) if self.pol.shape[i] > 1])
        except NameError:
            print 'running detection operations...'
            self.run()
            dim = len(self.pol.shape)
        # make sure parse_axes is legit
        if xaxis > dim or xaxis < 0:
            print 'parse_axes specifies axes out of bounds'
            raise TypeError
        # plotting process differs based on number of plots and plot dimensions
        # new strategy:  make y the second-to-last axis, x the last
        transpose_order = list(range(len(self.pol.shape)))
        if yaxis:
            transpose_order.pop(max(xaxis,yaxis))
            transpose_order.pop(min(xaxis,yaxis))
            transpose_order.append(yaxis)
        else:
            transpose_order.pop(xaxis)
        transpose_order.append(xaxis)
        temp_pol = self.pol.transpose(*transpose_order)
        # use temp_pol and order is [...,yaxis, xaxis]        
        # all that's left is figuring out how many plots there are!
        iterator_shape = temp_pol.shape[:-len(var_list)]
        #print transpose_order
        index_map=[]
        # we will make an exception if iterator_shape is empty--this can happen
        # if plot_dim is equal or greater to dim
        if dim - plot_dim >= 1:
            # if more than one plot, then we must loop through
            # present a series of 2d contours along the specified axis
            # make sure we're not going to make too many windows:
            num_plots = np.prod(iterator_shape)
            if num_plots > Measure.plot_lim: 
                print 'number of plots must be less than Measure.plot_lim = ', Measure.plot_lim
                raise TypeError
            for i in xrange(len(self.pol.shape)):
                if i not in var_list:
                    index_map.append(i)
        elif dim - plot_dim < 0:
            # unable to use:  plot dimensionality is higher than actual scan dimensions
            print 'requested dimensionality exceeds sample space'
            raise TypeError
        # declare the artist instance for plotting
        if isinstance(artist, fscolors_beta.plot_artist):
            # treatment for a predefined artist
            pass
        else:
            # make a new artist
            x_obj = self.scan_obj.axis_objs[xaxis]
            #xvar = x_obj.pulse_var
            #xlabel = x_obj.name
            if yaxis:
                y_obj = self.scan_obj.axis_objs[yaxis]
                artist = fscolors_beta.plot_artist(xaxis=x_obj, yaxis=y_obj)
            else:
                # 1D scans
                artist = fscolors_beta.plot_artist(xlabel=x_obj)

        # define a global max
        z_max, z_min = self.pol.max(), self.pol.min()
        try:
            output_folder = self.scan_obj.output_folder
        except AttributeError:
            # output folder is not defined
            output_folder = NISE_path + r'\data'
            print 'Could not find output folder.  Plots written to {0}'.format(
                output_folder)
        nstr = r''
        if dim - plot_dim >= 1:
            # many images, so generate its own output folder
            if folder_name:
                output_folder += r'\{0}'.format(folder_name)
            else:
                output_folder += r'\img ' + strftime("%Y.%m.%d %H-%M-%S")
                output_folder += r' - x={0}, y={1}'.format(xaxis, yaxis)
                mkdir_p(output_folder)
            # a loop to write plots (if necessary)
            for indices in np.ndindex(np.array(iterator_shape)):
                tstr = r''
                nstr = r'{0} '.format(indices)
                # generate labels for the plot (tstr) and for the filename (nstr):
                for i in xrange(len(indices)):
                    index = indices[i]
                    tstr += self.scan_obj.axis_objs[index_map[i]].name
                    tstr += ' = {0} '.format(self.scan_obj.axis_objs[
                                index_map[i]].points[index])
                    nstr += self.scan_obj.axis_objs[index_map[i]].pulse_var
                    nstr += ' = {0} '.format(self.scan_obj.axis_objs[
                                index_map[i]].points[index])
               # select the x,y and z data
                """
                i=0
                j=0
                # whittle array down to what we want?
                full_index = []
                while i < dim-plot_dim:
                    if j in var_list: 
                        # include all points
                        full_index.append(slice(shape[j]))
                    else:
                        full_index.append(indices[i])
                        i += 1
                    j += 1
                    if j > dim:
                        # failsafe
                        'indexing loop failed to terminate correctly'
                        raise RuntimeError
                        break
                del i,j
                z = self.pol[full_index]
                if xaxis < yaxis:
                    z = z.T
                """
                z = temp_pol[indices]
                x = self.scan_obj.axis_objs[xaxis].points.copy()
                #plt.figure()
                if plot_dim ==2:
                    y = self.scan_obj.axis_objs[yaxis].points.copy()
                    xyz = fscolors_beta.xyz(x,y,z, zmax=z_max)
                    if zoom is not None:
                        xyz.zoom(zoom)
                    artist.plot(xyz, alt_z='amp', contour=True)
                    # title has to be written before more objects are printed
                    # to keep focus on main subplot
                    plt.title(tstr)
                    artist.colorbar()
                    out.append(xyz)
                save_path = output_folder + '\\' + r''.join([nstr,'.png'])
                artist.savefig(fname = save_path)
        else:
            if folder_name:
                output_folder += r'\{0}'.format(folder_name)
            else:
                output_folder += r'\img ' + strftime("%Y.%m.%d %H-%M-%S")
                output_folder += r' - x={0}, y={1}'.format(xaxis, yaxis)
                mkdir_p(output_folder)
            indices = range(dim)
            # generate labels for the filename (nstr):
            for i in xrange(len(indices)):
                index = indices[i]
                nstr = self.scan_obj.axis_objs[index].pulse_var
            z = self.pol
            # meant to handle length-1 dimensions
            while len(z.shape) > plot_dim:
                z = z[0]
            # orient correctly
            if xaxis < yaxis:
                z = z.T
            x = self.scan_obj.axis_objs[xaxis].points
            #plt.figure()
            if plot_dim ==2:
                y = self.scan_obj.axis_objs[yaxis].points
                xyz = fscolors_beta.xyz(x,y,z, zmax=z_max)
                if zoom is not None:
                    xyz.zoom(zoom)
                artist.plot(xyz, alt_z='amp', contour=True)
                artist.colorbar()
            save_path = output_folder + '\\' + r''.join([nstr,'.png'])
            artist.savefig(fname = save_path)
            out.append(xyz)
        return out

class Mono:
    """
    fourier transform and (optional) apply a slit
    all frequencies are tabulated relative to the driven signal frame
    slitwidth in wavenumbers
    compute central frequency from phase matching
    """
    slitwidth = 120.0 # default slitwidth for mono (in wavenumbers)
    in_space = 't'
    out_space = 'f'
    # aliasing issue:  data is under-aliased, so translating frequency 
    # by a non-gridded value can introduce aliasing errors
    # problem can be minimized by adding trailing zeroes, effectively 
    # interpolating the fft data set
    alias_factor = 2

    @classmethod
    def run(cls, scan_obj, sig):
        tprime = np.arange(sig.shape[-1]) * scan_obj.timestep 
        if sig.shape[-1] == scan_obj.iprime:
            # why do I do this again?  I phase all signals to the arrival of
            # the last pulse, but why?
            #tprime -= scan_obj.early_buffer
            pass
        else: # sig has been stretched and put into absolute coordinates
            pass
        #fft_tprime = np.fft.fftfreq(sig.shape[-1], d=scan_obj.timestep)
        #tprime = np.fft.fftshift(np.fft.fftfreq(sig.shape[-1], d=np.abs(fft_tprime[1]-fft_tprime[0])))
        # apply rw to orient polarization to driven frame (where we 
        # traditionally put the mono)
        # need scan_obj.efp to determine the rw
        wm = scan_obj.get_color() * wn_to_omega # sometimes omega, sometimes wavenumbers
        out_shape = list(sig.shape)
        out_shape[-1] *= cls.alias_factor+1
        out = np.zeros(out_shape,dtype=np.complex64)
        tz_shape = list(sig.shape[-2:])
        tz_shape[-1] *= cls.alias_factor
        trailing_zeroes = np.zeros(tz_shape)
        #tprime2 = np.fft.fftfreq(out.shape[-1], 
        #                         d=np.abs(fft_tprime[1]-fft_tprime[0]))
        fft_tprime2 = np.fft.fftfreq(out.shape[-1], 
                                    d=scan_obj.timestep)
        for ind in np.ndindex(sig.shape[:-2]):
            wmi = wm[ind]
            outi = sig[ind] * np.exp(1j*wmi * tprime)
            outi = np.concatenate([outi, trailing_zeroes], axis=-1)
            outi = np.fft.fft(outi, axis=-1)
            outi = np.fft.fftshift(outi, axes=-1)
            out[ind] = outi
        if isinstance(cls.slitwidth, float):
            # convert to wavenumbers for mono functions
            var = (np.fft.fftshift(fft_tprime2) * 2 * np.pi / wn_to_omega)
            mono_func = np.zeros(var.shape)
            for i in range(len(var)):
                # new:  use triangle function (convolved slitwidths of width 
                # specified)
                if np.abs(var[i]) <= cls.slitwidth:
                    mono_func[i] = 1.0 * (cls.slitwidth - np.abs(var[i])) / cls.slitwidth
            out *= mono_func

        # redefine fourier and time axes after trailing zeroes are added
        #tprime = np.arange(out.shape[-1]+tz_shape[-1])*scan_obj.timestep
        """
        #plt.figure()
        #plt.plot(out[5,5,0])
        #plt.figure()
        #out = sig * np.exp(1j*wm[..., None, None]*wn_to_omega*tprime)
        # last axis is always the time dimension
        # to minimize aliasing, add trailing zeroes to signals
        # do this after mono because trailing zeros make the array
        # humungous--prone to memory errors
        tz_shape = list(sig.shape)
        tz_shape[-1] *= cls.alias_factor
        trailing_zeroes = np.zeros(tz_shape)
        out = np.concatenate([out, trailing_zeroes], axis=-1)
        out = np.fft.fft(out, axis=-1)
        out = np.fft.fftshift(out, axes=-1)
        fft_tprime = np.fft.fftfreq(out.shape[-1], d=scan_obj.timestep)
        tprime = np.fft.fftfreq(out.shape[-1], d=np.abs(fft_tprime[1]-fft_tprime[0]))
        """
        # multiply by timestep, since we are integrating
        # apply frequency apodization if slitwidth is specified
        return out 
    
    @classmethod
    def check_space(cls, space):
        if space == cls.in_space:
            return True
        else:
            return False

class SLD:
    """
    squared law detector
    """
    #delta = None
    in_space = None
    out_space = 'i' # intensity, or perhaps integrated

    @classmethod
    def run(cls, scan_obj, sig):
        if cls.in_space == 'f': # freq or time works for us
            # using parseval's theorem we can verify that 
            # dw = timestep / sig.shape[-1]
            delta = scan_obj.timestep / sig.shape[-1]
        elif cls.in_space == 't': #time
            delta = scan_obj.timestep
        else:
            print 'using a delta of 1'
            delta = 1.
        # how to decide whether to use dt or dw?
        # include timestep later?
        out = (np.abs(sig.sum(axis=-2))**2).sum(axis=-1) * delta 
        return out

    @classmethod
    def check_space(cls, space):
        """
        method for checking if space can agree with calculations
        """
        # spacing is determined by 
        if space in ['f', 't']: # freq or time works for us
            cls.in_space = space
            return True
        else:
            # cannot compute absolute signal without knowing the space we are in
            print 'cannot discern whether integration space is time or frequency'
            return False

class Scatter:
    """
    add e-fields to the output coherences
    """
    ratio = 100. # ratio of lo amplitude to max sig amp
    # either list of coefficients--one scalar for each pulse--or a 
    # int
    pulse = 0 
    in_space = 't'
    out_space = 't'
    chop = False
    #windowed = False

    @classmethod
    def run(cls, scan_obj, sig):
        #--------part 1:  retrieve efields--------#
        # adding conjugate to lo will result in a real field
        lo = scan_obj.efields().real
        #if not cls.windowed:
        #    lo = scan_obj.efields().real
        #else:
        #    # lo is not of length iprime anymore (necessarily)
        #    lo = scan_obj.efields(windowed=False).real
        #    # how to phase this with signal???
        #--------part 2:  scale efield amplitudes and add together--------#
        lo_max = lo.max()
        sig_max = np.abs(sig).sum(axis=-2).max()
        try: # first assume pulse is a list-like type
            for i in range(len(cls.pulse)): 
                lo[...,i,:] *= cls.pulse[i]
            lo = lo.sum(axis=-2)
        except TypeError: lo = lo[...,cls.pulse,:]
        lo *= sig_max / lo_max * cls.ratio
        #--------part 3:  add lo to output signal--------#
        if sig.shape[-1] == lo.shape[-1]:
            if cls.chop: # block the output so baseline can be established
                # write with reference to original sig array so size of 
                # outgroups axis is maintained, even though lo doesn't have 
                # this axis
                sig *= 0.
            sig += lo[...,None,:]
        else:
            print 'sig shape does not match lo shape; aborting heterodyne'
            return sig
        #print sig.shape
        return sig

    @classmethod
    def check_space(cls, space):
        """
        method for checking if space can agree with calculations
        """
        # spacing is determined by 
        if space == 't': # lets only worry about time for now
            return True
        elif space == 'f':
            print 'heterodyne only supports time-interference at this time'
            return False
        else:
            # cannot compute absolute signal without knowing the space we are in
            print 'cannot discern whether space is time or frequency'
            return False
        pass

class LO:
    """
    temporal interferrogram of signal with an electric field
    """
    ratio = 100. # ratio of lo amplitude to max sig amp
    pulse = 0 # index of the pulse used for the heterodyne
    in_space = 't'
    out_space = 't'
    chop = False    
    
    @classmethod
    def run(cls, scan_obj, sig):
        # adding conjugate to lo will result in a real field
        lo = scan_obj.efields().real()
        lo_max = lo.max()
        sig_max = np.abs(sig).sum(axis=-2).max()
        lo *= sig_max / lo_max * cls.ratio
        from scipy import convolve
        if sig.shape[-1] == lo.shape[-1]:
            if cls.chop: # block the output so baseline can be established
                # write with reference to original sig array so size of 
                # outgroups axis is maintained, even though lo doesn't have 
                # this axis
                sig *= 0.
            sig += lo[...,cls.pulse,None,:]
    
    @classmethod
    def check_space(cls, space):
        """
        method for checking if space can agree with calculations
        """
        # spacing is determined by 
        if space == 't': # lets only worry about time for now
            return True
        elif space == 'f':
            print 'heterodyne only supports time-interference at this time'
            return False
        else:
            # cannot compute absolute signal without knowing the space we are in
            print 'cannot discern whether space is time or frequency'
            return False
        pass

# don't include this yet
del LO