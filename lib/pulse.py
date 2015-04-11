from .misc import *

# factor used to convert FWHM to stdev for function definition
# gaussian defined such that FWHM,intensity scale = sigma * 2 * sqrt(ln(2))
# (i.e. FWHM, intensity scale ~ 1.67 * sigma, amplitude scale)
FWHM_to_sigma = 1./(2*np.sqrt(np.log(2)))

class Gauss_rwa:
    """
    returns an array of gaussian values with arguments of 
        input array vec, 
        mean mu, 
        standard deviation sigma
        amplitude amp
        frequency freq
        phase offset p
        and conjugate:  use -1 to select complex conjugate
    
    we are in rw approx, so if interaction is supposed to be negative, give
    frequency a negative sign
    """
    # go to values for arguments
    defaults = [1.0, 55., 0., 6900., 0.]
    # dictionary to associate array position to pulse attribute
    cols = {
        'A' : 0, # amplitude, a.u.
        's' : 1, # pulse FWHM (in fs)
        'd' : 2, # pulse center delay fs
        'w' : 3, # frequency in wavenumbers
        'p' : 4  # phase shift, in radians
    }
    # initial vars come from misc module, just as with scan module
    timestep = timestep
    early_buffer = early_buffer
    late_buffer = late_buffer
    
    @classmethod
    def pulse(cls, eparams, pm=None):
        """
            accepts a 2d array where the final index is the params for the 
            fields to be generated
        """
        # import if the size is right
        area  = eparams[:,0]
        sigma = eparams[:,1]
        mu    = eparams[:,2]
        freq  = eparams[:,3]
        p     = eparams[:,4]
        
        # proper unit conversions
        sigma *= FWHM_to_sigma
        freq *= wn_to_omega
        # redefine delays s.t. always relative to the first pulse listed
        offset = mu[0]
        # subtract off the value
        mu -= offset
    
        # normalize amplitdue to stdev
        y = area / (sigma*np.sqrt(2*np.pi))
        #print y, sigma, mu, freq, p
        # calculate t
        t = cls.get_t(mu)
        # incorporate complex conjugates if necessary
        if pm is None:
            cc = np.ones((eparams.shape[-1]))
        else:
            cc = np.sign(pm)
        x = rotor(cc[:,None]*(freq[:,None]*(t[None,:] - mu[:,None])+p[:,None]))
        x*= y[:,None] * np.exp(-(t[None,:] - mu[:,None])**2 / (2*sigma[:,None]**2) ) 
        return t, x
    
    @classmethod
    def get_t(cls, d):
        """
            returns the t array that is appropriate for the given pulse 
            parameters; needs the appropriate buffers/timestep as well 
        """
        t_min = min(d) - cls.early_buffer
        t_max = max(d) + cls.late_buffer
        t = np.arange(t_min, t_max, cls.timestep)
        return t
    
    @classmethod
    def export(cls):
        pass

    @classmethod
    def import_(cls):
        pass

