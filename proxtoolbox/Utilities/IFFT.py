# Inverse Fourier transform processor.  Takes data in the fourier domain
# that has been properly normalized and arranged (i.e. negative intensities
# resulting from Python's fft algorithm flipped to the propper sign) by the
# companion routine FFT.
# First we fftshift the data, then we take the negative sign of every other
# mode in each direction (x and y), and finally we
# multiply the bad boy by (sqrt(res))^2 to get the correct intensity scaling.

from numpy.random import rand
from numpy import exp, sqrt, log, tan, pi, floor, zeros, ones
from scipy.special import gammaln
import numpy as np

def IFFT(F):

    shape = F.shape;
    res=max(shape);
    F= np.fft.fftshift(F);

    if shape[0] ==  shape[1]:
        F[1:res:2,:]=-F[1:res:2,:]
        F[:,1:res:2]=-F[:,1:res:2];
        F=res*F;
        f= np.fft.ifft2(F);
    else:
        F = F.ravel('F'); #create 1-d array, colum-major order (matlab style), not really nice
        F[1:res:2] =-F[1:res:2];
        F = F.reshape(shape[1],shape[0]).T; #back to original shape
        F=np.sqrt(res)*F;
        f= np.fft.ifft(F,axis=0);
    return f;
