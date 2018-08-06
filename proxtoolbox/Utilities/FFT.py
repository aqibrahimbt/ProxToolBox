# Fourier transform processor.   FFT's data in the physical domain and
# corrects for the peculiarities of numpy's fft algorithm.
# First we fft the data, then we take the negative sign of every other
# mode in each direction (x and y), then we
# divide the bad boy by (sqrt(res))^2 to get the correct intensity scaling,
# finally we fftshift everything.

from numpy.random import rand
from numpy import exp, sqrt, log, tan, pi, floor, zeros, ones
from scipy.special import gammaln
import numpy as np

def FFT(f):

    shape = f.shape;
    res=max(shape);

    if shape[0] == shape[1]:
      F= np.fft.fft2(f);
      F=F/(res);
      F[1:res:2,:]=-F[1:res:2,:];
      F[:,1:res:2]=-F[:,1:res:2];
    else:
      F = np.fft.fft(f,axis=0);
      F = F.ravel('F'); #create 1-d array, colum-major order (matlab style), not really nice
      F[1:res:2] =-F[1:res:2];
      F = F.reshape(shape[1],shape[0]).T; #back to original shape
      F=F/np.sqrt(res);

    return np.fft.fftshift(F);

