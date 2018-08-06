#                  CDI_data_processor.m
#                written on Jan 10, 2012 by
#                     Russell Luke
#   Inst. Fuer Numerische und Angewandte Mathematik
#                Universitaet Goettingen
#
# DESCRIPTION:  Data processing subroutine called by
#               Loads data and prepares it for
#               phase retrieval algorithms.
#
# INPUT: all the input is set in Phase_in
#               data_file    = character string identifying the data
#                              file to load
#               step_up      = number of dyads to increase resolution.
#               noise        = numeric toggle for adding noise to original
#                              image
#               data_ball          = noise parameter
#               support_type = character string for type of
#                              support. Either 'autocorr', 'box' or 'autocorr_box'
#               tightness    = numeric parameter adjusts the tightness
#                              of the autocorrelation support constraint.
#
# OUTPUT:       u_0 = the algorithm initial guess (I use a bad initial guess on
#                                                   purpose)
#               S   = the object domain constraints.  This is either a 
#                        structure or an array.
#               Func_params   = image domain measurements and paramters.
#                                This is a data structure.
#
# USAGE: [S,M,u_0,true_object,success] =
#    Data_processor(data_file,program,step_up, noise,data_ball,support_type,tightness)
#
# Data loaded:  CDI data
# Nonstandard function calls:  none
#
#################################################################################

from scipy.io import loadmat
from matplotlib.pyplot import colorbar, show, imshow
from numpy import log10, log2, ceil, floor, argmax, unravel_index, zeros, where, nonzero, pi, sqrt, ones, float64
from numpy.fft import fftshift, ifft2
#from pyfftw.interfaces.scipy_fftpack import fftshift, ifft2
from numpy.linalg import norm
from numpy.random import rand

def CDI_data_processor(config):

    if config['distance'] =='near field':   #moved here from demo file in matlab
            config['fresnel_nr'] = 1*2*np.pi*config['Nx']
    else:
            config['fresnel_nr'] = 0   

    data_file = config['data_filename']
    data_ball = config['data_ball']


    f = loadmat('../InputData/Phase/CDI_intensity.mat')
    # diffraction pattern
    dp = f['intensity']['img'][0,0]

    im = imshow(log10(dp))
    colorbar(im)
    show(block=False)

    orig_res = dp.size # actual data size
    step_up = ceil(log2(config['Nx'])-log2(orig_res))
    workres = 2**(step_up)*2**(floor(log2(orig_res)))#desired array size

    N=int(workres)

    ## center data and use region of interest around center
    #central pixel
    #find max data value
    argmx = unravel_index(argmax(dp), dp.shape)

    Xi = argmx[0]+1
    Xj = argmx[1]+1

    #get desired roi:
    #necessary conversion
    Di = N/2-(Xi-1)
    Dj = N/2-(Xj-1)

    # roi around central pixel
    Ni = 2*(Xi+Di*(Di<0)-1)
    Nj = 2*(Xj+Dj*(Dj<0)-1)

    tmp = zeros((N,N))

    tmp[int(Di*(Di>0)):int(Di*(Di>0)+Ni),int(Dj*(Dj>0)):int(Dj*(Dj>0)+Nj)] = dp[int(Xi-Ni/2)-1:int(Xi+Ni/2-1),int(Xj-Nj/2)-1:int(Xj+Nj/2-1)]

    M=(fftshift(tmp))**(.5)
    # M=tmp.^(.5)
    ## define support
    DX = ceil(N/7)
    S = zeros((N,N))
    S[int(N/4+1+DX)-1:int(3/4*N-1-DX),int(N/4+1+DX)-1:int(3/4*N-1-DX)] = 1

    # config['data_ball=config['Nx*config['Ny*data_ball
    config['rt_data']=M
    # standard for the main program is that 
    # the data field is the magnitude SQUARED
    # in Luke.m this is changed to the magnitude.
    config['data']=M**2
    config['norm_rt_data']=norm(ifft2(M),'fro')
    config['data_zeros'] = where(M==0)
    config['support_idx'] = nonzero(S)
    config['product_space_dimension'] = 1

    # use the abs_illumination field to represent the 
    # support constraint.
    config['abs_illumination'] = S 
    config['supp_phase'] = []

    # initial guess
    config['u_0'] = 0.5*S*ones((N,N),float64)#S*rand(N,N)
    config['u_0'] = config['u_0']/norm(config['u_0'],'fro')*config['norm_rt_data'] 

    if config['fresnel_nr']*config['magn'] <=  2*pi*sqrt(config['Nx']*config['Ny']):
        config['use_farfield_formula'] = 1
        print('Using farfield formula.')
    else:
        config['use_farfield_formula'] = 0
        print('Using nearfield formula.')

