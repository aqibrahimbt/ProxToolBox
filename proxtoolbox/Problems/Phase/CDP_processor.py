from numpy.random import randn, random_sample
from numpy.linalg import norm
from numpy import sqrt, conj, tile, mean, exp, angle, trace
from numpy.fft import fft, ifft
import numpy as np
import time

def CDP_processor(config):
    # Implementation of the Wirtinger Flow (WF) algorithm presented in the paper 
    # "Phase Retrieval via Wirtinger Flow: Theory and Algorithms" 
    # by E. J. Candes, X. Li, and M. Soltanolkotabi
    # integrated into the ProxToolbox by 
    # Russell Luke, September 2016.

    # The input data are coded diffraction patterns about a random complex
    # valued image. 

    ## Make image
    n1 = config['Ny']
    n2 = config['Nx'] # for 1D signals, this will be 1
    x = randn(n1,n2) + 1j*randn(n1,n2)
    config['truth']=x
    config['norm_truth']=norm(x,'fro')
    config['truth_dim'] = x.shape


    ## Make masks and linear sampling operators

    L = config['product_space_dimension']                  # Number of masks 
    if n2==1:
        Masks = np.random.choice(np.array([1j, -1j, 1, -1]),(n1,L))
    elif n1==1:
        Masks = np.random.choice(np.array([1j, -1j, 1, -1]),(L,n2))
    else:
        Masks = zeros((n1,n2,L));  # Storage for L masks, each of dim n1 x n2
        # Sample phases: each symbol in alphabet {1, -1, i , -i} has equal prob. 
        for ll in range(L):
            Masks[:,:,ll] = np.random.choice(np.array([1j, -1j, 1, -1]),(n1,n2))

    # Sample magnitudes and make masks 
    temp = random_sample(Masks.shape) #works like rand but accepts tuple as argument
    Masks = Masks * ( (temp <= 0.2)*sqrt(3) + (temp > 0.2)/sqrt(2) )
    config['Masks'] = conj(Masks)
    # Saving the conjugate of the mask saves on computing the conjugate
    # every time the mapping A (below) is applied.

    if n2==1:
        # Make linear operators; A is forward map and At its scaled adjoint (At(Y)*numel(Y) is the adjoint)
        A = lambda I: fft(conj(Masks) * tile(I,[1, L]))  # Input is n x 1 signal, output is n x L array
        At = lambda Y: mean(Masks * ifft(Y), 1).reshape((n1,n2))           # Input is n x L array, output is n x 1 signal
    elif n1==1 :
        # Make linear operators; A is forward map and At its scaled adjoint (At(Y)*numel(Y) is the adjoint)
        A = lambda I: fft(conj(Masks) * tile(I,[L, 1]))  # Input is 1 x n signal, output is L x n array
        At = lambda Y: mean(Masks * ifft(Y), 0).reshape((n1,n2))            # Input is L x n array, output is 1 x n signal
    else:
        A = lambda I:  fft2(config['Masks'] * reshape(tile(I,[1, L]), (I.shape[0],I.shape[1], L)))  # Input is n1 x n2 image, output is n1 x n2 x L array
        At = lambda Y: mean(Masks * ifft2(Y), 2).reshape((n1,n2))                                           # Input is n1 x n2 X L array, output is n1 x n2 image

    # Data 
    Y = abs(A(x))
    config['rt_data']=Y
    Y=Y**2
    config['data']=Y
    config['norm_data']=sum(sum(Y))/Y.size
    normest = sqrt(config['norm_data']) # Estimate norm to scale eigenvector 
    config['norm_rt_data']=normest


    ## Initialization

    npower_iter = config['warmup_iter']
    z0 = randn(n1,n2)
    z0 = z0/norm(z0,'fro') # Initial guess 
    tic = time.time()                                     # Power iterations 
    for _ in range(npower_iter):
        z0 = At(Y*A(z0))
        z0 = z0/norm(z0)

    toc  = time.time()

    z = normest * z0                  # Apply scaling 
    if n2==1:
        Relerrs = norm(x - exp(-1j*angle(trace(x.T*z))) * z, 'fro')/norm(x,'fro')
        config['u_0'] = tile(z,[1,L])
    elif n1==1:
        Relerrs = norm(x - exp(-1j*angle(trace(z.T*x))) * z, 'fro')/norm(x,'fro')
        config['u_0'] = tilet(z,[L,1])
    else:
        Relerrs = norm(x - exp(-1j*angle(trace(x.T*z))) * z, 'fro')/norm(x,'fro')
        config['u_0']=reshape(tile(z,[1, L]), (z.shape[0], z.shape[1], L))

    print('Run time of initialization: %.2f  seconds', toc-tic)
    print('Relative error after initialization: %.2f', Relerrs)
    print('\n')