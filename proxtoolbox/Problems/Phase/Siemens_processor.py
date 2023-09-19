# data reader/processor

import numpy as np
from numpy.linalg import norm
from numpy.fft import fft2, ifft2
from numpy.random import rand
import proxtoolbox.Utilities as Utilities
from pathlib import Path
import os
import proxtoolbox.Utilities.GetData as GetData
import urllib.request
import tarfile

def Siemens_processor(config):

    data_ball = config['data_ball']

    noise = config['noise']
    snr = config['data_ball']

    my_file = Path("../InputData/Phase/Siemens_star_200px.mat")
    if not(my_file.is_file()):
        print("Phase input data is missing.") 
        if GetData.query_yes_no("Do you want to download the phase input data?"):
            urllib.request.urlretrieve("http://num.math.uni-goettingen.de/~r.luke/tmp/Phase.tar.gz","../InputData/Phase.tar.gz", reporthook=GetData.dlProgress)
            print("\nExtracting data...")
            tar = tarfile.open("../InputData/Phase.tar.gz", "r:gz")
            tar.extractall("../InputData/")
            tar.close()
    if not(my_file.is_file()):
            print('***************************************************************************************')
            print('* Input data still missing.  Please try automatic download again or manually download *') 
            print('*    http://num.math.uni-goettingen.de/data/Phase.tar.gz                              *')
            print('* Save and unpack the Phase.tar.gz datafile in the                                    *')
            print('*    ProxMatlab/InputData subdirectory                                                *')
            print('***************************************************************************************')

    print('Loading data file: Siemens_star_200px.mat')
    S = np.loadtxt('../InputData/Phase/Siemens_star_200px.mat')

    S=Utilities.ZeroPad(S)
    dim =np.shape(S)
    config['Nx'] = dim[1]
    config['Ny'] = dim[0]

    if config['object'] in ['real', 'nonnegative']:
        M=abs(fft2(S))
        # config['data_ball']=config['Nx']*config['Ny']*data_ball
        config['rt_data']=M
        # standard for the main program is that
        # the data field is the magnitude SQUARED
        # in Luke.m this is changed to the magnitude.
        config['data']=M**2
        config['norm_rt_data']=norm(S,'fro')
        config['norm_data']=config['norm_rt_data']**2
        config['data_zeros'] = np.where(M==0)
        # below, we make the support too small to increase the inconsistency
        # Stmp=zeros(size(S))
        # Stmp((m/2-m/4):(m/2+m/4),(n/2-n/4):(n/2+n/4))=S((m/2-m/4):(m/2+m/4),(n/2-n/4):(n/2+n/4))
        # S=Stmp
        config['support_idx'] = np.nonzero(S)

        # use the abs_illumination field to represent the
        # support constraint.
        config['amplitude'] = config['norm_rt_data']*S/(norm(S,'fro'))
        config['amplitude'] = config['amplitude']/norm(config['amplitude'],'fro')*config['norm_rt_data']
        config['supp_phase'] = np.nonzero(np.ones((config['Ny'],config['Nx'])))
        config['illumination_phase'] = np.nonzero(np.ones((config['Ny'],config['Nx'])))

    elif config['object'] == 'complex':
        # put some phase across S
        points = config['Nx']
        config['norm_rt_data']=norm(S,'fro')
        config['norm_data']=config['norm_rt_data']**2
        # use the amplitude field to represent the
        # support constraint.
        config['amplitude'] = config['norm_rt_data']*S/(norm(S,'fro'))
        # input.amplitude = input.amplitude/norm(input.amplitude,'fro')*input.norm_rt_data(1)
        G=np.zeros(S.shape)# Gaussian(points,10,[points/2+1,points/2+1])
        W=np.matmul(S,np.exp(1j*2*np.pi*G))
        M=abs(fft2(W))
        # config['data_ball']=config['Nx']*config['Ny']*data_ball
        config['rt_data']=M
        # standard for the main program is that
        # the data field is the magnitude SQUARED
        # in Luke.m this is changed to the magnitude.
        config['data']=M**2
        config['data_zeros'] = np.where(M==0)
        config['support_idx'] = np.nonzero(S)    

    elif config['object'] == 'phase':
        # put some phase across S
        points = config['Nx']
        # use the abs_illumination field to represent the
        # support constraint.
        config['abs_illumination'] = ones(S.shape)  
        G=Gaussian(points,10,[points/2+1,points/2+1])
        W=np.exp((1j*2*pi)*S*G)
        config['norm_rt_data']=norm(W,'fro')
        M=abs(fft2(W))
        # config['data_ball=config['Nx*config['Ny*data_ball
        config['rt_data']=M
        config['norm_data']=config['norm_rt_data']**2
        # standard for the main program is that
        # the data field is the magnitude SQUARED
        # in Luke.m this is changed to the magnitude.
        config['data']=M**2
        config['data_zeros'] = np.where(M==0)
        config['support_idx'] = np.nonzero(W)    
        config['supp_phase'] = S
        config['illumination_phase'] = S    

    # initial guess
    config['u_0'] = ifft2(M*np.exp(2*np.pi*1j*0.5*np.ones((dim[0],dim[1]))))# rand(dim[0],dim[1])))
    # config['u_0 = ifft2(M.*exp(2*pi*1i*Gaussian(N,N/2,[N/2+1,N/2+1]))).*S
    # config['u_0 = config['u_0/norm(config['u_0,'fro')*config['abs_illumination 

    config['product_space_dimension'] = 1
    config['truth'] = S
    config['truth_dim'] = np.shape(S)
    config['norm_truth']=norm(S,'fro')

    



