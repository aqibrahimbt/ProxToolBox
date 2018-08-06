'''
                  JWST_data_processor.m
                written on May 27, 2012 by
                     Russell Luke
   Inst. Fuer Numerische und Angewandte Mathematik
                Universitaet Goettingen

 DESCRIPTION:  

 INPUT: input = a data structure
 OUTPUT: input = a data structure

 USAGE: input = JWST_data_processor(input) 

 Data loaded:  
 Nonstandard function calls:  Rbin, ZeroPad, Mgen, Resize, OnesPad



 data reader/processor
'''

import numpy as np
from numpy import fromfile, exp, nonzero, zeros, pi, resize
from numpy.random import rand
from numpy.linalg import norm
from numpy.fft import fftshift
import proxtoolbox.Utilities as Utilities

#for loading data
from pathlib import Path
import proxtoolbox.Utilities.GetData as GetData
import urllib.request
import tarfile

def JWST_data_processor(config):

    my_file = Path("../InputData/Phase/pupil.pmod")
    if not(my_file.is_file()):
        print("Phase input data is missing.") 
        if GetData.query_yes_no("Do you want to download the phase input data?"):
            urllib.request.urlretrieve(" http://vaopt.math.uni-goettingen.de/data/Phase.tar.gz","../InputData/Phase.tar.gz", reporthook=GetData.dlProgress)
            print("\nExtracting data...")
            tar = tarfile.open("../InputData/Phase.tar.gz", "r:gz")
            tar.extractall("../InputData/Phase")
            tar.close()
    if not(my_file.is_file()):
            print('***************************************************************************************')
            print('* Input data still missing.  Please try automatic download again or manually download *') 
            print('*     http://vaopt.math.uni-goettingen.de/data/Phase.tar.gz                           *')
            print('* Save and unpack the Phase.tar.gz datafile in the                                    *')
            print('*    ProxMatlab/InputData subdirectory                                                *')
            print('***************************************************************************************')

    #moved here from JWST_in since if statements not possible in dictonary
    if 'distance' in config:
        if config['distance'] =='near field':
            config['fresnel_nr'] = 1*2*np.pi*config['Nx']
            config['use_farfield_formula'] = 0
        else:
            config['fresnel_nr'] = 0
            config['use_farfield_formula'] = 1   


    data_ball = config['data_ball'];

    newres = config['Nx'];
    noise = config['noise'];
    snr = config['data_ball'];

    print('Loading data')
    with open('../InputData/Phase/pupil.pmod','r') as fid:
    # read lower endian float <f
        Xi_A = np.fromfile(fid, dtype='<f');
        Xi_A = Xi_A.astype(np.float64)
        Xi_A = Xi_A.reshape((512,512)).T;

    diversity=3;

    with open('../InputData/Phase/phase_p37.pmod','r') as fid:
    # read lower endian float <f
        temp1 = np.fromfile(fid, dtype='<f');
        temp1 = temp1.astype(np.float64)
        temp1 = temp1.reshape(512,512).T;

    with open('../InputData/Phase/phase_m37.pmod','r') as fid:
    # read lower endian float <f
        temp2 = np.fromfile(fid, dtype='<f');
        temp2 = temp2.astype(np.float64)
        temp2 = temp2.reshape(512,512).T;

    defocus=(temp1-temp2)/2;
    theta=(temp1+temp2)/2;

    if newres != 512:
      defocus = Utilities.Resize(defocus,newres,newres);
      theta = Utilities.Resize(theta,newres,newres);
      Xi_A = Utilities.Resize(Xi_A,newres,newres);

    aberration = np.array( [np.zeros((newres,newres)), defocus, -defocus]);
    aberration = aberration.transpose((1,2,0));
 
   # aberration[:,:,1]=np.zeros((newres,newres));  # Note this order!!!!!
   # aberration[:,:,2]=defocus;  # Note this order!!!!!
   # aberration[:,:,3]=-defocus; # Note this order!!!!!

    true_object = np.multiply(Xi_A,exp(1j*(theta)));
    k= np.zeros((newres,newres,diversity));
    rt_k= np.zeros((newres,newres,diversity));

    for j in range(diversity):
      Pj = Xi_A * exp(1j*(aberration[:,:,j]+theta));
      k[:,:,j] = np.square(abs(Utilities.FFT(Pj)));

    
    epsilon=0;
    norm_rt_data=np.sqrt(sum(sum(Xi_A)));
    data_zeros= [[],[],[]];
    for i in range(diversity):
      rt_k[:,:,i]= np.sqrt(k[:,:,i]); # *newres ?
      # normalize the data so that the corresponding 
      # pupil function has amplitude 1 everywhere
      rt_k[:,:,i] = rt_k[:,:,i]/norm(Utilities.IFFT(rt_k[:,:,i]))*norm_rt_data;
      k[:,:,i] = np.square(rt_k[:,:,i]);
      if config['noise'] == 'Poisson':
          # Add Poisson noise according to Ning Lei:
          # f = alpha*4*4*pi/q/q*abs(cos(qx(iX))*sin(qy(iY))*(sin(q)/q-cos(q)));
          #        f2 = PoissonRan(f*f*2)/alpha/alpha/2;   # 2 is for I(-q)=I(q)
          #	sigma(iX, iY, iZ) = f/np.sqrt(2)/alpha/alpha/abs(ff)/abs(ff);
          # July 15, 2010:  Since this is meant to model photon counts (which
          # are integers) we add noise and then round down
          k[:,:,i] = k[:,:,i]/snr;

          for ii in range(newres):
              for jj in range(newres):
                  k[ii,jj,i]= np.random.poisson(k[ii,jj,i])*snr #use built in numpy possion instead of Utilities.PoissonRan(k[ii,jj,i])*snr;

          k[:,:,i]=np.round(k[:,:,i]);
          rt_k[:,:,i]= np.sqrt(k[:,:,i]);

      # fftshift and scale the data for use with fft:
      rt_k[:,:,i] = fftshift(rt_k[:,:,i])*newres;
      k[:,:,i] = fftshift(k[:,:,i])*newres;
      data_zeros[i] = rt_k[:,:,i].nonzero() #np.where(rt_k[:,:,i]==0)[1];
      #data_zeros[range(temp.size),i] = temp;

    config['rt_data'] = rt_k;
    config['data'] = k;
    config['norm_rt_data'] =norm_rt_data; #  this is correct since
                                     # is is calculated in the 
                                     # object domain
    config['data_zeros'] = data_zeros;
    config['supp_ampl'] = nonzero(Xi_A);
    config['indicator_ampl'] = zeros(Xi_A.shape);
    (config['indicator_ampl'])[config['supp_ampl']] = 1;

    config['product_space_dimension'] = diversity+1;
    # config.Nz = 1; # already set in the config file, but just in case

    config['abs_illumination'] = Xi_A;
    # Normalize the illumination so that it has the 
    # same energy as the data.
    config['abs_illumination'] = config['abs_illumination']/norm(config['abs_illumination'])*config['norm_rt_data'];#*config['norm_rt_data'][0]; 
    config['illumination_phase'] = aberration;

    # initial guess
    # config.u_0 = Xi_A.*exp(1i*2*pi*rand(newres));
    config['u_0']= zeros((newres,newres,config['product_space_dimension']),dtype = np.complex128);
    for j in range(config['product_space_dimension']):
        config['u_0'][:,:,j] = np.multiply(config['abs_illumination'],exp(1j*2*pi*0.5*np.ones(newres)));#rand(newres))); 
    config['truth'] = true_object;
    config['truth_dim'] = true_object.shape;
    config['norm_truth'] = norm(true_object);


