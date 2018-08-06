#                  ART_data_processor.m
#                written on May 27, 2012 by
#                     Russell Luke
#   Inst. Fuer Numerische und Angewandte Mathematik
#                Universitaet Goettingen
#
# DESCRIPTION:  
#
# INPUT: input = a data structure
# OUTPUT: input = a data structure
#
# USAGE: input = ART_data_processor(input) 
#
# Data loaded:  
# Nonstandard function calls:  Rbin, ZeroPad, Mgen, Resize, OnesPad
#
#################################################################################

import scipy.io
from numpy import sqrt, zeros, diag
#for loading data
from pathlib import Path
import proxtoolbox.Utilities.GetData as GetData
import urllib.request
import tarfile

def ART_data_processor(config):

    my_file = Path("../InputData/CT/ART_SheppLogan.mat")
    if not(my_file.is_file()):
        print("CT input data is missing.") 
        if GetData.query_yes_no("Do you want to download the CT input data?"):
            urllib.request.urlretrieve("http://vaopt.math.uni-goettingen.de/data/CT.tar.gz","../InputData/CT.tar.gz", reporthook=GetData.dlProgress)
            print("\nExtracting data...")
            tar = tarfile.open("../InputData/CT.tar.gz", "r:gz")
            tar.extractall("../InputData/CT/")
            tar.close()
    if not(my_file.is_file()):
            print('***************************************************************************************')
            print('* Input data still missing.  Please try automatic download again or manually download *') 
            print('*    http://vaopt.math.uni-goettingen.de/data/CT.tar.gz                               *')
            print('* Save and unpack the CT.tar.gz datafile in the                                       *')
            print('*    ProxMatlab/InputData subdirectory                                                *')
            print('***************************************************************************************')

    # load toy data
    print('Loading data file ART_SheppLogan.mat ')
    
    ART = scipy.io.loadmat('../InputData/CT/ART_SheppLogan.mat')

    config['Ny']=(ART['N'][0,0])**2
    config['inner_dimension']=ART['p'][0,0]
    config['outer_dimension']=ART['theta'].size
    config['Nx']=1 # config['inner_dimension']*config['outer_dimension']
    config['Nz']=1
    config['block_step']=ART['p'][0,0] # ceil(config['inner_dimension']/N) # p
    # the next is a generic scaling
    # that removes the dependence of the 
    # norms from the problem dimension. 
    # More specific normalizations are up to the user. 
    config['norm_data'] = sqrt(config['Ny'])
    config['A']=ART['A'].toarray() #otherwise A is csc sparse scipy matrix
    config['b']=ART['b_ex'].reshape(ART['b_ex'].size) #save as 1D array
    if 'rescale' in config:
        tmp=1/(diag(config['A']  @ config['A'].T)+1e-20)
        config['A']=(diag(tmp))@config['A']
        config['b']= config['b']*tmp
        #config['b']= config['b'].reshape(config['b'].size)*tmp
        #config['b'] = config['b'].reshape((config['b'].size,1))

    config['product_space_dimension']=config['block_step']*config['outer_dimension']
    config['u_0']=zeros((ART['N'][0,0]**2,config['product_space_dimension']))
    # config['product_space_dimension']=ART['b_ex'].size

