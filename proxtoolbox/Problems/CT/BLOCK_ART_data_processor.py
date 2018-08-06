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
# USAGE: input = BLOCK_ART_data_processor(input) 
#
# Data loaded:  
# Nonstandard function calls:  Rbin, ZeroPad, Mgen, Resize, OnesPad
#
#################################################################################

import numpy as np

#for loading data
from pathlib import Path
import proxtoolbox.Utilities.GetData as GetData
import urllib.request
import tarfile

import scipy.io
from numpy import ceil, zeros, sqrt


def BLOCK_ART_data_processor(config):

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
    if config['fanbeam'] != 'yes':
         print('Loading data file ART_SheppLogan.mat ')
         ART = scipy.io.loadmat('../InputData/CT/ART_SheppLogan.mat')
    else:
         print('Loading data file ART_SheppLogan_fanbeam.mat ')
         ART = scipy.io.loadmat('../InputData/CT/ART_SheppLogan_fanbeam.mat')

    config['Nx']=1
    config['Nz']=1
    config['Ny']=ART['A'].shape[1]
    config['inner_dimension']=int(ART['p'][0,0])
    config['outer_dimension']=ART['theta'].size
    config['full_product_space_dimension']=config['inner_dimension']*config['outer_dimension']

    # I rearrange the rows of the matrix A and the corresponding
    # measurements b_ex so that the blocks are contiguous.
    # The projectors downstream will always presume a contiguous 
    # block structure with possibly different block sizes.  How
    # the blocks are constructed is completely up to the user and 
    # will be defined in the data processor as per this example. 
    # For graphics, the user will also have to undo the block 
    # packing in a corresponding graphics module.  
    # This makes it easier to have blocks of different sizes and 
    # different structures since the projections are relatively dumb.
    # The block sizes are saved and passed to the projector
    # so that the projector knows how to expand the blocks
    # for the microprojections, but the projectors won't know the 
    # phyical ordering of the data. 

    #----------------------------------
    #
    # Rearrange the data into blocks.
    # This is not the only way to do it and 
    # in general will be data specific.  Since this
    # example is a parallel beam CT scan, I block my data
    # into parallel beams (rows of A) at the same angle with enough 
    # space between the beams (rows of A) so that no two beams will 
    # pass through the same pixel (i.e. so that the grouped rows 
    # are orthogonal).  This way the corresponding 
    # images WITHIN THE BLOCKS will be orthogonal and hence independent.
    #  The blocks will in general be of different sizes.
    #----------------------------------

    block_step= int(ceil(config['inner_dimension']/ART['N'])) # p


    tmp_A = ART['A'].toarray()
    A = ART['A'].toarray()
    tmp_b= np.copy(ART['b_ex'])
    # block_map keeps track of how many rows are in the
    # individual blocks
    config['product_space_dimension']=int(block_step*config['outer_dimension'])
    block_map=zeros((config['product_space_dimension']+1,1), dtype=np.int_).flatten() #we will later discard the last index, flatten to avoid varnings when entries are used as indices
    # the next vector tells us how to rearrange the
    # solution x to the reordered problem so that it
    # matches the original physical representation
    pointer_map=zeros(tmp_b.size);
    block_map[0]=0

    for k in range(config['outer_dimension']):
        block_index_skip=0
        for j in range(block_step):
            block_index= np.arange(j,config['inner_dimension'],block_step)
            block_size= block_index.size
            #print(k*config['inner_dimension']+block_index_skip+np.arange(1,block_size))
            l = k*config['inner_dimension']+block_index_skip         
            #print(np.arange(l,l+block_size))
            #print(tmp_A[l:(l+block_size),:])
            tmp_A[l:l+block_size,:]= A[k*config['inner_dimension']+block_index,:]
            tmp_b[l:l+block_size]=ART['b_ex'][k*config['inner_dimension']+block_index]
            pointer_map[l:l+block_size]=(k*config['inner_dimension']+block_index)
            block_index_skip=block_size + block_index_skip
            block_map[k*block_step+j+1]=block_size+block_map[k*block_step+j]

    block_map=block_map[:-1]
    # the next is a generic scaling
    # that removes the dependence of the
    # norms from the problem dimension. 
    # More specific normalizations are up to the user. 
    config['norm_data'] = sqrt(config['Ny'])

    config['A']=tmp_A
    config['b']=tmp_b
    config['block_map']=block_map
    config['pointer_map']=pointer_map

    config['u_0']=zeros((int(ART['N'][0,0]**2),config['product_space_dimension']))
    # input.product_space_dimension=length(b_ex);

