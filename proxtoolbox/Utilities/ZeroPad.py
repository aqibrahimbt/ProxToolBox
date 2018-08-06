#######################################################################
#                                                                    
#                             ZeroPad.m
#                            written by
#                           Russell Luke
#                      Nachwuchsforschergruppe
#                   Inst. fuer Num. u. Angew. Math
#                       Universitaet Goettingen
#                          August 12, 2001
# DESCRIPTION: Function for zero padding an array preparatory to
# using the FFT.  Automatically resizes the array to the next
# largest diad  by adding np.zeros symmetrically to the array and again
# symmetrically doubles the size of the array.  Takes one or two
# dimensional arrays.
#
# INPUT:
#              m = data array
#
# OUTPUT:      padded_m = the zero padded array ready for use with
#                         the FFT.
#
# USAGE:       padded_m = ZeroPad(m)
# Non-standard function calls: 
#
########################################################################
import numpy as np

def ZeroPad(m):

    dim = np.shape(m)
    orient = np.argmax(dim)
    major = dim[orient]
    minor = np.min(dim)
    tmp = np.log2(major+1)
    n = np.ceil(tmp)
    tmp = np.round(2**n-(2**tmp-1))
    if np.mod(tmp,2)==1:
      leftpad =  int((tmp+1)/2)
      rightpad = int((tmp-1)/2)
    else:
      leftpad =  int(tmp/2)
      rightpad = leftpad

    if orient==2:
      tmp = np.zeros((minor,leftpad))
      padded_m = np.concatenate((tmp,m),axis=1)
      tmp = np.zeros((minor,rightpad))
      padded_m = np.concatenate((padded_m,tmp),axis=1)
    else:
      tmp = np.zeros((leftpad,minor))
      padded_m = np.concatenate((tmp,m),axis=0)
      tmp = np.zeros((rightpad,minor))
      padded_m = np.concatenate((padded_m,tmp),axis=0)

    major = int(2**n)

    if minor!=1:
      tmp = np.log2(minor)
      tmp = np.round(2**n-2**tmp)
      if np.mod(tmp,2)==1:
        leftpad =  int((tmp+1)/2)
        rightpad = int((tmp-1)/2)
      else:
        leftpad =  int(tmp/2)
        rightpad = leftpad
      if orient==2:
        tmp = np.zeros((leftpad,major))
        padded_m = np.concatenate((tmp,padded_m),axis=0)
        tmp = np.zeros((rightpad,major))
        padded_m = np.concatenate((padded_m,tmp),axis=0)
      else:
        tmp = np.zeros((major,leftpad))
        padded_m = np.concatenate((tmp,padded_m),axis=1)
        tmp = np.zeros((major,rightpad))
        padded_m = np.concatenate((padded_m,tmp),axis=1)

    return padded_m

