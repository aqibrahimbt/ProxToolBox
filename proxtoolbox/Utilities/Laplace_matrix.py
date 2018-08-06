# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 09:22:48 2016

@author: rebecca
"""

from numpy import array, concatenate, ones
from scipy.sparse import csr_matrix, eye, kron, spdiags

def Laplace_matrix(n):
    """
    Discrete Laplace operator with Neumann boundary conditions
    
    This method calculates the Laplacian for a 2-dim. signal (usually an image) of size n x n.
    The resulting matrix has size (n*n)x(n*n).
    
    Parameters
    ----------
    n : int - Dimension of input signal of size n x n
    
    Returns
    -------
    L : sparse array - Discrete Laplace operator of size (n*n)x(n*n)
    """
    
    dd = concatenate([-1*ones((1,n)), concatenate(([[3]], 4*ones((1,n-2)), [[3]]), axis=1), -1*ones((1,n))])
    
    D = spdiags(dd, array([-1,0,1]), n,n)
    
    E = eye(n)
    ll = spdiags(ones((2,n*n)),array([-n,n]),n*n,n*n)
    L = kron(E,D) - ll
    
    L = csr_matrix(L)
    
    for jj in range(0,n):
        L[jj,jj] = L[jj,jj]-1
        L[n*n-(n-jj),n*n-(n-jj)] = L[n*n-(n-jj),n*n-(n-jj)]-1
    
    return L
