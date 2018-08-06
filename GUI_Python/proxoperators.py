# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 13:23:26 2014

@author: stefan
"""

import numpy

__all__ = ["ProxOperator","P_diag","P_parallel","magproj"]


def magproj(u, constr):
    """
    Projection operator onto a magnitude constraint
    
    Parameters
    ----------
    u : array_like
        The function to be projected onto constr (can be complex)
    constr : array_like
        A nonnegative array that is the magnitude constraint
    
    Returns
    -------
    array_like
        The projection
    """
    modsq_u = numpy.conj(u) * u;
    denom = modsq_u+3e-30;
    denom2 = numpy.sqrt(denom);
    r_eps = (modsq_u/denom2) - constr;
    dr_eps = (denom+3e-30)/(denom*denom2);
    return (1 - (dr_eps*r_eps)) * u;




class ProxOperator:
    """
    Generic interface for prox operators
    
    Methods
    -------
    work(u)
        Apply the operator
    """

    def __init__(self):
        raise NotImplementedError,"This is just an abstract interface"
    
    def work(self,u):
        """
        Applies a prox operator to some input data
        
        Parameters
        ----------
        u : array_like
            Input data for the operator
        
        Returns
        -------
        array_like
            Result of the operation
        """
        raise NotImplementedError,"This is just an abstract interface"





class P_diag(ProxOperator):
    """
    Projection onto the diagonal of a product space
    
    See Also
    --------
        ProxOperator : Generic interface for prox operators
    """
    
    def __init__(self,config):
        """
        Parameters
        ----------
        config : dict
            Dictionary containing the problem configuration. It must contain the
            following mappings:
                'Nx':int
                    Row-dim of the product space elements
                'Ny':int
                    Column-dim of the product space elements
                'Nz':int
                    Depth-dim of the product space elements
                'dim':int
                    Size of the product space
        """
        self.n = config['Nx']; self.m = config['Ny']; self.p = config['Nz'];
        self.K = config['dim'];
    
    def work(self,u):
        """
        See Also
        --------
            ProxOperator : Generic interface for prox operators
        """
        n = self.n; m = self.m; p = self.p; K = self.K;        
        
        if m == 1:
            tmp = numpy.sum(u,axis=1,dtype=u.dtype);
        elif n == 1:
            tmp = numpy.sum(u,axis=0,dtype=u.dtype);
        elif p == 1:
            tmp = numpy.zeros((n,m),dtype=u.dtype);
            for k in range(K):
                tmp += u[:,:,k];
        else:
            tmp = numpy.zeros((n,m,p),dtype=u.dtype);
            for k in range(K):
                tmp += u[:,:,:,k];
        
        tmp /= K;
        
        if m == 1:
            return numpy.dot(tmp,numpy.ones((1,K),dtype=u.dtype));
        elif n == 1:
            return numpy.dot(numpy.ones((K,1),dtype=u.dtype),tmp);
        elif p == 1:
            u_diag = numpy.empty((n,m,K),dtype=u.dtype);
            for k in range(K):
                u_diag[:,:,k] = tmp;
            return u_diag;
        else:
            u_diag = numpy.empty((n,m,p,K),dtype=u.dtype);
            for k in range(K):
                u_diag[:,:,:,k] = tmp;
            return u_diag;





class P_parallel(ProxOperator):
    """
    Projection onto the diagonal of a product space
    
    See Also
    --------
        ProxOperator : Generic interface for prox operators
    """
    
    def __init__(self,config):
        """
        Parameters
        ----------
        config : dict
            Dictionary containing the problem configuration. It must contain the
            following mappings:
                'Nx':int
                    Row-dim of the product space elements
                'Ny':int
                    Column-dim of the product space elements
                'Nz':int
                    Depth-dim of the product space elements
                'dim':int
                    Size of the product space
                'projs':sequence of ProxOperator
                    Sequence of prox operators to be used. (The classes, no instances)
        """
        self.n = config['Nx']; self.m = config['Ny']; self.p = config['Nz'];
        self.K = config['dim'];
        self.proj = [];
        for p in config['projectors']:
            self.proj.append(p(config));
        
    def work(self,u):
        """
        See Also
        --------
            ProxOperator : Generic interface for prox operators
        """
        n = self.n; m = self.m; p = self.p; K = self.K; proj = self.proj;

        if m == 1:
            u_parallel = numpy.empty((K,n),dtype=u.dtype);
            for k in range(K):
                u_parallel[k,:] = proj[k].work(u[k,:]);
        elif n == 1:
            u_parallel = numpy.empty((m,K),dtype=u.dtype);
            for k in range(K):
                u_parallel[:,k] = proj[k].work(u[:,k]);
        elif p == 1:
            u_parallel = numpy.empty((n,m,K),dtype=u.dtype);
            for k in range(K):
                u_parallel[:,:,k] = proj[k].work(u[:,:,k]);
        else:
            u_parallel = numpy.empty((n,m,p,K),dtype=u.dtype);
            for k in range(K):
                u_parallel[:,:,:,k] = proj[k].work(u[:,:,:,k]);
        
        return u_parallel;


