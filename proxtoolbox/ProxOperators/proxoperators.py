# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 13:23:26 2014

@author: stefan

The "ProxOperators"-module contains various specific operators that do the actual calculations within the ProxToolbox.
"""

import numpy as np
from numpy import conj, dot, empty, ones, sqrt, sum, zeros, exp, nonzero, log, tile, shape, real, zeros_like
#from pyfftw.interfaces.scipy_fftpack import fft2, ifft2
from numpy.fft import fft2, ifft2, ifft

__all__ = ["P_diag","P_parallel","magproj", "Approx_P_JWST_Poisson", "P_amp", "P_SP", "Approx_PM_Gaussian", "Approx_PM_Poisson", "P_S", "P_S_real", "P_sequential_hyperplane_odd", "P_sequential_hyperplane_even", "P_parallel_hyperplane", "P_block_parallel_hyperplane", "P_block_sequential_hyperplane", "Q_Heau", "P_Amod", "Approx_P_FreFra_Poisson"]



class ProxOperator:
    """
    Generic interface for prox operators
    """

    def __init__(self, config):
        """
        Initialization method for a concrete instance
        
        Parameters
        ----------
        config : dict - Parameters to configure the algorithm
        """
        raise NotImplementedError("This is just an abstract interface")
    
    def work(self,u):
        """
        Applies a prox operator to some input data
        
        Parameters
        ----------
        u : array_like - Input data for the operator
        
        Returns
        -------
        array_like - Result of the operation
        """
        raise NotImplementedError("This is just an abstract interface")


#@profile
def magproj(u, constr):
    """
    Projection operator onto a magnitude constraint

    Parameters
    ----------
    u : array_like - The function to be projected onto constr (can be complex)
    constr : array_like - A nonnegative array that is the magnitude constraint

    Returns
    -------
    array_like - The projection
    """
    """
    # naive implementation: should work now, roughly as fast as below
    mod_u = sqrt(u.real**2+u.imag**2)
    with np.errstate(divide='ignore', invalid='ignore'):
        proj = constr/mod_u
        proj[np.isnan(proj)] = 0 #then mod_u=0 and constr=0
    proj = proj*u
    index_inf = np.isinf(proj)
    proj[index_inf] = constr[index_inf] #then mod_u=0 and constr!=0
    return  proj
    """

    """  
    Inexact, but stable implementation of magnitude projection.  
    See LukeBurkeLyon, SIREV 2002
    """
    
    eps = 3e-30

    modsq_u = u.real**2+u.imag**2 #beaware: for U * conj(U) subsequent calculations are much slower since complex (more than double computation time)
    denom = modsq_u+eps
    denom2 = sqrt(denom)
    r_eps = (modsq_u/denom2) - constr
    dr_eps = (denom+eps)/(denom*denom2)

    return (1 - (dr_eps*r_eps)) * u

    
class P_diag(ProxOperator):
    """
    Projection onto the diagonal of a product space
    
    """
    
    def __init__(self,config):
        """
        Initialization
        
        Parameters
        ----------
        config : dict - Dictionary containing the problem configuration. It must contain the following mappings:
            
                'Nx': int
                    Row-dim of the product space elements
                'Ny': int
                    Column-dim of the product space elements
                'Nz': int
                    Depth-dim of the product space elements
                'dim': int
                    Size of the product space
        """
        self.n = config['Nx']; self.m = config['Ny']; self.p = config['Nz'];
        self.K = config['product_space_dimension'];
    
    def work(self,u):
        """
        Projects the input data onto the diagonal of the product space given by its dimensions.
        
        """
        n = self.n; m = self.m; p = self.p; K = self.K;        
        
        if m == 1:
            tmp = sum(u,axis=0,dtype=u.dtype)
        elif n == 1:
            tmp = sum(u,axis=1,dtype=u.dtype)
        elif p == 1:
            tmp = zeros((n,m),dtype=u.dtype)
            for k in range(K):
                tmp += u[:,:,k]
        else:
            tmp = zeros((n,m,p),dtype=u.dtype)
            for k in range(K):
                tmp += u[:,:,:,k]
        
        tmp /= K
        
        if m == 1:
            return ones((K,1),dtype=u.dtype)@ tmp.reshape(tmp.size,1)
        elif n == 1:
            return tmp.reshape(tmp.size,1)@ ones((1,K),dtype=u.dtype)
        elif p == 1:
            u_diag = empty((n,m,K),dtype=u.dtype)
            for k in range(K):
                u_diag[:,:,k] = tmp
            return u_diag
        else:
            u_diag = empty((n,m,p,K),dtype=u.dtype)
            for k in range(K):
                u_diag[:,:,:,k] = tmp
            return u_diag





class P_parallel(ProxOperator):
    """
    Projection onto the diagonal of a product space
    
    """
    
    def __init__(self,config):
        """
        Initialization
        
        Parameters
        ----------
        config : dict - Dictionary containing the problem configuration. It must contain the following mappings:
            
                'Nx': int
                    Row-dim of the product space elements
                'Ny': int
                    Column-dim of the product space elements
                'Nz': int
                    Depth-dim of the product space elements
                'dim': int
                    Size of the product space
                'projs': sequence of ProxOperator
                    Sequence of prox operators to be used. (The classes, no instances)
        """
        self.n = config['Nx']
        self.m = config['Ny']
        self.p = config['Nz']
        self.K = config['product_space_dimension']
        self.proj = []
        for p in config['projectors']:
            self.proj.append(p(config))
        
#        for p in config['projectors']:
#            self.proj.append(globals()[p](config))
#            pp = globals()[p]
#            self.proj.append(pp(config))

        
    def work(self,u):
        """
        Sequentially applies the projections of 'self.proj' to the input data.
        
        Parameters
        ----------
        u : array_like - Input
        
        Returns
        -------
        u_parallel : array_like - Projection
        """
        n = self.n; m = self.m; p = self.p; K = self.K; proj = self.proj

        if m == 1:
            u_parallel = empty((K,n),dtype=u.dtype)
            for k in range(K):
                u_parallel[k,:] = proj[k].work(u[k,:])
        elif n == 1:
            u_parallel = empty((m,K),dtype=u.dtype)
            for k in range(K):
                u_parallel[:,k] = proj[k].work(u[:,k])
        elif p == 1:
            u_parallel = empty((n,m,K),dtype=u.dtype)
            for k in range(K):
                u_parallel[:,:,k] = proj[k].work(u[:,:,k])
        else:
            u_parallel = empty((n,m,p,K),dtype=u.dtype)
            for k in range(K):
                u_parallel[:,:,:,k] = proj[k].work(u[:,:,:,k])
        
        return u_parallel


#original matlab comment
#                      Approx_PM_Poisson.m
#             written on Feb. 18, 2011 by 
#                   Russell Luke
#   Inst. Fuer Numerische und Angewandte Mathematik
#                Universitaet Gottingen
#
# DESCRIPTION:  Projection subroutine for projecting onto Fourier
#               magnitude constraints
#
# INPUT:        input = data structure
#                       .data_ball is the regularization parameter described in 
#                        D. R. Luke, Nonlinear Analysis 75 (2012) 1531–1546.
#               u = function in the physical domain to be projected
#
# OUTPUT:       
#               
# USAGE: u_epsilon = P_M(input,u)
#
# Nonstandard Matlab function calls:  MagProj

class Approx_P_JWST_Poisson(ProxOperator):

    def __init__(self,config):
        self.TOL2 = config['TOL2'];
        self.epsilon = config['data_ball'];
        self.indicator_ampl = config['indicator_ampl'];
        self.illumination_phase = config['illumination_phase'];
        self.data_zeros = config['data_zeros'];
        self.data_sq = config['data_sq'];
        self.data = config['data'];
        self.product_space_dimension = config['product_space_dimension'];
        self.abs_illumination = config['abs_illumination'];
    

        #calculate the following once (not every iteration) for better performance (speedup ~20# for 500 Iterations JWST)
        self.exp_illu = exp(1j*self.illumination_phase)*tile(self.indicator_ampl[...,None],(1,1,self.product_space_dimension-1))
        self.exp_illu_minus = exp(-1j*self.illumination_phase)*tile(self.indicator_ampl[...,None],(1,1,self.product_space_dimension-1))

    #@profile    
    def work(self,u):
        TOL2 = self.TOL2;
        epsilon = self.epsilon;
        indicator_ampl = self.indicator_ampl;
        illumination_phase = self.illumination_phase;
        data_zeros = self.data_zeros;
        data_sq = self.data_sq;
        data = self.data;
        product_space_dimension = self.product_space_dimension;
        abs_illumination = self.abs_illumination
        exp_illu = self.exp_illu
        exp_illu_minus = self.exp_illu_minus

        u_new = zeros(u.shape,u.dtype)
        

        for j in range(product_space_dimension-1):
            U = fft2( exp_illu[:,:,j] * u[:,:,j]);
            U_sq = U.real**2+U.imag**2 #slightly faster real(U * conj(U)) and be aware: for U * conj(U) subsequent calculations much slower since complex
            tmp = U_sq/data_sq[:,:,j];
            mask = ones(tmp.shape, np.bool);
            mask[data_zeros[j]] = 0;
            tmp[mask]=1;
            U_sq[mask]=0;
            IU= tmp==0;
            tmp[IU]=1;
            tmp=log(tmp)
            hU = sum(U_sq *tmp + data_sq[:,:,j] - U_sq)
            if hU>=epsilon+TOL2:
                U0 = magproj(U,data[:,:,j]); #argument order changed compared to matlab implementation!!!
                u_new[:,:,j] = exp_illu_minus[:,:,j] *ifft2(U0);
            else:
                u_new[:,:,j] = u[:,:,j]
            #else:
            # no change

        # now project onto the pupil constraint.
        # this is a qualitative constraint, so no 
        # noise is taken into account.
        j=product_space_dimension-1;
        u_new[:,:,j] = magproj(u[:,:,j],abs_illumination);#argument order changed compared to matlab implementation!!!
        
        return u_new;

#original matlab comment
#                      P_amp.m
#             written on Feb 3, 2012 by 
#                   Russell Luke
#   Inst. Fuer Numerische und Angewandte Mathematik
#                Universitaet Gottingen
#
# DESCRIPTION:  Projection subroutine for projecting onto
#               amplitude constraints

# INPUT:        input.amplitude =  nonnegative array
#               u = complex-valued array to be projected onto the ampltude constraints
#
# OUTPUT:       p_amp    = the projection 
#               
# USAGE: p_amp = P_amp(S,u)
#
# Nonstandard Matlab function calls:  MagProj

class P_amp(ProxOperator):

    def __init__(self,config):
        self.amplitude = config['amplitude'];

    def work(self,u):
        return magproj(u,self.amplitude); #argument order changed compared to matlab implementation!!!


#                      P_SP.m
#             written on May 23, 2002 by 
#                   Russell Luke
#   Inst. Fuer Numerische und Angewandte Mathematik
#                Universitaet Gottingen
#
# DESCRIPTION:  Projection subroutine for projecting onto
#               nonnegativity and support constraints
#
# INPUT:        input, a data structure with .supp_ampl a vector of indeces of the nonzero elements of the array
#               u = function in the physical domain to be projected
#
# OUTPUT:       p_SP    = the projection IN THE PHYSICAL (time) DOMAIN
#               
# USAGE: p_SP = P_SP(input,u)
#
# Nonstandard Matlab function calls:  

class P_SP(ProxOperator):

    def __init__(self,config):
        self.support_idx = config['support_idx']

    def work(self,u):
        p_SP=zeros(shape(u), np.float64)
        p_SP[self.support_idx] = np.maximum(real(u[self.support_idx]),0)
        return p_SP

#                      Approx_PM_Poisson.m
#             written on Feb. 18, 2011 by 
#                   Russell Luke
#   Inst. Fuer Numerische und Angewandte Mathematik
#                Universitaet Gottingen
#
# DESCRIPTION:  Projection subroutine for projecting onto Fourier
#               magnitude constraints
#
# INPUT:        Func_params = a data structure with .data = nonegative real FOURIER DOMAIN CONSTRAINT
#                             .data_ball is the regularization parameter described in 
#                             D. R. Luke, Nonlinear Analysis 75 (2012) 1531–1546.
#                            .TOL2 is an extra tolerance. 
#               u = function in the physical domain to be projected
#
# OUTPUT:       p_M    = the projection IN THE PHYSICAL (time) DOMAIN
#               phat_M = projection IN THE FOURIER DOMAIN
#               
# USAGE: u_epsilon = Approx_PM_Gaussian(Func_params,u)
#
# Nonstandard Matlab function calls:  MagProj

class Approx_PM_Gaussian(ProxOperator):

    def __init__(self,config):
        self.TOL2 = config['TOL2']
        self.b = config['data'] * config['data']
        self.epsilon = config['data_ball']
        self.data = config['data']

    def work(self,u):
        TOL2 = self.TOL2;
        b = self.b
        epsilon = self.epsilon;
        U = fft2(u)
        U0 = magproj(U, self.data)
        U0_sq = U0 * conj(U0)
        tmp = U0_sq - b
        h=sum(sum(tmp * conj(tmp)))
        if h>=epsilon+TOL2:
            return ifft2(U0)
        else:
            return u

#                      Approx_PM_Poisson.m
#             written on Feb. 18, 2011 by 
#                   Russell Luke
#   Inst. Fuer Numerische und Angewandte Mathematik
#                Universitaet Gottingen
#
# DESCRIPTION:  Projection subroutine for projecting onto Fourier
#               magnitude constraints
#
# INPUT:       Func_params = a data structure with 
#                             .data = nonegative real FOURIER DOMAIN CONSTRAINT
#                             .data_sq = the elementwise square of .data
#                             .data_ball is the regularization parameter described in 
#                             D. R. Luke, Nonlinear Analysis 75 (2012) 1531–1546.
#                            .TOL2 is an extra tolerance. 
#               u = function in the physical domain to be projected
#
# OUTPUT:       
#               
# USAGE: u_epsilon = P_M(Func_params,u)
#
# Nonstandard Matlab function calls:  MagProj

class Approx_PM_Poisson(ProxOperator):

    def __init__(self,config):
        self.TOL2 = config['TOL2']
        self.M = config['data']
        self.b = config['data_sq']
        self.Ib = config['data_zeros']
        self.epsilon = config['data_ball']

    
    def work(self,u):
        TOL2 = self.TOL2
        M = self.M
        b = self.b
        Ib = self.Ib
        epsilon = self.epsilon

        U = fft2(u)
        U_sq = U * conj(U)
        tmp = U_sq/b; tmp[Ib] = 1
        U_sq[Ib]=0
        IU = tmp==0
        tmp[IU]=1
        tmp=log(tmp)
        hU = sum(sum(U_sq * tmp + b - U_sq))
        if hU>=epsilon+TOL2:
            U0 = magproj(U,M)
            return ifft2(U0)
        else:
            return u

class P_S(ProxOperator):

    def __init__(self,config):
        self.support_idx = config['support_idx']

    def work(self,u):
        p_S=zeros(shape(u),dtype = np.complex128)
        p_S[self.support_idx] = u[self.support_idx]
        return p_S
        
class P_S_real(ProxOperator):

    def __init__(self,config):
        self.support_idx = config['support_idx']

    def work(self,u):
        p_S=zeros(shape(u),dtype = np.complex128)
        p_S[self.support_idx] = real(u[self.support_idx])
        return p_S

#  hyperplane projection subroutine of the
#  ART method  
#  for tomographic recostruction of a 
#  density profile.
# 
#  Russell Luke
#  Universitaet Goettingen
#  May 29, 2012
#
class P_sequential_hyperplane_odd(ProxOperator):

    def __init__(self,config):
        self.m = config['A'].shape[0]
        self.A = config['A']
        self.b = config['b']
        self.config = config #this is not a copy

    def work(self,u):
        K=np.mod(2*self.config['iter'],self.m)+1
        a=self.A[K,:]
        a=a.reshape((1,a.size)) #numpy doesn't differentiate between row and colum vectors for 1D
        aaT=(a @ a.T)+1e-20
        b=self.b[K]/aaT
        v=a @ u /aaT;
        v=u+a.T@(b-v);
        return v

class P_sequential_hyperplane_even(ProxOperator):

    def __init__(self,config):
        self.m = config['A'].shape[0]
        self.A = config['A']
        self.b = config['b']
        self.config = config #this is not a copy

    def work(self,u):
        K=np.mod(2*self.config['iter'],self.m)+1
        a=self.A[K,:]
        a=a.reshape((1,a.size)) #numpy doesn't differentiate between row and colum vectors for 1D
        aaT=(a @ a.T)+1e-20
        b=self.b[K]/aaT
        v=a @ u /aaT;
        v=u+a.T@(b-v);
        return v

#  hyperplane projection subroutine of the
#  Cimmino ART method  
#  for tomographic recostruction of a 
#  density profile.
# 
#  Russell Luke
#  Universitaet Goettingen
#  May 29, 2012
#

class P_parallel_hyperplane(ProxOperator):

    def __init__(self,config):
        self.A = config['A']
        self.b = config['b']

    def work(self,u):
        A = self.A
        s = A.shape
        v = zeros((s[1],s[0]))
        for K in range(s[0]):
            a=A[K,:]
            aaT=dot(a,a)+1e-20
            b=self.b[K]
            v_scal=dot(a,u[:,K])
            u[:, K]+a.T
            v[:,K]=u[:, K]+a.T*((b-v_scal)/aaT)
        return v

#  hyperplane projection subroutine of the
#  Cimmino ART method  
#  for tomographic recostruction of a 
#  density profile.
# 
#  Russell Luke
#  Universitaet Goettingen
#  May 29, 2012
#
class P_block_parallel_hyperplane(ProxOperator):

    def __init__(self,config):
        self.Ny = config['Ny']
        self.full_product_space_dimension = config['full_product_space_dimension']
        self.product_space_dimension = config['product_space_dimension']
        self.A = config['A']
        self.b = config['b']
        self.block_map = config['block_map']

    def work(self,u0):
        v0= zeros_like(u0)
        u=  zeros((self.Ny,self.full_product_space_dimension))

        # expand the blocks:
        for j in range (self.product_space_dimension-1):
            u[:,self.block_map[j]:self.block_map[j+1]]= u0[:,j].reshape(shape(u0)[0],1) @ ones((1,self.block_map[j+1]-self.block_map[j]), dtype=np.float64)
 
        u[:,self.block_map[self.product_space_dimension-1]:] = u0[:,-1].reshape(shape(u0)[0],1) @ ones((1,self.full_product_space_dimension-self.block_map[-1]), dtype=np.float64)

        # compute the micro-projections

        AAT= (self.A*self.A).sum(-1) + 1e-20        #this is the same as np.diagonal(self.A @ self.A.T)+1e-20, but now algorithm is 100 times faster
        AAT = AAT.reshape((AAT.size,1))
        b= self.b/AAT
        v= (self.A*u.T).sum(-1)      #this is the same as np.diagonal(self.A @ self.A.T)+1e-20, but now algorithm is 100 times faster
        v= v.reshape(v.size,1)/AAT
        v=u+self.A.T *( ones((self.Ny,1), dtype=np.float64) @ (b-v).T)

        # now average and condense the blocks
        for j in range(self.product_space_dimension-1):
          v0[:,j]= sum(v[:,self.block_map[j]:self.block_map[j+1]],1)/(self.block_map[j+1]-self.block_map[j])
        
        v0[:,-1]= sum(v[:,self.block_map[self.product_space_dimension-1]:],1)/(v.shape[1]-self.block_map[self.product_space_dimension-1]-1)
        # should probably iterate the above proceedure until 'convergence',
        # but one step is close enough if the blocks consist of orthogonal 
        # elemenets...

        return v0

#  hyperplane projection subroutine of the
#  Cimmino ART method  
#  for tomographic recostruction of a 
#  density profile.
# 
#  Russell Luke
#  Universitaet Goettingen
#  May 29, 2012
#
class P_block_sequential_hyperplane(ProxOperator):

    def __init__(self,config):
        self.Ny = config['Ny']
        self.full_product_space_dimension = config['full_product_space_dimension']
        self.product_space_dimension = config['product_space_dimension']
        self.A = config['A']
        self.b = config['b']
        self.block_map = config['block_map']

    def work(self,u0):
    # compute the sequential projections on the blocks
    # this could easily be parallelized since the blocks are
    # independent.  I only do one cycle of this, but should 
    # iterate until 'convergence'
        v0 = np.copy(u0)
        for j in range(self.product_space_dimension-1):
          for k in  range(self.block_map[j],self.block_map[j+1]):
              a=self.A[k,:]
              aaT=(a @ a.T)+1e-20
              b=self.b[k]/aaT
              v=a @ v0[:,j]/aaT
              v0[:,j]=v0[:,j]+a.T * (b-v)

        for k in range(self.block_map[self.product_space_dimension-1],self.full_product_space_dimension):
              a=self.A[k,:]
              aaT=(a @ a.T)+1e-20
              b=self.b[k]/aaT
              v=a @ v0[:,-1]/aaT
              v0[:,-1]=v0[:,-1]+a.T * (b-v) 

        # should probably iterate the above proceedure until 'convergence',
        # but one step is close enough if the blocks consist of orthogonal 
        # elemenets...
        return v0

#                      Q_Heau.m
#             written on MAY 31, 2006 by 
#           Russell Luke & Daniel Gempesaw
#             University of Delaware
#
# DESCRIPTION:  Q operator for the Heaugaseau-like algorithm as
#               described in Bauschke Combettes&Luke Jour. Approx. Thry., 2006
#
# INPUT:        x,y,z - all column vectors, not nec. real
#
# OUTPUT:       y_new = Q(x,y,z)
#               
# USAGE: y_new = Q_Heau(x,y,z)
#
# Nonstandard Matlab function calls:  none

from numpy import real, imag

class Q_Heau(ProxOperator):

    def work(x,y,z):

        xsi = real((y-z).T @ (x-y))+imag((y-z).T @ (x-y))
        eta = (x-y).T @ (x-y)
        nu =  (y-z).T  @ (y-z)
        rho = eta*nu-xsi**2

        if((rho<=0) and (xsi>=0)):
            return z
        elif((rho>0) and (xsi*nu>=rho)):
            return x+(1+xsi/nu)*(z-y)
        elif((rho>0) and (xsi*nu<rho)):
            return y+(nu/rho)*(xsi*(x-y)+nu*(z-y))


#                      P_Amod.m
#             written on June 1, 2017 by 
#                   Russell Luke
#   Inst. Fuer Numerische und Angewandte Mathematik
#                Universitaet Gottingen
#
# DESCRIPTION:  Projection subroutine for projecting onto
#               support constraints
#
# INPUT:        input.supp_ampl = OBJECT (TIME) DOMAIN CONSTRAINT:  indeces of the 
#               nonzero elements in the object domain.
#               u = function in the physical domain to be projected
#
# OUTPUT:       p_A    = the projection IN THE PHYSICAL (time) DOMAIN,
#                        leaves the phase alone on the support, resets
#                        amplitude to 1 and sets amplitude to 1 with 0
#                        phase outside support.
#               
# USAGE: p_A = P_Amod(input,u)
#
# Nonstandard Matlab function calls:  

class P_Amod(ProxOperator):

    def __init__(self,config):
        self.support_idx = config['support_idx']

    def work(self,u):
        p_A= np.ones_like(u)
        p_A[self.support_idx] = u[self.support_idx]/abs(u[self.support_idx])
        return(p_A)


#                Approx_P_FreFra_Poisson.m
#                written around Jan 23, 2012 by 
#                    Robert Hesse and Russell Luke
#   Inst. Fuer Numerische und Angewandte Mathematik
#                Universitaet Gottingen
#
# DESCRIPTION:  Projection subroutine for projecting onto Fourier
#               magnitude constraints.  This is an approximate
#               projector onto a ball around the data determined
#               by the Kullback-Leibler divergence, as appropriate
#               for Poisson noise.  The analysis of such 
#               approximate projections was studied in 
#               D.R.Luke, ``Local Linear Convergence of Approximate 
#               Projections onto Regularized Sets '', Nonlinear Analysis, 
#               75(2012):1531--1546.
#
# INPUT:        input, a data structure
#               f = function in the physical domain to be projected
#
# OUTPUT:       p_M    = the projection IN THE PHYSICAL (time) DOMAIN
#               phat_M = projection IN THE FOURIER DOMAIN
#               
# USAGE: p_M = Approx_P_FreFra_Poisson(M,u)
#
# Nonstandard Matlab function calls:  MagProj

class Approx_P_FreFra_Poisson(ProxOperator):

    def __init__(self,config):
        self.fresnel_nr = config['fresnel_nr']
        self.Nx = config['Nx']
        self.FT_conv_kernel = config['FT_conv_kernel']
        self.use_farfield_formula = config['use_farfield_formula']
        if self.use_farfield_formula:
            self.illumination = config['abs_illumination']
        self.magn = config['magn']
        self.data_sq = config['data_sq']
        self.data = config['data']
        self.data_zeros = config['data_zeros']
        if 'beam' in config:
            self.beam = config['beam']
        self.data_ball = config['data_ball']
        self.TOL2 = config['TOL2']

    def work(self,f):
     
        if self.use_farfield_formula:
            if  self.fresnel_nr>0:       
                fhat = -1j*self.fresnel_nr/(self.Nx**2*2*pi) * ifftshift(fft2(fftshift(f-self.illumination)))+ ifftshift(self.FT_conv_kernel)
            else:
                fhat = (-1j/self.Nx**2)*ifftshift(fft2(fftshift(f)))
        else:
            if hasattr(self,'beam'):
               fhat = ifft2(self.FT_conv_kernel*fft2(f*self.beam))/self.magn
            else:
               fhat = ifft2(self.FT_conv_kernel*fft2(f))/self.magn

        fhat_sq = fhat.real**2+fhat.imag**2
        tmp = fhat_sq/self.data_sq; tmp[self.data_zeros]=1;
        fhat_sq[self.data_zeros]=0
        Ifhat= tmp==0
        tmp[Ifhat]=1
        tmp=log(tmp)
        hfhat = sum(fhat_sq*tmp + self.data_sq - fhat_sq)
        if hfhat>=self.data_ball+self.TOL2:
            p_Mhat = magproj(fhat,self.data)
            if self.use_farfield_formula:
                if  self.fresnel_nr>0:   
                    p_M=1j/self.fresnel_nr*(self.Nx**2*2*pi) *ifftshift(ifft2(fftshift(p_Mhat)))
                else:
                    p_M=1j*(self.Nx**2) *ifftshift(ifft2(fftshift(p_Mhat)))
            else:
                if hasattr(self,'beam'):
                    p_M=ifft2(fft2(p_Mhat*self.magn)/self.FT_conv_kernel)/self.beam
                else:
                    p_M=ifft2(fft2(p_Mhat*self.magn)/self.FT_conv_kernel)
        else:
            p_M = f
       
        return(p_M)

