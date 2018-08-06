# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 13:08:06 2015

@author: rebecca
"""

from math import sqrt
from numpy import zeros
from scipy.linalg import norm
from .algorithms import Algorithm

class AP_expert(Algorithm):
    """
    Alternating Projections
    """
    
    def __init__(self, config):
        """
        Parameters
        ----------
        config : dict        
            Dictionary containing the problem configuration.
            It must contain the following mappings:
            
                prox1: ProxOperator
                    First ProxOperator (the class, no instance)
                prox2: ProxOperator
                    Second ProxOperator (the class, no instance)
                beta0: number
                    Starting relaxation parmater
                beta_max: number
                    Maximum relaxation parameter
                beta_switch: int
                    Iteration at which beta moves from beta0 -> beta_max
                norm_data: number
                    ?
                Nx: int
                    Row-dim of the product space elements
                Ny: int
                    Column-dim of the product space elements
                Nz: int
                    Depth-dim of the product space elements
                dim: int
                    Size of the product space
        """
        self.prox1 = config['proxoperators'][0](config); self.prox2 = config['proxoperators'][1](config);
        self.norm_data = config['norm_data'];
        self.Nx = config['Nx']; self.Ny = config['Ny']; self.Nz = config['Nz'];
        self.product_space_dimension = config['product_space_dimension'];
        self.iters = 0

        if 'truth' in config:
            self.truth = config['truth']
            self.truth_dim = config['truth_dim']
            self.norm_truth = config['norm_truth']
    
    def run(self, u, tol, maxiter):
        """
        Runs the algorithm for the specified input data
        """
        prox1 = self.prox1; prox2 = self.prox2;

        if u.ndim < 3:
            p = 1
            q = 1
        elif u.ndim == 3:
            p = u.shape[2]
            q = 1
        else:
            p = u.shape[2]
            q = u.shape[3]


        norm_data = self.norm_data;
        
        iters = self.iters
        change = zeros(maxiter+1);
        change[0] = 999;
        gap = change.copy();

        if  hasattr(self, 'truth'):
            Relerrs = change.copy()
        
        tmp1 = prox2.work(u);

        while iters < maxiter and change[iters] >= tol:
            iters += 1;
            
            tmp_u = prox1.work(tmp1);
            tmp1 = prox2.work(tmp_u);
            
            tmp_change = 0; tmp_gap = 0;

            if p==1 and q==1:
                tmp_change= (norm(u-tmp_u, 'fro')/norm_data)**2
                tmp_gap = (norm(tmp1-tmp_u,'fro')/norm_data)**2
                if hasattr(self, 'truth'):
                    if self.truth_dim[0] == 1:
                      z=tmp_u[0,:]
                    elif self.truth_dim[1]==1:
                      z=tmp_u[:,0]
                    else:
                      z=tmp_u
                    Relerrs[iter] = norm(self.truth - exp(-1j*angle(trace(self.truth.T*z))) * z, 'fro')/self.norm_truth

            elif q==1:
              for j in range(self.product_space_dimension):
                  tmp_change= tmp_change+(norm(u[:,:,j]-tmp_u[:,:,j], 'fro')/norm_data)**2
                  # compute (||P_SP_Mx-P_Mx||/norm_data)^2:
                  tmp_gap = tmp_gap+(norm(tmp1[:,:,j]-tmp_u[:,:,j],'fro')/norm_data)**2
            else:
              for j in range(self.product_space_dimension):
                  for k in range(self.Nz):
                      tmp_change= tmp_change+(norm(u[:,:,k,j]-tmp_u[:,:,k,j], 'fro')/norm_data)**2
                      # compute (||P_Sx-P_Mx||/norm_data)^2:
                      tmp_gap = tmp_gap+(norm(tmp1[:,:,k,j]-tmp_u[:,:,k,j],'fro')/(norm_data))**2

            change[iters] = sqrt(tmp_change);
            gap[iters] = sqrt(tmp_gap);
            
            u = tmp_u;
            
        tmp = prox1.work(u);
        tmp2 = prox2.work(u);


        if self.Nx == 1:
            u1 = tmp[:,0]
            u2 = tmp2[:,0]
        elif self.Ny == 1:
            u1 = tmp[0,:]
            u2 = tmp2[0,:]
        elif self.Nz == 1 and tmp.ndim > 2:
            u1 = tmp[:,:,0]
            u2 = tmp2[:,:,0]
        else:
            u1 = tmp
            u2 = tmp2
        change = change[1:iters+1];
        gap = gap[1:iters+1];    
        
        return {'u1': u1, 'u2': u2, 'iter': iters, 'change': change, 'gap': gap}
