# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 12:59:44 2015

@author: rebecca
"""

from math import exp, sqrt
from numpy import zeros
from scipy.linalg import norm
from .algorithms import Algorithm

class HPR(Algorithm):
    """
    Heinz-Patrick-Russell algorithm
    
    References
    ----------
    - TODO
    
    Bauschke, Combettes & Luke, Journal of the Optical Society of America,
    20(6):1025-1034, 2003
    """
    
    def __init__(self, config):
        """
        Parameters
        ----------
        config : dict        
            Dictionary containing the problem configuration. It must contain the
            following mappings:
            
                proj1: ProxOperator
                    First ProxOperator (the class, no instance)
                proj2: ProxOperator
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
        self.proj1 = config['proj1'](config);
        self.proj2 = config['proj2'](config);
        self.norm_data = config['norm_data'];
        self.beta0 = config['beta0'];
        self.beta_max = config['beta_max'];
        self.beta_switch = config['beta_switch'];
        self.Nx = config['Nx']; self.Ny = config['Ny']; self.Nz = config['Nz'];
        self.dim = config['dim'];
        self.iters = 0
    
    def run(self, u, tol, maxiter):
        """
        Runs the algorithm for the specified input data
        """
        
        ##### PREPROCESSING
        beta_max = self.beta_max;
        beta_switch = self.beta_switch;
        proj1 = self.proj1; proj2 = self.proj2;
        norm_data = self.norm_data;
        
        beta = self.beta0;
        iters = self.iters
        change = zeros(maxiter+1,dtype=u.dtype);
        change[0] = 999;
        gap = change.copy();
        
        tmp = proj2.work(u);
        
        ##### LOOP
        while iters < maxiter and change[iters] >= tol:
            a = exp((-iters/beta_switch)**3);
            beta = (a*self.beta0) + ((1-a)*beta_max);
            iters += 1;
            
            tmp2 = 2*proj2.work(u) - u + (beta-1)*tmp;
            tmp_u = (2*proj1.work(tmp2) - tmp2 + u + (1-beta)*tmp)/2;
            tmp = proj2.work(tmp_u);
            tmp3 = proj1.work(tmp);
            
            tmp_change = 0; tmp_gap = 0;
            if self.Ny == 1 or self.Nx == 1:
                tmp_change = (norm(u-tmp_u,'fro')/norm_data)**2;
                tmp_gap = (norm(tmp3-tmp2,'fro')/norm_data)**2;
            elif self.Nz == 1:
                for j in range(self.dim):
                    tmp_change += (norm(u[:,:,j]-tmp_u[:,:,j],'fro')/norm_data)**2;
                    tmp_gap += (norm(tmp3[:,:,j]-tmp2[:,:,j])/norm_data,'fro')**2;
            else:
                for j in range(self.dim):
                    for k in range(self.Nz):
                        tmp_change += (norm(u[:,:,k,j]-tmp_u[:,:,k,j],'fro')/norm_data)**2;
                        tmp_gap += (norm(tmp3[:,:,k,j]-tmp2[:,:,k,j],'fro')/norm_data)**2;
            
            change[iters] = sqrt(tmp_change);
            gap[iters] = sqrt(tmp_gap);
            
            u = tmp_u;
            
        
        ##### POSTPROCESSING
        tmp = proj1.work(u);
        tmp2 = proj2.work(u);
        
        if self.Ny == 1:
            u1 = tmp[:,1];
            u2 = tmp2[:,1];
        elif self.Nx == 1:
            u1 = tmp[1,:];
            u2 = tmp2[1,:];
        elif self.Nz == 1:
            u1 = tmp[:,:,1];
            u2 = tmp2[:,:,1];
        else:
            u1 = tmp;
            u2 = tmp2;
        change = change[1:iters+1];
        gap = gap[1:iters+1];    
        
        return u1, u2, iters, change, gap;
