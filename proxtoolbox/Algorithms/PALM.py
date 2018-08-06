# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 18:06:37 2015

@author: rebecca
"""

from numpy import zeros
from .algorithms import Algorithm

class PALM(Algorithm):
    """
    Proximal alternating (linearized) minimization algorithm
    
    .. seealso:: :class:`.Algorithm`
    """
    
    def __init__(self, config):
        """
        :param: config (dictionary)
        
        Dictionary containing the problem configuration. It must contain the
        following mappings:
        
            projector_sequence:sequence of Proxoperator
                 Sequence of prox operators to be used. (The classes, no instances)
            Nx:int
                Row-dim of the product space elements
            Ny:int
                Column-dim of the product space elements
            Nz:int
                Depth-dim of the product space elements
            dim:int
                Size of the product space
            norm_data:number
                ?
            ignore_error:boolean
                Whether to ignore errors
        """
        self.projs = [p(config) for p in config['projector_sequence']];
        
        self.ignore_error = config['ignore_error'];
        self.custom_error = PtychographyStats(config);
    
    
    def run(self, u, tol, maxiter):
        """
        Runs the algorithm for the specified input data
        """
        projs = self.projs;
        ignore_error = self.ignore_error;
        
        iters = 0;
        change = zeros(maxiter+1);
        change[0] = tol+1;        

        if ignore_error == False:
            num_custom_errors = len(self.custom_error.customerror(u,u))
            custom_errors = zeros((num_custom_errors,maxiter+1));
        
        tmp_u = u;
        
        while iters < maxiter and change[iters] > tol:
            iters += 1;
            
            for p in projs:
                tmp_u = p.work(tmp_u);
            
            change[iters] = self.custom_error.change(u,tmp_u);
            
            if ignore_error == False:
                custom_errors[:,iters] = self.custom_error.customerror(u,tmp_u);
                
            u = tmp_u;
        
        change = change[1:iters+1];
        custom_errors = numpy.empty((5,0)) if ignore_error else custom_errors[:,1:iters+1];
        
        return tmp_u, tmp_u.copy(), iters, change, custom_errors;
