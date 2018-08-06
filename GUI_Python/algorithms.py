# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 10:40:18 2014

@author: stefan
"""
import numpy, scipy.linalg, math

__all__ = ["Algorithm","RAAR","HPR","AP"]


class Algorithm:
    """
    Generic interface for ProxToolbox algorithms.
    """
    
    def __init__(self):
        raise NotImplementedError,"This is just an abstract interface"
        
    def run(self, u, tol, maxiter):
        """
        Runs the algorithm one some input data
        
        Parameters
        ----------
        u : array_like
            Input data
        tol : number
            Tolerance
        maxiter : int
            Maximum number of iterations
        
        Returns
        -------
        u1 : array_like
            Result
        u2 : array_like
            Result
        iters : int
            Number of iterations the algorithm performed
        change : array_like
            The percentage change in the norm
        gap : array_like
            Squared gap distance normalized by the magnitude constraint
        """
        raise NotImplementedError,"This is just an abstract interface"





class RAAR(Algorithm):
    """
    Relaxed Averaged Alternating Reflection algorithm
    
    See Also
    --------
        Algorithm : Generic interface for ProxToolbox algorithms
    """
    
    def __init__(self,config):
        """
        Parameters
        ----------
        config : dict
            Dictionary containing the problem configuration. It must contain the
            following mappings:
                proj1:ProxOperator
                    First ProxOperator (the class, no instance)
                proj2:ProxOperator
                    Second ProxOperator (the class, no instance)
                beta0:number
                    Starting relaxation parmater
                beta_max:number
                    Maximum relaxation parameter
                beta_switch:int
                    Iteration at which beta moves from beta0 -> beta_max
                norm_data:number
                    ?
                Nx:int
                    Row-dim of the product space elements
                Ny:int
                    Column-dim of the product space elements
                Nz:int
                    Depth-dim of the product space elements
                dim:int
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

    def run(self, u, tol, maxiter):
        
        beta0 = self.beta0; beta_max = self.beta_max; beta_switch = self.beta_switch;
        proj1 = self.proj1; proj2 = self.proj2;
        norm_data = self.norm_data;
        Nx = self.Nx; Ny = self.Ny; Nz = self.Nz;  dim = self.dim;
        
        beta = beta0;
        iters = 0;
        change = numpy.zeros(maxiter+1,dtype=u.dtype);
        change[0] = 999;
        gap = change.copy();
        
        tmp1 = 2*proj2.work(u) - u;
        
        while iters < maxiter and change[iters] >= tol:
            tmp = math.exp((-iters/beta_switch)**3);
            beta = (tmp*beta0) + ((1-tmp)*beta_max);
            iters += 1;
            
            tmp3 = proj1.work(tmp1);
            tmp_u = ((beta*(2*tmp3-tmp1)) + ((1-beta)*tmp1) + u)/2;
            tmp2 = proj2.work(tmp_u);
            
            tmp3 = proj1.work(tmp2);
            
            tmp_change = 0; tmp_gap = 0;
            if Ny == 1 or Nx == 1:
                tmp_change = (scipy.linalg.norm(u-tmp_u,'fro')/norm_data)**2;
                tmp_gap = (scipy.linalg.norm(tmp3-tmp2,'fro')/norm_data)**2;
            elif Nz == 1:
                for j in range(dim):
                    tmp_change += (scipy.linalg.norm(u[:,:,j]-tmp_u[:,:,j],'fro')/norm_data)**2;
                    tmp_gap += (scipy.linalg.norm(tmp3[:,:,j]-tmp2[:,:,j])/norm_data,'fro')**2;
            else:
                for j in range(dim):
                    for k in range(Nz):
                        tmp_change += (scipy.linalg.norm(u[:,:,k,j]-tmp_u[:,:,k,j],'fro')/norm_data)**2;
                        tmp_gap += (scipy.linalg.norm(tmp3[:,:,k,j]-tmp2[:,:,k,j],'fro')/norm_data)**2;
            
            change[iters] = math.sqrt(tmp_change);
            gap[iters] = math.sqrt(tmp_gap);
            
            u = tmp_u;
            tmp1 = (2*tmp2) - tmp_u;
            
        u = tmp2;
        tmp = proj1.work(u);
        tmp2 = proj2.work(u);
        if Ny == 1:
            u1 = tmp[:,1];
            u2 = tmp2[:,1];
        elif Nx == 1:
            u1 = tmp[1,:];
            u2 = tmp2[1,:];
        elif Nz == 1:
            u1 = tmp[:,:,1];
            u2 = tmp2[:,:,1];
        else:
            u1 = tmp;
            u2 = tmp2;
        change = change[1:iters+1];
        gap = gap[1:iters+1];    
        
        return u1, u2, iters, change, gap;




class HPR(Algorithm):
    """
    Heinz-Patrick-Russell algorithm
    
    See Also
    --------
        Algorithm : Generic interface for ProxToolbox algorithms
    
    References
    ----------
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
                proj1:ProxOperator
                    First ProxOperator (the class, no instance)
                proj2:ProxOperator
                    Second ProxOperator (the class, no instance)
                beta0:number
                    Starting relaxation parmater
                beta_max:number
                    Maximum relaxation parameter
                beta_switch:int
                    Iteration at which beta moves from beta0 -> beta_max
                norm_data:number
                    ?
                Nx:int
                    Row-dim of the product space elements
                Ny:int
                    Column-dim of the product space elements
                Nz:int
                    Depth-dim of the product space elements
                dim:int
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
    
    
    def run(self, u, tol, maxiter):
        beta0 = self.beta0; beta_max = self.beta_max; beta_switch = self.beta_switch;
        proj1 = self.proj1; proj2 = self.proj2;
        norm_data = self.norm_data;
        Nx = self.Nx; Ny = self.Ny; Nz = self.Nz;  dim = self.dim;
        
        beta = beta0;
        iters = 0;
        change = numpy.zeros(maxiter+1,dtype=u.dtype);
        change[0] = 999;
        gap = change.copy();
        
        tmp = proj2.work(u);
        
        while iters < maxiter and change[iters] >= tol:
            a = math.exp((-iters/beta_switch)**3);
            beta = (a*beta0) + ((1-a)*beta_max);
            iters += 1;
            
            tmp2 = 2*proj2.work(u) - u + (beta-1)*tmp;
            tmp_u = (2*proj1.work(tmp2) - tmp2 + u + (1-beta)*tmp)/2;
            tmp = proj2.work(tmp_u);
            tmp3 = proj1.work(tmp);
            
            tmp_change = 0; tmp_gap = 0;
            if Ny == 1 or Nx == 1:
                tmp_change = (scipy.linalg.norm(u-tmp_u,'fro')/norm_data)**2;
                tmp_gap = (scipy.linalg.norm(tmp3-tmp2,'fro')/norm_data)**2;
            elif Nz == 1:
                for j in range(dim):
                    tmp_change += (scipy.linalg.norm(u[:,:,j]-tmp_u[:,:,j],'fro')/norm_data)**2;
                    tmp_gap += (scipy.linalg.norm(tmp3[:,:,j]-tmp2[:,:,j])/norm_data,'fro')**2;
            else:
                for j in range(dim):
                    for k in range(Nz):
                        tmp_change += (scipy.linalg.norm(u[:,:,k,j]-tmp_u[:,:,k,j],'fro')/norm_data)**2;
                        tmp_gap += (scipy.linalg.norm(tmp3[:,:,k,j]-tmp2[:,:,k,j],'fro')/norm_data)**2;
            
            change[iters] = math.sqrt(tmp_change);
            gap[iters] = math.sqrt(tmp_gap);
            
            u = tmp_u;
            
        tmp = proj1.work(u);
        tmp2 = proj2.work(u);
        
        if Ny == 1:
            u1 = tmp[:,1];
            u2 = tmp2[:,1];
        elif Nx == 1:
            u1 = tmp[1,:];
            u2 = tmp2[1,:];
        elif Nz == 1:
            u1 = tmp[:,:,1];
            u2 = tmp2[:,:,1];
        else:
            u1 = tmp;
            u2 = tmp2;
        change = change[1:iters+1];
        gap = gap[1:iters+1];    
        
        return u1, u2, iters, change, gap;


class AP(Algorithm):
    """
    Alternating Projections
    
    See Also
    --------
        Algorithm : Generic interface for ProxToolbox algorithms
    """
    
    def __init__(self, config):
        """
        config : dict
            Dictionary containing the problem configuration. It must contain the
            following mappings:
                proj1:ProxOperator
                    First ProxOperator (the class, no instance)
                proj2:ProxOperator
                    Second ProxOperator (the class, no instance)
                beta0:number
                    Starting relaxation parmater
                beta_max:number
                    Maximum relaxation parameter
                beta_switch:int
                    Iteration at which beta moves from beta0 -> beta_max
                norm_data:number
                    ?
                Nx:int
                    Row-dim of the product space elements
                Ny:int
                    Column-dim of the product space elements
                Nz:int
                    Depth-dim of the product space elements
                dim:int
                    Size of the product space
        """
        self.proj1 = config['proj1'](config); self.proj2 = config['proj2'](config);
        self.norm_data = config['norm_data'];
        self.Nx = config['Nx']; self.Ny = config['Ny']; self.Nz = config['Nz'];
        self.dim = config['dim'];
    
    
    def run(self, u, tol, maxiter):
        proj1 = self.proj1; proj2 = self.proj2;
        norm_data = self.norm_data;
        Nx = self.Nx; Ny = self.Ny; Nz = self.Nz;  dim = self.dim;
        
        iters = 0;
        change = numpy.zeros(maxiter+1);
        change[0] = 999;
        gap = change.copy();
        
        tmp1 = proj2.work(u);
        
        while iters < maxiter and change[iters] >= tol:
            iters += 1;
            
            tmp_u = proj1.work(tmp1);
            tmp1 = proj2.work(tmp_u);
            
            tmp_change = 0; tmp_gap = 0;
            if Ny == 1 or Nx == 1:
                tmp_change = (scipy.linalg.norm(u-tmp_u,'fro')/norm_data)**2;
                tmp_gap = (scipy.linalg.norm(tmp1-tmp_u,'fro')/norm_data)**2;
            elif Nz == 1:
                for j in range(dim):
                    tmp_change += (scipy.linalg.norm(u[:,:,j]-tmp_u[:,:,j],'fro')/norm_data)**2;
                    tmp_gap += (scipy.linalg.norm(tmp1[:,:,j]-tmp_u[:,:,j])/norm_data,'fro')**2;
            else:
                for j in range(dim):
                    for k in range(Nz):
                        tmp_change += (scipy.linalg.norm(u[:,:,k,j]-tmp_u[:,:,k,j],'fro')/norm_data)**2;
                        tmp_gap += (scipy.linalg.norm(tmp1[:,:,k,j]-tmp_u[:,:,k,j],'fro')/norm_data)**2;
            change[iters] = math.sqrt(tmp_change);
            gap[iters] = math.sqrt(tmp_gap);
            
            u = tmp_u.copy();
            
        tmp = proj1.work(u);
        tmp2 = proj2.work(u);
        
        if Ny == 1:
            u1 = tmp[:,1];
            u2 = tmp2[:,1];
        elif Nx == 1:
            u1 = tmp[1,:];
            u2 = tmp2[1,:];
        elif Nz == 1:
            u1 = tmp[:,:,1];
            u2 = tmp2[:,:,1];
        else:
            u1 = tmp;
            u2 = tmp2;
        change = change[1:iters+1];
        gap = gap[1:iters+1];    
        
        return u1, u2, iters, change, gap;


