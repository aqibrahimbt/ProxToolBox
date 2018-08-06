# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 12:48:26 2015

@author: rebecca
"""

from numpy import zeros, angle, trace, exp, sqrt

from numpy.linalg import norm

from .algorithms import Algorithm

class RAAR_expert(Algorithm):
    """
    Relaxed Averaged Alternating Reflection algorithm
    """
    
    def __init__(self,config):
        """
        Parameters
        ----------
        config : dict        
                 Dictionary containing the problem configuration.
                 It must contain the following mappings:
            
                proxoperators: 2 ProxOperators
                    Tuple of ProxOperators (the class, no instance)
                beta_0: number
                    Starting relaxation parmater
                beta_max: number
                    Maximum relaxation parameter
                beta_switch: int
                    Iteration at which beta moves from beta_0 -> beta_max
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

        self.prox1 = config['proxoperators'][0](config)
        self.prox2 = config['proxoperators'][1](config)
        self.norm_data = config['norm_data']
        self.beta_0 = config['beta_0']
        self.beta_max = config['beta_max']
        self.beta_switch = config['beta_switch']
        self.Nx = config['Nx']; self.Ny = config['Ny']; self.Nz = config['Nz'];
        self.product_space_dimension = config['product_space_dimension']
        self.iter = 0

        if 'truth' in config:
            self.truth = config['truth']
            self.truth_dim = config['truth_dim']
            self.norm_truth = config['norm_truth']

    def run(self, u, tol, maxiter):
        """
        Runs the algorithm for the specified input data
        """
        
        ##### PREPROCESSING
        norm_data = self.norm_data
        
        beta = self.beta_0
        iter = self.iter
        if u.ndim < 3:
            p = 1
            q = 1
        elif u.ndim == 3:
            p = u.shape[2]
            q = 1
        else:
            p = u.shape[2]
            q = u.shape[3]

        change = zeros(maxiter+1,dtype=u.dtype)
        change[0] = 999
        gap = change.copy()
        shadow_change = change.copy()
        if  hasattr(self, 'truth'):
            Relerrs = change.copy()

        tmp1 = 2*self.prox2.work(u) - u
        shadow = self.prox1.work(u)

        
        
        ##### LOOP
        while iter < maxiter and change[iter] >= tol:

            tmp = exp((-(iter+1)/self.beta_switch)**3);
            beta = (tmp*self.beta_0) + ((1-tmp)*self.beta_max)
            
            iter += 1;
            
            tmp3 = self.prox1.work(tmp1)
            tmp_u = ((beta*(2*tmp3-tmp1)) + ((1-beta)*tmp1) + u)/2
            tmp2 = self.prox2.work(tmp_u)
            
            tmp3 = self.prox1.work(tmp2)
            
            tmp_change = 0; tmp_gap = 0; tmp_shadow=0;

            if p==1 and q==1:
              tmp_change= (norm(u-tmp_u, 'fro')/norm_data)**2
              tmp_shadow = (norm(tmp3-shadow,'fro')/norm_data)**2
              tmp_gap = (norm(tmp3-tmp2,'fro')/norm_data)**2

              if hasattr(self, 'truth'):
                  if self.truth_dim[0] == 1:
                      z=tmp3[0,:]
                  elif self.truth_dim[1] == 1:
                      z=tmp3[:,0]
                  else:
                      z=tmp3;
                  Relerrs[iter] = norm(self.truth - exp(-1j*angle(trace(self.truth.T*z))) * z, 'fro')/self.norm_truth

            elif q==1:
              for j in range(self.product_space_dimension):
                  tmp_change= tmp_change+ (norm(u[:,:,j]-tmp_u[:,:,j], 'fro')/norm_data)**2
                  # compute (||P_SP_Mx-P_Mx||/norm_data)^2:
                  tmp_gap = tmp_gap+(norm(tmp3[:,:,j]-tmp2[:,:,j],'fro')/norm_data)**2
                  tmp_shadow = tmp_shadow+(norm(tmp3[:,:,j]-shadow[:,:,j],'fro')/norm_data)**2

              if hasattr(self, 'truth'):
                 z=tmp3[:,:,0]
                 Relerrs[iter] = norm(self.truth - exp(-1j*angle(trace(self.truth.T*z))) * z, 'fro')/self.norm_truth

            else:
              for j in range(self.product_space_dimension):
                  for k in range(Nz):
                      tmp_change= tmp_change+ (norm(u[:,:,k,j]-tmp_u[:,:,k,j], 'fro')/norm_data)**2
                      # compute (||P_Sx-P_Mx||/norm_data)^2:
                      tmp_gap = tmp_gap+(norm(tmp3[:,:,k,j]-tmp2[:,:,k,j],'fro')/(norm_data))**2
                      tmp_shadow = tmp_shadow+(norm(tmp3[:,:,k,j]-shadow[:,:,k,j],'fro')/(norm_data))**2
   
            change[iter] = sqrt(tmp_change)
            gap[iter] = sqrt(tmp_gap)
            shadow_change[iter] = sqrt(tmp_shadow) # this is the Euclidean norm of the gap to
            # the unregularized set.  To monitor the Euclidean norm of the gap to the
            # regularized set is expensive to calculate, so we use this surrogate.
            # Since the stopping criteria is on the change in the iterates, this
            # does not matter.
            # graphics
            u = tmp_u
            tmp1 = (2*tmp2) - tmp_u
            
        
        ##### POSTPROCESSING
        u = tmp2;
        tmp = self.prox1.work(u);
        uu = tmp
        tmp2 = self.prox2.work(u);
        if self.Nx == 1:
            u1 = tmp[:,0];
            u2 = tmp2[:,0];
        elif self.Ny == 1:
            u1 = tmp[0,:];
            u2 = tmp2[0,:];
        elif self.Nz == 1 and tmp.ndim > 2:
            u1 = tmp[:,:,0]
            u2 = tmp2[:,:,0]
        else:
            u1 = tmp;
            u2 = tmp2;

        change = change[1:iter+1];
        gap = gap[1:iter+1];
        shadow_change = shadow_change[1:iter+1]

        output = {'u' : uu, 'u1': u1, 'u2': u2, 'iter': iter, 'change': change, 'gap': gap, 'shadow_change' : shadow_change}
        if hasattr(self, 'truth'):
            output['Relerrs']=Relerrs
        return output
