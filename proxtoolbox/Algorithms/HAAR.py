#                      HAAR.m
#             written on May 31, 2006 by
#          Russell Luke & Daniel Gempesaw
#               University of Delaware
#
# DESCRIPTION:  Heaugazeau-variant of RAAR, as proposed in 
#   H. H. Bauschke, P. L. Combettes and D. R. Luke,
#   `` A Strongly Convergent Reflection Method for Finding the 
#    Projection onto the Intersection of Two Closed Convex 
#    Sets in a Hilbert Space"  
#    Journal of Approximation Theory 141(1):63-69 (2006).
# 
#  The fixed point of the algorithm is in fact the PROJECTION
#  onto the intersection.  
#
# INPUT:  method_input, a data structure
#
# OUTPUT: method_output, a data structure with 
#               u = the algorithm fixed point
#               stats = [iter, change, gap]  where
#                     gap    = squared gap distance normalized by
#                              the magnitude constraint.
#                     change = the percentage change in the norm
#                              squared difference in the iterates
#                     iter   = the number of iterations the algorithm
#                              performed
#
# USAGE: method_output = HAAR(method_input)
#
# Nonstandard Matlab function calls:  method_input.Prox1 and .Prox2, RAAR, Q_Heau

from numpy import zeros, shape, sqrt
from scipy.linalg import norm
from .algorithms import Algorithm
from proxtoolbox import Algorithms
from proxtoolbox import ProxOperators


class HAAR(Algorithm):

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
        self.prox1 = config['proxoperators'][0](config)
        self.prox2 = config['proxoperators'][1](config)
        self.norm_data = config['norm_data']
        self.Nx = config['Nx']
        self.Ny = config['Ny']
        self.Nz = config['Nz']
        self.product_space_dimension = config['product_space_dimension']
        self.iter = 0

        if 'truth' in config:
            self.truth = config['truth']
            self.truth_dim = config['truth_dim']
            self.norm_truth = config['norm_truth']

        self.T='RAAR_expert'
        self.T_config = config.copy()
        self.T_config['MAXIT']=1
        self.T_config['beta_max']=.1
        self.T_config['beta_0']=.1
        self.T_config['anim']=0
        self.T = getattr(Algorithms, self.T)(self.T_config)

        self.diagnostic = 'diagnostic' in config
        self.y0 = config['u_0']
        self.Q_Heau = getattr(ProxOperators, "Q_Heau")
        self.mu = config['beta_0']

    def run(self, u, tol, maxiter):
        """
        Runs the algorithm for the specified input data
        """
        prox1 = self.prox1
        prox2 = self.prox2;
        norm_data = self.norm_data

        y0 = self.y0

        # preallocate the error monitors:

        #[m,n,p,q]
        s =shape(y0)
        dim_number = y0.ndim

        change = zeros(maxiter+1,dtype=u.dtype)
        change[0] = 999

        ######################################
        #  set up diagnostic arrays
        ######################################
        if self.diagnostic:
            gap = change.copy()
            shadow_change = change.copy()
        if  hasattr(self, 'truth'):
            Relerrs = change.copy()

        y=y0

        if dim_number==2 and (s[0] !=1 ) and(s[1] !=1):
            y0=y.reshape((s[0]*s[1],1))
        elif dim_number == 3:
            y0=zeros((s[0]*s[1]*s[2],1,1))
            for j in range(s[3]):
                ytmp=y[:,:,j].reshape((s[0]*s[1],1))
                y0[n*m*j:m*n*(j+1)]= ytmp
        elif dim_numner == 4:
            print('not ready for 4D arrays')
            return

        shadow = prox2.work(y)
        iter = 0

        while iter < maxiter and change[iter] >= tol:

            iter += 1
            self.iter=iter

            # next iterate
            # Ty= RAARy with beta=1 and MAXIT=1
            T_output = self.T.run(self.T_config['u_0'],self.T_config['TOL'],self.T_config['MAXIT'])
            Ty = T_output['u']

            if dim_number == 2:
                Ty = Ty.reshape((s[0]*s[1],1))
                tmp_y= y.reshape((s[0]*s[1],1))
            elif dim_number == 3:
                tmp_y = zeros((s[0]*s[1]*s[2],1,1))
                tmp_Ty=tmp_y
                for j in range(p):
                    ytmp = y[:,:,j].reshape((s[0]*s[1],1))
                    tmp_y[s[0]*s[1]*j:s[0]*s[1]*(j+1)] = ytmp
                    ytmp = Ty[:,:,j].reshape((s[0]*s[1],1))
                    tmp_Ty[s[0]*s[1]*j:s[0]*s[1]*(j+1)] = ytmp
                Ty=tmp_Ty
            elif dim_number == 4:
                print('not ready for 4D arrays')

            y_new = self.Q_Heau.work(y0,tmp_y,(1-self.mu)*tmp_y+self.mu*Ty)

            if dim_number == 2 and (s[0]!=1) and (s[1]!=1):
                y_new = y_new.reshape((s[0],s[1]))
            elif dim_number == 3:
                tmpy = zeros((s[0],s[1],s[2]))
                for j in range(p):
                    tmpy[:,:,j]= y_new[s[0]*s[1]*j:s[0]*s[1]*(j+1)].reshape((s[0],s[1]))
                y_new=tmpy
            elif dim_number == 4:
                print('not ready for 4D arrays')


            if self.diagnostic == True:
                # the next prox operations only gets used in the computation of
                # the size of the gap.  This extra step is not
                # required in alternating projections, which makes RAAR
                # and special cases thereof more expensive to monitor.
                # compute the normalized change in successive iterates:
                tmp = prox2.work(y_new)
                tmp2 = prox1.work(tmp)

            # compute the normalized change in successive iterates: 
            # change(iter) = sum(sum((feval('P_M',M,u)-tmp).^2))/norm_data;
            tmp_change=0
            tmp_gap=0
            tmp_shadow=0;

            if dim_number <= 2:
                tmp_change= (norm(y-y_new, 'fro')/norm_data)**2
                if self.diagnostic == True:
                    # For Douglas-Rachford,in general it is appropriate to monitor the
                    # SHADOWS of the iterates, since in the convex case these converge
                    # even for beta=1.
                    # (see Bauschke-Combettes-Luke, J. Approx. Theory, 2004)
                    tmp_shadow = (norm(tmp-shadow,'fro')/norm_data)**2
                    tmp_gap = (norm(tmp-tmp2,'fro')/norm_data)**2

                    if hasattr(self, 'truth'):
                      if self.truth_dim[0] == 1:
                          z=tmp[0,:]
                      elif self.truth_dim[1] == 1:
                          z=tmp[:,0]
                      else:
                          z=tmp;
                      Relerrs[iter] = norm(self.truth - exp(-1j*angle(trace(self.truth.T*z))) * z, 'fro')/self.norm_truth

            elif dim_number == 3:
                for j in range(self.product_space_dimension):
                    tmp_change= tmp_change+ (norm(y[:,:,j]-y_new[:,:,j], 'fro')/norm_data)**2
                    if self.diagnostic == True:
                        # compute (||P_SP_Mx-P_Mx||/norm_data)^2:
                        tmp_gap = tmp_gap+(norm(tmp[:,:,j]-tmp2[:,:,j],'fro')/norm_data)**2
                        tmp_shadow = tmp_shadow+(norm(tmp[:,:,j]-shadow[:,:,j],'fro')/norm_data)**2
                if self.diagnostic == True:
                    if hasattr(self, 'truth'):
                        z=tmp[:,:,0]
                        Relerrs[iter] = norm(self.truth - exp(-1j*angle(trace(self.truth.T*z))) * z, 'fro')/self.norm_truth


            change[iter] = sqrt(tmp_change)
            if self.diagnostic == True:
                gap[iter] = sqrt(tmp_gap)
                shadow_change[iter] = sqrt(tmp_shadow) # this is the Euclidean norm of the gap to
                # the unregularized set.  To monitor the Euclidean norm of the gap to the
                # regularized set is expensive to calculate, so we use this surrogate.
                # Since the stopping criteria is on the change in the iterates, this
                # does not matter.
                # graphics
            # update
            y=y_new
            self.T_config['u_0']=y

            if self.diagnostic == True:
                # For Douglas-Rachford,in general it is appropriate to monitor the
                # SHADOWS of the iterates, since in the convex case these converge
                # even for beta=1.
                # (see Bauschke-Combettes-Luke, J. Approx. Theory, 2004)
                shadow=tmp


        ##### POSTPROCESSING
        u = tmp2
        tmp1 = self.prox1.work(u)
        uu =tmp1
        tmp2 = self.prox2.work(u)

        if self.Nx == 1:
            u1 = tmp[:,0];
            u2 = tmp2[:,0];
        elif self.Ny == 1:
            u1 = tmp[0,:];
            u2 = tmp2[0,:];
        elif self.Nz == 1 and u1.ndim > 2:
            u1 = tmp[:,:,0]
            u2 = tmp2[:,:,0]

        change = change[1:iter+1]

        output = {'u' : u, 'u1': u1, 'u2': u2, 'iter': iter, 'change': change}

        if self.diagnostic == True:
            gap = gap[1:iter+1]
            shadow_change = shadow_change[1:iter+1]
            output['gap']=gap
            output['shadow_change']= shadow_change
            if hasattr(self, 'truth'):
                output['Relerrs']=Relerrs
        return output
