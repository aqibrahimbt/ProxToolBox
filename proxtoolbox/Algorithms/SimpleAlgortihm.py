#Simple Algortihm
#Analog to AlgortihmWrapper in ProxMatlab

from .algorithms import Algorithm

from numpy import zeros, angle, trace, exp, sqrt

from numpy.linalg import norm

class SimpleAlgorithm:

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
        self.Nx = config['Nx']; self.Ny = config['Ny']; self.Nz = config['Nz']
        self.product_space_dimension = config['product_space_dimension']
        self.iter = 0
        self.config = config

        if 'truth' in config:
            self.truth = config['truth']
            self.truth_dim = config['truth_dim']
            self.norm_truth = config['norm_truth']

        if 'diagnostic' in config:
            self.diagnostic = True
        else:
            self.diagnostic = False


    def run(self, u, tol, maxiter):
        """
        Runs the algorithm for the specified input data
        """

        ##### PREPROCESSING
        config = self.config
        norm_data = self.norm_data

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

        ######################################
        #  set up diagnostic arrays
        ######################################
        if self.diagnostic:
            gap = change.copy()
            shadow_change = change.copy()
        if  hasattr(self, 'truth'):
            Relerrs = change.copy()



        # Different algorithms will have different quatities that one 
        # should monitor.  Since we take a fixed point perspective, the main
        # thing to monitor is the change in the iterates.
        shadow = self.prox1.work(u)

        ##### LOOP
        while iter < maxiter and change[iter] >= tol:
            iter += 1
            config['iter'] =iter
            #makes call to algorithm:
            config['u']=u
            u_new = self.evaluate(u)

            if 'diagnostic' in self.config:
                # the next prox operation only gets used in the computation of
                # the size of the gap.  This extra step is not
                # required in alternating projections, which makes RAAR
                u2 = self.prox2.work(u_new)
                u1 = self.prox1.work(u2)

            # compute the normalized change in successive iterates: 
            # change(iter) = sum(sum((feval('P_M',M,u)-tmp).^2))/norm_data;
            tmp_change=0; tmp_gap=0; tmp_shadow=0;

            if p==1 and q==1:
                tmp_change= (norm(u-u_new, 'fro')/norm_data)**2
                if 'diagnostic' in self.config:
                    # For Douglas-Rachford,in general it is appropriate to monitor the
                    # SHADOWS of the iterates, since in the convex case these converge
                    # even for beta=1.
                    # (see Bauschke-Combettes-Luke, J. Approx. Theory, 2004)
                    tmp_shadow = (norm(u2-shadow,'fro')/norm_data)**2
                    tmp_gap = (norm(u1-u2,'fro')/norm_data)**2

                    if hasattr(self, 'truth'):
                      if self.truth_dim[0] == 1:
                          z=u1[0,:]
                      elif self.truth_dim[1] == 1:
                          z=u1[:,0]
                      else:
                          z=u1;
                      Relerrs[iter] = norm(self.truth - exp(-1j*angle(trace(self.truth.T*z))) * z, 'fro')/self.norm_truth

            elif q==1:
                for j in range(self.product_space_dimension):
                    tmp_change= tmp_change+ (norm(u[:,:,j]-u_new[:,:,j], 'fro')/norm_data)**2
                    if 'diagnostic' in self.config:
                        # compute (||P_SP_Mx-P_Mx||/norm_data)^2:
                        tmp_gap = tmp_gap+(norm(u1[:,:,j]-u2[:,:,j],'fro')/norm_data)**2
                        tmp_shadow = tmp_shadow+(norm(u2[:,:,j]-shadow[:,:,j],'fro')/norm_data)**2
                if 'diagnostic' in self.config:
                    if hasattr(self, 'truth'):
                     z=u1[:,:,0]
                     Relerrs[iter] = norm(self.truth - exp(-1j*angle(trace(self.truth.T*z))) * z, 'fro')/self.norm_truth

            else:
                if 'diagnostic' in self.config:
                    if hasattr(self, 'truth'):
                        Relerrs[iter] = 0
                for j in range(self.product_space_dimension):
                  for k in range(self.Nz):
                    tmp_change= tmp_change+ (norm(u[:,:,k,j]-u_new[:,:,k,j], 'fro')/norm_data)**2
                    if 'diagnostic' in self.config:
                        # compute (||P_Sx-P_Mx||/norm_data)^2:
                        tmp_gap = tmp_gap+(norm(u1[:,:,k,j]-u2[:,:,k,j],'fro')/(norm_data))**2
                        tmp_shadow = tmp_shadow+(norm(u2[:,:,k,j]-shadow[:,:,k,j],'fro')/(norm_data))**2
                if hasattr(self, 'truth') and (j==0):
                    Relerrs[iter] = Relerrs[iter] + norm(self.truth - exp(-1j*angle(trace(self.truth.T*u1[:,:,k,1]))) * u1[:,:,k,1], 'fro')/self.norm_truth


            change[iter] = sqrt(tmp_change)
            if 'diagnostic' in self.config:
                gap[iter] = sqrt(tmp_gap)
                shadow_change[iter] = sqrt(tmp_shadow) # this is the Euclidean norm of the gap to
                # the unregularized set.  To monitor the Euclidean norm of the gap to the
                # regularized set is expensive to calculate, so we use this surrogate.
                # Since the stopping criteria is on the change in the iterates, this
                # does not matter.
                # graphics

            #update
            u=u_new
            if 'diagnostic' in self.config:
                # For Douglas-Rachford,in general it is appropriate to monitor the
                # SHADOWS of the iterates, since in the convex case these converge
                # even for beta=1.
                # (see Bauschke-Combettes-Luke, J. Approx. Theory, 2004)
                shadow=u2

        ##### POSTPROCESSING
        u2 = self.prox2.work(u);
        u1 = self.prox1.work(u2);
        if self.Nx == 1:
            u1 = u1[:,0];
            u2 = u2[:,0];
        elif self.Ny == 1:
            u1 = u1[0,:];
            u2 = u2[0,:];
        elif self.Nz == 1 and u1.ndim > 2:
            u1 = u1[:,:,0]
            u2 = u2[:,:,0]

        change = change[1:iter+1];

        output = {'u' : u, 'u1': u1, 'u2': u2, 'iter': iter, 'change': change}
        if 'diagnostic' in self.config:
            gap = gap[1:iter+1]
            shadow_change = shadow_change[1:iter+1]
            output['gap']=gap
            output['shadow_change']= shadow_change
            if hasattr(self, 'truth'):
                output['Relerrs']=Relerrs
        return output
