#                      GRAAL.m
#             written on Aug 18, 2017 by
#                   Russell Luke
#   Inst. Fuer Numerische und Angewandte Mathematik
#                Universitaet Goettingen
#
# DESCRIPTION:  user-friendly version of Malitskyi's Golden-Ratio algorithm 
#
#
# INPUT:  u, an array, and self, a data structure
#
# OUTPUT: self, an appended data structure with 
#               xnew = the algorithm fixed point
#
# USAGE: [xnew, self] = GRAAL_simple(x, self)
#
# Nonstandard Matlab function calls:  self.Prox1 and .Prox2


from .SimpleAlgortihm import SimpleAlgorithm
from numpy.linalg import norm


class GRAAL(SimpleAlgorithm):

#    """
#   Golden Ratio Algorithm for the problem x = Tx
#
#    T is the operator
#    J is the energy which we want to check, like |x-x*|
#   x0 is the starting point
#
#    """


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
        self.Nx = config['Nx']
        self.Ny = config['Ny']
        self.Nz = config['Nz'];
        self.product_space_dimension = config['product_space_dimension']
        self.iter = 0

        if 'truth' in config:
            self.truth = config['truth']
            self.truth_dim = config['truth_dim']
            self.norm_truth = config['norm_truth']

        self.diagnostic = 'diagnostic' in config
        self.GRAAL_initialized = False

        self.phi = 1.4 if 'phi' not in config else config['phi']
        self.tau = 1. / self.phi + 1. / self.phi**2

        if 'inv_Lipschitz_const' not in config:
            self.inv_Lipschitz_const = 0.5
        else:
            self.inv_Lipschitz_const = config['inv_Lipschitz_const']

        self.th = 1

        self.config = config


    def evaluate(self, x):

        if(not(self.GRAAL_initialized)):
            self.xbar = x
            tmp1 = self.prox1.work(x)
            tmp2 = self.prox2.work(tmp1)
            self.Fx = x-tmp2 
            self.GRAAL_initialized= True
            x1=x
            return x1+999
        else:
            phi= self.phi
            la= self.inv_Lipschitz_const
            x1 = self.xbar - la * self.Fx
            tmp1 = self.prox1.work(x1)
            tmp2 = self.prox2.work(tmp1)
            Fx1 = x1-tmp2
            
            n1 = norm(x1 - x,'fro')**2
            n2 = norm(Fx1 - self.Fx, 'fro')**2

            la1 = min(self.tau * la, 0.25 * phi * self.th / la * (n1 / n2))

            self.xbar = ((phi - 1) * x1 + self.xbar) / phi
            self.th = phi * la1 / la
            self.inv_Lipschitz_const=la1
            self.Fx=Fx1
            
            return x1
