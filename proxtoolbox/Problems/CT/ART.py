#                      ART.m
#             written on May 23, 2012 by 
#                   Russell Luke
#   Inst. Fuer Numerische und Angewandte Mathematik
#                Universitaet Gottingen
#            Modified by:  Matthew Tam
#             University of Newcastle
#                 30th Oct 2013.
#
#
# This module sets the projectors based on the problem type
# within the problem family 'ART'
# This module sets the projectors based on the problem type
# For the problem family 'Custom' it doesn't really do anything

from proxtoolbox.Problems.problems import Problem
from proxtoolbox import ProxOperators
from proxtoolbox.ProxOperators.proxoperators import ProxOperator
from proxtoolbox.Problems.CT import Graphics
from proxtoolbox import Algorithms
from numpy.linalg import norm
from numpy import sqrt

class ART(Problem):
    """
    ART class
    """

    config = {
    }

    # for problem family 'Custom' we assume that the user has 
    # named the projectors in either the parameter_filename 
    # or the data_filename files (and, of course,  written 
    # the projectors)
    # apparently `exist' doesn't allow you to check elements of 
    # data structures...:(
    # if(~exist('input.product_space_dimension','var'))
    #    input.product_space_dimension = 1;
    # end
    # if you are only working with two sets but 
    # want to do averaged projections
    # (= alternating projections on the product space)
    # or RAAR on the product space (=swarming), then
    # you will want to change product_space_dimension=2
    # and adjust your input files and projectors accordingly. 
    # you could also do this within the data processor

    def __init__(self, new_config={}):
        """
        The initialization of a ART instance takes the default configuration
        and updates the parameters with the arguments in new_config.
        
        Parameters
        ----------
        new_config : dict, optional - Parameters to initialize the problem. If unspecified, the default config is used.
        """
        self.config.update(new_config)

        #call data processor, read data

        module = __import__(self.config['data_filename'])
        data_processor = getattr(module, self.config['data_filename'])
        data_processor(self.config)

        self.config['iter']=0
        self.config['TOL2'] = 1e-15
        self.config['proj_iter']=1

        self.config['proxoperators'] = []
        self.config['proxoperators'].append(getattr(ProxOperators, self.config['Prox1']))
        self.config['proxoperators'].append(getattr(ProxOperators, self.config['Prox2']))
        Prox1 = self.config['proxoperators'][0](self.config)
        u_1 = Prox1.work(self.config['u_0'])
        self.config['proj_iter']=2
        Prox2 = self.config['proxoperators'][1](self.config)
        u_2 = Prox2.work(u_1)
        # estimate the gap in the relevant metric
        tmp_gap=0
        if (self.config['Ny']==1) or (self.config['Nx']==1) or (self.config['Nz']==1):
            tmp_gap = (norm(u_1-u_2,'fro')/(self.config['norm_data']))**2
        else:
            for j in range(self.config['product_space_dimension']):
                # compute (||P_Sx-P_Mx||/norm_data)^2:
                tmp_gap = tmp_gap+(norm(u_1[:,j]-u_2[:,j],'fro')/(self.config['norm_data']))**2

        gap_0=sqrt(tmp_gap)
        # sets the set fattening to be a percentage of the
        # initial gap to the unfattened set with 
        # respect to the relevant metric (KL or L2), 
        # that percentage given by
        # input.data_ball input by the user.
        self.config['data_ball'] = self.config['data_ball']*gap_0
        # the second tolerance relative to the oder of 
        # magnitude of the metric
        self.config['TOL2'] = self.config['data_ball']*1e-15

        self.config['norm_data'] = self.config['norm_data'] #copy this is needed by AP
        self.algorithm = getattr(Algorithms, self.config['algorithm'])(self.config);


    def _presolve(self):
        """
        Prepares argument for actual solving routine
        """
        
    
    def _solve(self):
        """
        Runs the algorithm to solve the given sudoku problem
        """
        #algorithm = self.config['algorithm'](self.config)
        
        self.output = self.algorithm.run(self.config['u_0'],self.config['TOL'],self.config['MAXIT'])
        print('Iterations:' + str(self.output['iter']))

    
    def _postsolve(self):
        """
        Processes the solution and generates the output
        """
        
        
        
    def show(self):
        """
        Generates graphical output from the solution
        """
        
        print("Calculation time:")
        print(self.elapsed_time)
        graphics = getattr(Graphics,self.config['graphics_display'])
        graphics(self.config,self.output)

