# -*- coding: utf-8 -*-

from proxtoolbox.Problems.problems import Problem
from proxtoolbox import Algorithms
from proxtoolbox import ProxOperators
from proxtoolbox.ProxOperators.proxoperators import ProxOperator
from proxtoolbox.Problems.Phase import Graphics
from numpy.linalg import norm
import numpy as np
import h5py
from numpy import square, sqrt, nonzero, size


class Phase(Problem):
    """
    Phase Problem
    """
    config = {
    }
    
    def __init__(self, new_config={}):
        """
        The initialization of a Phase instance takes the default configuration
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

        # reshape and rename the data
        if 'data_sq' not in self.config:
            self.config['data_sq'] = self.config['data']
            self.config['data'] = self.config['rt_data']
            if('norm_data' in self.config):
                self.config['norm_data_sq']= self.config['norm_data']
            self.config['norm_data']=self.config['norm_rt_data']

            tmp = self.config['data'].shape
            if(tmp[0]==1 or tmp[1]==1):
                    self.config['data_sq'] = self.config['data_sq'].reshape((self.config['Nx'],self.config['Ny']))
                    #the prox algorithms work with the square root of the measurement:
                    self.config['data'] = self.config['data'].reshape((self.config['Nx'],self.config['Ny']))

            if 'Nz' not in self.config:
                self.config['Nz'] = 1


        #If method_config[formulation is does not exist, i.e. not specified in 
        #the *_in.m file, use the product space as the default.
        if 'formulation' in self.config:
            formulation = self.config['formulation']
        else:
            formulation = 'product space'

        # Set the projectors and inputs based on the types of constraints and 
        # experiments
        proxoperators = ['','','']

        if self.config['constraint'] == 'hybrid':
            proxoperators[0] = 'P_cP' # This will be problem specific
        elif self.config['constraint'] == 'support only':
            proxoperators[0] = 'P_S'
        elif self.config['constraint'] == 'real and support':
            proxoperators[0] ='P_S_real'
        elif self.config['constraint'] =='nonnegative and support':
            proxoperators[0] ='P_SP'
        elif self.config['constraint'] =='amplitude only':
            proxoperators[0] ='P_amp'
        elif self.config['constraint'] == 'phase on support':
            proxoperators[0] ='P_Amod'
        elif self.config['constraint'] =='minimum amplitude':
            proxoperators[0] = 'P_min_amp'
        elif self.config['constraint'] =='sparse':
            proxoperators[0] = 'not in yet' 
        elif self.config['constraint'] =='phaselift':
            proxoperators[0] = 'P_mean_SP'
        elif self.config['constraint'] =='phaselift2':
            proxoperators[0] ='P_liftM'
            proxoperators[2] ='Approx_PM_Poisson' # Patrick: This is just to monitor the change of phases!  

        if self.config['experiment'] in ['single diffraction', 'CDI']:
            if self.config['distance'] == 'far field':
                if self.config['constraint'] == 'phaselift':
                    proxoperators[1] = 'P_Rank1'
                elif self.config['constraint'] == 'phaselift2':
                    proxoperators[1] = 'P_rank1_SR'
                else:
                    proxoperators[1] = (
                        'Approx_PM_Poisson'
                        if self.config['noise'] == 'Poisson'
                        else 'Approx_PM_Gaussian'
                    )
            elif self.config['noise'] == 'Poisson':
                proxoperators[1]='Approx_P_FreFra_Poisson'
        elif self.config['experiment'] in ['Krueger', 'Near_field_cell_syn']:
            proxoperators[0]='P_Amod'
            proxoperators[1]='Approx_P_FreFra_Poisson'
        elif self.config['experiment'] in ('dict', 'dictyM103_stx6_600frames','xenopus','living_worm'):
            proxoperators[0]='P_Amod' # Not sure this is the appropriate prox operator for these
                                   # experiments...
            proxoperators[1]='Approx_P_FreFra_Poisson'

        elif self.config['experiment'] == 'diversity diffraction' and formulation == 'sequential':
            proxoperators[1] = 'Approx_P_RCAAR_JWST_Poisson'
            proxoperators[0] = proxoperators[1]
        elif self.config['experiment'] == 'JWST': 
            proxoperators[1] = 'Approx_P_JWST_Poisson'  
            proxoperators[2] = proxoperators[0]
            proxoperators[0] = 'P_diag'
        elif self.config['experiment'] == 'CDP':
            proxoperators[1] = 'P_CDP'
            proxoperators[2] = proxoperators[0]
            proxoperators[0] = 'P_diag'
        elif self.config['experiment'] == 'ptychography':
            proxoperators[1] = 'not in yet'
        elif self.config['experiment'] == 'complex':
            proxoperators[1] = 'not in yet'
        elif self.config['constraint'] == 'phaselift':
            proxoperators[1] ='P_PL_lowrank'

        self.config['proxoperators'] = [
            getattr(ProxOperators, prox) for prox in proxoperators if prox != ''
        ]
        # input.Proj1_input.F=F;  % is it any more expensive to pass everything
        # into the projectors rather than just a selection of data and
        # parameters?  If not, and we pass everything anyway, there is no need
        # to create a new structure element.

        if 'product_space_dimension' not in self.config:
            self.config['product_space_dimension'] = 1

        # set the animation program:
        self.config['animation']='Phase_animation'
        #
        # if you are only working with two sets but
        # want to do averaged projections
        # (= alternating projections on the product space)
        # or RAAR on the product space (=swarming), then
        # you will want to change product_space_dimension=2
        # and adjust your input files and projectors accordingly. 
        # you could also do this within the data processor

        self.config['TOL2'] = 1e-15

        #To estimate the gap in the sequential formulation, we build the
        # appropriate point in the product space. This allows for code reuse.
        # Note for sequential diversity diffraction, input.Proj1 is the "RCAAR"
        # version of the function.
        if formulation == 'sequential':
            for j in range(self.config['product_space_dimension']):
                self.config['proj_iter'] =j
                proj1 = self.config['proxoperators'][0](self.config)
                u_1[:,:,j]= proj1.work(self.config['u_0'])
                self.config['proj_iter'] = mod(j,config['product_space_dimension'])+1
                proj1 = self.config['proxoperators'][0](self.config)
                u_1[:,:,j]= proj1.work(self.config['u_0'])
            end
        else: #i.e. formulation=='product space'
            proj1 = self.config['proxoperators'][0](self.config)
            u_1 = proj1.work(self.config['u_0'])
            proj2 = self.config['proxoperators'][1](self.config)
            u_2 = proj2.work(u_1)

        # estimate the gap in the relevant metric
        if self.config['Nx'] ==1 or self.config['Ny']==1 :
            tmp_gap = square(norm(u_1-u_2)/self.config['norm_rt_data'])
        elif self.config['product_space_dimension'] == 1:
            tmp_gap = (norm(u_1-u_2)/self.config['norm_rt_data'])**2
        else:
            tmp_gap=0
            for j in range(self.config['product_space_dimension']):
                # compute (||P_Sx-P_Mx||/norm_data)^2:
                tmp_gap = tmp_gap+(norm(u_1[:,:,j]-u_2[:,:,j])/self.config['norm_rt_data'])**2

        gap_0=sqrt(tmp_gap)

        # sets the set fattening to be a percentage of the
        # initial gap to the unfattened set with 
        # respect to the relevant metric (KL or L2), 
        # that percentage given by
        # input.data_ball input by the user.
        self.config['data_ball']=self.config['data_ball']*gap_0
        # the second tolerance relative to the oder of 
        # magnitude of the metric
        self.config['TOL2'] = self.config['data_ball']*1e-15
        self.config['proxoperators']
        self.algorithm = getattr(Algorithms, self.config['algorithm'])(self.config)
        
    
    
    def _presolve(self):
        """
        Prepares argument for actual solving routine
        """
        
    
    def _solve(self):
        """
        Runs the algorithm to solve the given sudoku problem
        """
#        algorithm = self.config['algorithm'](self.config)
        
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

    def compare_to_matlab(self):
        """
        Routine to test and verify results by comparing to matlab
        Note that this is only for development and should not be used by a normal user
        For result to match u_0 should be chosen as np.multiply(config['abs_illumination'],exp(1j*2*pi*0.5*np.ones(newres)))']
        """

        print(self.config['proxoperators'])

        if self.config['experiment'] == 'JWST':
            if self.config['algorithm'] == 'RAAR':
                if self.config['MAXIT'] == 1:
                    f = h5py.File('Phase_test_data/u1_1.mat')
                elif self.config['MAXIT'] == 500 :
                    f = h5py.File('Phase_test_data/u1_500.mat')
                else:
                    print("No file available to compare to.")
                    return
            elif self.config['algorithm'] == 'AP':
                f = h5py.File('Phase_test_data/JWST_u1_ap_' + str(self.config['MAXIT']) + '.mat')

            u1 = f['u1'].value.view(np.complex)
        elif self.config['data_filename'] == 'Siemens_processor' and self.config['constraint'] == 'amplitude':
            f = h5py.File('Phase_test_data/siemens_amplitude_u1_' + str(self.config['MAXIT']) + '.mat')
            u1 = f['u1'].value.view(np.complex)
        elif self.config['data_filename'] == 'Siemens_processor' and self.config['constraint'] == 'nonnegative and support':
            f = h5py.File('Phase_test_data/siemens_nonneg_u1_' + str(self.config['MAXIT']) + '.mat')
            u1 = f['u1'].value.view(np.float64)
        elif self.config['data_filename'] == 'Siemens_processor' and self.config['constraint'] == 'real and support':
            f = h5py.File('Phase_test_data/siemens_real_u1_' + str(self.config['MAXIT']) + '.mat')
            u1 = f['u1'].value.view(np.float64)
        else:
            if self.config['algorithm'] == 'RAAR':
                if self.config['beta_0'] == 0.95:
                    if self.config['MAXIT'] == 1000 :
                        f = h5py.File('Phase_test_data/tasse_u1_1000.mat')
                    elif self.config['MAXIT'] == 20:
                        f = h5py.File('Phase_test_data/tasse_u1_20.mat')
                    elif self.config['MAXIT'] == 1:
                        f = h5py.File('Phase_test_data/tasse_u1_1.mat')
                    else:
                        print("No file available to compare to.")
                        return
                elif self.config['beta_0'] == 0.50:
                    f = h5py.File('Phase_test_data/tasse_u1_'+ str(self.config['MAXIT']) + '_beta_0_5.mat')
                else:
                    print("No file available to compare to.")
                    return
            elif self.config['algorithm'] == 'AP' and self.config['constraint'] == 'support only':
                        f = h5py.File('Phase_test_data/tasse_supp_u1_ap_' + str(self.config['MAXIT']) + '.mat')
            elif (
                self.config['algorithm'] in ['AP', 'AP_expert']
                and self.config['constraint'] == 'nonnegative and support'
            ):
                f = h5py.File('Phase_test_data/tasse_u1_ap_' + str(self.config['MAXIT']) + '.mat')

            u1 = f['u1'].value.view(np.float64)

        u1 =np.array(u1)
        u1 = u1.T
        print("Compare u1:")
        #print("Nonzero indices matlab:")
        #print(nonzero(u1))
        #print("Nonzero indices python:")
        #print(nonzero(self.output['u1']))
        print("Nonzero indices equal:")
        print(np.array_equal(nonzero(u1),nonzero(self.output['u1'])))
        #print("Nonzero values matlab:")
        #print(u1[nonzero(u1)])
        #print("Nonzero values python:")
        #print(self.output['u1'][nonzero(self.output['u1'])])
        #print("Difference at nonzero values:")
        #nonz = nonzero(u1)
        diff = u1 - self.output['u1']
        #print(diff[nonz])
        print("Maximum norm of difference:")
        print(np.amax(abs(diff)))
        print("Frobenius norm of difference:")
        print(norm(diff))
        print("Frobenius norm of matlab u1:")
        print(norm(u1))
        print("Frobenius norm of python u1:")
        print(norm(self.output['u1']))

    def compare_data_to_matlab(self):
        """
        Routine to test and verify results by comparing to matlab
        Note that this is only for development and should not be used by a normal user
        For result to match u_0 should be chosen as np.multiply(config['abs_illumination'],exp(1j*2*pi*0.5*np.ones(newres)))'] =
        """

        if self.config['data_filename'] == 'CDI_data_processor':
            f = h5py.File('Phase_test_data/CDI_data_processor_rt_data.mat')
            data_mat = f['rt_data'].value.view(np.float)
            data_py = self.config['rt_data']
        elif self.config['data_filename'] == 'Siemens_processor':
            f = h5py.File('Phase_test_data/Siemens_processor_truth.mat')
            data_mat = f['truth'].value.view(np.float)
            data_py = self.config['truth']

        data_mat =np.array(data_mat)
        data_mat = data_mat.T
        print("Compare:")
        print("Nonzero indices equal:")
        print(np.array_equal(nonzero(data_mat),nonzero(data_py)))
        diff = data_mat - data_py
        print("Maximum norm of difference:")
        print(np.amax(abs(diff)))
        print("Frobenius norm of difference:")
        print(norm(diff))
        print("Frobenius norm of matlab data:")
        print(norm(data_mat))
        print("Frobenius norm of python data:")
        print(norm(data_py))

        if self.config['data_filename'] == 'CDI_data_processor':
            f = h5py.File('Phase_test_data/CDI_data_processor_S.mat')
            S = f['S'].value.view(np.float)

            S =np.array(S)
            S = S.T
            print("Compare S:")
            print("Nonzero indices equal:")
            print(np.array_equal(nonzero(S),nonzero(self.config['abs_illumination'])))
            diff = S - self.config['abs_illumination']
            print("Maximum norm of difference:")
            print(np.amax(abs(diff)))
            print("Frobenius norm of difference:")
            print(norm(diff))
            print("Frobenius norm of matlab S:")
            print(norm(S))
            print("Frobenius norm of python S:")
            print(norm(self.config['abs_illumination']))
        
