# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 10:15:22 2014

@author: stefan
"""

from numpy import amax, argmax, array, float32, int32, zeros, zeros_like
from matplotlib import pyplot

from .problems import Problem
from proxtoolbox import Algorithms
from proxtoolbox import ProxOperators
from proxtoolbox.ProxOperators.proxoperators import ProxOperator

#__all__ = ["Sudoku"]


class ProjGiven(ProxOperator):
    """
    Projection onto the given entries in a Sudoku problem
    
    """
    
    def __init__(self,config):
        """
        Initialization method
        
        Parameters
        ----------
        config : dict - Dictionary containing the problem configuration. It must have the key 'sudoku' mapping to the input Sudoku.
        """
        self.given = config['sudoku'].copy()
        

    def work(self,u):
        """
        Applies the prox operator to the input data
        
        Parameters
        ----------
        u : array_like - Input data for the operator
        
        Returns
        -------
        array_like - Result of the operation
        """
        A = zeros((9,9,9),dtype=u.dtype)
        
        for x in range(9):
            for y in range(9):
                z = self.given[x,y]
                if z > 0:
                    A[x,y,z-1] = 1
                else:
                    A[x,y,argmax(u[x,y,:])] = 1
        
        return A


class ProjSquare(ProxOperator):
    """
    Projection onto the box constraints in a Sudoku problem
    
    """
    def __init__(self,config=None):
        """
        Initialization
        
        Parameters
        -----------
        config : dict, optional - Not used here
        """
        return
        
    def work(self,u):
        """
        Applies the prox operator to the input data
        
        Parameters
        ----------
        u : array_like - Input data for the operator
        
        Returns
        -------
        array_like - Result of the operation
        """
        Q = zeros((9,9,9),dtype=u.dtype)
        for z in range(9):
            for x in range(0,9,3):
                for y in range(0,9,3):
                    v = argmax(u[x:(x+3),y:(y+3),z],axis=0)
                    w = argmax(amax(u[x:(x+3),y:(y+3),z],axis=0))
                    
                    Q[x+v[w],y+w,z] = 1
        return Q


class ProjColumn(ProxOperator):
    """
    Projection onto the column constraints in a Sudoku problem
    
    """
    def __init__(self,config):
        """
        Initialization
        
        Parameters
        ----------
        config : dict, optional - Not used here
        """
        return
    
    def work(self,u):
        """
        Applies the prox operator to the input data
        
        Parameters
        ----------
        u : array_like - Input data for the operator
        
        Returns
        -------
        array_like - Result of the operation
        """
        C = zeros((9,9,9),dtype=u.dtype)
        for x in range(9):
            for z in range(9):
                y = argmax(u[x,:,z])
                C[x,y,z] = 1
        return C


class ProjRow(ProxOperator):
    """
    Projection onto the row constraints in a Sudoku problem
    
    """
    
    def __init__(self,config):
        """
        Initialization
        
        Parameters
        ----------
        config : dict, optional - Not used here
        """
        return
    
    def work(self,u):
        """
        Applies the prox operator to the input data
        
        Parameters
        ----------
        u : array_like - Input data for the operator
        
        Returns
        -------
        array_like - Result of the operation
        """
        R = zeros((9,9,9),dtype=u.dtype)
        for y in range(9):
            for z in range(9):
                x = argmax(u[:,y,z])
                R[x,y,z] = 1
        return R



class Sudoku(Problem):
    """
    Sudoku Problem
    
    The goal of a standard Sudoku puzzle is to fill a 9x9 array with the numbers
    from 1 to 9. Every row, every column and each of the 9 3x3 subblocks should
    contain each number only once.
    
    Starting point is a grid that already contains several numbers.
    Usually, there exists a unique solution.
    
    """
    config = {
        # This is the algorithm we use. RAAR and HPR will work.
        'algorithm':'RAAR',
        # RAAR requires 2 ProxOperators
        'proxoperators':('P_diag','P_parallel'),
        # P_parallel requires a sequence of projectors
        'projectors':('ProjRow','ProjColumn','ProjSquare','ProjGiven'),
        # Relaxation parameters for RAAR/HPR
        'beta_0':1,
        'beta_max':1,
        'beta_switch':1,
        # Any algorithm requires these
        'maxiter':2000,
        'tol':1e-9,
        # Dimension parameters
        # which are the same for every standard Sudoku
        'Nx':9,
        'Ny':9,
        'Nz':9,
        'product_space_dimension':4,
        'norm_data':81,
        # Just a random Sudoku. Not too easy, but no challenge for
        # the mighty ProxToolbox!
        'sudoku':((2,0,0,0,0,1,0,3,0),
                  (4,0,0,0,8,6,1,0,0),
                  (0,0,0,0,0,0,0,0,0),
                  (0,0,0,0,1,0,0,0,0),
                  (0,0,0,0,0,0,9,0,0),
                  (0,0,5,0,0,3,0,0,7),
                  (0,0,0,0,0,0,0,0,0),
                  (1,0,0,0,0,7,4,9,0),
                  (0,2,4,1,0,0,0,0,0)),
        'diagnostic': True
    }
    
    def __init__(self, new_config={}):
        """
        The initialization of a Sudoku instance takes the default configuration
        and updates the parameters with the arguments in new_config.
        
        Parameters
        ----------
        new_config : dict, optional - Parameters to initialize the problem. If unspecified, the default config is used.
        """
        self.config.update(new_config)
        
        add_config = self.config.copy()
        add_config['proxoperators'] = []
        for prox in self.config['proxoperators']:
            add_config['proxoperators'].append(getattr(ProxOperators, prox))
#        add_config.update({'prox1':getattr(ProxOperators, self.config['proj1']),
#                           'prox2':getattr(ProxOperators, self.config['proj2'])})
        
#        self.config['algorithm'] = getattr(Algorithms, self.config['algorithm'])
#        self.config['proj1'] = getattr(ProxOperators, self.config['proj1'])
#        self.config['proj2'] = getattr(ProxOperators, self.config['proj2'])
        
        add_config['projectors'] = []
        for p in self.config['projectors']:
            add_config['projectors'].append(globals()[p])
            
        add_config['sudoku'] = array(self.config['sudoku'],dtype=float32)
        
        self.algorithm = getattr(Algorithms, self.config['algorithm'])(add_config)
        
        self.config['sudoku'] = array(self.config['sudoku'],dtype=float32)
    
    
    def _presolve(self):
        """
        Prepares argument for actual solving routine
        """
        u = zeros((9,9,9,4),dtype=self.config['sudoku'].dtype)
        
        for x in range(9):
            for y in range(9):
                z = self.config['sudoku'][x,y]-1
                if z >= 0:
                    u[x,y,z,:] = 1
        
        self.u = u
    
    
    def _solve(self):
        """
        Runs the algorithm to solve the given sudoku problem
        """
#        algorithm = self.config['algorithm'](self.config)
        
        self.output = self.algorithm.run(self.u,self.config['tol'],self.config['maxiter'])
    
    
    def _postsolve(self):
        """
        Processes the solution and generates the output
        """
        solution = zeros_like(self.config['sudoku'])
        A = self.output['u1'][:,:,:,0]
        
        for x in range(9):
            for y in range(9):
                for z in range(9):
                    if A[x,y,z] > 0:
                        solution[x,y] = z+1
                        break
        
        self.solution = solution
        
        
    def show(self):
        """
        Generates graphical output from the solution
        """
        fig = pyplot.figure('Sudoku')
        
        # plot plain Sudoku puzzle
        ax = pyplot.subplot(2,2,1)
        ax.title.set_text('Given Sudoku')
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        table = ax.table(cellText=self.config['sudoku'].astype(int32),loc='center')
        for cell in table.properties()['child_artists']:
            cell.set_height(0.1)
            cell.set_width(0.1)
            txt = cell.get_text()
            if txt.get_text() == '0':
                txt.set_text('')
                
        # plot solution
        ax = pyplot.subplot(2,2,2)
        ax.title.set_text('Solution')
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        table = ax.table(cellText=self.solution.astype(int32),loc='center')
        for cell in table.properties()['child_artists']:
            cell.set_height(0.1)
            cell.set_width(0.1)
        
        # plot the change from one iterate to the next
        ax = pyplot.subplot(2,2,3)
        ax.xaxis.label.set_text('Iterations')
        ax.yaxis.label.set_text('Change')
        pyplot.semilogy(self.output['change'])
        
        # plot the gap
        ax = pyplot.subplot(2,2,4)
        ax.xaxis.label.set_text('Iterations')
        ax.yaxis.label.set_text('Gap')
        pyplot.semilogy(self.output['gap'])
        
        fig.tight_layout()
        pyplot.show()
