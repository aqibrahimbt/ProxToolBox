# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 10:15:22 2014

@author: stefan
"""

import numpy, matplotlib.pyplot as pyplot;

from problem import Problem
from algorithms import RAAR
from proxoperators import ProxOperator, P_diag, P_parallel

__all__ = ["Sudoku"];




class ProjGiven(ProxOperator):
    """
    Projection onto the given entries in a Sudoku problem
    
    See Also
    --------
        ProxOperator : Generic interface for prox operators
    """
    
    def __init__(self,config):
        """
        Parameters
        ----------
        config : dict
            Dictionary containing the problem configuration. It must have the
            key 'given_sudoku' mapping to the input Sudoku.
        """
        self.given = config['given_sudoku'].copy();
        

    def work(self,u):
        A = numpy.zeros((9,9,9),dtype=u.dtype);
        given = self.given;
        
        for x in range(9):
            for y in range(9):
                z = given[x,y];
                if z > 0:
                    A[x,y,z-1] = 1;
                else:
                    A[x,y,numpy.argmax(u[x,y,:])] = 1;
        
        return A;


class ProjSquare(ProxOperator):
    """
    Projection onto the box constraints in a Sudoku problem
    
    See Also
    --------
        ProxOperator : Generic interface for prox operators
    """
    def __init__(self,config=None):
        """
        Parameters
        ----------
        config : dict, optional
            Not used here.
        """
        return;
        
    def work(self,u):
        Q = numpy.zeros((9,9,9),dtype=u.dtype);
        for z in range(9):
            for x in range(0,9,3):
                for y in range(0,9,3):
                    v = numpy.argmax(u[x:(x+3),y:(y+3),z],axis=0);
                    w = numpy.argmax(numpy.amax(u[x:(x+3),y:(y+3),z],axis=0));
                    
                    Q[x+v[w],y+w,z] = 1;
        return Q;


class ProjColumn(ProxOperator):
    """
    Projection onto the column constraints in a Sudoku problem
    
    See Also
    --------
        ProxOperator : Generic interface for prox operators
    """
    def __init__(self,config):
        """
        Parameters
        ----------
        config : dict, optional
            Not used here.
        """
        return;
    
    def work(self,u):
        C = numpy.zeros((9,9,9),dtype=u.dtype);
        for x in range(9):
            for z in range(9):
                y = numpy.argmax(u[x,:,z]);
                C[x,y,z] = 1;
        return C;


class ProjRow(ProxOperator):
    """
    Projection onto the row constraints in a Sudoku problem
    
    See Also
    --------
        ProxOperator : Generic interface for prox operators
    """
    
    def __init__(self,config):
        """
        Parameters
        ----------
        config : dict, optional
            Not used here.
        """
        return;
    
    def work(self,u):
        R = numpy.zeros((9,9,9),dtype=u.dtype);
        for y in range(9):
            for z in range(9):
                x = numpy.argmax(u[:,y,z]);
                R[x,y,z] = 1;
        return R;



class Sudoku(Problem):
    
    default_config = {
        # This is the algorithm we use. RAAR and HPR will work.
        'algorithm':RAAR,
        # RAAR requires 2 ProxOperators
        'proj1':P_diag,
        'proj2':P_parallel,
        # P_parallel requires a sequence of projectors
        'projectors':(ProjRow,ProjColumn,ProjSquare,ProjGiven),
        # Relaxation parameters for RAAR/HPR
        'beta0':1,
        'beta_max':1,
        'beta_switch':1,
        # Any algorithm requires these
        'maxiter':2000,
        'tol':1e-9,
        # Just a random Sudoku. Not too easy, but no challenge for
        # the mighty ProxToolbox!
        'given_sudoku':numpy.array(((2,0,0,0,0,1,0,3,0),
                                    (4,0,0,0,8,6,1,0,0),
                                    (0,0,0,0,0,0,0,0,0),
                                    (0,0,0,0,1,0,0,0,0),
                                    (0,0,0,0,0,0,9,0,0),
                                    (0,0,5,0,0,3,0,0,7),
                                    (0,0,0,0,0,0,0,0,0),
                                    (1,0,0,0,0,7,4,9,0),
                                    (0,2,4,1,0,0,0,0,0)),dtype=numpy.float32)
    };
    
    def __init__(self, config=default_config):
        """
        Parameters
        ----------
        config : dict, optional
            Dictionary containing the problem configuration. If unspecified,
            Sudoku.default_config is used.
        """
        #config['projectors'] = (ProjRow,ProjColumn,ProjSquare,ProjGiven)
        # These are always the same for Sudoku, so they don't need to be
        # in config
        self.config = config.copy();
        self.config['Nx'] = self.config['Ny'] = self.config['Nz'] = 9;
        self.config['dim'] = 4; self.config['norm_data'] = 81;
        self.config['given_sudoku'] = config['given_sudoku'].copy();
        
        self.algorithm = self.config['algorithm'](self.config);
        self.sudoku = self.config['given_sudoku'];
    
    
    def _presolve(self):
        sudoku = self.sudoku;
        
        u = numpy.zeros((9,9,9,4),dtype=sudoku.dtype);
        
        for x in range(9):
            for y in range(9):
                z = sudoku[x,y]-1;
                if z >= 0:
                    u[x,y,z,:] = 1;
        
        self.u = u;
    
    
    def _solve(self):
        self.u1,self.u2,self.iters,self.change,self.gap = \
            self.algorithm.run(self.u,self.config['tol'],self.config['maxiter']);
    
    
    def _postsolve(self):
        solution = numpy.zeros_like(self.sudoku);
        A = self.u1[:,:,:,0];
        
        for x in range(9):
            for y in range(9):
                for z in range(9):
                    if A[x,y,z] > 0:
                        solution[x,y] = z+1;
                        break;
        
        self.solution = solution;
        
        
        fig = pyplot.figure('Sudoku');
        
        ax = pyplot.subplot(2,2,1);
        ax.title.set_text('Given Sudoku');
        ax.xaxis.set_visible(False);
        ax.yaxis.set_visible(False);
        table = ax.table(cellText=self.sudoku.astype(numpy.int32),loc='center');
        for cell in table.properties()['child_artists']:
            cell.set_height(0.1);
            cell.set_width(0.1);
            txt = cell.get_text();
            if txt.get_text() == '0':
                txt.set_text('');
        
        ax = pyplot.subplot(2,2,2);
        ax.title.set_text('Solution');
        ax.xaxis.set_visible(False);
        ax.yaxis.set_visible(False);
        table = ax.table(cellText=self.solution.astype(numpy.int32),loc='center');
        for cell in table.properties()['child_artists']:
            cell.set_height(0.1);
            cell.set_width(0.1);
        
        ax = pyplot.subplot(2,2,3);
        ax.xaxis.label.set_text('Iterations');
        ax.yaxis.label.set_text('change');
        pyplot.semilogy(self.change);
        
        ax = pyplot.subplot(2,2,4);
        ax.xaxis.label.set_text('Iterations');
        ax.yaxis.label.set_text('gap');
        pyplot.semilogy(self.gap);
        
        fig.show();
