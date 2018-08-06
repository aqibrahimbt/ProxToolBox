# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 13:49:37 2016

@author: rebecca

This class is a collection of tests to check integrity of various aspects of the Problems-module.
How to start:
open a terminal and ipython, then type 1. 'import Testing_Problems' and 2. 'run Testing_Problems'
"""

import unittest

from numpy import array, float32

from proxtoolbox.Algorithms import *
from proxtoolbox.Problems import *

from sudoku_in import new_config

class Test_Problems(unittest.TestCase):
    """
    Collection of various (unit) tests for the Problems-module.
    """      
    def test_instantiation(self):
        """
        Tests if after the initialization of an instance,
        each object is an instance of the respective Problem-class.
        """
        
        prob_sudoku = Sudoku()
        self.assertIsInstance(prob_sudoku, Sudoku)
        self.assertTrue("solve" in dir(prob_sudoku))
        
        prob_pty = Ptychography()
        self.assertTrue(isinstance(prob_pty, Ptychography))
        self.assertTrue("solve" in dir(prob_pty))
        
        
    def test_sudoku(self):
        """
        Tests the integrity of Problem 'Sudoku'
        by initializing a test data set, solving it and checking the result.
        The test data set is created such that there exists an unique solution.
        """
        
        config = {
            # This is the algorithm we use. RAAR and HPR will work.
            'algorithm':'RAAR',
            # RAAR requires 2 ProxOperators
            'proj1':'P_diag',
            'proj2':'P_parallel',
            # P_parallel requires a sequence of projectors
            'projectors':('ProjRow','ProjColumn','ProjSquare','ProjGiven'),
            # Relaxation parameters for RAAR/HPR
            'beta0':1,
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
            'dim':4,
            'norm_data':81,
            # Just a random Sudoku. Not too easy, but no challenge for
            # the mighty ProxToolbox!
            'sudoku':((0,4,5,7,0,0,0,6,1),
                      (0,0,0,8,0,5,4,0,0),
                      (0,0,7,6,4,0,0,0,2),
                      (5,0,0,9,6,4,0,0,0),
                      (6,0,8,5,0,1,0,3,0),
                      (0,7,0,0,0,0,0,9,0),
                      (3,0,2,4,0,0,1,0,0),
                      (9,0,0,1,0,0,0,0,6),
                      (0,0,0,0,0,3,8,0,0))
        }
        
        prob_sudoku = Sudoku(config)
        
        prob_sudoku.solve()
        self.assertEqual((prob_sudoku.config['sudoku'])[0,1], 4)
        self.assertEqual(prob_sudoku.solution[-1,-1], 9)
        
        self.assertTrue((prob_sudoku.solution == array(((8,4,5,7,3,2,9,6,1),
                                                        (2,6,9,8,1,5,4,7,3),
                                                        (1,3,7,6,4,9,5,8,2),
                                                        (5,2,3,9,6,4,7,1,8),
                                                        (6,9,8,5,7,1,2,3,4),
                                                        (4,7,1,3,2,8,6,9,5),
                                                        (3,8,2,4,9,6,1,5,7),
                                                        (9,5,4,1,8,7,3,2,6),
                                                        (7,1,6,2,5,3,8,4,9)),dtype=float32)).all())
                                                        
    def test_sudoku_configfile(self):
        """
        Tests the initialization of a Sudoku instance with a config provided in a file
        """
        prob_sudoku = Sudoku(new_config)
        
        prob_sudoku.solve()
        
        self.assertTrue((prob_sudoku.solution == array(((8,4,5,7,3,2,9,6,1),
                                                        (2,6,9,8,1,5,4,7,3),
                                                        (1,3,7,6,4,9,5,8,2),
                                                        (5,2,3,9,6,4,7,1,8),
                                                        (6,9,8,5,7,1,2,3,4),
                                                        (4,7,1,3,2,8,6,9,5),
                                                        (3,8,2,4,9,6,1,5,7),
                                                        (9,5,4,1,8,7,3,2,6),
                                                        (7,1,6,2,5,3,8,4,9)),dtype=float32)).all())
    
    
    def test_ptychography(self):
        """
        Tests the integrity of Problem 'Ptychography'
        """
        
        
        config = {
            'sim_data_type':'gaenseliesel',
            'Nx':64,
            'Ny':64,
            'Nz':1,
            'scan_stepsize':3.5e-7,
            'nx':25,
            'ny':25,
            'sample_area_center':(0.5,0.5),
            'noise':None,
            'switch_probemask':True,
            'probe_mask_gamma':1,
            'rmsfraction':0.5,
            'scan_type':'raster',
            'switch_object_support_constraint':True,
            'probe_guess_type':'circle',
            'object_guess_type':'random',
            'warmup':True,
            'warmup_alg':'PALM',
            'warmup_maxiter':1,
            'algorithm':'PALM',
            'maxiter':1, # This is usually 300. I reduced it to be able to test it.
            'tol':0.625,
            'ignore_error':False,
            'ptychography_prox':'Rodenburg',
            'blocking_switch':True,
            'blocking_scheme':'divide',
            'between_blocks_scheme':'sequential',
            'within_blocks_scheme':'parallel',
            'block_rows':2,
            'block_cols':2,
            'block_maxiter':1,
            'rodenburg_inner_it':1,
            'beta0':1,
            'beta_max':1,
            'beta_switch':30,
            'bs_mask':None,
            'bs_factor':1,
            'fmask':None,
            'overrelax':1
        }
        
        prob_pty = Ptychography(config)
        prob_pty.solve()
        
        self.assertTrue(prob_pty.u_final.dtype == 'complex128')
        # what to take for integrity-check?
            

if __name__ == "__main__":
    unittest.main()
