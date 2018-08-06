# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 10:52:12 2016

@author: rebecca

This class is a collection of tests to check integrity of functions in the Utilities-module.
How to start:
open a terminal, go to the directory containing this file,
start ipython, then type 1. 'import Testing_Utilities' and 2. 'run Testing_Problems'
"""

import unittest

from proxtoolbox import Utilities

class Test_Utilities(unittest.TestCase):  
    
    def test_Laplacian(self):
        
        Lap = Utilities.Laplace_matrix(3)
        
        self.assertTrue((Lap.toarray() == ([[2.,-1.,0.,-1.,0.,0.,0.,0.,0.],
                                            [-1.,3.,-1.,0.,-1.,0.,0.,0.,0.],
                                            [0.,-1.,2.,0.,0.,-1.,0.,0.,0.],
                                            [-1.,0.,0.,3.,-1.,0.,-1.,0.,0.],
                                            [0.,-1.,0.,-1.,4.,-1.,0.,-1.,0.],
                                            [0.,0.,-1.,0.,-1.,3.,0.,0.,-1.],
                                            [0.,0.,0.,-1.,0.,0.,2.,-1.,0.],
                                            [0.,0.,0.,0.,-1.,0.,-1.,3.,-1.],
                                            [0.,0.,0.,0.,0.,-1.,0.,-1.,2.]])).all())
        

if __name__ == "__main__":
    unittest.main()
