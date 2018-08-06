# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 15:12:47 2016

@author: rebecca

How to start:
open a terminal, go to the ProxPython directory,
start ipython, then type 'run TestSuite/TestSuite_all'
"""

import unittest

from TestSuite.Testing_Problems import Test_Problems
from TestSuite.Testing_Utilities import Test_Utilities

def Suite():
    """
    Test suite containing all tests implemented (if added)
    """
    suite = unittest.TestSuite()
    
    # makeSuite(): convenience function constructing a test suite that comprises all test cases in a test case class
    # passing the test cases to the test suite
    suite.addTest(unittest.makeSuite(Test_Problems))
    suite.addTest(unittest.makeSuite(Test_Utilities))
    
    return suite
    
    
if __name__ == "__main__":
    """
    Handler to provide a process for running the test suite
    """
    start = unittest.TextTestRunner()
    testsuite = Suite()
    start.run(testsuite)
