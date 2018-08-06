# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 10:40:18 2014

@author: stefan
"""

__all__ = ["Algorithm"]

class Algorithm:
    """
    Generic interface for ProxToolbox algorithms.
    It predefines the standard methods (names, input data).
    """
    
    def __init__(self, config):
        """
        Initialization method for a concrete instance.
        
        Parameters
        ----------
        config : dict
                 Parameters to configure the algorithm
                 
        """
        raise NotImplementedError ("This is just an abstract interface")
        
    def run(self, u, tol, maxiter):
        """
        Runs the algorithm for the specified input data

        Parameters
        ----------
        u : array_like - Input data
        tol : float - Tolerance
        maxiter : int - Maximum number of iterations
        
        Returns
        -------
        (u1, u2, iters, change, gap) : tuple consisting of
                * u1 : array_like - Result
                * u2 : array_like - Result
                * iters : int - Number of iterations the algorithm performed
                * change : array_like - Percentage change in norm
                * gap : array_like - Squared gap distance normalized by the magnitude constraint
        
        """
        raise NotImplementedError ("This is just an abstract interface")
