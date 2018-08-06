# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 11:36:27 2014

@author: stefan
"""

import time


__all__ = ["Problem"]

class Problem:
    """
    Generic interface for ProxToolbox problems
    """
    
    def __init__(self, config={}):
        """
        Initialization method
        
        Parameters
        ----------
        config : dict, optional
                 Dictionary containing problem parameters.
                 If unspecified, a default configuration will be used.
        """
        raise NotImplementedError("This is just an abstract interface")
    
    def _presolve(self):
        """
        Preprocessing data.
        This procedure is called before the problem will be solved.
        """
        raise NotImplementedError("This is just an abstract interface")

    def _solve(self):
        """
        This procedure solves the problem.
        """
        raise NotImplementedError("This is just an abstract interface")
        
    def _postsolve(self):
        """
        Postprocessing data.
        This procedure is called after the problem has been solved.
        """
        raise NotImplementedError("This is just an abstract interface")
        
    def solve(self):
        """
        Solves the problem, as specified by _presolve, _solve and _postsolve.
        
        This is NOT intended to be overloaded. Overload _presolve, _solve and
        _postsolve instead.
        """
        self._presolve();
        
        t = time.time();
        self._solve();
        self.elapsed_time = time.time() - t;
        
        self._postsolve();
        
    def show(self):
        """
        This procedure realizes textual and/or graphical output of the Problem's result.
        """
        raise NotImplementedError("This is just an abstract interface")

    def save(self):
        """
        This procedure stores input and output data.
        """
        raise NotImplementedError("This is just an abstract interface")
