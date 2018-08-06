# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 11:36:27 2014

@author: stefan
"""

__all__ = ["Problem"]


class Problem:
    """
    Interface for ProxToolbox problems
    """
    
    def __init__(self, config={}):
        """
        Parameters
        ----------
        config : dict, optional
            Dictionary containing problem parameters. If unspecified, a
            default configuration will be used.
        """
        raise NotImplementedError,"This is just an abstract interface"
    
    def _presolve(self):
        """
        This procedure is called before the problem will be solved.
        """
        raise NotImplementedError,"This is just an abstract interface"

    def _solve(self):
        """
        This procedure solves the problem.
        """
        raise NotImplementedError,"This is just an abstract interface"
        
    def _postsolve(self):
        """
        This procedure is called after the problem has been solved.
        """
        raise NotImplementedError,"This is just an abstract interface"
        
    def solve(self):
        """
        Solves the problem, as specified by _presolve, _solve and _postsolve.
        
        This is NOT intended to be overloaded. Overload _presolve, _solve and
        _postsolve instead.
        """
        self._presolve();
        self._solve();
        self._postsolve();
        
