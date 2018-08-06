# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 14:02:02 2016

@author: rebecca

The "ProxOperators"-module contains all proximal operators that are already implemented in the ProxToolbox.
"""

from .proxoperators import *

__all__ = ["P_diag","P_parallel","magproj", "Approx_P_JWST_Poisson", "P_amp", "P_SP", "Approx_PM_Gaussian", "Approx_PM_Poisson", "P_S", "P_S_real", "P_sequential_hyperplane_odd", "P_sequential_hyperplane_even","P_parallel_hyperplane", "P_block_parallel_hyperplane", "P_block_sequential_hyperplane", "Q_Heau", "P_Amod", "Approx_P_FreFra_Poisson"]
