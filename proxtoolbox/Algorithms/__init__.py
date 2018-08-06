# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 12:53:43 2015

@author: rebecca

The "Algorithms"-module contains all algorithms provided by the ProxToolbox.
"""

#from algorithms import Algorithm
from .AP import *
from .AP_expert import *
from .HPR import *
# from .PALM import *
from .RAAR import *
from .RAAR_expert import *
from .GRAAL import *
from .HAAR import *
from .QNAP import *

__all__ = ["AP","HPR","RAAR", "AP_expert", "GRAAL", "HAAR", "RAAR_expert", "QNAP"]
