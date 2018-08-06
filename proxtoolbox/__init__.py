# -*- coding: utf-8 -*-

""" 
  ProxToolbox
  ~~~~~~~~~~~

  ..

  :Copyright: (c) 2015 by `D. Russell Luke <http://num.math.uni-goettingen.de/~r.luke/>`_
  :License: BSD, see LICENSE for more details.
  :Contributors:
   * Robert Hesse, Inst. for Numerical and Applied Math, Universität Göttingen 
   * Pär Mattson, Inst. for Numerical and Applied Math, Universität Göttingen (Ptychography)
   * Rebecca Nahme, Inst. for Numerical and Applied Math, Universität Göttingen
   * Alexander Stalljahn, Inst. for Numerical and Applied Math, Universität Göttingen (Sudoku)
   * `Matthew Tam <http://num.math.uni-goettingen.de/~mtam/>`_, CARMA, University of Newcastle, Australia (Ptychography)
   * Robin Wilke, Inst. for X-Ray Physics, Univesität Göttingen (Ptychography and Phase)
   * Stefan Ziehe, Inst. for Computer Science, Universität Göttingen.
  :Funding:	 This has grown over the years and has been supported in part by:
  
       * NASA grant NGT5-66
       * Pacific Institute for Mathematical Sciences (PIMS)		
       * USA NSF Grant DMS-0712796
       * German DFG grant SFB-755 TPC2

"""

from .Algorithms import *
from .Problems import *
from .ProxOperators import *
from .Utilities import *

__all__ = ["Algorithms","Problems","ProxOperators","Utilities"]

