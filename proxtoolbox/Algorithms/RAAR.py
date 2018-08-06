#                    RAAR.m
#             written on Aug.18 , 2017 by
#                   Russell Luke
#   Inst. Fuer Numerische und Angewandte Mathematik
#
# DESCRIPTION:  User-friendly version of the Relaxed Averaged Alternating Reflection algorithm.  
# 		For background see:
#                D.R.Luke, Inverse Problems 21:37-50(2005)
#                D.R. Luke, SIAM J. Opt. 19(2): 714--739 (2008).

from .SimpleAlgortihm import SimpleAlgorithm
from numpy import exp

class RAAR(SimpleAlgorithm):

    def evaluate(self, u):
        iter = self.config['iter']+1 #add 1 to conform with matlab version (iter influences beta and we sometimes want to compare results between python and matlab), matlab counts starting at 1, python starts at 0
        beta0 = self.config['beta_0']
        beta_max = self.config['beta_max']
        beta_switch = self.config['beta_switch']
        tmp1 = 2*self.prox2.work(u)-u
        tmp2=self.prox1.work(tmp1)
        # update
        beta = exp((-iter/beta_switch)**3)*beta0+(1-exp((-iter/beta_switch)**3))*beta_max # unrelaxes as the
        unew = (beta*(2*tmp2-tmp1) + (1-beta)*tmp1 + u)/2
        return unew


