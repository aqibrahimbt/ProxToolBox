# Poisson random number generation with mean of xm
# April 6, 2005 Ning Lei
# see numerical recipe C++ version
# ref: y(i+1) = xm^i*exp(-xm)/factorial(i); see Leo Breiman Probability

from numpy.random import rand
from numpy import exp, sqrt, log, tan, pi, floor, zeros, ones
from scipy.special import gammaln
import numpy as np

def PoissonRan(xm):
  oldm = -1;
  if xm<12:
    if xm != oldm:
      oldm = xm;
      g = exp(-xm);

    em = -1;
    t = 1;

    em = em+1;
    t = t*rand(1);
    while t>g :
      em = em+1;
      t = t*rand(1);

  else:
    if xm != oldm:
      oldm = xm;
      sq = sqrt(2.0*xm);
      alxm = log(xm);
      g = xm*alxm-gammaln(xm+1);

    y = tan(pi*rand(1));
    em = sq*y+xm;

    while em < 0:
      y = tan(pi*rand(1));
      em = sq*y+xm;
  
    em = floor(em);
    t = 0.9*(1+y*y)*exp(em*alxm-gammaln(em+1)-g);

    while rand(1) > t:
      y = tan(pi*rand(1));
      em = sq*y+xm;

      while em < 0:
        y = tan(pi*rand(1));
        em = sq*y+xm;

      em = floor(em);
      t = 0.9*(1+y*y)*exp(em*alxm-gammaln(em+1)-g);
   
  return em;
