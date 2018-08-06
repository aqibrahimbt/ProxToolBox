# resizes arrays.  Takes an array A that is M by N and resizes it to Mnew by
# Nnew.
# usage: Anew = Resize(A, Mnew,Nnew)

from numpy.random import rand
from numpy import exp, sqrt, log, tan, pi, floor, zeros, ones
from scipy.special import gammaln
import numpy as np

def Resize(A, Mnew,Nnew):

    M = A.shape[0];
    N = A.shape[1];
    Anew = zeros((Mnew,Nnew));
    if Nnew <= N:
      Nstep = int(N/Nnew);
      Mstep = int(M/Mnew);

      countN = -1;
      for j in  range(0,N,Nstep):
        countN = countN+1;
        countM = -1;
        for i in range(0,M,Mstep):
          countM = countM+1;
          Anew[countM,countN] = A[i,j];
 
    else:
      Nstep = int(Nnew/N);
      Mstep = int(Mnew/M);
      temp=ones((Mstep,Nstep));

      countN = -1;
      for j in range(N):
        countN = countN+1;
        countM = -1;
        for i in range(M):
          countM = countM+1;
          Anew[i*Mstep:(i+1)*Mstep,j*Nstep:(j+1)*Nstep] = A[countM,countN]*temp;
    
    return Anew
 

