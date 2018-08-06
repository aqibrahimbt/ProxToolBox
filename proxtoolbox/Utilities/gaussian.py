# DESCRIPTION: Script for forming a 1-2D Gaussian of height 1 on a
#              Cartesian grid. 
# INPUT:  points       := number of points for gaussian
#                         representation (assumes square arrays)          
#         var          := a 1 or 2 dimensional array that gives the
#                         variance of the gaussian in the x and y 
#                         directions respectively.
#         center_index := a 1 or 2 dimensional array that gives the
#                         index coordinates of the center of the
#                         gaussian.
#                         Note that center index needs to be shifted by one compared to matlab
# OUTPUT: g, a (points)^n-dimensional gaussian of unit magnitude,
#         where n is the dimension and points is the discretization
#         resolution.
# USAGE:  g=Gaussian(points,var,center_index)
#
# Written by Russell Luke, July 17, 2001
# luke@math.uni-goettingen.de
#
from numpy import arange, pi, exp, zeros

def gaussian(points,var,center_index):

    dim = center_index.size
    L = 2
    x = arange(-L/2,L/2,L/points)
    if dim==1:
      centerx = x[center_index]
      g = exp(-((x-centerx)**2)*pi/var[1])
    elif dim==2:
        
      if( (var.size==1) and (dim==2)):
        var[1]=var[0]
        
      centerx=x[center_index[0]]
      centery=x[center_index[1]]
      
      g = zeros((points,points))
      
      for j in range(points):
        g[j,:]=exp(-(((x-centerx)**2)/var[0]+((x[j]-centery)**2)/var[1])*pi)
    
    return g