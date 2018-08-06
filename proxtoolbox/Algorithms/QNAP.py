#                      QNAP.m
#             written on May 23, 2017 by
#                   Russell Luke
#   Inst. Fuer Numerische und Angewandte Mathematik
#                Universitaet Gottingen
#
# DESCRIPTION:  Averaged prox mappings with a Quasi-Newton acceleration
#
# INPUT:  method_input, a data structure
#
# OUTPUT: method_output, a data structure with
#               u = the algorithm fixed point
#               stats = [iter, change, gap]  where
#                     gap    = squared gap distance normalized by
#                              the magnitude constraint.
#                     change = the percentage change in the norm
#                              squared difference in the iterates
#                     iter   = the number of iterations the algorithm
#                              performed
#
# USAGE: method_output = QNAP(method_input)
#
# Nonstandard Matlab function calls:  samsara, self.prox1 and .Prox2
#    It is assumed that the Prox1 is a projection onto the diagonal of the
#    product space, P_Diag.  Must have samsara, a reverse communication nonlinear
#    optimization package installed
#
#

from numpy import zeros, isreal, zeros_like
import numpy as np
from .algorithms import Algorithm
import sys
sys.path.append('../../samsara/python')
from samsara import Samsara

class QNAP(Algorithm):

    def __init__(self,config):

        self.iter = 0
        self.prox1 = config['proxoperators'][0](config)
        self.prox2 = config['proxoperators'][1](config)
        self.norm_data = config['norm_data']
        self.Nx = config['Nx']; self.Ny = config['Ny']; self.Nz = config['Nz']
        self.product_space_dimension = config['product_space_dimension']

        if 'truth' in config:
            self.truth = config['truth']
            self.truth_dim = config['truth_dim']
            self.norm_truth = config['norm_truth']

        if 'diagnostic' in config:
            self.diagnostic = True
        else:
            self.diagnostic = False

        self.samsara = Samsara()

    def run(self, u, tol, maxiter):

        """
        Runs the algorithm for the specified input data
        """

        ##### PREPROCESSING

        iter = self.iter
        prox1 = self.prox1; prox2 = self.prox2

        if u.ndim < 3:
            p = 1
            q = 1
        elif u.ndim == 3:
            p = u.shape[2]
            q = 1
        else:
            p = u.shape[2]
            q = u.shape[3]

        change = zeros(maxiter+1,dtype=u.dtype)
        change[0] = 999

        ######################################
        #  set up diagnostic arrays
        ######################################
        if self.diagnostic:
            gap = change.copy()
            shadow_change = change.copy()
        if  hasattr(self, 'truth'):
            Relerrs = change.copy()

        norm_data = self.norm_data

        tmp1 = self.prox2.work(u)
        
        gradfnew_vec = zeros_like(u[:,0])

        while iter < maxiter and change[iter] >= tol:
            iter += 1
            self.iter =iter
            tmp_u = prox1.work(tmp1)
            # make call to SAMSARA for acceleration step
            # since for the product space Prox1 is the
            # projection onto the diagonal, all arrays in
            # tmp_u are the same, so we can just send one
            # array to samsara.  Have to reshape the complex
            # array as a real-valued vector.  Use the gap as the
            # objective function, and the negative gradient tmp_u- u
        if (p==1) and (q==1):
            if (np.all(isreal(u))):
                if iter > 4:
                    uold_vec=u[:,0]
                    unew_vec=tmp_u[:,0]
                    gradfold_vec = gradfnew_vec
                    gradfnew_vec = (u[:,0]- tmp_u[:,0])
                    unew_vec, uold_vec, gap[iter-1], gradfold_vec, change[iter] = \
                        self.samsara.run(uold_vec, unew_vec,
                        gap[iter-2]*self.Nx*self.Ny,
                        gap[iter-1]*self.Nx*self.Ny,
                        gradfold_vec, gradfnew_vec)
                    for j in range(self.product_space_dimension):
                        tmp_u[:,j]=unew_vec
                    gap[iter-1]=gap[iter-1]/(self.Nx*self.Ny)
                else:
                    unew_vec=tmp_u[:,0]
                    uold_vec=u[:,0]
                    gradfnew_vec = (u[:,0]- tmp_u[:,0])
            else:
                if iter>3:
                    uold_vec=reshape(concatenate(real(u[:,0]), imag(u[:,0])),
                        self.Ny*self.Nx*2, 0);
                    unew_vec=reshape(concatenate(real(tmp_u[:,0]), imag(tmp_u[:,0])),
                        self.Ny*self.Nx*2, 0)
                    gradfold_vec = gradfnew_vec
                    gradfnew_vec = uold_vec-unew_vec
                    unew_vec,uold_vec,gap[iter-0],gradfold_vec, change[iter]= \
                        feval('samsara', uold_vec, unew_vec,
                        gap(iter-2)*self.Nx*self.Ny,
                        gap(iter-1)*self.Nx*self.Ny,
                        gradfold_vec, gradfnew_vec)
                    # now reshape unew_vec
                    tmp_u_vec = unew_vec[0:self.Ny*self.Nx-1]+1j*unew_vec[self.Ny*self.Nx:self.Ny*self.Nx*2-1]
                    for j in range(self.product_space_dimension):
                        tmp_u[:,j]=tmp_u_vec
                    gap[iter-1]=gap[iter-1]/(self.Nx*self.Ny)
                else:
                    uold_vec=reshape([real(u[:,0]), imag(u[:,0])],
                        self.Ny*self.Nx*2, 1)
                    unew_vec=reshape([real(tmp_u[:,0]), imag(tmp_u[:,0])],
                        self.Ny*self.Nx*2, 1)
                    gradfnew_vec = uold_vec-unew_vec

        elif (p!=1) and (q==1):
            if (np.all(isreal(u))):
                tmp2 = u[:,:,0]
                uold_vec=reshape(tmp2,self.Nx*self.Ny,1)
                tmp2 = tmp_u[:,:,0]
                unew_vec = reshape(tmp2,self.Nx*self.Ny,1);
                if iter<=3:
                    gradfnew_vec = uold_vec-unew_vec
                else:
                    gradfold_vec = gradfnew_vec
                    gradfnew_vec = uold_vec-unew_vec
                    unew_vec,uold_vec,gap[iter-1],gradfold_vec, change[iter]= \
                        feval('samsara', uold_vec, unew_vec, 
                        gap(iter-2)*self.Nx*self.Ny,
                        gap(iter-1)*self.Nx*self.Ny,
                        gradfold_vec, gradfnew_vec);
                    # now reshape and replace u 
                    gap[iter-1]=gap[iter-1]/(self.Nx*self.Ny);

                tmp2=reshape(unew_vec,self.Ny,self.Nx);
                for j in range(self.product_space_dimension):
                    tmp_u[:,:,j]=tmp2
                end
            else:
                tmp2=concatenate(real(u[:,:,0]),imag(u[:,:,0]))
                uold_vec=reshape(tmp2,self.Nx*self.Ny*2,1)
                tmp2=concatenate(real(tmp_u[:,:,0]), imag(tmp_u[:,:,0]))
                unew_vec = reshape(tmp2,self.Nx*self.Ny*2,1)
                if iter<=3:
                    gradfnew_vec = uold_vec-unew_vec
                else:
                    gradfold_vec = gradfnew_vec
                    gradfnew_vec = uold_vec-unew_vec
                    unew_vec,uold_vec,gap[iter-1],gradfold_vec, change[iter]= \
                        feval('samsara', uold_vec, unew_vec,
                        gap(iter-2)*self.Nx*self.Ny,
                        gap(iter-1)*self.Nx*self.Ny,
                        gradfold_vec, gradfnew_vec)
                    # now reshape and replace u 
                    gap[iter-1]=gap[iter-1]/(self.Nx*self.Ny)
                    tmp2=reshape(unew_vec,self.Ny,self.Nx*2)
                    unew=tmp2[:,1:self.Nx] +1j*tmp2[:,self.Nx+1:2*self.Nx]
                    for j in range(self.product_space_dimension):
                        tmp_u[:,:,j]=unew
        else: # product space of 3D arrays
            print('Cannot handle 3D arrays on the product space yet')

