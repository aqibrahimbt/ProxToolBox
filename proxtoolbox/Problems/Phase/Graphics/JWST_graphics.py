#                      JWST_graphics.m
#                  written on May 23, 2012 by 
#                         Russell Luke
#       Inst. Fuer Numerische und Angewandte Mathematik
#                    Universitaet Gottingen
#
#
# DESCRIPTION:  Script driver for viewing results from projection
#               algorithms on various toy problems
#
# INPUT:  
#              algorithm = character string for the algorithm used.
#         true_object = the original image
#                 u_0 = the initial guess
#                   u = the algorithm "fixed point"
#              change = the norm square of the change in the
#                              iterates
#              error  = squared set distance error at each
#                              iteration
#              noneg = the norm square of the nonnegativity/support
#                              constraint at each iteration
#              gap  = the norm square of the gap distance, that is
#                     the distance between the projections of the
#                     iterates to the sets
#
# OUTPUT:       graphics
# USAGE: JWST_graphics(algorithm_input,algorithm_output)
#
#############################################################

from matplotlib.pyplot import subplots, show
from numpy import real,angle

def JWST_graphics(config, output):
              
    algorithm=config['algorithm'];

    if 'beta_0' in config:
      beta0 = config['beta_0']
      beta_max = config['beta_max']
    elif config['algorithm'] =='ADMMPlus':
        beta0=config['stepsize']
        beta_max=beta0
    else:
        beta0=1
        beta_max=1

    u_0 = config['u_0']
    u = output['u1']
    u2 = output['u2']
    iter = output['iter']
    change = output['change']
    gap  = output['gap']

    f, ((ax1, ax2), (ax3, ax4)) = subplots(2, 2)
    im=ax1.imshow(abs(u),cmap='gray')
    f.colorbar(im, ax=ax1)
    ax1.set_title('best approximation amplitude - physical constraint satisfied')
    im = ax2.imshow(real(u), cmap='gray')
    f.colorbar(im, ax=ax2)
    ax2.set_title('best approximation phase - physical constraint satisfied')
    ax3.semilogy(change)
    ax3.set_xlabel('iteration')
    ax3.set_ylabel('$||x^{2k+2}-x^{2k}||$')
    ax4.semilogy(gap)
    ax4.set_xlabel('iteration')
    ax4.set_ylabel('$||x^{2k+1}-x^{2k}||$')

    g, ((bx1, bx2), (bx3, bx4)) = subplots(2, 2)
    im = bx1.imshow(abs(u),cmap='gray')
    f.colorbar(im, ax=bx1)
    bx1.set_title('best approximation amplitude - physical constraint satisfied')
    im = bx2.imshow(real(u),cmap='gray')
    f.colorbar(im, ax=bx2)
    bx2.set_title('best approximation phase - physical constraint satisfied')
    im = bx3.imshow(abs(u2),cmap='gray')
    f.colorbar(im, ax=bx3)
    bx3.set_title('best approximation amplitude - Fourier constraint satisfied')
    im = bx4.imshow(real(u2),cmap='gray')
    f.colorbar(im, ax=bx4)
    bx4.set_title('best approximation phase - Fourier constraint satisfied')

    show()
