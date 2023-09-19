#                      Phase_graphics.m
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
#              method = character string for the algorithm used.
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
# USAGE: Phase_graphics(config,output)
#
#############################################################

from matplotlib.pyplot import subplots, show
import numpy as np


def Phase_graphics(config, output):
              
    algortihm=config['algorithm']
    beta0 = config['beta_0']
    beta_max = config['beta_max']
    u_0 = config['u_0']
    if output['u1'].ndim == 2:
        u = output['u1']
        u2 = output['u2']
    else:
        u = output['u1'][:,:,0]
        u2 = output['u2'][:,:,0]
    iter = output['iter']
    change = output['change']

    time = output['time'] if 'time' in output else 999
    f, ((ax1, ax2), (ax3, ax4)) = subplots(2, 2)

    im=ax1.imshow(np.abs(u),cmap='gray')
    f.colorbar(im, ax=ax1)
    ax1.set_title('best approximation amplitude - physical constraint satisfied')

    im=ax2.imshow(np.real(u),cmap='gray')
    f.colorbar(im, ax=ax2)
    ax2.set_title('best approximation phase - physical constraint satisfied')

    im=ax3.imshow(np.abs(u2),cmap='gray')
    f.colorbar(im, ax=ax3)
    ax3.set_title('best approximation amplitude - Fourier constraint satisfied')

    im=ax4.imshow(np.real(u2),cmap='gray')
    f.colorbar(im, ax=ax4)
    ax4.set_title('best approximation amplitude - Fourier constraint satisfied')

    g, ((bx1, bx2), (bx3, bx4)) = subplots(2, 2)
    im=bx1.imshow(np.abs(u),cmap='gray')
    f.colorbar(im, ax=bx1)
    bx1.set_title('best approximation amplitude - physical constraint satisfied')
    im = bx2.imshow(np.real(u), cmap='gray')
    f.colorbar(im, ax=bx2)
    bx2.set_title('best approximation phase - physical constraint satisfied')
    bx3.semilogy(change)
    bx3.set_xlabel('iteration')
    bx3.set_ylabel('$||x^{2k+2}-x^{2k}||$')
    if 'diagnostic' in config:
        bx4.semilogy(output['gap'])
        bx4.set_xlabel('iteration')
        bx4.set_ylabel('$||x^{2k+1}-x^{2k}||$')


    show()
