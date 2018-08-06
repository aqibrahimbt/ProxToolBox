#                      ART_graphics.m
#                  written on 29. Mai, 2012 by
#                        Russell Luke
#                  Universität Göttingen
#
# DESCRIPTION:  Script driver for viewing results from projection
#               algorithms on various toy problems
#
# INPUT:  
#              method_input/output = data structures
#
# OUTPUT:       graphics
# USAGE: ART_graphics(method_input,method_output)
#
#############################################################
from matplotlib.pyplot import subplots, show
import numpy as np

def ART_graphics(config, output):

    algorithm=config['algorithm']
    beta0 = config['beta_0']
    beta_max = config['beta_max']
    u_0 = config['u_0']
    N= int(np.sqrt(config['Ny']))
    u = np.reshape(output['u1'],(N, N))
    u2 = np.reshape(output['u2'],(N,N))
    iter = output['iter']
    change = output['change']

    if 'time' in output:
        time = output['time']
    else:
        time=999



    f, ((ax1, ax2), (ax3, ax4)) = subplots(2, 2)

    im=ax1.imshow(u,cmap='gray')
    f.colorbar(im, ax=ax1)
    ax1.set_title('best approximation- physical dpmain')

    im=ax2.imshow(u2,cmap='gray')
    f.colorbar(im, ax=ax2)
    ax2.set_title('best approximation - data constraint')

    ax3.plot(change)
    ax3.set_yscale('log')
    ax3.set_xlabel('iteration' + ', time = ' + str(time) + 's')
    ax3.set_ylabel('log of iterate difference')
    ax3.set_title('Algorithm: '+ algorithm + ', ' + r'$\beta =$' + str(beta0) + ' to ' + str(beta_max))

    if 'diagnostic' in config:
        gap = output['gap']
        ax4.plot(gap)
        #ax4.set_yscale('log')
        ax4.set_xlabel('iteration' + ', time = ' + str(time) + 's')
        ax4.set_ylabel('log of gap distance')
        ax4.set_title('Algorithm: '+ algorithm + ', ' + r'$\beta =$' + str(beta0) + ' to ' + str(beta_max))
        
    
    show()
    
    '''    
    im=ax3.imshow(np.abs(u2),cmap='gray')
    f.colorbar(im, ax=ax3)
    ax3.set_title('best approximation amplitude - Fourier constraint satisfied')

    im=ax4.imshow(np.real(u2),cmap='gray')
    f.colorbar(im, ax=ax4)
    ax4.set_title('best approximation amplitude - Fourier constraint satisfied')
    show()


    figure(900);

        subplot(2,2,1); imagesc((u)); colormap gray; axis equal tight; colorbar; title('best approximation - physical domain'); drawnow; 
        subplot(2,2,2); imagesc((u2)); colormap gray; axis equal tight; colorbar; title('best approximation - data constraint'); drawnow; #caxis([4.85,5.35]);
        label = [ 'iteration', ', time = ',num2str(time), 's'];
        subplot(2,2,3);   semilogy(change),xlabel(label),ylabel(['log of iterate difference'])
        label = ['Algorithm: ',method, ', \beta=',num2str(beta0),' to ',num2str(beta_max)];
        title(label)
        if(any(strcmp('diagnostic', fieldnames(method_input))))        
            gap  = method_output.stats.gap;    
            label = [ 'iteration', ', time = ',num2str(time), 's'];
            subplot(2,2,4);   semilogy(gap),xlabel(label),ylabel(['log of gap distance'])
            label = ['Algorithm: ',method, ', \beta=',num2str(beta0),' to ',num2str(beta_max)];
            title(label)
        end


    success = 1;
    return
    '''
