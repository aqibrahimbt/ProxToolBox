# -*- coding: utf-8 -*-
"""
Created on Thu Dec 18 09:52:17 2014

@author: stefan

Just a small benchmarking script
"""

from matplotlib import pyplot
from ptychography import Ptychography_NTT_01_26210

import numpy, matplotlib


# Sort out some functions useless for our benchmark
class BenchmarkPtychography(Ptychography_NTT_01_26210):
    
    def _postsolve(self):
        return;


def go():
    config = BenchmarkPtychography.default_config.copy();
    
    config['blocking_switch'] = True;
    config['maxiter'] = 250;
    config['block_maxiter'] = 250;
    
    # This might still need some thinking
    config['tol'] = 1e-5;
    
    blocking_sizes = (2,4,8,16,32,64);
    blocking_schemes = [('sequential','none'),('parallel','none'),
                        ('sequential','sequential'),('parallel','sequential')];
    algorithms = ['PALM','PALMregPhiPtwise','Rodenburg','Thibault'];
    
    for i in range(6):
        config['block_rows'] = blocking_sizes[i];
        config['block_cols'] = blocking_sizes[i];
        
        for j in range(4):
            config['between_blocks_scheme'] = blocking_schemes[j][0];
            config['within_blocks_scheme'] = blocking_schemes[j][1];
            
            for k in range(4):
                config['ptychography_prox'] = algorithms[k];
                
                problem = BenchmarkPtychography(config);
                print 'Experiment', k + 4*j + 16*i + 1;
                print 'block_rows/block_cols:', blocking_sizes[i];
                print 'between_blocks_scheme: %s' % blocking_schemes[j][0];
                print 'within_blocks_scheme: %s'  % blocking_schemes[j][1];
                print 'ptychography_prox : %s'    % algorithms[k];
                problem.solve();
                print '';
                
                pyplot.imsave('%d%d%d_amplitude.png' % (i,j,k), \
                    numpy.abs(problem.u_final.obj));
                pyplot.imsave('%d%d%d_phase.png' % (i,j,k), \
                    numpy.angle(problem.u_final.obj));
                pyplot.imsave('%d%d%d_probe.png' % (i,j,k), \
                    matplotlib.colors.hsv_to_rgb(problem._im2hsv(problem.u_final.probe,1.0)));
                
                numpy.savez_compressed('%d%d%d_data' % (i,j,k), \
                    phi=problem.u_final.phi, \
                    obj=problem.u_final.obj, \
                    probe=problem.u_final.probe, \
                    change=problem.change, \
                    rms_err_obj=problem.custom_errors[0,:], \
                    rms_err_probe=problem.custom_errors[1,:], \
                    phi_change=problem.custom_errors[2,:], \
                    rfactor=problem.custom_errors[3,:], \
                    objective_value=problem.custom_errors[4,:], \
                    block_iters=problem.algorithm.block_iters);
                