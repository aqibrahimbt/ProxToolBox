# -*- coding: utf-8 -*-
"""
Created on Thu Dec 18 09:52:17 2014

@author: stefan

Just a small benchmarking script
"""

import numpy, matplotlib
import json, time

from matplotlib import pyplot
from Problems import Ptychography



# Sort out some functions that are useless for our benchmark
class BenchmarkPtychography(Ptychography):
    
    def _postsolve(self):
        return;


def go():
    config = BenchmarkPtychography.default_config.copy();
    
    config['blocking_switch'] = True;
    config['maxiter'] = 300;
    config['between_blocks_scheme'] = 'parallel';
    config['within_blocks_scheme'] = 'none';
    #config['overrelax'] = 15;
    
    config['tol'] = 1e-5;
    
    algorithms = ['PALM','PALMregPhiPtwise','Rodenburg','Thibault'];
    blocks = [(2,2), (4,4)];   

    times = {};

    for b in blocks:    
        for a in algorithms:
            config['ptychography_prox'] = a;
            config['block_rows'] = b[0];
            config['block_cols'] = b[1];
            problem = BenchmarkPtychography(config);
            print('Algorithm {}, {}x{} blocks'.format(a, b[0], b[1]));
        
            t = time.clock();
            problem.solve();
            t = time.clock() - t;
            times[a + '_parallel'] = t;
        
            pyplot.imsave('%s_parallel_amplitude.png' % a, \
                numpy.abs(problem.u_final.obj));
            pyplot.imsave('%s_parallel_phase.png' % a, \
                numpy.angle(problem.u_final.obj));
            pyplot.imsave('%s_parallel_probe.png' % a, \
            matplotlib.colors.hsv_to_rgb(problem._im2hsv(problem.u_final.probe,1.0)));
        
            numpy.savez_compressed('%s_parallel_data' % a, \
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
    
        config['blocking_switch'] = False;
    
    for a in algorithms:
        config['ptychography_prox'] = a;
        problem = BenchmarkPtychography(config);
        print('Algorithm {}, no blocking'.format(a));
        
        t = time.clock();
        problem.solve();
        t = time.clock() - t;
        times[a + '_sequential'] = t;
        
        pyplot.imsave('%s_sequential_amplitude.png' % a, \
            numpy.abs(problem.u_final.obj));
        pyplot.imsave('%s_sequential_phase.png' % a, \
            numpy.angle(problem.u_final.obj));
        pyplot.imsave('%s_sequential_probe.png' % a, \
        matplotlib.colors.hsv_to_rgb(problem._im2hsv(problem.u_final.probe,1.0)));
        
        numpy.savez_compressed('%s_sequential_data' % a, \
            phi=problem.u_final.phi, \
            obj=problem.u_final.obj, \
            probe=problem.u_final.probe, \
            change=problem.change, \
            rms_err_obj=problem.custom_errors[0,:], \
            rms_err_probe=problem.custom_errors[1,:], \
            phi_change=problem.custom_errors[2,:], \
            rfactor=problem.custom_errors[3,:], \
            objective_value=problem.custom_errors[4,:]);
    
    json.dump(times, open('time.json','w'));
        
                
