#                      ART_alternating_proj_in.m
#              written on October 6, 2011 by
#                        Russell Luke
#                 University of Goettingen
#
# DESCRIPTION:  parameter input file for main_ProxToolbox.m
#
##########################################################################
new_config={

## We start very general.
##
## What type of problem is being solved?  Classification 
## is according to the geometry:  'Affine', 'Phase', 
## 'Affine-sparsity',  'Custom'

'problem_family' : 'ART',

##==========================================
## Problem parameters
##==========================================
## What is the name of the data file?
'rescale' : 1,
'data_filename' : 'ART_data_processor',



## What type of object are we working with?
## Options are: 'phase', 'real', 'nonnegative', 'complex'
'object' : 'complex',

## What type of constraints do we have?
## Options are: 'support only', 'real and support', 'nonnegative and support',
##              'amplitude only', 'sparse real', 'sparse complex', and 'hybrid'
##              'convex'
'constraint' : 'convex',

## What type of measurements are we working with?
## Options are: 'single diffraction', 'diversity diffraction', 
##              'ptychography', 'complex', and 'diversity affine'
'experiment' : 'convex',


## What are the dimensions of the measurements?
# given in the ART_data_processor


##==========================================
##  Algorithm parameters
##==========================================
## Now set some algorithm parameters that the user should be 
## able to control (without too much damage)


# Point to appropriate projectors.

'Prox1' : 'P_parallel_hyperplane',# 'P_parallel_hyperplane', # 'P_sequential_hyperplane_odd', # projection onto support and nonnegativity constraint
'Prox2' : 'P_diag', #'P_Diag', #'P_sequential_hyperplane_even', # projection onto mask constraint

## Algorithm:
'algorithm' : 'AP', #'AP', 'Cimmino',  
'numruns' : 1, # the only time this parameter will
# be different than 1 is when we are
# benchmarking, that is, when algorithm performance statistics 
# are being generated for randomly generated problems/initial 
# values etc.  
## relaxaton parameters in RAAR, HPR and HAAR
'beta_0' : 0.75,  # starting relaxation prameter (only used with
# HAAR, HPR and RAAR)
'beta_max' : 0.75,   # maximum relaxation prameter (only used with
# HAAR, RAAR, and HPR)
'beta_switch' : 13, # iteration at which beta moves from beta_0 -> beta_max



## maximum number of iterations and tolerances
'MAXIT' : 20,
'TOL' : -1e-6,


## parameter for the data regularization
## need to discuss how/whether the user should
## put in information about the noise
'data_ball' : 1e-15,
# the above is the percentage of the gap
# between the measured data and the
# initial guess satisfying the
# qualitative constraints.  For a number
# very close to one, the gap is not expected
# to improve much.  For a number closer to 0
# the gap is expected to improve a lot.
# Ultimately the size of the gap depends
# on the inconsistency of the measurement model
# with the qualitative constraints.

##==========================================
## parameters for plotting and diagnostics
##==========================================
'diagnostic' : True, # to stop the diagnostics, just comment this field out.
'verbose' : 1, # options are 0 or 1
'graphics' : 1, # whether or not to display figures, options are 0 or 1.
                   # default is 1.
'anim' : 1,  # whether or not to disaply ``real time" reconstructions
                # options are 0=no, 1=yes, 2=make a movie
                # default is 1. Animation does not currently work for reordered arrays 
                
# below is just a hint of how one might make a movie of the reconstruction
# 'writerObj = VideoWriter('out.avi'), # Name it.
# 'writerObj.FrameRate = 25, # How many frames per second.
'animation' : 'ART_animation',
'graphics_display' : 'ART_graphics', # unless specified, a default 
                            # plotting subroutine will generate 
                            # the graphics.  Otherwise, the user
                            # can write their own plotting subroutine

##======================================================================
##  Technical/software specific parameters
##======================================================================
## Given the parameter values above, the following technical/algorithmic
## parameters are automatically set.  The user does not need to know 
## about these details, and so probably these parameters should be set in 
## a module one level below this one.  

}
