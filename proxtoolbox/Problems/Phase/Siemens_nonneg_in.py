#                   Siemens_amplitude_in
#              written on October 6, 2011 by
#                        Russell Luke
#                 University of Goettingen
#
# DESCRIPTION:  parameter input file for main_dvr.m
#
##########################################################################
## This is the input file that the user sees/modifies.  It should be simple, 
## avoid jargon or acronyms, and should be a model for a menu-driven GUI

new_config = {

    ## We start very general.
    ## What type of problem is being solved?  Classification 
    ## is according to the geometry:  Affine, Cone, Convex, 
    ## Phase, Affine-sparsity, Nonlinear-sparsity, Sudoku

    'problem_family' : 'Phase',

    ##==========================================
    ## Problem parameters
    ##==========================================
    ## What is the name of the data file?
    'data_filename' : 'Siemens_processor',


    ## What type of object are we working with?
    ## Options are: 'phase', 'real', 'nonnegative', 'complex'
    'object' : 'nonnegative',

    ## What type of constraints do we have?
    ## Options are: 'support only', 'real and support', 'nonnegative and support',
    ##              'amplitude only', 'minimum amplitude', 'sparse real', 'sparse complex', and 'hybrid'
    'constraint' : 'nonnegative and support',

    ## What type of measurements are we working with?
    ## Options are: 'single diffraction', 'diversity diffraction', 
    ##              'ptychography', and 'complex'
    'experiment' : 'single diffraction',

    ## Next we move to things that most of our users will know 
    ## better than we will.  Some of these may be overwritten in the 
    ## data processor file which the user will most likely write. 
    ## Are the measurements in the far field or near field?
    ## Options are: 'far field' or 'near field', 
    'distance' : 'far field',

    ## What are the dimensions of the measurements?
    'Nx' : 512,
    'Ny' : 512,
    'Nz':1, # only 2 dimensional arrays

    ## What are the noise characteristics (Poisson or Gaussian or none)?
    'noise':'Gaussian',

    ##==========================================
    ##  Algorithm parameters
    ##==========================================
    ## Now set some algorithm parameters that the user should be 
    ## able to control (without too much damage)

    ## Algorithm:
    'algorithm' : 'RAAR',#'RAAR', 'AP',  
    'numruns':1, # the only time this parameter will
    # be different than 1 is when we are
    # benchmarking...not something a normal user
    # would be doing.


    ## The following are parameters specific to RAAR, HPR, and HAAR that the 
    ## user should be able to set/modify.  Surely
    ## there will be other algorithm specific parameters that a user might 
    ## want to play with.  Don't know how best 
    ## to do this.  Thinking of a GUI interface, we could hard code all the 
    ##  parameters the user might encounter and have the menu options change
    ## depending on the value of the method field. 
    ## do different things depending on the chosen algorithm:
        ## maximum number of iterations and tolerances
        'MAXIT' : 500,
        'TOL' : 1e-5,

        ## relaxaton parameters in RAAR, HPR and HAAR
        'beta_0' : 1.0,                # starting relaxation prameter (only used with
        # HAAR, HPR and RAAR)
        'beta_max' :0.5,             # maximum relaxation prameter (only used with
        # HAAR, RAAR, and HPR)
        'beta_switch' : 50,           # iteration at which beta moves from beta_0 -> beta_max

        ## parameter for the data regularization 
        ## need to discuss how/whether the user should
        ## put in information about the noise
        'data_ball' : 0.1,  
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
    'diagnostic' : True, #to turn this off, just comment out
    'verbose' : 1, # options are 0 or 1
    'graphics' : 1, # whether or not to display figures, options are 0 or 1.
                       # default is 1.
    'anim' : 1,  # whether or not to disaply ``real time" reconstructions
                    # options are 0=no, 1=yes, 2=make a movie
                    # default is 1.
    'graphics_display' : 'Phase_graphics', # unless specified, a default 
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
