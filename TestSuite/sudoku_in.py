# -*- coding: utf-8 -*-

new_config = {
    # This is the algorithm we use. RAAR and HPR will work.
    'algorithm':'RAAR',
    # RAAR requires 2 ProxOperators
    'proj1':'P_diag',
    'proj2':'P_parallel',
    # P_parallel requires a sequence of projectors
    'projectors':('ProjRow','ProjColumn','ProjSquare','ProjGiven'),
    # Relaxation parameters for RAAR/HPR
    'beta0':1,
    'beta_max':1,
    'beta_switch':1,
    # Any algorithm requires these
    'maxiter':2000,
    'tol':1e-9,
    # Dimension parameters
    # which are the same for every standard Sudoku
    'Nx':9,
    'Ny':9,
    'Nz':9,
    'dim':4,
    'norm_data':81,
    # Just a random Sudoku. Not too easy, but no challenge for
    # the mighty ProxToolbox!
    'sudoku':((0,4,5,7,0,0,0,6,1),
              (0,0,0,8,0,5,4,0,0),
              (0,0,7,6,4,0,0,0,2),
              (5,0,0,9,6,4,0,0,0),
              (6,0,8,5,0,1,0,3,0),
              (0,7,0,0,0,0,0,9,0),
              (3,0,2,4,0,0,1,0,0),
              (9,0,0,1,0,0,0,0,6),
              (0,0,0,0,0,3,8,0,0))
}