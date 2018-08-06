
.. currentmodule:: proxtoolbox.Problems

Get started
-----------

The probably easiest, descriptive Problems instance implemented by the ProxToolbox is :class:`Sudoku`. It solves the well-known Sudoku puzzles.

.. contents::
   :local:


Run default Sudoku
++++++++++++++++++

The user needs to type the following lines into a terminal to import the Sudoku problem, initialize an instance with the default configuration and to start the solving process::

   from proxtoolbox.Problems import Sudoku
   prob_sudoku = Sudoku()
   prob_sudoku.solve()

Now, the user can access the solution for the specified Sudoku puzzle by::

   prob_sudoku.solution

To visualize the solution and to get some benchmarking results, the Problems class provides a separate method::

   prob_sudoku.show()


Adapt the configuration
+++++++++++++++++++++++

As already indicated, the Sudoku class provides a dictionary as default configuration, i.e.::

     config = {
        # This is the algorithm we use. RAAR and HPR will work.
        'algorithm':'RAAR',
        # RAAR requires 2 ProxOperators
        'proxoperators':('P_diag','P_parallel'),
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
        'normM':81,
        # Just a random Sudoku. Not too easy, but no challenge for
        # the mighty ProxToolbox!
        'sudoku':((2,0,0,0,0,1,0,3,0),
                  (4,0,0,0,8,6,1,0,0),
                  (0,0,0,0,0,0,0,0,0),
                  (0,0,0,0,1,0,0,0,0),
                  (0,0,0,0,0,0,9,0,0),
                  (0,0,5,0,0,3,0,0,7),
                  (0,0,0,0,0,0,0,0,0),
                  (1,0,0,0,0,7,4,9,0),
                  (0,2,4,1,0,0,0,0,0))
     }

Any Sudoku instance needs to contain a dictionary with the mappings listed above. Each parameter is assigned to a basic type, like a number, a boolean, a string or a collection of one of these types.

To initialize a Sudoku instance with a specific configuration, the user has to create a dictionary which contains mappings corresponding to the list above. This dictionary does not need to redefine every mapping but only those being different to the default values.

Assume the new configuration called "new_config" is stored in file "sudoku_in.py". Then the following commands initialize a Sudoku instance with this configuration::

   # from proxtoolbox.Problems import Sudoku
   from sudoku_in import new_config
   prob_sudoku = Sudoku(new_config)

Note that it is recommended to create a new Sudoku instance for every configuration.

