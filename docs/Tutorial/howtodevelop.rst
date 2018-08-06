
.. currentmodule:: proxtoolbox

Behind the scenes
-----------------

There is more behind the routines which the user can see at first glance.
Have a closer look at:

.. contents::
   :local:


Problem interface
+++++++++++++++++

The class :class:`Sudoku<Problems.Sudoku>` inherits from the interface :class:`Problem<Problems.problems.Problem>` which determines some basic methods Sudoku should have. Most of them still needs to be implemented, i.e.
   * __init__
   * _presolve
   * _solve
   * _postsolve
   * show
   * save


Dependencies Problem-Algorithm-ProxOperators
++++++++++++++++++++++++++++++++++++++++++++

The solving process of a Sudoku puzzle is based on an algorithm and several proximal operators - as intended by the structure of the ProxToolbox. They are all set in the configuration, remember::

     config = {
        # This is the algorithm we use. RAAR and HPR will work.
        'algorithm':'RAAR',
        # RAAR requires 2 ProxOperators
        'proxoperators':('P_diag','P_parallel'),
        # P_parallel requires a sequence of projectors
        'projectors':('ProjRow','ProjColumn','ProjSquare','ProjGiven'),
        [..]
     }

The algorithm :class:`RAAR<Algorithms.RAAR>` is contained in the Problems module and inherits from the corresponding interface :class:`Problem<Problems.problems.Problem>`.

RAAR requires two proximal operators. The chosen operators, :class:`P_diag<ProxOperators.P_diag>` and :class:`P_parallel<ProxOperators.P_parallel>`, belong to the ProxOperators module with interface :class:`ProxOperator<ProxOperators.proxoperators.ProxOperator>`.

Furthermore, P_parallel processes a sequence of projections which also inherit from the ProxOperators module's interface. Since these four operators, 'ProjRow','ProjColumn','ProjSquare' and 'ProjGiven', work explicitly on the Sudoku data, they are implemented in the sudoku-file and do not belong to the ProxOperators module.

