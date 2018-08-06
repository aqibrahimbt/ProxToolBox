import sys
sys.path.append('../proxtoolbox/Problems/CT')
sys.path.append('..')
from proxtoolbox.Problems import Sudoku

sudoku = Sudoku()
sudoku.solve()
sudoku.show()


