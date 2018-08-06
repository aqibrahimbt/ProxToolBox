import os
import json
import yaml
import numpy
from sudoku import Sudoku, ProjRow, ProjColumn, ProjSquare, ProjGiven
from algorithms import *
from proxoperators import *
from problem import Problem

with open('data_in.json') as data_file:
	config = yaml.safe_load(data_file)
if os.path.isfile('edit_sudoku.json'):
	# use your own sudoku from the gui or the defined sudoku below
	with open('edit_sudoku.json') as data_file:
		s = yaml.safe_load(data_file)
		
	config.update({'given_sudoku':numpy.array(((s['a_1'],s['a_2'],s['a_3'],s['b_1'],s['b_2'],s['b_3'],s['c_1'],s['c_2'],s['c_3']),
											(s['a_4'],s['a_5'],s['a_6'],s['b_4'],s['b_5'],s['b_6'],s['c_4'],s['c_5'],s['c_6']),
											(s['a_7'],s['a_8'],s['a_9'],s['b_7'],s['b_8'],s['b_9'],s['c_7'],s['c_8'],s['c_9']),
											(s['d_1'],s['d_2'],s['d_3'],s['e_1'],s['e_2'],s['e_3'],s['f_1'],s['f_2'],s['f_3']),
											(s['d_4'],s['d_5'],s['d_6'],s['e_4'],s['e_5'],s['e_6'],s['f_4'],s['f_5'],s['f_6']),
											(s['d_7'],s['d_8'],s['d_9'],s['e_7'],s['e_8'],s['e_9'],s['f_7'],s['f_8'],s['f_9']),
											(s['g_1'],s['g_2'],s['g_3'],s['h_1'],s['h_2'],s['h_3'],s['i_1'],s['i_2'],s['i_3']),
											(s['g_4'],s['g_5'],s['g_6'],s['h_4'],s['h_5'],s['h_6'],s['i_4'],s['i_5'],s['i_6']),
											(s['g_7'],s['g_8'],s['g_9'],s['h_7'],s['h_8'],s['h_9'],s['i_7'],s['i_8'],s['i_9'])),dtype=numpy.float32),
											})

else:
	config.update({'given_sudoku':numpy.array(((2,0,0,0,0,1,0,3,0),
							(4,0,0,0,8,6,1,0,0),
							(0,0,0,0,0,0,0,0,0),
							(0,0,0,0,1,0,0,0,0),
							(0,0,0,0,0,0,9,0,0),
							(0,0,5,0,0,3,0,0,7),
							(0,0,0,0,0,0,0,0,0),
							(1,0,0,0,0,7,4,9,0),
							(0,2,4,1,0,0,0,0,0)),dtype=numpy.float32),
							})
if config['algorithm'] == 'RAAR':
	config['algorithm'] = RAAR
elif config['algorithm'] == 'HPR':
	config['algorithm'] = HPR
elif config['algorithm'] == 'AP':
	config['algorithm'] = AP
if config['proj1'] == 'P_parallel':
	config['proj1'] = P_parallel
if config['proj2'] == 'P_diag':
	config['proj2'] = P_diag
config['projectors'] = (ProjRow,ProjColumn,ProjSquare,ProjGiven)
config['tol'] = float(config['tol'])

# set up a ptychography class element and start algorithm
sudo = Sudoku(config)
sudo.solve()
raw_input('Press Enter to exit!')