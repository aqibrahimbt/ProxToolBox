# -*- coding: utf-8 -*-
#!/usr/bin/python
import os
import json
import yaml
import numpy
from ptychography import Ptychography, PALM
from algorithms import *
from proxoperators import *
from problem import Problem

# assign properties from json file
# be careful about the type of each property
with open('data_in.json') as data_file:
	config = yaml.safe_load(data_file)
config['scan_stepsize'] = float(config['scan_stepsize'])
if config['noise'] == 'none':
	config['noise'] = None
config['switch_probemask'] = config['switch_probemask'] == 'true'
config['probe_mask_gamma'] = float(config['probe_mask_gamma'])
config['switch_object_support_constraint'] = (
	config['switch_object_support_constraint'] == 'true'
)
if config['warmup'] == 'true':
	config['warmup'] = True
	if config['warmup_alg'] == 'PALM':
		config['warmup_alg'] = PALM
else:
	config['warmup'] = False
if config['algorithm'] == 'PALM':
	config['algorithm'] = PALM
config['tol'] = float(config['tol'])
config['ignore_error'] = config['ignore_error'] == 'true'
config['blocking_switch'] = config['blocking_switch'] == 'true'
if config['bs_mask'] == 'none':
	config['bs_mask'] = None
if config['fmask'] == 'none':
	config['fmask'] = None

# set up a ptychography class element and start algorithm
ptychography = Ptychography(config)
ptychography.solve()
raw_input('Press Enter to exit!')
