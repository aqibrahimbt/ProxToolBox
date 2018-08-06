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
if config['switch_probemask'] == 'true':
	config['switch_probemask'] = True
else:
	config['switch_probemask'] = False
config['probe_mask_gamma'] = float(config['probe_mask_gamma'])
if config['switch_object_support_constraint'] == 'true':
	config['switch_object_support_constraint'] = True
else:
	config['switch_object_support_constraint'] = False
if config['warmup'] == 'true':
	config['warmup'] = True
	if config['warmup_alg'] == 'PALM':
		config['warmup_alg'] = PALM
else:
	config['warmup'] = False
if config['algorithm'] == 'PALM':
	config['algorithm'] = PALM
config['tol'] = float(config['tol'])
if config['ignore_error'] == 'true':
	config['ignore_error'] = True
else:
	config['ignore_error'] = False
if config['blocking_switch'] == 'true':
	config['blocking_switch'] = True
else:
	config['blocking_switch'] = False
if config['bs_mask'] == 'none':
	config['bs_mask'] = None
if config['fmask'] == 'none':
	config['fmask'] = None

# set up a ptychography class element and start algorithm
ptychography = Ptychography(config)
ptychography.solve()
raw_input('Press Enter to exit!')
