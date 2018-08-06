# -*- coding: utf-8 -*-
#!/usr/bin/python
import sys, os
import numpy
import json
from PyQt4 import QtGui, QtCore, uic

class MainWindow(QtGui.QMainWindow):
	def __init__(self, parent=None):
		super(MainWindow, self).__init__(parent)
		uic.loadUi("main_ptychography.ui", self)
		
		# Connection to button
		self.connect(self.startButton,QtCore.SIGNAL("clicked()"), self.start)
		self.connect(self.warm_upButton,QtCore.SIGNAL("clicked()"), self.warm_up)
		self.connect(self.blockingButton,QtCore.SIGNAL("clicked()"), self.blocking)
		self.connect(self.input_dataButton,QtCore.SIGNAL("clicked()"), self.input_data)
	
	# start further guis for additional values
	def warm_up(self):
		os.system("./warm_up_gui.py")
		
	def blocking(self):
		os.system("./blocking_gui.py")
		
	def input_data(self):
		os.system("./input_data.py")
	
	# get values from gui and start algorithm
	def start(self):
		config = {
		'algorithm':str(self.algorithm.currentText()),
		'maxiter':int(self.maxiter.toPlainText()),
		'ptychography_prox':str(self.ptychography_prox.currentText()),
		'probe_guess_type':str(self.probe_guess_type.currentText()),
		'object_guess_type':str(self.object_guess_type.currentText()),
		'tol':float(self.tol.toPlainText()),
		'noise':str(self.noise.currentText()),
		'switch_probemask':str(self.switch_probemask.currentText()),
		'probe_mask_gamma':float(self.probe_mask_gamma.toPlainText()),
		'rmsfraction':float(self.rmsfraction.toPlainText()),
		'scan_type':str(self.scan_type.currentText()),
		'switch_object_support_constraint':str(self.switch_object_support_constraint.currentText()),
		'ignore_error':str(self.ignore_error.currentText()),
		'bs_mask':str(self.bs_mask.currentText()),
		'bs_factor':int(self.bs_factor.toPlainText()),
		'fmask':str(self.fmask.currentText()),
		'overrelax':int(self.overrelax.toPlainText()),}
		
		if os.path.isfile('data_input.json'):
			with open('data_input.json') as data_file:
				data_input = json.load(data_file)
				config.update(data_input)
			os.remove('data_input.json')
			
		if os.path.isfile('data_warmup.json'):
			with open('data_warmup.json') as data_file:
				data_warmup = json.load(data_file)
				config.update(data_warmup)
			os.remove('data_warmup.json')
		else:
			config.update({'warmup':'false',})
		
		if os.path.isfile('data_blocking.json'):
			with open('data_blocking.json') as data_file:
				data_blocking = json.load(data_file)
				config.update(data_blocking)
			os.remove('data_blocking.json')
		else:
			config.update({'blocking_switch':'false',})
		
		with open('data_in.json','w') as outfile:
			json.dump(config, outfile)
		
		
		if str(self.algorithm.currentText()) == 'RAAR' or str(self.algorithm.currentText()) == 'HPR':
			os.system("relaxtion_parameters_in_raar_hpr_haar.py")
		
		self.close()
		
		os.system("./start_ptychography.py")
		
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = MainWindow()
	window.show()
	app.exec_()
