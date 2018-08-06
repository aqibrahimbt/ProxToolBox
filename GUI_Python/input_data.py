# -*- coding: utf-8 -*-
#!/usr/bin/python
import sys, os
from PyQt4 import QtGui, QtCore, uic
import json

class MainWindow(QtGui.QMainWindow):
	def __init__(self, parent=None):
		super(MainWindow, self).__init__(parent)
		uic.loadUi("input_data.ui", self)

		# Connection to button
		self.connect(self.okButton,QtCore.SIGNAL("clicked()"), self.ok)

	def ok(self):
		# get values from gui and store them into a json file
		config = {
		'sim_data_type':str(self.sim_data_type.toPlainText()),
		'Nx':int(self.Nx.toPlainText()),
		'Ny':int(self.Ny.toPlainText()),
		'Nz':int(self.Nz.toPlainText()),
		'scan_stepsize':float(self.scan_stepsize.toPlainText()),
		'nx':int(self.nx.toPlainText()),
		'ny':int(self.ny.toPlainText()),
		'sample_area_center':str(self.sample_area_center.toPlainText()),};
		
		with open('data_input.json','wb') as outfile:
			json.dump(config, outfile)
				
		self.close()
		
		
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = MainWindow()
	window.show()
	app.exec_()
