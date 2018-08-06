import sys, os
from PyQt4 import QtGui, QtCore, uic
import json

class MainWindow(QtGui.QMainWindow):
	def __init__(self, parent=None):
		super(MainWindow, self).__init__(parent)
		uic.loadUi("relaxtion_parameters_in_raar_hpr_haar.ui", self)

		# Connection to button
		self.connect(self.okButton,QtCore.SIGNAL("clicked()"), self.ok)

	def ok(self):
		# get values from gui and store them into a json file
		with open('data_in.json') as data_file:
			config = json.load(data_file)
				
		config.update({
		'beta0':float(self.beta_0.toPlainText()),
		'beta_max':float(self.beta_max.toPlainText()),
		'beta_switch':float(self.beta_switch.toPlainText()),});
		
		with open('data_in.json','wb') as outfile:
			json.dump(config, outfile)
				
		self.close()
		
		
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = MainWindow()
	window.show()
	app.exec_()