import sys, os
import numpy
import json
from PyQt4 import QtGui, QtCore, uic

class MainWindow(QtGui.QMainWindow):
	def __init__(self, parent=None):
		super(MainWindow, self).__init__(parent)
		uic.loadUi("main_sudoku.ui", self)
		
		# Connection to button
		self.connect(self.startButton,QtCore.SIGNAL("clicked()"), self.start)
		self.connect(self.edit_sudoku,QtCore.SIGNAL("clicked()"), self.sudoku)
	
	# start gui to define own sudoku
	def sudoku(self):
		os.system("edit_sudoku.py")
	
	# get values from gui and start algorithm
	def start(self):
		config = {
		'algorithm':str(self.algorithm.currentText()),
		'proj1':str(self.proj1.toPlainText()),
		'proj2':str(self.proj2.toPlainText()),
		'maxiter':int(self.maxiter.toPlainText()),
		'tol':float(self.tol.toPlainText()),};
		
		with open('data_in.json','w') as outfile:
			json.dump(config, outfile)
		
		if str(self.algorithm.currentText()) == 'RAAR' or str(self.algorithm.currentText()) == 'HPR':
			os.system("relaxtion_parameters_in_raar_hpr_haar.py")
		
		self.close()
		
		os.system("start_sudoku.py")

		
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = MainWindow()
	window.show()
	app.exec_()