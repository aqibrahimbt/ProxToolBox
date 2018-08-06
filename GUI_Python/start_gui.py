# -*- coding: utf-8 -*-
#!/usr/bin/python
import sys, os
from PyQt4 import QtGui, QtCore, uic

class MainWindow(QtGui.QMainWindow):
	def __init__(self, parent=None):
		super(MainWindow, self).__init__(parent)
		uic.loadUi("start_gui.ui", self)
		
		# Connection to button
		self.connect(self.okButton,QtCore.SIGNAL("clicked()"), self.ok)

	def ok(self):
		# start main ptychography or sudoku gui
		self.close()
		str = self.problem_family.currentText()
		if str == "Ptychography":
			os.system("./main_ptychography.py")
		elif str == "Sudoku":
			os.system("./main_sudoku.py")
		
		
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = MainWindow()
	window.show()
	app.exec_()
