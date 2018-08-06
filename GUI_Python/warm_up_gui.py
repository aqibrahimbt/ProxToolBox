# -*- coding: utf-8 -*-
#!/usr/bin/python
import sys, os
from PyQt4 import QtGui, QtCore, uic
import json

class MainWindow(QtGui.QMainWindow):
	def __init__(self, parent=None):
		super(MainWindow, self).__init__(parent)
		uic.loadUi("warm_up_gui.ui", self)
		
		# Connection to button
		self.connect(self.okButton,QtCore.SIGNAL("clicked()"), self.ok)

	def ok(self):
		# get values from gui and store them into a json file
		if str(self.warmup.currentText()) == 'true':
			config = {
			'warmup':str(self.warmup.currentText()),
			'warmup_alg':str(self.warmup_alg.currentText()),
			'warmup_maxiter':int(self.warmup_maxiter.toPlainText()),};
		
		else:
			config = {'warmup':str(self.warmup.currentText()),}

		with open('data_warmup.json','wb') as outfile:
			json.dump(config, outfile)
				
		self.close()
		
		
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = MainWindow()
	window.show()
	app.exec_()
