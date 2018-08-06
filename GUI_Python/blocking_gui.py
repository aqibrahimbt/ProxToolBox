# -*- coding: utf-8 -*-
#!/usr/bin/python
import sys, os
from PyQt4 import QtGui, QtCore, uic
import json

class MainWindow(QtGui.QMainWindow):
	def __init__(self, parent=None):
		super(MainWindow, self).__init__(parent)
		uic.loadUi("blocking_strategies_gui.ui", self)

		# Connection to button
		self.connect(self.okButton,QtCore.SIGNAL("clicked()"), self.ok)

	def ok(self):
		# get values from gui and store them into a json file
		if str(self.blocking_switch.currentText()) == 'true':
			config = {
			'blocking_switch':str(self.blocking_switch.currentText()),
			'blocking_scheme':str(self.blocking_scheme.currentText()),
			'between_blocks_scheme':str(self.between_blocks_scheme.currentText()),
			'within_blocks_scheme':str(self.within_blocks_scheme.currentText()),
			'block_maxiter':int(self.block_maxiter.toPlainText()),
			'block_rows':int(self.block_rows.toPlainText()),
			'block_cols':int(self.block_cols.toPlainText()),};
		else:
			config = {'blocking_switch':str(self.blocking_switch.currentText()),}
		
		with open('data_blocking.json','wb') as outfile:
			json.dump(config, outfile)
				
		self.close()
		
		
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = MainWindow()
	window.show()
	app.exec_()
