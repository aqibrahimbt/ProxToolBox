import sys, os
from PyQt4 import QtGui, QtCore, uic

class MainWindow(QtGui.QMainWindow):
	def __init__(self, parent=None):
		super(MainWindow, self).__init__(parent)
		uic.loadUi("start_gui.ui", self)

		self.connect(self.okButton,QtCore.SIGNAL("clicked()"), self.ok)

	def ok(self):
		self.close()
		str = self.start_gui.currentText()
		if str == "Ptychography":
			os.system("main_ptychography.py")
		elif str == "Sudoku":
			os.system("main_sudoku.py")
		
		
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = MainWindow()
	window.show()
	app.exec_()
	
'''
python C:\Python34\Lib\site-packages\PyQt4\uic\pyuic.py problem_family_gui.ui -o problem_family_gui.py -x
'''