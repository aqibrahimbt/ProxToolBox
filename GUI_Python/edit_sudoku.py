import sys, os
from PyQt4 import QtGui, QtCore, uic
import json

class MainWindow(QtGui.QMainWindow):
	def __init__(self, parent=None):
		super(MainWindow, self).__init__(parent)
		uic.loadUi("edit_sudoku.ui", self)

		# Connection to button
		self.connect(self.okButton,QtCore.SIGNAL("clicked()"), self.ok)

	def ok(self):
		# get vaules from each cell and define the sudoku
		a_1 = int(self.a_1.toPlainText())
		a_2 = int(self.a_2.toPlainText())
		a_3 = int(self.a_3.toPlainText())
		a_4 = int(self.a_4.toPlainText())
		a_5 = int(self.a_5.toPlainText())
		a_6 = int(self.a_6.toPlainText())
		a_7 = int(self.a_7.toPlainText())
		a_8 = int(self.a_8.toPlainText())
		a_9 = int(self.a_9.toPlainText())
		
		b_1 = int(self.b_1.toPlainText())
		b_2 = int(self.b_2.toPlainText())
		b_3 = int(self.b_3.toPlainText())
		b_4 = int(self.b_4.toPlainText())
		b_5 = int(self.b_5.toPlainText())
		b_6 = int(self.b_6.toPlainText())
		b_7 = int(self.b_7.toPlainText())
		b_8 = int(self.b_8.toPlainText())
		b_9 = int(self.b_9.toPlainText())
		
		c_1 = int(self.c_1.toPlainText())
		c_2 = int(self.c_2.toPlainText())
		c_3 = int(self.c_3.toPlainText())
		c_4 = int(self.c_4.toPlainText())
		c_5 = int(self.c_5.toPlainText())
		c_6 = int(self.c_6.toPlainText())
		c_7 = int(self.c_7.toPlainText())
		c_8 = int(self.c_8.toPlainText())
		c_9 = int(self.c_9.toPlainText())
		
		d_1 = int(self.d_1.toPlainText())
		d_2 = int(self.d_2.toPlainText())
		d_3 = int(self.d_3.toPlainText())
		d_4 = int(self.d_4.toPlainText())
		d_5 = int(self.d_5.toPlainText())
		d_6 = int(self.d_6.toPlainText())
		d_7 = int(self.d_7.toPlainText())
		d_8 = int(self.d_8.toPlainText())
		d_9 = int(self.d_9.toPlainText())
		
		e_1 = int(self.e_1.toPlainText())
		e_2 = int(self.e_2.toPlainText())
		e_3 = int(self.e_3.toPlainText())
		e_4 = int(self.e_4.toPlainText())
		e_5 = int(self.e_5.toPlainText())
		e_6 = int(self.e_6.toPlainText())
		e_7 = int(self.e_7.toPlainText())
		e_8 = int(self.e_8.toPlainText())
		e_9 = int(self.e_9.toPlainText())
		
		f_1 = int(self.f_1.toPlainText())
		f_2 = int(self.f_2.toPlainText())
		f_3 = int(self.f_3.toPlainText())
		f_4 = int(self.f_4.toPlainText())
		f_5 = int(self.f_5.toPlainText())
		f_6 = int(self.f_6.toPlainText())
		f_7 = int(self.f_7.toPlainText())
		f_8 = int(self.f_8.toPlainText())
		f_9 = int(self.f_9.toPlainText())
		
		g_1 = int(self.g_1.toPlainText())
		g_2 = int(self.g_2.toPlainText())
		g_3 = int(self.g_3.toPlainText())
		g_4 = int(self.g_4.toPlainText())
		g_5 = int(self.g_5.toPlainText())
		g_6 = int(self.g_6.toPlainText())
		g_7 = int(self.g_7.toPlainText())
		g_8 = int(self.g_8.toPlainText())
		g_9 = int(self.g_9.toPlainText())
		
		h_1 = int(self.h_1.toPlainText())
		h_2 = int(self.h_2.toPlainText())
		h_3 = int(self.h_3.toPlainText())
		h_4 = int(self.h_4.toPlainText())
		h_5 = int(self.h_5.toPlainText())
		h_6 = int(self.h_6.toPlainText())
		h_7 = int(self.h_7.toPlainText())
		h_8 = int(self.h_8.toPlainText())
		h_9 = int(self.h_9.toPlainText())
		
		i_1 = int(self.i_1.toPlainText())
		i_2 = int(self.i_2.toPlainText())
		i_3 = int(self.i_3.toPlainText())
		i_4 = int(self.i_4.toPlainText())
		i_5 = int(self.i_5.toPlainText())
		i_6 = int(self.i_6.toPlainText())
		i_7 = int(self.i_7.toPlainText())
		i_8 = int(self.i_8.toPlainText())
		i_9 = int(self.i_9.toPlainText())
		
		sudoku = {
			'a_1':a_1,'a_2':a_2,'a_3':a_3,'a_4':a_4,'a_5':a_5,'a_6':a_6,'a_7':a_7,'a_8':a_8,'a_9':a_9,
			'b_1':b_1,'b_2':b_2,'b_3':b_3,'b_4':b_4,'b_5':b_5,'b_6':b_6,'b_7':b_7,'b_8':b_8,'b_9':b_9,
			'c_1':c_1,'c_2':c_2,'c_3':c_3,'c_4':c_4,'c_5':c_5,'c_6':c_6,'c_7':c_7,'c_8':c_8,'c_9':c_9,
			'd_1':d_1,'d_2':d_2,'d_3':d_3,'d_4':d_4,'d_5':d_5,'d_6':d_6,'d_7':d_7,'d_8':d_8,'d_9':d_9,
			'e_1':e_1,'e_2':e_2,'e_3':e_3,'e_4':e_4,'e_5':e_5,'e_6':e_6,'e_7':e_7,'e_8':e_8,'e_9':e_9,
			'f_1':f_1,'f_2':f_2,'f_3':f_3,'f_4':f_4,'f_5':f_5,'f_6':f_6,'f_7':f_7,'f_8':f_8,'f_9':f_9,
			'g_1':g_1,'g_2':g_2,'g_3':g_3,'g_4':g_4,'g_5':g_5,'g_6':g_6,'g_7':g_7,'g_8':g_8,'g_9':g_9,
			'h_1':h_1,'h_2':h_2,'h_3':h_3,'h_4':h_4,'h_5':h_5,'h_6':h_6,'h_7':h_7,'h_8':h_8,'h_9':h_9,
			'i_1':i_1,'i_2':i_2,'i_3':i_3,'i_4':i_4,'i_5':i_5,'i_6':i_6,'i_7':i_7,'i_8':i_8,'i_9':i_9,
		};
		
		with open('edit_sudoku.json','wb') as outfile:
			json.dump(sudoku, outfile)
				
		self.close()
		
		
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = MainWindow()
	window.show()
	app.exec_()