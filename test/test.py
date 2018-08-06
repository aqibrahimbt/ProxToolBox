import sys

from PyQt5.QtCore import QRect
from PyQt5.QtWidgets import (QWidget, QToolTip, QPushButton, QApplication, QDesktopWidget, QMessageBox, QAction, qApp,
                             QMainWindow, QMenu, QAbstractButton, QLabel, QComboBox)
from PyQt5.QtGui import (QFont, QIcon, QPainter, QPixmap)


class MainFrame(QMainWindow):

    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):

        QToolTip.setFont(QFont('SansSerif', 10))

        # Menu bar
        menubar = self.menuBar()

        # Main Options
        fileMenu = menubar.addMenu('&Start')
        helpMenu = menubar.addMenu('&Help')

        # Sub Menu for File Menu
        problem = QMenu('Problem', self)
        ct = QAction('CT', self)
        phase = QAction('Phase',  self)
        problem.addAction(ct)
        problem.addAction(phase)
        fileMenu.addMenu(problem)

        algorithm = QMenu('Algorithm', self)
        ap_default = QAction('AP (Default)', self)
        ap_expert = QAction('AP (Expert)', self)
        graal = QAction('GRAAL', self)
        haar = QAction('HAAR', self)
        hpr = QAction('HPR', self)
        palm = QAction('PALM', self)
        qnap = QAction('QNAP', self)
        raar = QAction('RAAR (Default)', self)
        raar_expert = QAction('RAAR (Expert)', self)
        algorithm.addAction(ap_default)
        algorithm.addAction(ap_expert)
        algorithm.addAction(graal)
        algorithm.addAction(haar)
        algorithm.addAction(hpr)
        algorithm.addAction(palm)
        algorithm.addAction(qnap)
        algorithm.addAction(raar)
        algorithm.addAction(raar_expert)
        fileMenu.addMenu(algorithm)

        export = QAction('Export', self)
        fileMenu.addAction(export)

        # Sub Menu for Help Menu
        gettingStarted = QAction('Getting Started', self)
        apiReference = QAction('API Reference', self)
        whatisnew = QAction('What is new in Prox Toolbox', self)
        support = QAction('Support Center', self)
        feedback = QAction('Submit feedback', self)
        helpMenu.addAction(gettingStarted)
        helpMenu.addAction(apiReference)
        helpMenu.addAction(whatisnew)
        helpMenu.addAction(support)
        helpMenu.addAction(feedback)

        # Exit application from Menu or by pressing Ctrl + Q
        exitAct = QAction(QIcon('cancel.png'), '&Exit', self)
        exitAct.setShortcut('Ctrl+Q')
        exitAct.setStatusTip('Exit application')
        exitAct.triggered.connect(qApp.quit)

        #self.toolbar = self.addToolBar('Exit')
        #self.toolbar.addAction(exitAct)

        # # Buttons
        # btn = QPushButton('Algorithm', self)
        # btn.setToolTip('Click to select <b>Algorithm</b>')
        # btn.resize(btn.sizeHint())
        # btn.move(50, 50)
        #
        # qbtn = QPushButton('Quit', self)
        # qbtn.clicked.connect(self.clickMethod)
        # qbtn.resize(qbtn.sizeHint())
        # qbtn.move(200, 50)

        newbtn = PicButton(QPixmap('images/paper.png'), self)
        newbtn.move(100, 100)
        newbtn.resize(60, 60)
        newbtnlbl = QLabel('New', self)
        newbtnlbl.move(115, 155)

        openbtn = PicButton(QPixmap('images/folder.png'), self)
        openbtn.move(250, 100)
        openbtn.resize(60, 60)
        openlbl = QLabel('Open', self)
        openlbl.move(265, 155)

        examplebtn = PicButton(QPixmap('images/tutorial.png'), self)
        examplebtn.move(400, 100)
        examplebtn.resize(60, 60)
        examplebtnlbl = QLabel('Example', self)
        examplebtnlbl.move(405, 155)

        tutorialbtn = PicButton(QPixmap('images/006-notes.png'), self)
        tutorialbtn.move(100, 250)
        tutorialbtn.resize(60, 60)
        tutorialbtnlbl = QLabel('Tutorials', self)
        tutorialbtnlbl.move(105, 305)

        docbtn = PicButton(QPixmap('images/website.png'), self)
        docbtn.move(250, 250)
        docbtn.resize(60, 60)
        docbtnlbl = QLabel('Documentation', self)
        docbtnlbl.move(233, 305)

        exitbtn = PicButton(QPixmap('images/exit-2.png'), self)
        exitbtn.move(405, 250)
        exitbtn.resize(60, 60)
        exitbtnlbl = QLabel('Quit', self)
        exitbtnlbl.move(412, 305)

        exitbtn.clicked.connect(self.clickMethod)

        self.setAutoFillBackground(True)
        self.center()
        self.setWindowTitle('Prox ToolBox')
        self.setFixedSize(570, 450)
        self.setWindowIcon(QIcon('images/domain.png'))
        self.show()

    # Centers the screen to the middle of the desktop
    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

    # Displays the message box. Can be used for something else later
    def clickMethod(self):
        quitResponse = QMessageBox.question(self, 'Achtung!', "Are you sure you want to Quit?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if quitResponse == QMessageBox.Yes:
            sys.exit()
        if quitResponse == QMessageBox.No:
            print('No clicked.')


class PicButton(QAbstractButton):
    def __init__(self, pixmap, parent=None):
        super(PicButton, self).__init__(parent)
        self.pixmap = pixmap

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.drawPixmap(event.rect(), self.pixmap)

    def sizeHint(self):
        return self.pixmap.size()


class ComputationFrame(QMainWindow):

    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):

        centralWidget = QWidget(self)
        self.setCentralWidget(centralWidget)

        # Create combobox and add items.
        self.comboBox = QComboBox(centralWidget)
        self.comboBox.setGeometry(QRect(40, 40, 491, 31))
        self.comboBox.setObjectName(("comboBox"))
        self.comboBox.addItem("PyQt")
        self.comboBox.addItem("Qt")
        self.comboBox.addItem("Python")
        self.comboBox.addItem("Example")

        self.setAutoFillBackground(True)
        self.center()
        self.setWindowTitle('Prox ToolBox')
        self.setFixedSize(570, 450)
        self.setWindowIcon(QIcon('images/domain.png'))
        self.show()

    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())


if __name__ == '__main__':

    app = QApplication(sys.argv)
    ex = ComputationFrame()
    sys.exit(app.exec_())