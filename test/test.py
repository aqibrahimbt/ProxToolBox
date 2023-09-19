import sys

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *


class MainFrame(QMainWindow):

    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):

        QToolTip.setFont(QFont('SansSerif', 10))

        # Menu bar
        menubar = self.menuBar()

        # Main Options
        fileMenu = menubar.addMenu('&File')
        helpMenu = menubar.addMenu('&Help')

        # Sub Menu for File Menu
        new = QAction('New', self)
        open = QAction('Open', self)
        export = QAction('Export', self)
        fileMenu.addAction(new)
        fileMenu.addAction(open)
        fileMenu.addAction(export)

        # Sub Menu for Help Menu
        gettingStarted = QAction('Getting Started', self)
        examples = QAction('Examples', self)
        tutorials = QAction('Tutorials', self)
        apiReference = QAction('API Reference', self)
        whatisnew = QAction('What is new in Prox Toolbox', self)
        support = QAction('Support Center', self)
        feedback = QAction('Submit feedback', self)
        about = QAction('About', self)
        helpMenu.addAction(gettingStarted)
        helpMenu.addAction(examples)
        helpMenu.addAction(tutorials)
        helpMenu.addAction(apiReference)
        helpMenu.addAction(whatisnew)
        helpMenu.addAction(support)
        helpMenu.addAction(feedback)
        helpMenu.addAction(about)

        # Exit application from Menu or by pressing Ctrl + Q
        exitAct = QAction(QIcon('cancel.png'), '&Exit', self)
        exitAct.setShortcut('Ctrl+Q')
        exitAct.setStatusTip('Exit application')
        exitAct.triggered.connect(qApp.quit)

        newbtn = PicButton(QPixmap('images/new.png'), self)
        newbtn.move(100, 100)
        newbtn.resize(60, 60)
        newbtnlbl = QLabel('New', self)
        newbtnlbl.move(115, 155)

        newbtn.clicked.connect(self.on_pushButton_clicked)

        openbtn = PicButton(QPixmap('images/open.png'), self)
        openbtn.move(250, 100)
        openbtn.resize(60, 60)
        openlbl = QLabel('Open', self)
        openlbl.move(265, 155)

        examplebtn = PicButton(QPixmap('images/example.png'), self)
        examplebtn.move(400, 100)
        examplebtn.resize(60, 60)
        examplebtnlbl = QLabel('Example', self)
        examplebtnlbl.move(405, 155)

        tutorialbtn = PicButton(QPixmap('images/tutorials.png'), self)
        tutorialbtn.move(100, 250)
        tutorialbtn.resize(60, 60)
        tutorialbtnlbl = QLabel('Tutorials', self)
        tutorialbtnlbl.move(105, 305)

        docbtn = PicButton(QPixmap('images/documentation.png'), self)
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

    # Displays the message box.
    def clickMethod(self):
        quitResponse = QMessageBox.question(
            self, 'Achtung!', "Are you sure you want to Quit?",
                                    QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if quitResponse == QMessageBox.Yes:
            sys.exit()
        if quitResponse == QMessageBox.No:
            print('No clicked.')

    def on_pushButton_clicked(self):
        self.w = ComputationFrame()
        self.w.show()


class PicButton(QAbstractButton):
    def __init__(self, pixmap, parent=None):
        super(PicButton, self).__init__(parent)
        self.pixmap = pixmap

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.drawPixmap(event.rect(), self.pixmap)

    def sizeHint(self):
        return self.pixmap.size()


class ComputationFrame(QWidget):

    def __init__(self):
        super().__init__()

        self.items = {
            'Phase': ['AP', 'RAAR', 'RAAR_expert'],
            'ART': ['AP', 'RAAR', 'HAAR', 'QNAP', 'Cimmino']
        }
        self.initUI()

    def initUI(self):

        grid = QGridLayout()
        grid.setSpacing(5)

        # Create ComboBoxes
        self.comboBox = QComboBox(self)
        self.comboBox.addItem("Phase")
        self.comboBox.addItem("ART")
        self.comboBoxlbl = QLabel('MODULE')
        self.comboBox.activated[str].connect(self.onActivated)
        grid.addWidget(self.comboBoxlbl, 1, 0)
        grid.addWidget(self.comboBox, 1, 1)

        self.comboBox1 = QComboBox(self)
        self.comboBox1lbl = QLabel('ALGORITHM')
        grid.addWidget(self.comboBox1lbl, 2, 0)
        grid.addWidget(self.comboBox1, 2, 1, 1, 1)

        self.file = QLabel('FILE')
        self.filebtn = QPushButton("UPLOAD FILE")
        self.filebtn.clicked.connect(self.openFileNameDialog)
        grid.addWidget(self.file, 3, 0)
        grid.addWidget(self.filebtn, 3, 1)

        self.constraint = QLabel('CONSTRAINT')
        self.constraintEdit = QComboBox(self)
        self.constraintEdit.addItem("Convex")
        self.constraintEdit.addItem("support only")
        self.constraintEdit.addItem("real and support")
        self.constraintEdit.addItem("nonnegative and support")
        self.constraintEdit.addItem("amplitude only")
        self.constraintEdit.addItem("sparse real")
        self.constraintEdit.addItem("sparse complex")
        grid.addWidget(self.constraint, 4, 0)
        grid.addWidget(self.constraintEdit, 4, 1)

        self.beta_0 = QLabel('BETA_0')
        self.beta_0Edit = QLineEdit()
        self.beta_0Edit.setPlaceholderText(' starting relaxation parameter')
        grid.addWidget(self.beta_0, 5, 0)
        grid.addWidget(self.beta_0Edit, 5, 1)

        self.beta_max = QLabel('BETA_MAX')
        self.beta_maxEdit = QLineEdit()
        self.beta_maxEdit.setPlaceholderText(' maximum relaxation parameter')
        grid.addWidget(self.beta_max, 6, 0)
        grid.addWidget(self.beta_maxEdit, 6, 1)

        self.beta_switch = QLabel('BETA_SWITCH')
        self.beta_switchEdit = QLineEdit()
        self.beta_switchEdit.setPlaceholderText(' iteration at which beta moves from beta_0 -> beta_max')
        grid.addWidget(self.beta_switch, 7, 0)
        grid.addWidget(self.beta_switchEdit, 7, 1)

        self.max_iter = QLabel('MAX ITER.')
        self.max_iterEdit = QLineEdit()
        self.max_iterEdit.setPlaceholderText(' maximum number of iterations')
        grid.addWidget(self.max_iter, 8, 0)
        grid.addWidget(self.max_iterEdit, 8, 1)

        self.tol = QLabel('TOLERANCE')
        self.tolEdit = QLineEdit()
        self.tolEdit.setPlaceholderText(' maximum tolerance')
        grid.addWidget(self.tol, 9, 0)
        grid.addWidget(self.tolEdit, 9, 1)

        # Buttons
        self.submitbtn = QPushButton('SUBMIT')
        self.submitbtn.resize(self.submitbtn.sizeHint())
        self.submitbtn.setStyleSheet("background-color: #17a2b8")
        grid.addWidget(self.submitbtn, 10, 1)
        self.submitbtn.clicked.connect(self.on_click)

        self.setLayout(grid)
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

    def onActivated(self, text):
        print(self.items[text])
        self.comboBox1.clear()
        self.comboBox1.addItems(self.items[text])

    def openFileNameDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(
            self,"QFileDialog.getOpenFileName()",
            "", "All Files (*);;Python Files (*.py)", options=options)
        if fileName:
            print(fileName)

    def saveFileDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(
            self, "QFileDialog.getSaveFileName()"
            , "", "All Files (*);;Text Files (*.txt)", options=options)
        if fileName:
            print(fileName)

    @pyqtSlot()
    def on_click(self):
        comboBox = self.comboBox.currentText()
        comboBox1 = self.comboBox1.currentText()
        constraintEdit = self.constraintEdit.currentText()
        beta_0Edit = self.beta_0Edit.text()
        beta_maxEdit = self.beta_maxEdit.text()
        beta_switchEdit = self.beta_switchEdit.text()
        max_iterEdit = self.max_iterEdit.text()
        tolEdit = self.tolEdit.text()

        print(
            f'{comboBox} {comboBox1} {constraintEdit} {beta_0Edit} {beta_maxEdit} {beta_switchEdit} {max_iterEdit} {tolEdit}'
        )


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = MainFrame()
    sys.exit(app.exec_())
