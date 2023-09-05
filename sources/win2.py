from PyQt5 import QtCore, QtGui, QtWidgets, uic


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        uic.loadUi('win2.ui', self)
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Course project"))
