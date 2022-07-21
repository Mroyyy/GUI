from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import QDir, QEventLoop
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import QFileDialog, QTextEdit, QAction, QDialog, QSizePolicy, QGridLayout, QLabel

from matplotlib.pyplot import title


import make_dashboard
from make_dashboard import *


class Ui_SecondWindow(object):
    def setupOutput(self, SecondWindow):
        SecondWindow.setObjectName("SecondWindow")
        SecondWindow.resize(800, 600)
        self.OutputWindow = QtWidgets.QWidget(SecondWindow)
        self.OutputWindow.setObjectName("SecondWindow")

        self.label = QtWidgets.QLabel(self.OutputWindow)
        self.label.setGeometry(QtCore.QRect(0, 70, 881, 461))
        self.label.setText("")

        #self.label.setPixmap(QtGui.QPixmap("output_dir/fig1.png"))
        self.label.setPixmap(QtGui.QPixmap("/home/gallegolab/Desktop/GUI/TFM/fig1.png"))
        self.label.setObjectName("label")

        self.pushButton = QtWidgets.QPushButton(self.OutputWindow)
        self.pushButton.setGeometry(QtCore.QRect(700, 520, 89, 25))
        self.pushButton.setObjectName("pushButton")

        SecondWindow.setCentralWidget(self.OutputWindow)
        self.menubar = QtWidgets.QMenuBar(SecondWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 22))
        self.menubar.setObjectName("menubar")
        SecondWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(SecondWindow)
        self.statusbar.setObjectName("statusbar")
        SecondWindow.setStatusBar(self.statusbar)

        self.retranslateUi(SecondWindow)
        QtCore.QMetaObject.connectSlotsByName(SecondWindow)




    def retranslateUi(self, SecondWindow):
        _translate = QtCore.QCoreApplication.translate
        SecondWindow.setWindowTitle(_translate("SecondWindow", "SecondWindow"))
        self.pushButton.setText(_translate("SecondWindow", "Next"))











