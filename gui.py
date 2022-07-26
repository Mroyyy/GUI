# Pyqt5 import functions
from PyQt5 import QtCore, QtGui, QtWidgets, QtWebEngineWidgets
from PyQt5.QtCore import QDir, QEventLoop, Qt, QUrl
from PyQt5.QtWidgets import QFileDialog, QTextEdit, QAction, QDialog, QSizePolicy, QGridLayout, QSpacerItem, \
    QScrollArea, QWidget, QVBoxLayout, QLabel, QSizeGrip, QPlainTextEdit, QPushButton

##########################
# Ferran's script functions
##########################

import Bio.SeqIO as IO

import plotly
from matplotlib.backends.backend_template import FigureCanvas
from plotly.subplots import make_subplots
from pygments.lexers import go


from bin.blast import *
from bin.PDB_retriever import *
import bin.config as cfg
import logging as l
from pathlib import Path, PurePosixPath
import os
import shutil

from bin.dashboard.dashboard_functions import read_DFI_csvs, read_hng_files, read_compsite_files
from bin.graphical_summary import StructuReport
from bin.extract_flexible_residues import extract_residue_list
from bin.process_predicted_model import *
from matplotlib import pyplot as plt
import fnmatch
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor, CellExecutionError
from bin.custom_top import make_rb_list, make_composite, write_custom_topology
import pandas as pd

import os
from bin.custom_top import  write_custom_topology
from pathlib import Path

# File management/OS
from pathlib import Path, PurePosixPath


#Plotting
# import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd

## DASH
from dash import Dash, dcc, html, Input, Output, State
import subprocess

#Custom topology
from bin.custom_top import make_rb_list

import plotly.io as pio
import plotly.express as px


class Ui_HelpWindow(object):
    '''
    Set up Help Window to provide a step by step usage
    Shortcut: Ctrl + H
    '''
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("HelpWindow")
        MainWindow.resize(1150, 610)
        MainWindow.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(40, 0, 1091, 91))
        font = QtGui.QFont()
        font.setFamily("Chandas")
        font.setPointSize(14)
        font.setBold(True)
        font.setWeight(75)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(9, 82, 1116, 423))
        font = QtGui.QFont()
        font.setFamily("Chandas")
        self.label.setFont(font)
        self.label.setObjectName("label")

        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1150, 22))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("HelpWindow", "HelpWindow"))
        self.label_2.setText(_translate("HelpWindow", "IMPORTANT CONSIDERATIONS"))
        self.label.setText(_translate("HelpWindow", "<html><head/><body><p align=\"justify\"><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; color:#24292f; background-color:#ffffff;\">This tool automatically combines structures derived from </span><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; font-weight:600; color:#24292f; background-color:#ffffff;\">experimental</span><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; color:#24292f; background-color:#ffffff;\"> (Protein Data Bank) and </span><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; font-weight:600; color:#24292f; background-color:#ffffff;\">computational</span><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; color:#24292f; background-color:#ffffff;\"> methods (AlphaFold and RoseTTaFold) for building an </span></p><p align=\"justify\"><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; color:#24292f; background-color:#ffffff;\">atomic structure composite to be used as input for the </span><a href=\"https://integrativemodeling.org/\"><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; text-decoration: underline; color:#000000; background-color:#ffffff;\">Integrative modeling Platform</span></a><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; color:#24292f; background-color:#ffffff;\"> in the form of a Topology File, maximizing the coverage of high resolution structures </span></p><p align=\"justify\"><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; color:#24292f; background-color:#ffffff;\">with respect to a target sequence, and providing flexibility predictions. This tool will increase the accuracy of the models generated by IMP.</span></p><p align=\"justify\"><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; color:#24292f; background-color:#ffffff;\">This program uses external programs to work, so you need to have them installed:</span></p><p align=\"justify\"><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; color:#24292f; background-color:#ffffff;\">- </span><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; font-weight:600; color:#24292f; background-color:#ffffff;\">BLAST (mandatory)</span></p><p align=\"justify\"><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; color:#24292f; background-color:#ffffff;\">- AlphaFold</span></p><p align=\"justify\"><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; color:#24292f; background-color:#ffffff;\">- RoseTTaFold</span></p><p align=\"justify\"><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; color:#24292f; background-color:#ffffff;\">Also, in order to install the Python packages you will need </span><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; font-weight:600; color:#24292f; background-color:#ffffff;\">Conda</span><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; color:#24292f; background-color:#ffffff;\"> and </span><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; font-weight:600; color:#24292f; background-color:#ffffff;\">pip</span><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; color:#24292f; background-color:#ffffff;\">.</span></p><p align=\"justify\"><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; color:#24292f; background-color:#ffffff;\">Positional arguments:</span></p><p align=\"justify\"><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; color:#24292f; background-color:#ffffff;\">- FASTA input sequence in FASTA format. </span><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; font-weight:600; color:#24292f;\">It needs to have the same file name and header e. g. the file SEC3.fasta starts with the line &quot;&gt;SEC3&quot;</span></p><p align=\"justify\"><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; color:#24292f; background-color:#ffffff;\">- Output directory to store the retrieved PDBs</span></p><p align=\"justify\"><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; color:#24292f; background-color:#ffffff;\">Optional arguments:</span></p><p align=\"justify\"><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; color:#24292f; background-color:#ffffff;\">- AlphaFold2 model, RoseTTaFold model and PAE JSON (from AF-EBI server) in PDB format</span></p><p align=\"justify\"><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; color:#24292f; background-color:#ffffff;\">- Custom Templates: A custom experimentally solved PDB provided by the user</span></p><p align=\"justify\"><span style=\" font-family:\'-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif,Apple Color Emoji,Segoe UI Emoji\'; font-size:10pt; color:#24292f; background-color:#ffffff;\">- Run AlphaFold or RoseTTaFold: Send a batch script using SLURM or run locally</span></p></body></html>"))


class Ui_MainWindow(object):
    def openWindow(self):
        '''
        Creates connection between input window and output window (second window)
        '''
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_SecondWindow()
        self.ui.setupOutput(self.window)
        self.window.show()

    def show_help(self):
        '''
        Opens Ui_HelpWindow
        '''
        self.help = QtWidgets.QMainWindow()
        self.hlp = Ui_HelpWindow()
        self.hlp.setupUi(self.help)
        self.help.show()

    def setupUi(self, MainWindow):
        '''
        Set up Main Window where files and output are selected
        '''
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1069, 614)
        MainWindow.setSizeIncrement(QtCore.QSize(0, 0))
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.verticalLayout_10 = QtWidgets.QVBoxLayout()
        self.verticalLayout_10.setObjectName("verticalLayout_10")
        self.label_9 = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Chandas")
        font.setPointSize(12)
        font.setBold(True)
        font.setWeight(75)
        self.label_9.setFont(font)
        self.label_9.setObjectName("label_9")
        self.verticalLayout_10.addWidget(self.label_9)
        self.gridLayout_3.addLayout(self.verticalLayout_10, 0, 2, 1, 8)
        self.verticalLayout_7 = QtWidgets.QVBoxLayout()
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.toolButton_4 = QtWidgets.QToolButton(self.centralwidget)
        self.toolButton_4.setObjectName("toolButton_4")
        self.verticalLayout_7.addWidget(self.toolButton_4)
        self.toolButton_5 = QtWidgets.QToolButton(self.centralwidget)
        self.toolButton_5.setObjectName("toolButton_5")
        self.verticalLayout_7.addWidget(self.toolButton_5)
        self.toolButton_6 = QtWidgets.QToolButton(self.centralwidget)
        self.toolButton_6.setObjectName("toolButton_6")
        self.verticalLayout_7.addWidget(self.toolButton_6)
        self.gridLayout_3.addLayout(self.verticalLayout_7, 2, 13, 1, 1)
        self.gridLayout_2 = QtWidgets.QGridLayout()
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.label_7 = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Chandas")
        self.label_7.setFont(font)
        self.label_7.setObjectName("label_7")
        self.gridLayout_2.addWidget(self.label_7, 0, 0, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Chandas")
        self.label_6.setFont(font)
        self.label_6.setObjectName("label_6")
        self.gridLayout_2.addWidget(self.label_6, 1, 0, 1, 1)
        self.label_5 = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Chandas")
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")
        self.gridLayout_2.addWidget(self.label_5, 2, 0, 1, 1)
        self.gridLayout_3.addLayout(self.gridLayout_2, 2, 9, 1, 2)
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.label_4 = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Chandas")
        self.label_4.setFont(font)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 4, 0, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Chandas")
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 2, 0, 1, 1)
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Chandas")
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 3, 0, 1, 1)
        self.label = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Chandas")
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 1, 0, 1, 1)
        self.gridLayout_3.addLayout(self.gridLayout, 2, 1, 1, 4)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.lineEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit.setSizeIncrement(QtCore.QSize(0, 0))
        font = QtGui.QFont()
        font.setFamily("Chandas")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.lineEdit.setFont(font)
        self.lineEdit.setObjectName("lineEdit")
        self.verticalLayout.addWidget(self.lineEdit)
        self.horizontalLayout_2.addLayout(self.verticalLayout)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.pushButton = QtWidgets.QPushButton(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Chandas")
        self.pushButton.setFont(font)
        self.pushButton.setObjectName("pushButton")
        self.horizontalLayout.addWidget(self.pushButton)
        self.horizontalLayout_2.addLayout(self.horizontalLayout)
        self.gridLayout_3.addLayout(self.horizontalLayout_2, 1, 0, 1, 14)
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.label_8 = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Chandas")
        font.setItalic(True)
        self.label_8.setFont(font)
        self.label_8.setObjectName("label_8")
        self.verticalLayout_3.addWidget(self.label_8)
        self.gridLayout_3.addLayout(self.verticalLayout_3, 6, 2, 1, 1)
        self.verticalLayout_4 = QtWidgets.QVBoxLayout()
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.lineEdit_2 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_2.setObjectName("lineEdit_2")
        self.verticalLayout_4.addWidget(self.lineEdit_2)
        self.lineEdit_3 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_3.setObjectName("lineEdit_3")
        self.verticalLayout_4.addWidget(self.lineEdit_3)
        self.lineEdit_4 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_4.setObjectName("lineEdit_4")
        self.verticalLayout_4.addWidget(self.lineEdit_4)
        self.lineEdit_5 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_5.setObjectName("lineEdit_5")
        self.verticalLayout_4.addWidget(self.lineEdit_5)
        self.gridLayout_3.addLayout(self.verticalLayout_4, 2, 5, 1, 2)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.pushButton_2 = QtWidgets.QPushButton(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Chandas")
        self.pushButton_2.setFont(font)
        self.pushButton_2.setObjectName("pushButton_2")
        self.horizontalLayout_3.addWidget(self.pushButton_2)
        self.gridLayout_3.addLayout(self.horizontalLayout_3, 4, 12, 4, 3)
        self.verticalLayout_5 = QtWidgets.QVBoxLayout()
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.lineEdit_6 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_6.setObjectName("lineEdit_6")
        self.verticalLayout_5.addWidget(self.lineEdit_6)
        self.lineEdit_7 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_7.setObjectName("lineEdit_7")
        self.verticalLayout_5.addWidget(self.lineEdit_7)
        self.lineEdit_8 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_8.setObjectName("lineEdit_8")
        self.verticalLayout_5.addWidget(self.lineEdit_8)
        self.gridLayout_3.addLayout(self.verticalLayout_5, 2, 11, 1, 2)
        self.verticalLayout_6 = QtWidgets.QVBoxLayout()
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.toolButton_2 = QtWidgets.QToolButton(self.centralwidget)
        self.toolButton_2.setObjectName("toolButton_2")
        self.verticalLayout_6.addWidget(self.toolButton_2)
        self.toolButton = QtWidgets.QToolButton(self.centralwidget)
        self.toolButton.setObjectName("toolButton")
        self.verticalLayout_6.addWidget(self.toolButton)
        self.toolButton_3 = QtWidgets.QToolButton(self.centralwidget)
        self.toolButton_3.setObjectName("toolButton_3")
        self.verticalLayout_6.addWidget(self.toolButton_3)
        self.toolButton_7 = QtWidgets.QToolButton(self.centralwidget)
        self.toolButton_7.setObjectName("toolButton_7")
        self.verticalLayout_6.addWidget(self.toolButton_7)
        self.gridLayout_3.addLayout(self.verticalLayout_6, 2, 7, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1069, 22))
        self.menubar.setObjectName("menubar")
        self.menuHelp = QtWidgets.QMenu(self.menubar)
        self.menuHelp.setObjectName("menuHelp")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionHelp = QtWidgets.QAction(MainWindow)
        self.actionHelp.setObjectName("actionHelp")
        self.actionHelp.triggered.connect(self.show_help)
        self.actionHelp.setShortcut("Ctrl+H")


        self.menuHelp.addAction(self.actionHelp)
        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.label_9.setText(
            _translate("MainWindow", "Automated structural information retrieval for Integrative Modeling"))
        self.toolButton_4.setText(_translate("MainWindow", "...")) # Custom templates
        self.toolButton_5.setText(_translate("MainWindow", "...")) # Run RosettaFold
        self.toolButton_6.setText(_translate("MainWindow", "...")) # RosettaFold model
        self.toolButton_7.setText(_translate("MainWindow", "...")) # AlphaFold PAE JSON
        self.label_7.setText(_translate("MainWindow", "Custom Templates*"))
        self.label_6.setText(_translate("MainWindow", "Run RosettaFold*"))
        self.label_5.setText(_translate("MainWindow", "RosettaFold (PDB format)*"))
        self.label_4.setText(_translate("MainWindow", "Alphafold PAE JSON file*"))
        self.label_2.setText(_translate("MainWindow", "Alphafold model*"))
        self.label_3.setText(_translate("MainWindow", "Run AlphaFold*"))
        self.label.setText(_translate("MainWindow", "Output directory"))
        self.toolButton_2.setText(_translate("MainWindow", "...")) # Output directory
        self.toolButton.setText(_translate("MainWindow", "...")) # Alphamodel
        self.toolButton_3.setText(_translate("MainWindow", "...")) # Run AlphaFold
        self.lineEdit.setText(_translate("MainWindow", "Select sequence file and output directory"))
        self.pushButton.setText(_translate("MainWindow", "Select file"))
        self.label_8.setText(_translate("MainWindow", "Arguments with * are not required (optional)"))
        self.pushButton_2.setText(_translate("MainWindow", "START"))
        self.menuHelp.setTitle(_translate("MainWindow", "Help"))
        self.actionHelp.setText(_translate("MainWindow", "Help"))
        self.actionHelp.setShortcut("Ctrl+H")

        self.toolButton_2.clicked.connect(self.create_dir)
        self.pushButton_2.clicked.connect(self.run) # START BUTTON
        # click button and call function
        self.pushButton.clicked.connect(self.get_file)

    def get_file(self):
        '''
        Function that gets filename of fasta, its directory, and the sequence itself
        '''
        global source
        global filename
        global query_name
        global fasta
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        name, _ = QFileDialog.getOpenFileName(None, 'QFileDialog.getOpenFileName()', "")
        filename = name.split("/")[-1]  # get filename from path
        source = name
        file = open(name, 'r')
        query_name, fasta = list(file)
        # query_name = print(list(file))
        query_name = query_name.replace(">", "").strip()

        return self.lineEdit.setText(source), filename, query_name, fasta

    def create_dir(self):
        '''
        Creates output directory
        '''
        global output_dir
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        name, _ = QFileDialog.getSaveFileName(None, 'QFileDialog.getSaveFileName()', "")
        output_dir = name

        return self.lineEdit_2.setText(output_dir)

    def run(self):
        '''
        Main program/function
        :return: Atomic structure composite
        '''
        verbose = True

        ### Initializing the LOG system ###
        global logdir
        global plot_dir
        logdir = os.path.join(output_dir, query_name, "LOG", "")
        Path(logdir).mkdir(parents=True, exist_ok=True)
        l.basicConfig(format="%(levelname)s:%(message)s",
                      filename=os.path.join(logdir, f"{query_name}.log"),
                      level = l.DEBUG)


        l.debug("...STARTING...\n")

        # If verbose is set, the LOG file is also printed in STDOUT
        if verbose:
            l.getLogger().addHandler(l.StreamHandler())

        record_dict = IO.to_dict(IO.parse(source, "fasta"))
        if len(record_dict.keys()) > 1:
            raise NotImplemented("This program only works with one sequence at a time")
        if len(record_dict.keys()) == 1:
            if list(record_dict.keys())[0] != query_name:
                raise NameError(f"""Please, make sure your filename and fasta identifier 
                coincide. filename: {query_name} / ID name: {record_dict.keys()}""")

            query_length = len(record_dict[query_name].seq)
            l.info(f"Query length: {query_length}")

        l.info(f"Create folders:")
        blast_dir = os.path.join(output_dir, query_name, "BLAST", "")
        pdb_dir = os.path.join(output_dir, query_name, "PDB", "")
        fasta_dir = os.path.join(output_dir, query_name, "FASTA", "")
        report_dir = os.path.join(output_dir, query_name, "REPORT", "")
        hinges_dir = os.path.join(output_dir, query_name, "HINGES", "")
        IMP_dir = os.path.join(output_dir, query_name, "IMP", "")
        plot_dir = os.path.join(output_dir, query_name, "PLOTS", "")
        l.info(f"{pdb_dir}, {fasta_dir}, {report_dir}, {hinges_dir}, {IMP_dir} .")

        Path(blast_dir).mkdir(parents=True, exist_ok=True)
        Path(pdb_dir).mkdir(parents=True, exist_ok=True)
        Path(fasta_dir).mkdir(parents=True, exist_ok=True)
        Path(report_dir).mkdir(parents=True, exist_ok=True)
        Path(hinges_dir).mkdir(parents=True, exist_ok=True)
        Path(IMP_dir).mkdir(parents=True, exist_ok=True)
        Path(plot_dir).mkdir(parents=True, exist_ok=True)


        ## 1. Check if the input sequence is already in the PDB

        l.info("### BLAST ###")

        l.info(f"The BLAST output will be stored in:{blast_dir}")

        # Run BLAST
        # from bin.blast import query as q_blast
        # q_blast = query_name
        outblast = run_blast_local(source, blast_dir)

        # Catch exact matches
        exact_matches = exact_match_retriever(outblast)
        l.info(f""" The target sequence has close homologs in the PDB with 
            code/s: {exact_matches.keys()}""")

        structures_for_query = []
        l.info(f"Structures for query: {structures_for_query}")

        # Retrieve from the PDB
        l.info("Retrieving structures from the pdb")
        if exact_matches:
            retrieve_pdb_info(exact_matches, pdb_dir, fasta_dir)
            # Check lengths of the actual PDB Chains and store them accordingly
            for file in os.listdir(pdb_dir):
                current = os.path.join(pdb_dir, file)
                if os.path.isfile(current):
                    l.info(f"File being processed: {file}")
                    identifier = file.split(".")[0].upper()

                    # Make the directory for the chains
                    chain_dir = os.path.join(pdb_dir, "CHAINS", "")
                    l.info(f"Making directory for the chains at {chain_dir}")
                    Path(chain_dir).mkdir(parents=True, exist_ok=True)

                    # Extract the desired chain
                    l.info(f"Extracting the chain")
                    splitter = ChainSplitter(mmcif=True, out_dir=chain_dir)
                    chain_path = splitter.make_pdb(os.path.join(pdb_dir, file),
                                                   exact_matches[identifier], overwrite=True)

                    # Store partial matches (<95% of the query length)
                    pdb_len = check_PDB_len(chain_path, exact_matches[identifier])
                    l.info(f"""Length of the template {PurePosixPath(chain_path).name}: 
                                        {pdb_len}""")

                    l.info(f"PDB_LEN: {pdb_len} . QUERY_LEN: {query_length}")
                    if pdb_len > 10 and pdb_len < (0.95 * query_length):
                        l.info(f"""{PurePosixPath(chain_path).name} has length {pdb_len}, 
                        it will be stored as a partial match""")
                        newpath = os.path.join(pdb_dir, "partial",
                                               f"{PurePosixPath(chain_path).name}")
                        try:
                            l.info(f"MOVING {chain_path} TO {newpath}")
                            shutil.move(chain_path, newpath)
                            structures_for_query.append(newpath)
                        except Exception:
                            directory = os.path.join(pdb_dir, "partial", "")
                            l.info(f"\"{directory}\" does not exist, it will be created")
                            Path(directory).mkdir(parents=True, exist_ok=True)
                            shutil.move(chain_path, newpath)
                            structures_for_query.append(newpath)
                    if pdb_len > 10 and pdb_len > (0.95 * query_length):
                        l.info(f"""{PurePosixPath(chain_path).name} has length {pdb_len}, 
                        it will be stored as a full-length match""")
                        newpath = os.path.join(pdb_dir, "total",
                                               f"{PurePosixPath(chain_path).name}")
                        try:
                            shutil.move(chain_path, newpath)
                            structures_for_query.append(newpath)
                        except Exception:
                            directory = os.path.join(pdb_dir, "total", "")
                            l.info(f"\"{directory}\" does not exist, it will be created")
                            Path(directory).mkdir(parents=True, exist_ok=True)
                            shutil.move(chain_path, newpath)
                            structures_for_query.append(newpath)



        l.info(f"CONFIDENT FILES: {structures_for_query}")
        nrow = len(structures_for_query)

        l.info("### HINGE DETECTION and DFI ###")
        for structure in structures_for_query:
            reporter = StructuReport(structure)
            # Get coverage of the structure; fasta variable
            coverage_df = reporter.get_coverage(source, save_csv=True, outdir=report_dir)
            # Get hinges and save the .hng files
            hinges = reporter.get_hinges(alpha_range=cfg.PACKMANconfig["alpha_range"],
                                         save_hng=True, outdir=hinges_dir)
            # Get DFI
            dfi_df = reporter.get_dfi_coverage(reference_fasta=source, save_csv=True, outdir=report_dir)

        ## Write Topology file ##

        rigid_bodies = make_rb_list(structures_for_query, source)

        composite_rb = make_composite(rigid_bodies)

        # Export the composite coverage in .csv
        out_path = os.path.join(report_dir, "COVERAGE", f"{query_name}_composite_coverage.csv")
        i = 0
        for rb in composite_rb:
            coverage = rb.get_coverage()
            if i == 0:
                composite_coverage = coverage
                i += 1
            else:
                if i == 1:
                    merged_left = pd.merge(left=composite_coverage, right=coverage,
                                           how="left", left_on="ResID", right_on="ResID")
                    i += 1
                else:
                    merged_left = pd.merge(left=merged_left, right=coverage,
                                           how="left", left_on="ResID", right_on="ResID")
                    i += 1

        if i == 1:
            composite_coverage.to_csv(out_path, encoding='utf-8',
                                      index=False, float_format='%.3f')
        elif i > 1:
            merged_left.to_csv(out_path, encoding='utf-8',
                               index=False, float_format='%.3f')

        # Convert to list and sort by the ones who start earlier in the sequence
        composite_rb.sort(key=lambda x: x.residue_range[0])

        # Write the topology file;
        write_custom_topology(os.path.join(IMP_dir, f"{query_name}.topology"), composite_rb)
        os.rmdir("obsolete")
        os.remove("DCI_pymol_output.txt")

        ui.openWindow()


class Ui_SecondWindow(object):
    def show_plot(self):
        '''
        Connection function: show coverage plot
        Calls update_graph() at the beginning in order to save the image and then show it
        '''

        self.update_graph()
        self.plot_out = QtWidgets.QMainWindow()
        self.pl = Ui_plot()
        self.pl.setup(self.plot_out)
        self.plot_out.show()

    def show_Secondplot(self):
        '''
        Connection function: show hinges and flexibility plot
        Calls update_Secondgraph() at the beginning in order to save the image and then show it
        '''
        self.update_Secondgraph()
        self.plot_Secondout = QtWidgets.QMainWindow()
        self.plSec = Ui_Secondplot()
        self.plSec.setup(self.plot_Secondout)
        self.plot_Secondout.show()

    def update_graph(self):
        '''
        Function that creates coverage plot, reads data from output directory and makes graph
        :return: coverage plot
        '''
        global fig1
        global png1_html
        global query_name
        i = 0
        df_list = []
        structure_list = []
        for child in Path(os.path.join(output_dir, query_name, "REPORT", "COVERAGE")).iterdir():
            if child.is_file() and "composite" not in str(child):
                i += 1
                df = pd.read_csv(child)
                df_list.append(df)
                structure_list.append(child)

        fig1 = make_subplots(rows=i, cols=1, shared_xaxes=True, x_title="Residue position")

        i = 1
        for df in df_list:
            fig1.append_trace(go.Scatter(
                x=df[df.columns[0]],  # ResID
                y=df[df.columns[1]],
                fill='tozeroy',
                name=str(structure_list[i - 1])
            ), row=i, col=1)
            i += 1

        # fig1.update_layout(height=400, width=1000, title_text="Coverage")
        fig1.update_layout(title_text="Coverage")
        fig1.update_yaxes(showgrid=False, range=[0, 1], showticklabels=False)

        # function that saves image as .png
        png1 = pio.write_image(fig1, f"{plot_dir}/coverage_plot.png", scale=1, width=1400, height=700)
        png1_html = pio.write_html(fig1, f"{plot_dir}/coverage_plot.html")

        return fig1

    def update_Secondgraph(self):
        '''
        Function that reads data from output directory and predicts hinges and dynamic flexibility
        :return: DFI profiles + Predicted hinges
        '''
        dfi_dict = read_DFI_csvs(os.path.join(output_dir, query_name, "REPORT", "DFI"))
        hng_dict = read_hng_files(os.path.join(output_dir, query_name, "HINGES"))
        dfi_files = [dfi_file for dfi_file in dfi_dict.keys() if
                     "AF_DFI" not in str(dfi_file.stem) and "RF_DFI" not in str(dfi_file.stem)]

        fig2 = make_subplots(
            rows=len(dfi_files), cols=1, shared_xaxes=True,
            x_title="Residue position"
        )

        i = 1
        for dfi_file in dfi_files:
            df = dfi_dict[dfi_file]
            fig2.append_trace(go.Scatter(
                x=df[df.columns[0]],  # resIDs
                y=df[df.columns[1]],  # pctdfi
                name=str(dfi_file)
            ), row=i, col=1)
            j = 1
            for hng_file in hng_dict.keys():
                if str(PurePosixPath(dfi_file).stem)[0:-13] == str(PurePosixPath(hng_file).stem):
                    for hinge in hng_dict[hng_file]:
                        fig2.add_vrect(
                            x0=hinge.split(':')[0],
                            x1=hinge.split(':')[1],
                            annotation_text=f"H{j}", annotation_position="top left",
                            fillcolor="#52BE80", opacity=0.2,
                            layer="below", line_width=0,
                            row=i, col=1)
                        j += 1
            i += 1
        fig2.update_layout(title_text="DFI profiles + Predicted hinges",
                           margin_pad=10, barmode="group", legend=dict(orientation="h", y=-0.35))
        fig2.update_yaxes(showgrid=False, range=[0, 1], nticks=2)

        png2 = pio.write_image(fig2, f"{plot_dir}/hinges_prediction.png", scale=1, width=1400,
                               height=700)
        png2_html = pio.write_html(fig2, f"{plot_dir}/hinges_prediction.html")

        return fig2

    def openThirdWindow(self):
        '''
        Creates connection between input window and output window (second window)
        '''
        self.Thirdwindow2 = QtWidgets.QMainWindow()
        self.ui = Ui_ThirdWindow()
        self.ui.setupOutput(self.Thirdwindow2)
        self.Thirdwindow2.show()


    def setupOutput(self, SecondWindow):
        '''
        set up second window, when introduced input fasta and output directory, window that will be opened to show
        results
        '''
        SecondWindow.setObjectName("SecondWindow")
        SecondWindow.resize(1059, 686)
        self.OutputWindow = QtWidgets.QWidget(SecondWindow)
        self.OutputWindow.setObjectName("OutputWindow")
        self.label = QtWidgets.QLabel(self.OutputWindow)
        self.label.setGeometry(QtCore.QRect(40, 80, 971, 61))
        font = QtGui.QFont()
        font.setFamily("Chandas")
        font.setPointSize(10)
        self.label.setFont(font)
        self.label.setWordWrap(True)
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(self.OutputWindow)
        self.label_2.setGeometry(QtCore.QRect(40, 40, 161, 31))
        font = QtGui.QFont()
        font.setFamily("Chandas")
        font.setPointSize(16)
        font.setBold(True)
        font.setWeight(75)
        self.label_2.setFont(font)
        self.label_2.setWordWrap(True)
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(self.OutputWindow)
        self.label_3.setGeometry(QtCore.QRect(40, 170, 341, 31))
        font = QtGui.QFont()
        font.setFamily("Chandas")
        font.setPointSize(16)
        font.setBold(True)
        font.setWeight(75)
        self.label_3.setFont(font)
        self.label_3.setWordWrap(True)
        self.label_3.setObjectName("label_3")
        self.label_4 = QtWidgets.QLabel(self.OutputWindow)
        self.label_4.setGeometry(QtCore.QRect(40, 220, 971, 331))
        font = QtGui.QFont()
        font.setFamily("Chandas")
        font.setPointSize(10)
        self.label_4.setFont(font)
        self.label_4.setWordWrap(True)
        self.label_4.setObjectName("label_4")
        self.pushButton = QtWidgets.QPushButton(self.OutputWindow)
        self.pushButton.setGeometry(QtCore.QRect(920, 150, 89, 25))
        font = QtGui.QFont()
        font.setFamily("Chandas")
        font.setPointSize(10)
        self.pushButton.setFont(font)
        self.pushButton.setObjectName("pushButton")
        self.pushButton_2 = QtWidgets.QPushButton(self.OutputWindow)
        self.pushButton_2.setGeometry(QtCore.QRect(930, 510, 89, 25))
        font = QtGui.QFont()
        font.setFamily("Chandas")
        font.setPointSize(10)
        self.pushButton_2.setFont(font)
        self.pushButton_2.setObjectName("pushButton_2")
        SecondWindow.setCentralWidget(self.OutputWindow)
        self.menubar = QtWidgets.QMenuBar(SecondWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1059, 22))
        self.menubar.setObjectName("menubar")
        self.menuPage_1 = QtWidgets.QMenu(self.menubar)
        self.menuPage_1.setObjectName("menuPage_1")
        SecondWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(SecondWindow)
        self.statusbar.setObjectName("statusbar")
        SecondWindow.setStatusBar(self.statusbar)
        self.actionCoverage = QtWidgets.QAction(SecondWindow)
        self.actionCoverage.setObjectName("actionCoverage")
        self.actionHinges_and_Flexibility = QtWidgets.QAction(SecondWindow)
        self.actionHinges_and_Flexibility.setObjectName("actionHinges_and_Flexibility")
        self.actionComposite_and_Topology_File = QtWidgets.QAction(SecondWindow)
        self.actionComposite_and_Topology_File.setObjectName("actionComposite_and_Topology_File")
        self.actionCustom_hinges = QtWidgets.QAction(SecondWindow)
        self.actionCustom_hinges.setObjectName("actionCustom_hinges")
        self.menuPage_1.addAction(self.actionCoverage)
        self.menuPage_1.addAction(self.actionHinges_and_Flexibility)
        self.menuPage_1.addAction(self.actionComposite_and_Topology_File)
        self.menuPage_1.addAction(self.actionCustom_hinges)
        self.menubar.addAction(self.menuPage_1.menuAction())

        self.pushButton_3 = QtWidgets.QPushButton(self.OutputWindow)
        self.pushButton_3.setFont(font)
        self.pushButton_3.setObjectName("pushButton_3")
        self.pushButton_3.setGeometry(QtCore.QRect(810, 540, 89, 25))

        self.retranslateUi(SecondWindow)
        QtCore.QMetaObject.connectSlotsByName(SecondWindow)

    def retranslateUi(self, SecondWindow):
        _translate = QtCore.QCoreApplication.translate
        SecondWindow.setWindowTitle(_translate("SecondWindow", "SecondWindow"))
        self.label.setText(_translate("SecondWindow",
                                      "<html><head/><body><p align=\"justify\">In this graph,you can see which parts of the reference FASTA sequence are covered by structure. This structures come from either the <a href=\"https://www.rcsb.org/\"><span style=\" text-decoration: underline; color:#0000ff;\">Protein Data Bank</span></a>, <a href=\"https://www.deepmind.com/blog/alphafold-a-solution-to-a-50-year-old-grand-challenge-in-biology\"><span style=\" text-decoration: underline; color:#0000ff;\">AlphaFold</span></a> models or <a href=\"https://www.ipd.uw.edu/2021/07/rosettafold-accurate-protein-structure-prediction-accessible-to-all/\"><span style=\" text-decoration: underline; color:#0000ff;\">RoseTTaFold</span></a> models.<br/></p></body></html>"))
        self.label_2.setText(_translate("SecondWindow", "COVERAGE"))
        self.label_3.setText(_translate("SecondWindow", "HINGES AND FLEXIBILITY"))
        self.label_4.setText(_translate("SecondWindow",
                                        "<html><head/><body><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\">Flexibility is an important feature of proteins, since they need to move to perform their function and interact with their substrates. In the following section, we provide you with two types of flexibility prediction: the Dynamic Flexibility Index and Hinge Prediction. The overlap of these two measures might be helpful for you, in case you wanted to modify the final topology file with some of these hinges.</span></p><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; font-style:italic; color:#323232;\">Dynamic Flexibility Index</span><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\"><br/>This is per-residue index indicating the contribution of each residue to the overall flexibility of the protein. It uses a method based in an Elastic Network Model (ENM), which is a more lightweight (but less precise, obviously) alternative to Molecular Dynamics. for more info, </span><a href=\"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3673471/\"><span style=\" text-decoration: underline; color:#0000ff;\">here</span></a><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\"> is the original paper.</span></p><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; font-style:italic; color:#323232;\">Hinge Prediction</span><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\"><br/>Hinges are the regions of the protein that allow it to move and change conformations. Using </span><a href=\"https://academic.oup.com/bioinformaticsadvances/article/2/1/vbac007/6525212?login=true\"><span style=\" text-decoration: underline; color:#0000ff;\">this tool</span></a><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\"> we provide you with some suggested hinge regions. Note that this information is only available for experimental structures. This is due to the use of ENM, it is not designed to work with predicted models that might contain important artifacts, and, in this case, that are split into the highly confidently predicted regions.</span></p><p><a name=\"modebar-338c4d\"/><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\"><br/></span></p></body></html>"))
        self.pushButton.setText(_translate("SecondWindow", "Show plot"))
        self.pushButton.clicked.connect(self.show_plot)

        self.pushButton_2.setText(_translate("SecondWindow", "Show plot"))
        self.pushButton_2.clicked.connect(self.show_Secondplot)

        self.pushButton_3.setText(_translate("SecondWindow", "Next"))
        self.pushButton_3.clicked.connect(self.openThirdWindow)
        self.pushButton_3.clicked.connect(SecondWindow.hide)

        self.menuPage_1.setTitle(_translate("SecondWindow", ""))
        self.actionCoverage.setText(_translate("SecondWindow", "Coverage"))
        self.actionHinges_and_Flexibility.setText(_translate("SecondWindow", "Hinges and Flexibility"))
        self.actionComposite_and_Topology_File.setText(_translate("SecondWindow", "Composite and Topology File"))
        self.actionCustom_hinges.setText(_translate("SecondWindow", "Custom hinges"))


class Ui_ThirdWindow(object):
    def update_dropdown(self):
        structure_list = []
        for child in Path(os.path.join(output_dir, query_name, "REPORT", "COVERAGE")).iterdir():
            if child.is_file() and "composite" not in str(child):
                structure_list.append(str(child))

        comp_dict = read_compsite_files(os.path.join(output_dir, query_name, "REPORT", "COVERAGE"))

        split_path = os.path.join(output_dir, query_name).split("/")
        out_name = split_path[-1]

        filename = os.path.join(output_dir, query_name, "REPORT", "COVERAGE", f"{out_name}_composite_coverage.csv")
        df = pd.read_csv(filename)
        initial_structure_list = []
        for str1 in df.columns[1:]:
            id1 = str(os.path.basename(str1))
            for str2 in structure_list:
                id2 = str(os.path.basename(str2))
                print(f"COMPARING: {id1} and {id2}")
                if id1[0:-4] == id2[0:-13]:
                    initial_structure_list.append(str2)
                    print(f"ADDED {str2}")

        return structure_list, initial_structure_list

    def update_Thirdgraph(self):
        options_chosen = [file for file in Path(os.path.join(output_dir, query_name, "REPORT", "COVERAGE")).iterdir()]

        if len(options_chosen) == 0:
            return None

        if options_chosen is None:
            fig = {}
            return fig

        i = 0
        df_list = []
        structure_list = []
        for file in options_chosen:
            if "composite" not in str(file):
                i += 1
                print(file)
                df = pd.read_csv(file)
                df_list.append(df)
                structure_list.append(file)

        fig4 = make_subplots(rows=i, cols=1, shared_xaxes=True, x_title="Residue position")

        i = 1
        for df in df_list:
            fig4.append_trace(go.Scatter(
                x=df[df.columns[0]],  # ResID
                y=df[df.columns[1]],
                fill='tozeroy',
                name=str(structure_list[i - 1])
            ), row=i, col=1)
            i += 1

        # fig4.update_layout(height=400, width=1000, title_text="Coverage")
        fig4.update_layout(title_text="Coverage")
        fig4.update_yaxes(showgrid=False, range=[0, 1], showticklabels=False)

        png4 = pio.write_image(fig4, f"{plot_dir}/structure_plot.png", scale=1, width=1400, height=700)
        png4_html = pio.write_html(fig4, f"{plot_dir}/structure_plot.html")


        return fig4

    def show_Compositeplot(self):
        '''
        Connection function: show coverage plot
        Calls update_graph() at the beginning in order to save the image and then show it
        '''
        self.update_dropdown()
        self.update_Thirdgraph()
        self.plot_out = QtWidgets.QMainWindow()
        self.pl = Ui_Thirdplot()
        self.pl.setup(self.plot_out)
        self.plot_out.show()

    def close_third(self, main_w):
        main_w.hide()
        self.returnsecond = QtWidgets.QMainWindow()
        self.s = Ui_SecondWindow()
        self.s.setupOutput(self.returnsecond)
        self.returnsecond.show()

    def search(self):
        '''
        Retrieve structures we want to include in our composite
        '''
        global textboxValue
        global textbox
        textboxValue = self.textEdit.toPlainText()
        textbox = self.textEdit_2.toPlainText()
        if textboxValue != '' or textbox != '':
            self.onclick_topology()
        else:
            self.warning_message = QtWidgets.QMainWindow()
            self.warning_message.resize(488, 272)
            self.label_message = QtWidgets.QLabel(self.warning_message)
            self.label_message.setGeometry(QtCore.QRect(40, 90, 401, 61))
            font = QtGui.QFont()
            font.setFamily("Chandas")
            font.setPointSize(14)
            font.setBold(True)
            font.setWeight(75)
            self.label_message.setFont(font)
            self.label_message.setText("WARNING: NO HINGES INTRODUCED")
            self.warning_message.show()


    def onclick_topology(self):
        '''
        Retrieve fragments and hinges selected and create custom topology file
        '''
        str_hinges_input = textbox
        selected_fragments = textboxValue
        output_directory = output_dir
        structure_list = []
        print(str_hinges_input)
        print(selected_fragments)
        try:
            for child in Path(os.path.join(output_dir, query_name, "PDB", "total")).iterdir():
                if child.is_file() and "composite" not in str(child):
                    if os.path.basename(child)[0:-4] == selected_fragments:
                        structure_list.append(child)
        except:
            pass

        try:
            for child in Path(os.path.join(output_dir, query_name, "PDB", "partial")).iterdir():
                if child.is_file() and "composite" not in str(child):
                    if os.path.basename(child)[0:-4] == selected_fragments:
                        structure_list.append(child)


        except:
            pass

        try:
            for child in Path(os.path.join(output_dir, query_name, "PDB", "CHAINS")).iterdir():
                if child.is_file() and "composite" not in str(child):
                    if os.path.basename(child)[0:-4] == selected_fragments:
                        structure_list.append(child)
        except:
            pass

        try:
            for child in Path(os.path.join(output_dir, query_name, "ALPHAFOLD", "DOMAINS")).iterdir():
                if child.is_file() and "confident" not in str(child) and "domains" not in str(child):
                    if str(os.path.basename(child)[0:-4]) == str(selected_fragments):
                        structure_list.append(child)
        except:
            pass

        try:
            for child in Path(os.path.join(output_dir, query_name, "ROSETTAFOLD", "DOMAINS")).iterdir():
                if child.is_file() and "confident" not in str(child) and "domains" not in str(child):
                    if str(os.path.basename(child)[0:-4]) == str(selected_fragments):
                        structure_list.append(child)
        except:
            pass

        fasta = "input_fasta/" + str(os.path.basename(query_name)) + ".fasta"

        rigid_bodies = make_rb_list(structure_list, fasta)

        ## incorporate the hinges
        hinges_list = [hinge for hinge in str_hinges_input.split(",")]
        hinges_list = [tuple(i.split(':')) for i in hinges_list]
        hinges_list = [(int(i[0]), int(i[1])) for i in hinges_list]

        final_rigid_bodies = []

        for rb in rigid_bodies:
            print(f"RB: {rb.pdb_fn}")
            split_rb = rb.split_rb_hinges(hinges_list)
            print(f"RB SPLIT: {[rb.residue_range for rb in split_rb]}")
            final_rigid_bodies = final_rigid_bodies + split_rb

        print(f"INITIAL = {len(rigid_bodies)}, FINAL = {len(final_rigid_bodies)}")

        final_rigid_bodies.sort(key=lambda x: x.residue_range[0])
        str_out = str(os.path.join(output_dir, query_name))
        out_name = str_out.split("/")[-1]
        # Write the topology file
        write_custom_topology(os.path.join(output_dir, query_name, "IMP", f"{out_name}_custom.topology"),
                              final_rigid_bodies)

        return f"Topology file created with:{[str(rb.pdb_fn) for rb in final_rigid_bodies]}"

    def setupOutput(self, ThirdWindow):
        '''
        set up second window, when introduced input fasta and output directory, window that will be opened to show
        results
        '''
        ThirdWindow.setObjectName("ThirdWindow")
        ThirdWindow.resize(1066, 686)
        ThirdWindow.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.SecondOutputWindow = QtWidgets.QWidget(ThirdWindow)
        self.SecondOutputWindow.setObjectName("SecondOutputWindow")
        self.label = QtWidgets.QLabel(self.SecondOutputWindow)
        self.label.setGeometry(QtCore.QRect(40, 120, 971, 131))
        font = QtGui.QFont()
        font.setFamily("Chandas")
        font.setPointSize(10)
        self.label.setFont(font)
        self.label.setWordWrap(True)
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(self.SecondOutputWindow)
        self.label_2.setGeometry(QtCore.QRect(40, 80, 431, 31))
        font = QtGui.QFont()
        font.setFamily("Chandas")
        font.setPointSize(16)
        font.setBold(True)
        font.setWeight(75)
        self.label_2.setFont(font)
        self.label_2.setWordWrap(True)
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(self.SecondOutputWindow)
        self.label_3.setGeometry(QtCore.QRect(40, 290, 341, 31))
        font = QtGui.QFont()
        font.setFamily("Chandas")
        font.setPointSize(16)
        font.setBold(True)
        font.setWeight(75)
        self.label_3.setFont(font)
        self.label_3.setWordWrap(True)
        self.label_3.setObjectName("label_3")
        self.label_4 = QtWidgets.QLabel(self.SecondOutputWindow)
        self.label_4.setGeometry(QtCore.QRect(40, 330, 971, 191))
        font = QtGui.QFont()
        font.setFamily("Chandas")
        font.setPointSize(10)
        self.label_4.setFont(font)
        self.label_4.setWordWrap(True)
        self.label_4.setObjectName("label_4")
        self.pushButton = QtWidgets.QPushButton(self.SecondOutputWindow)
        self.pushButton.setGeometry(QtCore.QRect(930, 270, 89, 25))
        font = QtGui.QFont()
        font.setFamily("Chandas")
        font.setPointSize(10)
        self.pushButton.setFont(font)
        self.pushButton.setObjectName("pushButton")
        self.pushButton_2 = QtWidgets.QPushButton(self.SecondOutputWindow)
        self.pushButton_2.setGeometry(QtCore.QRect(810, 540, 211, 25))
        font = QtGui.QFont()
        font.setFamily("Chandas")
        font.setPointSize(10)
        self.pushButton_2.setFont(font)
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton_3 = QtWidgets.QPushButton(self.SecondOutputWindow, clicked = lambda: self.close_third(ThirdWindow))
        self.pushButton_3.setGeometry(QtCore.QRect(888, 600, 101, 31))
        font = QtGui.QFont()
        font.setFamily("Chandas")
        font.setPointSize(10)
        self.pushButton_3.setFont(font)
        self.pushButton_3.setObjectName("pushButton")
        self.textEdit = QtWidgets.QTextEdit(self.SecondOutputWindow)
        self.textEdit.setGeometry(QtCore.QRect(250, 500, 201, 31))
        self.textEdit.setObjectName("textEdit")
        self.textEdit_2 = QtWidgets.QTextEdit(self.SecondOutputWindow)
        self.textEdit_2.setGeometry(QtCore.QRect(250, 560, 201, 31))
        self.textEdit_2.setObjectName("textEdit_2")
        font2 = QtGui.QFont()
        font2.setFamily("Open Sans")
        font2.setPointSize(10)
        self.label_5 = QtWidgets.QLabel(self.SecondOutputWindow)
        self.label_5.setGeometry(QtCore.QRect(90, 560, 151, 31))
        self.label_5.setFont(font2)
        self.label_5.setObjectName("label")
        self.label_6 = QtWidgets.QLabel(self.SecondOutputWindow)
        self.label_6.setGeometry(QtCore.QRect(90, 500, 151, 31))
        self.label_6.setFont(font2)
        self.label_6.setObjectName("label")
        ThirdWindow.setCentralWidget(self.SecondOutputWindow)
        self.menubar = QtWidgets.QMenuBar(ThirdWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1059, 22))
        self.menubar.setObjectName("menubar")
        self.menuPage_1 = QtWidgets.QMenu(self.menubar)
        self.menuPage_1.setObjectName("menuPage_1")
        ThirdWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(ThirdWindow)
        self.statusbar.setObjectName("statusbar")
        ThirdWindow.setStatusBar(self.statusbar)
        self.actionCoverage = QtWidgets.QAction(ThirdWindow)
        self.actionCoverage.setObjectName("actionCoverage")
        self.actionHinges_and_Flexibility = QtWidgets.QAction(ThirdWindow)
        self.actionHinges_and_Flexibility.setObjectName("actionHinges_and_Flexibility")
        self.actionComposite_and_Topology_File = QtWidgets.QAction(ThirdWindow)
        self.actionComposite_and_Topology_File.setObjectName("actionComposite_and_Topology_File")
        self.actionCustom_hinges = QtWidgets.QAction(ThirdWindow)
        self.actionCustom_hinges.setObjectName("actionCustom_hinges")
        self.menuPage_1.addAction(self.actionCoverage)
        self.menuPage_1.addAction(self.actionHinges_and_Flexibility)
        self.menuPage_1.addAction(self.actionComposite_and_Topology_File)
        self.menuPage_1.addAction(self.actionCustom_hinges)
        self.menubar.addAction(self.menuPage_1.menuAction())

        self.retranslateUi(ThirdWindow)
        QtCore.QMetaObject.connectSlotsByName(ThirdWindow)


    def retranslateUi(self, ThirdWindow):
        _translate = QtCore.QCoreApplication.translate
        ThirdWindow.setWindowTitle(_translate("ThirdWindow", "ThirdWindow"))
        self.label.setText(_translate("ThirdWindow",
                                      "<html><head/><body><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\">From all the structures retrieved by the program and provided by the user, the program generates this composite, trying to cover as much of the reference sequence as possible, avoiding overlaps.</span></p><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\">This composite can be used to automatically build an </span><a href=\"https://integrativemodeling.org/2.5.0/doc/ref/classIMP_1_1pmi_1_1topology_1_1TopologyReader.html\"><span style=\" text-decoration: underline; color:#0000ff;\">IMP topology file</span></a></p><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\">If you want to generate another topology file, with your custom choice of fragments, simply select the fragments you want to include and click the button &quot;Create Topology File&quot;. it will save yhe output inj the folder IMP of your selected output folder.</span></p></body></html>"))
        self.label_2.setText(_translate("ThirdWindow", "COMPOSITE AND TOPOLOGY FILE"))
        self.label_3.setText(_translate("ThirdWindow", "CUSTOM HINGES"))
        self.label_4.setText(_translate("ThirdWindow",
                                        "<html><head/><body><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\">In this section you can introduce hinge regions. Hinge regions are those regions of the protein that bend, allowing the movement of the more rigid domains, which is essential for the interaciton of the proteins with other biomolecules.</span></p><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; font-weight:600; color:#323232;\">How are hinges encoded?</span><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\"><br/>Let\'s imagine you have a proterin of 200 amino acids. The DFI and PACKMAN hinge prediction are indicating a putative hinge region between positions 50 and 100 and another one between 120 and 130.</span></p><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\">In the box, you will need to introduce the hinges with the following format: 50:100,120:130</span></p><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\">The program will split the structures, in the topology file, according to the hinges introduced.</span></p><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;"))
        self.label_5.setText(_translate("MainWindow", "Selected fragments:"))
        self.label_6.setText(_translate("MainWindow", "Selected structures:"))
        self.pushButton.setText(_translate("ThirdWindow", "Show plot"))
        self.pushButton.clicked.connect(self.show_Compositeplot)

        self.pushButton_2.setText(_translate("ThirdWindow", "Generate IMP Topology File"))
        self.pushButton_2.clicked.connect(self.search)
        self.pushButton_3.setText(_translate("MainWindow", "Return"))


        self.menuPage_1.setTitle(_translate("ThirdWindow", ""))
        self.actionCoverage.setText(_translate("ThirdWindow", "Coverage"))
        self.actionHinges_and_Flexibility.setText(_translate("ThirdWindow", "Hinges and Flexibility"))
        self.actionComposite_and_Topology_File.setText(_translate("ThirdWindow", "Composite and Topology File"))
        self.actionCustom_hinges.setText(_translate("ThirdWindow", "Custom hinges"))


class Ui_plot(object):
    '''
    set up window to show image of plot
    '''

    def setup(self, plotWindow):
        plotWindow.resize(1400, 733)
        self.plotwindow = QtWidgets.QWidget(plotWindow)
        self.plotwindow.setObjectName("plotwindow")
        self.m_output = QtWebEngineWidgets.QWebEngineView(self.plotwindow)
        self.m_output.setGeometry(0, -10, 1391, 721)
        url = QtCore.QUrl.fromLocalFile(f"{plot_dir}/coverage_plot.html")
        self.m_output.load(url)
        plotWindow.setCentralWidget(self.plotwindow)
        self.menubar = QtWidgets.QMenuBar(plotWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1066, 22))
        self.menubar.setObjectName("menubar")
        plotWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(plotWindow)
        self.statusbar.setObjectName("statusbar")
        plotWindow.setStatusBar(self.statusbar)

        self.retranslateUi(plotWindow)
        QtCore.QMetaObject.connectSlotsByName(plotWindow)

    def retranslateUi(self, plotWindow):
        _translate = QtCore.QCoreApplication.translate
        plotWindow.setWindowTitle(_translate("plotWindow", "plotWindow"))


class Ui_Secondplot(object):
    '''
    set up another window to show predicted hinges and flexibility plot
    '''
    def setup(self, plotWindow):
        plotWindow.resize(1400, 733)
        self.plotwindow = QtWidgets.QWidget(plotWindow)
        self.plotwindow.setObjectName("plotwindow")
        self.s_output = QtWebEngineWidgets.QWebEngineView(self.plotwindow)
        self.s_output.setGeometry(0, -10, 1391, 721)
        url = QtCore.QUrl.fromLocalFile(f"{plot_dir}/hinges_prediction.html")
        self.s_output.load(url)
        plotWindow.setCentralWidget(self.plotwindow)
        self.menubar = QtWidgets.QMenuBar(plotWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1066, 22))
        self.menubar.setObjectName("menubar")
        plotWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(plotWindow)
        self.statusbar.setObjectName("statusbar")
        plotWindow.setStatusBar(self.statusbar)

        self.retranslateUi(plotWindow)
        QtCore.QMetaObject.connectSlotsByName(plotWindow)

    def retranslateUi(self, plotWindow):
        _translate = QtCore.QCoreApplication.translate
        plotWindow.setWindowTitle(_translate("plotWindow", "plotWindow"))


class Ui_Thirdplot(object):
    '''
    set up another window to show predicted hinges and flexibility plot
    '''
    def setup(self, plotWindow):
        plotWindow.resize(1400, 733)
        self.plotwindow = QtWidgets.QWidget(plotWindow)
        self.plotwindow.setObjectName("plotwindow")
        self.m_output = QtWebEngineWidgets.QWebEngineView(self.plotwindow)
        self.m_output.setGeometry(0, -10, 1391, 721)
        url = QtCore.QUrl.fromLocalFile(f"{plot_dir}/structure_plot.html")
        self.m_output.load(url)
        plotWindow.setCentralWidget(self.plotwindow)
        self.menubar = QtWidgets.QMenuBar(plotWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1066, 22))
        self.menubar.setObjectName("menubar")
        plotWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(plotWindow)
        self.statusbar.setObjectName("statusbar")
        plotWindow.setStatusBar(self.statusbar)

        self.retranslateUi(plotWindow)
        QtCore.QMetaObject.connectSlotsByName(plotWindow)

    def retranslateUi(self, plotWindow):
        _translate = QtCore.QCoreApplication.translate
        plotWindow.setWindowTitle(_translate("plotWindow", "plotWindow"))







if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    # GUI pops up
    MainWindow.show()
    #loop = QEventLoop()
    #loop.exec()  # waits

    sys.exit(app.exec_())

