import plotly
from PyQt5 import QtCore, QtGui, QtWidgets, QtWebEngineWidgets
from PyQt5.QtCore import QDir, QEventLoop, Qt, QUrl
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import QFileDialog, QTextEdit, QAction, QDialog, QSizePolicy, QGridLayout, QSpacerItem, \
    QScrollArea, QWidget, QVBoxLayout, QLabel, QSizeGrip, QPlainTextEdit, QPushButton

##########################
# Ferran's script functions
##########################

import Bio.SeqIO as IO
# from custom_parser import parser
from PySide2.QtWebEngineWidgets import QWebEngineView
from matplotlib.backends.backend_template import FigureCanvas
from plotly.subplots import make_subplots
from pygments.lexers import go

import SecondWindow
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



class Ui_MainWindow(object):
    def openWindow(self):
        '''
        Creates connection between input window and output window (second window)
        '''
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_SecondWindow()
        self.ui.setupOutput(self.window)
        self.window.show()

        #view.show()

    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(800, 600)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")

        self.formLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.formLayoutWidget.setGeometry(QtCore.QRect(0, 0, 2, 2))
        self.formLayoutWidget.setObjectName("formLayoutWidget")
        self.formLayout = QtWidgets.QFormLayout(self.formLayoutWidget)
        self.formLayout.setContentsMargins(0, 0, 0, 0)
        self.formLayout.setObjectName("formLayout")

        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(40, 40, 461, 41))
        self.label.setObjectName("label")
        self.pushButton = QtWidgets.QPushButton(self.centralwidget) # START program
        self.pushButton.setGeometry(QtCore.QRect(660, 40, 71, 61))
        self.pushButton.setObjectName("pushButton")
        self.pushButton_2 = QtWidgets.QPushButton(self.centralwidget) # Select fasta button
        self.pushButton_2.setGeometry(QtCore.QRect(550, 40, 71, 61))
        self.pushButton_2.setObjectName("pushButton_2")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(40, 120, 161, 31))
        self.label_2.setObjectName("label_2")


        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(40, 240, 181, 31))
        self.label_3.setObjectName("label_3")
        self.label_4 = QtWidgets.QLabel(self.centralwidget)
        self.label_4.setGeometry(QtCore.QRect(420, 120, 171, 31))
        self.label_4.setObjectName("label_4")
        self.label_5 = QtWidgets.QLabel(self.centralwidget)
        self.label_5.setGeometry(QtCore.QRect(420, 170, 161, 31))
        self.label_5.setObjectName("label_5")
        self.label_6 = QtWidgets.QLabel(self.centralwidget)
        self.label_6.setGeometry(QtCore.QRect(40, 180, 121, 31))
        self.label_6.setObjectName("label_6")
        self.label_7 = QtWidgets.QLabel(self.centralwidget)
        self.label_7.setGeometry(QtCore.QRect(40, 300, 121, 31))
        self.label_7.setObjectName("label_7")
        self.label_8 = QtWidgets.QLabel(self.centralwidget)
        self.label_8.setGeometry(QtCore.QRect(420, 230, 121, 31))
        self.label_8.setObjectName("label_8")
        self.lineEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit.setGeometry(QtCore.QRect(220, 120, 113, 25))
        self.lineEdit.setObjectName("lineEdit")
        self.lineEdit_2 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_2.setGeometry(QtCore.QRect(490, 230, 113, 25))
        self.lineEdit_2.setObjectName("lineEdit_2")
        self.lineEdit_3 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_3.setGeometry(QtCore.QRect(160, 180, 113, 25))
        self.lineEdit_3.setObjectName("lineEdit_3")
        self.lineEdit_4 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_4.setGeometry(QtCore.QRect(230, 240, 113, 25))
        self.lineEdit_4.setObjectName("lineEdit_4")
        self.lineEdit_5 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_5.setGeometry(QtCore.QRect(170, 300, 113, 25))
        self.lineEdit_5.setObjectName("lineEdit_5")
        self.lineEdit_6 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_6.setGeometry(QtCore.QRect(560, 180, 113, 25))
        self.lineEdit_6.setObjectName("lineEdit_6")
        self.lineEdit_7 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_7.setGeometry(QtCore.QRect(600, 120, 113, 25))
        self.lineEdit_7.setObjectName("lineEdit_7")
        self.lineEdit_8 = QtWidgets.QLineEdit(self.centralwidget) # Line where fasta path is printed
        self.lineEdit_8.setGeometry(QtCore.QRect(30, 44, 491, 31))
        self.lineEdit_8.setObjectName("lineEdit")
        self.toolButton = QtWidgets.QToolButton(self.centralwidget)
        self.toolButton.setGeometry(QtCore.QRect(310, 120, 26, 24))
        self.toolButton.setObjectName("toolButton")

        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 22))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)


    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.lineEdit_8.setText(_translate("MainWindow", "Select sequence file and output directory"))
        # self.label.setText(_translate("MainWindow", "Select sequence file and output directory"))
        self.pushButton.setText(_translate("MainWindow", "Select File"))
        self.label_2.setText(_translate("MainWindow", "Output directory"))
        self.label_3.setText(_translate("MainWindow", "RosettaFold (PDB format)*"))
        self.label_4.setText(_translate("MainWindow", "AlphaFold PAE JSON file*"))
        self.label_5.setText(_translate("MainWindow", "Custom templates*"))
        self.label_6.setText(_translate("MainWindow", "Run AlphaFold*"))
        self.label_7.setText(_translate("MainWindow", "Run RosettaFold*"))
        self.label_8.setText(_translate("MainWindow", "AlphaFold (PDB format)*"))
        self.toolButton.setText(_translate("MainWindow", "..."))
        self.toolButton.clicked.connect(self.create_dir)


        self.pushButton.setText(_translate("MainWindow", "START"))
        self.pushButton.clicked.connect(self.run)
        # click button and call function
        self.pushButton_2.setText(_translate("MainWindow", "Select file"))
        self.pushButton_2.clicked.connect(self.get_file)


    def get_file(self):
        '''
        Function that gets filename of fasta, its directory, and the sequence it self
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
        #query_name = print(list(file))
        query_name = query_name.replace(">", "").strip()

        return self.lineEdit_8.setText(source), filename, query_name, fasta

    def create_dir(self):
        '''
        Creates output directory
        '''
        global output_dir
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        name, _ = QFileDialog.getSaveFileName(None, 'QFileDialog.getSaveFileName()', "")
        output_dir = name

        return self.lineEdit.setText(output_dir)


    def run(self):
        # Set fasta path and outdir and wait until those variables are set
        verbose = True


        ### Initializing the LOG system ###

        logdir = os.path.join(output_dir, query_name, "LOG", "")
        Path(logdir).mkdir(parents=True, exist_ok=True)

        l.basicConfig(format="%(levelname)s:%(message)s",
                      filename=os.path.join(logdir, f"{query_name}.log"),
                      level=l.DEBUG)

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
        l.info(f"{pdb_dir}, {fasta_dir}, {report_dir}, {hinges_dir}, {IMP_dir} .")

        Path(blast_dir).mkdir(parents=True, exist_ok=True)
        Path(pdb_dir).mkdir(parents=True, exist_ok=True)
        Path(fasta_dir).mkdir(parents=True, exist_ok=True)
        Path(report_dir).mkdir(parents=True, exist_ok=True)
        Path(hinges_dir).mkdir(parents=True, exist_ok=True)
        Path(IMP_dir).mkdir(parents=True, exist_ok=True)

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

        '''
        # Don't forget the user's templates!
        l.info("Checking user's templates")
        if parser.custom_templates:
            structures_for_query.append(parser.custom_templates)


        ### ALPHAFOLD & PAE

        # If you want to use your AF model and PAE file:
        if parser.alphamodel:
            l.info(f"custom AF model detected: {parser.alphamodel}")
            af_dir = os.path.join(output_dir, query_name, "ALPHAFOLD", "")
            Path(af_dir).mkdir(parents=True, exist_ok=True)
            AF_server_model = parser.alphamodel
            shutil.copy(AF_server_model, os.path.join(af_dir,
                                                      PurePosixPath(AF_server_model).name))
        if parser.PAE_json:
            l.info(f"custom PAE matrix detected: {parser.PAE_json}")
            PAE_dir = os.path.join(af_dir, "PAE", "")
            Path(PAE_dir).mkdir(parents=True, exist_ok=True)
            PAE_json = parser.PAE_json

        from bin.utilities import submit_AF_to_SLURM, submit_RF_to_SLURM

        if parser.run_alphafold:
            l.info(f"Submitting AF2 SLURM batch script")
            af_dir = os.path.join(output_dir, query_name, "ALPHAFOLD", "")
            Path(af_dir).mkdir(parents=True, exist_ok=True)
            submit_AF_to_SLURM(fasta, af_dir, workload_manager="sbatch",
                               dummy_dir=".", max_jobs_in_queue=None)

        ### ROSETTAFOLD

        # If you want to use your RF model:
        if parser.rosettamodel:
            l.info(f"custom RF model detected: {parser.rosettamodel}")
            rf_dir = os.path.join(output_dir, query_name, "ROSETTAFOLD", "")
            custom_rf_dir = os.path.join(rf_dir, "CUSTOM")
            Path(rf_dir).mkdir(parents=True, exist_ok=True)
            Path(custom_rf_dir).mkdir(parents=True, exist_ok=True)
            RF_custom_model = parser.rosettamodel
            shutil.copy(RF_custom_model, os.path.join(custom_rf_dir, PurePosixPath(RF_custom_model).name))
        if parser.run_rosettafold:
            l.info(f"Submitting RF SLURM batch script")
            rf_dir = os.path.join(output_dir, query_name, "ROSETTAFOLD", "")
            Path(rf_dir).mkdir(parents=True, exist_ok=True)
            submit_RF_to_SLURM(fasta, rf_dir, workload_manager="sbatch", dummy_dir=".", max_jobs_in_queue=None)

        ### Extract confident regions

        if (parser.alphamodel and parser.PAE_json) or (parser.run_alphafold):
            # Setting up the parameters for the PHENIX library
            master_phil = iotbx.phil.parse(master_phil_str)
            params = master_phil.extract()
            master_phil.format(python_object=params).show(out=sys.stdout)
            p = params.process_predicted_model
            p.domain_size = cfg.CCTBXconfig["AF2_domain_size"]
            p.maximum_rmsd = cfg.CCTBXconfig["AF2_maximum_rmsd"]
            p.b_value_field_is = 'lddt'

            from iotbx.data_manager import DataManager

            dm = DataManager()
            dm.set_overwrite(True)

            l.info("Extracting AF2 high confidence domains")
            domains_dir = os.path.join(af_dir, "DOMAINS", "")
            Path(domains_dir).mkdir(parents=True, exist_ok=True)
            l.info(f"Domains will be stored in:{domains_dir}")

            af_conficent_regions = []
            for filename in os.listdir(af_dir):
                if os.path.isfile(os.path.join(af_dir, filename)):
                    l.info(f"Processing file: {filename}")
                    print("\nProcessing and splitting model into domains")

                    m = dm.get_model(os.path.join(af_dir, filename))
                    pae_matrix = pae_matrix = parse_json_PAE(PAE_json)
                    model_info = process_predicted_model(m, params, pae_matrix)

                    chainid_list = model_info.chainid_list
                    print("Segments found: %s" % (" ".join(chainid_list)))

                    mmm = model_info.model.as_map_model_manager()

                    # Write all the domains in one file
                    mmm.write_model(os.path.join(domains_dir,
                                                 f"{PurePosixPath(filename).stem}_domains.pdb"))

                    # Write different domains in different files
                    for chainid in chainid_list:
                        selection_string = "chain %s" % (chainid)
                        ph = model_info.model.get_hierarchy()
                        asc1 = ph.atom_selection_cache()
                        sel = asc1.selection(selection_string)
                        m1 = model_info.model.select(sel)
                        filepath = os.path.join(domains_dir,
                                                f"{PurePosixPath(filename).stem}_{chainid}_AF.pdb")
                        dm.write_model_file(m1, filepath)
                        structures_for_query.append(filepath)

                    conf_domains = extract_residue_list(os.path.join(domains_dir,
                                                                     f"{PurePosixPath(filename).stem}_domains.pdb"),
                                                        domains_dir)
                    l.info(f"Residue list of confident domains: {conf_domains}")

        # For Rosettafold models
        if os.path.exists(os.path.join(output_dir, query_name, "ROSETTAFOLD", "")):
            rf_dir = os.path.join(output_dir, query_name, "ROSETTAFOLD", "")
            rf_models_dir = os.path.join(rf_dir, "model", "")
            # Setting up the parameters for the PHENIX library
            master_phil = iotbx.phil.parse(master_phil_str)
            params = master_phil.extract()
            master_phil.format(python_object=params).show(out=sys.stdout)
            p = params.process_predicted_model
            p.domain_size = cfg.CCTBXconfig["RF_domain_size"]
            p.maximum_rmsd = cfg.CCTBXconfig["RF_maximum_rmsd"]
            p.b_value_field_is = 'rmsd'

            from iotbx.data_manager import DataManager

            dm = DataManager()
            dm.set_overwrite(True)

            l.info("Extracting RoseTTaFold high confidence domains")
            domains_dir = os.path.join(rf_dir, "DOMAINS", "")
            Path(domains_dir).mkdir(parents=True, exist_ok=True)
            l.info(f"Domains will be stored in:{domains_dir}")
            abs_rf_dir = os.path.abspath(rf_models_dir)
            abs_custom_rf_dir = os.path.abspath(rf_models_dir)
            rf_conficent_regions = []
            for filename in os.listdir(abs_rf_dir):
                if os.path.isfile(os.path.join(abs_rf_dir, filename)) and \
                        (fnmatch.fnmatch(filename, "model_*.crderr.pdb") or \
                         fnmatch.fnmatch(filename, "model_*-crderr.pdb")):
                    newname = filename.split(".")
                    l.info(f"LIST NEWNAME: {newname}")
                    noext = newname[0:-1]
                    noext = "-".join(noext)
                    ext = newname[-1]
                    newname = noext + "." + ext
                    l.info(f"NEWNAME: {newname}")
                    filename = os.path.join(abs_rf_dir, filename)
                    newname = os.path.join(abs_rf_dir, newname)
                    os.rename(filename, newname)

                    l.info(f"Processing file: {newname}")
                    print("\nProcessing and splitting model into domains")

                    m = dm.get_model(newname)
                    model_info = process_predicted_model(m, params)

                    chainid_list = model_info.chainid_list
                    print("Segments found: %s" % (" ".join(chainid_list)))

                    mmm = model_info.model.as_map_model_manager()

                    # Write all the domains in one file
                    mmm.write_model(os.path.join(domains_dir,
                                                 f"{PurePosixPath(newname).stem}_domains.pdb"))

                    # Write different domains in different files
                    for chainid in chainid_list:
                        selection_string = "chain %s" % (chainid)
                        ph = model_info.model.get_hierarchy()
                        asc1 = ph.atom_selection_cache()
                        sel = asc1.selection(selection_string)
                        m1 = model_info.model.select(sel)
                        filepath = os.path.join(domains_dir,
                                                f"{PurePosixPath(newname).stem}_{chainid}_RF.pdb")
                        dm.write_model_file(m1, filepath)
                        structures_for_query.append(filepath)

                    conf_domains = extract_residue_list(os.path.join(domains_dir,
                                                                     f"{PurePosixPath(newname).stem}_domains.pdb"),
                                                        domains_dir)

                    conf_domains = extract_residue_list(os.path.join(domains_dir,
                                                                     f"{PurePosixPath(newname).stem}_domains.pdb"),
                                                        domains_dir)
                    l.info(f"Residue list of confident domains: {conf_domains}")
        '''
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

        # Exprt the composite coverage in .csv
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

        MainWindow.hide()
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
        png1 = pio.write_image(fig1, f"{output_dir}/{query_name}/coverage_plot.png", scale=1, width=1400, height=700)

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

        png2 = pio.write_image(fig2, f"{output_dir}/{query_name}/hinges_prediction.png", scale=1, width=1400, height=700)

        return fig2


    def openThirdWindow(self):
        '''
        Creates connection between input window and output window (second window)
        '''
        #SecondWindow.hide()
        self.Thirdwindow = QtWidgets.QMainWindow()
        self.ui = Ui_ThirdWindow()
        self.ui.setupOutput(self.Thirdwindow)
        self.Thirdwindow.show()

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
        self.label.setText(_translate("SecondWindow","<html><head/><body><p align=\"justify\">In this graph,you can see which parts of the reference FASTA sequence are covered by structure. This structures come from either the <a href=\"https://www.rcsb.org/\"><span style=\" text-decoration: underline; color:#0000ff;\">Protein Data Bank</span></a>, <a href=\"https://www.deepmind.com/blog/alphafold-a-solution-to-a-50-year-old-grand-challenge-in-biology\"><span style=\" text-decoration: underline; color:#0000ff;\">AlphaFold</span></a> models or <a href=\"https://www.ipd.uw.edu/2021/07/rosettafold-accurate-protein-structure-prediction-accessible-to-all/\"><span style=\" text-decoration: underline; color:#0000ff;\">RoseTTaFold</span></a> models.<br/></p></body></html>"))
        self.label_2.setText(_translate("SecondWindow", "COVERAGE"))
        self.label_3.setText(_translate("SecondWindow", "HINGES AND FLEXIBILITY"))
        self.label_4.setText(_translate("SecondWindow","<html><head/><body><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\">Flexibility is an important feature of proteins, since they need to move to perform their function and interact with their substrates. In the following section, we provide you with two types of flexibility prediction: the Dynamic Flexibility Index and Hinge Prediction. The overlap of these two measures might be helpful for you, in case you wanted to modify the final topology file with some of these hinges.</span></p><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; font-style:italic; color:#323232;\">Dynamic Flexibility Index</span><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\"><br/>This is per-residue index indicating the contribution of each residue to the overall flexibility of the protein. It uses a method based in an Elastic Network Model (ENM), which is a more lightweight (but less precise, obviously) alternative to Molecular Dynamics. for more info, </span><a href=\"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3673471/\"><span style=\" text-decoration: underline; color:#0000ff;\">here</span></a><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\"> is the original paper.</span></p><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; font-style:italic; color:#323232;\">Hinge Prediction</span><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\"><br/>Hinges are the regions of the protein that allow it to move and change conformations. Using </span><a href=\"https://academic.oup.com/bioinformaticsadvances/article/2/1/vbac007/6525212?login=true\"><span style=\" text-decoration: underline; color:#0000ff;\">this tool</span></a><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\"> we provide you with some suggested hinge regions. Note that this information is only available for experimental structures. This is due to the use of ENM, it is not designed to work with predicted models that might contain important artifacts, and, in this case, that are split into the highly confidently predicted regions.</span></p><p><a name=\"modebar-338c4d\"/><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\"><br/></span></p></body></html>"))
        self.pushButton.setText(_translate("SecondWindow", "Show plot"))
        self.pushButton.clicked.connect(self.show_plot)

        self.pushButton_2.setText(_translate("SecondWindow", "Show plot"))
        self.pushButton_2.clicked.connect(self.show_Secondplot)

        self.pushButton_3.setText(_translate("SecondWindow", "Next"))
        self.pushButton_3.clicked.connect(self.openThirdWindow)


        self.menuPage_1.setTitle(_translate("SecondWindow", "Page 1"))
        self.actionCoverage.setText(_translate("SecondWindow", "Coverage"))
        self.actionHinges_and_Flexibility.setText(_translate("SecondWindow", "Hinges and Flexibility"))
        self.actionComposite_and_Topology_File.setText(_translate("SecondWindow", "Composite and Topology File"))
        self.actionCustom_hinges.setText(_translate("SecondWindow", "Custom hinges"))


class Ui_plot(object):
    '''
    set up window to show image of plot
    '''
    def setup(self, plotWindow):
        plotWindow.resize(1400, 733)
        self.plotwindow = QtWidgets.QWidget(plotWindow)
        self.plotwindow.setObjectName("plotwindow")
        self.label = QtWidgets.QLabel(self.plotwindow)
        self.label.setGeometry(QtCore.QRect(0, -10, 1391, 721))
        self.label.setText("")
        self.label.setPixmap(QtGui.QPixmap(f"{output_dir}/{query_name}/coverage_plot.png"))
        self.label.setObjectName("label")
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
        self.label = QtWidgets.QLabel(self.plotwindow)
        self.label.setGeometry(QtCore.QRect(0, -10, 1391, 721))
        self.label.setText("")
        self.label.setPixmap(QtGui.QPixmap(f"{output_dir}/{query_name}/hinges_prediction.png"))
        self.label.setObjectName("label")
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
        self.label = QtWidgets.QLabel(self.plotwindow)
        self.label.setGeometry(QtCore.QRect(0, -10, 1391, 721))
        self.label.setText("")
        self.label.setPixmap(QtGui.QPixmap(f"{output_dir}/{query_name}/structure_plot.png"))
        self.label.setObjectName("label")
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

        png4 = pio.write_image(fig4, f"{output_dir}/{query_name}/structure_plot.png", scale=1, width=1400, height=700)

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

    def onclick_topology(selected_fragments, n_clicks, output_dir, str_hinges_input):
        structure_list = []
        clicks = n_clicks
        try:
            for child in Path(os.path.join(output_dir, "PDB", "total")).iterdir():
                if child.is_file() and "composite" not in str(child):
                    for name in selected_fragments:
                        if os.path.basename(child)[0:-4] == os.path.basename(name)[0:-13]:
                            structure_list.append(child)
        except:
            pass

        try:
            for child in Path(os.path.join(output_dir, "PDB", "partial")).iterdir():
                if child.is_file() and "composite" not in str(child):
                    for name in selected_fragments:
                        if os.path.basename(child)[0:-4] == os.path.basename(name)[0:-13]:
                            structure_list.append(child)
        except:
            pass

        try:
            for child in Path(os.path.join(output_dir, "PDB", "CHAINS")).iterdir():
                if child.is_file() and "composite" not in str(child):
                    for name in selected_fragments:
                        if os.path.basename(child)[0:-4] == os.path.basename(name)[0:-13]:
                            structure_list.append(child)
        except:
            pass

        try:
            for child in Path(os.path.join(output_dir, "ALPHAFOLD", "DOMAINS")).iterdir():
                if child.is_file() and "confident" not in str(child) and "domains" not in str(child):
                    for name in selected_fragments:
                        if str(os.path.basename(child)[0:-4]) == str(os.path.basename(name)[0:-13]):
                            structure_list.append(child)
        except:
            pass

        try:
            for child in Path(os.path.join(output_dir, "ROSETTAFOLD", "DOMAINS")).iterdir():
                if child.is_file() and "confident" not in str(child) and "domains" not in str(child):
                    for name in selected_fragments:
                        if str(os.path.basename(child)[0:-4]) == str(os.path.basename(name)[0:-13]):
                            structure_list.append(child)
        except:
            pass

        fasta = "input_fasta/" + str(os.path.basename(output_dir)) + ".fasta"

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
        str_out = str(output_dir)
        out_name = str_out.split("/")[-1]
        # Write the topology file
        write_custom_topology(os.path.join(output_dir, "IMP", f"{out_name}_custom.topology"), final_rigid_bodies)

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
        self.textEdit = QtWidgets.QTextEdit(self.SecondOutputWindow)
        self.textEdit.setGeometry(QtCore.QRect(250, 500, 201, 31))
        self.textEdit.setObjectName("textEdit")
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
                                      "<html><head/><body><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\">From all the structures retrieved by the program and provided by the user, the program generates this composite, trying to cover as much of the reference sequence as possible, avoiding overlaps.</span></p><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\">This comopsite can be used to automatically build an </span><a href=\"https://integrativemodeling.org/2.5.0/doc/ref/classIMP_1_1pmi_1_1topology_1_1TopologyReader.html\"><span style=\" text-decoration: underline; color:#0000ff;\">IMP topology file</span></a></p><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\">If you want to generate another topology file, with your custom choice of fragments, simply select the fragments woy want to include and click the button &quot;Create Topology File&quot;. it will save yhe output inj the folder IMP of your selected output folder.</span></p></body></html>"))
        self.label_2.setText(_translate("ThirdWindow", "COMPOSITE AND TOPOLOGY FILE"))
        self.label_3.setText(_translate("ThirdWindow", "CUSTOM HINGES"))
        self.label_4.setText(_translate("ThirdWindow",
                                        "<html><head/><body><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\">In this section you can introduce hinge regions. Hinge regions are those regions of the protein that bend, allowing the movement of the more rigid domains, which is essential for the interaciton of the proteins with other biomolecules.</span></p><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; font-weight:600; color:#323232;\">How are hinges encoded?</span><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\"><br/>Let\'s imagine you have a proterin of 200 amino acids. The DFI and PACKMAN hinge prediction are indicating a putative hinge region between positions 50 and 100 and another one between 120 and 130.</span></p><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\">In the box, you will need to introduce the hinges with the following format: 50:100,120:130</span></p><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\">The program will split the structures, in the topology file, according to the hinges introduced.</span></p><p align=\"justify\"><span style=\" font-family:\'Open Sans,HelveticaNeue,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:15px; color:#323232;\">Topology file created with:</span></p></body></html>"))
        self.pushButton.setText(_translate("ThirdWindow", "Show plot"))
        self.pushButton.clicked.connect(self.show_Compositeplot)

        self.pushButton_2.setText(_translate("ThirdWindow", "Generate IMP Topology File"))
        self.menuPage_1.setTitle(_translate("ThirdWindow", "Page 1"))
        self.actionCoverage.setText(_translate("ThirdWindow", "Coverage"))
        self.actionHinges_and_Flexibility.setText(_translate("ThirdWindow", "Hinges and Flexibility"))
        self.actionComposite_and_Topology_File.setText(_translate("ThirdWindow", "Composite and Topology File"))
        self.actionCustom_hinges.setText(_translate("ThirdWindow", "Custom hinges"))



if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    # GUI pops up
    MainWindow.show()
    loop = QEventLoop()
    loop.exec()  # waits



    sys.exit(app.exec_())


