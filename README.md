# Graphical User Interface 
Here is the source code for an application about **fpegenaute Master Thesis for the MSc in Bioinformatics for Health Sciences** using Python and Pyqt5 library.
The idea of his project (["Automated structural information retrieval for Integrative Modeling"](https://github.com/fpegenaute/TFM/blob/main/README.md#automated-structural-information-retrieval-for-integrative-modeling)) was to create a tool that combines information from already existing and experimental data in order to predict the dynamics of proteins.
Basically, we gather information from PDB (Protein Data Bank) and AlphaFold or RoseTTaFold and use it as an input for IMP ([Integrative Modeling Platform](https://integrativemodeling.org/)), which gives us an hybrid approach integrating data from diverse biochemical and biophysical experiments.

## Installation and Requirements
The Graphical User Interface, as the original program, uses external programs to work:
 - **BLAST** [Install BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
 - **AlphaFold**, non-docker set up (optional) [Install AlphaFold](https://github.com/kalininalab/alphafold_non_docker)
 - **RoseTTaFold** (optional) [Install RoseTTaFold](https://github.com/RosettaCommons/RoseTTAFold)

Also, to install the Python packages you will need Conda and [pip](https://pip.pypa.io/en/stable/installation/).

## Set up the program
 1.- Open a Terminal
 
 2.- Download this repository
 ```ruby
 git clone https://github.com/Mroyyy/GUI.git
 ```
 
 3.- **Create a conda environment** to avoid dependency issues
```ruby
conda create --name your_env python=3.8
```

4.- Install the **Python packages** needed by the program in the working directory
```ruby
cd path/to/working/directory
pip install -r requirements.txt
```
- **What is recommended**, if you have conda installed, you can create the environment with all the dependencies as:
```ruby
conda env create -f environment.yaml --name your_env
```

5.- **Configure the program**. Open the file config.py inside the bin/ folder, and in the line 6 change:
```ruby
"blastdb" : "/path/BLAST/database"
```

 
 ## Practical example: SEC3
 Sec3 is an essential protein which participates in exocytosis, whose full structure is still unkown. Its sequence in FASTA format is on "input_fasta" directory. **All sequences you want to study have to be in the same directory as the program to work**

There are instructions in a Help Window, but to start you will only need two arguments:

 - **Input sequence** in FASTA format
 - **Output directory** to store the retrieved PDBs

1.- Activate your environment!
```ruby
conda activate your_env
```
2.- Run GUI
```ruby
python3 gui.py
```
