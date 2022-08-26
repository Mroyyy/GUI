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

![helpred](https://user-images.githubusercontent.com/93529369/186919664-9dc2c591-ce28-41bb-9f7e-2d7f819590be.png)


![help](https://user-images.githubusercontent.com/93529369/186918736-843cd346-d96f-4d55-a026-d0e7ed606eda.png)

1.- **Activate your environment**
```ruby
conda activate your_env
```
2.- **Run GUI**
```ruby
python3 gui.py
```

![Screenshot from 2022-08-26 16-01-57](https://user-images.githubusercontent.com/93529369/186921187-c3aa6067-9ca0-4528-a8c1-28381872e9f4.png)

3.- **Output**
 - **Coverage**: Which parts of the FASTA are covered
 - **Hinges and flexibilility**: Hinges are regions of the protein that allow it to move and change conformations. In this plot you can observe its prediction
 - **Composite and Topology File**:  From all the structures retrieved by the program, a composite is generated, trying to cover as much of the reference sequence as possible, avoiding overlaps.
 - **Custom hinges**: In this section, hinge regions can be introduced, allowing the movement of the more rigid domains.
   An example is shown below:
   
![hinges](https://user-images.githubusercontent.com/93529369/186927248-be52fb78-ad76-4ecc-9eb0-9fa04acbe816.png)

**This custom topology file will be the input file for IMP:**

![top](https://user-images.githubusercontent.com/93529369/186927563-05ea9067-d66c-4e62-bc57-bdabb61e031a.png)

**IMP Output**

![imp](https://user-images.githubusercontent.com/93529369/186927966-b44bdaa4-2981-4cb3-8804-f36e063e3c44.png)



At this pont, the output directory should contain the files and directories described bellow:
```
SEC3/
 BLAST/
  SEC3_blast.out
 FASTA/
  5LG4.fa
  5YFP.fa
 HINGES/
  5lg4_B.hng
  5lg4_B_packman_output.txt
  5yfp_A.hng
  5yfp_A_packman_output.txt
 IMP/
  SEC3.topology
  *SEC3_custom.topology*
 LOG/
  SEC3.log
 PDB/
  CHAINS/
  partial/
   5lg4_B.pdb
   5yfp_A.pdb
  5lg4.cif
  5yfp.cif
 PLOTS/
  coverage_plot.html
  coverage_plot.png
  hinges_prediction.html
  hinges_prediction.png
  structure_plot.html
  structure_plot.png
 REPORT/
  COVERAGE/
   5lg4_B_coverage.csv
   5yfp_A_coverage.csv
   SEC3_composite_coverage.csv
  DFI/
   5lg4_B_DFI_coverage.csv
   5yfp_A_DFI_coverage.csv
```

## FAQS
**Why in the main window there are more labels than used?** The idea is to implement AlphaFold and RoseTTaFold models (or to run the server directly) and obtain the most complete structure as possible. If it was the case in the Composite and Topology File we could also observe and select fragments from AlphaFold model for example.
**I can only run it on Linux?** Currently you have to run the GUI from the terminal, but the future idea is to obtain an .exe one file with all packages and needed files.
