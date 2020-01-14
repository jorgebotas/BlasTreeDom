# **_BlastTreeDom_**


#### Given query FASTA file(s) and genBank file(s):
* Parse genBank file(s) CDS protein sequences
* Perform blastp analysis
* Compute Neighbour-Joining phylogenetic tree(s) using MUSCLE
* Extract ProSite domains of each of the hit sequences from blastp

#### Required packages
* Bio
* blast+
* muscle
* numpy
* pandas


## Required package installation

[Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) or [MiniConda](https://docs.conda.io/en/latest/miniconda.html) installation is highly recommended

The following installation will be performed using conda from the command line

#### [biopython](https://anaconda.org/anaconda/biopython)
`conda install -c anaconda biopython`
#### blast
* Via [Anaconda](https://anaconda.org/bioconda/blast)
* Via [NCBI website](https://www.ncbi.nlm.nih.gov/books/NBK279671/)
#### [muscle](https://anaconda.org/bioconda/muscle)
`conda install -c bioconda muscle`
` conda install -c bioconda/label/cf201901 muscle `
#### [numpy](https://anaconda.org/anaconda/numpy)
`conda install -c anaconda numpy`
#### [pandas](https://anaconda.org/anaconda/pandas)
`conda install -c anaconda pandas`


## USAGE

### User-friendly
After cloning and cd into **_BlasTreeDom_** repository:
* Run `python main.py -ui` from the command line

Following steps will be prompted on the Terminal

### One-shot command
Run `python main.py [arguments]`
As arguments the following will be expected:
* A query FASTA file or a directory containing several files:\
  `-query query.fasta`
  
* Subject sequence againt which to perform the analysis:
  - A genBank file or directory containing several files from which to parse the CDS protein sequences\
 `-genBank genBank`
  - Alternatively a multifasta file\
  `-multifasta subject.fasta`
  
* If no blast analysis is required, unaligned sequences in a FASTA file can be provided\
`-unaligned unaligned.fasta`
