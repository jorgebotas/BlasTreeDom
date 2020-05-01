# **_BlastTreeDom_**


#### Given query FASTA file(s) and genBank file(s):
* Parse genBank file(s) CDS protein sequences and info
* Perform blastp analysis
* Compute Neighbour-Joining phylogenetic tree(s) using MUSCLE
* Extract ProSite domains of each of the hit sequences from blastp and present  
* Graph blastp and ProSite-domain output   

#### Required packages
* biopython  
* blast  
* muscle  


## Required package installation

[Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) or [MiniConda](https://docs.conda.io/en/latest/miniconda.html) installation is highly recommended  

The following installation will be performed using conda from the command line  

#### [biopython](https://anaconda.org/anaconda/biopython)  
`conda install -c anaconda biopython`  

#### [blast](https://www.ncbi.nlm.nih.gov/books/NBK279671/)  
* Via [Anaconda](https://anaconda.org/bioconda/blast)  
  `conda install -c bioconda blast`   
* Via [NCBI website](https://www.ncbi.nlm.nih.gov/books/NBK279671/)  

#### [muscle](https://anaconda.org/bioconda/muscle)  
`conda install -c bioconda muscle`  
  
  
## USAGE  

### User-friendly  
After cloning and _cd_ into **_BlasTreeDom_** repository:  
* Run `python main.py -ui` from the command line  

Following steps will be prompted on the Terminal  

### One-shot command  
Run `python main.py [arguments]`  
  
As arguments the following will be expected:  
* Query FASTA file(s) or directory containing several files:  
  `-query query(.fasta)`  
  
* Subject sequence against which to perform the analysis:  
  - A genBank file or directory containing several files from which to parse the CDS protein sequences and info  
 `-genBank genBank`  
  - Alternatively, provide multifasta file(s)  
  `-multifasta subject(.fasta)`  
    - If a database from sequences have already been computed, it can be provided additionally as:  
      `-database database`
   

## Output  

Results will be stored in a directory named by date and time of command execution.  
Within this directory the following can be found: a _log_ file, a _database_ folder, merged _query_ and _subject_ FASTA files, _blast_ and _merged_ output files, and a separate folder for each input query   

### blast  

* A tsv file containing all results  
* Graphical representation of output:  

![](images/blast.png)  


### Alignment and N-J tree  

* Unaligned and aligned sequences will be stored in each query directory  
* Neighbor-Joining tree computed using MUSCLE  

### ProSite domains  

* A tsv file containing information on: sequence accession number (id) and name, accession, description, pattern, start and end location of found ProSite domains  
* A txt file for each query and hit sequence containing the domains information plus additional text parsed from prosite.doc  
* Graphical representation of found domains for each protein sequence:  

![](images/domains_amplified.png)  

### Merged results  

* A tsv file containing _blast_output_, genBank parsed fields and extracted ProSite domain names to create an integrated ouput file
