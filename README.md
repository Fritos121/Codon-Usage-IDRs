# Codon-Usage-IDRs

Codon-Usage-IDRs is a repository storing the pipeline used in my basic analysis of the conservation of codon usage bias in Intrinsically Disordered Regions.

## Packages Used

A number of python packages were used in this pipeline:
  1. Biopython (v1.72)
  1. Bioservices (v1.7.4)
  1. Numpy (v1.15.4)
  1. Pandas (v0.24.2)
  1. Matplotlib (v3.0.2)


## Other Software Used

VSL2 - Used to predict protein disorder (coming to this repo soon)

[Jupyter Notebook](https://jupyter.org/install) - used to visualize and analyze data


## Databases Accessed

[PANTHER](http://pantherdb.org/) was used to retrieve orthologs and the family multiple sequence alignment.

[UniProt](https://www.uniprot.org/) was used to retrieve each ortholog’s taxonomy IDs and EMBL coding sequence accession numbers.

[European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/home) was used to retrieve the orthologs’ coding sequences.

[National Center for Biotechnology Information (NCBI)](https://www.ncbi.nlm.nih.gov/) was used to retrieve all protein coding sequences for an organism.


## Usage

All programs in the pipeline are command-line compatible. Each program's description, required input, and outputs are provided below.

* <ins>get_ortholog_info.py</ins>:

  * <ins>Description</ins>: Retrieve orthologs and their respecitve family multiple sequence alignments for each protein.
  
  * <ins>Input Files</ins>: CSV file with one or more sets of UniProt IDs, PANTHER Family IDs, and Taxonomy IDs, in that order.
  
  * <ins>Output</ins>: A text file of multiple sequence alignments in Fasta format.

* <ins>get_species_info_by_uid.py</ins>:

  * <ins>Description</ins>: Retreive taxonomy ID's and coding sequences. Filter out any protein with incomplete information.
  
  * <ins>Input Files</ins>: A text file of multiple sequence alignments in Fasta format where UniProtKB IDs are the headers.
  
  * <ins>Output</ins>: 1) CSV file containg partial links to the NCBI FTP server, and the corresponding taxonomy ID. 2) Fasta file of all DNA coding sequences from input list. 3) Error Log text file recording which proteins were removed.

* <ins>get_assembly_seqs_from_ftp.py</ins>:

  * <ins>Description</ins>: Download file from NCBI FTP that contains all the protein coding sequences in an organism.
  
  * <ins>Input Files</ins>: CSV file containg partial links to the NCBI FTP server, and the corresponding taxonomy ID.
  
  * <ins>Output</ins>: Fasta file containing all protein coding sequences for an organism.
  
  
* <ins>get_distributions.sh</ins>:

  * <ins>Description</ins>: Create CSV file(s) containing the codon distributions for: 1) Each organism in directory 2) Each gene in each Fasta file in subdirectories (-m option). Uses codon_dist_from_fasta.py for codon distribution calculation.
  
  * <ins>Input Files</ins>: Fasta file containing all protein coding sequences for an organism.
  
  * <ins>Output</ins>: CSV file containing codon distribution, formatted as "codon,amino acid"





Below is a visualization of the pipeline.
![Pipeline Flowchart](pipeline.png)
