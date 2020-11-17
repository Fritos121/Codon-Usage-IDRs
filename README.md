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

get_ortholog_info.py:

  Description: 
  
  Required Input:
  
  Output:

get_species_info_by_uid.py


Below is a visualization of the pipeline.
![Pipeline Flowchart](pipeline.png)
