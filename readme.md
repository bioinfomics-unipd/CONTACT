# CONTACT 
## COmparisoN beTween mAgs CharacTeristics and genomes of isolates

This code was developed as part of a project carried out during the course of Microbial Metagenomics (Molecular Biology master degree) at the University of Padova. The project was supervised by Prof. Stefano Campanaro, Dr. Maria Silvia Morlino, Dr. Edoardo Bizzotto, and Dr. Gabriele Ghiotto. 
The code was developed by Anna Bortolato, Alessia Prestanti, Sara Hosseini.

The script was developed in order to verify the coherency between the characteristics of the MAGs and those reported for NCBI genomes belonging to the same taxonomic group, with particular reference to the discrepancy in terms of genome size and GC content taking into account the MAG’s completeness and level of contamination. 

## Requirements

The following python libraries are required in order to run the script : `pandas`, `numpy`, `sys`, `seaborn`, `matplotlib.pyplot`.

## How to generate input files

The input files the user has to provide in order to run the script are:
 
1. The output file of `CheckM`, that is a tabular file in a text format 
2. A tabular file in `txt` format with MAGs classification from gtdb and NCBI database (output of the script `gtdb_to_ncbi_majority_vote.py`)
3. The NCBI taxonomy database for prokaryotes in `txt` format. To generate this table, the summary from NCBI database must be downloaded (`assembly_summary.txt`). Using `ete3` toolkit *taxids* are converted in *taxranks*. Then, the complete taxonomy for each organism is generated with `ncbi.get_taxid_translator` function.

## Information on the input requireds and how to use the software

In order to use our script the command line must be organized as in the following example:

`python3 contact.py checkM_output.txt gtdb_NCBI_conversion.txt NCBI_prokaryotes_taxonomy.txt taxon min_number_genomes`

The first and the second terms are fixed and they represent the Python version and the name of the script, then we have 3 input files in a tabular `txt` format:

- `checkM_output.txt` is the output of `CheckM` including a series of features for each MAG such as the _genome size_, _number of scaffolds_, _GC content_, _completeness_, _contamination_ etc.;
- `gtdb_NCBI_conversion.txt` is a tabular file with the Genome Taxonomy Database (gtdb) assignment for each MAG and its conversion to NCBI taxonomy;
- `NCBI_prokaryotes_taxonomy.txt` is a tabular file with information for each isolate in the NCBI database including _genome size_ (column 6), _GC content_ (column 7) and _taxonomy_ from root to species level (column 23);
- `taxon` is the taxonomic level defined by the user in order to perform the analysis: it could be chosen among _genus_, _family_ or _species_ in this implementation of the software;
- `min_number_genomes` is the minimum number of NCBI reference genomes associated with a certain taxonomic level among _genus_, _family_ or _species_ that the user wants to consider to perform the analysis, if not specified it is automatically set to 5.

## Output description

The software generates in the working directory a tabular output (`csv` format) in which each row represents a specific MAG, identified by an ID and associated with its taxonomy (genus, family or species), reported in the second column. Other features reported in columns are MAG’s completeness, genome size and contamination and both the average and the standard deviation of genome size and GC content of NCBI reference genomes belonging to the same taxonomic group. Two additional columns show the distance normalized by the standard deviation between the average genome size of the reference genomes and MAG’s size and the distance between the average GC content and the MAG’s GC content.  

The graphical output is saved automatically in the working directory in `pdf` format and it is reporting for each MAG the distance from the average genome size, measured in number of standard deviations, as a function of the MAG’s completeness. Each dot identifies a specific genome and the distance from 0 represents the number of standard deviations from the average genome size of prokaryotes from NCBI according to a specific taxonomic group. Colors identify different levels of contamination.
