# CONTACT 
## COmparisoN beTween mAgs CharacTeristics and genomes of isolates

This code was developed as part of a project carried out during the course of Microbial Metagenomics (Molecular Biology master degree) at the University of Padova. The project was supervised by Prof. Stefano Campanaro, Dr. Maria Silvia Morlino, Dr. Edoardo Bizzotto, and Dr. Gabriele Ghiotto. 
The code was developed by Anna Bortolato, Alessia Prestanti, Sara Hosseini.

The script was developed in order to verify the coherency between the characteristics of the MAGs and those reported for NCBI genomes belonging to the same taxonomic group, with particular reference to the discrepancy in terms of genome size and GC content taking into account the MAG’s completeness and level of contamination. 

## Requirements

The following python libraries are required in order to run the script : `pandas`, `numpy`, `sys`, `seaborn`, `matplotlib.pyplot`.

## How to generate input files

The input files the user has to provide in order to run the script are:
 
1. The output file of `CheckM`, that is a tabular file in a text format, the header is the first row and the first column (MAG Id) is set as index 
2. A tabular file in `txt` format with MAGs classification from gtdb and NCBI database; this file is generated with the script `gtdb_to_ncbi_majority_vote.py`, the header is the first row and the first column is set as index
3. The NCBI taxonomy database for prokaryotes in `txt` format. To generate this table, the summary from NCBI database must be downloaded (`assembly_summary.txt`). Using `ete3` toolkit *taxids* are converted in *taxranks*. Then, the complete taxonomy for each organism is generated with `ncbi.get_taxid_translator` function; this table should have no header.

## Information on the input requireds and how to use the software

In order to use our script the command line must be organized as in the following example:

`python3 contact.py checkM_output.txt gtdb_NCBI_conversion.txt NCBI_prokaryotes_taxonomy.txt taxon min_number_genomes`

The first and the second terms are fixed and they represent the Python version and the name of the script, then we have 3 input files in a tabular `txt` format:

- `checkM_output.txt` is the output of `CheckM` including a series of features for each MAG such as the _genome size_, _number of scaffolds_, _GC content_, _completeness_, _contamination_ etc.;
- `gtdb_NCBI_conversion.txt` contains the taxonomic assignment for each MAG from Genome Taxonomy Database (gtdb) and its conversion to NCBI taxonomy;
- `NCBI_prokaryotes_taxonomy.txt` contains genome characteristics for each prokaryotic isolate deposited in NCBI, including _genome size_ (column 6) and _GC content_ (column 7), and its _taxonomy_ from root to species level (column 23);
- `taxon` is the taxonomic level defined by the user in order to perform the analysis: it could be chosen among _genus_, _family_ or _species_ in this implementation of the software;
- `min_number_genomes` is the minimum number of NCBI reference genomes associated with a certain taxonomic level among _genus_, _family_ or _species_ that the user wants to consider to perform the analysis, if not specified it is automatically set to 5.

## Output description

The software generates in the working directory a tabular output in `csv` format in which each row represents a specific MAG, identified by an ID and associated with its taxonomy (genus, family or species), which is reported in the second column. Other features reported in columns are MAG’s completeness, genome size and contamination and both the average and the standard deviation of genome size and GC content of NCBI reference genomes belonging to the same taxonomic group. Users can find columns showing the distance (difference in size) normalized by the standard deviation between the average genome size of the reference genomes and MAG’s size and the distance between the average GC content and the MAG’s GC content. Two additional columns report the expected genome size for each MAG and the expected missing or exceeding genomic portion.

Two graphical outputs are saved automatically in the working directory in pdf format; one of them is reporting for each MAG the distance from the average genome size, measured in number of standard deviations, as a function of the MAG’s completeness. Each dot identifies a specific genome and the distance from 0 represents the number of standard deviations from the average genome size of prokaryotes from NCBI according to a specific taxonomic group. Colors identify different levels of contamination. 
The second graph shows the completeness in the y axis as a function of the missing or exceeding genomic portion, which can explain the difference in size between MAGs and the average genome size. 

