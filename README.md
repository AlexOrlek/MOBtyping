#Introduction

MOB typing is a method to classify plasmids. It is based on the detection of 'relaxase' proteins, which are invovled in plasmid mobility. To MOB type plasmids from bacteria of the Gammaproteobacteria group, plasmids are BLAST searched against 6 relaxase proteins; each represents a protein family and corresponds to one of 6 MOB types (MOBC, MOBF, MOBH, MOBP, MOBQ, MOBV). Standard non-iterative BLAST has been commonly used for MOB typing. However, this approach is likely to lack power because relaxase proteins - even within the same relaxase protein family - are diverse. This repository provides a curated plasmid database and scripts for detecting relaxases / assigning MOB types, using more powerful (iterative) PSI-BLAST searching. PSI-BLAST harnesses position-specific information about conservation of relaxase protein residues amongst plasmids in the database. The bioinformatic steps can be summarised as follows: one or more untyped plasmids (provided by the user) are appended to the existing database, and then 6 MOB-relaxase proteins are queried against the updated plasmid database using PSI-BLAST.


#Installation

To install scripts and download the plasmid database, navigate to an empty directory and run the following code:

```

git clone https://github.com/AlexOrlek/MOBtyping.git                                                                                                                                                      

```

#Requirements

The scripts have been tested using Ubuntu Version 16.04 and Python Version 2.7. Up to 6 threads will be used to run code. To use the scripts, one or more untyped plasmids must be added to the '/untyped_plasmids' directory in either FASTA format or Genbank format. The file name extensions must be either '.fa' or '.gb'. If retrieving a Genbank file from NCBI, download the 'full' file which includes the nucleotide sequence. A Genbank file retrieved from NCBI is likely to include annotations indicating CDS (coding) features - this information will be used when selecting the best PSI-BLAST hit at a locus (a locus will be defined based on the position of the CDS feature). If a FASTA file is provided instead, CDS feature annotations can be provided separately in a .tsv file with first 4 columns containing the following information: plasmid accession, CDS nucleotide start position, CDS nucleotide end position, CDS description (e.g. protein product), plasmid accession length (bp). If no CDS feature information is provided, the code will run using a more simplistic method for selecting the best PSI-BLAST hit at a locus (based on overlapping position of hits).


#Usage

A parent python script (parentscript.py) calls all other necessary scripts using a custom function that invokes the subprocess module. The parent script takes arguments from the command line:

1. the name (excluding file extension) of the FASTA or Genbank file containing the untyped plasmid(s).  
2. the number of PSI-BLAST iterations. We recommend 14 iterations, and at least 5 iterations.

Therefore, if 'untypedplasmids.gb' has been added to the '/untyped_plasmids' directory, the parent script is called from the command line as follows:

```

python parentscript.py untypedplasmids 14

```

A brief explanation of each script invoked by the parent script is given below:

**plasmidtranslate.py**  untyped plasmid(s) are translated in all 6 frames, and the resulting protein sequences are added to the plasmid database (in FASTA format) to create the file <name>_plasmidproteins.fa.

**cdsfeatures.py**  if available, CDS feature information on the untyped plasmids is extracted and appended to information on plasmids in the origian database to create the file pickle.CDSfeatures_<name>

**makeblastdb.py**  a protein BLAST database is created from <name>_plasmidproteins.fa.

**psiblast.py**  MOB queries are searched against the BLAST database, running through iterations up until iteration 14 (or custom maximum iteration).

**additerationnumber.py**  hits produced from each MOB query are combined per iteration; in the file containing all hits for the maximum iteration, the iteration at which each hit was retrieved is recoreded.

**blastfilter.py**  where multiple hits align at the same locus, a best hit is selected, favouring hits retrieved at lower iteration, and having a higher coverage/percent identity score.

**mobtyping.py**  the best hit(s) associated with each plasmid accession are used to create a file recording the MOB type for each plasmid (previously untyped plasmids are included at the bottom).



#Example

After cloning the repository, you will find a Genbank file of a plasmid called pNUC (accession KU852461.1) which is not present in the plasmid database, having been added to NCBI since the original database was compiled. To MOB type this plasmid, run the following code from the command line:

```

python parentscript.py pNUC 14

```

Navigate to the '/output_final' directory and select the file 'FINAL_OUTPUT_mobtyped_plasmids_1_14.tsv'. At the bottom of the file you will find that pNUC has been typed as MOBQ. To explain the '1_14' suffix: '1' indicates that the file has been generated using the first column of evalues specified in the 'mobtyping_evalues' file (by default this is the only set of evalues used); '14' is the maximum number of iterations run.



#Outputs

The '/output_intermediate' directory stores intermediate files. By default all intermediate files are deleted (see bottom of blastfilter.py script).

The '/output_final' directory contains 3 .tsv files. The 'FINAL_OUTPUT...' file shows [plasmid accession names of plasmids - including the previously untyped plasmid(s) - and associated MOB types. The other files show the PSI-BLAST hits immediately prior to best hit selection ('all_hits') and after best hit selection ('best_hits').


#Caveats

The scripts are designed for typing plasmids of the Gammaproteobacteria group, which are thought to belong to 6 MOB types. However, other researchers have used a further 2 relaxase queries to type plasmids from outside Gammaproteobacteria - scripts would need to be modified if this is the aim.

Increasing database size inflates evalues, potentially influencing the overall MOB typing results. The scripts are intended to be used for typing a small number of untyped plasmids, such that the size of the plasmid database is not changed considerably after appending the untyped plasmids.
