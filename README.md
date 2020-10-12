UPDATE: I now consider the MOBtyping tool to be deprecated (Oct 2020). I would intead recommend using the MacSyFinder program with the CONJscan module to identify conjugative systems and assign MOB types. For more information on MacSyFinder/CONJscan, see [Cury et al. 2020](https://pubmed.ncbi.nlm.nih.gov/31584169/).<br>


# Introduction

MOB typing is a method of plasmid classification. It is based on the detection of 'relaxase' proteins, which are involved in plasmid mobility. To MOB type plasmids from bacteria of the Gammaproteobacteria group, 6 relaxase proteins are used as queries to BLAST search a plasmid dataset, in order to detect relaxase homologs. Each relaxase query represents a protein family and corresponds to one of 6 MOB types (MOBC, MOBF, MOBH, MOBP, MOBQ, MOBV). Standard non-iterative BLAST has been commonly used for MOB typing. However, this approach is likely to lack power because relaxase proteins - even within the same relaxase protein family - are diverse. This repository provides scripts for detecting relaxases / assigning MOB types, using more powerful (iterative) PSI-BLAST searching. PSI-BLAST harnesses position-specific information about conservation of relaxase protein residues amongst plasmids in a database. A curated plasmid database is provided in an accompanying [Figshare repository](https://figshare.com/s/18de8bdcbba47dbaba41). The bioinformatic steps can be summarised as follows: one or more untyped plasmids (provided by the user) are appended to the curated plasmid database; 6 MOB-relaxase proteins are then queried against the updated plasmid database using PSI-BLAST.


# Installation

To install scripts, navigate to an empty directory and run the following code:

```

git clone https://github.com/AlexOrlek/MOBtyping.git

```

# Software dependencies and input data requirements

The scripts have been tested using Ubuntu Version 16.04 and Python Version 2.7. 
[Command-line NCBI-BLAST](https://www.ncbi.nlm.nih.gov/books/NBK52640/#_chapter1_Installation_) will need to be installed

To use the scripts, a database of plasmids (as protein sequences) must be added to the '/plasmiddatabase' directory - you can download the curated plasmid database ('translatedproteinseq.fa') provided in the Figshare repository (linked above). In addition, one or more untyped plasmids must be added to the '/untyped_plasmids' directory in either FASTA format or Genbank format. The file name extensions must be either '.fa' or '.gb'. If retrieving a Genbank file from NCBI, download the 'full' file which includes the nucleotide sequence. A Genbank file retrieved from NCBI is likely to include annotations indicating CDS (coding) features - this information will be used when selecting the best PSI-BLAST hit at a locus (a locus will be defined based on the position of the CDS feature). If untyped plasmids are instead provided as a FASTA file, CDS feature annotations can be provided separately in a .tsv file with the first 4 columns containing the following information: plasmid accession, CDS nucleotide start position, CDS nucleotide end position, CDS description (e.g. protein product), plasmid accession length (bp). If no CDS feature information is provided, the code will run using a more simplistic method for selecting the best PSI-BLAST hit at a locus (a locus will be defined based on alignments overlapping at the same position).


# Usage

A parent python script (parentscript.py) calls all other necessary scripts using a custom function that invokes the Python [subprocess module](https://docs.python.org/2/library/subprocess.html). The parent script takes arguments from the command line:

1. the name (excluding file extension) of the FASTA or Genbank file containing the untyped plasmid(s).  
2. the number of PSI-BLAST iterations. We recommend 14 iterations.
3. the number of computer threads to allocate

Therefore, if 'untypedplasmids.gb' has been added to the '/untyped_plasmids' directory, the parent script is called as follows, running 14 PSI_BLAST iterations and allocating 10 threads:

```

python parentscript.py untypedplasmids 14 10

```

A brief explanation of each script invoked by the parent script is given below:

**plasmidtranslate.py**  untyped plasmid(s) are translated in all 6 frames, and the resulting protein sequences are added to the original plasmid database ('translatedproteinseq.fa') to create the file 'name'_plasmidproteins.fa.

**cdsfeatures.py**  if available, information on the CDS features encoded on the untyped plasmids is extracted and appended to existing information on plasmids in the original database, to create the file pickle.CDSfeatures_'name'

**makeblastdb.py**  a protein BLAST database is created from 'name'_plasmidproteins.fa.

**psiblast.py**  MOB queries are searched against the BLAST database, running through iterations up until iteration 14 (or a custom maximum iteration).

**additerationnumber.py**  hits produced from each MOB query are combined per iteration; in the file containing all hits for the maximum iteration, the iteration at which each hit was retrieved is recorded.

**blastfilter.py**  where multiple hits align at the same locus, a best hit is selected, favouring hits retrieved at lower iteration, and having a higher coverage/percent identity score.

**mobtyping.py**  the best hit(s) associated with each plasmid accession are used to create a file recording the MOB type for each plasmid. This script outputs two 'FINAL_OUTPUT' files (one includes MOB typing results for previously untyped plasmids only, the other also includes database plasmid typing).



# Example

After cloning the repository, you will find a Genbank file of a plasmid called pNUC (accession KU852461.1) which is not present in the plasmid database, having been added to NCBI since the original database was compiled. For details see Olivia et al. 2017 (article provided in '/untyped_plasmids' directory). Make sure you have also downloaded the plasmid database and added it to the '/plasmiddatabase' directory, as described above. To MOB type the untyped pNUC plasmid, run the following code from the command line:

```

python parentscript.py pNUC 5 10

```

Navigate to the '/output_final' directory and select the file 'FINAL_OUTPUT_mobtyped_plasmids_1_5.tsv'. You will find that pNUC has been typed as MOBQ. To explain the '1_5' suffix: '1' indicates that the file has been generated using the first column of evalues specified in the 'mobtyping_evalues' file (by default this is the only set of evalues that can be selected); '5' is the maximum number of iterations run (in this example, 5 iterations have been used for speed).



# Outputs

The '/output_intermediate' directory stores intermediate files. By default all intermediate files are removed when no longer needed (see bottom of blastfilter.py script).

The '/output_final' directory contains 4 .tsv files. The two 'FINAL_OUTPUT...' files generated by mobtyping.py (see above) and two other files showing the PSI-BLAST hits immediately prior to best hit selection ('all_hits') and after best hit selection ('best_hits').


# Caveats

The scripts are designed for MOB typing plasmids of the Gammaproteobacteria group, which are thought to belong to 6 MOB types. Likewise, the curated plasmid database (linked above) includes only Gammaproteobacterial plasmids (of the Enterobacteriaceae family). Background information and details on how evalues were optimised can be found in an associated [research paper](http://www.sciencedirect.com/science/article/pii/S0147619X16301032). A further 2 relaxase queries are required to MOB type plasmids from outside the Gammaproteobacteria group [(Guglielmini et al., 2011)](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002222). If the aim is to type plasmids from a broad taxonomic range of bacteria, scripts would need to be modified accordingly, in order to incorporate the additional 2 MOB qeuries, and you may want to re-optimise evalue threhsolds. The reference plasmid database would need to be expanded to include non-Gammaproteobacterial plamsids, such that all 8 relaxase protein families would be represented in the database.

BLAST search evalues are designed to reduce false-positive alignments. Increasing database size inflates evalues, making a given alignment less likely to pass the evalue theshold, and accounting for the fact that more false-positive alignments will be detected when a larger database is searched. Therefore, you should be aware that database size will influence typing results and if a large number of untyped plamids are appended to the database this will reduce the power of the PSI-BLAST search.
