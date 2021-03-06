OVERVIEW
========

This folder contains Python code used to determine significant transcription factor (TF) binding in regions around the transcription start sites (TSS) of genes in particular clusters. The program works by finding the proportion of all genes with TF binding in a region surrounding a TSS and determining if the number of genes with TF binding in a particular cluster is significantly larger than would be expected by chance. This script was used to produce results in the paper “Inference of RNA Polymerase II Transcription Dynamics from Chromatin Immunoprecipitation Time Course Data,” 
by Ciira wa Maina et al. published in [PLOS Computational Biology](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003598).


REQUIREMENTS
============
The programs require Python 2.7 or later and the following python libraries

1. numpy >= 1.6.1
2. sys
3. os
4. string

Also  intersectBed from BEDtools must be installed. Follow these [installation instructions](http://bedtools.readthedocs.org/en/latest/content/installation.html).

INSTALLATION
============
Download the following files in this folder and place them in a folder of your choice.

1. TFbind.py

In the same folder include the following subfolders

1. Clusters/, This folder contained text files named cluster\_X.txt where X runs from 1 to num\_clust 
2. LocalBeds/, All bed files generated by the program are stored here


In the same folder include the following files

1. refGene.txt, this is the refseq file corresponding  to the genome used in the study. For example if using hg_19 download it [here](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refGene.txt.gz)
2. The bedfile with the binding sites
3. A list of genes (allGenes.txt is included here for human genes)


To run the program open an Ipython shell and type

	run TFbind.py [gene_file] [interval] [num_clust] [bedfile] [p-value]

or type

	python TFbind.py [gene_file] [interval] [num_clust] [bedfile] [p-value]

where

	The inputs to the program are

	1) gene_file:  A list of all genes in the genome
	2) Interval: The region around the TSS to consider in base pairs
	3) num_clust: The number of clusters. There must be a corresponding file in Clusters/cluster_X.txt (where X runs from 1 to num_clust) which lists genes in the cluster. 
	4) bedfile: A bedfile which contains the TF binding sites
	5) p-value: The level of significance say 0.05


EXAMPLE
=======

To reproduce the analysis in the paper shown in table 9 for ERalpha follow the following steps:

1. Download the binding sites from the cistrome database [here](http://cistrome.org/NR_Cistrome/Cistrome/ChIP_seq/Human_ER_seq/Human_MCF-7_PolII_17b-E2-1hr_Stunnenberg.bed)
2. Download the cluster files included in the Cluster folder in this repository.


Run the analysis by typing 

	python TFbind.py allGenes.txt 20000 12 Human_MCF-7_ESR1_17b-E2-1hr_Stunnenberg.bed 0.05


The output on the command line is 

	Generating BED file for all genes...
	Done generating BED file.
	Generating BED file for individual clusters...
	Done generating cluster BED files.
	Determining cluster binding proportions and significance...
	1 27 0.000454
	2 31 0.002659
	3 11 0.12291
	4 20 0.007266
	5 15 0.173521
	6 27 0.002978
	7 10 0.691283
	8 32 0.001042
	9 18 0.010105
	10 30 1.9e-05
	11 5 0.765658
	12 19 0.258238
	Done!

For each of the twelve clusters, we get the number of genes with TF binding sites near the TSS and the corresponding p-value. 


	
