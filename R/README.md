# emilio
Collection of scripts and utilities for gene co-expression network analysis

# dependencies


# We will use as example the files 
GSE153345_TIS_counts.txt: count data 
run_selector_moni.txt: metadata from the sra selector, we need to identify biological replicates
metadata.txt: table relating samples with biological conditions

# Explanation
We download the raw count matrix from SRA (accession/bioproject:xxxxxxxxx) wich complies a timeline for several biological conditions
(see experimental desing).
* Data reading and filtering  ok
* Preliminary analysis (PCA)  ok
* Network construction (pearson correlation calculation and formating network) ok 
* Module detection (mcl) ok
* Functional enrichment
* Heatmap plotting ok
* Network statistics

