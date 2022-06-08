# BF528

This repository consists of various scripts that were generated for BF528 Applications in Translational Bioinformatics class projects. The idea was to replicate the results of the given papers and report the findings. Since it was a group project, I assigned a particular role for each project which kept on rotating to get the experience of all roles. I was responsible for replicating the results for a particular portion of the project. 

## Project 1

Paper: Marisa et al. Gene Expression Classification of Colon Cancer into Molecular Subtypes: Characterization, Validation, and Prognostic Value. PLoS Medicine, May 2013. PMID: 23700391

### Role: Biologist



## Project 2

Paper: O’Meara et al. Transcriptional Reversion of Cardiac Myocyte Fate During Mammalian Cardiac Regeneration. Circ Res. Feb 2015. PMID: 25477501

### Role: Analyst 

## Project 3

Paper: Wang, Charles, Binsheng Gong, Pierre R. Bushel, Jean Thierry-Mieg, Danielle Thierry-Mieg, Joshua Xu, Hong Fang, et al. 2014. “A comprehensive study design reveals treatment- and transcript abundance–dependent concordance between RNA-seq and microarray data” Nature Biotechnology 32 (9): 926–32. PMID: 4243706

### Role: Data Curator

#### STAR.qsub

This script aligns the paired end reads to the reference genome where input is in the form of .fastq files and output as .bam files with alignment statistics.

* Dependencies: STAR Aligner

* Input: There are 3 positional arguments out of which 1st and 2nd take the 1st and 2nd read of paired end respectively and GENOMEDIR in the script specifies the path of reference genome to which the reads are aligned.

* Execution:`qsub STAR.qsub SRR1177XXX_1.fastq.gz SRR1177XXX_2.fastq.gz`(where XXX depends on the sample number)

* Output: The 3rd positional argument defines the name of the output file i.e. star_output in which the outputs are stored as .bam files.

### multiqc.qsub

This script reports the summary statistics using the FastQC and STAR alignment results generated earlier.

* Dependencies: multiqc

* Input: STAR alignment output (.bam files) and fastq files

* Excecution: `qsub multiqc.qsub`

* Output: Multiqc output is in the form of an HTML report which is displayed in the browser (multiqc_report.html)


## Project 4

Paper: Baron, Maayan, Adrian Veres, Samuel L. Wolock, Aubrey L. Faust, Renaud Gaujoux, Amedeo Vetere, Jennifer Hyoje Ryu, et al. 2016. “A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-Cell Population Structure.” Cell Systems 3 (4): 346–60.e4. PMID: 27667365

### Role: Programmer

#### programmer.R

* This script is for filtering, normalization and clustering of genes.

* Input: quants_mat.gz which is output from ALevin

* Dependencies: R

* Execution: It is recommended to run this script in Rstudio by installing the following packages: Seurat, biomaRt, GenomicFeatures, tidyverse, Matrix, tximport, SeqGSEA, fishpond, EnsDb.Hsapiens.v79

   Alternatively to run on command line:

            module load R/4.0.2

            Rscript programmer.R

* Output: panc_cells.rda file

