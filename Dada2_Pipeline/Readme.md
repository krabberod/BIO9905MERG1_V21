# Dada2 pipeline.

This folder contains the necessary files and data to run the DADA2 pipeline in R. For this course we will mainly use Rstudio when running , which you already should have installed R. ALso make sure that you have installed the necessary packages.

- Make a working directory for the DADA2 pipeline
- Download the Rscript called [DADA2_pipeline.R](DADA2_pipeline.R)
- Make as subdirectory called *fastq* in the working directroy
- Download the fastq files to the subdirectory
- Open the script in Rstudio

**The fastq files have already been processed:**
- The samples have been demultiplexed, i.e. split into individual per-sample fastq files.
- Primers and adapters and have been removed.
- Paired-end data has been matched in the same order (i.e. sequences in the *R1* and *R2* fastq files arein the same order)


The script is based on the DADA2 tutorial https://benjjneb.github.io/dada2/tutorial.html
