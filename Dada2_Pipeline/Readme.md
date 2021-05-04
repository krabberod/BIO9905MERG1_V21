# Dada2 pipeline.

This folder contains the necessary files and data to run the DADA2 pipeline in R. For this course we will mainly use Rstudio when running , which you already should have installed R. Also make sure that you have installed the necessary packages.

- The [introduction Lecture](./../Lectures_and_groups/DADA2_lecture.pdf)

### Hands on run-through:
- Make a working directory for the DADA2 pipeline on your computer
- Optional: Make an R-project in your working directory
- Download the R script either from Rstudio (type in the console):
```download.file("https://raw.githubusercontent.com/krabberod/BIO9905MERG1_V21/main/Dada2_Pipeline/DADA2_pipeline.R", "DADA2_pipeline.R")```
- Or the boring way [DADA2_pipeline.R](https://raw.githubusercontent.com/krabberod/BIO9905MERG1_V21/main/Dada2_Pipeline/DADA2_pipeline.R)
- Open the script in Rstudio
- Make as subdirectory called *fastq* in the working directroy
- Download the fastq files to the subdirectory (link is in the script)

**The fastq files have already been processed:**
- The samples have been demultiplexed, i.e. split into individual per-sample fastq files.
- Primers and adapters and have been removed.
- Paired-end data has been matched in the same order (i.e. sequences in the *R1* and *R2* fastq files arein the same order)
- A link can be found in the DAD2 script

- If you didn't manage to finish the pipeline during hands-on an RData image of the full pipeline has been added to the repository containing all the necessary results from the pipeline: [dada2.RData](dada2.RData)

This can be downloaded and the opened i Rstudio.
load("dada2.RData")

The script is based on the DADA2 tutorial https://benjjneb.github.io/dada2/tutorial.html
