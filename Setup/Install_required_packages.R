# Please make sure that you have the required R packages installed.
# You can get a list of the package already installed on your computer s by executing

installed.packages()

# in R alternatively you can just run the installation for each package
# to make sure that you have the latest version.
# Run these installation commands line-by-line in R (or Rstudio)
# and answer yesif you are asked to update any previously installed pakages:


install.packages("readr")     # To read and write files
install.packages("readxl")    # To read excel files
install.packages("dplyr")     # To manipulate dataframes
install.packages("tibble")    # To work with data frames
install.packages("tidyr")     # To work with data frames
install.packages("stringr")   # To manipulate strings
install.packages("ggplot2")   # To do plots

install.packages("Matrix")	# Manipulation of large matrices
install.packages("igraph")	# Network analysis
install.packages("standardize")
install.packages("vegan")
#install.packages("cowplot")

install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager")
BiocManager::install(c("dada2", "phyloseq","Biostrings"))
BiocManager::install("microbiome")
#BiocManager::install("SpiecEasi") #Network construction

library(devtools)
install_github("zdk123/SpiecEasi")
