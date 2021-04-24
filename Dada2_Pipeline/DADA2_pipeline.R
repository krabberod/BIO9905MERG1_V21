#### Load libraries ####

library("dada2")
library("phyloseq")
library("Biostrings")
library("ggplot2")
library("dplyr")
library("tidyr")
library("tibble")
library("readxl")
library("readr")
library("stringr")
library("kableExtra")
library("tidyverse")

#### Prepare Directories ####
# Check your current working directory:
getwd()

# If you want to set a new working directory:
# setwd("~/path_to_my_directory/DADA2_pipeline")


# Define the name of directories to use.
fastq_dir <- "fastq"  # fastq directory with the samples that will be used
database_dir <- "databases/"  # folder with the PR2 database https://github.com/vaulot/metabarcodes_tutorials/tree/master/databases

filtered_dir <- "fastq_filtered/"  # for the fastq-files after filtering
qual_dir <- "qual_pdf/"  # quality scores plots
dada2_dir <- "dada2_results/"  # dada2 results
blast_dir <- "blast/"  # blast2 results

# Create the directories
dir.create(filtered_dir)
dir.create(qual_dir)
dir.create((dada2_dir))
dir.create(blast_dir)

#### Primers
# The primers in this example are for commonly used for 18S

primer_set_fwd = c("CCAGCAGCCGCGGTAATTCC", "CCAGCACCCGCGGTAATTCC", "CCAGCAGCTGCGGTAATTCC",
                   "CCAGCACCTGCGGTAATTCC")
primer_set_rev = c("ACTTTCGTTCTTGATYRATGA")
primer_length_fwd <- str_length(primer_set_fwd[1])
primer_length_rev <- str_length(primer_set_rev[1])


#### Examine fastq files
# get a list of all fastq files in the fastq" directory and separate R1 and R2
fns <- sort(list.files(fastq_dir, full.names = T))
fns <- fns[str_detect(basename(fns), ".fastq.gz")]
fns_R1 <- fns[str_detect(basename(fns), "R1")]
fns_R2 <- fns[str_detect(basename(fns), "R2")]

# Extract sample names, assuming filenames have format: 18S_SAMPLENAME_XXX.fastq.gz
sample.names <- str_split(basename(fns_R1), pattern = "_", simplify = TRUE)
sample.names <- sample.names[, 2]

#### Make a dataframe with the number of sequences in each file ####
# NB! Any sequence of length 0 will cause the loop to crash
# These should be removed in advance

df <- data.frame()

for (i in 1:length(fns_R1)) {

  # use the Biosstrings function fastq.geometry
  geom <- fastq.geometry(fns_R1[i])

  # extract the information on number of sequences and file name
  df_one_row <- data.frame(n_seq = geom[1], file_name = basename(fns_R1[i]))

  # add one line to data frame
  df <- bind_rows(df, df_one_row)
}

# knitr::kable(df) # to make html table.
View(df)

# If you want to write the table to your working directory remove the hashtag and use:
# write.table(df, file = 'n_seq.txt', sep='\t', row.names = FALSE, na='',quote=FALSE)


# plot the histogram with number of sequences
# The plot for the example data looks kind of uninformative, why?

ggplot(df, aes(x = n_seq)) +
  geom_histogram(alpha = 0.5, position = "identity", binwidth = 15000)

hist(df$n_seq, breaks = 10)

#### Plot Quality for each fastq file ####
for (i in 1:length(fns)) {

  # Use dada2 function to plot quality
  p1 <- plotQualityProfile(fns[i])

  # Only plot on screen for first 2 files
  if (i <= 2) {
    print(p1)
  }

  # save the file as a pdf file 
  p1_file <- paste0(qual_dir, basename(fns[i]), ".qual.pdf")

  ggsave(plot = p1, filename = p1_file, device = "pdf", width = 15, height = 15,
         scale = 1, units = "cm")
}



####
# Prepare the outputnames for filtered reads:
filt_R1 <- str_c(filtered_dir, sample.names, "_R1_filt.fastq")
filt_R2 <- str_c(filtered_dir, sample.names, "_R2_filt.fastq")


####FILTER WITH CUTADAPT (not executed inside R)
# Cutadapt is the preferred way of filtering sequences, because it will trim
# primers allowing some mismatches and ambiguities. It is not implemented
# in R at the moment. But you should have learned how to do prior to using this pipeline. 
# Why is this not an optimal solution?


# Alternatively install the following package and do it in R (not tested, as of april 2021)
# library(devtools)
# install_github("omicsCore/SEQprocess")
# cutadpat(fq1, output.dir, adpat.seq="insert primer sequence here", m=1, mc.cores=1, run.cmd=TRUE)


#### SIMPLE FILTER BY LENGTH OF PRIMER
# This is is a workaround if you can't filter by any other means. However, it is highly 
# advisable to remove the primers with for instance cutadatp. 

out <- filterAndTrim(fns_R1, filt_R1, fns_R2, filt_R2, truncLen = c(250, 200),
                     trimLeft = c(primer_length_fwd, primer_length_rev), maxN = 0,
                     maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
                     compress = FALSE, multithread = FALSE)

#### DADA2 ####
# If your setup allows running multiple threads set multithread = TRUE
# Windows: multithread = FALSE
# Can take a couple of mminutes: 

err_R1 <- learnErrors(filt_R1, multithread = TRUE)
plotErrors(err_R1, nominalQ = TRUE)

err_R2 <- learnErrors(filt_R2, multithread = T)
plotErrors(err_R2, nominalQ = TRUE)


#### Dereplicate the reads ####

derep_R1 <- derepFastq(filt_R1, verbose = FALSE)
derep_R2 <- derepFastq(filt_R2, verbose = FALSE)

# Name the derep-class objects by the sample names
names(derep_R1) <- sample.names
names(derep_R2) <- sample.names

####  Sequence-variant inference algorithm on the dereplicated data ####
# If your computer can run multiple threads set multithread = TRUE

ptm <- proc.time()
dada_R1 <- dada(derep_R1, err = err_R1, multithread = TRUE, pool = FALSE)
dada_mulit<-proc.time() - ptm
dada_R2 <- dada(derep_R2, err = err_R2, multithread = T, pool = FALSE)

# Viewing the first entry in each of the dada objects
dada_R1[[1]]
dada_R2[[1]]

#### Merge Sequences
mergers <- mergePairs(dada_R1, derep_R1, dada_R2, derep_R2, verbose = TRUE)
head(mergers[[1]])

####   ####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Make a transposed version of seqtab to make it similar to data in mothur
t_seqtab <- t(seqtab) # the function t() is a simple transposing of the matrix
table(nchar(getSequences(seqtab)))
plot(table(nchar(getSequences(seqtab)))) #simple plot of length distribution


#### Remove chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE,
                                    verbose = TRUE)

# Compute % of non chimeras
paste0("% of non chimeras : ", sum(seqtab.nochim)/sum(seqtab) * 100)
paste0("total number of sequences : ", sum(seqtab.nochim))

# What are chimeras, and why do we remove them? How is it done in Dada2?

#### Track number of reads at each step
getN <- function(x) sum(getUniques(x)) # example of a function in R

track <- cbind(out, sapply(dada_R1, getN), sapply(mergers, getN), rowSums(seqtab),
               rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names

#View the output
track

#### Transforming and saving the OTU sequences

seqtab.nochim_trans <- as.data.frame(t(seqtab.nochim)) %>% rownames_to_column(var = "sequence") %>%
  rowid_to_column(var = "OTUNumber") %>% mutate(OTUNumber = sprintf("OTU_%05d",
                                                                    OTUNumber)) %>% mutate(sequence = str_replace_all(sequence, "(-|\\.)", ""))

df <- seqtab.nochim_trans
seq_out <- Biostrings::DNAStringSet(df$sequence)
names(seq_out) <- df$OTUNumber
seq_out

Biostrings::writeXStringSet(seq_out, str_c(dada2_dir, "OTU_no_taxonomy.fasta"),
                            compress = FALSE, width = 20000)
#### Assinging taxonomy
# The PR2 database can be found here:
# https://pr2-database.org/
# The PR2 database has preformatted files suitable for dada2  for both 18S and 16S
#
# This step depends on the kind of taxonomic assignment that will be used later
# The PR2 database is a curated quality database for protissts, with 8 taxonomic ranks

PR2_tax_levels <- c("Kingdom", "Supergroup", "Division", "Class",
                    "Order", "Family",
                    "Genus", "Species")

pr2_file <- paste0("databases/pr2_version_4.13.0_18S_dada2.fasta.gz")

# OBS! The next step takes a long time. ~45-60 min on a medium fast computer...
# So in case we are running late skip this next command. If we have time
# start the process (i.e. remove hashtags), and have some coffee/lunch.

# ptm <- proc.time()
# taxa <- assignTaxonomy(seqtab.nochim, refFasta = pr2_file, taxLevels = PR2_tax_levels,
#                       minBoot = 0, outputBootstraps = TRUE, verbose = TRUE)
# taxo_assign<-proc.time() - ptm

### Runtime result:
# user   system  elapsed
# 6612.075  504.603 2576.869

# Single objects can be saved to a file with the saveRDS() function.
# saveRDS(taxa, str_c(dada2_dir, "dada2.taxa.rds"))
# I have prepared a taxonomy file that I can put on github, if necessary.
# taxa <- readRDS(str_c(dada2_dir, "taxa.rds"))
# Seqtab.nochim_trans <- read.RDS(str_c("seqtab.nochim_trans.rds"))

# Export information in tab or comma separated files
# Tab:
write_tsv(as_tibble(taxa$tax), file = str_c(dada2_dir, "taxa.txt"))

# Csv:
#write.csv(taxa$tax, file = str_c(dada2_dir, "taxa.txt"))

#### Appending taxonomy and boot to the sequence table ####
taxa_tax <- as.data.frame(taxa$tax)
taxa_boot <- as.data.frame(taxa$boot) %>% rename_all(funs(str_c(., "_boot")))
seqtab.nochim_trans <- taxa_tax %>% bind_cols(taxa_boot) %>% bind_cols(seqtab.nochim_trans)


#### Filter for 18S ####
# Define a minimum bootstrap value for filtering
# Think before applying the cut-off! What is the benefits of removing
# OTUs with low support? What are the drawbacks?
bootstrap_min <- 80

# Remove OTU with annotation below the bootstrap value
seqtab.nochim_18S <- seqtab.nochim_trans %>% dplyr::filter(Supergroup_boot >=
                                                             bootstrap_min)
seqtab.nochim_18S_lowsupport<- seqtab.nochim_trans %>% dplyr::filter(Supergroup_boot <=
                                                                       bootstrap_min)
write_tsv(seqtab.nochim_18S, str_c(dada2_dir, "OTU_table.tsv"))
write_tsv(seqtab.nochim_18S_lowsupport, str_c(dada2_dir, "OTU_table_lowsupport.tsv"))

####  Write FASTA file for BLAST analysis with taxonomy ####
# Blasting is an alternative to RDP classifier:

df <- seqtab.nochim_trans
seq_out <- Biostrings::DNAStringSet(df$sequence)

names(seq_out) <- str_c(df$OTUNumber, df$Supergroup, df$Division, df$Class,
                        df$Order, df$Family, df$Genus, df$Species, df$Species_boot1, sep = "|")

Biostrings::writeXStringSet(seq_out, str_c(blast_dir, "OTU.fasta"), compress = FALSE,
                            width = 20000)

#### EXTRA Example for blast on the cluster Saga ####
##!/bin/sh
##SBATCH --job-name=blastn
##SBATCH --account= #add your own project
##SBATCH --output=slurm-%j.base
##SBATCH --cpus-per-task=16
##SBATCH --time=100:00:00
##SBATCH --mem-per-cpu=6G
#
#module purge
#module load BLAST+/2.8.1-intel-2018b
#
#FASTA=OTU.fasta
#BLAST_TSV=OTU.blast.tsv
#DB=/cluster/shared/databases/blast/latest/nt
#
#
#OUT_FMT="6 qseqid sseqid sacc stitle sscinames staxids sskingdoms sblastnames pident slen length mismatch gapopen qstart qend sstart send evalue bitscore"
#
#blastn -max_target_seqs 100 -evalue 1.00e-10 -query $FASTA -out $BLAST_TSV -db "$DB" -outfmt "$OUT_FMT" -num_threads 16
###############

#### Make Phyloseq Object ####
samdf <- data.frame(sample_name = sample.names)
rownames(samdf) <- sample.names

rownames(seqtab.nochim_18S)<-seqtab.nochim_18S$OTUNumber


OTU <- seqtab.nochim_18S %>% select_if(is.numeric) %>%
  select(-contains("_boot")) %>% as.matrix() %>% otu_table(taxa_are_rows = TRUE)

TAX <- seqtab.nochim_18S %>% select(Kingdom:Species) %>%
  as.matrix() %>% tax_table()

ps_dada2 <- phyloseq(OTU, sample_data(samdf), TAX)

### Saving and loading data ####
# You can save selected objects:
saveRDS(ps_dada2, str_c(dada2_dir, "phyloseq.rds"))
# ps_dada2<-readRDS(str_c(dada2_dir, "phyloseq.rds"))

# Or save the entire workspace:
# save.image("dada2.RData")

# Which can be loaded with:
# Can be loaded with
# load("dada2.RData"")
