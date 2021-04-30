## BIO9905MERG1_V21
## Ramiro Logares (ICM, CSIC, Barcelona)
## 2021

######################
## Community ecology
######################

# Install packages (in case you didn't before)
install.packages("vegan")          # Installs the community ecology package Vegan with hundreds of functions
install.packages("readr")          # To read and write files
install.packages("readxl")         # To read excel files
install.packages("dplyr")          # To manipulate dataframes
install.packages("tibble")         # To work with data frames
install.packages("tidyr")          # To work with data frames
install.packages("stringr")        # To manipulate strings
install.packages("ggplot2")        # To do plots
install.packages("kableExtra")     # Nice table formatting with knitr
install.packages("tidyverse")      # Collection of packages for data science
install.packages("spaa")           # Installs the ecological package spaa
install.packages("compositions")   # To work with compositional data
install.packages("zCompositions")  # To work with compositional data
install.packages("devtools")       # Developer tools
install.packages("mixOmics")       # Multivariate methods
install.packages("ape")            # Phylogenetic tools
install.packages("recluster")      # Clustering tools
install.packages("dendextend")     # To work with dendrograms
install.packages("corrplot")       # makes nice correlation plots
install.packages("RcmdrMisc")      # diverse tools

# We install other packages with another method
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install(version = "3.10")
BiocManager::install(c("dada2", "phyloseq","Biostrings", "PCAtools"))

# We install other packages with another method
devtools::install_github("pr2database/pr2database") # Installs directly from github resources that are not in R repos
devtools::install_github("GuillemSalazar/EcolUtils") # Installs other tools for ecological analyses
devtools::install_github('fawda123/ggord')

#Load packages
library(vegan)
library(spaa)
library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(readxl)
library(readr)
library(stringr)
library(kableExtra) 
library(tidyverse)
#library("pr2database")
library(compositions)
library(zCompositions)
library(PCAtools)
library(mixOmics)
library(ape)
library(recluster)
library(dendextend)
library(corrplot)
library(RcmdrMisc)
library(devtools)
library(ggord)
library(ggplot2)
# Perhaps other packages may be needed, and we'll install them as we go


## Starting community ecology analyses

# Read dada2 otuput
otu.tab<-read_tsv("https://raw.githubusercontent.com/krabberod/BIO9905MERG1_V21/main/Dada2_Pipeline/dada2_results/OTU_table.tsv")
dim(otu.tab) # 2107   26
#Let's reorder the table
otu.tab<-otu.tab[,c(17,19:26,1:16,18)]
#We assign to rownames the OTU names
otu.tab <- column_to_rownames(otu.tab, var = "OTUNumber") # %>% as_tibble()
rownames(otu.tab) # Let's check the names
dim(otu.tab) # 2107   25
otu.tab.simple<-otu.tab[,1:8] # We'll need this table for community ecology analyses
#We transpose the table, as this is how Vegan likes it
otu.tab.simple<-t(otu.tab.simple)
otu.tab.simple[1:5,1:5]

#         OTU_00001 OTU_00002 OTU_00004 OTU_00005 OTU_00006
# BL040126      4996     12348     11426         0      3958
# BL040419       739       684        97     16605      4702
# BL040719         0         0       166         0       806
# BL041019        78        74         0       184       286
# BL050120     30697     12885      5417         0      3739

# Richness estimations
richness<-estimateR(otu.tab.simple)

#            BL040126   BL040419   BL040719   BL041019   BL050120  BL050413  BL050705   BL051004
# S.obs    642.000000 499.000000 263.000000 414.000000 942.000000 89.000000 69.000000 227.000000
# S.chao1  642.000000 499.000000 263.000000 414.000000 943.250000 89.000000 69.000000 227.000000
# se.chao1   0.000000   0.000000   0.000000   0.000000   1.621617  0.000000  0.000000   0.000000
# S.ACE    642.000000 499.000000 263.000000 414.000000 943.399653 89.000000 69.000000 227.000000
# se.ACE     7.091415   9.573887   4.198497   6.011262  10.391615  2.539574  2.797514   2.604638

# Above we have the estimators Chao and ACE as well as the species number.

# Rarefaction

#Let's calculate the number of reads per sample
rowSums(otu.tab.simple)

# BL040126 BL040419 BL040719 BL041019 BL050120 BL050413 BL050705 BL051004 
#   182462    66827    55896    39672   189636    29053    10771    87192 

rarecurve (otu.tab.simple, step=100, xlab= "Number of reads", ylab="Richness", col="blue")

#Accumulation curves
accum.curve<-specaccum(otu.tab.simple, method="collector")
plot(accum.curve)

#Evenness
plot(colSums(otu.tab.simple),log="y",xlab="Rank", ylab="Abundance", pch=19, cex=0.5, col="blue")

#Fitting rank-abundance distribution models to the data
mod<-radfit(otu.tab.simple)
plot(mod)

mod.all<-radfit(colSums(otu.tab.simple))
plot(mod.all)

#Fitting data to the Preston model
preston<-prestonfit(colSums(otu.tab.simple))
preston.dist<-prestondistr(colSums(otu.tab.simple))
plot(preston)
lines(preston.dist, line.col="blue3")

## Extrapolated richness
veiledspec(preston)
# Extrapolated     Observed       Veiled 
# 2113.475329  2107.000000     6.475329 

veiledspec(preston.dist)
# Extrapolated     Observed       Veiled 
# 2113.236021  2107.000000     6.236021 


#Shannon H index (considers richness and evenness)

H<-diversity(otu.tab.simple, index="shannon")
# BL040126 BL040419 BL040719 BL041019 BL050120 BL050413 BL050705 BL051004 
# 5.049747 4.185494 4.627698 5.236017 4.849669 3.698185 3.406164 4.358232 

plot(H, type="l", col="blue")

#Pielou's index of evenness (range 0-1, 1 = maximum evenness)
# J=H/Hmax
# J=Shannon (H) / log(S=species richness)

J=H/log(rowSums(otu.tab.simple>0))
#  BL040126  BL040419  BL040719  BL041019  BL050120  BL050413  BL050705  BL051004 
# 0.7811398 0.6737098 0.8305043 0.8689236 0.7081871 0.8238995 0.8044587 0.8033681 

# Inverse Simpson's D index (richness+evenness. Larger values, larger diversity)
inv.simpson<-diversity(otu.tab.simple, "invsimpson")
plot(inv.simpson, type="l", col="blue")

# BL040126  BL040419  BL040719  BL041019  BL050120  BL050413  BL050705  BL051004 
# 54.13768  13.15796  48.69382 107.30411  26.16040  25.27907  19.71550  37.93128 

#################
## Beta diversity
#################

#We rarefy all samples to the same sequencing depth, to reduce biases
min(rowSums(otu.tab.simple)) # We calculate the sample with the minimum amount of reads
# [1] 10771

otu.tab.simple.ss<-rrarefy(otu.tab.simple, 10771) #Samples are rarefied to 10771 reads per sample
rowSums(otu.tab.simple.ss) # We check the number of reads per sample
# BL040126 BL040419 BL040719 BL041019 BL050120 BL050413 BL050705 BL051004 
#  10771    10771    10771    10771    10771    10771    10771    10771 

#Check the dimensions of the tables
dim(otu.tab.simple) 
# [1]    8 2107
dim(otu.tab.simple.ss)
# [1]    8 2107

#Tables have the same size, but, after removing reads, several OTUs are left with zero abundances
length(which(colSums(otu.tab.simple)==0))
# [1] 0 #No OTU has an abundance sum that is 0, as expected

length(which(colSums(otu.tab.simple.ss)==0))
#  [1] 273 # A total of 273 OTUs were found in the rarefied table with zero abundance. Let's corroborate

which(colSums(otu.tab.simple.ss)==0) # Show the OTUs and the position in the table that have 0 abundance
# A small subsample of them
# OTU_00814 OTU_01076 OTU_01077 OTU_01232 OTU_01242 
#    772      1020      1021      1166      1176     

otu.tab.simple[,772] # This gives the abundance of the OTU_00814 across the different samples in the table that is NOT subsampled
# BL040126 BL040419 BL040719 BL041019 BL050120 BL050413 BL050705 BL051004 
#    0        0        0        0       88        0        0        0 

otu.tab.simple.ss[,772] # # This gives the abundance of the OTU_00814 across the different samples in the table that IS subsampled
# BL040126 BL040419 BL040719 BL041019 BL050120 BL050413 BL050705 BL051004
# 0        0        0        0        0        0        0        0 

otu.tab.simple.ss.nozero<-otu.tab.simple.ss[,-(which(colSums(otu.tab.simple.ss)==0))] # Removes OTUs with cero abundance
length(which(colSums(otu.tab.simple.ss.nozero)==0)) # Check that no cero abundance OTUs are left
# [1] 0 # correct
# Let's check dimensions
dim(otu.tab.simple.ss)
# [1]    8 2107
dim(otu.tab.simple.ss.nozero)
# [1]    8 1834
# 2107-1834 = 273 , This is the number of OTUs that we expected to be removed.

## Compositional data analyses
# We install packages to work with compositional data
install.packages("compositions")
install.packages("zCompositions")
library(compositions)
library(zCompositions)

otu.tab.simple.gbm<-cmultRepl(t(otu.tab.simple), output = "p-counts")  # replace zeros (problems with log calculations) with pseudo-counts
# No. corrected values:  12246 
otu.tab.simple.gbm[1:5,1:5] # We have a look to the replaced values
#               BL040126 BL040419    BL040719    BL041019     BL050120
# OTU_00001  4996.000000      739   0.9810100  78.0000000 30697.000000
# OTU_00002 12348.000000      684   0.9744656  74.0000000 12885.000000
# OTU_00004 11426.000000       97 166.0000000   0.6851427  5417.000000
# OTU_00005     3.229364    16605   0.9892938 184.0000000     3.356335
# OTU_00006  3958.000000     4702 806.0000000 286.0000000  3739.000000

# centered log-ratio (clr) transformation
otu.tab.simple.gbm.clr<-clr(otu.tab.simple.gbm) # We apply a centered log-ratio (clr) transformation 
otu.tab.simple.gbm.clr[1:5,1:5] #Values now look different than counts.
# clr values indicate how OTUs behave relative to the per-sample average
#            BL040126   BL040419   BL040719    BL041019  BL050120
# OTU_00001  3.016847 1.10575200 -5.5187186 -1.14283710  4.832374
# OTU_00002  5.034361 2.14106914 -4.4127548 -0.08282368  5.076930
# OTU_00004  4.818212 0.04927624  0.5865531 -4.90356291  4.071863
# OTU_00005 -2.082582 6.46259197 -3.2656311  1.96006859 -2.044017
# OTU_00006  3.237162 3.40941121  1.6457517  0.60965979  3.180241


# Distance metrics
# We calculate the Bray Curtis dissimilarities for the rarefied dataset
otu.tab.simple.ss.nozero.bray<-vegdist(otu.tab.simple.ss.nozero, method="bray")
as.matrix(otu.tab.simple.ss.nozero.bray)[1:5,1:5]
#            BL040126  BL040419  BL040719  BL041019  BL050120
# BL040126 0.0000000 0.8087457 0.9264692 0.8720639 0.5661498
# BL040419 0.8087457 0.0000000 0.9017733 0.8754062 0.8352985
# BL040719 0.9264692 0.9017733 0.0000000 0.7490484 0.9118002
# BL041019 0.8720639 0.8754062 0.7490484 0.0000000 0.8183084
# BL050120 0.5661498 0.8352985 0.9118002 0.8183084 0.0000000


# We calculate the Euclidean distance based on the clr data (also known as Aitchison distance)
otu.tab.simple.gbm.clr.euclidean<-dist(t(otu.tab.simple.gbm.clr), method = "euclidean") 
as.matrix(otu.tab.simple.gbm.clr.euclidean)[1:5,1:5]
#          BL040126 BL040419 BL040719 BL041019 BL050120
# BL040126   0.0000 162.3852 153.9187 178.3993 176.3504
# BL040419 162.3852   0.0000 142.0086 164.0183 198.7059
# BL040719 153.9187 142.0086   0.0000 134.3082 182.8134
# BL041019 178.3993 164.0183 134.3082   0.0000 205.9812
# BL050120 176.3504 198.7059 182.8134 205.9812   0.0000

#Let's compare the distance matrices
identical(rownames(as.matrix(otu.tab.simple.ss.nozero.bray)),rownames(as.matrix(otu.tab.simple.gbm.clr.euclidean))) # Check matrices have same order
# [1] TRUE

plot(otu.tab.simple.ss.nozero.bray, otu.tab.simple.gbm.clr.euclidean, pch=19, xlab="Bray Curtis", ylab="Aitchison") # We generate a simple x-y plot
lm<-lm(otu.tab.simple.gbm.clr.euclidean~otu.tab.simple.ss.nozero.bray) # Fitting a linear model (regression)
abline(lm, col="red")

mantel(otu.tab.simple.ss.nozero.bray, otu.tab.simple.gbm.clr.euclidean) # The correlation between distance matrices is tested with a Mantel test
# Mantel statistic based on Pearson's product-moment correlation 

# Call:
# mantel(xdis = otu.tab.simple.ss.nozero.bray, ydis = otu.tab.simple.gbm.clr.euclidean) 

# Mantel statistic r: 0.06774 
#       Significance: 0.346  # Correlations between distances matrices is not significant

# Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
# 0.280 0.355 0.419 0.462 
# Permutation: free
# Number of permutations: 999

#Ordination and clustering

#PCA

# We install PCAtools
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('PCAtools')

library(PCAtools)

#PCA rarefied table
otu.tab.simple.ss.nozero.pca<-pca(t(otu.tab.simple.ss.nozero), scale=FALSE) # Runs de PCA
biplot(otu.tab.simple.ss.nozero.pca, showLoadings = T, lab=rownames(otu.tab.simple.ss.nozero)) # Plots de PCA
screeplot(otu.tab.simple.ss.nozero.pca, axisLabSize = 18, titleLabSize = 22) # We plot the percentage of variance explained by each axis

#We install mixOmics
install.packages("mixOmics") # We change the package mixOmics, as PCAtools had some issues with clr tables
library(mixOmics)
#PCA clr table (calculated with the vegan "rda" function, as pca from PCAtools gives errors
otu.tab.simple.gbm.clr.pca<-pca(otu.tab.simple.gbm.clr, scale=FALSE, ncomp=6) # NB: the pca used here is from "mixOmics" while the pca above is from "PCAtools"
plotVar(otu.tab.simple.gbm.clr.pca)
plot(otu.tab.simple.gbm.clr.pca)


#PCoA

#PCoA not included in Vegan, so we install the APE package
install.packages("ape")
library(ape)
# Calculates the PCoA of the rarefied table using Bray Curtis distances
otu.tab.simple.ss.nozero.bray.pcoa<-pcoa(otu.tab.simple.ss.nozero.bray) # Calculates PCoA
barplot(otu.tab.simple.ss.nozero.bray.pcoa$values$Relative_eig[1:10]) # Plot Eigenvalues (amount of variance explained by Axes)
biplot.pcoa(otu.tab.simple.ss.nozero.bray.pcoa)

# Calculates the PCoA of the clr table using Euclidean distances
otu.tab.simple.gbm.clr.euclidean.pcoa<-pcoa(otu.tab.simple.gbm.clr.euclidean) # Calculates PCoA
barplot(otu.tab.simple.gbm.clr.euclidean.pcoa$values$Relative_eig[1:10]) # Plot Eigenvalues (amount of variance explained by Axes)
biplot.pcoa(otu.tab.simple.gbm.clr.euclidean.pcoa)

#NMDS
# We will define the function NMDS.scree() that automatically performs a NMDS for 1-7 dimensions 
# and plots the number of dimensions vs. stress

set.seed(666) # We include this value to make results reproducible
NMDS.scree <- function(x) { # x is the name of the distance matrix
  plot(rep(1, 7), replicate(7, metaMDS(x, autotransform = F, k = 1)$stress), xlim = c(1, 7),ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")
  for (i in 1:7) {
    points(rep(i + 1,7),replicate(7, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

# Using the function to determine the optimal number of dimensions
# Using the rarefied table
NMDS.scree(otu.tab.simple.ss.nozero.bray)
# Using the clr table
NMDS.scree(otu.tab.simple.gbm.clr.euclidean)

# We calculate NMDS for k(dimensions)=2
# Rarefied table (we use the dataframe to have access to sample and OTU names)
otu.tab.simple.ss.nozero.bray.nmds<-metaMDS(otu.tab.simple.ss.nozero, k=2, trymax=100, trace=F, autotransform = F, distance="bray")
stressplot(otu.tab.simple.ss.nozero.bray.nmds) # Make stressplot

# clr table (we use the dataframe to have access to sample and OTU names)
otu.tab.simple.gbm.clr.euclidean.nmds<-metaMDS(t(as.data.frame(otu.tab.simple.gbm.clr)), k=2, trymax=100, trace=F, autotransform = F, distance="euclidean")
stressplot(otu.tab.simple.gbm.clr.euclidean.nmds) # Make stressplot

# Simple plotting
# Rarefied table 
plot(otu.tab.simple.ss.nozero.bray.nmds, display="sites", type="n")
points(otu.tab.simple.ss.nozero.bray.nmds, display = "sites", col = "red", pch=19)
text(otu.tab.simple.ss.nozero.bray.nmds, display ="sites")

# clr table 
plot(otu.tab.simple.gbm.clr.euclidean.nmds, display="sites", type="n")
points(otu.tab.simple.gbm.clr.euclidean.nmds, display = "sites", col = "red", pch=19)
text(otu.tab.simple.gbm.clr.euclidean.nmds, display ="sites")

# Let's make nicer plots
# We define seasons for samples
seasons<-c("Winter","Spring","Summer","Autumn","Winter","Spring","Summer","Autumn")
months<-c("January","April","July","October","January","April","July","October")

library(ggplot2) # Generates nice plots
library(ggrepel) # Adds in to ggplot

# Rarefied table 
# We generate a table of nmds scores and other features
otu.tab.simple.ss.nozero.bray.nmds.scores<-as.data.frame(scores(otu.tab.simple.ss.nozero.bray.nmds))
otu.tab.simple.ss.nozero.bray.nmds.scores$seasons<-seasons
otu.tab.simple.ss.nozero.bray.nmds.scores$months<-months
otu.tab.simple.ss.nozero.bray.nmds.scores$samples<-rownames(otu.tab.simple.ss.nozero.bray.nmds.scores)

#                 NMDS1       NMDS2 seasons  months  samples
# BL040126 -0.192087931 -0.34552707  Winter January BL040126
# BL040419  0.163687487  0.01138097  Spring   April BL040419
# BL040719 -0.293448084  0.13565597  Summer    July BL040719
# BL041019 -0.284857321  0.13150682  Autumn October BL041019
# BL050120 -0.209189049 -0.34417159  Winter January BL050120
# BL050413 -0.009003643  0.52375809  Spring   April BL050413
# BL050705  0.652757387 -0.08086158  Summer    July BL050705
# BL051004  0.172141153 -0.03174161  Autumn October BL051004


# Create the plot
p <- ggplot(otu.tab.simple.ss.nozero.bray.nmds.scores) +
  geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour = seasons), size=3) +
  coord_fixed()+## need aspect ratio of 1!
  geom_text_repel(box.padding = 0.5, aes(x = NMDS1, y = NMDS2, label = samples),
            size = 3)
p

# clr table 
# We generate a table of nmds scores and other features
otu.tab.simple.gbm.clr.euclidean.nmds.scores<-as.data.frame(scores(otu.tab.simple.gbm.clr.euclidean.nmds))
otu.tab.simple.gbm.clr.euclidean.nmds.scores$seasons<-seasons
otu.tab.simple.gbm.clr.euclidean.nmds.scores$months<-months
otu.tab.simple.gbm.clr.euclidean.nmds.scores$samples<-rownames(otu.tab.simple.gbm.clr.euclidean.nmds.scores)

#                 NMDS1      NMDS2 seasons  months  samples
# BL040126  -77.250023 -44.586424  Winter January BL040126
# BL040419    8.589223 -54.666069  Spring   April BL040419
# BL040719   -3.408569  13.911075  Summer    July BL040719
# BL041019   23.754238  59.015021  Autumn October BL041019
# BL050120 -143.136183  31.763211  Winter January BL050120
# BL050413   73.779344   9.957496  Spring   April BL050413
# BL050705   91.926317  -1.829225  Summer    July BL050705
# BL051004   25.745653 -13.565085  Autumn October BL051004

p <- ggplot(otu.tab.simple.gbm.clr.euclidean.nmds.scores) +
  geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour = seasons), size=3) +
  coord_fixed()+## need aspect ratio of 1!
  geom_text_repel(box.padding = 0.5, aes(x = NMDS1, y = NMDS2, label = samples),
                  size = 3)
p


#Clustering of samples

# Allows determining the similarity between samples as well as the organization of samples in groups.
# Hierarchical clustering: samples will be organized in ranks according to their similarity and all samples will be included in a large group
# Unweighted Pair-Group Method Using Arithmetic Averages (UPGMA): This linkage method will link samples by considering their distance
#   to a subgroup arithmetic average. This is a method widely used in ecology 

install.packages("recluster")
library("recluster")

#UPGMA

# Rarefied dataset 
# We generate 100 trees by resampling and then, we use the consensus
otu.tab.simple.ss.nozero.bray.upgma<-recluster.cons(otu.tab.simple.ss.nozero.bray, tr=100, p=0.5, method="average") 
plot(otu.tab.simple.ss.nozero.bray.upgma$cons) # plot consensus tree
# We'll calculate bootstrap support values (0: bad - 100: perfect)
# This allows us to know how well supported is the branching pattern
otu.tab.simple.ss.nozero.bray.upgma.boot<-recluster.boot(otu.tab.simple.ss.nozero.bray.upgma$cons, otu.tab.simple.ss.nozero, 
                                                         tr=100, p=0.5, method="average", boot=1000, level=1)
recluster.plot(otu.tab.simple.ss.nozero.bray.upgma$cons, otu.tab.simple.ss.nozero.bray.upgma.boot) # We add bootstrap values to the branching pattern

#clr transformed dataset
# We generate 100 trees by resampling and then, we build the consensus
otu.tab.simple.gbm.clr.euclidean.upgma<-recluster.cons(otu.tab.simple.gbm.clr.euclidean, tr=100, p=0.5, method="average") 
plot(otu.tab.simple.gbm.clr.euclidean.upgma$cons) # plot consensus tree
# We'll calculate bootstrap support values (0: bad - 100: perfect)
otu.tab.simple.gbm.clr.euclidean.upgma.boot<-recluster.boot(otu.tab.simple.gbm.clr.euclidean.upgma$cons, t(otu.tab.simple.gbm.clr), 
                                                         tr=100, p=0.5, method="average", boot=100, level=1)
recluster.plot(otu.tab.simple.gbm.clr.euclidean.upgma$cons, otu.tab.simple.gbm.clr.euclidean.upgma.boot) # We add bootstrap values to the branching pattern

#Let's compare both dendrograms using tanglegrams

install.packages("dendextend")
library(dendextend)

dendlist(as.dendrogram(otu.tab.simple.ss.nozero.bray.upgma$cons), as.dendrogram(otu.tab.simple.gbm.clr.euclidean.upgma$cons)) %>%
  untangle(method = "step1side") %>% # Find the best alignment layout
  tanglegram(cex_main=0.7, cex_sub=1, lwd=2.0, main_left="rarefied", main_right="clr transformed",cex_main_left=2, lab.cex=1.5, edge.lwd=2)

#Analyses using environmental variation
# We aim at investigating the environmental variation that may explain community variance.
# Read environmental table
bbmo.metadata.course<-read_tsv("https://raw.githubusercontent.com/krabberod/BIO9905MERG1_V21/main/community.ecology/bbmo.metadata.course.tsv", col_names = T)
bbmo.metadata.course<-as.data.frame(bbmo.metadata.course)
rownames(bbmo.metadata.course)<-bbmo.metadata.course[,1]
bbmo.metadata.course<-bbmo.metadata.course[,-1]

#                            BL040126 BL040419 BL040719 BL041019 BL050120 BL050413 BL050705 BL051004
# ENV_Temp                         14     12.6       24     19.2       13       13       24     21.5
# ENV_SECCHI                       14        6       24       12       19       18       22       17
# ENV_SAL_CTD                    37.9     35.9     36.9     37.5       37     37.7    37.35     35.1
# ENV_CHL_total                   1.1      1.4      0.4      0.3      0.5        2      0.1      0.6
# ENV_PO4                         0.2      0.2      0.1      0.1      0.2      0.3      0.2      0.2
# ENV_NH4                         0.3      1.5        1      0.5      1.1      2.1      1.4      1.5
# ENV_NO2                         0.3      0.4      0.2      0.1      0.2      0.4      0.1      0.1
# ENV_NO3                         1.5      2.5      0.1      0.4      1.1      3.3      0.2      2.4
# ENV_SI                          1.8      6.1      1.4      1.4      2.6      3.4      1.8      1.6
# ENV_BACTERIA                 854356  1046779  1654834  1083724   582655   788163  1127596   885144
# ENV_SYNECHOS                   5927     1411    38741  30915.5     8253     4169    24823    33866
# ENV_Micromonas                 9258     1424      203      730     4414     1543      505      573
# ENV_PNF_tot                   11451     2266     1228     2811     5853     2506     1699     2052
# ENV_HNF_tot                     329     1793     1357      822      420      669     1528      837
# ENV_Day_length_Hours_light      9.8    13.51    14.81    10.94     9.61     13.2    15.12    11.67
# Month                        01_jan   04_apr   07_jul   10_oct   01_jan   04_apr   07_jul   10_oct
# Season                          win      spr      sum      aut      win      spr      sum      aut
# Season_corr                     win      spr      sum      aut      win      spr      sum      aut
# Year                           2004     2004     2004     2004     2005     2005     2005     2005

#Double check samples are correct
identical(colnames(otu.tab.simple.gbm.clr),colnames(bbmo.metadata.course))
# [1] TRUE #Both tables have the same names
#We transform variables 1:15 using z-scores to have comparable ranges of variation
bbmo.metadata.course.15vars<-bbmo.metadata.course[1:15,] #We select continuous variables
bbmo.metadata.course.15vars[]<- lapply(bbmo.metadata.course.15vars, as.character) #We transform the datatype to characters
bbmo.metadata.course.15vars[]<- lapply(bbmo.metadata.course.15vars, as.numeric) #We transform to numeric
#lapply : applies a function to a list object
bbmo.metadata.course.15vars.zscores<-scale(t(bbmo.metadata.course.15vars), center = T, scale = T)
bbmo.metadata.course.15vars.zscores[,1:5]

#            ENV_Temp  ENV_SECCHI ENV_SAL_CTD ENV_CHL_total    ENV_PO4
# BL040126 -0.7223777 -0.43425521  1.02225526     0.4644927  0.1950474
# BL040419 -0.9985084 -1.82387188 -1.06132234     0.9289853  0.1950474
# BL040719  1.2499845  1.30276563 -0.01953354    -0.6193235 -1.3653316
# BL041019  0.3032507 -0.78165938  0.60553974    -0.7741544 -1.3653316
# BL050120 -0.9196139  0.43425521  0.08464534    -0.4644927  0.1950474
# BL050413 -0.9196139  0.26055313  0.81389750     1.8579706  1.7554264
# BL050705  1.2499845  0.95536146  0.44927142    -1.0838162  0.1950474
# BL051004  0.7568940  0.08685104 -1.89475338    -0.3096618  0.1950474

#Let's check the correlation in environmental variables
install.packages("corrplot") # makes nice correlation plots
install.packages("RcmdrMisc") # diverse tools
library("corrplot")
library("RcmdrMisc")

#We calculate correlations and p-values
env.corr.signif.adjust<-rcorr.adjust(as.matrix(bbmo.metadata.course.15vars.zscores)) # The p-values are corrected for multiple inference using Holm's method (see p.adjust).
#More info in: https://en.wikipedia.org/wiki/Multiple_comparisons_problem
#Holm corrected values for multiple comparisons
env.corr.signif.r<-env.corr.signif.adjust$R$r
env.corr.signif.p<-env.corr.signif.adjust$P
# We edit the objetc to replace any "<" by "0" using the function "gsub"
env.corr.signif.p<-gsub("<","0", env.corr.signif.p)
# We modify the object to be numeric datatype. #NB: the transformation is done so the matrix of p values can be read as numeric!
env.corr.signif.p <- apply(env.corr.signif.p, 2 ,as.numeric)
# We plot the correlation plot
corrplot(env.corr.signif.r , type="upper", order="hclust", p.mat = env.corr.signif.p, sig.level = 0.05,
         insig = "pch", hclust.method = c("average"), tl.cex= 0.8, tl.col="black", diag=F)

#Fitting environmental variables to ordinations
#   envfit will fit the environmental variables to the NMDS ordination as vectors

#Rarefied table
otu.tab.simple.ss.nozero.bray.nmds.envfit<-envfit(otu.tab.simple.ss.nozero.bray.nmds,bbmo.metadata.course.15vars.zscores, permu=999)
# ***VECTORS
#                              NMDS1    NMDS2     r2   Pr(>r)  
# ENV_Temp                    0.95531  0.29562 0.0939  0.784  
# ENV_SECCHI                  0.76214  0.64742 0.0088  0.983  
# ENV_SAL_CTD                -0.95959  0.28139 0.0737  0.815  
# ENV_CHL_total              -0.17996  0.98367 0.2149  0.557  
# ENV_PO4                     0.86025  0.50987 0.1955  0.599  
# ENV_NH4                     0.60519  0.79608 0.6083  0.086 .
# ENV_NO2                    -0.41451  0.91004 0.0993  0.763  
# ENV_NO3                     0.23214  0.97268 0.1340  0.718  
# ENV_SI                      0.73053  0.68288 0.0655  0.838  
# ENV_BACTERIA               -0.11815  0.99300 0.0809  0.825  
# ENV_SYNECHOS               -0.17166  0.98516 0.0147  0.976  
# ENV_Micromonas             -0.41872 -0.90811 0.5020  0.182  
# ENV_PNF_tot                -0.44058 -0.89771 0.5355  0.164  
# ENV_HNF_tot                 0.88769  0.46043 0.3467  0.348  
# ENV_Day_length_Hours_light  0.66598  0.74597 0.5677  0.115  
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  Permutation: free
#  Number of permutations: 999
# The two last columns indicate the squared correlation coefficient and the associated p-value
# We plot the vectors of the significant correlations
plot(otu.tab.simple.ss.nozero.bray.nmds, type="t", display="sites")
plot(otu.tab.simple.ss.nozero.bray.nmds.envfit) # We plot all vectors
plot(otu.tab.simple.ss.nozero.bray.nmds.envfit, p.max=0.1) # We plot all vectors that are significant, what happens?

#clr table
otu.tab.simple.gbm.clr.euclidean.nmds.envfit<-envfit(otu.tab.simple.gbm.clr.euclidean.nmds,bbmo.metadata.course.15vars.zscores, permu=999)
# ***VECTORS
#                               NMDS1    NMDS2    r2  Pr(>r)
# ENV_Temp                    0.64026  0.76816 0.2842  0.459
# ENV_SECCHI                  0.09551  0.99543 0.1629  0.646
# ENV_SAL_CTD                -0.09439  0.99553 0.1041  0.749
# ENV_CHL_total               0.10611 -0.99435 0.2019  0.572
# ENV_PO4                     0.19861 -0.98008 0.1624  0.662
# ENV_NH4                     0.89853 -0.43890 0.2838  0.451
# ENV_NO2                    -0.07233 -0.99738 0.2818  0.454
# ENV_NO3                     0.13336 -0.99107 0.2223  0.532
# ENV_SI                      0.02947 -0.99957 0.2719  0.460
# ENV_BACTERIA                0.96856  0.24880 0.1424  0.687
# ENV_SYNECHOS                0.30439  0.95255 0.2973  0.459
# ENV_Micromonas             -0.67561 -0.73726 0.5987  0.104
# ENV_PNF_tot                -0.73873 -0.67400 0.5464  0.148
# ENV_HNF_tot                 0.68687 -0.72678 0.3742  0.315
# ENV_Day_length_Hours_light  0.94169 -0.33648 0.5634  0.149
# Permutation: free
# Number of permutations: 999

plot(otu.tab.simple.gbm.clr.euclidean.nmds, type="t", display="sites")
plot(otu.tab.simple.gbm.clr.euclidean.nmds.envfit) # We plot all vectors
plot(otu.tab.simple.gbm.clr.euclidean.nmds.envfit, p.max=0.1) # We plot all vectors that are significant, what happens?

#Constrained Ordination
# dbRDA (distance-based redundancy analyses)
# Selection of the most important variables for dbRDA

#Rarefaction table
mod0.rarefaction<-capscale(otu.tab.simple.ss.nozero.bray~1, as.data.frame(bbmo.metadata.course.15vars.zscores)) # model containing only species matrix and intercept
mod1.rarefaction<-capscale(otu.tab.simple.ss.nozero.bray~ ., as.data.frame(bbmo.metadata.course.15vars.zscores)) # # model including all variables from env matrix (the dot after tilde (~) means ALL!)
ordistep(mod0.rarefaction, scope = formula(mod1.rarefaction), perm.max = 1000, direction="forward")

# Start: otu.tab.simple.ss.nozero.bray ~ 1 
#                              Df    AIC      F Pr(>F)  
# + ENV_PNF_tot                 1 9.2535 1.3702  0.050 *
# + ENV_Day_length_Hours_light  1 9.1702 1.4474  0.055 .
# + ENV_Micromonas              1 9.2311 1.3909  0.055 .
# + ENV_BACTERIA                1 9.4129 1.2248  0.110  
# + ENV_Temp                    1 9.3168 1.3121  0.170  
# + ENV_PO4                     1 9.4892 1.1562  0.195  
# + ENV_HNF_tot                 1 9.4548 1.1870  0.240  
# + ENV_SYNECHOS                1 9.4613 1.1812  0.255  
# + ENV_NH4                     1 9.5305 1.1193  0.285  
# + ENV_NO3                     1 9.5767 1.0784  0.350  
# + ENV_NO2                     1 9.6558 1.0087  0.380  
# + ENV_CHL_total               1 9.5684 1.0857  0.385  
# + ENV_SAL_CTD                 1 9.6782 0.9891  0.485  
# + ENV_SI                      1 9.7718 0.9078  0.590  
# + ENV_SECCHI                  1 9.8076 0.8770  0.675  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Step: otu.tab.simple.ss.nozero.bray ~ ENV_PNF_tot 
#                              Df     AIC      F Pr(>F)
# + ENV_SAL_CTD                 1  9.5007 1.2248  0.225
# + ENV_PO4                     1  9.4897 1.2333  0.235
# + ENV_CHL_total               1  9.5584 1.1800  0.285
# + ENV_NO3                     1  9.6004 1.1477  0.285
# + ENV_NH4                     1  9.6031 1.1456  0.330
# + ENV_NO2                     1  9.7120 1.0625  0.370
# + ENV_Temp                    1  9.7212 1.0555  0.380
# + ENV_SYNECHOS                1  9.7690 1.0195  0.465
# + ENV_SI                      1  9.7931 1.0013  0.600
# + ENV_BACTERIA                1  9.8945 0.9258  0.620
# + ENV_HNF_tot                 1  9.9316 0.8983  0.645
# + ENV_SECCHI                  1  9.9584 0.8786  0.675
# + ENV_Day_length_Hours_light  1 10.0502 0.8115  0.720
# + ENV_Micromonas              1 10.0363 0.8216  0.745

# Call: capscale(formula = otu.tab.simple.ss.nozero.bray ~ ENV_PNF_tot, data = as.data.frame(bbmo.metadata.course.15vars.zscores))
# NB: the variables in this model are the ones that were selected.

# Inertia Proportion Rank
# Total          2.7072     1.0000     
# Constrained    0.5033     0.1859    1
# Unconstrained  2.2039     0.8141    6
# Inertia is squared Bray distance 

# Eigenvalues for constrained axes:
#   CAP1 
# 0.5033 

# Eigenvalues for unconstrained axes:
#  MDS1   MDS2   MDS3   MDS4   MDS5   MDS6 
# 0.5958 0.4353 0.3912 0.2791 0.2778 0.2246 

#clr table
mod0.clr<-capscale(otu.tab.simple.gbm.clr.euclidean~1, as.data.frame(bbmo.metadata.course.15vars.zscores)) # model containing only species matrix and intercept
mod1.clr<-capscale(otu.tab.simple.gbm.clr.euclidean~ ., as.data.frame(bbmo.metadata.course.15vars.zscores)) # # model including all variables from env matrix (the dot after tilde (~) means ALL!)
ordistep(mod0.clr, scope = formula(mod1.clr), perm.max = 1000, direction="forward")

# Start: otu.tab.simple.gbm.clr.euclidean ~ 1 
# Df    AIC      F Pr(>F)  
# + ENV_Day_length_Hours_light  1 76.928 2.0796  0.025 *
# + ENV_Micromonas              1 77.021 1.9868  0.060 .
# + ENV_PNF_tot                 1 77.049 1.9591  0.065 .
# + ENV_NH4                     1 77.529 1.4954  0.110  
# + ENV_HNF_tot                 1 77.774 1.2687  0.190  
# + ENV_Temp                    1 78.088 0.9896  0.400  
# + ENV_SECCHI                  1 78.249 0.8497  0.500  
# + ENV_BACTERIA                1 78.185 0.9050  0.535  
# + ENV_SI                      1 78.382 0.7368  0.675  
# + ENV_SYNECHOS                1 78.334 0.7778  0.685  
# + ENV_PO4                     1 78.364 0.7525  0.695  
# + ENV_SAL_CTD                 1 78.523 0.6196  0.860  
# + ENV_NO2                     1 78.604 0.5528  0.940  
# + ENV_NO3                     1 78.680 0.4909  0.980  
# + ENV_CHL_total               1 78.756 0.4298  0.995  
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Step: otu.tab.simple.gbm.clr.euclidean ~ ENV_Day_length_Hours_light 

# Df    AIC      F Pr(>F)
# + ENV_BACTERIA    1 77.310 1.1210  0.365
# + ENV_Micromonas  1 77.406 1.0484  0.390
# + ENV_SECCHI      1 77.418 1.0391  0.470
# + ENV_PNF_tot     1 77.477 0.9945  0.475
# + ENV_HNF_tot     1 77.568 0.9267  0.565
# + ENV_PO4         1 77.596 0.9060  0.630
# + ENV_NH4         1 77.634 0.8781  0.640
# + ENV_SI          1 77.620 0.8880  0.650
# + ENV_NO3         1 77.918 0.6733  0.735
# + ENV_SAL_CTD     1 77.857 0.7165  0.760
# + ENV_SYNECHOS    1 77.855 0.7183  0.795
# + ENV_NO2         1 77.943 0.6556  0.825
# + ENV_Temp        1 78.013 0.6063  0.870
# + ENV_CHL_total   1 78.151 0.5100  0.920

# Call: capscale(formula = otu.tab.simple.gbm.clr.euclidean ~ ENV_Day_length_Hours_light, data = as.data.frame(bbmo.metadata.course.15vars.zscores))

# Inertia Proportion Rank
# Total         1.400e+04  1.000e+00     
# Constrained   3.605e+03  2.574e-01    1
# Unconstrained 1.040e+04  7.426e-01    6
# Inertia is mean squared Euclidean distance 

# Eigenvalues for constrained axes:
#  CAP1 
#  3605 

# Eigenvalues for unconstrained axes:
#  MDS1 MDS2 MDS3 MDS4 MDS5 MDS6 
#  3253 2380 1780 1498  931  558 


#Generate the ordination
# We will use two more variables that we know they are important in for this dataset

#We install ggord for nicer plots
library(devtools)
install_github('fawda123/ggord')
library(ggord)
library(ggplot2)

#rarefied table
ggord(dbrda(formula = otu.tab.simple.ss.nozero.bray ~ ENV_PNF_tot+ENV_Day_length_Hours_light+ENV_Temp, data = as.data.frame(bbmo.metadata.course.15vars.zscores)))
screeplot(dbrda(formula = otu.tab.simple.ss.nozero.bray ~ ENV_PNF_tot+ENV_Day_length_Hours_light+ENV_Temp, data = as.data.frame(bbmo.metadata.course.15vars.zscores)))
dbrda(formula = otu.tab.simple.ss.nozero.bray ~ ENV_PNF_tot+ENV_Day_length_Hours_light+ENV_Temp, data = as.data.frame(bbmo.metadata.course.15vars.zscores))

# Call: dbrda(formula = otu.tab.simple.ss.nozero.bray ~ ENV_PNF_tot + ENV_Day_length_Hours_light + ENV_Temp, data =
#               as.data.frame(bbmo.metadata.course.15vars.zscores))
#                Inertia Proportion Rank
# Total          2.7072     1.0000     
# Constrained    1.1945     0.4412    3   # Community variation constrained by the used variables
# Unconstrained  1.5127     0.5588    4
# Inertia is squared Bray distance  # Inertia = variance in species abundances

# Eigenvalues for constrained axes:
#   dbRDA1 dbRDA2 dbRDA3 
#   0.5701 0.3699 0.2545 

# Eigenvalues for unconstrained axes:
#     MDS1   MDS2   MDS3   MDS4 
#  0.5818 0.3960 0.2997 0.2353 

#clr table
ggord(dbrda(formula = otu.tab.simple.gbm.clr.euclidean ~ ENV_PNF_tot+ENV_Day_length_Hours_light+ENV_Temp, data = as.data.frame(bbmo.metadata.course.15vars.zscores)))
screeplot(dbrda(formula = otu.tab.simple.gbm.clr.euclidean ~ ENV_PNF_tot+ENV_Day_length_Hours_light+ENV_Temp, data = as.data.frame(bbmo.metadata.course.15vars.zscores)))
dbrda(formula = otu.tab.simple.gbm.clr.euclidean ~ ENV_PNF_tot+ENV_Day_length_Hours_light+ENV_Temp, data = as.data.frame(bbmo.metadata.course.15vars.zscores))

# Call: dbrda(formula = otu.tab.simple.gbm.clr.euclidean ~ ENV_PNF_tot + ENV_Day_length_Hours_light + ENV_Temp, data =
#              as.data.frame(bbmo.metadata.course.15vars.zscores))

#                  Inertia Proportion Rank
# Total         14004.823      1.000     
# Constrained    6399.625      0.457    3   # Community variation constrained by the used variables
# Unconstrained  7605.199      0.543    4
# Inertia is mean squared Euclidean distance 

# Eigenvalues for constrained axes:
#  dbRDA1 dbRDA2 dbRDA3 
#   3848   1491   1061 

# Eigenvalues for unconstrained axes:
#     MDS1   MDS2   MDS3   MDS4 
#   3033.4 2224.1 1504.9  842.9 

#Constrained Correspondence Analyses (CCA)
# Selection of the most important variables
# rarefaction table
mod0.cca.rarefaction<-cca(otu.tab.simple.ss.nozero~1, as.data.frame(bbmo.metadata.course.15vars.zscores)) # model containing only species matrix and intercept
mod1.cca.rarefaction<-cca(otu.tab.simple.ss.nozero~ ., as.data.frame(bbmo.metadata.course.15vars.zscores)) # # model including all variables from env matrix (the dot after tilde (~) means ALL!)
ordistep(mod0.cca.rarefaction, scope = formula(mod1.cca.rarefaction), perm.max = 1000, direction="forward")

# Call: cca(formula = otu.tab.simple.ss.nozero ~ ENV_Day_length_Hours_light, data = as.data.frame(bbmo.metadata.course.15vars.zscores)) # Best model

# Inertia Proportion Rank
# Total          4.1915     1.0000     
# Constrained    0.6904     0.1647    1   # Community variation constrained by the used variables
# Unconstrained  3.5012     0.8353    6
# Inertia is scaled Chi-square 

# Eigenvalues for constrained axes:
#   CCA1 
# 0.6904 

# Eigenvalues for unconstrained axes:
# CA1    CA2    CA3    CA4    CA5    CA6 
# 0.7987 0.6988 0.6026 0.5395 0.5271 0.3345 

ggord(cca(formula = otu.tab.simple.ss.nozero ~ ENV_PNF_tot+ENV_Day_length_Hours_light+ENV_Temp, data = as.data.frame(bbmo.metadata.course.15vars.zscores)), obslab=T, addsize=0.6)
screeplot(cca(formula = otu.tab.simple.ss.nozero ~ ENV_PNF_tot+ENV_Day_length_Hours_light+ENV_Temp, data = as.data.frame(bbmo.metadata.course.15vars.zscores)))
cca(formula = otu.tab.simple.ss.nozero ~ ENV_PNF_tot+ENV_Day_length_Hours_light+ENV_Temp, data = as.data.frame(bbmo.metadata.course.15vars.zscores))

# Call: cca(formula = otu.tab.simple.ss.nozero ~ ENV_PNF_tot + ENV_Day_length_Hours_light + ENV_Temp, data =
#            as.data.frame(bbmo.metadata.course.15vars.zscores))

# Inertia Proportion Rank
# Total          4.1915     1.0000     
# Constrained    1.7827     0.4253    3   # Community variation constrained by the used variables
# Unconstrained  2.4089     0.5747    4
# Inertia is scaled Chi-square 

# Eigenvalues for constrained axes:
#   CCA1   CCA2   CCA3 
# 0.7137 0.5997 0.4692 

# Eigenvalues for unconstrained axes:
#  CA1    CA2    CA3    CA4 
# 0.7667 0.6562 0.5585 0.4275 

#PERMANOVA

#rarefied table
permanova.rarefaction<-adonis(otu.tab.simple.ss.nozero~ENV_PNF_tot+ENV_Day_length_Hours_light+ENV_Temp, data = as.data.frame(bbmo.metadata.course.15vars.zscores), 
                              method="bray", permutations=9999)
# Call:
# adonis(formula = otu.tab.simple.ss.nozero ~ ENV_PNF_tot + ENV_Day_length_Hours_light + ENV_Temp, data = as.data.frame(bbmo.metadata.course.15vars.zscores), permutations = 9999, method = "bray") 

# Permutation: free
# Number of permutations: 9999
# Terms added sequentially (first to last)

#                            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# ENV_PNF_tot                 1   0.50330 0.50330 1.33088 0.18591 0.0844 . # PNF is close to significance, explaining ca. 18% of the variance
# ENV_Day_length_Hours_light  1   0.30774 0.30774 0.81378 0.11368 0.7359  
# ENV_Temp                    1   0.38345 0.38345 1.01397 0.14164 0.4271  
# Residuals                   4   1.51267 0.37817         0.55877         
# Total                       7   2.70716                 1.00000         
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#
#clr-transformed table
permanova.clr<-adonis(t(otu.tab.simple.gbm.clr)~ENV_PNF_tot+ENV_Day_length_Hours_light+ENV_Temp, data = as.data.frame(bbmo.metadata.course.15vars.zscores), 
                      method="euclidean", permutations=9999)
# Call:
#  adonis(formula = t(otu.tab.simple.gbm.clr) ~ ENV_PNF_tot + ENV_Day_length_Hours_light +      ENV_Temp, data = as.data.frame(bbmo.metadata.course.15vars.zscores),      permutations = 9999, method = "euclidean") 

# Permutation: free
# Number of permutations: 9999
# Terms added sequentially (first to last)

#                             Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)
# ENV_PNF_tot                 1     24130 24130.3 1.81307 0.24614 0.1085
# ENV_Day_length_Hours_light  1     13180 13180.0 0.99030 0.13444 0.3917
# ENV_Temp                    1      7487  7487.1 0.56255 0.07637 0.9368
# Residuals                   4     53236 13309.1         0.54304       
# Total                       7     98034                 1.00000 


#Tutorial
#Follow all the steps indicated above with the provided data.

#Files
# All content in https://github.com/krabberod/BIO9905MERG1_V21/tree/main/community.ecology

# OTU table raw
otu.tab<-read_tsv("https://raw.githubusercontent.com/krabberod/BIO9905MERG1_V21/main/Dada2_Pipeline/dada2_results/OTU_table.tsv")
dim(otu.tab) # 2107   26
#Let's reorder the table
otu.tab<-otu.tab[,c(17,19:26,1:16,18)]
#We assign to rownames the OTU names
otu.tab <- column_to_rownames(otu.tab, var = "OTUNumber") # %>% as_tibble()
dim(otu.tab) # 2107   25 <- Dimensions of the table
otu.tab.simple<-otu.tab[,1:8] # We'll need this table for community ecology analyses
#We transpose the table, as this is how Vegan likes it
otu.tab.simple<-t(otu.tab.simple)
otu.tab.simple # => ready to use <= 


# OTU table rarefied 
otu.tab.simple.ss.nozero<-read_tsv("https://raw.githubusercontent.com/krabberod/BIO9905MERG1_V21/main/community.ecology/otu.tab.simple.ss.nozero.tsv", col_names = T)
otu.tab.simple.ss.nozero<-as.data.frame(otu.tab.simple.ss.nozero) # transform to dataframe
rownames(otu.tab.simple.ss.nozero)<-otu.tab.simple.ss.nozero[,1] # fix row names
otu.tab.simple.ss.nozero<-otu.tab.simple.ss.nozero[,-1] # fix row names
otu.tab.simple.ss.nozero # => ready to use <=

# Bray Curtis distance matrix
otu.tab.simple.ss.nozero.bray<-read_tsv("https://raw.githubusercontent.com/krabberod/BIO9905MERG1_V21/main/community.ecology/otu.tab.simple.ss.nozero.mat.tsv", col_names = T)
otu.tab.simple.ss.nozero.bray<-as.dist(otu.tab.simple.ss.nozero.bray) # Transform to a distance object
otu.tab.simple.ss.nozero.bray # <- ready to use

# OTU table centered log-ratio transformed
otu.tab.simple.gbm.clr<-read_tsv("https://raw.githubusercontent.com/krabberod/BIO9905MERG1_V21/main/community.ecology/otu.tab.simple.gbm.clr.tsv", col_names = T)
otu.tab.simple.gbm.clr<-as.data.frame(otu.tab.simple.gbm.clr) # transform to dataframe
rownames(otu.tab.simple.gbm.clr)<-otu.tab.simple.gbm.clr[,1]# fix row names
otu.tab.simple.gbm.clr<-otu.tab.simple.gbm.clr[,-1] # fix row names
otu.tab.simple.gbm.clr # => ready to use <=

# Aitchison distance matrix
otu.tab.simple.gbm.clr.euclidean<-read_tsv("https://raw.githubusercontent.com/krabberod/BIO9905MERG1_V21/main/community.ecology/otu.tab.simple.gbm.clr.mat.tsv", col_names = T)
otu.tab.simple.gbm.clr.euclidean<-as.dist(otu.tab.simple.gbm.clr.euclidean) # Transform to a distance object
otu.tab.simple.gbm.clr.euclidean # => ready to use <=

# Environmental table
bbmo.metadata.course<-read_tsv("https://raw.githubusercontent.com/krabberod/BIO9905MERG1_V21/main/community.ecology/bbmo.metadata.course.tsv", col_names = T)
bbmo.metadata.course<-as.data.frame(bbmo.metadata.course)
rownames(bbmo.metadata.course)<-bbmo.metadata.course[,1]
bbmo.metadata.course<-bbmo.metadata.course[,-1]
bbmo.metadata.course # <- Ready to use, table to z-score transformed
#We transform variables 1:15 using z-scores to have comparable ranges of variation
bbmo.metadata.course.15vars<-bbmo.metadata.course[1:15,] #We select continuous variables
bbmo.metadata.course.15vars[]<- lapply(bbmo.metadata.course.15vars, as.character) #We transform the datatype to characters
bbmo.metadata.course.15vars[]<- lapply(bbmo.metadata.course.15vars, as.numeric) #We transform to numeric
bbmo.metadata.course.15vars.zscores<-scale(t(bbmo.metadata.course.15vars), center = T, scale = T) #zscore transform
bbmo.metadata.course.15vars.zscores # => ready to use <=

# END





















