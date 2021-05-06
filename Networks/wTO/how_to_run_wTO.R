library(tidyverse)
library(vegan)
library(wTO)
library(magrittr)
library(readr)
library(ptm)
source("filterCoOccuringwTO.R")

###  read OTU table that will be used for network construction ### 

table <- read.csv("./../Dada2_Pipeline/dada2_results/OTU_table.tsv",sep = "\t") 
#row.names(table)<-table$OTUNumber # Set the name of the OTUs
row.names(table)<-paste0(table$OTUNumber,"_",table$Species)


otu.table<-table[,c(19:26)]

#  To filter only the prevalent OTUs this is what we do: 
otu.table.pa<- decostand(otu.table, method = "pa")
otu.table.pa.pr<- otu.table.pa[rowSums(otu.table.pa)>=2,] 
otus.prev<- otu.table[rownames(otu.table) %in% rownames(otu.table.pa.pr),]

# wTO requires a data frame to work on
otus.prev<-as.data.frame(otus.prev)


###  To determine the autocorrelation of each OTU ###
# This is important when working with temporal series as it determines the parameter "lag" 
# that will be used in the wTO.complete function. 
# If most of the OTUs present an autocorrelation of 2, then lag=2.

par(mfrow = c(3,3))
for ( i in 1:nrow(otus.prev)){
  acf(t(otus.prev[i,]))
}

###  Network calculation  ### 
ptm_network <- proc.time()

otus.prev.complete = wTO.Complete(k = 1, n = 250, Data = otus.prev, Overlap = row.names(otus.prev), method = 's' , 
                              method_resampling = 'BlockBootstrap', lag = 2, pvalmethod = 'BH', savecor = TRUE, 
                              expected.diff = 0.2, plot = T)
ptm_network <- proc.time() - ptm_network
ptm_network


saveRDS (otus.prev.complete, "otus.prev.complete.prev.2.rds")
# These are the paramenters that we used, but they can be changed: 
# k=threads to be used
# n=amount of replicates
# method=s:Spearman, p=Pearson
# method_resampling= BlockBootstrap. This is the method to be used when working with temporal series. It can be used: Bootstrap or Reshuffle 
# lag=2; time dependency. When using blockBootstraping method
# pvalmethod= BH p-value correction method 
# savecor = TRUE/FALSE: if need to save the correlation
# expected.diff=expected difference between the resampled values and the ωi,j
# plot = TRUE/FALSE. If TRUE, the diagnosis plots are generated and let you know how well the analysis has been done

### Filtering associations ###

otus.prev.compl.filtered <- filterCooccuring(otus.prev.complete$wTO, otus.prev) 
otus.prev.compl.filtered.0.05 <- subset(otus.prev.compl.filtered, otus.prev.compl.filtered$pval_sig < 0.05)    
write.table(otus.prev.compl.filtered.0.001, file = "wTO_network.txt", append = FALSE, quote = F, sep = "\t", dec = ".", row.names = F, col.names = TRUE)


### To export the resulting network without filtering associations ###
wTO.export(otus.prev.complete, "wTO_network.txt", sign = TRUE, padj = 0.05) # This is similar to: 

otus.wTO.filtered <- subset(otus.prev.complete$wTO, otus.prev.complete$wTO$Padj_sig<0.05) 
write.table(otus.wTO.filtered, file = "wTO_network.txt", append = FALSE, quote = F, sep = "\t", dec = ".", row.names = F, col.names = TRUE)

