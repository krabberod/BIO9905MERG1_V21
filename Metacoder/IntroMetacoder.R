###########METACODER SCRIPT FOR BIOMERG9905

# Install metacoder

install.packages("metacoder")

#you can also install the development version for the newest functions
install.packages("devtools")
devtools::install_github("grunwaldlab/metacoder")

#Read libraries
library(metacoder)
library(ggplot2)
library(taxa)
library(readr)
library(dplyr)
library(tidyverse)



# Read your files, containing OTU table and Metadata
#make sure that the taxonomy is attached as a column in your OTU table

otu_data<-read_tsv("rarotutable_small.txt")
print(otu_data)

sample_data<- read_tsv("sampledata.txt",
                       col_types = "cccccccccccccccccccc")
                      #the c's have to add up to the number of columns 
                      #you have in your sample data
print(sample_data)

#check that the taxonomy is attached at the end
head(otu_data$taxonomy, 10)

                 

#Parsing taxa (this sorts your taxa using the package taxa into a 
#taxmap object)
obj <- parse_tax_data(otu_data,
                      class_cols = "taxonomy",
                      class_sep = ";", # What each taxon is separated by
                      class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                      class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))

#View your taxmap object:

print(obj)

#Check if it returns ok taxon names (e.g. "Fungi", "Basidiomycota")
head(taxa:::taxon_names(obj)) 

#rename the tax_data
names(obj$data) <- "read_abundance"
#to view the entire object
print(obj)


#to make heat trees we need to calculate taxon abundance
obj$data$tax_abund <- calc_taxon_abund(obj, "read_abundance",
                                       cols = sample_data$SampleID)

#For calculating total read abundance per taxa - this gives you a list of the
#sum of all taxa across all samples
obj$data$tax_abund$total <- rowSums(obj$data$tax_abund[, -1]) 
# -1 = taxon_id column

#The dataset may contain uwanted names/taxa, these can be filtered out
#For removing names=unidentified

taxa_to_remove <- c("unidentified")
obj <- taxa::filter_taxa(obj, taxon_names %in% taxa_to_remove, invert = TRUE, subtaxa = FALSE, supertaxa = FALSE)

#Intitial tree, no layout specifications:
heat_tree(obj, node_label = taxon_names, node_size = total, 
          node_color = total)


#Heat tree with more appropriate settings:
#Here we have set taxon_ranks="o", which is order level, 
#we clean up the names, and we use 
#the read abundance (total) to make the edges and nodes
obj %>%
  taxa::filter_taxa(taxon_ranks == "o", supertaxa = TRUE, reassign_obs = FALSE) %>%
  mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = "")) %>%
  taxa::filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$"), reassign_obs = FALSE) %>%
  heat_tree(node_label = taxon_names,
            node_size = total,
            node_color = total,
            #node_color_range = c("gray",  "#85CEBF","#018571", "#A78959", "#A6611A"),
            #node_color_trans= "ln area",
            layout = "da", initial_layout = "re")

#if you want to change layout, check 
#?heat_tree()


######## Plotting heat trees for different sample types:

#In this dataset we have samples from roots and soil.
#Let's plot these individually:

view(sample_data) 
#Sample types are Bistorta, Potentilla and Soil

#Calculate read abundances for the different sample types 
#Here you have to specify groups= 

obj$data$type_abund <- calc_taxon_abund(obj, "read_abundance",
                                        cols = sample_data$SampleID,
                                        groups = sample_data$Sample)

#Make a heat tree with only taxa present in Bistorta vivipara roots
# filter_taxa is used to sort out only taxa present in Bistorta with at least 1 read:

obj %>%
  taxa::filter_taxa(Bistorta>0) %>% #if you wish to sort out very rare OTUs (i.e. less than 50 reads), you can set (Bistorta>50)
  # taxa:: needed because of phyloseq::filter_taxa
  taxa::filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names)) %>% 
  # remove "odd" taxa
  taxa::filter_taxa(taxon_ranks == "o", supertaxa = TRUE) %>% 
  # to change taxonomic resolutions, 
  #set taxon_ranks to "k","p","c","o","g","s"
  heat_tree(node_label = taxon_names,
            node_size = Bistorta, 
            node_color = Bistorta, 
            #node_color_range= quantative_palette(),
            #node_color_range=diverging_palette(),
            #node_color_range = c("gray",  "pink",  "blue"),
            #node_color_trans= "area",
            layout = "da", initial_layout = "re", 
            title = "Taxa in Bistorta vivipara")

#Plotting the taxa present with more than 50 reads in Potentilla erecta
obj %>%
  taxa::filter_taxa(Potentilla>50) %>% # if you wish to sort out very rare OTUs (i.e. less than 50 reads) # taxa:: needed because of phyloseq::filter_taxa
  taxa::filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names)) %>% # remove "odd" taxa
  taxa::filter_taxa(taxon_ranks == "o", supertaxa = TRUE) %>% # to change taxonomic resolutions, 
  #set taxon_ranks to "k","p","c","o","g","s"
  heat_tree(node_label = taxon_names,
            node_size = Potentilla, 
            node_color = Potentilla, 
            layout = "da", initial_layout = "re", 
            title = "Taxa in Potentilla erecta")


obj %>%
  taxa::filter_taxa(Soil>0) %>%  # taxa:: needed because of phyloseq::filter_taxa
  taxa::filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names)) %>% # remove "odd" taxa
  taxa::filter_taxa(taxon_ranks == "o", supertaxa = TRUE) %>% # to change taxonomic resolutions, 
  #set taxon_ranks to "k","p","c","o","g","s"
  heat_tree(node_label = taxon_names,
            node_size = Soil, 
            node_color = Soil, 
            node_color_range = c("gray",  "#85CEBF","#018571", "#A78959", "#A6611A"),
            #node_color_trans= "area",
            layout = "da", initial_layout = "re", 
            title = "Taxa in soil")


##############Compare two groups######################
#Here we compare taxa present in Plant roots and soil using Wilcoxon Rank Sum 
#test and the function compare_groups(). 
#Columns should be specified to SampleID and groups must be specified.

obj$data$diff_table <- compare_groups(obj, data = "tax_abund",
                                      cols = sample_data$SampleID,
                                      groups = sample_data$Type)

View(obj$data$diff_table)

#We then need to correct for multiple comparison, in this case we use "false discovery rate"
#but others can be specified, see:
#?p.adjust()

obj <- mutate_obs(obj, "diff_table",
                  wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))

#then set the non-significant p-values to zero to aid visualization

obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0

#Plotting the tree. Change taxon_ranks=="", for "c" (class),"f" (family), "o" (order), "g"(genus)  "s" (species)
set.seed(2)
obj %>%
  taxa::filter_taxa(taxon_ranks == "o", supertaxa = TRUE, reassign_obs = FALSE) %>%
  mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = "")) %>%
  taxa::filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$")) %>%
  heat_tree(node_label = cleaned_names,
            node_size = total, # number of OTUs
            node_color = log2_median_ratio, # difference between groups
            node_color_interval = c(-3, 3), # symmetric interval
            node_color_range = c("#D55E00", "gray", "#009E73"), # should use diverging colors
            node_size_axis_label = "read abundance",
            node_color_axis_label = "Log 2 ratio of median counts",
            layout = "da", initial_layout = "re", # good layout for large trees
            title = "Plant roots vs soil samples")




########Comparing more than two treatments###################

#Here we compare three different sample types: soil, Potentilla ercta roots and Bistorta vivipara
#roots
#(This is the same tax_abund calculation as earlier in the script, 
#somehow it gets removed and has to be redone after the previous 
#comparison.)
obj$data$tax_abund <- calc_taxon_abund(obj, "read_abundance",
                                       cols = sample_data$SampleID)
obj$data$tax_abund$total <- rowSums(obj$data$tax_abund[, -1])


#Make a diff_table using Wilcoxon Rank Sum test, specifiying groups = 
#This may give you warnings, if there are no overlapping taxa between groups. You can ignore this


obj$data$diff_table <- compare_groups(obj, data = "tax_abund",
                                      cols = sample_data$SampleID,
                                      groups = sample_data$Sample)



#Look at the diff_table
View(obj$data$diff_table)


#We then need to correct for multiple comparisons and 
#set non-significant differences to zero as above

obj <- mutate_obs(obj, "diff_table",
                  wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))

obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0

view(obj$data$diff_table)

#Plotting the tree. 
#Change taxon_ranks=="", for "c" (class),"f" (famili), "o" (order), "g"(genus)  "s" (species)
obj %>%
  taxa::filter_taxa(taxon_ranks == "o", supertaxa = TRUE, reassign_obs = FALSE) %>%
  mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = "")) %>%
  taxa::filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$"), reassign_obs = FALSE) %>%
  heat_tree_matrix(data = "diff_table",
                   node_label = cleaned_names,
                   node_size = total, # read abundance
                   node_color = log2_median_ratio, # difference between groups
                   node_color_trans = "linear",
                   node_color_interval = c(-3, 3), # symmetric interval
                   edge_color_interval = c(-3, 3), # symmetric interval
                   node_color_range = diverging_palette(), # diverging colors
                   node_size_axis_label = "read abundance",
                   node_color_axis_label = "Log 2 ratio of median counts",
                   layout = "da", initial_layout = "re",
                   key_size = 0.67,
                   seed = 2)


#Let's try on genus level:
obj %>%
  taxa::filter_taxa(taxon_ranks == "g", supertaxa = TRUE, reassign_obs = FALSE) %>%
  mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = "")) %>%
  taxa::filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$"), reassign_obs = FALSE) %>%
  heat_tree_matrix(data = "diff_table",
                   node_label = cleaned_names,
                   node_size = total, # number of OTUs
                   node_color = log2_median_ratio, # difference between groups
                   node_color_trans = "linear",
                   node_color_interval = c(-3, 3), # symmetric interval
                   edge_color_interval = c(-3, 3), # symmetric interval
                   node_color_range = diverging_palette(), # diverging colors
                   node_size_axis_label = "read abundance",
                   node_color_axis_label = "Log 2 ratio of median counts",
                   layout = "da", initial_layout = "re",
                   key_size = 0.67,
                   seed = 2)


