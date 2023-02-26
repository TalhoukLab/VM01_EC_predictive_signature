## Required packages
library(metagenomeSeq)
library(vegan)
library(biomformat)
library(phyloseq)
library(tidyr)
library(microbiome)
library(ggpubr)
library(upstartr)
library(dplyr)
library(OTUtable)
library(picante)
library(dplyr)
library(reshape2)
library(mia)
library(miaViz)
library(pheatmap)
library(ggtree)
library(gsubfn)
library(dplyr)

## Making phyloseq objects 

cohorts <- c("Angel", "Antonio", "Gressel", "Tsementzi", "Walsh")
for (cohort in cohorts) {
  print(cohort)
  raw_table <- read_biom(file.path(paste0('~/Desktop/R&R/vaginalMicrobiome/Tsementzi_pipeline/results/', cohort, '_cohort/filtered_sortedUnique_otu.biom')))
  raw_otus <- as.data.frame(as.matrix(biom_data(raw_table)))
  
  # Remove empty sample 
  to_keep <- which((colSums(raw_otus))>1)
  raw_otus <- raw_otus[, to_keep]

  # CSS normalization
  metaSeqObject <- newMRexperiment(raw_otus) 
  CSS <- cumNorm(metaSeqObject, p=cumNormStat(metaSeqObject))
  outs_CSS = data.frame(MRcounts(CSS, norm=TRUE, log=T))

  # read in tax data 
  tax_table <- read.delim(file.path(paste0('~/Desktop/R&R/vaginalMicrobiome/Tsementzi_pipeline/results/', cohort, '_cohort/seqs_sorted_rep_set_tax_assignments.txt')), header = FALSE, sep = "\t")
  tax_table <- separate(data = tax_table, col = V2 , into = c("kingdom", "phylum", "class", "order", "family", "genus"), sep = "\\;")
  rownames(tax_table) <- tax_table$V1
  tax_table[is.na(tax_table)] <- ""
  tax_table = subset(tax_table, select = -c(V3, V1))
  tax_table$kingdom<-gsub("k__","",as.character(tax_table$kingdom))
  tax_table$phylum <- gsub("p__", "", gsub("\\[", "",  gsub("\\]", "" ,as.character(tax_table$phylum))))
  tax_table$class <- gsub("c__", "", gsub("\\[", "",  gsub("\\]", "" ,as.character(tax_table$class))))
  tax_table$order <- gsub("o__", "", gsub("\\[", "",  gsub("\\]", "" ,as.character(tax_table$order))))
  tax_table$family <- gsub("f__", "", gsub("\\[", "",  gsub("\\]", "" ,as.character(tax_table$family))))
  tax_table$genus <- gsub("g__", "", gsub("\\[", "",  gsub("\\]", "" ,as.character(tax_table$genus))))
  tax_table <- mutate_all(tax_table, .funs=tolower)
  
  phylo_tree <- read_tree(file.path(paste0('~/Desktop/R&R/vaginalMicrobiome/Tsementzi_pipeline/results/', cohort, '_cohort/rep_set.tre')))
  
  feature_table <- read.delim(file.path(paste0('~/Desktop/R&R/vaginalMicrobiome/00-helperfiles/', cohort, 'FT.csv')), header = TRUE, sep = ",")
  rownames(feature_table) <- feature_table$sraID
  feature_table$histology[feature_table$histology == "ACH"] <- "EC"
  feature_table$histology <- as.factor(feature_table$histology)
  feature_table$histology <- relevel(feature_table$histology,"Benign")
  
  otus_table <- otu_table(outs_CSS, taxa_are_rows = TRUE)
  tax_table_phy = tax_table(as.matrix(tax_table))
  tree_phy <- phy_tree(phylo_tree)
  samples = sample_data(feature_table)
  
  # Making phyloseq object without tree
  #phyloseq_pre <- phyloseq(otus_table, tax_table_phy, samples)
  #phyloseq_pre_raw <- phyloseq(otu_table(raw_otus, taxa_are_rows = T), tax_table_phy, samples)
  # Filter phyloseq object without tree
  #phyloseq_obj <- subset_taxa(phyloseq_pre, kingdom!="Unclassified")
  #phyloseq_obj_raw <- subset_taxa(phyloseq_pre_raw, kingdom!="Unclassified")
  #assign(paste0(cohort, "_Tsementziphyloseq"), phyloseq_pre,.GlobalEnv)
  #assign(paste0(cohort, "_Tsementziphyloseq_raw"),phyloseq_pre_raw,.GlobalEnv)
  
  # Making phyloseq object with tree
  phyloseq_pre_tree <- phyloseq(otus_table, tax_table_phy, samples, tree_phy)
  phyloseq_pre_tree_raw <- phyloseq(otu_table(raw_otus, taxa_are_rows = T), tax_table_phy, samples, tree_phy)
  
  assign(paste0(cohort, "_Tsementziphyloseq_tree"),phyloseq_pre_tree,.GlobalEnv)
  assign(paste0(cohort, "_Tsementziphyloseq_tree_raw"),phyloseq_pre_tree_raw,.GlobalEnv)
}
