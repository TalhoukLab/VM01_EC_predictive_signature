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
library(pheatmap)
library(stringr)

cohorts <- c("Antonio", "Antonio_train", "Antonio_test", "Chao", "Gressel", "Tsementzi", "Walsh", "Walsh_train", "Walsh_test")
for (cohort in cohorts) {
  print(cohort)
  raw_table <- read_biom(file.path(paste0('../VM01_reproducibility_replicability/Antonio_Walsh_pipeline/results/', cohort, '/test_paired.biom')))
  raw_otus <- as.data.frame(as.matrix(biom_data(raw_table)))
  to_keep <- which((colSums(raw_otus))>1)
  otus_filtered_abun_nonEmpty <- raw_otus[, to_keep]

  # read in tax data 
  tax_table <- read.delim(file.path(paste0('../VM01_reproducibility_replicability/Antonio_Walsh_pipeline/results/', cohort, '/test_paired.probs.taxonomy')), header = FALSE, sep = "\t")
  tax_table$V2 <- gsub("\\s*\\([^\\)]+\\)","",tax_table$V2)
  
  tax_table <- separate(data = tax_table, col = V2 , into = c("kingdom", "phylum", "class", "order", "family", "genus"), sep = ";")
  rownames(tax_table) <- tax_table$V1
  tax_table[is.na(tax_table)] <- " "
  tax_table = subset(tax_table, select = -c(V1))
  tax_table$kingdom<-gsub("d:","",as.character(tax_table$kingdom))
  tax_table$phylum <- gsub("p:", "", gsub("\\[", "",  gsub("\\]", "" ,as.character(tax_table$phylum))))
  tax_table$class <- gsub("c:", "", gsub("\\[", "",  gsub("\\]", "" ,as.character(tax_table$class))))
  tax_table$order <- gsub("o:", "", gsub("\\[", "",  gsub("\\]", "" ,as.character(tax_table$order))))
  tax_table$family <- gsub("f:", "", gsub("\\[", "",  gsub("\\]", "" ,as.character(tax_table$family))))
  tax_table$genus <- gsub("g:", "", gsub("\\[", "",  gsub("\\]", "" ,as.character(tax_table$genus))))
  tax_table <- mutate_all(tax_table, .funs=tolower)
  
  phylo_tree <- read_tree(file.path(paste0('../VM01_reproducibility_replicability/Antonio_Walsh_pipeline/results/', cohort, '/test_paired.tree')))
  rooted_tree <- phangorn::midpoint(phylo_tree)
  feature_table <- read.delim(file.path(paste0('../VM01_reproducibility_replicability/00-helperfiles/', cohort, 'FT.csv')), header = TRUE, sep = ",")
  rownames(feature_table) <- feature_table$sraID
  feature_table$histology[feature_table$histology == "ACH"] <- "EC"
  feature_table$histology <- as.factor(feature_table$histology)
  feature_table$histology <- relevel(feature_table$histology,"Benign")
  otus_table <- otu_table(otus_filtered_abun_nonEmpty, taxa_are_rows = TRUE)
  tax_table_phy = tax_table(as.matrix(tax_table))
  tree_phy <- phy_tree(rooted_tree)
  samples = sample_data(feature_table)

  phyloseq_obj <- phyloseq(otus_table, tax_table_phy, samples, tree_phy)
  ps.rarefied = rarefy_even_depth(phyloseq_obj, rngseed=1, sample.size=0.9*min(sample_sums(phyloseq_obj)), replace=F)
  assign(paste0(cohort, "_Antoniophyloseq_tree"),ps.rarefied,.GlobalEnv)
  assign(paste0(cohort, "_Antoniophyloseq_tree_raw"),phyloseq_obj,.GlobalEnv)
  assign(paste0(cohort, "_Antoniophyloseq_tree_ML"),ps.rarefied,.GlobalEnv)
}
  