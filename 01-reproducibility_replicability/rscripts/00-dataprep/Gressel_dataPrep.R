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


cohorts <- c("Antonio", "Chao", "Gressel", "Gressel_train", "Gressel_test", "Tsementzi", "Walsh")
for (cohort in cohorts) {
  print(cohort)
  raw_table <- read_biom(file.path(paste0('../VM01_reproducibility_replicability/Gressel_pipeline/results/', cohort, '/feature-table.biom')))
  raw_otus <- as.data.frame(as.matrix(biom_data(raw_table)))
  to_keep <- which((colSums(raw_otus))>1)
  raw_otus <- raw_otus[, to_keep]
  tax_table <- read.delim(file.path(paste0('../VM01_reproducibility_replicability/Gressel_pipeline/results/', cohort, '/taxonomy.tsv')), header = TRUE, sep = "\t")
  tax_table <- separate(data = tax_table, col = Taxon , into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = "\\;")
  rownames(tax_table) <- tax_table$Feature.ID
  tax_table = subset(tax_table, select = -c(Feature.ID, Confidence))
  tax_table$kingdom<-gsub("k__","",as.character(tax_table$kingdom))
  tax_table$phylum <- gsub("p__", "", gsub("\\[", "",  gsub("\\]", "" ,as.character(tax_table$phylum))))
  tax_table$class <- gsub("c__", "", gsub("\\[", "",  gsub("\\]", "" ,as.character(tax_table$class))))
  tax_table$order <- gsub("o__", "", gsub("\\[", "",  gsub("\\]", "" ,as.character(tax_table$order))))
  tax_table$family <- gsub("f__", "", gsub("\\[", "",  gsub("\\]", "" ,as.character(tax_table$family))))
  tax_table$genus <- gsub("g__", "", gsub("\\[", "",  gsub("\\]", "" ,as.character(tax_table$genus))))
  tax_table[is.na(tax_table)] <- " "
  tax_table <- mutate_all(tax_table, .funs=tolower)
  
  feature_table <- read.delim(file.path(paste0('../VM01_reproducibility_replicability/00-helperfiles/', cohort, 'FT.csv')), header = TRUE, sep = ",")
  rownames(feature_table) <- feature_table$sraID
  feature_table$histology[feature_table$histology == "ACH"] <- "EC"
  feature_table$histology <- as.factor(feature_table$histology)
  feature_table$histology <- relevel(feature_table$histology,"Benign")
  otus_table <- otu_table(raw_otus, taxa_are_rows = TRUE)
  tax_table_phy = tax_table(as.matrix(tax_table))
  samples = sample_data(feature_table)

  # Making phyloseq object without tree
  phyloseq_pre <- phyloseq(otus_table, tax_table_phy, samples)
  assign(paste0(cohort, "_Gresselphyloseq_tree_raw"),phyloseq_pre,.GlobalEnv)
  assign(paste0(cohort, "_Gresselphyloseq_tree"),phyloseq_pre,.GlobalEnv)
  assign(paste0(cohort, "_Gresselphyloseq_tree_ML"),phyloseq_pre,.GlobalEnv)
  
  
}

