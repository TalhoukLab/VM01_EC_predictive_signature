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
library(gsubfn)
library(dplyr)
library(DrImpute)
library(msa)

library(seqinr)

dataPrep <- function(raw_otus, raw_tax, phylo_tree, feature_table){
  rownames(raw_otus) <- raw_otus$X.OTU.ID
  raw_otus = subset(raw_otus, select = -c(X.OTU.ID))
  otus <- mutate_all(raw_otus,function(x) as.numeric(as.character(x)))
  names(otus) = gsub(pattern = "_", replacement = "", x = names(otus))
  rownames(otus) <- toupper(rownames(otus))
  raw_tax <- raw_tax[-1, ]
  raw_tax <-  as.data.frame(tidyr::separate_wider_delim(raw_tax, cols = V1, delim = "\t", names = c("OTUID", "domain")))
  colnames(raw_tax) <- c("OTUID", "domain", "phylum", "class", "order", "family", "genus", "species")
  raw_tax <- raw_tax[!raw_tax$phylum=="", ]
  raw_tax$domain<-gsub("d__","",  raw_tax$domain)
  raw_tax$phylum<-gsub("p__","",  raw_tax$phylum)
  raw_tax$class<-gsub("c__","",  raw_tax$class)
  raw_tax$order<-gsub("o__","",  raw_tax$order)
  raw_tax$family<-gsub("f__","",  raw_tax$family)
  raw_tax$genus<-gsub("g__","",  raw_tax$genus)
  raw_tax$species<-gsub("s__","",  raw_tax$species)
  rownames(raw_tax) <- raw_tax$OTUID
  raw_tax <- dplyr::select(raw_tax, select = -c(OTUID))  
  feature_table$histology[feature_table$histology == "ACH"] <- "EC"
  rownames(feature_table) <- feature_table$sraID
  otus_table <- otu_table(otus, taxa_are_rows = TRUE)
  rooted_tree <- phangorn::midpoint(phylo_tree)
  tree_phy <- phy_tree(rooted_tree)
  taxa_names(tree_phy) <- toupper(taxa_names(tree_phy))
  tax_table_phy = tax_table(as.matrix(raw_tax))
  samples = sample_data(feature_table)
  phyloseq_pre <- phyloseq(otus_table, tax_table_phy, samples, tree_phy)
  phyloseq_pre_un <- subset_taxa(phyloseq_pre, is.na(phylum)==FALSE)
  phyloseq_obj_ne <- prune_samples(sample_sums(phyloseq_pre_un) > 1, phyloseq_pre_un)
  
  
  phylo_raw <- phyloseq(otu_table(otus, taxa_are_rows=TRUE), tax_table(phyloseq_obj_ne), 
                        sample_data(phyloseq_obj_ne), phy_tree(phyloseq_obj_ne))
  
  phylo_raw <- subset_taxa(phylo_raw, phylum!="Cyanobacteria/Chloroplast")
  
  otus <- abundances(phylo_raw)
  to_keep <- which((rowSums(otus))>2)
  otus <- otus[to_keep, ]
  otus <- otus + 1
  otu_table(phylo_raw) <- otu_table(otus, taxa_are_rows = TRUE)
  otu.tss <- log10(data.frame((apply(otus, 2, function(x){x/sum(x)}))))
  phylo_tss <- phyloseq(otu_table(otu.tss, taxa_are_rows=TRUE), tax_table(phylo_raw), 
                        sample_data(phylo_raw), phy_tree(phylo_raw))
  
  res <- list(raw = phylo_raw, tss = phylo_tss)
  return(res)
}
