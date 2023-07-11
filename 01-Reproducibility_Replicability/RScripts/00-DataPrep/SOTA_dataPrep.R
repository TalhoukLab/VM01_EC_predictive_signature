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
library(rRDP)
library(rRDPData)
library(seqinr)

cohorts <- c("Antonio", "Chao", "Gressel", "Tsementzi", "Walsh")
for(cohort in cohorts){
  raw_otus <- read.delim(file.path(paste0('~/Desktop/thesis/vaginalMicrobiome/01-Reproducibility_Replicability/SOTA_pipeline/results/', cohort, '_cohort/all.otutab.txt')))
  rownames(raw_otus) <- raw_otus$X.OTU.ID
  raw_otus = subset(raw_otus, select = -c(X.OTU.ID))
  otus <- mutate_all(raw_otus,function(x) as.numeric(as.character(x)))
  names(otus) = gsub(pattern = "_", replacement = "", x = names(otus))
  rep_seqs <- readDNAStringSet(file.path(paste0('~/Desktop/thesis/vaginalMicrobiome/01-Reproducibility_Replicability/SOTA_pipeline/results/', cohort, '_cohort/all.otus.fasta')))
  names(rep_seqs) <- sapply(strsplit(names(rep_seqs), ";"), "[", 1)
  pred <- predict(rdp(), rep_seqs)
  
  pred_copy <- pred
  pred_copy <- pred_copy[is.na(pred_copy$phylum)==FALSE, ]
  pred_copy <- pred_copy%>%
    mutate(genus= ifelse(is.na(genus), family, genus)) %>%
    mutate(genus= ifelse(is.na(genus), order, genus)) %>%
    mutate(genus= ifelse(is.na(genus), class, genus)) %>%
    mutate(genus= ifelse(is.na(genus), phylum ,genus))


  unique_genus <- unique(pred_copy$genus)
  otu_filter <- otus
  #otu_filter <- filter_taxa(otu_filter, 0.1, 5)
  #otu_filter <- otu_filte[, colSums(otu_filter != 0) > 0]

  otu_copy <- otu_filter
  pred_copy1 <- pred_copy
  for(genus_use in unique_genus){
    print(genus_use)
    otu_ids <- rownames(subset(pred_copy, genus==genus_use))
    if(length(otu_ids)>1){
      seqs_msa <- rep_seqs[otu_ids]
      msa_align <- msa(seqs_msa, "ClustalOmega", type = "dna")
      msa_align2 <- msaConvert(msa_align, type="seqinr::alignment")
      d <- as.matrix(100*(1-(dist.alignment(msa_align2, "similarity")^2)))
      d[upper.tri(d)] <- NA
      d_mat <- melt(d)
      d_matFil <- as.matrix(subset(d_mat, value >= 97))
      if(nrow(d_matFil)>1){
        j <- 1
        dict_list <- list()
        dict_list[1] <- "NULL VALUE"
        for(i in 1:nrow(d_matFil)){
          index1 = which(sapply(dict_list, function(x) d_matFil[i, 'Var1'] %in% x))
          index2 = which(sapply(dict_list, function(x) d_matFil[i, 'Var2'] %in% x))
          if((length(index1)==1 && length(index2)==1) && (index1!=index2)){
            dict_list[[index1]] <- c(dict_list[[index1]], dict_list[[index2]])
            dict_list <- dict_list[-index2]
            j <- j -1
            idex_use <- which(sapply(dict_list, function(x) d_matFil[i, 'Var1'] %in% x))
            dict_list[[idex_use]] <- append(dict_list[[idex_use]], d_matFil[i, 'Var1'])
            dict_list[[idex_use]] <- append(dict_list[[idex_use]], d_matFil[i, 'Var2'])
          } else if(length(index1)==0 && length(index2)==0){
            dict_list[j] <- list(c(d_matFil[i, 'Var1'], d_matFil[i, 'Var2']))
            j <- j+1
          } else if(length(index2)==1 && length(index1)==0){
            dict_list[[index2]] <- append(dict_list[[index2]], d_matFil[i, 'Var1'])
          } else {
            dict_list[[index1]] <- append(dict_list[[index1]], d_matFil[i, 'Var2'])
          }
        }
        for(i in 1:length(dict_list)){
          to_merge <- unique(unlist(dict_list[i]))
          otu_copy <- aggregate(otu_copy, list(Group=replace(rownames(otu_copy),rownames(otu_copy) %in%
                                                               to_merge, to_merge[1])), sum)
          rownames(otu_copy) <- otu_copy$Group
          otu_copy = subset(otu_copy, select = -c(Group))
          to_remove <- to_merge[-1]
          pred_copy1 <- pred_copy1[!(row.names(pred_copy1) %in% to_remove),]
          pred_copy1[to_merge[1], ]$genus <- paste0(pred_copy1[to_merge[1], ]$genus, "_merged_", i)
        }
        to_merge <- unlist(dict_list)
        remaining <- which(otu_ids %in% to_merge ==FALSE)
        if(length(remaining)>0){
          for(i in 1:length(remaining)){
            pred_copy1[otu_ids[remaining[i]], ]$genus <- paste0(pred_copy1[otu_ids[remaining[i]], ]$genus, "_unmerged_", i)
          }
        }

      }
    }

  }
  #otu_copy <- otus
  metaSeqObject <- newMRexperiment(otu_copy) 
  CSS <- cumNorm(metaSeqObject, p=cumNormStat(metaSeqObject))
  outs_CSS = data.frame(MRcounts(CSS, norm=TRUE, log=F))
  #outs_CSS_lg <- log(outs_CSS + 1)
 
  feature_table <- read.delim(file.path(paste0('~/Desktop/thesis/vaginalMicrobiome/01-Reproducibility_Replicability/00-helperfiles/', cohort, 'FT.csv')), header = TRUE, sep = ",")
  #ids <- c(read.delim("~/Desktop/TAW_train.txt", sep = "\t", header = FALSE))
  #feature_table <- vm.metadata[ids$V1,]
  feature_table$histology[feature_table$histology == "ACH"] <- "EC"
  rownames(feature_table) <- feature_table$sraID
  otus_table <- otu_table(otu_copy, taxa_are_rows = TRUE)
  #pred_copy1 <- pred
  tax_table_phy = tax_table(as.matrix(pred_copy1))
  samples = sample_data(feature_table)
  phyloseq_pre <- phyloseq(otus_table, tax_table_phy, samples)
  phyloseq_pre_un <- subset_taxa(phyloseq_pre, is.na(phylum)==FALSE)
  #phyloseq_pre_un1 <- subset_taxa(phyloseq_pre_un, phylum!="cyanobacteria/chloroplast")
  phyloseq_obj_ne <- prune_samples(sample_sums(phyloseq_pre_un) > 1, phyloseq_pre_un)
  #otu_table_clr <- as.data.frame(otu_table(phyloseq_obj_ne))
  #imputed_clr <- DrImpute(as.matrix(otu_table_clr), ks = 2)
  #colnames(imputed_clr) <- colnames(otu_table_clr)
  #phylo_imputed <- phyloseq(otu_table(imputed_clr, taxa_are_rows=TRUE), tax_table(phyloseq_obj_ne), 
  #                          sample_data(phyloseq_obj_ne)) 
  #assign(paste0(cohort, "_InHousephyloseq_tree"), phylo_imputed,.GlobalEnv)
  phylo_raw <- phyloseq(otu_table(otu_copy, taxa_are_rows=TRUE), tax_table(phyloseq_obj_ne), 
                       sample_data(phyloseq_obj_ne))
  
  assign(paste0(cohort, "_SOTAphyloseq_tree_raw"),phylo_raw ,.GlobalEnv)
}

