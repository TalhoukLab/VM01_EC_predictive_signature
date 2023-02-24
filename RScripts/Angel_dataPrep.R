cohorts <- c("Angel", "Antonio", "Gressel", "Tsementzi", "Walsh")
for(cohort in cohorts){
  print(cohort)
  raw_table <- read.delim(file.path(paste0('~/Desktop/R&R/vaginalMicrobiome/Angel_pipeline/results/', cohort, '_cohort'), 'zotutab_raw.txt'), sep = '\t')
  raw_otus <- as.data.frame(as.matrix((raw_table)))
  rownames(raw_otus) <- raw_otus$X.OTU.ID
  raw_otus = subset(raw_otus, select = -c(X.OTU.ID))
  otus <- mutate_all(raw_otus,function(x) as.numeric(as.character(x)))
  names(otus) = gsub(pattern = "_", replacement = "", x = names(otus))
  
  #metaSeqObject <- newMRexperiment(otus) 
  #CSS <- cumNorm(metaSeqObject, p=cumNormStat(metaSeqObject))
  #outs_CSS = data.frame(MRcounts(CSS, norm=TRUE, log=T))
  
  tax_table <- read.delim(file.path(paste0('~/Desktop/R&R/vaginalMicrobiome/Angel_pipeline/results/', cohort, '_cohort'), 'sintax.txt'), header = FALSE, sep = "\t")
  tax_table$V2 <- gsub("\\s*\\([^\\)]+\\)","",as.character(tax_table$V2))
  tax_table <- separate(data = tax_table, col = V2 , into = c("kingdom", "phylum", "class", "order", "family", "genus"), sep = "\\,")
  rownames(tax_table) <- tax_table$V1
  tax_table[is.na(tax_table)] <- " "
  tax_table = subset(tax_table, select = -c(V3, V4, V1))
  tax_table$kingdom<-gsub("d:","",as.character(tax_table$kingdom))
  tax_table$phylum <- gsub("p:", "", gsub("\\[", "",  gsub("\\]", "" ,as.character(tax_table$phylum))))
  tax_table$class <- gsub("c:", "", gsub("\\[", "",  gsub("\\]", "" ,as.character(tax_table$class))))
  tax_table$order <- gsub("o:", "", gsub("\\[", "",  gsub("\\]", "" ,as.character(tax_table$order))))
  tax_table$family <- gsub("f:", "", gsub("\\[", "",  gsub("\\]", "" ,as.character(tax_table$family))))
  tax_table$genus <- gsub("g:", "", gsub("\\[", "",  gsub("\\]", "" ,as.character(tax_table$genus))))
  tax_table <- mutate_all(tax_table, .funs=tolower)
  feature_table <- read.delim(paste0('~/Desktop/R&R/vaginalMicrobiome/00-helperfiles/', cohort,'FT.csv'), header = TRUE, sep = ",")
  feature_table$histology[feature_table$histology == "ACH"] <- "EC"
  rownames(feature_table) <- feature_table$sraID
  phylo_tree <- read_tree(file.path(paste0('~/Desktop/R&R/vaginalMicrobiome/Angel_pipeline/results/', cohort, '_cohort'),'zotus.tree'))
  otus_table <- otu_table(otus, taxa_are_rows = TRUE)
  tax_table_phy = tax_table(as.matrix(tax_table))
  samples = sample_data(feature_table)
  tree_phy <- phy_tree(phylo_tree)
  ## Phyloseq object without tree 
  phyloseq_pre <- phyloseq(otus_table, tax_table_phy, samples, tree_phy)
  phyloseq_pre_un <- subset_taxa(phyloseq_pre, kingdom!="Unclassified")
  phyloseq_obj_ne <- prune_samples(sample_sums(phyloseq_pre_un) > 1, phyloseq_pre_un)
  phyloseq_com <- microbiome::transform(phyloseq_obj_ne, 'compositional')
  otu_table <- data.frame(otu_table(phyloseq_com)*100)
  otu_table[otu_table == 0] <- 0.001
  phyloseq_obj <- phyloseq(otu_table(otu_table,  taxa_are_rows = TRUE), tax_table(phyloseq_com), sample_data(phyloseq_com), phy_tree(phyloseq_com))
  phyloseq_lg <- microbiome::transform(phyloseq_obj, 'log10')
  
  assign(paste0(cohort, "_Angelphyloseq_tree"), phyloseq_lg,.GlobalEnv)
  assign(paste0(cohort, "_Angelphyloseq_tree_raw"), phyloseq_obj_ne,.GlobalEnv)
}
