## This script will format the input for Lefse
pipeline <- "Tsementzi"
levels <- c("OTU")
cohorts <- c("Antonio","Chao", "Gressel", "Tsementzi", "Walsh")
for(cohort in cohorts){
  for(level in levels){
    ## OTU level
    if(level=="OTU"){
      rel_phyloseq <- microbiome::transform(eval(parse(text=paste0(cohort, "_", pipeline, "phyloseq_tree_raw"))), "compositional")
      data <- data.frame(otu_table(rel_phyloseq))
      taxa_names <- data.frame(tax_table(rel_phyloseq))
      taxa_names$merged <- paste(rownames(taxa_names), taxa_names$kingdom, taxa_names$phylum, taxa_names$class,
                                 taxa_names$order, taxa_names$family, taxa_names$genus, sep = "|")
      rownames(data) <- taxa_names$merged
    } else {
      agg_taxa <- aggregate_taxa(eval(parse(text = paste0(cohort, "_", pipeline, "phyloseq_tree_raw"))), level = level)
      data <- data.frame(otu_table(agg_taxa))
      taxa_names <- data.frame(tax_table(agg_taxa))
      rownames(data) <- taxa_names$unique
    }
    flip_data <- data.frame(t(data), check.names = F)
    if(cohort=="Antonio"| cohort=="Walsh"){
      FT <- as.data.frame(sample_data(eval(parse(text = paste0(cohort, "_", pipeline, "phyloseq_tree_raw")))))
      #FT <- data.frame(FT) %>% 
      #  mutate(Obesity = if_else(BMI >=30, "Obese", "Non-obese"))
      FT$pH <- as.numeric(FT$pH)
      FT$BMI <- as.numeric(FT$BMI)
      covariates = FT[, c('Age', 'histology', 'pH', 'BMI', 'sraID')]
      data_full <- cbind(covariates, flip_data)
    } else if ( cohort=="Tsementzi" ){
      FT <- as.data.frame(sample_data(eval(parse(text = paste0(cohort, "_", pipeline, "phyloseq_tree_raw")))))
      #FT <- data.frame(FT) %>% 
      # mutate(Obesity = if_else(BMI >=30, "Obese", "Non-obese"))
      FT$pH <- as.numeric(FT$pH)
      FT$BMI <- as.numeric(FT$BMI)
      covariates = FT[, c('age', 'histology', 'pH', 'BMI', 'sraID')]
      #covariates = FT[, c('histology', 'sraID')]
      data_full <- cbind(covariates, flip_data)
    } else{
      FT <- as.data.frame(sample_data(eval(parse(text = paste0(cohort, "_", pipeline, "phyloseq_tree_raw")))))
      covariates = FT[, c('histology', 'sraID')]
      data_full <- cbind(covariates, flip_data)
    }
    data_full <- na.omit(data_full)
    data_full1 <- data.frame(t(data_full), check.names = F)
    write.table(data_full1, paste0("~/Desktop/lefse_", level, "_", pipeline, "pipeline_", cohort, ".txt"), quote=FALSE, sep = "\t", col.names = FALSE)
  }
}

