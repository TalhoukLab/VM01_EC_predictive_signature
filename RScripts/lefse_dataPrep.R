## This script will format the input for Lefse
pipeline <- "Angel"
levels <- c("OTU", "genus", "family", "order")
cohorts <- c("Angel", "Antonio", "Gressel", "Tsementzi", "Walsh")
for(cohort in cohorts){
  for(level in levels){
    ## OTU level
    if(level=="OTU"){
      data <- data.frame(otu_table(eval(parse(text = paste0(cohort, "_", pipeline, "phyloseq_tree")))))
      taxa_names <- data.frame(tax_table(eval(parse(text = paste0(cohort, "_", pipeline, "phyloseq_tree")))))
      taxa_names$merged <- paste(rownames(taxa_names), taxa_names$kingdom, taxa_names$phylum, taxa_names$class,
                                 taxa_names$order, taxa_names$family, taxa_names$genus, sep = "|")
      rownames(data) <- taxa_names$merged
    } else {
      agg_taxa <- aggregate_taxa(eval(parse(text = paste0(cohort, "_", pipeline, "phyloseq_tree"))), level = level)
      data <- data.frame(otu_table(agg_taxa))
      taxa_names <- data.frame(tax_table(agg_taxa))
      rownames(data) <- taxa_names$unique
    }
    flip_data <- data.frame(t(data), check.names = F)
    if(cohort=="Antonio" | cohort=="Tsementzi" | cohort=="Walsh"){
      FT <- as.data.frame(sample_data(eval(parse(text = paste0(cohort, "_", pipeline, "phyloseq_tree")))))
      FT <- data.frame(FT) %>% 
        mutate(Obesity = if_else(BMI >=30, "Obese", "Non-obese"))
      FT$pH <- as.numeric(FT$pH)
      covariates = FT[, c('Age', 'histology', 'pH', 'Obesity', 'sraID')]
      data_full <- cbind(covariates, flip_data)
    } else{
      FT <- as.data.frame(sample_data(eval(parse(text = paste0(cohort, "_", pipeline, "phyloseq_tree")))))
      covariates = FT[, c('histology', 'sraID')]
      data_full <- cbind(covariates, flip_data)
    }
    data_full <- na.omit(data_full)
    data_full1 <- data.frame(t(data_full), check.names = F)
    write.table(data_full1, paste0("~/Desktop/lefse_", level, "_", pipeline, "pipeline_", cohort, ".txt"), quote=FALSE, sep = "\t", col.names = FALSE)
  }
}

