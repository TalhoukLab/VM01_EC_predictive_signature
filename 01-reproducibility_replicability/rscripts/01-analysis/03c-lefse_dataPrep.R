## This script will format the input for Lefse
pipeline <- "Chao"
levels <- c("phylum", "class", "order", "family", "genus")
cohorts <- c("Antonio","Chao", "Gressel", "Tsementzi", "Walsh")
for(cohort in cohorts){
  for(level in levels){
    agg_taxa  <- subset_taxa(eval(parse(text = paste0(cohort, "_", pipeline, 'phyloseq_tree_raw'))), 
                                c(eval(parse(text=level))!=" " & eval(parse(text=level)) !=""))
    agg_taxa <- aggregate_taxa(agg_taxa, level = level)
    data <- data.frame(otu_table(agg_taxa))
    taxa_names <- data.frame(agg_taxa@tax_table)
    rownames(data) <- taxa_names$unique
    flip_data <- data.frame(t(data), check.names = F)
    if(cohort=="Antonio"){
      FT <- as.data.frame(sample_data(eval(parse(text = paste0(cohort, "_", pipeline, "phyloseq_tree_raw")))))
      FT$age <- as.numeric(FT$age)
      FT$BMI <- as.numeric(FT$BMI)
      FT$pHRecoded <- as.numeric(as.factor(FT$pHRecoded))
      FT$menopausal.status <- as.numeric(as.factor(FT$menopausal.status))
      covariates = FT[, c('histology', 'sraID', 'pHRecoded', 'BMI', 'age', 'menopausal.status')]
      data_full <- cbind(covariates, flip_data)
    }
    if(cohort=="Walsh"){
      FT <- as.data.frame(sample_data(eval(parse(text = paste0(cohort, "_", pipeline, "phyloseq_tree_raw")))))
      FT$age <- as.numeric(FT$age)
      FT$BMI <- as.numeric(FT$BMI)
      FT$pHRecoded <- as.numeric(as.factor(FT$pHRecoded))
      FT$menopausal.status <- as.numeric(as.factor(FT$menopausal.status))
      FT$ethnicityRecoded <- as.numeric(as.factor(FT$ethnicityRecoded))
      covariates = FT[, c('histology', 'sraID', 'age','pHRecoded', 'BMI','menopausal.status', 'ethnicityRecoded' )]
      data_full <- cbind(covariates, flip_data)
    } 
    if(cohort == "Tsementzi" ){
      FT <- as.data.frame(sample_data(eval(parse(text = paste0(cohort, "_", pipeline, "phyloseq_tree_raw")))))
      FT$age <- as.numeric(FT$age)
      FT$BMI <- as.numeric(FT$BMI)
      FT$pHRecoded <- as.numeric(as.factor(FT$pHRecoded))
      FT$ethnicityRecode <- as.numeric(as.factor(FT$ethnicityRecode))
      covariates = FT[, c('histology','sraID', 'age','pHRecoded', 'BMI', 'ethnicityRecode')]
      data_full <- cbind(covariates, flip_data)
    } 
    if(cohort == "Chao"){
      FT <- as.data.frame(sample_data(eval(parse(text = paste0(cohort, "_", pipeline, "phyloseq_tree_raw")))))
      covariates = FT[, c('histology', 'sraID', 'age')]
      FT$age <- as.numeric(FT$age)
      data_full <- cbind(covariates, flip_data)
    }
    if(cohort == "Gressel"){
      FT <- as.data.frame(sample_data(eval(parse(text = paste0(cohort, "_", pipeline, "phyloseq_tree_raw")))))
      covariates = FT[, c('histology', 'sraID')]
      data_full <- cbind(covariates, flip_data)
    }
    data_full <- na.omit(data_full)
    data_full1 <- data.frame(t(data_full), check.names = F)
    write.table(data_full1, paste0("../vaginalMicrobiome/01-Reproducibility_Replicability/Results/03-DET/Chao_pipeline/dataPrep/lefse_", level, "_", pipeline, "pipeline_", cohort, ".txt"), quote=FALSE, sep = "\t", col.names = FALSE)
  }
}

