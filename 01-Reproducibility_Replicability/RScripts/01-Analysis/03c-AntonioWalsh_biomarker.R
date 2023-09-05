## This script reproduces DET analysis for Antonio and Walsh. 
# Identified DET at phylum, family, genus, and species 
# Filter rare taxa prevalent at less than 10% of samples or taxa with a maximum proportion less than 0.2% 
# Fit a generalized linear model with over-dispersed Poisson distribution to the count data 
# Estimated library sizes (sequencing depth) using the geometric mean of paired rations (GMPR) normlaized method 
# log of GMPR size factors were used in the Poisson model as an offset to account for variable library sizes 
# Data was winsorized at 97% upper quatile to reduce the possible impact of outliers on parameter estimates before fitting the model

library(GUniFrac)
library(data.table)
library(MatrixGenerics)

cohorts <- c("Chao", "Antonio", "Gressel", "Tsementzi", "Walsh")
levels <- c("phylum", "class", "order", "family", "genus")

for(cohort in cohorts){
  for(level in levels){
    print(cohort)
    print(level)
    ## Collapse at level of interest
    if(level=="OTU"){
      phylo_use <- eval(parse(text = paste0(cohort, "_", "Antonio", "phyloseq_tree_raw")))
      agg_phylo <- phylo_use
    } else {
      phylo_use <- eval(parse(text = paste0(cohort, "_", "Antonio", "phyloseq_tree_raw")))
      agg_phylo <- aggregate_taxa(phylo_use, level = level)
    }
    
    ## Winsorize (adding even if no outliers)
    count_data <- as.data.frame(t(as.matrix(agg_phylo@otu_table@.Data)))
    to_keep <- which((colSums(count_data))>0)
    count_data_nonEmpty <- count_data[, to_keep]
    dep <- colSums(count_data_nonEmpty)
    
    count_data.p <- t(t(count_data_nonEmpty)/dep)
    count_data.p <- apply(count_data.p, 1, function(x) {
      cutoff <- quantile(x, 0.97)
      x[x >= cutoff] <- cutoff
      x
    }
    )
    count_data_nonEmpty <- t(round(count_data.p * dep))
    to_keep <- which((colSums(count_data_nonEmpty))>0)
    count_data_nonEmpty <- count_data_nonEmpty[, to_keep]
    ## Filter rare taxa 
    prop <- t(t(count_data_nonEmpty) / colSums(count_data_nonEmpty))
    prop <- prop[rowSums(prop!=0) > 0.1 * ncol(prop), ,drop=FALSE]	
    count_data_nonEmpty <- count_data_nonEmpty[rownames(prop), , drop=FALSE]
    
    prop <- prop[rowMaxs(prop) > 0.002, , drop=FALSE]	
    count_data_nonEmpty <- count_data_nonEmpty[rownames(prop), , drop=FALSE]
    
    ## GMPR 
    count_data_mat <- t(data.frame(count_data_nonEmpty))
    dep <- data.frame(GMPR(count_data_mat))
    
    ## Poisson model 
    count_data_fit <- data.frame(t(count_data_mat))
    sraID <- data.frame(rownames(count_data_fit))
    rownames(count_data_fit) <- NULL
    count_data_fit <- cbind(count_data_fit, dep = dep$GMPR.count_data_mat., sraID = sraID$rownames.count_data_fit.)
    sample_data <- data.frame(sample_data(agg_phylo))
    merged <- merge(count_data_fit, sample_data, by = "sraID")
    if(cohort == "Antonio") {
      merged <- subset(merged, select = -c(sraID, cohort, Stage, pH, X))
      merged$BMI <- as.numeric(merged$BMI)
      merged$age <- as.numeric(merged$age)
      merged$pH <- as.numeric(as.factor(merged$pHRecoded))
      merged$histology <- as.numeric(merged$histology)
      merged$menopausal.status <- as.numeric(as.factor(merged$menopausal.status))
      model <- glm(histology ~ ., data = merged, family = poisson)
    }
    if(cohort == "Walsh"){
      merged <- subset(merged, select = -c(sraID, cohort, stage, ethnicity, pH, X))
      merged$BMI <- as.numeric(merged$BMI)
      merged$age <- as.numeric(merged$age)
      merged$pHRecoded <- as.numeric(as.factor(merged$pHRecoded))
      merged$histology <- as.numeric(merged$histology)
      merged$menopausal.status <- as.numeric(as.factor(merged$menopausal.status))
      merged$ethnicityRecoded <- as.numeric(as.factor(merged$ethnicityRecoded))
      model <- glm(histology ~ ., data = merged, family = poisson)
    }
    if(cohort == "Tsementzi"){
      merged <- subset(merged, select = -c(sraID, cohort, stage, menopausal.status, pH, X))
      merged$BMI <- as.numeric(merged$BMI)
      merged$age <- as.numeric(merged$age)
      merged$pH <- as.numeric(as.factor(merged$pHRecoded))
      merged$histology <- as.numeric(merged$histology)
      merged$ethnicity <- as.numeric(as.factor(merged$ethnicity))
      model <- glm(histology ~ ., data = merged, family = poisson)
    }
    if(cohort == "Chao"){
      merged <- subset(merged, select = -c(sraID, cohort, stage,menopausal.status))
      merged$age <- as.numeric(merged$age)
      merged$histology <- as.numeric(merged$histology)
      model <- glm(histology ~ ., data = merged, family = poisson)
    }
    if(cohort == "Gressel"){
      merged <- subset(merged, select = -c(sraID, cohort, X))
      merged$histology <- as.numeric(merged$histology)
      model <- glm(histology ~ ., data = merged, family = poisson)
    }

    print(summary(model))
    results <- data.frame(summary(model)$coefficients)
    results$q <- p.adjust(results$Pr...z..,method="BH")
    #results_write <- results[results$Pr...z.. <= 0.05, ]
    write.csv(results, paste0("~/Desktop/", cohort, "_level_", level, ".csv"), row.names = TRUE)
    #if(nrow(results_write)>0){
     # write.csv(results_write, paste0("~/Desktop/", cohort, "_level_", level, ".csv"), row.names = FALSE)
    #}
  }
}

