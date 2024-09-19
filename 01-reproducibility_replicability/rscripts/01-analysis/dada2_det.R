## Gressel pipeline DET 
library(ANCOMBC)
library(tidyverse)
library(DescTools)
limit <- 2
cohorts <- c("Antonio", "Chao", "Gressel", "Tsementzi", "Walsh")
level <- "species"
pipelines <- c("dada2")
for(cohort in cohorts){
  for(pipeline in pipelines){
    phylo_use_raw <- eval(parse(text = paste0(cohort, "_", pipeline, "phyloseq_tree_raw")))
    if(cohort == "Antonio"){
      out = ancombc(data = phylo_use_raw, assay_name = "counts", 
                    tax_level =  level, 
                    formula = "menopausal.status + age + BMI + pHRecoded + histology", 
                    p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
                    group = "histology", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
                    n_cl = 1, verbose = TRUE)
    }
    if(cohort == "Chao"){
      out = ancombc(data = phylo_use_raw, assay_name = "counts", 
                    tax_level =  level, 
                    formula = "age + histology", 
                    p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
                    group = "histology", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
                    n_cl = 1, verbose = TRUE)
    }
    if(cohort == "Gressel"){
      out = ancombc(data = phylo_use_raw, assay_name = "counts", 
                    tax_level =  level, 
                    formula = "histology", 
                    p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
                    group = "histology", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
                    n_cl = 1, verbose = TRUE)
    }
    if(cohort == "Tsementzi"){
      out = ancombc(data = phylo_use_raw, assay_name = "counts", 
                    tax_level =  level, 
                    formula = "BMI + age + pHRecoded + histology + ethnicityRecode", 
                    p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
                    group = "histology", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
                    n_cl = 1, verbose = TRUE)
      
    }
    if(cohort == "Walsh"){
      sample_data <- data.frame(sample_data(phylo_use_raw))
      sample_data <- na.omit(sample_data)
      sample_data$age <- as.numeric(sample_data$age)
      sample_data$BMI <- as.numeric(sample_data$BMI)
      sample_data$menopausal.status <- as.factor(sample_data$menopausal.status)
      phylo_use_raw <- phyloseq(otu_table(phylo_use_raw), tax_table(phylo_use_raw), sample_data(sample_data))
      out = ancombc(data = phylo_use_raw, assay_name = "counts", 
                    tax_level =  level, 
                    formula = "menopausal.status + BMI + age + pHRecoded + histology + ethnicityRecoded", 
                    p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
                    group = "histology", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
                    n_cl = 1, verbose = TRUE)
    }
    assign(paste0(cohort, "_", pipeline, "_ancombc"), out,.GlobalEnv)
    res_ancombc <- merge(out$res$lfc, out$res$diff_abn, by = "taxon", suffixes = c("_lfc","_diff_abun") )
    res_ancombc$cohort <- cohort 
    res_ancombc <- res_ancombc %>% 
                            filter_at(vars(ends_with('_diff_abun')),  any_vars(. == TRUE))
    res = out$res
    tab_W <- res$W
    library(dplyr)
    tab_diff <- res$diff_abn
    
    tab_diff_longer <- tab_diff %>% 
      pivot_longer(-taxon)%>% 
      filter(value == TRUE) 
    
    lfc_longer <- res$lfc %>% 
      pivot_longer(-taxon) 
    
    diff_lfc <- merge(tab_diff_longer, lfc_longer, 
                      by = c("taxon", "name"), suffixes = c("_lfc","_diff_abun") )

    diff_lfc$expression <- with(diff_lfc, 
                            ifelse(value_diff_abun>0, 'positive', 'negative'))
    
    diff_lfc$cohort <- cohort
    
    W_diff <- merge(tab_W, tab_diff, by = "taxon")
   # W_diff <- W_diff %>% 
  #    filter_at(vars(ends_with('.x')),  any_vars(abs(.) >=limit))
    tab_diff <- tab_diff %>% filter(taxon %in% W_diff$taxon)
    assign(paste0(cohort, "_", pipeline, "_ancombc_res"), res_ancombc,.GlobalEnv)
    assign(paste0(cohort, "_", pipeline, "_networkplot"), diff_lfc,.GlobalEnv)
    
  }
}
