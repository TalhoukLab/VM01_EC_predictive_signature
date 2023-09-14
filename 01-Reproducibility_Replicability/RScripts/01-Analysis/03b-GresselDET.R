## Gressel pipeline DET 
library(ANCOMBC)
library(tidyverse)

cohorts <- c("Antonio", "Chao", "Gressel", "Tsementzi", "Walsh")
levels <- c("phylum", "class", "order", "family", "genus")
for(cohort in cohorts){
  for(level in levels){
    phylo_use_raw <- eval(parse(text = paste0(cohort, "_Gresselphyloseq_tree")))
    if(cohort == "Antonio"){
      out = ancombc(data = phylo_use_raw, assay_name = "counts", 
                    if(level == "OTU")  tax_level <- NULL else  tax_level =  level, 
                    formula = "menopausal.status + age + BMI + pHRecoded + histology", 
                    p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
                    group = "histology", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
                    n_cl = 1, verbose = TRUE)
    }
    if(cohort == "Chao"){
      out = ancombc(data = phylo_use_raw, assay_name = "counts", 
                    if(level == "OTU")  tax_level <- NULL else  tax_level =  level, 
                    formula = "age + histology", 
                    p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
                    group = "histology", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
                    n_cl = 1, verbose = TRUE)
    }
    if(cohort == "Gressel"){
      out = ancombc(data = phylo_use_raw, assay_name = "counts", 
                    if(level == "OTU")  tax_level <- NULL else  tax_level =  level, 
                    formula = "histology", 
                    p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
                    group = "histology", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
                    n_cl = 1, verbose = TRUE)
    }
    if(cohort == "Tsementzi"){
      out = ancombc(data = phylo_use_raw, assay_name = "counts", 
                    if(level == "OTU")  tax_level <- NULL else  tax_level =  level, 
                    formula = "BMI + age + pHRecoded + histology + ethnicityRecode", 
                    p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
                    group = "histology", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
                    n_cl = 1, verbose = TRUE)
      
    }
    if(cohort == "Walsh "){
      out = ancombc(data = phylo_use_raw, assay_name = "counts", 
                    if(level == "OTU")  tax_level <- NULL else  tax_level =  level, 
                    formula = "menopausal.status + BMI + age + pHRecoded + histology + ethnicityRecoded", 
                    p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
                    group = "histology", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
                    n_cl = 1, verbose = TRUE)
    }
    assign(paste0(cohort, "_", level, "_ancombc"), out,.GlobalEnv)
    res_ancombc <- as.data.frame(cbind(taxon = out$res$diff_abn$taxon, diff_abun = out$res$diff_abn$histologyEC, lefc = out$res$lfc$histologyEC, 
                                       se = out$res$se$histologyEC, W = out$res$W$histologyEC))
    rownames(res_ancombc) <-  gsub(paste0(level,":"), '', out$res$lfc$taxon)
    res_ancombc_fil <- res_ancombc[(res_ancombc$W>10 & res_ancombc$diff_abun==TRUE), ]
    res_ancombc_fil <- res_ancombc_fil[res_ancombc_fil$taxon %like% level, ]
    assign(paste0(cohort, "_", level, "_ancombc_res"), res_ancombc,.GlobalEnv)
    if(nrow(res_ancombc_fil)>0){
      write.csv(res_ancombc_fil, paste0("~/Desktop/", cohort, "_level_", level, ".csv"), row.names = TRUE)
    }
  }
}