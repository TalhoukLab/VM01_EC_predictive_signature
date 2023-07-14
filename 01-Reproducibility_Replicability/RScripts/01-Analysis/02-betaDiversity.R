library(patchwork)
library(microViz)
## Chao beta diversity - only use unweighted uniFrac distance with an PCA plot. 
cohorts <- c("Antonio","Chao", "Gressel", "Tsementzi", "Walsh")

mul_var_analysis <- function(cohort, dist_mat, phylo_to_use){
  print(phylo_to_use)
  if(cohort == "Chao"){
    test_mul <- adonis2(dist_mat ~ sample_data(phylo_to_use)$histology + 
                          sample_data(phylo_to_use)$age, 
                        permutations = 1000, na.action = "na.omit")
  } 
  if(cohort == "Antonio" || cohort == "Tsementzi" || cohort == "Walsh"){
    test_mul <- adonis2(dist_mat ~ sample_data(phylo_to_use)$histology + 
                          sample_data(phylo_to_use)$age + 
                          sample_data(phylo_to_use)$BMI + 
                          sample_data(phylo_to_use)$pHRecoded,
                        permutations = 1000, na.action = "na.omit")
  }
  if(cohort == "Gressel"){
    test_mul <- adonis2(dist_mat ~ sample_data(phylo_to_use)$histology,
                        permutations = 1000, na.action = "na.omit")
  }
  return(test_mul)
}


sink("~/Desktop/Chao_qual.txt")
for(cohort in cohorts){
  print(cohort)
  phylo_use <- eval(parse(text = paste0(cohort, "_Chaophyloseq_tree")))
  dist_obj <- phylo_use %>%
    tax_transform("identity", rank = "unique") %>%
    dist_calc("unifrac") 
  
  p1 <- phylo_use %>%
    tax_transform("identity", rank = "unique") %>%
    dist_calc("unifrac") %>%
    ord_calc("PCA")%>%
    ord_plot(axes = c(1, 2), color = "histology", size = 2.5) +
    scale_color_manual(values=c("yellowgreen", "tomato3"))
  
  mul_var_res <- mul_var_analysis(cohort = cohort, dist_mat = dist_obj@dist, phylo_to_use = phylo_use)
  assign(paste0(cohort, "_Chao_beta"), p1,.GlobalEnv)
  print(mul_var_res)
}
sink()

## Antonio beta diversity - only use unweighted unifrac using NMDS ordination 
sink("~/Desktop/Antonio_qual.txt")
for(cohort in cohorts){
  phylo_use <- eval(parse(text = paste0(cohort, "_Antoniophyloseq_tree")))
  p1 <- phylo_use %>%
    tax_transform("identity", rank = "unique") %>%
    dist_calc("unifrac") %>% 
    ord_calc("NMDS")%>%
    ord_plot(axes = c(1, 2), color = "histology", size = 2.5) +
    scale_color_manual(values=c("yellowgreen", "tomato3"))
  
  dist_obj <- phylo_use %>%
    tax_transform("identity", rank = "unique") %>%
    dist_calc("unifrac") 
 
  mul_var_res <- mul_var_analysis(cohort = cohort, dist_mat = dist_obj@dist,phylo_to_use = phylo_use)
  assign(paste0(cohort, "_Antonio_beta"), p1,.GlobalEnv)
  print(mul_var_res)
}
sink()
  
## Tsementzi beta diversity - use bray and jaccard using NMDS ordination 
sink("~/Desktop/Tsementzi_qual.txt")
for(cohort in cohorts){
  phylo_use <- eval(parse(text = paste0(cohort, "_Tsementziphyloseq_tree")))
  phylo_use <- prune_samples(sample_sums(phylo_use)> 0, phylo_use)
  p1 <- phylo_use %>%
    tax_transform("identity", rank = "unique") %>%
    dist_calc("bray") %>% 
    ord_calc("NMDS")%>%
    ord_plot(axes = c(1, 2), color = "histology", size = 2.5) +
    scale_color_manual(values=c("yellowgreen", "tomato3"))
  
  dist_obj <- phylo_use %>%
    tax_transform("identity", rank = "unique") %>%
    dist_calc("bray") 
  
  mul_var_res <- mul_var_analysis(cohort = cohort, dist_mat = dist_obj@dist, phylo_to_use = phylo_use)
  assign(paste0(cohort, "_Tsementzi_beta"), p1,.GlobalEnv)
  print(mul_var_res)
}
sink()

## Walsh beta diversity - use 
sink("~/Desktop/Walsh_qual.txt")
for(cohort in cohorts){
  phylo_use <- eval(parse(text = paste0(cohort, "_Antoniophyloseq_tree")))
  phylo_use <- prune_samples(sample_sums(phylo_use)> 0, phylo_use)
  p1 <- phylo_use %>%
    tax_transform("identity", rank = "unique") %>%
    dist_calc("wunifrac") %>% 
    ord_calc("PCoA")%>%
    ord_plot(axes = c(1, 2), color = "histology", size = 2.5) +
    scale_color_manual(values=c("yellowgreen", "tomato3"))
  
  dist_obj <- phylo_use %>%
    tax_transform("identity", rank = "unique") %>%
    dist_calc("wunifrac")
  
  mul_var_res <- mul_var_analysis(cohort = cohort, dist_mat = dist_obj@dist, phylo_to_use = phylo_use)
  assign(paste0(cohort, "_Walsh_beta"), p1,.GlobalEnv)
  print(mul_var_res)
}
sink()

## SOTA pipeline
sink("~/Desktop/SOTA_qual.txt")
for(cohort in cohorts){
  phylo_use <- eval(parse(text = paste0(cohort, "_SOTAphyloseq_tree_raw")))
  phylo_use <- prune_samples(sample_sums(phylo_use)> 0, phylo_use)
  p1 <- phylo_use %>%
    tax_fix() %>%
    tax_transform("clr", rank = "genus") %>%
    ord_calc(method = "PCA") %>%
    ord_plot(axes = c(1, 2), color = "histology", size = 2.5) +
    scale_color_manual(values=c("yellowgreen", "tomato3"))
  
  dist_obj <- phylo_use %>%
    tax_fix() %>%
    tax_transform("identity", rank = "genus") %>%
    dist_calc("aitchison")
  
  mul_var_res <- mul_var_analysis(cohort = cohort, dist_mat = dist_obj@dist, phylo_to_use = phylo_use)
  assign(paste0(cohort, "_SOTA_beta"), p1,.GlobalEnv)
  print(mul_var_res)
}
sink()
chao_pipeline_plots <- (Antonio_Chao_beta  + theme(legend.position = "none") |
                          Chao_Chao_beta + theme(legend.position = "none") | 
                           Gressel_Chao_beta + theme(legend.position = "none") |
                           Tsementzi_Chao_beta  + theme(legend.position = "none") |
                           Walsh_Chao_beta+  theme(legend.position = "none"))

antonio_pipeline_plots <- (Antonio_Antonio_beta  +  ggtitle("Antonio") +theme(legend.position = "none")|
                            Chao_Antonio_beta + ggtitle("Chao") + theme(legend.position = "none")| 
                             Gressel_Antonio_beta +ggtitle("Gressel") + theme(legend.position = "none") |
                             Tsementzi_Antonio_beta + ggtitle("Tsementzi") + theme(legend.position = "none") |
                             Walsh_Antonio_beta + ggtitle("Walsh") + theme(legend.position = "none"))

tsementzi_pipeline_plots <- (Antonio_Tsementzi_beta + theme(legend.position = "none") |
                               Chao_Tsementzi_beta + theme(legend.position = "none")| 
                               Gressel_Tsementzi_beta  + theme(legend.position = "none") |
                               Tsementzi_Tsementzi_beta + theme(legend.position = "none") |
                               Walsh_Tsementzi_beta +
                               theme(legend.title = element_text(size=14),
                                     legend.text = element_text(size=13)) +
                               guides(color = guide_legend(override.aes = list(size = 4))))

walsh_pipeline_plots <- (Antonio_Walsh_beta  + theme(legend.position = "none") |
                           Chao_Walsh_beta  + theme(legend.position = "none")| 
                            Gressel_Walsh_beta  + theme(legend.position = "none") |
                            Tsementzi_Walsh_beta + theme(legend.position = "none") |
                            Walsh_Walsh_beta + theme(legend.position = "none"))

SOTA_pipeline_plots <- (Antonio_SOTA_beta  + theme(legend.position = "none") |
                                Chao_SOTA_beta  + theme(legend.position = "none") | 
                              Gressel_SOTA_beta  + theme(legend.position = "none") |
                              Tsementzi_SOTA_beta + theme(legend.position = "none") |
                              Walsh_SOTA_beta + theme(legend.position = "none"))

p <- (antonio_pipeline_plots / chao_pipeline_plots/ tsementzi_pipeline_plots / walsh_pipeline_plots / SOTA_pipeline_plots)
png(paste0("../vaginalMicrobiome/01-Reproducibility_Replicability/Results/betadiversity.png"), width = 6500, height = 3500, res = 300)
print(p)
dev.off()
