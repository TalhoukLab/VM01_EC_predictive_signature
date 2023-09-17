library(patchwork)
library(microViz)
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


source(file = "../vaginalMicrobiome/01-Reproducibility_Replicability/RScripts/00-DataPrep/Antonio_dataPrep.R")
source(file = "../vaginalMicrobiome/01-Reproducibility_Replicability/RScripts/00-DataPrep/Chao_dataPrep.R")
source(file = "../vaginalMicrobiome/01-Reproducibility_Replicability/RScripts/00-DataPrep/Gressel_dataPrep.R")
source(file = "../vaginalMicrobiome/01-Reproducibility_Replicability/RScripts/00-DataPrep/Tsementzi_dataPrep.R")
source(file = "../vaginalMicrobiome/01-Reproducibility_Replicability/RScripts/00-DataPrep/SOTA_dataPrep.R")

## Chao beta diversity - only use unweighted uniFrac distance with an PCA plot. 
cohorts <- c("Antonio","Chao", "Gressel", "Tsementzi", "Walsh")

mul_var_analysis <- function(cohort, dist_mat, phylo_to_use){
  print(phylo_to_use)
  if(cohort == "Chao"){
    test_mul <- adonis2(dist_mat ~ sample_data(phylo_to_use)$histology + 
                          sample_data(phylo_to_use)$age, 
                        permutations = 100, na.action = "na.omit", by = "margin")
  } 
  if(cohort == "Antonio"){
    test_mul <- adonis2(dist_mat ~ sample_data(phylo_to_use)$histology + 
                          sample_data(phylo_to_use)$age + 
                          sample_data(phylo_to_use)$BMI + 
                          sample_data(phylo_to_use)$pHRecoded,
                        permutations = 100, na.action = "na.omit", by = "margin")
  } 
  if(cohort == "Tsementzi"){
    test_mul <- adonis2(dist_mat ~ sample_data(phylo_to_use)$histology + 
                          sample_data(phylo_to_use)$age + 
                          sample_data(phylo_to_use)$BMI + 
                          as.factor(sample_data(phylo_to_use)$pHRecoded) + 
                          as.factor(sample_data(phylo_to_use)$ethnicity),
                        permutations = 100, na.action = "na.omit", by = "margin")
  } 
  if(cohort == 'Walsh'){
    test_mul <- adonis2(dist_mat ~ sample_data(phylo_to_use)$histology + 
                          sample_data(phylo_to_use)$age + 
                          sample_data(phylo_to_use)$BMI + 
                          as.factor(sample_data(phylo_to_use)$pHRecoded) + 
                          as.factor(sample_data(phylo_to_use)$ethnicityRecoded),
                        permutations = 100, na.action = "na.omit", by = "margin")
  }
  if(cohort == "Gressel"){
    test_mul <- adonis2(dist_mat ~ sample_data(phylo_to_use)$histology,
                        permutations = 100, na.action = "na.omit", by = "margin")
  }
  return(test_mul)
}


sink("../vaginalMicrobiome/01-Reproducibility_Replicability/Results/02-BetaDiversity/Chao_qual_margin.txt")
for(cohort in cohorts){
  print(cohort)
  phylo_use <- eval(parse(text = paste0(cohort, "_Chaophyloseq_tree_raw")))
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
sink("../vaginalMicrobiome/01-Reproducibility_Replicability/Results/02-BetaDiversity/Antonio_qual.txt")
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
sink("../vaginalMicrobiome/01-Reproducibility_Replicability/Results/02-BetaDiversity/Tsementzi_qual.txt")
for(cohort in cohorts){
  phylo_use <- eval(parse(text = paste0(cohort, "_Tsementziphyloseq_tree_raw")))
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
sink("../vaginalMicrobiome/01-Reproducibility_Replicability/Results/02-BetaDiversity/Walsh_qual.txt")
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
sink("../vaginalMicrobiome/01-Reproducibility_Replicability/Results/02-BetaDiversity/SOTA_qual_margin.txt")
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
                          Walsh_Chao_beta+  theme(legend.position = "none") |
                          Tsementzi_Chao_beta  + theme(legend.position = "none") |
                           Gressel_Chao_beta + theme(legend.position = "none") |
                          Chao_Chao_beta + theme(legend.position = "none"))

antonio_pipeline_plots <- (Antonio_Antonio_beta  +  ggtitle("Antonio") +theme(legend.position = "none")|
                             Walsh_Antonio_beta + ggtitle("Walsh") + theme(legend.position = "none") |
                             Tsementzi_Antonio_beta + ggtitle("Tsementzi") + theme(legend.position = "none") |
                             Gressel_Antonio_beta +ggtitle("Gressel") + theme(legend.position = "none") |
                             Chao_Antonio_beta + ggtitle("Chao") + theme(legend.position = "none"))

tsementzi_pipeline_plots <- (Antonio_Tsementzi_beta + theme(legend.position = "none") |
                               Walsh_Tsementzi_beta + theme(legend.position = "none")| 
                               Tsementzi_Tsementzi_beta  + theme(legend.position = "none") |
                               Gressel_Tsementzi_beta + theme(legend.position = "none") |
                               Chao_Tsementzi_beta +
                               theme(legend.title = element_text(size=14),
                                     legend.text = element_text(size=13)) +
                               guides(color = guide_legend(override.aes = list(size = 4))))

walsh_pipeline_plots <- (Antonio_Walsh_beta  + theme(legend.position = "none") |
                           Walsh_Walsh_beta  + theme(legend.position = "none")| 
                           Tsementzi_Walsh_beta + theme(legend.position = "none") |
                           Gressel_Walsh_beta + theme(legend.position = "none") |
                           Chao_Walsh_beta + theme(legend.position = "none"))

SOTA_pipeline_plots <- (Antonio_SOTA_beta  + theme(legend.position = "none", plot.background = element_rect(fill = rgb(0.96, 0.96, 0.96, alpha = 0.6), 
                                                                                                            color = rgb(0.96, 0.96, 0.96, alpha = 0.6))) | 
                          Walsh_SOTA_beta + theme(legend.position = "none",  plot.background = element_rect(fill = rgb(0.96, 0.96, 0.96, alpha = 0.6), 
                                                                                                            color = rgb(0.96, 0.96, 0.96, alpha = 0.6))) | 
                          Tsementzi_SOTA_beta + theme(legend.position = "none",  plot.background = element_rect(fill = rgb(0.96, 0.96, 0.96, alpha = 0.6), 
                                                                                                                color = rgb(0.96, 0.96, 0.96, alpha = 0.6))) |
                          Gressel_SOTA_beta + theme(legend.position = "none",  plot.background = element_rect(fill = rgb(0.96, 0.96, 0.96, alpha = 0.6), 
                                                                                                              color = rgb(0.96, 0.96, 0.96, alpha = 0.6))) |
                          Chao_SOTA_beta + theme(legend.position = "none",  plot.background = element_rect(fill = rgb(0.96, 0.96, 0.96, alpha = 0.6), 
                                                                                                           color = rgb(0.96, 0.96, 0.96, alpha = 0.6))))

p <- (antonio_pipeline_plots / walsh_pipeline_plots/ tsementzi_pipeline_plots / chao_pipeline_plots / SOTA_pipeline_plots)
png(paste0("../vaginalMicrobiome/01-Reproducibility_Replicability/Results/betadiversity.png"), width = 6500, height = 3500, res = 300)
print(p)
dev.off()

library(ggplot2)
library(latex2exp)


pdf("../vaginalMicrobiome/01-Reproducibility_Replicability/Results/02-BetaDiversity/margin.pdf",  width=16, height=5)


df <- read.csv("../vaginalMicrobiome/01-Reproducibility_Replicability/Results/02-BetaDiversity/betaDiversity_R2_margin.csv", header = TRUE, sep = ",")
df$pipeline <- paste0(df$pipeline, "_pipeline")
df$R2 <- round(df$R2, 2)
df$covariate <- factor(df$covariate, levels = c("histology", "BMI", "pH", "age", "ethnicity"))
df$cohort <- factor(df$cohort, levels = c("Antonio", "Walsh", "Tsementzi", "Gressel", "Chao"))
df$pipeline <- factor(df$pipeline, levels = c("Antonio_pipeline", "Walsh_pipeline", "Tsementzi_pipeline", "Gressel_pipeline", "Chao_pipeline", "SOTA_pipeline"))
ggplot(df, aes(cohort, covariate, fill= R2)) +  geom_tile(aes(fill = R2)) + 
  geom_tile(data = df, fill="transparent", color = ifelse(df$sig=="Sig", 'black', NA), size = 0.3) +
  geom_text(aes(label = R2), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "red", name = unname(TeX(c("$R^2$")))) + facet_wrap(~pipeline,  ncol=5) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title=element_text(size=16),
        strip.text.x = element_text(size = 12),
        legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))
dev.off()

