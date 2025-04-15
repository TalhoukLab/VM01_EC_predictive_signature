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


source(file = "../VM01_reproducibility_replicability/RScripts/00-DataPrep/Antonio_dataPrep.R")
source(file = "../VM01_reproducibility_replicability/RScripts/00-DataPrep/Chao_dataPrep.R")
source(file = "../VM01_reproducibility_replicability/RScripts/00-DataPrep/Gressel_dataPrep.R")
source(file = "../VM01_reproducibility_replicability/RScripts/00-DataPrep/Tsementzi_dataPrep.R")
source(file = "../VM01_reproducibility_replicability/RScripts/00-DataPrep/dada2_dataPrep.R")

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
                          sample_data(phylo_to_use)$BMI ,
                        permutations = 100, na.action = "na.omit", by = "margin")
  } 
  if(cohort == "Tsementzi"){
    test_mul <- adonis2(dist_mat ~ sample_data(phylo_to_use)$histology + 
                          sample_data(phylo_to_use)$age + 
                          sample_data(phylo_to_use)$bmi + 
                          as.factor(sample_data(phylo_to_use)$ethnicityRecode),
                        permutations = 100, na.action = "na.omit", by = "margin")
  } 
  if(cohort == 'Walsh'){
    test_mul <- adonis2(dist_mat ~ sample_data(phylo_to_use)$histology + 
                          sample_data(phylo_to_use)$age + 
                          as.numeric(sample_data(phylo_to_use)$bmi) + 
                          as.factor(sample_data(phylo_to_use)$ethnicityRecode), 
                        permutations = 100, na.action = "na.omit", by = "margin")
  }
  if(cohort == "Gressel"){
    test_mul <- adonis2(dist_mat ~ sample_data(phylo_to_use)$histology,
                        permutations = 100, na.action = "na.omit", by = "margin")
  }
  return(test_mul)
}


sink("~/Desktop/Chao_qual_margin.txt")
for(cohort in cohorts){
  print(cohort)
  phylo_use <- eval(parse(text = paste0(cohort, "_Chaophyloseq_tree_raw")))
  
  dist = phyloseq::distance(phylo_use,
                            method = "unifrac")
  ordi = ordinate(phylo_use, distance="unifrac", method="PCoA")
  p1 <- plot_ordination(phylo_use, ordi, "samples", color="histology") +
    geom_point(size = 2.5) +
    scale_color_manual(values=c("yellowgreen", "tomato3"))
  mul_var_res <- mul_var_analysis(cohort = cohort, dist_mat = dist, phylo_to_use = phylo_use)
  assign(paste0(cohort, "_Chao_beta"), p1,.GlobalEnv)
  print(mul_var_res)
}

sink()

## Antonio beta diversity - only use unweighted unifrac using NMDS ordination 
sink("~/Desktop/Antonio_qual_margin.txt")
for(cohort in cohorts){
  phylo_use <- eval(parse(text = paste0(cohort, "_Antoniophyloseq_tree")))
  dist = phyloseq::distance(phylo_use,
                            method = "unifrac")
  ordi = ordinate(phylo_use, distance="unifrac", method="NMDS")
  p1 <- plot_ordination(phylo_use, ordi, "samples", color="histology") +
    geom_point(size = 2.5) +
    scale_color_manual(values=c("yellowgreen", "tomato3"))
  mul_var_res <- mul_var_analysis(cohort = cohort, dist_mat = dist, phylo_to_use = phylo_use)
  assign(paste0(cohort, "_Antonio_beta"), p1,.GlobalEnv)
  print(mul_var_res)
}
sink()
  
## Tsementzi beta diversity - use bray and jaccard using NMDS ordination 
sink("~/Desktop/Tsementzi_qual_margin.txt")
for(cohort in cohorts){
  phylo_use <- eval(parse(text = paste0(cohort, "_Tsementziphyloseq_tree")))
  phylo_use <- prune_samples(sample_sums(phylo_use)> 0, phylo_use)
  dist = phyloseq::distance(phylo_use,
                            method = "bray")
  ordi = ordinate(phylo_use, distance="bray", method="NMDS")
  p1 <- plot_ordination(phylo_use, ordi, "samples", color="histology") +
    geom_point(size = 2.5) +
    scale_color_manual(values=c("yellowgreen", "tomato3"))
  mul_var_res <- mul_var_analysis(cohort = cohort, dist_mat = dist, phylo_to_use = phylo_use)
  assign(paste0(cohort, "_Tsementzi_beta"), p1,.GlobalEnv)
  print(mul_var_res)
}
sink()

## Walsh beta diversity - use 
sink("~/Desktop/Walsh_qual_margin.txt")
for(cohort in cohorts){
  phylo_use <- eval(parse(text = paste0(cohort, "_Antoniophyloseq_tree")))
  phylo_use <- prune_samples(sample_sums(phylo_use)> 0, phylo_use)
  dist = phyloseq::distance(phylo_use,
                            method = "wunifrac")
  ordi = ordinate(phylo_use, distance="wunifrac", method="PCoA")
  p1 <- plot_ordination(phylo_use, ordi, "samples", color="histology") +
    geom_point(size = 2.5) +
    scale_color_manual(values=c("yellowgreen", "tomato3"))
  
  mul_var_res <- mul_var_analysis(cohort = cohort, dist_mat = dist, phylo_to_use = phylo_use)
  assign(paste0(cohort, "_Walsh_beta"), p1,.GlobalEnv)
  print(mul_var_res)
}
sink()

## SOTA pipeline
library(phyloseq)
library(ggplot2)
library(microViz)
sink("~/Desktop/DADA2_qual_margin.txt")
for(cohort in cohorts){
  phylo_use <- eval(parse(text = paste0(cohort, "_dada2phyloseq_tree")))
  phylo_use <- microbiome::transform(phylo_use, transform = "log")
  dist = phyloseq::distance(phylo_use,
                     method = "wunifrac")
  ordi = ordinate(phylo_use, distance="wunifrac", method="PCoA")
  p1 <- plot_ordination(phylo_use, ordi, "samples", color="histology") +
    geom_point(size = 2.5) +
    scale_color_manual(values=c("yellowgreen", "tomato3"))
  mul_var_res <- mul_var_analysis(cohort = cohort, dist_mat = dist, phylo_to_use = phylo_use)
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

SOTA_pipeline_plots <- (Antonio_SOTA_beta  + theme(legend.position = "none", 
                                                   plot.background = element_rect(fill = rgb(0.96, 0.96, 0.96, 
                                                                                             alpha = 0.6), 
                                                                                  color = rgb(0.96, 0.96, 0.96, 
                                                                                              alpha = 0.6))) | 
                          Walsh_SOTA_beta + theme(legend.position = "none", 
                                                  plot.background = element_rect(fill = rgb(0.96, 0.96, 0.96, alpha = 0.6), 
                                                                                                            color = rgb(0.96, 0.96, 0.96, 
                                                                                                                        alpha = 0.6))) | 
                          Tsementzi_SOTA_beta + theme(legend.position = "none",  plot.background = element_rect(fill = rgb(0.96, 0.96, 0.96, alpha = 0.6), 
                                                                                                                color = rgb(0.96, 0.96, 0.96, alpha = 0.6))) |
                          Gressel_SOTA_beta + theme(legend.position = "none",  plot.background = element_rect(fill = rgb(0.96, 0.96, 0.96, alpha = 0.6), 
                                                                                                              color = rgb(0.96, 0.96, 0.96, alpha = 0.6))) |
                          Chao_SOTA_beta + theme(legend.position = "none",  plot.background = element_rect(fill = rgb(0.96, 0.96, 0.96, alpha = 0.6), 
                                                                                                           color = rgb(0.96, 0.96, 0.96, alpha = 0.6))))

p <- (antonio_pipeline_plots / walsh_pipeline_plots/ tsementzi_pipeline_plots / chao_pipeline_plots / SOTA_pipeline_plots)
png(paste0("~/Desktop/betadiversity.png"), width = 6500, height = 3500, res = 300)
print(p)
dev.off()

library(ggplot2)
library(latex2exp)


pdf("~/Desktop/margin1.pdf",  width=25, height=8)


df <- read.csv("../VM01_EC_predictive_signature/01-reproducibility_replicability/results/02-betadiversity/betaDiversity_R2_margin.csv", header = TRUE, sep = ",")
df$pipeline <- paste0(df$pipeline, "_pipeline")
df$R2 <- as.numeric(df$R2)
df$p.value <- as.numeric(df$pvalue)
df$R2 <- round(df$R2, 2)
df$covariate <- factor(df$covariate, labels = c("Health condition","BMI", "pH", "Age", "Ethnicity"), levels = c("histology", "BMI", "pH", "age", "ethnicity"))
df$cohort <- factor(df$cohort, levels = c("Antonio", "Walsh", "Tsementzi", "Gressel", "Chao"))
df$pipeline <- factor(df$pipeline, levels = c("Antonio_pipeline", "Walsh_pipeline", "Tsementzi_pipeline", "Gressel_pipeline", "Chao_pipeline", "DADA2_pipeline"),
                      labels = c(unname(TeX(c("Antonio_pipeline\nUniFrac"))), unname(TeX(c("Walsh_pipeline\nwUniFrac"))), 
                                 unname(TeX(c("Tsementzi_pipeline\nBray-Curtis"))), unname(TeX(c("Gressel_pipeline"))), 
                                 unname(TeX(c("Chao_pipeline\nUniFrac"))), unname(TeX(c("DADA2_pipeline\nwUniFrac")))))
df_dada <- df[df$pipeline=="DADA2_pipeline\nwUniFrac", ]
beta <- ggplot(df, aes(cohort, covariate, fill= R2)) +  geom_tile(aes(fill = R2)) + facet_grid(. ~ pipeline) + 
  geom_tile(data = df, fill="transparent", color = ifelse(df$sig=="Sig", 'black', NA), size = 0.5) +
  geom_text(aes(label = R2), color = "black", size = 6.0) +
  scale_fill_gradient(low = "white", high = "red", name = unname(TeX(c("$R^2$"))))  + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 28),
        axis.text.y = element_text(size = 28),
        axis.title=element_text(size=29),
        strip.text.x = element_text(size = 28),
        legend.key.size = unit(2, 'cm'),
        legend.key.height = unit(2, 'cm'),
        legend.key.width = unit(2, 'cm'),
        legend.position = "bottom",
        legend.title = element_text(size=29),
        legend.text = element_text(size=25), plot.title = element_text(size = 29, hjust = 0.5, face = "bold")) + ggtitle("Beta diversity")
#print(beta)

pdf("~/Desktop/figure1.pdf",  width=30, height=19)

pw <- (bp1 / beta )
pw + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 23, face = "bold"))
dev.off()

