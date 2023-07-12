## This script is for beta diversity.
library(patchwork)
library(microViz)
## Chao beta diversity - only use unweighted uniFrac distance with an PCA plot. 
cohorts <- c("Antonio","Chao", "Gressel", "Tsementzi", "Walsh")
for(cohort in cohorts){
  phylo_use <- eval(parse(text = paste0(cohort, "_Chaophyloseq_tree")))
  p1 <- phylo_use %>%
    tax_transform("identity", rank = "unique") %>%
    dist_calc("unifrac") %>% 
    ord_calc("PCA")%>%
    ord_plot(axes = c(1, 2), color = "histology", size = 2.5) +
    scale_color_manual(values=c("yellowgreen", "tomato3"))
  assign(paste0(cohort, "_Chao_beta"), p1,.GlobalEnv)
}


## Antonio beta diversity - only use unweighted unifrac using NMDS ordination 
for(cohort in cohorts){
  phylo_use <- eval(parse(text = paste0(cohort, "_Antoniophyloseq_tree")))
  p1 <- phylo_use %>%
    tax_transform("identity", rank = "unique") %>%
    dist_calc("unifrac") %>% 
    ord_calc("NMDS")%>%
    ord_plot(axes = c(1, 2), color = "histology", size = 2.5) +
    scale_color_manual(values=c("yellowgreen", "tomato3"))
  assign(paste0(cohort, "_Antonio_beta"), p1,.GlobalEnv)
}
  
## Tsementzi beta diversity - use bray and jaccard using NMDS ordination 
for(cohort in cohorts){
  phylo_use <- eval(parse(text = paste0(cohort, "_Tsementziphyloseq_tree")))
  phylo_use <- prune_samples(sample_sums(phylo_use)> 0, phylo_use)
  p1 <- phylo_use %>%
    tax_transform("identity", rank = "unique") %>%
    dist_calc("bray") %>% 
    ord_calc("NMDS")%>%
    ord_plot(axes = c(1, 2), color = "histology", size = 2.5) +
    scale_color_manual(values=c("yellowgreen", "tomato3"))
  assign(paste0(cohort, "_Tsementzi_beta"), p1,.GlobalEnv)
}


## Walsh beta diversity - use 
for(cohort in cohorts){
  phylo_use <- eval(parse(text = paste0(cohort, "_Antoniophyloseq_tree")))
  phylo_use <- prune_samples(sample_sums(phylo_use)> 0, phylo_use)
  p1 <- phylo_use %>%
    tax_transform("identity", rank = "unique") %>%
    dist_calc("wunifrac") %>% 
    ord_calc("PCoA")%>%
    ord_plot(axes = c(1, 2), color = "histology", size = 2.5) +
    scale_color_manual(values=c("yellowgreen", "tomato3"))
  assign(paste0(cohort, "_Walsh_beta"), p1,.GlobalEnv)
}

## SOTA pipeline
for(cohort in cohorts){
  phylo_use <- eval(parse(text = paste0(cohort, "_SOTAphyloseq_tree_raw")))
  phylo_use <- prune_samples(sample_sums(phylo_use)> 0, phylo_use)
  p1 <- phylo_use %>%
    tax_fix() %>%
    tax_transform("clr", rank = "genus") %>%
    ord_calc(method = "PCA") %>%
    ord_plot(axes = c(1, 2), color = "histology", size = 2.5) +
    scale_color_manual(values=c("yellowgreen", "tomato3"))
  assign(paste0(cohort, "_SOTA_beta"), p1,.GlobalEnv)
}

chao_pipeline_plots <- (Chao_Chao_beta+ggtitle("Chao") + theme(legend.position = "none") | 
                           Antonio_Chao_beta+ggtitle("Antonio")  + theme(legend.position = "none") |
                           Gressel_Chao_beta+ggtitle("Gressel") + theme(legend.position = "none") |
                           Tsementzi_Chao_beta + ggtitle("Tsementzi") + theme(legend.position = "none") |
                           Walsh_Chao_beta+ggtitle("Walsh") +  theme(legend.position = "none"))

antonio_pipeline_plots <- (Chao_Antonio_beta + theme(legend.position = "none")| 
                             Antonio_Antonio_beta  + theme(legend.position = "none")|
                             Gressel_Antonio_beta + theme(legend.position = "none") |
                             Tsementzi_Antonio_beta + theme(legend.position = "none") |
                             Walsh_Antonio_beta + theme(legend.position = "none"))

tsementzi_pipeline_plots <- (Chao_Tsementzi_beta + theme(legend.position = "none")| 
                               Antonio_Tsementzi_beta + theme(legend.position = "none") |
                               Gressel_Tsementzi_beta  + theme(legend.position = "none") |
                               Tsementzi_Tsementzi_beta + theme(legend.position = "none") |
                               Walsh_Tsementzi_beta +
                               theme(legend.title = element_text(size=14),
                                     legend.text = element_text(size=13)) +
                               guides(color = guide_legend(override.aes = list(size = 4))))

walsh_pipeline_plots <- (Chao_Walsh_beta  + theme(legend.position = "none")| 
                            Antonio_Walsh_beta  + theme(legend.position = "none") |
                            Gressel_Walsh_beta  + theme(legend.position = "none") |
                            Tsementzi_Walsh_beta + theme(legend.position = "none") |
                            Walsh_Walsh_beta + theme(legend.position = "none"))

SOTA_pipeline_plots <- (Chao_SOTA_beta  + theme(legend.position = "none") | 
                              Antonio_SOTA_beta  + theme(legend.position = "none") |
                              Gressel_SOTA_beta  + theme(legend.position = "none") |
                              Tsementzi_SOTA_beta + theme(legend.position = "none") |
                              Walsh_SOTA_beta + theme(legend.position = "none"))

p <- (antonio_pipeline_plots / chao_pipeline_plots/ tsementzi_pipeline_plots / walsh_pipeline_plots / SOTA_pipeline_plots)
png(paste0("~/Desktop/betadiversity.png"), width = 1500, height = 900)
print(p)
dev.off()
