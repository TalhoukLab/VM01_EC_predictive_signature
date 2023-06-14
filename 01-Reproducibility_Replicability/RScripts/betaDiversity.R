## This script is for beta diversity.
library(patchwork)
## Angel beta diversity - only use unweighted uniFrac distance with an PCA plot. 
cohorts <- c("Angel", "Antonio", "Gressel", "Tsementzi", "Walsh")
for(cohort in cohorts){
  phylo_use <- eval(parse(text = paste0(cohort, "_Angelphyloseq_tree")))
  p1 <- phylo_use %>%
    tax_transform("identity", rank = "unique") %>%
    dist_calc("unifrac") %>% 
    ord_calc("PCA")%>%
    ord_plot(axes = c(1, 2), color = "histology", size = 2.5) 
  assign(paste0(cohort, "_Angel_beta"), p1,.GlobalEnv)
}


## Antonio beta diversity - only use unweighted unifrac using NMDS ordination 
for(cohort in cohorts){
  phylo_use <- eval(parse(text = paste0(cohort, "_Antoniophyloseq_tree")))
  p1 <- phylo_use %>%
    tax_transform("identity", rank = "unique") %>%
    dist_calc("unifrac") %>% 
    ord_calc("NMDS")%>%
    ord_plot(axes = c(1, 2), color = "histology", size = 2.5) 
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
    ord_plot(axes = c(1, 2), color = "histology", size = 2.5) 
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
    ord_plot(axes = c(1, 2), color = "histology", size = 2.5) 
  assign(paste0(cohort, "_Walsh_beta"), p1,.GlobalEnv)
}


angel_pipeline_plots <- (Angel_Angel_beta+ggtitle("Angel") + theme(legend.position = "none")  + scale_color_manual(values=c("khaki", "yellowgreen", "tomato3"))| 
                           Antonio_Angel_beta+ggtitle("Antonio")  + theme(legend.position = "none") + scale_color_manual(values=c("khaki", "yellowgreen", "tomato3"))|
                           Gressel_Angel_beta+ggtitle("Gressel") + theme(legend.position = "none")  + scale_color_manual(values=c("yellowgreen","tomato3")) |
                           Tsementzi_Angel_beta + ggtitle("Tsementzi") + theme(legend.position = "none")+ scale_color_manual(values=c("yellowgreen","tomato3")) |
                           Walsh_Angel_beta+ggtitle("Walsh") +  theme(legend.position = "none") +scale_color_manual(values=c("khaki", "yellowgreen", "tomato3")))

antonio_pipeline_plots <- (Angel_Antonio_beta + theme(legend.position = "none")  + scale_color_manual(values=c("yellowgreen", "khaki", "tomato3"))| 
                             Antonio_Antonio_beta  + theme(legend.position = "none") + scale_color_manual(values=c("yellowgreen", "khaki", "tomato3"))|
                             Gressel_Antonio_beta + theme(legend.position = "none")  + scale_color_manual(values=c("yellowgreen","tomato3")) |
                             Tsementzi_Antonio_beta + theme(legend.position = "none")  + scale_color_manual(values=c("yellowgreen","tomato3")) |
                             Walsh_Antonio_beta + theme(legend.position = "none") + scale_color_manual(values=c("yellowgreen", "khaki", "tomato3")))

tsementzi_pipeline_plots <- (Angel_Tsementzi_beta + theme(legend.position = "none")  + scale_color_manual(values=c("yellowgreen", "khaki", "tomato3"))| 
                               Antonio_Tsementzi_beta + theme(legend.position = "none")  + scale_color_manual(values=c("yellowgreen", "khaki", "tomato3"))|
                               Gressel_Tsementzi_beta  + theme(legend.position = "none") + scale_color_manual(values=c("yellowgreen","tomato3")) |
                               Tsementzi_Tsementzi_beta + theme(legend.position = "none") + scale_color_manual(values=c("yellowgreen","tomato3")) |
                               Walsh_Tsementzi_beta + scale_color_manual(values=c("yellowgreen", "khaki", "tomato3")) +
                               theme(legend.title = element_text(size=14),
                                     legend.text = element_text(size=13)) +
                               guides(color = guide_legend(override.aes = list(size = 4))))

walsh_pipeline_plots <- (Angel_Walsh_beta  + theme(legend.position = "none") + scale_color_manual(values=c("yellowgreen", "khaki", "tomato3"))| 
                            Antonio_Walsh_beta  + theme(legend.position = "none") + scale_color_manual(values=c("yellowgreen", "khaki", "tomato3"))|
                            Gressel_Walsh_beta  + theme(legend.position = "none") + scale_color_manual(values=c("yellowgreen","tomato3")) |
                            Tsementzi_Walsh_beta + theme(legend.position = "none") + scale_color_manual(values=c("yellowgreen","tomato3")) |
                            Walsh_Walsh_beta + theme(legend.position = "none") + scale_color_manual(values=c("yellowgreen", "khaki", "tomato3")))

inhouse_pipeline_plots <- (Angel_InHouse_beta  + theme(legend.position = "none") + scale_color_manual(values=c("khaki", "yellowgreen", "tomato3"))| 
                              Antonio_InHouse_beta  + theme(legend.position = "none") + scale_color_manual(values=c("khaki", "yellowgreen", "tomato3"))|
                              Gressel_InHouse_beta  + theme(legend.position = "none") + scale_color_manual(values=c("yellowgreen","tomato3")) |
                              Tsementzi_InHouse_beta + theme(legend.position = "none") + scale_color_manual(values=c("yellowgreen","tomato3")) |
                              Walsh_InHouse_beta + theme(legend.position = "none") +scale_color_manual(values=c("khaki", "yellowgreen", "tomato3")))

p <- (angel_pipeline_plots / antonio_pipeline_plots / tsementzi_pipeline_plots / walsh_pipeline_plots / inhouse_pipeline_plots)
png(paste0("~/Desktop/betadiversity.png"), width = 1500, height = 900)
print(p)
dev.off()
