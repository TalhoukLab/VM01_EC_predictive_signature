library(dplyr)
all_cohorts <- c("Antonio", "Walsh", "Tsementzi", "Gressel", "Chao")
all_levels <- c("phylum", "class", "order", "family", "genus")
source(file = "../VM01_reproducibility_replicability/RScripts/01-Analysis/03d-SOTA_biomarker.R")

## Unable to reproduce or replicate results from Antonio and Walsh et al. and Gressel et al. 

## Tsementzi and Chao pipelines 
pipeline <- "Chao_pipeline"
for(cohort in all_cohorts){
  for(level in all_levels){
    full_results <- read_delim(paste0("../vaginalMicrobiome/01-Reproducibility_Replicability/Results/03-DET/Chao_pipeline/results/", level, "_Chaopipeline_", cohort, ".res"),
                               col_names = F)
    full_results <- na.omit(full_results)
    full_results <- full_results %>%
                        dplyr::select(-c(X2, X4)) %>%
                        dplyr::filter(!(X1 %in% c("pHRecoded", "age", 
                                                  "BMI", "ethnicityRecoded", "menopausal_status"))) %>%
                        dplyr::filter(X5 <= 0.05)
    
    colnames(full_results) <- c("taxa", "class", "p-value")
    full_results$cohort <- cohort 
    full_results$expression <- with(full_results, 
                                ifelse(class=="EC", 'High in EC', 'Low in EC'))
    full_results$pipeline <- pipeline
    full_results <- full_results %>% 
                  dplyr::select(c(taxa, expression, cohort, pipeline))
    assign(paste0(cohort, "_", level, "_Chao_df"), full_results,.GlobalEnv)
    
  }
}


for(level in all_levels){
  Tsementzi_rbind <- rbind(eval(parse(text = paste0("Antonio_", level, "_Tsementzi_df"))),
                           eval(parse(text = paste0("Walsh_", level, "_Tsementzi_df"))),
                           eval(parse(text = paste0("Tsementzi_", level, "_Tsementzi_df"))),
                           eval(parse(text = paste0("Gressel_", level, "_Tsementzi_df"))),
                           eval(parse(text = paste0("Chao_", level, "_Tsementzi_df"))))
  
  
  Chao_rbind <- rbind(eval(parse(text = paste0("Antonio_", level, "_Chao_df"))),
                         eval(parse(text = paste0("Walsh_", level, "_Chao_df"))),
                         eval(parse(text = paste0("Tsementzi_", level, "_Chao_df"))),
                         eval(parse(text = paste0("Gressel_", level, "_Chao_df"))),
                         eval(parse(text = paste0("Chao_", level, "_Chao_df"))))
  
  SOTA_rbind <- rbind(eval(parse(text = paste0("Antonio_", level, "_deseq2_df"))),
                      eval(parse(text = paste0("Walsh_", level, "_deseq2_df"))),
                      eval(parse(text = paste0("Tsementzi_", level, "_deseq2_df"))),
                      eval(parse(text = paste0("Gressel_", level, "_deseq2_df"))),
                      eval(parse(text = paste0("Chao_", level, "_deseq2_df"))))
  all <- rbind(Tsementzi_rbind,
               Chao_rbind, SOTA_rbind)
  all$taxa <- tolower(all$taxa)
  all$cohort <- factor(all$cohort, levels = c("Antonio", "Walsh", "Tsementzi", "Gressel", "Chao"))
  all$pipeline <- factor(all$pipeline, levels = c("Tsementzi_pipeline", "Gressel_pipeline", "Chao_pipeline", "SOTA_pipeline"))
  
  all$taxa <- factor(all$taxa)
  all$expression <- as.factor(all$expression)
  
  library(forcats)
  all <- all %>%
    arrange(cohort, pipeline) %>%   # rearrange the df in the order we want
    mutate(taxa = factor(taxa, unique(taxa))) 
  a <-ggplot(all, aes(cohort, taxa, fill= expression)) +
    geom_tile(color = "white",
              lwd = 1,
              linetype = 1, alpha = 0.9) + 
    geom_rect(data = subset(all,pipeline == 'SOTA_pipeline'),aes(fill = pipeline),
              xmin = -Inf,xmax = Inf,
              ymin = -Inf,ymax = Inf,alpha = 0.09)+  geom_tile(color = "white",
                                                               lwd = 1,
                                                               linetype = 1, alpha = 0.9)+
    scale_fill_manual(values=c("tomato3", "yellowgreen", "whitesmoke")) + 
    ggtitle(paste0("Diffrentially expressed taxa at ",  level, " level")) + theme_bw() + facet_wrap(~pipeline)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 12),
          axis.title=element_text(size=16),
          strip.text.x = element_text(size = 12),
          legend.key.size = unit(1, 'cm'),
          legend.key.height = unit(1, 'cm'),
          legend.key.width = unit(1, 'cm'),
          legend.title = element_text(size=10),
          legend.text = element_text(size=10), strip.background =element_rect(fill="white")) 
  pdf(paste0("../vaginalMicrobiome/01-Reproducibility_Replicability/Results/03-DET/", level, ".pdf"),  width=10, height=12)
  
  a
  dev.off()
}
