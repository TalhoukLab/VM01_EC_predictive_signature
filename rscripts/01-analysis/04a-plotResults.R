library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

all_scores_cross <- read.csv("../VM01_reproducibility_replicability/Results/04-MLRep/crossLab.csv", sep=",", header = TRUE)
all_scores_cross$type <- "cross lab"
all_scores_within <- read.csv("../VM01_reproducibility_replicability/Results/04-MLRep/withinLab.csv", sep=",", header = TRUE)
all_scores_within$type <- "within lab" 
all_scores_within$training_cohort <- all_scores_within$training_cohort %>% str_replace("_.*", "")
all_scores_within$testing_cohort <- all_scores_within$testing_cohort %>% str_replace("_.*", "")
all_scores <- rbind(all_scores_cross, all_scores_within)

all_scores$level = factor(all_scores$level, levels=c("phylum", "class", "order", "family", "genus"), 
                                labels=c("phylum", "class", "order", "family", "genus")) 
all_scores$testing_cohort = factor(all_scores$testing_cohort, levels=c("Antonio","Gressel", "Chao", 
                                                                         "Tsementzi", "Walsh"), 
                                    labels=c("Antonio (N=22)","Gressel (N=29)", "Chao (N=30)", 
                                             "Tsementzi (N=36)", "Walsh (N=149)"))
all_scores$training_cohort <- factor(all_scores$training_cohort, levels=c("Antonio","Gressel", "Chao", 
                                                                          "Tsementzi", "Walsh"), 
                                     labels = c(levels=c("Antonio","Gressel", "Chao", 
                                                         "Tsementzi", "Walsh")))
all_scores$pipeline = factor(all_scores$pipeline, levels=c("Antonio", "Tsementzi", "Gressel", "Chao", "SOTA"),
                             labels=c("Antonio_Walsh_pipeline", "Tsementzi_pipeline", "Gressel_pipeline", "Chao_pipeline", "SOTA_pipeline"))

pdf("~/Desktop/ggplot1.pdf",  width=12, height=9)

all_scores %>%
  filter(pipeline != "SOTA") %>% 
  filter(level == "genus") %>%
  ggplot(aes(x = training_cohort, y=V1, group = pipeline, colour = pipeline, alpha = pipeline)) +
  geom_line(data = filter(all_scores, type =="cross lab" & level == "genus"), size = 2.0) +
  geom_line(data = filter(all_scores, pipeline == "SOTA_pipeline" & type =="cross" & level == "genus"), size = 2.3) +
  geom_point(aes(shape = type, size = type, alpha = pipeline))+
  geom_point(data = filter(all_scores, type == "within lab" & level == "genus"), aes(shape = type, size = type))+
  scale_shape_manual(values=c(16, 17))+
  scale_size_manual(values=c(2.6, 2.8)) +
  scale_alpha_manual(values=c(0.4, 0.4, 0.4, 0.4, 1)) +
  facet_grid(level ~ testing_cohort) +  ylim(0,1) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =18),
        axis.text.y = element_text(size = 18),
        axis.title=element_text(size=21),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        strip.background =element_rect(fill="white"),
        legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(size=21, face = "bold"),
        legend.text = element_text(size=21),
        legend.position="bottom", 
        legend.box="vertical") + ylab("AUC")
dev.off()

