library(dplyr)
library(tidyr)
library(ggplot2)

all_scores_cross <- read.csv("~/Desktop/crossLab.csv", sep=",", header = TRUE)
all_scores_cross$type <- "cross lab"
all_scores_within <- read.csv("~/Desktop/withinLab.csv", sep=",", header = TRUE)
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

pdf("~/Desktop/ggplot.pdf",  width=12, height=9)


all_scores %>%
  filter(pipeline != "SOTA") %>% 
  ggplot(aes(x = training_cohort, y=V1, group = pipeline, colour = pipeline, alpha = pipeline)) +
  geom_line(data = filter(all_scores, type =="cross lab"), size = 1.5) +
  geom_line(data = filter(all_scores, pipeline == "SOTA_pipeline" & type =="cross"), size = 1.8) +
  geom_point(aes(shape = type, size = type, alpha = pipeline))+
  geom_point(data = filter(all_scores, type == "within lab"), aes(shape = type, size = type))+
  scale_shape_manual(values=c(16, 17))+
  scale_size_manual(values=c(2.0, 2.2)) +
  scale_alpha_manual(values=c(0.4, 0.4, 0.4, 0.4, 1)) +
  facet_grid(level ~ testing_cohort) +  ylim(0,1) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =12),
        axis.text.y = element_text(size = 12),
        axis.title=element_text(size=15),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14)) + ylab("AUC")
dev.off()

