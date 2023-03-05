library(dplyr)
library(tidyr)
library(ggplot2)

all_scores <- read.csv("~/Desktop/thesis/vaginalMicrobiome/01-Reproducibility_Replicability/ML-replicability/all_modesScores.csv", sep=",", header = TRUE)
all_scores$level = factor(all_scores$level, levels=c("phylum", "class", "order", "family", "genus"), labels=c("phylum", "class", "order", "family", "genus")) 
all_scores$training_cohort = factor(all_scores$training_cohort, levels=c("Antonio (N=15)", "Gressel (N=19)", "Angel (N=23)", "Tsementzi (N=29)", "Walsh (N=99)"), 
                                    labels=c("Antonio (N=15)", "Gressel (N=19)", "Angel (N=23)", "Tsementzi (N=29)", "Walsh (N=99)"))
all_scores$pipeline = factor(all_scores$pipeline, levels=c("Angel", "Antonio", "Gressel", "Tsementzi"),
                             labels=c("Angel pipeline", "Antonio_Walsh pipeline", "Gressel pipeline", "Tsementzi pipeline"))


ggplot(all_scores, aes(x=level, y=AUC, group=testing_cohort)) +
  geom_line(aes(color=testing_cohort), alpha=0.6, size=2.0)+ 
  facet_grid(pipeline ~ training_cohort) +
  geom_point(aes(colour=testing_cohort), alpha=0.6, size=2.0) + 
  ylim(0,1) + ggtitle("Testing AUC across datasets") +
  theme(legend.position="bottom", axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=0.3))+
  ylab("AUC") 

