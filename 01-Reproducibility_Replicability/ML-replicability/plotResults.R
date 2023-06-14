library(dplyr)
library(tidyr)
library(ggplot2)

all_scores <- read.csv("~/Desktop/thesis/vaginalMicrobiome/01-Reproducibility_Replicability/ML-replicability/all_scores.csv", sep=",", header = TRUE)
all_scores$level = factor(all_scores$level, levels=c("phylum", "class", "order", "family", "genus"), labels=c("phylum", "class", "order", "family", "genus")) 
#all_scores$training_cohort = factor(all_scores$training_cohort, levels=c("Antonio (N=15)", "Gressel (N=19)", "Angel (N=23)", "Tsementzi (N=29)", "Walsh (N=99)"), 
#                                    labels=c("Antonio (N=15)", "Gressel (N=19)", "Angel (N=23)", "Tsementzi (N=29)", "Walsh (N=99)"))
#all_scores$pipeline = factor(all_scores$pipeline, levels=c("Angel", "Antonio", "Gressel", "Tsementzi"),
#                             labels=c("Angel pipeline", "Antonio_Walsh pipeline", "Gressel pipeline", "Tsementzi pipeline"))

subset_data <- subset(all_scores, level=="genus")
pdf("~/Desktop/ggplot.pdf",  width=24, height=24)
a <- ggplot(subset_data, aes(x = training_cohort, y=AUC)) +
  geom_boxplot() + 
  geom_jitter(aes(colour=testing_cohort), size=10.0, alpha =0.9, width = 0.2) + 
  facet_grid(pipeline ~ .) +  ylim(0,1)+
  theme(axis.text.x = element_text(size=34),
        axis.title=element_text(size=37),
        axis.text.y = element_text(size=37),
        strip.text = element_text(size = 37),
        legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(2, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=32), #change legend title font size
        legend.text = element_text(size=32))
print(a)
dev.off()

ggplot(all_scores, aes(x=level, y=AUC, group=testing_cohort)) +
  geom_line(aes(color=testing_cohort), alpha=0.6, size=2.0)+ 
  facet_grid(pipeline ~ training_cohort) +
  geom_point(aes(colour=testing_cohort), alpha=0.6, size=2.0) + 
  ylim(0,1) + ggtitle("Testing AUC across datasets") +
  theme(legend.position="bottom", axis.text.x = element_text(angle = 45, vjust = 0.3, hjust=0.3))+
  ylab("AUC") 

