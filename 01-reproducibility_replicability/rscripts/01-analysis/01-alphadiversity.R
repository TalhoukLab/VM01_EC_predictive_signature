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

cohorts <- c("Antonio", "Chao", "Gressel", "Tsementzi", "Walsh")
pipelines <- c("Antonio", "Chao", "Tsementzi", "Gressel", "dada2")

my_Shannon <- function(x){
  x <- x[x > 0]
  S <- length(x)
  p <- x/sum(x)
  (-sum(p * log(p)))
}

sink("~/Desktop/LM.txt")
for (pipeline in pipelines){
  for(cohort in cohorts){
    print(paste0("cohort_pipeline", cohort, '_', pipeline))
    phylo <- eval(parse(text = paste0(cohort, '_', pipeline, 'phyloseq_tree')))
    otu <- abundances(phylo)
    to_keep <- which((colSums(otu))!=0)
    otu <- otu[, to_keep]
    ev <- apply(otu, 2, function(x) {
      my_Shannon(x)
    })
    ps1.meta <- meta(phylo)
    ps1.meta <- ps1.meta[to_keep,]
    ps1.meta$shannon <- as.data.frame(ev)$ev
    assign(paste0(cohort, "_", pipeline, "alphaDiversity"),ps1.meta,.GlobalEnv)
    if(cohort=="Antonio" || cohort=="Walsh"){
      test.sig <- lm(ps1.meta$shannon ~ as.factor(ps1.meta$histology) +
                        as.numeric(ps1.meta$BMI) + as.factor(ps1.meta$pHRecoded) +
                        as.factor(ps1.meta$menopausal.status))
    }
    if(cohort == "Tsementzi"){
      test.sig <- lm(ps1.meta$shannon ~ as.factor(ps1.meta$histology) +
                       as.numeric(ps1.meta$BMI) + as.factor(ps1.meta$pHRecoded))
    }
    if(cohort == "Chao"){
      test.sig <- lm(ps1.meta$shannon ~ as.factor(ps1.meta$histology) +
                       as.factor(ps1.meta$menopausal.status))
    }
    if(cohort == "Gressel"){
      test.sig <- lm(ps1.meta$shannon ~ as.factor(ps1.meta$histology))
    }
    print(summary(test.sig))
    assign(paste0(cohort, "_", pipeline, "alphaDiversity_sig"),test.sig,.GlobalEnv)
  }
  
  
  all_cohorts <- plyr::rbind.fill(eval(parse(text = paste0('Antonio_', pipeline, 'alphaDiversity'))),
                       eval(parse(text = paste0('Chao_', pipeline, 'alphaDiversity'))),
                       eval(parse(text = paste0('Gressel_', pipeline, 'alphaDiversity'))),
                       eval(parse(text = paste0('Tsementzi_', pipeline, 'alphaDiversity'))),
                       eval(parse(text = paste0('Walsh_', pipeline, 'alphaDiversity'))))
  
  pathology <- levels(as.factor(all_cohorts$histology))
  pathology.pairs <- combn(seq_along(pathology), 2, simplify = FALSE, FUN = function(i)pathology[i])
  
  all_cohorts_long <- gather(all_cohorts, metric, value, shannon, factor_key=TRUE)
  all_cohorts_long$log_val <- log(all_cohorts_long$value)
  all_cohorts_long$pipeline <- pipeline
  all_cohorts_long$cohort <- as.factor(all_cohorts_long$cohort)
  assign(paste0(pipeline, "alphaDiversity"),all_cohorts_long,.GlobalEnv)
  if(pipeline == "Gressel"){
    anno_df = compare_means(value ~ histology, group.by = c("cohort", "metric"), data = all_cohorts_long, 
                            method = "kruskal.test", p.adjust.method = "BH")
  } else {
    anno_df = compare_means(value ~ histology, group.by = c("cohort"), data = all_cohorts_long, 
                            method = "t.test", p.adjust.method = "BH")
  }
  #print(pipeline)
  print(anno_df)
}
sink()

all_pipelines <- rbind(AntonioalphaDiversity,
                       ChaoalphaDiversity,
                       GresselalphaDiversity,
                       TsementzialphaDiversity,
                       dada2alphaDiversity)
all_pipelines$pipeline[all_pipelines$pipeline == "Antonio"] <- "Antonio_Walsh"
all_pipelines$pipeline[all_pipelines$pipeline == "dada2"] <- "DADA2"
all_pipelines$pipeline <- paste0(all_pipelines$pipeline, "_pipeline")
all_pipelines$pipeline <- factor(all_pipelines$pipeline, levels = c("Antonio_Walsh_pipeline",
                                                                    "Tsementzi_pipeline",
                                                                    "Gressel_pipeline", 
                                                                    "Chao_pipeline",
                                                                    "DADA2_pipeline"))
all_pipelines$cohort <- factor(all_pipelines$cohort, levels = c("Antonio", "Walsh", "Tsementzi", "Gressel", "Chao"))
pdf("~/Desktop/alpha.pdf",  width=25, height=8)

bp1 <- ggplot(all_pipelines, aes(x=cohort, y=value, fill = histology)) +
  geom_rect(data = subset(all_pipelines,pipeline == 'DADA2_pipeline'),aes(fill = "whitesmoke"),
            xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf,alpha = 0.12) +
  geom_boxplot(aes(fill=histology), outlier.shape = NA) + 
  scale_fill_manual(values=c("#7CAE00", "#F8766D", "whitesmoke"), ) + 
  theme_bw()+ ggtitle("Alpha diversity") +
  facet_wrap(~pipeline, ncol=5, scales = "fixed") + ylab("Shannon index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title=element_text(size=21),
        strip.text.x = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.position="bottom", 
        legend.title = element_text(size=21, face = "bold"),
        legend.text = element_text(size=10), strip.background =element_rect(fill="white")) 
bp1
dev.off()