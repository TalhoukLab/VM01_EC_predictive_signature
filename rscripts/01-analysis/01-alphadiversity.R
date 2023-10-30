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

cohorts <- c("Antonio", "Chao", "Gressel", "Tsementzi", "Walsh")
pipelines <- c("Antonio", "Chao", "Tsementzi", "Gressel", "SOTA")

my_Shannon <- function(x){
  # Ignore zeroes
  x <- x[x > 0]
  # Species richness (number of species)
  S <- length(x)
  # Relative abundances
  p <- x/sum(x)
  
  # Normlainzed Shannon index
 # ((-sum(p * log(p)))/log(S))

    (-sum(p * log(p)))
  
}

sink("../vaginalMicrobiome/01-Reproducibility_Replicability/Results/01-AlphaDiversity/LM.txt")
for (pipeline in pipelines){
  for(cohort in cohorts){
    print(paste0("cohort_pipeline", cohort, '_', pipeline))
    phylo <- eval(parse(text = paste0(cohort, '_', pipeline, 'phyloseq_tree_raw')))
    #tab <- microbiome::alpha(phylo, index = c("shannon", "chao1", "pielou"))
    otu <- abundances(phylo)
    to_keep <- which((colSums(otu))!=0)
    otu <- otu[, to_keep]
    
    #metaSeqObject <- newMRexperiment(otu) 
    #CSS <- cumNorm(metaSeqObject, p=cumNormStat(metaSeqObject))
    #outs_CSS = data.frame(MRcounts(CSS, norm=TRUE, log=T))
    
    ev <- apply(otu, 2, function(x) {
      my_Shannon(x)
    })
    ps1.meta <- meta(phylo)
    #ps1.meta <- subset(ps1.meta, select = c("cohort", "sraID", "histology"))
    ps1.meta <- ps1.meta[to_keep,]
    ps1.meta$shannon <- as.data.frame(ev)$ev
    #tab$sraID <- ps1.meta$sraID
    #tab$histology <- ps1.meta$histology
    #tab$cohort <- ps1.meta$cohort
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
  assign(paste0(pipeline, "alphaDiversity"),all_cohorts_long,.GlobalEnv)
  anno_df = compare_means(value ~ histology, group.by = c("metric", "cohort"), data = all_cohorts_long, method = "kruskal.test")
  #print(pipeline)
  #print(anno_df)
  
}
sink()

all_pipelines <- rbind(AntonioalphaDiversity,
                       ChaoalphaDiversity,
                       GresselalphaDiversity,
                       TsementzialphaDiversity,
                       SOTAalphaDiversity)
all_pipelines$pipeline[all_pipelines$pipeline == "Antonio"] <- "Antonio_Walsh"
all_pipelines$pipeline <- paste0(all_pipelines$pipeline, "_pipeline")
all_pipelines$pipeline <- factor(all_pipelines$pipeline, levels = c("Antonio_Walsh_pipeline",
                                                                    "Tsementzi_pipeline",
                                                                    "Gressel_pipeline", 
                                                                    "Chao_pipeline",
                                                                    "SOTA_pipeline"))

all_pipelines$cohort <- factor(all_pipelines$cohort, levels = c("Antonio", "Walsh", "Tsementzi", "Gressel", "Chao"))

bp1 <- ggplot(all_cohorts_long, aes(x=cohort, y=value, fill = histology)) +
  geom_rect(data = subset(all_cohorts_long,pipeline == 'SOTA_pipeline'),aes(fill = pipeline),
            xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf,alpha = 0.12) +
  geom_boxplot(aes(fill=histology)) + 
  scale_fill_manual(values=c("yellowgreen", "tomato3", "whitesmoke")) + 
  theme_bw()+
  facet_wrap(~pipeline, ncol=5, scales = "fixed") + ylab("Shannon index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title=element_text(size=16),
        strip.text.x = element_text(size = 12),
        legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10), strip.background =element_rect(fill="white")) 
bp1
