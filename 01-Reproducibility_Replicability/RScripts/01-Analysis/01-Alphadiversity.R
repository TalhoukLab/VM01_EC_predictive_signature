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


source(file = "../vaginalMicrobiome/01-Reproducibility_Replicability/RScripts/00-DataPrep/Angel_dataPrep.R")
source(file = "../vaginalMicrobiome/01-Reproducibility_Replicability/RScripts/00-DataPrep/Antonio_dataPrep.R")
source(file = "../vaginalMicrobiome/01-Reproducibility_Replicability/RScripts/00-DataPrep/Gressel_dataPrep.R")
source(file = "../vaginalMicrobiome/01-Reproducibility_Replicability/RScripts/00-DataPrep/Tsementzi_dataPrep.R")

cohorts <- c("Antonio", "Angel", "Gressel", "Tsementzi", "Walsh")
pipelines <- c("Angel", "Antonio", "Tsementzi", "Gressel", "InHouse")

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

for (pipeline in pipelines){
  for(cohort in cohorts){
    phylo <- eval(parse(text = paste0(cohort, '_', pipeline, 'phyloseq_tree_raw')))
    #tab <- microbiome::alpha(phylo, index = c("shannon", "chao1", "pielou"))
    otu <- abundances(phylo)
    to_keep <- which((colSums(otu))>1)
    otu <- otu[, to_keep]
    metaSeqObject <- newMRexperiment(otu) 
    CSS <- cumNorm(metaSeqObject, p=cumNormStat(metaSeqObject))
    outs_CSS = data.frame(MRcounts(CSS, norm=TRUE, log=T))
    ev <- apply(outs_CSS, 2, function(x) {
      my_Shannon(x)
    })
    ps1.meta <- meta(phylo)
    ps1.meta <- subset(ps1.meta, select = c("cohort", "sraID", "histology"))
    ps1.meta <- ps1.meta[to_keep,]
    ps1.meta$shannon <- as.data.frame(ev)$ev
    #tab$sraID <- ps1.meta$sraID
    #tab$histology <- ps1.meta$histology
    #tab$cohort <- ps1.meta$cohort
    assign(paste0(cohort, "_", pipeline, "alphaDiversity"),ps1.meta,.GlobalEnv)
  }
  
  all_cohorts <- rbind(eval(parse(text = paste0('Antonio_', pipeline, 'alphaDiversity'))),
                       eval(parse(text = paste0('Angel_', pipeline, 'alphaDiversity'))),
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
  print(pipeline)
  print(anno_df)
  
}


all_pipelines <- rbind(AngelalphaDiversity,
                       AntonioalphaDiversity,
                       GresselalphaDiversity,
                       TsementzialphaDiversity,
                       InHousealphaDiversity)
bp1 <- ggplot(all_pipelines, aes(x=cohort, y=value, fill = histology)) +
  geom_boxplot(aes(fill=histology)) + 
  scale_fill_manual(values=c("yellowgreen", "tomato3")) + 
  facet_wrap(~pipeline, ncol=5, scales = "fixed") + ylab("Shannon index")+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title=element_text(size=24),
        strip.text.x = element_text(size = 23),
        legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(size=16),   legend.text = element_text(size=16)) 
bp1
