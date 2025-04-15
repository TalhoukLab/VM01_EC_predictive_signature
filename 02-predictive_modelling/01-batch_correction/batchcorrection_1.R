library(microViz)
source(here::here("~/Desktop/thesis/VM01_EC_predictive_signature/02-predictive_modelling/01-batch_correction/dataPrep.R"), encoding = "UTF-8")
source(here::here("~/Desktop/thesis/VM02_meta_analysis/src/batchcorrection_functions.R"), encoding = "UTF-8")


testing <- c("Antonio_new")

for(tester in testing){
  raw_otus <- read.delim(file.path(paste0('~/Desktop/thesis/VM01_reproducibility_replicability/dada2_pipeline/results/testing/', tester, "/all.otutab.txt")))
  rep_seqs <- readDNAStringSet(file.path(paste0('~/Desktop/thesis/VM01_reproducibility_replicability/dada2_pipeline/results/testing/', tester, '/all.otus.fasta')))
  raw_tax <- data.frame(read.delim(file.path(paste0('~/Desktop/thesis/VM01_reproducibility_replicability/dada2_pipeline/results/testing/', tester,"/all_mergedmajority_taxonomy.tsv")), 
                                   sep = ";", header = FALSE))
  phylo_tree <- read_tree(file.path(paste0('~/Desktop/thesis/VM01_reproducibility_replicability/dada2_pipeline/results/testing/', tester, "/otus.tree")))
  feature_table <- read.delim(file.path(paste0('~/Desktop/thesis/VM02_meta_analysis/src/results/all_cohortsFT.csv')), header = TRUE, sep = ",")
  feature_table <- feature_table %>% filter(cohort != "Antonio")
  res <- dataPrep(raw_otus = raw_otus,raw_tax =  raw_tax, phylo_tree = phylo_tree, feature_table = feature_table)
  assign(paste0(tester, "_res"),res ,.GlobalEnv)
}

#Batch correction function 
library(PLSDAbatch)
library(mixOmics)
library(sva)
library(ConQuR)
library(doParallel)


phylo_obj <- Antonio_new_res
tax_level <- "genus"
phylo_obj_filtered <- phyloseq::filter_taxa(phylo_obj$raw, function(x) sum(x > 0) > (0.05*length(x)), TRUE)

vm.metadata <- data.frame(sample_data(phylo_obj_filtered))
vm.batch <- as.factor(vm.metadata$cohort)
batchid <- as.factor(as.numeric(vm.batch)-1)
vm.trt <- as.factor(vm.metadata$histology)
vm.primers <- as.factor(vm.metadata$Primers)
names(vm.batch) <- names(vm.trt) <- rownames(vm.metadata)
vm.mod <- model.matrix( ~ vm.trt)
rownames(vm.mod) <- names(vm.trt)

if(tax_level == "ASV"){
  counts.aggregated.raw <- t(data.frame(otu_table(phylo_obj_filtered)))
} else {
    phylo_obj_agg <- tax_glom(phylo_obj_filtered, taxrank = tax_level, 
                         NArm = TRUE, bad_empty = c(NA, "", " ", "\t"))
    idx <- which(rank_names(phylo_obj_agg)==tax_level)
    taxa_names(phylo_obj_agg) <- paste0(data.frame(tax_table(phylo_obj_agg))[,idx-3], "_",
                                    data.frame(tax_table(phylo_obj_agg))[,idx-2], "_", data.frame(tax_table(phylo_obj_agg))[,idx-1], "_",
                                    data.frame(tax_table(phylo_obj_agg))[,idx])
    counts.aggregated.raw <- t(data.frame(otu_table(phylo_obj_agg)))
}
phylo_obj_agg %>%
  tax_transform("identity", rank = NA) %>%
  dist_calc("aitchison") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "cohort", shape = "histology",  size = 3) 


counts.aggregated.raw.pse <- counts.aggregated.raw + 0
counts.aggregated.raw.pse.tss <- data.frame(t(apply(counts.aggregated.raw.pse, 1, function(x){x/sum(x)})))

phylo_obj_tss <- phyloseq(otu_table(counts.aggregated.raw.pse.tss, taxa_are_rows = FALSE), tax_table(phylo_obj_agg),
                          sample_data(phylo_obj_agg), phy_tree(phylo_obj_agg))
raw_fig <- phylo_obj_tss %>%
  tax_transform("log", rank = NA) %>%
  dist_calc(data, dist = "euclidean") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "cohort", shape = "histology",  size = 12, plot_samples = TRUE) +
  scale_shape_manual(values = c(1, 16))  +
  scale_color_brewer(palette = "Spectral")+
  ggside::geom_xsideboxplot(aes(fill = histology, y = histology), orientation = "y", outlier.size = 0.5) +
  ggside::geom_ysideboxplot(aes(fill = histology, x = histology), orientation = "x", outlier.size = 0.5) + 
  ggside::geom_xsideboxplot(aes(fill = cohort, y = cohort), orientation = "y", outlier.size = 0.5) +
  ggside::geom_ysideboxplot(aes(fill = cohort, x = cohort), orientation = "x", outlier.size = 0.5) + 
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 32),
                                      axis.text.y = element_text(size = 32),
                                      axis.title=element_text(size=33),
                                      strip.text.x = element_text(size = 32),
                                      legend.key.size = unit(1, 'cm'),
                                      legend.key.height = unit(1, 'cm'),
                                      legend.key.width = unit(1, 'cm'),
                                      legend.position="bottom", 
                                      legend.title = element_text(size=33, face = "bold"),
                                      legend.text = element_text(size=25)) 

covar_mat <- cbind(as.factor(vm.metadata$histology), as.factor(vm.metadata$Primers))

## combatseq
corrected_combatseq <- t(ComBat_seq(t(counts.aggregated.raw), batch = vm.batch, group = vm.trt))
phylo_obj_agg_combatSeq <- phyloseq(otu_table(corrected_combatseq, taxa_are_rows = FALSE), tax_table(phylo_obj_agg),
                          sample_data(phylo_obj_agg), phy_tree(phylo_obj_agg))
phylo_obj_agg_combatSeq %>%
  tax_transform("identity", rank = NA) %>%
  dist_calc("aitchison") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "cohort", shape = "histology",  size = 3)

## combat
corrected_combat <- combat(data = counts.aggregated.raw.pse.tss, batch = vm.batch, mod = vm.mod)
phylo_obj_agg_combat <- phyloseq(otu_table(corrected_combat, taxa_are_rows = FALSE), tax_table(phylo_obj_agg),
                                    sample_data(phylo_obj_agg), phy_tree(phylo_obj_agg))
combat_fig <- phylo_obj_agg_combat %>%
  dist_calc(data, dist = "euclidean") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "cohort", shape = "histology",  size =12)+
  scale_shape_manual(values = c(1, 16)) +
  scale_color_brewer(palette = "Spectral")+
  ggside::geom_xsideboxplot(aes(fill = histology, y = histology), orientation = "y", outlier.size = 0.5) +
  ggside::geom_ysideboxplot(aes(fill = histology, x = histology), orientation = "x", outlier.size = 0.5) + 
  ggside::geom_xsideboxplot(aes(fill = cohort, y = cohort), orientation = "y", outlier.size = 0.5) +
  ggside::geom_ysideboxplot(aes(fill = cohort, x = cohort), orientation = "x", outlier.size = 0.5) + 
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 32),
                                      axis.text.y = element_text(size = 32),
                                      axis.title=element_text(size=33),
                                      strip.text.x = element_text(size = 32),
                                      legend.key.size = unit(1, 'cm'),
                                      legend.key.height = unit(1, 'cm'),
                                      legend.key.width = unit(1, 'cm'),
                                      legend.position="bottom", 
                                      legend.title = element_text(size=33, face = "bold"),
                                      legend.text = element_text(size=25)) 



## LIMMA
corrected_limma <- limma(data = counts.aggregated.raw.pse.tss, batch = vm.batch, mod = vm.mod)
phylo_obj_agg_limma <- phyloseq(otu_table(corrected_limma, taxa_are_rows = FALSE), tax_table(phylo_obj_agg),
                                 sample_data(phylo_obj_agg), phy_tree(phylo_obj_agg))
phylo_obj_agg_limma %>%
  ord_calc(method = "PCA") %>%
  ord_plot(color = "cohort", shape = "histology",  size = 3)

## PLSDA
corrected_plsda <- myplsda(data=counts.aggregated.raw.pse.tss, trt = vm.trt, batch = vm.batch)
phylo_obj_agg_plsda <- phyloseq(otu_table(corrected_plsda, taxa_are_rows = FALSE), tax_table(phylo_obj_agg),
                                sample_data(phylo_obj_agg), phy_tree(phylo_obj_agg))
phylo_obj_agg_plsda %>%
  ord_calc(method = "PCA") %>%
  ord_plot(color = "cohort", shape = "histology",  size = 3)


## SPLSDA
corrected_splsda <- mysplsda(data=counts.aggregated.raw.pse.tss, trt = vm.trt, batch = vm.batch)
phylo_obj_agg_splsda <- phyloseq(otu_table(corrected_splsda, taxa_are_rows = FALSE), tax_table(phylo_obj_agg),
                                sample_data(phylo_obj_agg), phy_tree(phylo_obj_agg))
phylo_obj_agg_splsda %>%
  ord_calc(method = "PCA") %>%
  ord_plot(color = "cohort", shape = "histology",  size = 3)

## Conqur
corrected_conqur = ConQuR(tax_tab = counts.aggregated.raw.pse, batchid = batchid, 
                          covariates = vm.trt, batch_ref = 3,
                         logistic_lasso = T, quantile_type = "lasso", interplt = T, num_core = 3)
corrected_conqur.tss <- data.frame(t(apply(corrected_conqur, 1, function(x){x/sum(x)})))

phylo_obj_agg_conqur <- phyloseq(otu_table(corrected_conqur, taxa_are_rows = FALSE), tax_table(phylo_obj_agg),
                                    sample_data(phylo_obj_agg), phy_tree(phylo_obj_agg))
phylo_obj_agg_conqur %>%
  tax_transform("identity", rank = NA) %>%
  dist_calc("aitchison") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "cohort", shape = "histology",  size = 3)

#corrected_bmc <- bmc(data = counts.aggregated.raw.pse.tss, batch = vm.batch)
#corrected_pn <- pn(level = tax_level, tester = tester)
#corrected_pn <- data.frame(t(apply(corrected_pn, 1, function(x){x/sum(x)})))

res = list(raw = counts.aggregated.raw,
           raw.tss = counts.aggregated.raw.pse.tss,
          combat = corrected_combat,
          combatseq = corrected_combatseq,
          limma = corrected_limma,
          plsda = corrected_plsda,
          splsda = corrected_splsda,
          conqur = corrected_conqur.tss,
          #bmc = corrected_bmc,
           #pn = corrected_pn,
           labels = vm.trt,
           batches = vm.batch)

library(patchwork)
pdf("~/Desktop/batchCorrection.pdf", width = 34, height = 18)
(raw_fig | combat_fig)
dev.off()
