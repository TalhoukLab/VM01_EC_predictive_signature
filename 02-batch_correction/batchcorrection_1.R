library(microViz)
source(here::here("~/Desktop/thesis/VM02_meta_analysis/src/dataPrep.R"), encoding = "UTF-8")
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

batchcorrect <- function(phylo_obj, tax_level, tester){
  phylo_obj_filtered <- phyloseq::filter_taxa(phylo_obj, function(x) sum(x > 0) > (0.05*length(x)), TRUE)

  
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
    # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
    ord_calc(method = "PCoA") %>%
    ord_plot(color = "cohort", shape = "histology",  size = 3) 
  
  
  counts.aggregated.raw.pse <- counts.aggregated.raw + 0
  counts.aggregated.raw.pse.tss <- data.frame(t(apply(counts.aggregated.raw.pse, 1, function(x){x/sum(x)})))
  
  phylo_obj_tss <- phyloseq(otu_table(counts.aggregated.raw.pse.tss, taxa_are_rows = FALSE), tax_table(phylo_obj_agg),
                            sample_data(phylo_obj_agg), phy_tree(phylo_obj_agg))
  phylo_obj_tss %>%
    tax_transform("log", rank = NA) %>%
    # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
    ord_calc(method = "PCA") %>%
    ord_plot(color = "cohort", shape = "histology",  size = 3, plot_samples = TRUE)
  
  covar_mat <- cbind(as.factor(vm.metadata$histology), as.factor(vm.metadata$Primers))
  
  corrected_combatseq <- t(ComBat_seq(t(counts.aggregated.raw), batch = vm.batch, group = vm.trt))
  phylo_obj_agg_combatSeq <- phyloseq(otu_table(corrected_combatseq, taxa_are_rows = FALSE), tax_table(phylo_obj_agg),
                            sample_data(phylo_obj_agg), phy_tree(phylo_obj_agg))
  phylo_obj_agg_combatSeq %>%
    tax_transform("identity", rank = NA) %>%
    dist_calc("aitchison") %>%
    # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
    ord_calc(method = "PCoA") %>%
    ord_plot(color = "cohort", shape = "histology",  size = 3)
  
  corrected_combat <- combat(data = counts.aggregated.raw.pse.tss, batch = vm.batch, mod = vm.mod)
  phylo_obj_agg_combat <- phyloseq(otu_table(corrected_combat, taxa_are_rows = FALSE), tax_table(phylo_obj_agg),
                                      sample_data(phylo_obj_agg), phy_tree(phylo_obj_agg))
  phylo_obj_agg_combat %>%
    dist_calc(data, dist = "euclidean") %>%
    # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
    ord_calc(method = "PCoA") %>%
    ord_plot(color = "cohort", shape = "histology",  size = 3)
  
  corrected_limma <- limma(data = counts.aggregated.raw.pse.tss, batch = vm.batch, mod = vm.mod)
  phylo_obj_agg_limma <- phyloseq(otu_table(corrected_limma, taxa_are_rows = FALSE), tax_table(phylo_obj_agg),
                                   sample_data(phylo_obj_agg), phy_tree(phylo_obj_agg))
  phylo_obj_agg_limma %>%
    # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
    ord_calc(method = "PCA") %>%
    ord_plot(color = "cohort", shape = "histology",  size = 3)
  corrected_plsda <- myplsda(data=counts.aggregated.raw.pse.tss, trt = vm.trt, batch = vm.batch)
  phylo_obj_agg_plsda <- phyloseq(otu_table(corrected_plsda, taxa_are_rows = FALSE), tax_table(phylo_obj_agg),
                                  sample_data(phylo_obj_agg), phy_tree(phylo_obj_agg))
  phylo_obj_agg_plsda %>%
    # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
    ord_calc(method = "PCA") %>%
    ord_plot(color = "cohort", shape = "histology",  size = 3)
  corrected_splsda <- mysplsda(data=counts.aggregated.raw.pse.tss, trt = vm.trt, batch = vm.batch)
  phylo_obj_agg_splsda <- phyloseq(otu_table(corrected_splsda, taxa_are_rows = FALSE), tax_table(phylo_obj_agg),
                                  sample_data(phylo_obj_agg), phy_tree(phylo_obj_agg))
  phylo_obj_agg_splsda %>%
    # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
    ord_calc(method = "PCA") %>%
    ord_plot(color = "cohort", shape = "histology",  size = 3)
  
  corrected_conqur = ConQuR(tax_tab = counts.aggregated.raw.pse, batchid = batchid, 
                            covariates = vm.trt, batch_ref = 3,
                           logistic_lasso = T, quantile_type = "lasso", interplt = T, num_core = 3)
  corrected_conqur.tss <- data.frame(t(apply(corrected_conqur, 1, function(x){x/sum(x)})))

  phylo_obj_agg_conqur <- phyloseq(otu_table(corrected_conqur, taxa_are_rows = FALSE), tax_table(phylo_obj_agg),
                                      sample_data(phylo_obj_agg), phy_tree(phylo_obj_agg))
  phylo_obj_agg_conqur %>%
    tax_transform("identity", rank = NA) %>%
    dist_calc("aitchison") %>%
    # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
    ord_calc(method = "PCoA") %>%
    ord_plot(color = "cohort", shape = "histology",  size = 3)
  #corrected_bmc <- bmc(data = counts.aggregated.raw.pse.tss, batch = vm.batch)
  corrected_pn <- pn(level = tax_level, tester = tester)
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
  return(res)
}