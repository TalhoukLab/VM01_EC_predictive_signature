## Chao microbiome only stacked 
library(yardstick)
library(themis)
library(tidymodels)
library(stacks)
library(data.table)

source(here::here("~/Desktop/thesis/VM02_meta_analysis/src/dataPrep.R"), encoding = "UTF-8")
cohorts <- c("Antonio_train", "Antonio_test", "Antonio", "Antonio1", "Chao_train", "Chao_test", "Chao",
             "Gressel_train", "Gressel_test", "Gressel", "Gressel_forward", "Gressel_forward_train", "Gressel_forward_test", 
             "Tsementzi_train", "Tsementzi_test", "Tsementzi",
             "Walsh_train", "Walsh_test", "Walsh")
for(tester in cohorts){
  raw_otus <- read.delim(file.path(paste0('~/Desktop/thesis/VM01_reproducibility_replicability/dada2_pipeline/results/', tester, "/all.otutab.txt")))
  rep_seqs <- readDNAStringSet(file.path(paste0('~/Desktop/thesis/VM01_reproducibility_replicability/dada2_pipeline/results/', tester, '/all.otus.fasta')))
  raw_tax <- data.frame(read.delim(file.path(paste0('~/Desktop/thesis/VM01_reproducibility_replicability/dada2_pipeline/results/', tester,"/all_mergedmajority_taxonomy.tsv")), 
                                   sep = ";", header = FALSE))
  phylo_tree <- read_tree(file.path(paste0('~/Desktop/thesis/VM01_reproducibility_replicability/dada2_pipeline/results/', tester, "/otus.tree")))
  feature_table <- read.delim(file.path(paste0('~/Desktop/thesis/VM01_reproducibility_replicability/00-helperfiles/', tester, 'FT.csv')), header = TRUE, sep = ",")
  #feature_table <- feature_table %>% filter(cohort != tester)
  res <- dataPrep(raw_otus, raw_tax, phylo_tree, feature_table)
  assign(paste0(tester, "_solo_res"),res ,.GlobalEnv)
}

microbiome_my_tax_glom <- function(phylo_obj, level){
  aggregated <- tax_glom(phylo_obj, taxrank = level, 
                         NArm = TRUE, bad_empty = c(NA, "", " ", "\t"))
  idx <- which(rank_names(aggregated)==level)
  taxa_names(aggregated) <- paste0(data.frame(tax_table(aggregated))[,idx-2], "_",
                                   data.frame(tax_table(aggregated))[,idx-1], "_", 
                                   data.frame(tax_table(aggregated))[,idx])
  counts.aggregated.raw <- t(data.frame(otu_table(aggregated)))
  counts.aggregated.raw <- janitor::clean_names(counts.aggregated.raw)
  counts.aggregated.raw.pse <- counts.aggregated.raw + 1
  counts.aggregated.raw.pse.tss <- data.frame(t(apply(counts.aggregated.raw.pse, 1, function(x){x/sum(x)})))
  counts.aggregated.raw.pse.tss.log <- log10(counts.aggregated.raw.pse.tss)
  counts.aggregated.raw.pse.tss.log <- cbind(counts.aggregated.raw.pse.tss.log,
                                             "labels"=sample_data(phylo_obj)$histology)
  res <- list(raw = counts.aggregated.raw, tss = counts.aggregated.raw.pse.tss.log)
  return(res)
}

microbiome_fncols_antonio <- function(data, cname) {
  add <-cname[!cname%in%names(data)]
  if(length(add)!=0) data[add] <- 0
  data
}


microbiome_prep_test <- function(data, trained_data){
  data <- data.frame(janitor::clean_names(data))
  otu_test <- microbiome_fncols_antonio(data, setdiff(colnames(trained_data), colnames(data)))
  otu_test <- base::subset(otu_test, select = colnames(trained_data))
  shannon_list <- otu_test$shannon
  if("shannon" %in% colnames(otu_test)){
    otu_test <- subset(otu_test, select = -c(shannon, labels) )
  } else {
    otu_test <- otu_test
  }
  otu_test <- otu_test + 1
  otu_test.pse.tss <- data.frame(t(apply(otu_test, 1, function(x){x/sum(x)})))
  otu_test.pse.tss.log <- log10(otu_test.pse.tss)
  if("shannon" %in% colnames(trained_data)){
    otu_test.pse.tss.log$shannon <- shannon_list
  }
  return(otu_test.pse.tss.log)
}

microbiome_get_metrics_indi <- function(data, model, level, labels, cohort){
  preds <- microbiome_prep_test(data = data, trained_data = model$pre$actions$recipe$recipe$template) %>%
    bind_cols(predict(model,
                      ., type = "prob"))%>%
    bind_cols("labels"=labels)
  
  preds$labels <- factor(preds$labels,  levels = c("Benign", "EC"))
  preds <- preds %>%
    mutate(predicted = case_when(.pred_Benign > .pred_EC ~ "Benign",
                                 .pred_EC > .pred_Benign ~ "EC"))
  preds$predicted <- factor(preds$predicted, levels = c("Benign", "EC"))
  metric <-  metric_set(f_meas, precision, recall, specificity, npv)
  metrics_my <- metric(preds, truth = labels, estimate = predicted, event_level = "second", estimator = "binary")
  metrics_my$type <- level
  res <- list(met = metrics_my, predictions = preds)
  return(res)
}

microbiome_get_metrics_indi_tr <- function(data, model, level, labels, cohort){
  preds <- predict(model, data, type = "prob") %>%
    bind_cols("labels"=labels)
  
  preds$labels <- factor(preds$labels,  levels = c("Benign", "EC"))
  preds <- preds %>%
    mutate(predicted = case_when(.pred_Benign > .pred_EC ~ "Benign",
                                 .pred_EC > .pred_Benign ~ "EC"))
  preds$predicted <- factor(preds$predicted, levels = c("Benign", "EC"))
  metric <-  metric_set(f_meas, precision, recall, specificity, npv)
  metrics_my <- metric(preds, truth = labels, estimate = predicted, event_level = "second", estimator = "binary")
  metrics_my$type <- level
  res <- list(met = metrics_my, predictions = preds)
  return(res)
}

library(doMC)
registerDoMC(cores = 5)

shannon_test_chao <- as.data.frame(microbiome::alpha(Antonio1_solo_res$raw, 
                                                     index = c("diversity_shannon")))
test_fitered_chao <- phyloseq::filter_taxa(Antonio1_solo_res$raw, function(x) sum(x > 0) > (0.05*length(x)), TRUE)
test_obj_agg_chao <- tax_glom(test_fitered_chao, taxrank = "genus", 
                              NArm = TRUE, bad_empty = c(NA, "", " ", "\t"))
idx <- which(rank_names(test_obj_agg_chao)=="genus")
taxa_names(test_obj_agg_chao) <- paste0(data.frame(tax_table(test_obj_agg_chao))[,idx-2], "_",
                                        data.frame(tax_table(test_obj_agg_chao))[,idx-1], "_", 
                                        data.frame(tax_table(test_obj_agg_chao))[,idx])
test_data_chao <- data.frame(t(data.frame(otu_table(test_obj_agg_chao))))
test_data_chao$shannon <- (shannon_test_chao$diversity_shannon)

## Chao 
microbiome_chao_rf_grid <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Chao_microbiome_rf_grid.rds")
microbiome_chao_rf_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Chao_microbiome_rf_fit.rds")
microbiome_chao_nnet_grid <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Chao_microbiome_nnet_grid.rds")
microbiome_chao_nnet_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Chao_microbiome_nnet_fit.rds")
microbiome_chao_xgb_grid <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Chao_microbiome_xgb_grid.rds")
microbiome_chao_xgb_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Chao_microbiome_xgb_fit.rds")

## Gressel 
microbiome_gressel_rf_grid <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Gressel_microbiome_rf_grid.rds")
microbiome_gressel_rf_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Gressel_microbiome_rf_fit.rds")
microbiome_gressel_nnet_grid <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Gressel_microbiome_nnet_grid.rds")
microbiome_gressel_nnet_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Gressel_microbiome_nnet_fit.rds")
microbiome_gressel_xgb_grid <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Gressel_microbiome_xgb_grid.rds")
microbiome_gressel_xgb_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Gressel_microbiome_xgb_fit.rds")

## Tsementzi 
microbiome_tsementzi_rf_grid <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Tsementzi_microbiome_rf_grid.rds")
microbiome_tsementzi_rf_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Tsementzi_microbiome_rf_fit.rds")
microbiome_tsementzi_nnet_grid <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Tsementzi_microbiome_nnet_grid.rds")
microbiome_tsementzi_nnet_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Tsementzi_microbiome_nnet_fit.rds")
microbiome_tsementzi_xgb_grid <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Tsementzi_microbiome_xgb_grid.rds")
microbiome_tsementzi_xgb_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Tsementzi_microbiome_xgb_fit.rds")

## Walsh
microbiome_walsh_rf_grid <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Walsh_microbiome_rf_grid.rds")
microbiome_walsh_rf_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Walsh_microbiome_rf_fit.rds")
microbiome_walsh_nnet_grid <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Walsh_microbiome_nnet_grid.rds")
microbiome_walsh_nnet_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Walsh_microbiome_nnet_fit.rds")
microbiome_walsh_xgb_grid <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Walsh_microbiome_xgb_grid.rds")
microbiome_walsh_xgb_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Walsh_microbiome_xgb_fit.rds")

microbiome_get_metrics <- function(data, col){
  data <- data %>% mutate(class = case_when((col) >=0.5 ~ "EC",
                                            (col) < 0.5 ~ "Benign"))
  data$class <- factor(data$class, levels = c("Benign", "EC"))
  data$col <- col
  metric <- metric_set(npv, recall, precision, specificity, f_meas)
  metrics_my <- metric(data, truth = labels, estimate = class, event_level = "second")
  metrics_my[nrow(metrics_my) + 1,] <-  roc_auc(data = data, truth = labels, col, event_level = "second")
  return(metrics_my)
}

get_results_microbiome <- function(cohort, type){
  if(type == "training"){
    rf_results <- microbiome_get_metrics_indi_tr(data = eval(parse(text = paste0("microbiome_", cohort, "_rf_fit")))$pre$actions$recipe$recipe$template, 
                                   model = eval(parse(text = paste0("microbiome_", cohort, "_rf_fit"))), 
                                   level = "genus", labels = eval(parse(text = paste0("microbiome_", cohort, "_rf_fit")))$pre$mold$outcomes$labels, cohort = cohort)
    rf_results$met$model <- "rf"
    nnet_results<- microbiome_get_metrics_indi_tr(data = eval(parse(text = paste0("microbiome_", cohort, "_rf_fit")))$pre$mold$predictors, 
                                    model = eval(parse(text = paste0("microbiome_", cohort, "_nnet_fit"))), 
                                    level = "genus", labels = eval(parse(text = paste0("microbiome_", cohort, "_rf_fit")))$pre$mold$outcomes$labels, cohort = cohort)
    nnet_results$met$model <- "nnet"
    xgb_results <- microbiome_get_metrics_indi_tr(data = eval(parse(text = paste0("microbiome_", cohort, "_rf_fit")))$pre$mold$predictors, 
                                    model = eval(parse(text = paste0("microbiome_", cohort, "_xgb_fit"))), 
                                    level = "genus", labels = eval(parse(text = paste0("microbiome_", cohort, "_rf_fit")))$pre$mold$outcomes$labels, cohort = cohort)
    xgb_results$met$model <- "xgb"
  }
  if(type == "testing"){
    rf_results <- microbiome_get_metrics_indi(data = test_data_chao, model = eval(parse(text = paste0("microbiome_", cohort, "_rf_fit"))), level = "genus", labels = sample_data(Antonio1_solo_res$raw)$histology, cohort = cohort)
    rf_results$met$model <- "rf"
    nnet_results <- microbiome_get_metrics_indi(data = test_data_chao, model = eval(parse(text = paste0("microbiome_", cohort, "_nnet_fit"))), level = "genus", labels = sample_data(Antonio1_solo_res$raw)$histology, cohort = cohort)
    nnet_results$met$model <- "nnet"
    xgb_results <- microbiome_get_metrics_indi(data = test_data_chao, model = eval(parse(text = paste0("microbiome_", cohort, "_xgb_fit"))), level = "genus", labels = sample_data(Antonio1_solo_res$raw)$histology, cohort = cohort)
    xgb_results$met$model <- "xgb"
  }
  
  rf_pred <- rf_results$predictions %>% 
    dplyr::select(c(.pred_EC)) %>%
    dplyr::rename(!!paste0(cohort, "_rf") := .pred_EC)
  
  nnet_pred <- nnet_results$predictions %>%
    dplyr::select(c(.pred_EC)) %>%
    dplyr::rename(!!paste0(cohort, "_nnet") := .pred_EC)
  
  xgb_pred <- xgb_results$predictions %>%
    dplyr::select(c(.pred_EC)) %>%
    dplyr::rename(!!paste0(cohort, "_xgb") := .pred_EC)
  combined <- cbind(rf_pred,nnet_pred, xgb_pred)
  if(type == "training"){
    combined <- combined %>% 
      bind_cols("labels" = factor(eval(parse(text = paste0("microbiome_", cohort, "_rf_fit")))$pre$mold$outcomes$labels, levels = c("Benign", "EC")))
    combined$ensemble <- rowSums(combined[,c(1:3)])/3 
    combined$rf_nnet <- rowSums(combined[,c(1:2)])/2 
    combined$rf_xgb <- rowSums(combined[,c(1,3)])/2 
    combined$nnet_xgb <- rowSums(combined[,c(2:3)])/2 
  } else {
    combined <- combined %>%
      bind_cols("labels" = factor(sample_data(Antonio1_solo_res$raw)$histology, levels = c("Benign", "EC")),
                "sraID" = sample_data(Antonio1_solo_res$raw)$sraID)
    combined$ensemble <- rowSums(combined[,c(1:3)])/3
    combined$rf_nnet <- rowSums(combined[,c(1:2)])/2
    combined$rf_xgb <- rowSums(combined[,c(1,3)])/2
    combined$nnet_xgb <- rowSums(combined[,c(2:3)])/2
    combined <- combined %>% relocate(sraID, .after = nnet_xgb)
    
  }
  combined <- combined %>% relocate(labels, .after = nnet_xgb) %>%
                  dplyr::rename(!!paste0(cohort, "_ensemble") := ensemble,
                                !!paste0(cohort, "_rf_nnet") := rf_nnet,
                                !!paste0(cohort, "_rf_xgb") := rf_xgb,
                                !!paste0(cohort, "_nnet_xgb") := nnet_xgb)
  all_metrics <- lapply(combined[, c(1:7)], microbiome_get_metrics, data = combined)
  
  all_metrics_list <- Map(cbind, all_metrics, type = names(all_metrics))
  all_metrics_list <- data.frame(rbindlist(all_metrics_list))
  all <- all_metrics_list %>% select(c(.metric, .estimate, type))
  all_wide <- all %>% pivot_wider(names_from = .metric, values_from = .estimate)
  return(list(met = all_wide, pred = combined))
}

chao_microbiome_training <- get_results_microbiome(cohort = "chao",type = "training")
chao_microbiome_testing <- get_results_microbiome(cohort = "chao",type = "testing")

tsementzi_microbiome_training <- get_results_microbiome(cohort = "tsementzi",type = "training")
tsementzi_microbiome_testing <- get_results_microbiome(cohort = "tsementzi",type = "testing")

walsh_microbiome_training <- get_results_microbiome(cohort = "walsh",type = "training")
walsh_microbiome_testing <- get_results_microbiome(cohort = "walsh",type = "testing")

gressel_microbiome_training <- get_results_microbiome(cohort = "gressel",type = "training")
gressel_microbiome_testing <- get_results_microbiome(cohort = "gressel",type = "testing")

all_microbiome <- chao_microbiome_testing$pred %>%
  left_join(tsementzi_microbiome_testing$pred, by='sraID', copy = FALSE) %>% 
  left_join(walsh_microbiome_testing$pred, by = "sraID", copy = FALSE) %>%
  left_join(gressel_microbiome_testing$pred, by = "sraID", copy = FALSE)


all_microbiome <- all_microbiome[, c('sraID', 'labels.x', 'chao_rf', 'tsementzi_rf',  'tsementzi_xgb', 'walsh_rf', 'walsh_xgb', 'gressel_nnet')]

all_microbiome$ensemble <- rowSums(all_microbiome[,c(3:8)])/6

all_microbiome <- all_microbiome %>% mutate(ensemble_pred = case_when((ensemble) >=0.5 ~ "EC",
                                                (ensemble) < 0.5 ~ "Benign"))
all_microbiome$ensemble_pred <- factor(all_microbiome$ensemble_pred, levels = c("Benign", "EC"))
all_microbiome$labels <- factor(sample_data(Antonio1_solo_res$raw)$histology, levels = c("Benign", "EC"))

metric <- metric_set(f_meas, precision, recall, specificity, npv)
metrics_my_microbiome <- metric(all_microbiome, truth = labels, estimate = ensemble_pred, event_level = "second")
metrics_my_microbiome[nrow(metrics_my_microbiome) + 1,] <-  roc_auc(all_microbiome, truth = labels, ensemble, event_level = "second")
write.table(all_microbiome, "~/Desktop/microbiome_only.csv", quote = F, sep = ",", row.names = F)


