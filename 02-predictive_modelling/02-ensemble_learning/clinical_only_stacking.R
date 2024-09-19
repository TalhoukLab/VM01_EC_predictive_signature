## Chao microbiome only stacked 
library(yardstick)
library(themis)
library(tidymodels)
library(stacks)

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

my_tax_glom <- function(phylo_obj, level){
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

prep_test_chao <- function(data, trained_data){
  data <- data.frame(janitor::clean_names(data))
  trained_data <- subset(trained_data, select = -c(histology) )
  otu_test <- base::subset(data, select = colnames(trained_data))
  return(otu_test)
}

get_metrics_indi <- function(data, model, level, labels, cohort){
  preds <- prep_test_chao(data = data, trained_data = model$pre$actions$recipe$recipe$template) %>%
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

get_metrics_base <- function(data, model, level, labels){
  preds <- prep_test(data = data, trained_data = model$train) %>%
    bind_cols(predict(model,
                      ., type = "prob", members=TRUE))%>%
    bind_cols("labels"=labels)
  preds$labels <- factor(preds$labels,  levels = c("Benign", "EC"))
  preds <- preds %>%
    mutate(predicted = case_when(.pred_Benign > .pred_EC ~ "Benign",
                                 .pred_EC > .pred_Benign ~ "EC"))
  preds$predicted <- factor(preds$predicted, levels = c("Benign", "EC"))
  metric <-  metric_set(f_meas, precision, recall, specificity, npv)
  metrics_my <- metric(preds, truth = labels, estimate = predicted, event_level = "second", estimator = "binary")
  metrics_my$type <- level
  plot <- autoplot( model, type = "weights")
  res <- list(met = metrics_my, plot = plot, predictions = preds)
  return(res)
}

library(doMC)
registerDoMC(cores = 5)

shannon_test <- as.data.frame(microbiome::alpha(Antonio1_solo_res$raw, 
                                                index = c("diversity_shannon")))
test_fitered <- phyloseq::filter_taxa(Antonio1_solo_res$raw, function(x) sum(x > 0) > (0.05*length(x)), TRUE)
test_data <- my_tax_glom(test_fitered, level = "genus")$tss
test_data$shannon <- (shannon_test$diversity_shannon)
test_data$age <- (sample_data(Antonio1_solo_res$raw)$age)
test_data$bmi <- as.numeric(sample_data(Antonio1_solo_res$raw)$BMI)
test_data$pHRecoded <- factor(sample_data(Antonio1_solo_res$raw)$pHRecoded)
test_data$ethnicity <- factor("white", levels = c("white"), labels = c("white"))

test_data_tmp <- mice(test_data[, -189],m=1,maxit=50,meth='pmm',seed=123)
test_data_complete <- complete(test_data_tmp,1)
test_data_complete_sub <- test_data_complete %>% select(c(age, bmi, pHRecoded, ethnicity))


## Chao 
chao_rf_grid <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Chao_age_rf_grid.rds")
chao_rf_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Chao_age_rf_fit.rds")
chao_nnet_grid <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Chao_age_nnet_grid.rds")
chao_nnet_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Chao_age_nnet_fit.rds")
chao_xgb_grid <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Chao_age_xgb_grid.rds")
chao_xgb_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Chao_age_xgb_fit.rds")

## Tsementzi 
tsementzi_rf_grid <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Tsementzi_all_clinical_rf_grid.rds")
tsementzi_rf_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Tsementzi_all_clinical_rf_fit.rds")
tsementzi_nnet_grid <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Tsementzi_all_clinical_nnet_grid.rds")
tsementzi_nnet_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Tsementzi_all_clinical_nnet_fit.rds")
tsementzi_xgb_grid <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Tsementzi_all_clinical_xgb_grid.rds")
tsementzi_xgb_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Tsementzi_all_clinical_xgb_fit.rds")

## Walsh
walsh_rf_grid <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_all_clinical_rf_grid.rds")
walsh_rf_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_all_clinical_rf_fit.rds")
walsh_nnet_grid <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_all_clinical_nnet_grid.rds")
walsh_nnet_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_all_clinical_nnet_fit.rds")
walsh_xgb_grid <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_all_clinical_xgb_grid.rds")
walsh_xgb_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_all_clinical_xgb_fit.rds")

get_results_microbiome <- function(cohort){
  # stacks_model <- stacks() %>%
  #   add_candidates(candidates = eval(parse(text = paste0(cohort, "_rf_grid"))), name = "rf") %>%
  #   add_candidates(candidates = eval(parse(text = paste0(cohort, "_nnet_grid"))), name = "nnet") %>%
  #   add_candidates(candidates = eval(parse(text = paste0(cohort, "_xgb_grid"))), name = "bt") %>%
  #   blend_predictions() %>%
  #   fit_members()
  # stacks_results <- get_metrics_base(data = test_data_chao, model = stacks_model, level = "genus", labels = sample_data(test_obj_agg_chao)$histology)
  # stacks_results$met$model <- "stacks"
  # 
  # stacks_model_fmeas <- stacks() %>%
  #   add_candidates(candidates = eval(parse(text = paste0(cohort, "_rf_grid"))), name = "rf") %>%
  #   add_candidates(candidates = eval(parse(text = paste0(cohort, "_nnet_grid"))), name = "nnet") %>%
  #   add_candidates(candidates = eval(parse(text = paste0(cohort, "_xgb_grid"))), name = "bt") %>%
  #   blend_predictions(metric = metric_set(f_meas)) %>%
  #   fit_members()
  # stacks_results_fmeas <- get_metrics_base(data = test_data_chao, model = stacks_model_fmeas, level = "genus", labels = sample_data(test_obj_agg_chao)$histology)
  # stacks_results_fmeas$met$model <- "stacks_fmeas"

  rf_results <- get_metrics_indi(data = test_data_complete_sub, model = eval(parse(text = paste0(cohort, "_rf_fit"))), level = "genus", labels = sample_data(Antonio1_solo_res$raw)$histology, cohort = cohort)
  rf_results$met$model <- "rf"
  nnet_results <- get_metrics_indi(data = test_data_complete_sub, model = eval(parse(text = paste0(cohort, "_nnet_fit"))), level = "genus", labels = sample_data(Antonio1_solo_res$raw)$histology, cohort = cohort)
  nnet_results$met$model <- "nnet"
  xgb_results <- get_metrics_indi(data = test_data_complete_sub, model = eval(parse(text = paste0(cohort, "_xgb_fit"))), level = "genus", labels = sample_data(Antonio1_solo_res$raw)$histology, cohort = cohort)
  xgb_results$met$model <- "xgb"
  
  rf_pred <- rf_results$predictions %>% 
    dplyr::select(c(.pred_EC)) %>%
    dplyr::rename(!!paste0(cohort, "_rf") := .pred_EC)
  
  nnet_pred <- nnet_results$predictions %>%
    dplyr::select(c(.pred_EC)) %>%
    dplyr::rename(!!paste0(cohort, "_nnet") := .pred_EC)
  
  xgb_pred <- xgb_results$predictions %>%
    dplyr::select(c(.pred_EC)) %>%
    dplyr::rename(!!paste0(cohort, "_xgb") := .pred_EC) %>%
    bind_cols("labels" = factor(sample_data(Antonio1_solo_res$raw)$histology, levels = c("Benign", "EC")),
              "sraID" = sample_data(Antonio1_solo_res$raw)$sraID)
  combined <- cbind(rf_pred,nnet_pred, xgb_pred)
  combined$ECprob_voting <- rowSums(combined[,-c(4:5)])/3
  combined <- combined %>% mutate(voting_pred = case_when((ECprob_voting) >=0.5 ~ "EC",
                                                          (ECprob_voting) < 0.5 ~ "Benign"))
  combined$voting_pred <- factor(combined$voting_pred, levels = c("Benign", "EC"))
  metric <- metric_set(f_meas, precision, recall, specificity, npv)
  metrics_my <- metric(combined, truth = labels, estimate = voting_pred, event_level = "second")
  metrics_my$model <- "ensemble"
  metrics_my$type <- "none"
  
  all <- rbind(rf_results$met, nnet_results$met, xgb_results$met, metrics_my)
  all <- all %>% select(c(.metric, .estimate, model))
  all_wide <- all %>% pivot_wider(names_from = .metric, values_from = .estimate)
  return(list(met = all_wide, pred = combined))
}

chao_microbiome_only <- get_results_microbiome(cohort = "chao")
tsementzi_microbiome_only <- get_results_microbiome(cohort = "tsementzi")
walsh_microbiome_only <- get_results_microbiome(cohort = "walsh")

all <- chao_microbiome_only$pred %>%
  left_join(tsementzi_microbiome_only$pred, by='sraID', copy = FALSE) %>% 
  left_join(walsh_microbiome_only$pred, by = "sraID", copy = FALSE)

all <- all[!endsWith(names(all), '.y')]
all <- all[!endsWith(names(all), '.x')]
all$labels <- factor(sample_data(Antonio1_solo_res$raw)$histology, levels = c("Benign", "EC"))
all <-  all %>% relocate(sraID, .before = labels)      
all <- all %>% select(-c(ECprob_voting, voting_pred))
all <- all[!endsWith(names(all), 'xgb')]

all$ECprob_voting <- rowSums(all[,-c(7:8)])/6

all <- all %>% mutate(voting_pred = case_when((ECprob_voting) >=0.5 ~ "EC",
                                              (ECprob_voting) < 0.5 ~ "Benign"))
all$voting_pred <- factor(all$voting_pred, levels = c("Benign", "EC"))
#predictions$final_predprob <- predictions$ECvotes/4
metric <- metric_set(f_meas, precision, recall, specificity, npv)
metrics_my <- metric(all, truth = labels, estimate = voting_pred, event_level = "second")
metrics_my[nrow(metrics_my) + 1,] <-  roc_auc(all, truth = labels, ECprob_voting, event_level = "second")


