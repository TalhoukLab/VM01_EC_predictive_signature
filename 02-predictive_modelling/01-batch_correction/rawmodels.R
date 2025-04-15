library(yardstick)
library(themis)
library(tidymodels)

source(here::here("~/Desktop/thesis/VM02_meta_analysis/src/dataPrep.R"), encoding = "UTF-8")
cohorts <- c("Antonio1")
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
  taxa_names(aggregated) <- paste0(data.frame(tax_table(aggregated))[,idx-3], "_",
                                   data.frame(tax_table(aggregated))[,idx-2], "_",  data.frame(tax_table(aggregated))[,idx-1], "_", 
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

fncols_antonio <- function(data, cname) {
  add <-cname[!cname%in%names(data)]
  if(length(add)!=0) data[add] <- 0
  data
}


prep_test <- function(data, trained_data){
  data <- data.frame(janitor::clean_names(data))
  otu_test <- fncols_antonio(data, setdiff(colnames(trained_data), colnames(data)))
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

get_metrics_indi <- function(data, model, level, labels){
  preds <- prep_test(data = data, trained_data = model$pre$actions$recipe$recipe$template) %>%
    bind_cols(predict(model,
                      ., type = "prob"))%>%
    bind_cols("labels"=labels)
  preds$labels <- factor(preds$labels,  levels = c("Benign", "EC"))
  preds <- preds %>%
    mutate(predicted = case_when(.pred_Benign > .pred_EC ~ "Benign",
                                 .pred_EC > .pred_Benign ~ "EC"))
  preds$predicted <- factor(preds$predicted, levels = c("Benign", "EC"))
  metric <- metric_set(yardstick::f_meas, yardstick::precision, yardstick::recall, yardstick::specificity, yardstick::npv,yardstick::ppv)
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
  metric <- metric_set(yardstick::f_meas, yardstick::precision, yardstick::recall, yardstick::specificity, yardstick::npv,yardstick::ppv)
  metrics_my <- metric(preds, truth = labels, estimate = predicted, event_level = "second", estimator = "binary")
  metrics_my$type <- level
  plot <- autoplot( model, type = "weights")
  res <- list(met = metrics_my, plot = plot, predictions = preds)
  return(res)
}


shannon_test_chao <- as.data.frame(microbiome::alpha(Antonio1_solo_res$raw, 
                                                     index = c("diversity_shannon")))
test_fitered_chao <- phyloseq::filter_taxa(Antonio1_solo_res$raw, function(x) sum(x > 0) > (0.05*length(x)), TRUE)
test_obj_agg_chao <- tax_glom(test_fitered_chao, taxrank = "genus", 
                              NArm = TRUE, bad_empty = c(NA, "", " ", "\t"))
idx <- which(rank_names(test_obj_agg_chao)=="genus")
taxa_names(test_obj_agg_chao) <- paste0(data.frame(tax_table(test_obj_agg_chao))[,idx-3], "_",
                                        data.frame(tax_table(test_obj_agg_chao))[,idx-2], "_", data.frame(tax_table(test_obj_agg_chao))[,idx-1], "_",
                                        data.frame(tax_table(test_obj_agg_chao))[,idx])
test_data_chao <- data.frame(t(data.frame(otu_table(test_obj_agg_chao))))
test_data_chao$shannon <- (shannon_test_chao$diversity_shannon)

raw_data <- phylo_obj_tss
training_raw_data<- data.frame(otu_table(raw_data))
training_raw_data$labels <- factor(sample_data(raw_data)$histology, levels = c("Benign", "EC"))

shannon_train <- as.data.frame(microbiome::alpha(raw_data, 
                                                 index = c("diversity_shannon")))
training_raw_data$shannon <- (shannon_train$diversity_shannon)



neigh <- ifelse(min(table(training_raw_data$labels))>25, 5, min(table(training_raw_data$labels))-1)
training_raw_data$labels <- as.factor(training_raw_data$labels)
training_raw_data <- janitor::clean_names(training_raw_data)
training_raw_data_smote <- themis::smote(training_raw_data, var = "labels", k = neigh, over_ratio = 1)
set.seed(123)
folds <- rsample::vfold_cv(training_raw_data_smote, strata = "labels", v = 5, repeats = 3)
## remove correlation filter for base models 
data_rec <- recipe(labels ~ ., data = training_raw_data_smote) %>%
  step_normalize(shannon)
metric <- metric_set(yardstick::f_meas, yardstick::precision, yardstick::recall, yardstick::specificity, yardstick::roc_auc, yardstick::ppv, yardstick::npv, yardstick::roc_auc)
ctrl_res <- control_grid(event_level = "second", save_pred = TRUE, save_workflow = TRUE, verbose = TRUE)

## Random forest - set up
rand_forest_spec <- rand_forest(mtry = tune(),
                                min_n = tune(),
                                trees = tune()) %>%
  set_mode("classification") %>%
  set_engine("ranger", seed = 123, importance = "permutation", splitrule = "hellinger")
rand_forest_wflow <- workflow() %>%
  add_recipe(data_rec)%>%   
  add_model(rand_forest_spec) 

library(doMC)
registerDoMC(cores = 5)

## Grid tuning for random forest
rf_grid <- expand.grid(mtry = seq(150, 300, 10), min_n = seq(1, 15, 2), trees = seq(1500, 2500, 100))
rand_forest_res_grid <- tune_grid(object = rand_forest_wflow, 
                                  resamples = folds,  grid = rf_grid,
                                  metrics = metric,
                                  control = ctrl_res) 
final_rf_wf <- rand_forest_wflow %>% 
  finalize_workflow(select_best(rand_forest_res_grid, metric = "f_meas"))
final_fit_rf <- final_rf_wf %>%
  fit(training_raw_data_smote)
rf_results_chao_grid <- get_metrics_indi(data = test_data_chao, model = final_fit_rf, level = "genus", labels = data.frame(sample_data(test_obj_agg_chao))$histology)
saveRDS(rand_forest_res_grid, "~/Desktop/thesis/VM01_EC_predictive_signature/02-predictive_modelling/01-batch_correction/raw_rf_grid.rds")
saveRDS(final_fit_rf, "~/Desktop/thesis/VM01_EC_predictive_signature/02-predictive_modelling/01-batch_correction/raw_rf_fit.rds")


nnet_spec <- mlp(hidden_units = tune(), penalty = tune(), epochs = tune()) %>%
  set_mode("classification") %>%
  set_engine("nnet", MaxNWts = 10000)
nnet_wflow <- workflow() %>%
  add_recipe(data_rec)%>%
  add_model(nnet_spec) 

## Grid tuning for nnet
nnet_grid <- expand.grid(hidden_units = seq(10, 300, 10), penalty = seq(0, 1, 0.1), epochs = seq(5, 11, 2))

nnet_res <-tune_grid(object = nnet_wflow, metrics = metric, resamples = folds, 
                     grid = nnet_grid, control = ctrl_res)

final_nnet_wf <- nnet_wflow %>% 
  finalize_workflow(select_best(nnet_res, metric = "f_meas" ))
final_fit_nnet <- final_nnet_wf %>%
  fit(training_raw_data_smote) 
nnet_results_chao_pre <- get_metrics_indi(data = test_data_chao, model = final_fit_nnet, level = "genus", labels = data.frame(sample_data(test_obj_agg_chao))$histology)
saveRDS(nnet_res, "~/Desktop/thesis/VM01_EC_predictive_signature/02-predictive_modelling/01-batch_correction/raw_nnet_grid.rds")
saveRDS(final_fit_nnet, "~/Desktop/thesis/VM01_EC_predictive_signature/02-predictive_modelling/01-batch_correction/raw_nnet_fit.rds")



bt_spec <- boost_tree(mtry = tune(), trees = tune(), min_n = tune(), 
                      tree_depth = tune(), learn_rate = tune(), 
                      loss_reduction = tune(), sample_size = tune(), 
                      stop_iter = tune()) %>%
  set_mode("classification") %>%
  set_engine("xgboost")
bt_wflow <- workflow() %>%
  add_recipe(data_rec)%>%
  add_model(bt_spec)

## Grid tuning for xgb
xgb_grid <- expand.grid(mtry = seq(100, 300, 10), trees = 1500, min_n = seq(2, 10, 3), learn_rate = 0.3, 
                        loss_reduction = seq(0, 1, 0.5), stop_iter = Inf,
                        tree_depth = seq(2, 15, 3), sample_size = c(0.5, 1))

bt_res <-tune_grid(object = bt_wflow, metrics = metric, resamples = folds, 
                   grid = xgb_grid, control = ctrl_res)
final_bt_wf <- bt_wflow %>% 
  finalize_workflow(select_best(bt_res, metric = "f_meas" ))
final_fit_bt <- final_bt_wf %>%
  fit(training_raw_data_smote) 
bt_results_chao_pre <- get_metrics_indi(data = test_data_chao, model = final_fit_bt, level = "genus", labels = data.frame(sample_data(test_obj_agg_chao))$histology)
saveRDS(nnet_res, "~/Desktop/thesis/VM01_EC_predictive_signature/02-predictive_modelling/01-batch_correction/raw_xgb_grid.rds")
saveRDS(final_fit_nnet, "~/Desktop/thesis/VM01_EC_predictive_signature/02-predictive_modelling/01-batch_correction/raw_xgb_fit.rds")
