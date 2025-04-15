library(yardstick)
library(themis)
library(tidymodels)
library(mice)

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

prep_test <- function(data, trained_data){
  data <- data.frame(janitor::clean_names(data))
  #data$age <- scale(data$age, center = T, scale = T)[, 1]
  #data$bmi <- scale(data$bmi, center = T, scale = T)[, 1]
  
  return(data)
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

## MICE imputation 
shannon_train <- as.data.frame(microbiome::alpha(Walsh_solo_res$raw, 
                                                 index = c("diversity_shannon")))
train_fitered <- phyloseq::filter_taxa(Walsh_solo_res$raw, function(x) sum(x > 0) > (0.05*length(x)), TRUE)
training_data <- my_tax_glom(train_fitered, level = "genus")$tss
training_data$shannon <- (shannon_train$diversity_shannon)
training_data$age <- (sample_data(Walsh_solo_res$raw)$age)
training_data$bmi <- as.numeric(sample_data(Walsh_solo_res$raw)$bmi)
#training_data$pHRecoded <- factor(sample_data(Walsh_solo_res$raw)$pHRecoded)
training_data$ethnicity <- factor(sample_data(Walsh_solo_res$raw)$ethnicityRecoded)

training_data_tmp <- mice(training_data[, -227],m=1,maxit=50,meth='pmm',seed=123)
training_data_complete <- complete(training_data_tmp,1)
training_data_complete_sub <- training_data_complete %>% select(c(age, bmi, ethnicity))
training_data_complete_sub$histology <- factor(sample_data(Walsh_solo_res$raw)$histology, levels = c("Benign", "EC"))


shannon_test <- as.data.frame(microbiome::alpha(Antonio1_solo_res$raw, 
                                                index = c("diversity_shannon")))
test_fitered <- phyloseq::filter_taxa(Antonio1_solo_res$raw, function(x) sum(x > 0) > (0.05*length(x)), TRUE)
test_data <- my_tax_glom(test_fitered, level = "genus")$tss
test_data$shannon <- (shannon_test$diversity_shannon)
test_data$age <- (sample_data(Antonio1_solo_res$raw)$age)
test_data$bmi <- as.numeric(sample_data(Antonio1_solo_res$raw)$BMI)
#test_data$pHRecoded <- factor(sample_data(Antonio1_solo_res$raw)$pHRecoded)
test_data$ethnicity <- factor("white", levels = c("white"), labels = c("white"))

test_data_tmp <- mice(test_data[, -189],m=1,maxit=50,meth='pmm',seed=123)
test_data_complete <- complete(test_data_tmp,1)
test_data_complete_sub <- test_data_complete %>% select(c(age, bmi, ethnicity))


neigh <- ifelse(min(table(training_data_complete_sub$histology))>25, 5, min(table(training_data_complete_sub$histology))-1)
training_data_complete_sub$histology <- as.factor(training_data_complete_sub$histology)
training_data_complete_sub <- janitor::clean_names(training_data_complete_sub)
ttraining_data_complete_sub_smote <- themis::smotenc(training_data_complete_sub, var = "histology", k = neigh, over_ratio = 1)

set.seed(123)
folds <- rsample::vfold_cv(ttraining_data_complete_sub_smote, strata = "histology", v = 5, repeats = 3)
## remove correlation filter for base models 
data_rec <- recipe(histology ~ ., data = ttraining_data_complete_sub_smote) %>%
  step_normalize(age, bmi)
metric <- metric_set(yardstick::f_meas, yardstick::precision, yardstick::recall, yardstick::specificity, yardstick::roc_auc, yardstick::ppv, yardstick::npv)
ctrl_res <- control_grid(event_level = "second", save_pred = TRUE, save_workflow = TRUE, verbose = TRUE)

## Random forest - set up
rand_forest_spec <- rand_forest(mtry = tune(),
                                min_n = tune(),
                                trees = tune()) %>%
  set_mode("classification") %>%
  set_engine("randomForest", seed = 123)
rand_forest_wflow <- workflow() %>%
  add_recipe(data_rec)%>%   
  add_model(rand_forest_spec) 

library(doMC)
registerDoMC(cores = 5)
## Random length tuning just to get a sense of grid 

rf_grid <- expand.grid(mtry = seq(1, 4, 1), min_n = seq(2, 10, 1), trees = seq(1500, 2500, 100))
rand_forest_res_grid <- tune_grid(object = rand_forest_wflow, 
                                  resamples = folds,  grid = rf_grid,
                                  metrics = metric,
                                  control = ctrl_res) 
final_rf_wf <- rand_forest_wflow %>% 
  finalize_workflow(select_best(rand_forest_res_grid, metric = "f_meas"))
final_fit_rf <- final_rf_wf %>%
  fit(ttraining_data_complete_sub_smote)
rf_results_chao_grid <- get_metrics_indi(data = test_data_complete_sub, model = final_fit_rf, 
                                         level = "genus", labels = sample_data(Antonio1_solo_res$raw)$histology)
saveRDS(rand_forest_res_grid, "~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_all_clinical_no_ph_rf_grid.rds")
saveRDS(final_fit_rf, "~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_all_clinical_no_ph_rf_fit.rds")


## NNET
nnet_spec <- mlp(hidden_units = tune(), penalty = tune(), epochs = tune()) %>%
  set_mode("classification") %>%
  set_engine("nnet", MaxNWts = 10000, seed = 123)
nnet_wflow <- workflow() %>%
  add_recipe(data_rec)%>%
  add_model(nnet_spec) 

##Grid tuning for nnet 
nnet_grid <- expand.grid(hidden_units = seq(1, 30, 1), penalty = seq(0, 1, 0.1), epochs = seq(5, 11, 2))
nnet_res <-tune_grid(object = nnet_wflow, metrics = metric, resamples = folds, 
                     grid = nnet_grid, control = ctrl_res)

set.seed(123)
final_nnet_wf <- nnet_wflow %>% 
  finalize_workflow(select_best(nnet_res, metric = "f_meas" ))
final_fit_nnet <- final_nnet_wf %>%
  fit(ttraining_data_complete_sub_smote) 
nnet_results_chao_grid <- get_metrics_indi(data = test_data_complete_sub, model = final_fit_nnet, 
                                           level = "genus", labels = sample_data(Antonio1_solo_res$raw)$histology)
saveRDS(nnet_res, "~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_all_clinical_no_ph_nnet_grid.rds")
saveRDS(final_fit_nnet, "~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_all_clinical_no_ph_nnet_fit.rds")


## XGBOOST
bt_spec <- boost_tree(mtry = tune(), trees = tune(), min_n = tune(), 
                      tree_depth = tune(), learn_rate = tune(), 
                      loss_reduction = tune(), sample_size = tune(), 
                      stop_iter = tune()) %>%
  set_mode("classification") %>%
  set_engine("xgboost")

xgb_rec <- recipe(histology ~ ., data = ttraining_data_complete_sub_smote) %>%
  step_normalize(age, bmi) %>%
  step_dummy(all_nominal_predictors(), one_hot = TRUE)

bt_wflow <- workflow() %>%
  add_recipe(xgb_rec) %>%
  add_model(bt_spec)

##Grid tuning for xgboost
xgb_grid <- expand.grid(mtry = seq(1, 4, 1), trees = 1500, min_n = seq(2, 10, 1), learn_rate = 0.3, 
                        loss_reduction = seq(0, 1, 0.5), stop_iter = Inf,
                        tree_depth = seq(2, 15, 3), sample_size = c(0.5, 1))
xgb_res <-tune_grid(object = bt_wflow, metrics = metric, resamples = folds, 
                    grid = xgb_grid, control = ctrl_res)

final_bt_wf <- bt_wflow %>% 
  finalize_workflow(select_best(xgb_res, metric = "f_meas" ))
final_fit_bt <- final_bt_wf %>%
  fit(ttraining_data_complete_sub_smote) 
xgb_results_chao_grid <- get_metrics_indi(data = test_data_complete_sub, model = final_fit_bt, 
                                          level = "genus", labels = sample_data(Antonio1_solo_res$raw)$histology)
saveRDS(xgb_res, "~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_all_clinical_no_ph_xgb_grid.rds")
saveRDS(final_fit_bt, "~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_all_clinical_no_ph_xgb_fit.rds")
