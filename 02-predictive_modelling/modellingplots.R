library(yardstick)
library(ggplot2)
library(patchwork)

chao_FT <- read.csv("~/Desktop/thesis/VM01_reproducibility_replicability/00-helperfiles/Antonio1FT.csv")


all <- full_join(all, chao_FT, "sraID")
test_all <- all %>% relocate(sraID, .before = predicted1) %>% 
                    relocate(age, .after = sraID) %>% relocate(histology, .after = age) %>% 
                    relocate(labels, .after = histology) %>% relocate(ECprob_voting, .after = labels) %>% 
                    relocate(voting_pred, .after = ECprob_voting)

test_all_post <- test_all %>% filter(menopausal.status == "Post-menopausal" | menopausal.status == "Peri-menopausal")
metric <- metric_set(f_meas, precision, recall, specificity, npv)
metrics_my <- metric(test_all_post, truth = labels, estimate = voting_pred, event_level = "second")

test_all_pre <- test_all %>% filter(menopausal.status == "Pre-menopausal")
metric <- metric_set(f_meas, precision, recall, specificity, npv)
metrics_my <- metric(test_all_pre, truth = labels, estimate = voting_pred, event_level = "second")


## Null model
all_combined$null_pred <- factor(c("EC"), levels = c("Benign", "EC"))
all_combined$ECprob_null <- 0.5
metric <- metric_set(f_meas, precision, recall, specificity, npv)
metrics_my_null <- metric(all_combined, truth = labels, estimate = null_pred, event_level = "second")
metrics_my_null[nrow(metrics_my_null) + 1,] <-  roc_auc(all_combined, truth = labels, ECprob_null, event_level = "second")
null_model_cm <- conf_mat(all_combined, truth = labels, estimate = null_pred)
null_model_cm_plot <- autoplot(null_model_cm, type = "heatmap") +
  scale_fill_gradient(low="white",high = "skyblue") + ggtitle("Null model") + theme(plot.title = element_text(hjust = 0.5))


## clinical model 
metrics_my_voting <- metric(all_combined, truth = labels, estimate = clinical_class, event_level = "second")
metrics_my_voting[nrow(metrics_my_voting) + 1,] <-  roc_auc(all_combined, truth = labels, clinical_prob, event_level = "second")
voting_model_cm <- conf_mat(all_combined, truth = labels, estimate = clinical_class)
voting_model_cm_plot <- autoplot(voting_model_cm, type = "heatmap") +
  scale_fill_gradient(low="white",high = "skyblue") + ggtitle("Clinical model") + theme(plot.title = element_text(hjust = 0.5))

## microbiome model 
metrics_my_voting <- metric(all_combined, truth = labels, estimate = micro_class, event_level = "second")
metrics_my_voting[nrow(metrics_my_voting) + 1,] <-  roc_auc(all_combined, truth = labels, micro_prob, event_level = "second")
voting_model_cm <- conf_mat(all_combined, truth = labels, estimate = micro_class)
voting_model_cm_plot1 <- autoplot(voting_model_cm, type = "heatmap") +
  scale_fill_gradient(low="white",high = "skyblue") + ggtitle("Microbiome model") + theme(plot.title = element_text(hjust = 0.5))

## all 
metrics_my_voting <- metric(all_combined, truth = labels, estimate = clinical_micro_class, event_level = "second")
metrics_my_voting[nrow(metrics_my_voting) + 1,] <-  roc_auc(all_combined, truth = labels, clinial_micro_prob, event_level = "second")
voting_model_cm <- conf_mat(all_combined, truth = labels, estimate = clinical_micro_class)
voting_model_cm_plot2 <- autoplot(voting_model_cm, type = "heatmap") +
  scale_fill_gradient(low="white",high = "skyblue") + ggtitle("Clinical + Microbiome model") + theme(plot.title = element_text(hjust = 0.5))

library(patchwork)
null_model_cm_plot | voting_model_cm_plot | voting_model_cm_plot1 | voting_model_cm_plot2

##combaseq 
combatseq_model <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/results/combat_rf_sat_model.rds")
combatseq_results <- get_metrics_indi_all(data = test_data_chao, model = combatseq_model, 
                                    level = "genus", labels = sample_data(test_obj_agg_chao)$histology)
predictions_all$combat_pred <- combatseq_results$predictions$predicted
predictions_all$EC_prob_combat <- combatseq_results$predictions$.pred_EC
combatseq_model_cm <- conf_mat(predictions_all, truth = labels, estimate = combat_pred)
combatseq_model_cm_plot <- autoplot(combatseq_model_cm, type = "heatmap") +
  scale_fill_gradient(low="white",high = "skyblue") + ggtitle("Batch corrected") + theme(plot.title = element_text(hjust = 0.5))

## Raw model 
raw_model <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/results/raw_rf_model.rds")
raw_results <- get_metrics_indi_all(data = test_data_chao, model = raw_model, 
                                                  level = "genus", labels = sample_data(test_obj_agg_chao)$histology)
predictions_all$raw_pred <- raw_results$predictions$predicted
predictions_all$EC_prob_raw <- raw_results$predictions$.pred_EC
raw_model_cm <- conf_mat(predictions_all, truth = labels, estimate = raw_pred)
raw_model_cm_plot <- autoplot(raw_model_cm, type = "heatmap") +
  scale_fill_gradient(low="white",high = "skyblue") + ggtitle("Non-batch corrected") + theme(plot.title = element_text(hjust = 0.5))


library(patchwork)
null_model_cm_plot  | raw_model_cm_plot | combatseq_model_cm_plot |voting_model_cm_plot 
