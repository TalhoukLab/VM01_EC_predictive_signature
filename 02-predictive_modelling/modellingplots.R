library(yardstick)
library(ggplot2)
library(patchwork)

Antonio_FT <- read.csv("00-helperfiles/Antonio1FT.csv")
raw_models <- read.csv("VM01_EC_predictive_signature/02-predictive_modelling/raw.csv")
combat_models <- read.csv("VM01_EC_predictive_signature/02-predictive_modelling/combat.csv")
microbiome_only_models <- read.csv("VM01_EC_predictive_signature/02-predictive_modelling/microbiome_only.csv")
clinical_with_ph <- read.csv("M01_EC_predictive_signature/02-predictive_modelling/clinical_with_ph.csv")
clinical_without_ph <- read.csv("VM01_EC_predictive_signature/02-predictive_modelling/clinical_without_ph.csv")
clinical_with_ph_micro <- read.csv("VM01_EC_predictive_signature/02-predictive_modelling/clinical_with_ph_micro.csv")
clinical_without_ph_micro <- read.csv("VM01_EC_predictive_signature/02-predictive_modelling/clinical_no_ph_micro.csv")


## raw 
all_results <- raw_models %>% 
                    select(raw_rf, sraID) %>% 
                    dplyr::rename("raw_prob" = raw_rf) %>% 
                    dplyr::mutate("raw_pred" = ifelse(raw_prob >= 0.5, "EC", "Benign")) %>%
                    full_join(Antonio_FT, "sraID") 
all_results$labels <- as.factor(all_results$histology)
all_results$raw_pred <- factor(all_results$raw_pred, levels = c("Benign", "EC"))

## combat batch corrected 
all_results <- combat_models %>%
                select(combat_nnet, sraID) %>%
                dplyr::rename("combat_prob" = combat_nnet) %>%
                dplyr::mutate("combat_pred" = ifelse(combat_prob >= 0.5, "EC", "Benign")) %>%
                full_join(all_results, "sraID") 
all_results$combat_pred <- factor(all_results$combat_pred, levels = c("Benign", "EC"))

## Microbiome 
all_results <- microbiome_only_models %>%
  select(ensemble, ensemble_pred, sraID) %>%
  dplyr::rename("microbiome_only_prob" = ensemble,
                "microbiome_only_pred" = ensemble_pred) %>%
  full_join(all_results, "sraID")
all_results$microbiome_only_pred <- factor(all_results$microbiome_only_pred, levels = c("Benign", "EC"))

## clinical with pH
all_results <- clinical_with_ph %>%
  select(ensemble, ensemble_pred, sraID) %>%
  dplyr::rename("clinical_with_ph_prob" = ensemble,
                "clinical_with_ph_pred" = ensemble_pred) %>%
  full_join(all_results, "sraID")
all_results$clinical_with_ph_pred <- factor(all_results$clinical_with_ph_pred, levels = c("Benign", "EC"))

## clinical without pH
all_results <- clinical_without_ph %>%
  select(ensemble, ensemble_pred, sraID) %>%
  dplyr::rename("clinical_without_ph_prob" = ensemble,
                "clinical_without_ph_pred" = ensemble_pred) %>%
  full_join(all_results, "sraID")
all_results$clinical_without_ph_pred <- factor(all_results$clinical_without_ph_pred, levels = c("Benign", "EC"))

## Microbiome + clinical with ph
all_results <- clinical_with_ph_micro %>%
                select(ensemble, ensemble_pred, sraID) %>%
                dplyr::rename("micro_clinical_with_ph_prob" = ensemble,
                              "micro_clinical_with_ph_pred" = ensemble_pred) %>%
                full_join(all_results, "sraID")
all_results$micro_clinical_with_ph_pred <- factor(all_results$micro_clinical_with_ph_pred, levels = c("Benign", "EC"))


## Microbiome + clinical without ph
all_results <- clinical_without_ph_micro %>%
  select(ensemble, ensemble_pred, sraID) %>%
  dplyr::rename("micro_clinical_without_ph_prob" = ensemble,
                "micro_clinical_without_ph_pred" = ensemble_pred) %>%
  full_join(all_results, "sraID")
all_results$micro_clinical_without_ph_pred <- factor(all_results$micro_clinical_without_ph_pred, levels = c("Benign", "EC"))


all_results <- all_results %>% relocate(sraID, .before = microbiome_only_prob) %>% 
                    relocate(age, .after = sraID) %>% relocate(histology, .after = age) %>% 
                    relocate(labels, .after = histology) 




## Null model
all_results$null_pred <- factor(c("EC"), levels = c("Benign", "EC"))
all_results$null_prob <- 0.5
metric <- metric_set(f_meas, precision, recall, specificity, npv)
metrics_my_null <- metric(all_results, truth = labels, estimate = null_pred, event_level = "second")
metrics_my_null[nrow(metrics_my_null) + 1,] <-  roc_auc(all_results, truth = labels, null_prob, event_level = "second")

null_model_cm_table <- table("actual" = all_results$labels, "predicted" = all_results$null_pred)
null_model_cm_tibble <- as_tibble(null_model_cm_table)
null_model_cm_tibble$refere <- factor(null_model_cm_tibble$actual, levels = c("EC", "Benign"))
#null_model_cm_tibble$predicted <- factor(null_model_cm_tibble$predicted, levels = c("EC", "Benign"))
null_plot <- ggplot(data =  null_model_cm_tibble, mapping = aes(x = actual, y = predicted)) +
  geom_tile(aes(fill = n)) +
  geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1, size = 15) +
  scale_fill_gradient(low = "white", high = "deepskyblue2") +
  theme_bw() + theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 30), axis.text = element_text(size = 30))


## Raw model
metrics_raw <- metric(all_results, truth = labels, estimate = raw_pred, event_level = "second")
metrics_raw[nrow(metrics_raw) + 1,] <-  roc_auc(all_results, truth = labels, raw_prob, event_level = "second")
raw_model_cm <- conf_mat(all_results, truth = labels, estimate = raw_pred)
raw_model_cm_table <- table("actual" = all_results$labels, "predicted" = all_results$raw_pred)
raw_model_cm_tibble <- as_tibble(raw_model_cm_table)
raw_model_cm_tibble$actual <- factor(raw_model_cm_tibble$actual, levels = c("EC", "Benign"))
#raw_model_cm_tibble$predicted <- factor(raw_model_cm_tibble$predicted, levels = c("EC", "Benign"))
raw_plot <- ggplot(data =  raw_model_cm_tibble, mapping = aes(x = actual, y = predicted)) +
  geom_tile(aes(fill = n)) +
  geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1, size = 15) +
  scale_fill_gradient(low = "white", high = "deepskyblue2") +
  theme_bw() + theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 30), axis.text = element_text(size = 30))

## combat model
metrics_combat <- metric(all_results, truth = labels, estimate = combat_pred, event_level = "second")
metrics_combat[nrow(metrics_combat) + 1,] <-  roc_auc(all_results, truth = labels, combat_prob, event_level = "second")
combat_model_cm_table <- table("actual" = all_results$labels, "predicted" = all_results$combat_pred)
combat_model_cm_tibble <- as_tibble(combat_model_cm_table)
combat_model_cm_tibble$actual <- factor(combat_model_cm_tibble$actual, levels = c("EC", "Benign"))
#combat_model_cm_tibble$predicted <- factor(combat_model_cm_tibble$predicted, levels = c("EC", "Benign"))
combat_plot <- ggplot(data =  combat_model_cm_tibble, mapping = aes(x = actual, y = predicted)) +
  geom_tile(aes(fill = n)) +
  geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1, size = 15) +
  scale_fill_gradient(low = "white", high = "deepskyblue2") +
  theme_bw() + theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 30), axis.text = element_text(size = 30))

## microbiome model 
metrics_microbiome <- metric(all_results, truth = labels, estimate = microbiome_only_pred, event_level = "second")
metrics_microbiome[nrow(metrics_microbiome) + 1,] <-  roc_auc(all_results, truth = labels, microbiome_only_prob, event_level = "second")
microbiome_model_cm_table <- table("actual" = all_results$labels, "predicted" = all_results$microbiome_only_pred)
microbiome_model_cm_tibble <- as_tibble(microbiome_model_cm_table)
microbiome_model_cm_tibble$actual <- factor(microbiome_model_cm_tibble$actual, levels = c("EC", "Benign"))
#microbiome_model_cm_tibble$predicted <- factor(microbiome_model_cm_tibble$predicted, levels = c("EC", "Benign"))
microbiome_only_plot <- ggplot(data =  microbiome_model_cm_tibble, mapping = aes(x = actual, y = predicted)) +
  geom_tile(aes(fill = n)) +
  geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1, size = 15) +
  scale_fill_gradient(low = "white", high = "deepskyblue2") +
  theme_bw() + theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 30), axis.text = element_text(size = 30))


## clinical with ph model 
metrics_clinical <- metric(all_results, truth = labels, estimate = clinical_with_ph_pred, event_level = "second")
metrics_clinical[nrow(metrics_clinical) + 1,] <-  roc_auc(all_results, truth = labels, clinical_with_ph_prob, event_level = "second")
clinical_model_cm_table <- table("actual" = all_results$labels, "predicted" = all_results$clinical_with_ph_pred)
clinical_model_cm_tibble <- as_tibble(clinical_model_cm_table)
clinical_model_cm_tibble$actual <- factor(clinical_model_cm_tibble$actual, levels = c("EC", "Benign"))
#clinical_model_cm_tibble$predicted <- factor(clinical_model_cm_tibble$predicted, levels = c("EC", "Benign"))
clinical_with_ph_plot <- ggplot(data =  clinical_model_cm_tibble, mapping = aes(x = actual, y = predicted)) +
  geom_tile(aes(fill = n)) +
  geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1, size = 15) +
  scale_fill_gradient(low = "white", high = "deepskyblue2") +
  theme_bw() + theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 30), axis.text = element_text(size = 30))

## clinical without ph model 
metrics_clinical_no_ph <- metric(all_results, truth = labels, estimate = clinical_without_ph_pred, event_level = "second")
metrics_clinical_no_ph[nrow(metrics_clinical_no_ph) + 1,] <-  roc_auc(all_results, truth = labels, clinical_without_ph_prob, event_level = "second")
clinical_no_ph_model_cm_table <- table("actual" = all_results$labels, "predicted" = all_results$clinical_without_ph_pred)
clinical_no_ph_model_cm_table <- as_tibble(clinical_no_ph_model_cm_table)
clinical_no_ph_model_cm_table$actual <- factor(clinical_no_ph_model_cm_table$actual, levels = c("EC", "Benign"))
#clinical_model_cm_tibble$predicted <- factor(clinical_model_cm_tibble$predicted, levels = c("EC", "Benign"))
clinical_no_ph_plot <- ggplot(data =  clinical_no_ph_model_cm_table, mapping = aes(x = actual, y = predicted)) +
  geom_tile(aes(fill = n)) +
  geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1, size = 15) +
  scale_fill_gradient(low = "white", high = "deepskyblue2") +
  theme_bw() + theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 30), axis.text = element_text(size = 30))



## all 
metrics_micro_clinical <- metric(all_results, truth = labels, estimate = micro_clinical_with_ph_pred, event_level = "second")
metrics_micro_clinical[nrow(metrics_micro_clinical) + 1,] <-  roc_auc(all_results, truth = labels, micro_clinical_with_ph_prob, event_level = "second")

micro_clinical_model_cm_table <- table("actual" = all_results$labels, "predicted" = all_results$micro_clinical_with_ph_pred)
micro_clinical_model_cm_tibble <- as_tibble(micro_clinical_model_cm_table)
micro_clinical_model_cm_tibble$actual <- factor(micro_clinical_model_cm_tibble$actual, levels = c("EC", "Benign"))
#micro_clinical_model_cm_tibble$predicted <- factor(micro_clinical_model_cm_tibble$predicted, levels = c("EC", "Benign"))
micro_clinical_with_ph_plot <- ggplot(data =  micro_clinical_model_cm_tibble, mapping = aes(x = actual, y = predicted)) +
  geom_tile(aes(fill = n)) +
  geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1, size = 15) +
  scale_fill_gradient(low = "white", high = "deepskyblue2") +
  theme_bw() + theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 30), axis.text = element_text(size = 30))

# micro + clinical without pH
metrics_micro_clinical_no_ph <- metric(all_results, truth = labels, estimate = micro_clinical_without_ph_pred, event_level = "second")
metrics_micro_clinical_no_ph[nrow(metrics_micro_clinical_no_ph) + 1,] <-  roc_auc(all_results, truth = labels, micro_clinical_without_ph_prob, event_level = "second")

micro_clinical_no_ph_model_cm_table <- table("actual" = all_results$labels, "predicted" = all_results$micro_clinical_without_ph_pred)
micro_clinical_no_ph_model_cm_tibble <- as_tibble(micro_clinical_no_ph_model_cm_table)
micro_clinical_no_ph_model_cm_tibble$actual <- factor(micro_clinical_no_ph_model_cm_tibble$actual, levels = c("EC", "Benign"))
#micro_clinical_model_cm_tibble$predicted <- factor(micro_clinical_model_cm_tibble$predicted, levels = c("EC", "Benign"))
micro_clinical_without_ph_plot <- ggplot(data =  micro_clinical_no_ph_model_cm_tibble, mapping = aes(x = actual, y = predicted)) +
  geom_tile(aes(fill = n)) +
  geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1, size = 15) +
  scale_fill_gradient(low = "white", high = "deepskyblue2") +
  theme_bw() + theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 30), axis.text = element_text(size = 30))



library(patchwork)
pdf("~/Desktop/confu_mat1.pdf", width =40, height = 6)
null_model_cm_plot  | raw_model_cm_plot | combat_model_cm_plot | clinical_model_cm_plot | microbiome_model_cm_plot | micro_clinical_model_cm_plot 
dev.off()

combined_plot <- ggarrange(null_plot, 
                           raw_plot,combat_plot, 
                           microbiome_only_plot, clinical_no_ph_plot,
                           micro_clinical_without_ph_plot,
                          nrow = 1,
                          ncol = 8) 
combined_plot
dev.off()
## plot auc 

plot(pROC::smooth(pROC::roc(all_results$labels, all_results$raw_prob)), col = "#66c2a5", legacy.axes = T)
plot((pROC::roc(all_results$labels, all_results$combat_prob)), col = "#fc8d62", legacy.axes = T, add = T)
plot(pROC::smooth(pROC::roc(all_results$labels, all_results$clinical_prob)), col = "#8da0cb", legacy.axes = T, add = T)
plot(pROC::smooth(pROC::roc(all_results$labels, all_results$microbiome_prob)), col = "#e78ac3", legacy.axes = T, add = T)
plot(pROC::smooth(pROC::roc(all_results$labels, all_results$micro_clinical_prob)), col = "#a6d854", legacy.axes = T, add = T)
legend("bottomright",
       legend=c("Raw model", "Combat model", "Clinical model", "Microbiome model", "Clinical + microbiome model"),
       col=c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854"),
       lty=1,
       lwd=c(2,2))


## generate CI 
library(epiR)
## tp, fp, fn, tn
data_raw <- as.table(matrix(c(12,10,0,0), nrow = 2, byrow = TRUE))
rval_raw <- epi.tests(data_raw, conf.level = 0.95)
rval_raw$detail$model <- "raw"
raw_auc <- pROC::auc(all_results$labels, all_results$raw_prob)
raw_auc_ci <- pROC::ci.auc(raw_auc)
rval_raw$detail[nrow(rval_raw$detail) + 1,] <-  c("ROC_AUC",  raw_auc, raw_auc_ci[1], raw_auc_ci[2], "raw")


data_combat <- as.table(matrix(c(12,10,0,0), nrow = 2, byrow = TRUE))
rval_combat <- epi.tests(data_combat, conf.level = 0.95)
rval_combat$detail$model <- "combat"
combat_auc <- pROC::auc(all_results$labels, all_results$combat_prob)
combat_auc_ci <- pROC::ci.auc(combat_auc)
rval_combat$detail[nrow(rval_combat$detail) + 1,] <-  c("ROC_AUC",  combat_auc, combat_auc_ci[1], combat_auc_ci[2], "combat")

data_clinical <- as.table(matrix(c(7, 1, 5, 9), nrow = 2, byrow=T))
rval_clinical <- epi.tests(data_clinical, conf.level = 0.95)
rval_clinical$detail$model <- "clinical"
clinical_auc <- pROC::auc(all_results$labels, all_results$clinical_with_ph_prob)
clinical_auc_ci <- pROC::ci.auc(clinical_auc)
rval_clinical$detail[nrow(rval_clinical$detail) + 1,] <-  c("ROC_AUC",  clinical_auc, clinical_auc_ci[1], clinical_auc_ci[2], "clinical")

data_clinical_no_ph <- as.table(matrix(c(8, 1, 4, 9), nrow = 2, byrow=T))
rval_clinical_no_ph <- epi.tests(data_clinical_no_ph, conf.level = 0.95)
rval_clinical_no_ph$detail$model <- "clinical no ph"
clinical_no_ph_auc <- pROC::auc(all_results$labels, all_results$clinical_without_ph_prob)
clinical_no_ph_auc_ci <- pROC::ci.auc(clinical_no_ph_auc)
rval_clinical_no_ph$detail[nrow(rval_clinical_no_ph$detail) + 1,] <-  c("ROC_AUC",  clinical_no_ph_auc, clinical_no_ph_auc_ci[1], clinical_no_ph_auc_ci[2], "clinical no ph")

data_microbiome <- as.table(matrix(c(8, 3, 4, 7), nrow = 2, byrow=T))
rval_microbiome <- epi.tests(data_microbiome, conf.level = 0.95)
rval_microbiome$detail$model <- "microbiome"
microbiome_auc <- pROC::auc(all_results$labels, all_results$microbiome_only_prob)
microbiome_auc_ci <- pROC::ci.auc(microbiome_auc)
rval_microbiome$detail[nrow(rval_microbiome$detail) + 1,] <-  c("ROC_AUC",  microbiome_auc, microbiome_auc_ci[1], microbiome_auc_ci[2], "microbiome")

data_micro_clinical <- as.table(matrix(c(12, 4, 0, 6), nrow = 2, byrow=T))
rval_micro_clinical <- epi.tests(data_micro_clinical, conf.level = 0.95)
rval_micro_clinical$detail$model <- "clinical_microbiome"
micro_clinical_auc <- pROC::auc(all_results$labels, all_results$micro_clinical_with_ph_prob)
micro_clinical_auc_ci <- pROC::ci.auc(micro_clinical_auc)
rval_micro_clinical$detail[nrow(rval_micro_clinical$detail) + 1,] <-  c("ROC_AUC",  micro_clinical_auc, micro_clinical_auc_ci[1], micro_clinical_auc_ci[2], "clinical_microbiome")

data_micro_clinical_no_ph <- as.table(matrix(c(12, 3, 0, 7), nrow = 2, byrow=T))
rval_micro_clinical_no_ph <- epi.tests(data_micro_clinical_no_ph, conf.level = 0.95)
rval_micro_clinical_no_ph$detail$model <- "clinical_no_ph_microbiome"
micro_clinical_no_ph_auc <- pROC::auc(all_results$labels, all_results$micro_clinical_without_ph_prob)
micro_clinical_no_ph_auc_ci <- pROC::ci.auc(micro_clinical_no_ph_auc)
rval_micro_clinical_no_ph$detail[nrow(rval_micro_clinical_no_ph$detail) + 1,] <-  c("ROC_AUC",  micro_clinical_no_ph_auc, micro_clinical_auc_ci[1], micro_clinical_no_ph_auc_ci[2], "clinical_no_ph_microbiome")


all_results_st <- rbind(rval_raw$detail, rval_combat$detail, 
                        rval_clinical$detail, rval_clinical_no_ph$detail,
                        rval_microbiome$detail, rval_micro_clinical$detail, rval_micro_clinical_no_ph$detail)
all_results_st_fil <- all_results_st %>% filter(statistic == "se" |
                                                  statistic == "sp" | statistic == "pv.pos" |
                                                  statistic == "pv.neg" | statistic == "ROC_AUC")
all_results_st_fil$statistic[all_results_st_fil$statistic == "se" ] <- "sensitivity"
all_results_st_fil$statistic[all_results_st_fil$statistic == "sp" ] <- "specificity"
all_results_st_fil$statistic[all_results_st_fil$statistic == "pv.pos" ] <- "PPV"
all_results_st_fil$statistic[all_results_st_fil$statistic == "pv.neg" ] <- "NPV"
all_results_st_fil$est <- round(as.numeric(all_results_st_fil$est), 3)
all_results_st_fil$upper <- round(as.numeric(all_results_st_fil$upper), 3)
all_results_st_fil$lower <- round(as.numeric(all_results_st_fil$lower), 3)
all_results_st_fil$meas <- paste0(all_results_st_fil$est, " [", all_results_st_fil$lower, ", ", all_results_st_fil$upper, "]")
all_results_st_fil$metric <- all_results_st_fil$statistic

all_results_st_fil_sel <- all_results_st_fil %>% select(metric, model, meas)
all_results_st_fil_wide <- pivot_wider(all_results_st_fil_sel, names_from = model, values_from = meas)
rempsyc::nice_table(all_results_st_fil_wide)
print(all_results_st_fil_wide, preview = "docx")
flextable::save_as_docx(all_results_st_fil_wide, path = "~/Desktop/test.docx")
