#Variable importance 

## Micorbiome models 
#'chao_rf', 'tsementzi_rf',  'tsementzi_xgb', 'walsh_rf', 'walsh_xgb', 'gressel_nnet'
microbiome_chao_rf_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Chao_microbiome_rf_fit.rds")
microbiome_chao_rf_fit$name <- "Chao - RF"
microbiome_tsementzi_xgb_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Tsementzi_microbiome_xgb_fit.rds")
microbiome_tsementzi_xgb_fit$name <- "Tsementzi - XGB"
microbiome_walsh_rf_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Walsh_microbiome_rf_fit.rds")
microbiome_walsh_rf_fit$name <- "Walsh - RF"
microbiome_walsh_xgb_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Walsh_microbiome_xgb_fit.rds")
microbiome_walsh_xgb_fit$name <- "Walsh - XGB"
microbiome_gressel_nnet_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Gressel_microbiome_nnet_fit.rds")
microbiome_gressel_nnet_fit$name <- "Gressel - NNET"

microbiome_models <- list("Walsh - RF" = microbiome_walsh_rf_fit, "Walsh - XGB" = microbiome_walsh_xgb_fit, 
                          "Tsementzi - XGB" = microbiome_tsementzi_xgb_fit,
                          "Gressel - NNET" = microbiome_gressel_nnet_fit,  "Chao - RF" = microbiome_chao_rf_fit)
t <- lapply(microbiome_models, function(x) vip(x$fit$fit, num_features = 10L) + 
              theme(axis.title = element_text(size = 30), plot.title = element_text(size = 40),
                    axis.text = element_text(size = 30))+  ggtitle(x$name))

names(t) <- names(microbiome_models)
lapply(names(t), 
       function(x) ggsave(filename=paste("~/Desktop/paper1/VIP/", x, "microbiome_only.jpeg"), plot=t[[x]],  width = 20, height = 10))


## Clinical models 
#'chao_xgb', 'tsementzi_xgb', 'walsh_rf', 'walsh_nnet'
#clinical_chao_xgb_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Chao_age_xgb_fit.rds")
clinical_tsementzi_xgb_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Tsementzi_all_clinical_xgb_fit.rds")
clinical_tsementzi_xgb_fit$name <- "Tsementzi - XGB"
clinical_walsh_rf_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_all_clinical_rf_fit.rds")
clinical_walsh_rf_fit$name <- "Walsh - RF"
clinical_walsh_nnet_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_all_clinical_nnet_fit.rds")
clinical_walsh_nnet_fit$name <- "Walsh - NNET"

clinical_models <- list("Tsementzi - XGB" =clinical_tsementzi_xgb_fit,
                        "Walsh - RF" = clinical_walsh_rf_fit, "Walsh - NNET" =clinical_walsh_nnet_fit)
t1 <- lapply(clinical_models, function(x) vip(x$fit$fit, num_features = 10L) + 
              theme(axis.title = element_text(size = 12),
                    axis.text = element_text(size = 12)) +  ggtitle(x$name))
names(t1) <- names(clinical_models)
lapply(names(t1), 
       function(x) ggsave(filename=paste("~/Desktop/paper1/VIP/", x, "clinical_with_ph.jpeg"), plot=t1[[x]]))

## Clinical models no pH
#'chao_xgb', 'tsementzi_xgb', 'walsh_rf', 'walsh_nnet'
#clinical_chao_xgb_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Chao_age_xgb_fit.rds")
clinical_no_ph_tsementzi_xgb_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Tsementzi_all_clinical_no_ph_xgb_fit.rds")
clinical_no_ph_tsementzi_xgb_fit$name <- "Tsementzi - XGB"
clinical_no_ph_walsh_rf_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_all_clinical_no_ph_rf_fit.rds")
clinical_no_ph_walsh_rf_fit$name <- "Walsh - RF"
clinical_no_ph_walsh_nnet_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_all_clinical_no_ph_nnet_fit.rds")
clinical_no_ph_walsh_nnet_fit$name <- "Walsh -NNET"

clinical_no_ph_models <- list( "Tsementzi - XGB" = clinical_no_ph_tsementzi_xgb_fit,
                               "Walsh - RF" = clinical_no_ph_walsh_rf_fit, 
                               "Walsh -NNET" =clinical_no_ph_walsh_nnet_fit)
t2 <- lapply(clinical_no_ph_models, function(x) vip(x$fit$fit, num_features = 10L) + 
               theme(axis.title = element_text(size = 12),
                     axis.text = element_text(size = 12))+  ggtitle(x$name))
names(t2) <- names(clinical_no_ph_models)
lapply(names(t2), 
       function(x) ggsave(filename=paste("~/Desktop/paper1/VIP/", x, "clinical_without_ph.jpeg"), plot=t2[[x]]))


## Microbiome and clinical without pH
#'chao_xgb', 'chao_rf', 'tsementzi_rf',  'walsh_rf', 'walsh_xgb', 'gressel_nnet'

## Chao 
chao_rf_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Chao_clinical_micro_rf_fit.rds")
chao_rf_fit$name <- "Chao - RF"
chao_xgb_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Chao_clinical_micro_xgb_fit.rds")
chao_xgb_fit$name <- "Chao - XGB"

## Gressel 
gressel_nnet_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Gressel_microbiome_nnet_fit.rds")
gressel_nnet_fit$name <- "Gressel - NNET"

## Tsementzi 
tsementzi_rf_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Tsementzi_clinical_no_ph_micro_rf_fit.rds")
tsementzi_rf_fit$name <- "Tsementzi - RF"

## Walsh
walsh_rf_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_clinical_no_phmicro_rf_fit.rds")
walsh_rf_fit$name <- "Walsh - RF"
#walsh_xgb_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_clinical_no_ph.rds")
#walsh_xgb_fit$name <- "Walsh - XGB"

clinical_mciro_models <- list("Walsh - RF" = walsh_rf_fit, "Walsh - XGB" = walsh_xgb_fit,
                              "Tsementzi - RF" = tsementzi_rf_fit, "Gressel - NNET" = gressel_nnet_fit, 
                              "Chao - RF" = chao_rf_fit, "Chao - XGB" = chao_xgb_fit)
t3 <- lapply(clinical_mciro_models, function(x) vip(x$fit$fit, num_features = 10L) +
               theme(axis.title = element_text(size = 30), 
                     plot.title = element_text(size = 40),axis.text = element_text(size = 30))+  ggtitle(x$name))  
names(t3) <- names(clinical_mciro_models)
lapply(names(t3), 
       function(x) ggsave(filename=paste("~/Desktop/paper1/VIP/", x, "clinical_micro_with_ph.jpeg"), plot=t3[[x]], width = 20, height = 10))


## Microbiome and clinical 
#'chao_xgb', 'chao_rf', 'tsementzi_rf',  'walsh_rf', 'walsh_xgb', 'gressel_nnet'

## Chao 
chao_rf_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Chao_clinical_micro_rf_fit.rds")
chao_rf_fit$name <- "Chao - RF"
chao_xgb_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Chao_clinical_micro_xgb_fit.rds")
chao_xgb_fit$name <- "Chao - XGB"

## Gressel 
gressel_nnet_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Gressel_microbiome_nnet_fit.rds")
gressel_nnet_fit$name <- "Gressel - NNET"

## Tsementzi 
tsementzi_rf_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Tsementzi_clinical_micro_rf_fit.rds")
tsementzi_rf_fit$name <- "Tsementzi - RF"

## Walsh
walsh_rf_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_clinical_micro_rf_fit.rds")
walsh_rf_fit$name <- "Walsh - RF"
walsh_xgb_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_clinical_micro_xgb_fit.rds")
walsh_xgb_fit$name <- "Walsh - XGB"

clinical_mciro_models <- list("Walsh - RF" = walsh_rf_fit, "Walsh - XGB" = walsh_xgb_fit,
                              "Tsementzi - RF" = tsementzi_rf_fit, "Gressel - NNET" = gressel_nnet_fit, 
                              "Chao - RF" = chao_rf_fit, "Chao - XGB" = chao_xgb_fit)
t3 <- lapply(clinical_mciro_models, function(x) vip(x$fit$fit, num_features = 10L) +
               theme(axis.title = element_text(size = 30), 
                     plot.title = element_text(size = 40),axis.text = element_text(size = 30))+  ggtitle(x$name))  
names(t3) <- names(clinical_mciro_models)
lapply(names(t3), 
       function(x) ggsave(filename=paste("~/Desktop/paper1/VIP/", x, "clinical_micro_with_ph.jpeg"), plot=t3[[x]], width = 20, height = 10))

## Microbiome and clinical with no pH 
#'chao_xgb', 'chao_rf', 'tsementzi_rf',  'walsh_rf', 'walsh_xgb', 'gressel_nnet'

## Chao 
chao_rf_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Chao_clinical_micro_rf_fit.rds")
chao_rf_fit$name <- "Chao - RF"
chao_xgb_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Chao_clinical_micro_xgb_fit.rds")
chao_xgb_fit$name <- "Chao - XGB"

## Gressel 
gressel_nnet_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Gressel_microbiome_nnet_fit.rds")
gressel_nnet_fit$name <- "Gressel - NNET"

## Tsementzi 
tsementzi_rf_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/Tsementzi_clinical_no_ph_micro_rf_fit.rds")
tsementzi_rf_fit$name <- "Tsementzi - RF"

## Walsh
walsh_rf_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_clinical_no_phmicro_rf_fit.rds")
walsh_rf_fit$name <- "Walsh - RF"
walsh_xgb_fit <- readRDS("~/Desktop/thesis/VM02_meta_analysis/src/updated/updated1/Walsh_clinical_no_ph_micro_xgb_fit.rds")
walsh_xgb_fit$name <- "Walsh - XGB"

clinical_no_ph_mciro_models <- list("Walsh - RF" = walsh_rf_fit, "Walsh - XGB" = walsh_xgb_fit,
                              "Tsementzi - RF" = tsementzi_rf_fit, "Gressel - NNET" = gressel_nnet_fit, 
                              "Chao - RF" = chao_rf_fit, "Chao - XGB" = chao_xgb_fit)
t4 <- lapply(clinical_no_ph_mciro_models, function(x) vip(x$fit$fit, num_features = 10L) +
                     theme(axis.title = element_text(size = 30), 
                           plot.title = element_text(size = 40),axis.text = element_text(size = 30))+  ggtitle(x$name))  
names(t4) <- names(clinical_no_ph_mciro_models)
lapply(names(t4), 
       function(x) ggsave(filename=paste("~/Desktop/paper1/VIP/", x, "clinical_no_ph_micro.jpeg"), plot=t4[[x]], width = 20, height = 10))
