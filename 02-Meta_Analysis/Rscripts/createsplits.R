library(caret)

createPartitions <- function(cohort){
  FT <- read.delim(file.path(paste0('../vaginalMicrobiome/01-reproducibility_replicability/00-helperfiles/', cohort, 'FT.csv')), header = TRUE, sep = ",")
  FT$histology[FT$histology == "ACH"] <- "EC"
  rownames(FT) <- FT$sraID
  train_idx <- createDataPartition(FT$histology, p = 0.7, times = 1, list = FALSE)
  train <- FT[train_idx,]
  rest <- FT[-train_idx,]
  test_idx <- createDataPartition(rest$histology, p = 0.5, times = 1, list = FALSE)
  test <- rest[test_idx,]
  validation <-  rest[-test_idx,]
  splits <- list(train=train, test=test, validation=validation)  
  write.table(train$sraID, file.path(paste0('../vaginalMicrobiome/02-meta_analysis/00-helperfiles/', cohort, '_train_sraID.txt')), quote = F, sep = ",", row.names = F, col.names = F)
  write.csv(train, file.path(paste0('../vaginalMicrobiome/02-meta_analysis/00-helperfiles/', cohort, '_train.txt')), quote = FALSE)
  write.csv(test, file.path(paste0('../vaginalMicrobiome/02-meta_analysis/00-helperfiles/', cohort, '_test.txt')), quote = FALSE)
  write.table(test$sraID, file.path(paste0('../vaginalMicrobiome/02-meta_analysis/00-helperfiles/', cohort, '_test_sraID.txt')), quote = F, sep = ",", row.names = F, col.names = F)
  write.csv(validation, file.path(paste0('../vaginalMicrobiome/02-meta_analysis/00-helperfiles/', cohort, '_validation.txt')), quote  = FALSE)
  write.table(validation$sraID, file.path(paste0('../vaginalMicrobiome/02-meta_analysis/00-helperfiles/', cohort, '_validation_sraID.txt')), quote = F, sep = ",", row.names = F, col.names = F)
  return(splits)
}

antonio_splits <- createPartitions("Antonio")
walsh_splits <- createPartitions("Walsh")
tsementzi_splits <- createPartitions("Tsementzi")
gressel_splits <- createPartitions("Gressel")
chao_splits <- createPartitions("Chao")

