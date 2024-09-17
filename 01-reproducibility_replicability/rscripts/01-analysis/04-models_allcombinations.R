library(xgboost)
library(caret)
library(MLmetrics)
library(tibble)
library(stringr)

## Source all data prep files 
source(file = "../vaginalMicrobiome/01-Reproducibility_Replicability/RScripts/00-DataPrep/Antonio_dataPrep.R")
source(file = "../vaginalMicrobiome/01-Reproducibility_Replicability/RScripts/00-DataPrep/Chao_dataPrep.R")
source(file = "../vaginalMicrobiome/01-Reproducibility_Replicability/RScripts/00-DataPrep/Gressel_dataPrep.R")
source(file = "../vaginalMicrobiome/01-Reproducibility_Replicability/RScripts/00-DataPrep/Tsementzi_dataPrep.R")
source(file = "../vaginalMicrobiome/01-Reproducibility_Replicability/RScripts/00-DataPrep/SOTA_dataPrep.R")

## If the train or test set does not have a particular taxa present, 
## we manually add it in with an abundance of  0 (depending on the preprocessing step of the pipeline)
fncols_angel <- function(data, cname) {
  add <-cname[!cname%in%names(data)]
  if(length(add)!=0) data[add] <- log10(0.001)
  data
}

fncols_antonio <- function(data, cname) {
  add <-cname[!cname%in%names(data)]
  if(length(add)!=0) data[add] <- 0
  data
}

all_cohorts <- c("Antonio", "Chao", "Gressel", "Tsementzi", "Walsh")


## ML replicability was tested using the counts as returned by the pipeline (usually saved as *_tree)

mlRep <- function(training_cohort, pipeline, level, type){
  train_phylo  <- subset_taxa(eval(parse(text = paste0(training_cohort, "_", pipeline, 'phyloseq_tree'))), 
                              c(eval(parse(text=level))!=" " & eval(parse(text=level)) !=""))
  tab_train <- microbiome::alpha(train_phylo, 
                                 index = c("observed", "chao1", "diversity_shannon"))
  Train_level_phy <- aggregate_taxa(train_phylo, level = level)
  otu_train <- as.data.frame(t(as.matrix(Train_level_phy@otu_table@.Data)))
  otu_train <- merge(otu_train, tab_train, by = "row.names")
  rownames(otu_train) <- otu_train$Row.names
  otu_train <- subset(otu_train, select = -c(Row.names))

  if(type == "cross"){
    testing_cohorts <- setdiff(all_cohorts, training_cohort)
  } else {
    testing_cohorts <- c(training_cohort)
    testing_cohorts <- testing_cohorts %>% str_replace("_.*", "_test")
    
  }
  for(testing_cohort in testing_cohorts){
    test_phylo  <- subset_taxa(eval(parse(text = paste0(testing_cohort, "_", pipeline, 'phyloseq_tree'))),  
                               c(eval(parse(text=level))!=" " & eval(parse(text=level)) != ""))
    tab_test <- microbiome::alpha(test_phylo, 
                                  index = c("observed", "chao1", "diversity_shannon"))
    Test_level_phy <- aggregate_taxa(test_phylo, level = level)
    otu_test <- as.data.frame(t(as.matrix(Test_level_phy@otu_table@.Data)))
    otu_test <- merge(otu_test, tab_test, by = "row.names")
    rownames(otu_test) <- otu_test$Row.names
    otu_test <- subset(otu_test, select = -c(Row.names))
    if(pipeline == "Chao"){
      otu_train <- fncols_angel(otu_train, setdiff(names(otu_test), names(otu_train)))
    } else {
      otu_train <- fncols_antonio(otu_train, setdiff(names(otu_test), names(otu_train)))
    }
    y_ref <- Test_level_phy@sam_data$histology
    y_ref_reco <- factor(y_ref, labels=c("B", "C"))
    assign(paste0(testing_cohort, "_", pipeline, "_testing_data"),otu_test,.GlobalEnv)
    assign(paste0(testing_cohort, "_", pipeline, "_testing_data_hist"),y_ref_reco,.GlobalEnv)
  }
  
  labels <- Train_level_phy@sam_data$histology
  y <- factor(labels, labels=c("B", "C"))
  set.seed(123)
  otu_train <- cbind(otu_train, y)
  
  train_control = trainControl(method = "cv", number = 5, search = "grid",     
                               summaryFunction = prSummary, classProbs = TRUE)
  
  gbmGrid <-  expand.grid(max_depth = c(3, 5, 7, 9, 11), 
                          nrounds = (1:10)*5,    # number of trees
                          # default values below
                          eta = c(0.1, 0.2, 0.3),
                          gamma = 0,
                          subsample = c(0.5, 1),
                          min_child_weight = 1,
                          colsample_bytree = 0.6)
  
  model = train(y~., data = otu_train, method = "xgbTree", na.action = na.pass, 
                trControl = train_control, tuneGrid = gbmGrid,
                metric="AUC",
                objective="binary:logistic",
                nthread = 1, verbosity = 0)
  
  all.SC = array(0, dim =  c(length(testing_cohorts), 2))
  rownames(all.SC) <- testing_cohorts
  for(testing_cohort in testing_cohorts){
    otu_test_use <- eval(parse(text = paste0(testing_cohort, "_", pipeline, '_testing_data')))
    y_ref_recoded <- eval(parse(text = paste0(testing_cohort, "_", pipeline, "_testing_data_hist")))
    if(pipeline == "Chao"){
      otu_test_use <- fncols_angel(otu_test_use, setdiff(names(otu_train), names(otu_test_use)))
    } else {
      otu_test_use <- fncols_antonio(otu_test_use, setdiff(names(otu_train), names(otu_test_use)))
    }
    
    otu_test_use[is.na(otu_test_use)] <- 0
    y_pred <- predict(model, otu_test_use, type="prob")[,2]
    y_pred_MC <- predict(model, otu_test_use, type="raw")
    all.SC[testing_cohort, 1] <- PRAUC(y_pred, y_ref_recoded)
    all.SC[testing_cohort, 2] <- mean(y_pred_MC != y_ref_recoded)
  }
  return(all.SC)
}


all_pipelines <- c("Antonio", "Chao", "Gressel", "Tsementzi", "SOTA")
all_levels <- c("phylum", "order", "class", "family", "genus")
  
columns = c("AUC","MC", "testing_cohort","training_cohort", "pipeline", "level", "type") 
df = data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(df) = columns

for(training_cohort in all_cohorts){
  for(pipeline in all_pipelines){
    for(level in all_levels){
      scores <- mlRep(training_cohort = training_cohort, pipeline = pipeline, level = level, type = "cross")
      test <- as.data.frame(scores)
      test$testing_cohort <- rownames(test)
      rownames(test) <- NULL
      test$training_cohort <- training_cohort
      test$pipeline <- pipeline
      test$level <- level
      df <- rbind(df, test)
    }
  }
}

write.csv(df, "../vaginalMicrobiome/01-Reproducibility_Replicability/Results/04-MLRep/crossLab.csv", row.names = FALSE)

columns = c("AUC","MC", "testing_cohort","training_cohort", "pipeline", "level", "type") 
df = data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(df) = columns

within_training_cohorts <- c("Antonio_train", "Chao_train", "Gressel_train", "Tsementzi_train", "Walsh_train")

for(training_cohort in within_training_cohorts){
    for(level in all_levels){
      pipelines <- c(training_cohort %>% str_replace("_.*", ""), "SOTA")
      if(training_cohort == "Walsh_train"){
        pipelines <- c("Antonio", "SOTA")
      }
      for(pipeline in pipelines){
        scores <- mlRep(training_cohort = training_cohort, pipeline = pipeline, level = level, type = "within")
        test <- as.data.frame(scores)
        test$testing_cohort <- rownames(test)
        rownames(test) <- NULL
        test$training_cohort <- training_cohort
        test$pipeline <- pipeline
        test$level <- level
        df <- rbind(df, test)
      }
    }
}
write.csv(df, "../vaginalMicrobiome/01-Reproducibility_Replicability/Results/04-MLRep/withinLab.csv", row.names = FALSE)

