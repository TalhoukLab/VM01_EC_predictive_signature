library(xgboost)
library(caret)
library(MLmetrics)
library(tibble)

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

## Raw 
levels <- c("phylum", "class","order", "family", "genus")
#pipelines in Angel, Antonio_Walsh, Gressel, Tsementzi, InHouse
pipeline <-c("Angel")
#training cohorts in Angel, Antonio, Gressel, Tsementzi, Walsh
training_cohort <- c("Angel")
all_cohorts <- c("Angel", "Antonio", "Gressel", "Tsementzi", "Walsh")
testing_cohorts <- setdiff(all_cohorts,training_cohort)

all.MC = vector(mode="list", length=(length(testing_cohorts)*length(levels)))
all.SC = array(0, dim =  c(length(testing_cohorts), length(levels)))
colnames(all.SC) <- levels
rownames(all.SC) <- testing_cohorts
for(testing_cohort in testing_cohorts){
  for(level in levels){
    if(testing_cohort == paste0(training_cohort, "_test")){
      training_cohort <- paste0(training_cohort, "_train")
    }
    print(paste0("tesing cohort: ", testing_cohort, " level: ", level, " training cohort:", training_cohort))
    train_phylo  <- subset_taxa(eval(parse(text = paste0(training_cohort, "_", pipeline, 'phyloseq_tree_raw'))), c(eval(parse(text=level))!=" " & eval(parse(text=level)) !=""))
    test_phylo  <- subset_taxa(eval(parse(text = paste0(testing_cohort, "_", pipeline, 'phyloseq_tree_raw'))),  c(eval(parse(text=level))!=" " & eval(parse(text=level)) != ""))
    
    #tab_train <- microbiome::alpha(train_phylo, 
    #                                index = c("observed", "chao1", "diversity_gini_simpson", "diversity_shannon"))
    #tab_test <- microbiome::alpha(test_phylo, 
    #                                index = c("observed", "chao1", "diversity_gini_simpson", "diversity_shannon"))
    
    Train_level_phy <- aggregate_taxa(train_phylo, level = level)
    Test_level_phy <- aggregate_taxa(test_phylo, level = level)
    
    otu_train_temp <- as.data.frame(t(as.matrix(Train_level_phy@otu_table@.Data)))
    otu_test_temp <- as.data.frame(t(as.matrix(Test_level_phy@otu_table@.Data)))
    
    if(pipeline == "Angel"){
      otu_test <- fncols_angel(otu_test_temp, setdiff(names(otu_train_temp), names(otu_test_temp)))
      otu_train <- fncols_angel(otu_train_temp, setdiff(names(otu_test_temp), names(otu_train_temp)))
    } else {
      otu_test <- fncols_antonio(otu_test_temp, setdiff(names(otu_train_temp), names(otu_test_temp)))
      otu_train <- fncols_antonio(otu_train_temp, setdiff(names(otu_test_temp), names(otu_train_temp)))
    }
    #otu_train <- merge(otu_train, tab_train, by = "row.names")
    #rownames(otu_train) <- otu_train$Row.names
    #otu_train <- subset(otu_train, select = -c(Row.names))
    #otu_test <- merge(otu_test, tab_test, by = "row.names")
    #rownames(otu_test) <- otu_test$Row.names
    #otu_test <- subset(otu_test, select = -c(Row.names))
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
    all.MC[[paste0(testing_cohort, "_", level)]] <- model 
    otu_test[is.na(otu_test)] <- 0
    y_pred <- predict(model, otu_test, type="prob")[,2]
    y_ref <- Test_level_phy@sam_data$histology
    y_ref_reco <- factor(y_ref, labels=c("B", "C"))
    all.SC[testing_cohort, level] <- PRAUC(y_pred, y_ref_reco)
  }
}

data_wide <- as.data.frame(all.SC)
data_wide$training_cohort <- training_cohort
data_wide$pipeline <- pipeline
data_wide <- tibble::rownames_to_column(data_wide, "testing_cohort")
data_long <- gather(data_wide, level, AUC, order:genus, factor_key=TRUE)
write.csv(data_long, paste0("~/Desktop/train_", training_cohort, "_pipeline_", pipeline, ".csv"), row.names = FALSE)

