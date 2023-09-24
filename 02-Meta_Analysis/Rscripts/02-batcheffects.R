phylo_use <- eval(parse(text = paste0(cohort, "_SOTAphyloseq_tree_raw")))
phylo_use <- prune_samples(sample_sums(phylo_use)> 0, phylo_use)

## CLR - transformation 
phylo_clr  <- combat_corrected_test %>%
  tax_fix() %>%
  ord_calc(method = "PCA") %>%
  ord_plot(axes = c(1, 2), color = "histology", size = 2.5)

## Combat
library(ConQuR)
library(doParallel)
library(bapred)

data <- data.frame(t(otu_table(phylo_clr)))
FT <- sample_data(phylo_use)
FT_var <- data.frame(FT) %>% mutate(region =
                                      case_when(cohort == "Walsh" ~ "V3-V5", 
                                                cohort == "Antonio" ~ "V3-V5",
                                                cohort == "Gressel" ~ "V4",
                                                cohort == "Chao" ~ "V3-V4",
                                                cohort == "Tsementzi" ~ "V4"))

FT_var$histology <- as.factor(FT_var$histology)
FT_var$region <- as.factor(FT_var$region)

mat_df <- as.matrix(data)
y_num <- as.factor(as.numeric(FT_var$histology))
batch_num <- as.factor(as.numeric(as.factor(FT_var$cohort)))

test <- ba(mat_df, y_num, batch_num, method = "meancenter")

data <- data.frame(t(otu_table(phylo_clr)))
colnames(test$xadj) <- colnames(data)
rownames(test$xadj) <- rownames(data)

combat_corrected_test <- phyloseq(otu_table(test$xadj, taxa_are_rows=FALSE), 
                                  tax_table(phylo_clr), 
                                  sample_data(FT_var))

