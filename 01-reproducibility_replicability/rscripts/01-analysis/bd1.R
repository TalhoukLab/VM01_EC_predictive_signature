library(phyloseq)
library(ggplot2)
library(microViz)

mul_var_analysis <- function(dist, cohort, raw, normal){
  if(dist == "unifrac"){
    phylo_to_use <-  raw
  } else {
    phylo_to_use <- normal
    phylo_to_use <- microbiome::transform(phylo_to_use, transform = "log")
  }
  dist_mat = phyloseq::distance(phylo_to_use, method = dist)
 # dist_mat <- as.matrix(dist_mat)
  
  if(cohort == "Chao"){
    test_mul <- adonis2(dist_mat ~ sample_data(phylo_to_use)$histology + 
                          sample_data(phylo_to_use)$age,
                        permutations = 100, na.action = "na.omit", by = "margin")
  } 
  if(cohort == "Antonio"){
    test_mul <- adonis2(dist_mat ~ sample_data(phylo_to_use)$histology + 
                          sample_data(phylo_to_use)$age + 
                          sample_data(phylo_to_use)$BMI ,
                        permutations = 100, na.action = "na.omit", by = "margin")
  } 
  if(cohort == "Tsementzi"){
    test_mul <- adonis2(dist_mat ~ sample_data(phylo_to_use)$histology + 
                          sample_data(phylo_to_use)$age + 
                          sample_data(phylo_to_use)$bmi + 
                          as.factor(sample_data(phylo_to_use)$ethnicityRecode),
                        permutations = 100, na.action = "na.omit", by = "margin")
  } 
  if(cohort == 'Walsh'){
    test_mul <- adonis2(dist_mat ~ sample_data(phylo_to_use)$histology + 
                          sample_data(phylo_to_use)$age + 
                          as.numeric(sample_data(phylo_to_use)$bmi) + 
                          as.factor(sample_data(phylo_to_use)$ethnicityRecode), 
                        permutations = 100, na.action = "na.omit", by = "margin")
  }
  if(cohort == "Gressel"){
    test_mul <- adonis2(dist_mat ~ sample_data(phylo_to_use)$histology,
                        permutations = 100, na.action = "na.omit", by = "margin")
  }
  return(test_mul)
}

cohorts <- c("Antonio", "Chao", "Gressel", "Tsementzi", "Walsh")
dists <- c("bray", "jaccard", "unifrac", "wunifrac", "jsd")
colnames_list <- c("R2_bray", "p-value_bray", "R2_jaccard", "p-value_jaccard", "R2_unifrac", "p-value_unifrac", "R2_wunifrac", "p-value_wunifrac", "R2_jsd", "p-value_jsd")


Walsh_mulvar <- data.frame(base::sapply(base::sapply(dists, mul_var_analysis, cohort = "Walsh", raw = Walsh_dada2phyloseq_tree_raw, normal = Walsh_dada2phyloseq_tree)[c(3, 5), ], cbind))
colnames(Walsh_mulvar) <- colnames_list
Walsh_mulvar$covariate <- c("Health condition", "Age", "BMI", "Ethnicity", "Residual", "total")
Walsh_mulvar$cohort <- "Walsh"

Tsementzi_mulvar <- data.frame(base::sapply(base::sapply(dists, mul_var_analysis, cohort = "Tsementzi", raw = Tsementzi_dada2phyloseq_tree_raw, normal = Tsementzi_dada2phyloseq_tree)[c(3, 5), ], cbind))
colnames(Tsementzi_mulvar) <- colnames_list
Tsementzi_mulvar$covariate <- c("Health condition", "Age", "BMI", "Ethnicity", "Residual", "total")
Tsementzi_mulvar$cohort  <- "Tsementzi"

Gressel_mulvar <- data.frame(base::sapply(base::sapply(dists, mul_var_analysis, cohort = "Gressel", raw = Gressel_dada2phyloseq_tree_raw, normal = Gressel_dada2phyloseq_tree)[c(3, 5), ], cbind))
colnames(Gressel_mulvar) <- colnames_list
Gressel_mulvar$covariate <- c("Health condition", "Residual", "total")
Gressel_mulvar$cohort <- "Gressel"

Chao_mulvar <- data.frame(base::sapply(base::sapply(dists, mul_var_analysis, cohort = "Chao", raw = Chao_dada2phyloseq_tree_raw, normal = Chao_dada2phyloseq_tree)[c(3, 5), ], cbind))
colnames(Chao_mulvar) <- colnames_list
Chao_mulvar$covariate <- c("Health condition", "Age", "Residual", "total")
Chao_mulvar$cohort <- "Chao"

Antonio_mulvar <- data.frame(base::sapply(base::sapply(dists, mul_var_analysis, cohort = "Antonio", raw = Antonio_dada2phyloseq_tree_raw, normal = Antonio_dada2phyloseq_tree)[c(3, 5), ], cbind))
colnames(Antonio_mulvar) <-  colnames_list
Antonio_mulvar$covariate <- c("Health condition", "Age", "BMI", "Residual", "total")
Antonio_mulvar$cohort <- "Antonio"

all <- rbind(Walsh_mulvar, Tsementzi_mulvar, Gressel_mulvar, Chao_mulvar, Antonio_mulvar)
all_long <- pivot_longer(all, "R2_bray":"p-value_jsd")
all_long <- separate(all_long, col = name, into = c("type", "distance"), sep = "_")
all_long_wide <-pivot_wider(all_long, id_cols = c("covariate", "cohort", "distance"), names_from = "type", values_from = "value")
all_long_wide <- unique(all_long_wide)

all_long_wide <- all_long_wide %>% dplyr::filter(covariate %in% c("Health condition", "BMI", "Age", "Ethnicity"))
all_long_wide$R2 <- round(all_long_wide$R2, 2)
all_long_wide$sig <- ifelse(all_long_wide$`p-value` <=0.05, "Sig", "NS")
all_long_wide$covariate <- factor(all_long_wide$covariate, levels = c("Health condition", "BMI", "Age", "Ethnicity"))
all_long_wide$cohort <- factor(all_long_wide$cohort, levels = c("Antonio", "Walsh", "Tsementzi", "Gressel", "Chao"))
all_long_wide$distance <- factor(all_long_wide$distance, levels = c("bray", "jaccard", "unifrac", "wunifrac", "jsd"),
                      labels = c("Bray", "Jaccard", "UniFrac", "wUniFrac", "JSD"))
library(ggplot2)
library(latex2exp)
pdf("~/Desktop/dada2_all.pdf",  width=25, height=8)
beta <- ggplot(all_long_wide, aes(cohort, covariate, fill= R2)) +  geom_tile(aes(fill = R2)) + facet_grid(. ~ distance) + 
  geom_tile(data = all_long_wide, fill="transparent", color = ifelse(all_long_wide$sig=="Sig", 'black', NA), size = 0.5) +
  geom_text(aes(label = R2), color = "black", size = 6.0) +
  scale_fill_gradient(low = "white", high = "red", name = unname(TeX(c("$R^2$"))))  + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title=element_text(size=21),
        strip.text.x = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.position = "bottom",
        legend.title = element_text(size=21),
        legend.text = element_text(size=10), plot.title = element_text(size = 23, hjust = 0.5, face = "bold")) + ggtitle("DADA2 - Beta diversity")
print(beta)
dev.off()


boxplot_bd <- function(dist, cohort, raw, normal){
  if(dist == "unifrac"){
    phylo_to_use <- raw
  } else {
    phylo_to_use <- normal
  }
  dist_mat = phyloseq::distance(phylo_to_use, method = dist)
  dist_mat <- as.matrix(dist_mat)
  
  sub_dist <- list()
  groups_all <- factor(sample_data(phylo_to_use)$histology)
  for(group in levels(groups_all)) { 
    row_group <- which(groups_all == group)
    sample_group <- sample_names(phylo_to_use)[row_group]
    sub_dist[[group]] <- dist_mat[ sample_group, sample_group]
    sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
  }
  
  groups <- melt(sub_dist)
  df.dist <- groups[complete.cases(groups), ]
  #df.dist$L1 <- factor(df.dist$L1, levels=names(sub_dist))
  df.dist$cohort <- cohort
  df.dist$dist <- dist
  return(df.dist)
}

colnames_list <- c("Var 1", "Var 2", "Bray", "Jaccard", "UniFrac", "wUniFrac", "JSD", "histology", "cohort")
Walsh_bp_bd <- data.frame(base::sapply(t(base::sapply(dists, boxplot_bd, cohort = "Walsh", raw = Walsh_dada2phyloseq_tree_raw, normal = Walsh_dada2phyloseq_tree))[1:5, ], rbind)) %>%
                    dplyr::select(c(X1, X6, X11:X16, X21))
colnames(Walsh_bp_bd) <- colnames_list

Tsementzi_bp_bd <- data.frame(base::sapply(t(base::sapply(dists, boxplot_bd, cohort = "Tsementzi", raw = Tsementzi_dada2phyloseq_tree_raw, normal =Tsementzi_dada2phyloseq_tree ))[1:5, ], rbind)) %>%
                        dplyr::select(c(X1, X6, X11:X16, X21))
colnames(Tsementzi_bp_bd) <- colnames_list

Gressel_bp_bd <- data.frame(base::sapply(t(base::sapply(dists, boxplot_bd, cohort = "Gressel", raw = Gressel_dada2phyloseq_tree_raw, normal = Gressel_dada2phyloseq_tree))[1:5, ], rbind)) %>%
  dplyr::select(c(X1, X6, X11:X16, X21))
colnames(Gressel_bp_bd) <- colnames_list

Chao_bp_bd <- data.frame(base::sapply(t(base::sapply(dists, boxplot_bd, cohort = "Chao", raw = Chao_dada2phyloseq_tree_raw, normal = Chao_dada2phyloseq_tree))[1:5, ], rbind)) %>%
  dplyr::select(c(X1, X6, X11:X16, X21))
colnames(Chao_bp_bd) <- colnames_list

Antonio_bp_bd <- data.frame(base::sapply(t(base::sapply(dists, boxplot_bd, cohort = "Antonio", raw = Antonio_dada2phyloseq_tree_raw, normal = Antonio_dada2phyloseq_tree))[1:5, ], rbind)) %>%
  dplyr::select(c(X1, X6, X11:X16, X21))
colnames(Antonio_bp_bd) <- colnames_list

all_bp_bd <- rbind(Walsh_bp_bd, Tsementzi_bp_bd, Gressel_bp_bd, Chao_bp_bd, Antonio_bp_bd)
all_bp_bd_long <- pivot_longer(all_bp_bd, cols = "Bray":"JSD", names_to = "distance")
all_bp_bd_long$value <- as.numeric(all_bp_bd_long$value)
all_bp_bd_long$distance <- factor(all_bp_bd_long$distance, levels = c("Bray", "Jaccard", "UniFrac", "wUniFrac", "JSD"),
                                 labels = c("Bray", "Jaccard", "UniFrac", "wUniFrac", "JSD"))

pdf("~/Desktop/dada2_all_bp.pdf",  width=15, height=5)
ggplot(all_bp_bd_long, aes(x=cohort, y=value, fill=histology)) +
  geom_boxplot(aes(fill=histology), outlier.shape = NA) + 
  scale_fill_manual(values=c("yellowgreen", "tomato3")) + 
 facet_grid(. ~ distance) +
  ylab("Distance") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
dev.off()





