cohorts <- c("Antonio", "Angel", "Gressel", "Tsementzi", "Walsh")

i <- 1
for(cohort in cohorts){
  tab <- microbiome::alpha(eval(parse(text = paste0(cohort, '_', 'Tsementziphyloseq_tree'))), index = c("chao1", "evenness_pielou", "diversity_shannon"))
  ps1.meta <- meta(eval(parse(text = paste0(cohort, '_', 'Tsementziphyloseq_tree'))))
  ps1.meta <- subset(ps1.meta, select = c("cohort", "sraID", "histology"))
  ps1.meta$chao1 <- tab$chao1
  ps1.meta$diversity_shannon <- tab$diversity_shannon
  ps1.meta$evenness_pielou <- tab$evenness_pielou
  
  # phylogenetic diversity
  table_otu <- as.data.frame(otu_table(eval(parse(text = paste0(cohort, '_', 'Tsementziphyloseq_tree')))))
  tree_phy <- phy_tree(eval(parse(text = paste0(cohort, '_', 'Tsementziphyloseq_tree'))))
  phy_data <- pd(t(table_otu), tree_phy, include.root=F)
  
  ps1.meta$PD <- phy_data$PD

  assign(paste0(cohorts[i], "_TsementzialphaDiversity"),ps1.meta,.GlobalEnv)
  
  i <- i + 1
}

all_cohorts <- rbind(Antonio_TsementzialphaDiversity,
                     Angel_TsementzialphaDiversity,
                     Gressel_TsementzialphaDiversity,
                     Tsementzi_TsementzialphaDiversity,
                     Walsh_TsementzialphaDiversity)

pathology <- levels(as.factor(all_cohorts$histology))
pathology.pairs <- combn(seq_along(pathology), 2, simplify = FALSE, FUN = function(i)pathology[i])

all_cohorts_long <- gather(all_cohorts, metric, value, chao1:PD, factor_key=TRUE)
all_cohorts_long$log_val <- log(all_cohorts_long$value)

anno_df = compare_means(log_val ~ histology, group.by = c("metric", "cohort"), data = all_cohorts_long, method = "kruskal.test")

bp <- ggplot(all_cohorts_long, aes(x=cohort, y=log_val, fill = histology)) + 
  geom_boxplot(aes(fill=histology)) + 
  scale_fill_manual(values=c("#32CB46", "#FF3333")) +
  facet_wrap(~metric,ncol = 4, scales="free_y") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
bp
