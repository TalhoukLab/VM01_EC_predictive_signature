cohorts <- c("Antonio", "Angel", "Gressel", "Tsementzi", "Walsh")
pipelines <- c("Angel", "Antonio", "Tsementzi", "Gressel")

for (pipeline in pipelines){
  for(cohort in cohorts){
    if(pipeline == "Angel"){
      phylo <- eval(parse(text = paste0(cohort, '_', pipeline, 'phyloseq_tree_raw')))
      phylo_log <- microbiome::transform(phylo, 'log10p')
      tab <- microbiome::alpha(phylo_log, index = c("observed", "chao1", "diversity_gini_simpson", "diversity_shannon"))
      ps1.meta <- meta(eval(parse(text = paste0(cohort, '_', pipeline, 'phyloseq_tree'))))
    } else if(pipeline == "Antonio"){
      phylo <- eval(parse(text = paste0(cohort, '_', pipeline, 'phyloseq_tree')))
      tab <- microbiome::alpha(phylo, index = c("observed", "diversity_shannon"))
      ps1.meta <- meta(eval(parse(text = paste0(cohort, '_', pipeline, 'phyloseq_tree'))))
    } else if(pipeline == "Tsementzi"){
      tab <- microbiome::alpha(eval(parse(text = paste0(cohort, '_', pipeline, 'phyloseq_tree'))), index = c("chao1", "evenness_pielou", "diversity_shannon"))
      table_otu <- as.data.frame(otu_table(eval(parse(text = paste0(cohort, '_', 'Tsementziphyloseq_tree')))))
      tree_phy <- phy_tree(eval(parse(text = paste0(cohort, '_', 'Tsementziphyloseq_tree'))))
      phy_data <- pd(t(table_otu), tree_phy, include.root=F)
      
      tab$PD <- phy_data$PD
      ps1.meta <- meta(eval(parse(text = paste0(cohort, '_', pipeline, 'phyloseq_tree'))))
      
    } else {
      #phylo <- eval(parse(text = paste0(cohort, '_', pipeline, 'phyloseq_raw')))
      #phylo_log <- microbiome::transform(phylo, 'log10p')
      tab <- microbiome::alpha(phylo, index = c("chao1", "diversity_shannon", "diversity_fisher"))
      ps1.meta <- meta(eval(parse(text = paste0(cohort, '_', pipeline, 'phyloseq_raw'))))
    }
    ps1.meta <- subset(ps1.meta, select = c("cohort", "sraID", "histology"))
    tab$sraID <- ps1.meta$chao1
    tab$histology <- ps1.meta$histology
    tab$cohort <- ps1.meta$cohort
    assign(paste0(cohort, "_", pipeline, "alphaDiversity"),tab,.GlobalEnv)
  }
  
  all_cohorts <- rbind(eval(parse(text = paste0('Antonio_', pipeline, 'alphaDiversity'))),
                       eval(parse(text = paste0('Angel_', pipeline, 'alphaDiversity'))),
                       eval(parse(text = paste0('Gressel_', pipeline, 'alphaDiversity'))),
                       eval(parse(text = paste0('Tsementzi_', pipeline, 'alphaDiversity'))),
                       eval(parse(text = paste0('Walsh_', pipeline, 'alphaDiversity'))))
  
  pathology <- levels(as.factor(all_cohorts$histology))
  pathology.pairs <- combn(seq_along(pathology), 2, simplify = FALSE, FUN = function(i)pathology[i])
  
  all_cohorts_long <- gather(all_cohorts, metric, value, observed:diversity_shannon, factor_key=TRUE)
  all_cohorts_long$log_val <- log(all_cohorts_long$value)
  
  anno_df = compare_means(log_val ~ histology, group.by = c("metric", "cohort"), data = all_cohorts_long, method = "kruskal.test")
  
  bp <- ggplot(all_cohorts_long, aes(x=cohort, y=log_val, fill = histology)) + 
    geom_boxplot(aes(fill=histology)) + 
    scale_fill_manual(values=c("#32CB46", "#FF3333")) +
    facet_wrap(~metric,ncol = 4, scales="free_y") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  bp
  
}




