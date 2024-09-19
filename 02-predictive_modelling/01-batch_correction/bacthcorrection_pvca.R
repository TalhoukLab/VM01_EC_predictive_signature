
methods <- c("combat", "combatseq", "limma", "plsda", "splsda", "conqur", "bmc", "pn")
levels <- c("phylum", "class", "order", "family", "genus", "species")
testing <- c("Antonio", "Chao", "Gressel", "Tsementzi", "Walsh")

pdf(paste0("~/Desktop/no_primers_new.pdf"))

for(tester in testing){
  for(level in levels){
    for(method in methods){
      corrected <- readRDS(paste0("~/Desktop/thesis/VM02_meta_analysis/src/results/", 
                                  tester, "_", level, "corrected.rds"))
      corrected_counts <- corrected[[method]]
      if(ncol(corrected_counts)<50){
        next
      }
      df <- data.frame(cbind("histology"=corrected[['labels']],
                       "batches"=corrected[['batches']]))
      df <- df[match(rownames(corrected_counts), rownames(df)),]
      adf <- Biobase::AnnotatedDataFrame(df)
      eds <- Biobase::ExpressionSet(t(corrected_counts), phenoData = adf)
      batch.factors <- c("histology", "batches")
      pvcaObj <- pvcaBatchAssess(eds, batch.factors, 0.1) 
      pvca_df <- data.frame(cbind("labels"=pvcaObj$label, "variance"=pvcaObj$dat[1,]))
      pvca_df$variance <- as.numeric(pvca_df$variance)
      print(ggplot(data=pvca_df, aes(x=labels, y=variance)) +
        geom_bar(stat="identity") + coord_flip() + ggtitle(paste0(tester, "_", level, "_", method)))
    }
  }
}
dev.off()

