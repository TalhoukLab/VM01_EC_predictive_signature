## This script is for Angel beta diversity. They only use unweighted uniFrac distance with an ordination plot. 

cohorts <- c("Angel", "Antonio", "Gressel", " Tsementzi", "Walsh")

for(cohort in cohorts){
  dist <- phyloseq::distance(phylo_use, method = "unifrac")
  ps_ord <- ordinate(phylo_use, distance = dist)
  png(paste0("~/Desktop/R&R/vaginalMicrobiome/Angel_pipeline/Rplots/", cohort, "_uunifrac_Angelpipeline_histolog.png"))
  plot_ordination(eval(parse(text = paste0(cohort, '_', 'Angelphyloseq_tree'))), 
                  ps_ord, type = "samples", color = "histology", 
                  title = paste0(cohort, " using Angel pipeline based on ", "histology"))
  
  dev.off()
  
}
  