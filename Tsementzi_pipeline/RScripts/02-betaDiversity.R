## Script for beta diversity metrics Bray-Curtis and Jaccard for Tsementzi pipeline 

cohort <- "Tsementzi" ## Change to either Angel, Antonio, Gressel, Tsementzi, Walsh
var <- "histology" ## Change to either histology, age, pH, BMI (age, pH, BMI are only presnet for Antonio, Walsh and Tsementzi)

to_keep <- which((colSums(otu_table(eval(parse(text = paste0(cohort, '_', 'Tsementziphyloseq_tree'))))))>1)
otus <- otu_table(eval(parse(text = paste0(cohort, '_', 'Tsementziphyloseq_tree'))))[, to_keep]
FT <- as.data.frame(sample_data(eval(parse(text = paste0(cohort, '_', 'Tsementziphyloseq_tree')))))
FT <- data.frame(FT) %>% 
    mutate(Obesity = if_else(BMI >=30, "Obese", "Non-obese"))
FT$pH <- as.numeric(FT$pH)
new_phylo <- phyloseq(otus, sample_data(FT), 
                      tax_table(eval(parse(text = paste0(cohort, '_', 'Tsementziphyloseq_tree')))), 
                      phy_tree(eval(parse(text = paste0(cohort, '_', 'Tsementziphyloseq_tree')))))
beta_dist <- vegdist(t(otu_table(new_phylo)), index = "jaccard")

ps_ord <- ordinate(new_phylo, method = "NMDS", distance = "jaccard")
png(paste0("~/Desktop/R&R/vaginalMicrobiome/Tsementzi_pipeline/Rplots/", cohort, "_Jaccard_Tsementzipipeline_", var, ".png"))
plot_ordination(new_phylo, ps_ord, type = "samples", color = var, 
                title = paste0(cohort, " using Tsementzi pipeline based on ", var), )
dev.off()

##Univariate 
test <- adonis2(beta_dist ~ eval(parse(text = paste0("sample_data(new_phylo)$", var))), permutations = 100000, na.action = "na.omit")
test

## Multivariate analysis - no interactions 
test_mul <- adonis2(beta_dist ~ sample_data(new_phylo)$histology + sample_data(new_phylo)$age + sample_data(new_phylo)$Obesity + sample_data(new_phylo)$pH,
                    permutations = 100000, na.action = "na.omit")
test_mul
 
## Interactions 
test_mul_int <- adonis2(beta_dist ~ sample_data(new_phylo)$histology + sample_data(new_phylo)$age + 
                          sample_data(new_phylo)$Obesity + sample_data(new_phylo)$pH +
                        sample_data(new_phylo)$age*sample_data(new_phylo)$histology +
                        sample_data(new_phylo)$Obesity*sample_data(new_phylo)$histology +
                        sample_data(new_phylo)$pH*sample_data(new_phylo)$histology,
                    permutations = 100000, na.action = "na.omit")
test_mul_int


