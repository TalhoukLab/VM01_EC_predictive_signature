## This script will format the input for Lefse

agg_taxa <- aggregate_taxa(Tsementzi_Tsementziphyloseq_tree_raw, level = "family")

data <- data.frame(otu_table(Tsementzi_Tsementziphyloseq_tree_raw))
taxa_names <- data.frame(tax_table(Tsementzi_Tsementziphyloseq_tree_raw))
taxa_names$merged <- paste(rownames(taxa_names), taxa_names$kingdom, taxa_names$phylum, taxa_names$class,
                           taxa_names$order, taxa_names$family, taxa_names$genus, sep = "|")
rownames(data) <- taxa_names$merged
flip_data <- data.frame(t(data), check.names = F)


FT <- as.data.frame(sample_data(Tsementzi_Tsementziphyloseq_tree_raw))
FT <- data.frame(FT) %>% 
  mutate(Obesity = if_else(BMI >=30, "Obese", "Non-obese"))
FT$pH <- as.numeric(FT$pH)


covariates = FT[, c('age', 'histology', 'pH', 'sraID')]
data_full <- cbind(covariates, flip_data)
data_full <- na.omit(data_full)
data_full1 <- data.frame(t(data_full), check.names = F)

write.table(data_full1, "~/Desktop/lefse.txt", quote=FALSE, sep = "\t", col.names = FALSE)
