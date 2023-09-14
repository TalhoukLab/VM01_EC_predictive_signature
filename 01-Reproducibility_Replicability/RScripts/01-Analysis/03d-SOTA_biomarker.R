## Biomarker validation
library(DESeq2)
library(ALDEx2)
library(ANCOMBC)
library(microViz)
library(tidyverse)
library(VennDiagram)
library(DescTools)

all_cohorts <- c("Antonio", "Chao", "Gressel", "Tsementzi", "Walsh")
all_levels <- c("phylum", "class", "order", "family", "genus")
pipeline <- "SOTA_pipeline"
for(level in all_levels){
  for(cohort in all_cohorts){
    print(paste0(cohort, "_", "level"))
    phylo_use_raw <- eval(parse(text = paste0(cohort, "_SOTAphyloseq_tree_raw")))
    
    phylo_use <- phylo_use_raw %>%
        tax_fix() %>%
        tax_transform("identity", rank = level)
    counts <- data.matrix(as.data.frame(otu_table(phylo_use)))
    counts <- counts + 1
    if(cohort == "Chao"){
      print("Chao_DESeq2")
      sampleData <- data.frame(sample_data(phylo_use))
      sampleData$histology <- as.factor(sampleData$histology)
      sampleData_scaled <- sampleData %>% mutate_if(is.numeric, scale)
      dds <- DESeqDataSetFromMatrix(countData = counts,
                                    colData = sampleData_scaled,
                                    design= ~  age + histology)
    } else if(cohort == "Gressel"){
      sampleData <- data.frame(sample_data(phylo_use))
      sampleData$histology <- as.factor(sampleData$histology)
      
      print("Gressel_DESeq2")
      dds <- DESeqDataSetFromMatrix(countData = counts,
                                    colData = data.frame(sample_data(phylo_use)),
                                    design= ~ histology)
    } else if(cohort == "Tsementzi"){
      sampleData <- data.frame(sample_data(phylo_use))
      sampleData$histology <- as.factor(sampleData$histology)
      
      sampleData_scaled <- sampleData %>% mutate_if(is.numeric, scale)
      sampleData_scaled <- na.omit(sampleData_scaled)
      sampleData <- na.omit(sampleData)
      phylo_scaled <- phyloseq(otu_table(counts, taxa_are_rows=T), sample_data(sampleData_scaled), 
                               tax_table(phylo_use))
      
      print("Tsementzi_DESeq2")
      dds <- DESeqDataSetFromMatrix(countData = data.matrix(as.data.frame(otu_table(phylo_scaled))),
                                    colData =as.data.frame(sample_data(phylo_scaled)),
                                    design= ~  BMI + age + ethnicityRecode + pHRecoded + histology)
      
    } else if(cohort == "Antonio"){
      sampleData <- data.frame(sample_data(phylo_use))
      sampleData$histology <- as.factor(sampleData$histology)
      
      sampleData$Stage[is.na(sampleData$Stage)] <- "Not applicable"
      sampleData$menopausal.status <- factor(sampleData$menopausal.status)
      sampleData$pHRecoded <- factor(sampleData$pHRecoded)
      sampleData_scaled <- sampleData %>% mutate_if(is.numeric, scale)
      sampleData_scaled <- na.omit(sampleData_scaled)
      
      
      phylo_scaled <- phyloseq(otu_table(counts, taxa_are_rows=T), sample_data(sampleData_scaled), 
                               tax_table(phylo_use))
      
      #sampleData <- data.frame(sample_data(phylo_use))
      #sampleData[is.na(sampleData)] <- "not collected"
      print("Antonio_DESeq2")
      dds <- DESeqDataSetFromMatrix(countData = data.matrix(as.data.frame(otu_table(phylo_scaled))),
                                    colData =as.data.frame(sample_data(phylo_scaled)),
                                    design= ~  menopausal.status + age + BMI + pHRecoded + histology)
      
      
      } else {
        sampleData <- data.frame(sample_data(phylo_use))
        sampleData$histology <- as.factor(sampleData$histology)
        
        sampleData$BMI <- as.numeric(sampleData$BMI)
        sampleData_scaled <- sampleData %>% mutate_if(is.numeric, scale)
        sampleData_scaled <- na.omit(sampleData_scaled)
        sampleData <- na.omit(sampleData)
        phylo_scaled <- phyloseq(otu_table(counts, taxa_are_rows=T), sample_data(sampleData_scaled), 
                                 tax_table(phylo_use))
        
        
        print("Walsh_DESeq2")
        dds <- DESeqDataSetFromMatrix(countData = data.matrix(as.data.frame(otu_table(phylo_scaled))),
                                      colData =as.data.frame(sample_data(phylo_scaled)),
                                      design= ~  menopausal.status + BMI + age + ethnicityRecoded + pHRecoded + histology)
        
    }
    dds <- DESeq(dds)
    assign(paste0(cohort, "_", level, "_deseq2"), dds,.GlobalEnv)
    
    res_deseq2 <- as.data.frame(results(dds, contrast=c("histology","EC","Benign")))
    res_deseq2_fil <- subset(res_deseq2, padj <= 0.05)
    res_deseq2_fil <- subset(res_deseq2_fil, abs(log2FoldChange) >= 3)
    res_deseq2_fil$taxa <- rownames(res_deseq2_fil)
    
    if(nrow(res_deseq2_fil)>0){
      res_deseq2_fil$cohort <- cohort
      res_deseq2_fil$pipeline <- pipeline
      res_deseq2_fil$direction <- with(res_deseq2_fil, ifelse(log2FoldChange>0, 'Positive', 'Negative'))
      res_deseq2_fil <- res_deseq2_fil %>% 
        dplyr::select(c(taxa, cohort, pipeline, direction))
    }
    
    assign(paste0(cohort, "_", level, "_deseq2_df"), res_deseq2_fil,.GlobalEnv)
    
  }
}

 
 