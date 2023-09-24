library(PLSDAbatch)
library(mixOmics)

cohort <- AWT_SOTAphyloseq_tree_raw
vm.count <- t(data.frame(otu_table(cohort)))
dim(vm.count)
vm.metadata <- data.frame(sample_data(cohort))
vm.batch <- as.factor(vm.metadata$cohort)
vm.trt <- as.factor(vm.metadata$histology)
names(vm.batch) <- names(vm.trt) <- rownames(vm.metadata)

vm.filter.res <- PreFL(data = vm.count)
vm.filter <- vm.filter.res$data.filter
dim(vm.filter)

vm.filter.res$zero.prob
sum(vm.filter==0)/(nrow(vm.filter) * ncol(vm.filter))

vm.clr <- logratio.transfo(X = vm.filter, logratio = "CLR", offset = 1)
class(vm.clr) = "matrix"

vm.pca.before <- pca(vm.clr, ncomp = 3, scale = TRUE)
Scatter_Density(object = vm.pca.before, batch = vm.batch,
                trt = vm.trt, title = "raw", trt.legend.title = "histology")


vm.clr.s <- scale(vm.clr, center = TRUE, scale = TRUE)
vm.clr.ss <- scale(t(vm.clr.s), center = TRUE, scale = TRUE)

vm.anno_col <- data.frame(Batch = vm.batch, Treatment = vm.trt)
vm.anno_colors <- list(Batch = color.mixo(seq_len(length(unique(vm.batch)))),
                       Treatment = pb_color(seq_len(2)))
names(vm.anno_colors$Batch) = levels(vm.batch)
names(vm.anno_colors$Treatment) = levels(vm.trt)

pheatmap::pheatmap((vm.clr.ss), 
                   cluster_rows = FALSE, 
                   fontsize_row = 4, 
                   fontsize_col = 6,
                   fontsize = 8,
                   clustering_distance_rows = 'euclidean',
                   clustering_method = 'ward.D',
                   treeheight_row = 30,
                   annotation_col = vm.anno_col,
                   annotation_colors = vm.anno_colors,
                   border_color = 'NA',
                   main = 'VM data - Scaled')


vm.trt.tune <- plsda(X = vm.clr, Y = vm.trt, ncomp = 5)
vm.trt.tune$prop_expl_var

vm.batch.tune <- PLSDA_batch(X = vm.clr, 
                             Y.trt = vm.trt, 
                             Y.bat = vm.batch,
                             ncomp.trt = 1, ncomp.bat = 10, balance = FALSE)
vm.batch.tune$explained_variance.bat #5
sum(vm.batch.tune$explained_variance.bat$Y[seq_len(5)])

vm.PLSDA_batch.res <- PLSDA_batch(X = vm.clr, 
                                  Y.trt = vm.trt, 
                                  Y.bat = vm.batch,
                                  ncomp.trt = 1, ncomp.bat = 5, balance = FALSE)
vm.PLSDA_batch <- vm.PLSDA_batch.res$X.nobatch

set.seed(109)
vm.test.keepX = c(seq(1, 10, 1), seq(20, 300, 10), 
                  seq(300, 938, 50), 938)
vm.trt.tune.v <- tune.splsda(X = vm.clr, Y = vm.trt, 
                             ncomp = 1,
                             test.keepX = vm.test.keepX, 
                             validation = 'Mfold', folds = 4, 
                             nrepeat = 50, cpus = 30)
vm.trt.tune.v$choice.keepX #140


vm.batch.tune <- PLSDA_batch(X = vm.clr, 
                             Y.trt = vm.trt, Y.bat = vm.batch,
                             ncomp.trt = 1, keepX.trt = 150,
                             ncomp.bat = 5, balance = FALSE)
vm.batch.tune$explained_variance.bat 


vm.sPLSDA_batch.res <- PLSDA_batch(X = vm.clr, 
                                   Y.trt = vm.trt, 
                                   Y.bat = vm.batch,
                                   ncomp.trt = 1, 
                                   keepX.trt = 150,
                                   ncomp.bat = 5, balance = FALSE)
vm.sPLSDA_batch <- vm.sPLSDA_batch.res$X.nobatch


vm.pca.before <- pca(vm.clr, ncomp = 3, scale = TRUE)
vm.pca.PLSDA_batch <- pca(vm.PLSDA_batch, ncomp = 3, scale = TRUE)
vm.pca.sPLSDA_batch <- pca(vm.sPLSDA_batch, ncomp = 3, scale = TRUE)


vm.pca.before.plot <- Scatter_Density(object = vm.pca.before, 
                                      batch = vm.batch, legend.cex = 1,
                                      trt = vm.trt, title.cex = 3,density.lwd = 0.5,
                                      title = 'Before correction')
vm.pca.PLSDA_batch.plot <- Scatter_Density(object = vm.pca.PLSDA_batch, 
                                           batch = vm.batch, 
                                           trt = vm.trt, 
                                           title = 'PLSDA-batch')
vm.pca.sPLSDA_batch.plot <- Scatter_Density(object = vm.pca.sPLSDA_batch, 
                                            batch = vm.batch, legend.cex = 1,
                                            trt = vm.trt, title.cex = 3,density.lwd = 0.5,
                                            title = 'sPLSDA-batch')
#library(vegan3d)
#ordirgl(vm.pca.sPLSDA_batch, display = "sites", 
#        col= rainbow(5)[as.numeric(vm.batch)])

library(pvca)
pct_threshold <- 0.6
df <- vm.metadata
df$pHRecoded <- as.factor(df$pHRecoded)
df$menopausal.status <- as.factor(df$menopausal.status)
df$cohort <- as.factor(df$cohort)
df$histology <- as.factor(df$histology)
df <- df %>% dplyr::select(c("cohort", "histology", "pHRecoded",
                             "menopausal.status", "BMI"))
batch.factors <- c("cohort", "histology", "pHRecoded",
                   "menopausal.status", "BMI")

adf <- Biobase::AnnotatedDataFrame(df)
eds <- Biobase::ExpressionSet(t(vm.clr), phenoData = adf)

pvcaObj <- pvcaBatchAssess (eds, batch.factors, pct_threshold) 

bp <- barplot(pvcaObj$dat, horiz = T, 
              xlab = "Weighted average proportion variance", xlim= c(0,1.1),
              col = c("blue"), las=2, main="PVCA estimation bar chart - CLR")
axis(2, at = bp, labels = pvcaObj$label, cex.axis = 1.2, las=2)
values = pvcaObj$dat
new_values = round(values , 3)
text(pvcaObj$dat,bp, labels = new_values, pos=4, cex = 1) 

