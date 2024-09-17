library(dada2)
library(DECIPHER)
library(dplyr)
library(tibble)
library(phyloseq)
library(speedyseq)
library(msa)
library(phangorn)
library(Biostrings)
require(ShortRead)

args = commandArgs(trailingOnly=TRUE)
cohort <- args[1]
print(cohort)
workspace <- args[2]
path <- paste0(workspace, "/01b-filtering")
resPath <- paste0(workspace, "/denoised")

fnFs <- sort(list.files(path, pattern="_1_paired.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2_paired.fastq", full.names = TRUE))

if(cohort=="mock"){
  sample.names <- sapply(strsplit(basename(fnFs), '^[^_]*?[_][^_]*?(*SKIP)(*F)|_', perl=TRUE),`[`,1)
} else {
  sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
}


#plotQualityProfile(fnFs[1:5])
#plotQualityProfile(fnRs[1:5])

filtFs <- file.path(resPath, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(resPath, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

if(cohort=="Angel" | cohort == "Angel_train" | cohort == "Angel_test"){
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(50, 55), truncLen=c(230, 200),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) 
} else if (cohort=="Antonio" | cohort == "Antonio_train" | cohort == "Antonio_test") {
   out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(58, 66), truncLen=c(220,150),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) 
} else if (cohort=="Gressel" | cohort == "Gressel_train" | cohort == "Gressel_test") {
   out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(19, 20), truncLen=c(220,170),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) 
} else if (cohort=="Tsementzi" | cohort == "Tsementzi_train" | cohort == "Tsementzi_test") {
   out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(19, 20), truncLen=c(220,170),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) 
} else if (cohort=="Walsh" | cohort == "Walsh_train" | cohort == "Walsh_test") {
   out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(51, 47), truncLen=c(220,190),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) 
} 

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
#plotErrors(errR, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, justConcatenate=TRUE)

seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
seqtab.nochim_useASV <- seqtab.nochim
ASV_towrite <- colnames(seqtab.nochim_useASV)
ASV_seqs_to_write <- rownames(seqtab.nochim_useASV)
colnames_replace <- paste0("ASV_", seq(1:length(colnames(seqtab.nochim_useASV))))
colnames(seqtab.nochim_useASV) <- colnames_replace
write.csv(seqtab.nochim_useASV, paste0(workspace, "/02-ASVs.csv"), row.names=TRUE)
write.csv(ASV_seqs_to_write, paste0(workspace, "/02-ASVsSeqs.csv"), row.names=TRUE)


for(i in seq(nrow(seqtab.nochim))){
  samples <- rownames(seqtab.nochim)  
  ids <- paste0(samples[i], "_", seq(rowSums(seqtab.nochim)[i]))
  seqs.re_replicated <- rep(colnames(seqtab.nochim), times=seqtab.nochim[i, ])
  writeFasta(object = ShortRead(sread = DNAStringSet(seqs.re_replicated),
                id = BStringSet(ids)), file = paste0(workspace, "/02-ASVsamples/", samples[i],".fasta"), width = 20000)
}
