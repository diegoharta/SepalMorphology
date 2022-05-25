#!/usr/local/R-3.6.0/bin/Rscript


library(tximport)
library(rhdf5)
library(DESeq2)
library(apeglm)


# Files, folders and parameters
args <- commandArgs(trailingOnly=TRUE)
Metadatafile <- args[1]
AnnotationsTrasncToGene <- args[2]
ResultsFolder <- args[3]

metadata <- read.table(Metadatafile, h=T)

system_out <- system(paste("zcat ", AnnotationsTrasncToGene," | tail -n +2 | sort -u", sep =""), intern=T)
gene.to.tx <- read.table(text=system_out, h=F)
tx2gene <- as.data.frame(cbind(as.character(gene.to.tx[,2]), as.character(gene.to.tx[,1])))
#head(tx2gene)

system_out <- system(paste("ls ", ResultsFolder,"/*/abundance.h5", sep =""), intern=T)
samplefiles <- read.table(text=system_out, h=F)
samplefiles <- as.character(samplefiles[,1])
names(samplefiles) <- metadata$ID
print(samplefiles)


tx.imported <- tximport(files=samplefiles, type="kallisto", tx2gene=tx2gene, ignoreTxVersion=FALSE,txIn=TRUE,txOut=FALSE)
write.table(tx.imported$counts, file = paste(ResultsFolder,"/Counts_kallisto_tximport2.tab", sep =""), quote = F, sep="\t", col.names = NA, row.names = TRUE)
write.table(tx.imported$abundance, file = paste(ResultsFolder,"/TPM_kallisto_tximport2.tab", sep =""), quote = F, sep="\t", col.names = NA, row.names = TRUE)
write.table(tx.imported$length, file = paste(ResultsFolder,"/Length_kallisto_tximport2.tab", sep =""), quote = F, sep="\t", col.names = NA, row.names = TRUE)


coldataAll<-read.table('colData_60Samples.txt', h=T, stringsAsFactors=T)

Data <- tximport(files=samplefiles, type="kallisto", tx2gene=tx2gene, ignoreTxVersion=FALSE,txIn=TRUE,txOut=FALSE)
dds <- DESeqDataSetFromTximport(txi = Data, colData=coldataAll, design=~ temperature)
test <- DESeq(dds, test='Wald')
results.s <- lfcShrink(test, coef=resultsNames(test)[2], type="apeglm")
results <- results(test)
shortlab<-"60Samples"
write.table(results.s, file = paste(ResultsFolder, "/", shortlab, "_Shrink.txt", sep=""), quote = F, sep="\t", col.names = NA, row.names = TRUE)
write.table(results, file = paste(ResultsFolder, "/", shortlab, ".txt", sep=""), quote = F, sep="\t", col.names = NA, row.names = TRUE)

results.All=results(test)
results.All=results.All[order(results.All$padj),]

pdf(paste(ResultsFolder,"/volcanoPlot_DiffExp_TPM_60Samples.pdf", sep =""), width=10, height=5)
par(mfrow=c(1,1))
plot(results.All$log2FoldChange, -log10(results.All$padj), xlab='log2FoldChange', ylab='-log10 FDR',
     main='60Samples_TPM', pch=20, xlim=c(-5, 5), ylim=c(0, 200))
dev.off()


