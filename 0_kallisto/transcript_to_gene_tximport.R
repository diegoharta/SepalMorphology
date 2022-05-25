#!/usr/local/R-3.6.0/bin/Rscript


library(tximport)
library(rhdf5)


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
write.table(tx.imported$counts, file = paste(ResultsFolder,"/Counts_kallisto_tximport.tab", sep =""), quote = F, sep="\t", col.names = NA, row.names = TRUE)
write.table(tx.imported$abundance, file = paste(ResultsFolder,"/TPM_kallisto_tximport.tab", sep =""), quote = F, sep="\t", col.names = NA, row.names = TRUE)
write.table(tx.imported$length, file = paste(ResultsFolder,"/Length_kallisto_tximport.tab", sep =""), quote = F, sep="\t", col.names = NA, row.names = TRUE)










