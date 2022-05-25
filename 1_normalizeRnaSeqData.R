###############################################
### Sepal morphology and gene expression
### Author: Diego A. Hartasanchez
### Created: 25/08/2021
### Last edit: 20/04/2022
### - R version used: 4.0.3
### Pipeline:
### - Data processing: normalization, quality control, filtering
### 


#########################################
##### Packages to be installed ##########
#########################################

#install.packages("devtools")
#install.packages("FactoMineR")
#install.packages("ggfortify")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR")
#install.packages("edgeR")
#install.packages("viridis")
#install.packages("factoextra")
#install.packages("GGally")
#install.packages("moments")
#install.packages("tidyverse")

####################
###Load libraries ##

library("devtools")
#devtools::install_github("jaredhuling/jcolors")
library("FactoMineR")
library(ggfortify)
library(ggplot2)
library(edgeR)
library(viridis)
library("factoextra")
library(RColorBrewer)
library(GGally)
library(grid)
library(moments)
library(tidyverse)
library(gridExtra)
library(lattice)
library(jcolors)
library(cowplot)

### Outside this pipeline:
### Pseudoalignement to the TAIR10 reference genome was performed with kallisto
### Transcript to gene transformation was performed with tximport

## Set working directory
setwd('~/Documents/Academico/Proyectos/SepalGeneNetworks_2020/SepalMorphology_GitHubRepository/')
## Load data. Gene counts from 30 sepal samples 
dataCountsDEF <- read.table('Data/Counts_kallisto_tximport_DEF.tab', h=T, stringsAsFactors=F, row.names = 1)
dim(dataCountsDEF)
head(dataCountsDEF)

### Analysis with DiCoExpress revealed three outlier samples

### Removal of outlier samples F10, F12 and F15 (columns 21, 23 and 25)
### and reordering of columns
dataCountsDEF_27ord <- dataCountsDEF[,c(4,5,6,7,8,9,10,1,2,3,
                                        16,17,18,19,20,11,12,13,14,15,
                                        27,28,29,30,22,24,26)]
head(dataCountsDEF_27ord)

#change sample names
colnames(dataCountsDEF_27ord) <- c("D1","D2","D3","D4","D6","D7","D8","D10","D11","D12",
                                   "E2","E3","E4","E6","E9","E10","E14","E15","E16","E17",
                                   "F3","F5","F7","F9","F11","F14","F17")

dataCounts <- dataCountsDEF_27ord

############################
### NORMALIZATION WITH EdgeR 
############################

x.norm.edger <- DGEList(counts=dataCounts)

# Remove genes with low expression
# Keep genes with at least 5 copies per million in at least 14 samples 
# 14 samples are chosen since they represent more than 50% of the remaining 27 sampless in each sample
keep=rowSums(cpm(x.norm.edger) >= 5) >=  14
fcounts <- x.norm.edger[keep,]
dim(fcounts)
### This filtering criteria reduces our dataset to 14085 genes to analyze

#normalization with TMM
fcounts$samples$lib.size=colSums(fcounts$counts)
normCounts <- calcNormFactors(fcounts, method="TMM")
tmm <- normCounts$samples$norm.factors
N <- normCounts$samples$lib.size
f <- tmm * N/mean(tmm * N) 

## normalized read counts are obtained by dividing raw read counts by these re-scaled normalization factors
NormCounts <- scale(normCounts$counts, center=FALSE, scale=f)

## Log2 + 1 transformation of NormCounts
NormCounts_log2 <- log2(NormCounts+1)
dim(NormCounts_log2)

## Write Normalized Counts for 27 samples and list of gene names to Results folder
write.table(NormCounts_log2,file="Results/normCounts_log2_Sepal_DEF_14085.txt",row.names = TRUE, quote = FALSE)
#write.table(rownames(NormCounts_log2),file = "Results/geneNames_Sepal_DEF_14085.txt",quote=FALSE,row.names = FALSE,col.names = FALSE)


###################################################
### PLOTS gene expression density, boxplots and PCA
###################################################
# -------------------------------------------------------------------------
##Plot density of gene expression and boxplots for all samples 
colorSamples=viridis(27)
pdf('Plots/Supp_geneExpressionPerSample.pdf', width=8, height=6)
plot(density(NormCounts_log2[,1]),xlim=c(3,18),ylim=c(0,0.3),col=colorSamples[1])
for(i in 2:27){
  lines(density(NormCounts_log2[,i]),xlim=c(0,10000),col=colorSamples[i])
}
boxplot(NormCounts_log2,cex=0.4,axes=FALSE,col=colorSamples)
axis(1,cex.axis=0.8,labels = colnames(NormCounts_log2), at =seq(1,27,1),las=2)
axis(2,cex.axis=0.8)
dev.off()

#Calculate PCA for gene expression and check for Plant effect
tNormCounts <- t(NormCounts_log2)
project.pca <- prcomp(tNormCounts, scale. = TRUE)

#Plant Data specifies which samples correspond to which plant (D, E and F)
plantData <- read.table('Data/sampleDataType_Sepals_DEF.txt',header=TRUE,row.names=1)
dim(plantData)

pca_geneExp <- fviz_pca_ind(project.pca, habillage=plantData$Plant,
                     addEllipses=TRUE, ellipse.level=0.90,repel=TRUE,
                   #  xlab="PC1 (56.5%)", ylab="PC2 (22.3%)",
                     title=NULL,legend="none"
) + 
  theme(legend.position = "none",
        plot.title = element_blank())


#Plot PCA
pdf('Plots/Supp_Figure3_PCAgeneExp.pdf', width=6, height=4)
  pca_geneExp
dev.off()


