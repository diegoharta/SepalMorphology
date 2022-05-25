###############################################
### Sepal morphology and gene expression
### Author: Diego A. Hartasanchez
### Created: 01/03/2022
### Last edit: 
### - R version used: 4.0.3
### Pipeline:
### - Generate Figure 4B after WGCNA output
### 

## Load libraries
library('ggpubr')
library('network')
library('sna')
library('igraph')
library('intergraph')


setwd('~/Documents/Academico/Proyectos/SepalGeneNetworks_2020/SepalMorphology_GitHubRepository/')
ngs <- read.table("Results/normCounts_log2_Sepal_DEF_14085.txt",h=T)
shape <- read.table("Results/shapeParams_7_pipe.txt",h=T)

moduleColors <- read.table("Data/GeneratedData/moduleColors_DEF_27samp_7param_14085gene_50min_.35thre_3_50.txt",h=F)
modCol_pipe3 <- cbind(rownames(ngs),moduleColors)
colnames(modCol_pipe3) <- c("Gene","Module")
write.table(modCol_pipe3, file="Results/geneModuleColors.txt",row.names = FALSE,col.names = TRUE, quote = FALSE, sep="\t")

tf_unique <- read.table("Data/tf_list_family_unique.txt",h=T)

gtm_magenta <- read.table("Data/GeneratedData/geneTraitModule_DEF_27samp_7param_14085gene_50min_.35thre_3_Width_inModule_magenta.txt",h=T)
gtm_magenta$absTraitSig <- abs(gtm_magenta$traitSignificance)
geneList <- gtm_magenta$gene

ngs_magenta <- ngs[rownames(ngs) %in% geneList,]
equis <- seq(1,27,1)

#Read list of CW genes
cwGenes_ngs <- read.table("Data/cw_genes_in_14085.txt",h=F,sep=";")
gtm_magenta_cw <- gtm_magenta[gtm_magenta$gene %in% cwGenes_ngs,]
#dim(gtm_magenta_cw)
ngs_magenta_cw <- ngs[rownames(ngs) %in% gtm_magenta_cw$gene,]
#dim(ngs_magenta_cw)

gtm_magenta_tf <- gtm_magenta[gtm_magenta$gene %in% tf_unique$Gene_ID,]
ngs_magenta_tf <- ngs[rownames(ngs) %in% gtm_magenta_tf$gene,]
ngs_magenta_tf[1,]
ngs_magenta_tf[4,]

pdf(file = "Plots/Figure4A&B_magentaExpressionAcrossSamples_width.pdf", wi = 12, he = 5)
layout(matrix(c(1,2,3,3), nr=2, byrow=T),widths = c(5,7), heights=c(4,.3))

par(mar=c(5,6,1,2))
#par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(abs(gtm_magenta$traitSignificance),gtm_magenta$moduleMembership,
     col="mediumorchid2",pch=20,ylab="Module membership",
     xlab="abs(trait significance) for Width ",cex=1.8,cex.lab=1.8,cex.axis=1.2,las=1)
abline(h=.7,col="black",lty=2)
# abline(v=.3,col="black",lty=2)
points(abs(gtm_magenta$traitSignificance),gtm_magenta$moduleMembership,
       col="mediumorchid2",pch=20,cex=1.8)
points(abs(gtm_magenta_cw$traitSignificance),gtm_magenta_cw$moduleMembership,
       col="deepskyblue1",pch=1,cex=2,lwd=2)
points(abs(gtm_magenta_tf[1,]$traitSignificance),gtm_magenta_tf[1,]$moduleMembership,
       col="firebrick1",pch=1,cex=3,lwd=3)
points(abs(gtm_magenta_tf[4,]$traitSignificance),gtm_magenta_tf[4,]$moduleMembership,
       col="goldenrod2",pch=1,cex=3,lwd=3)

par(mar=c(5,5,1,5))
adjustcolor("magenta", alpha = .2)
log( gtm_magenta$moduleMembership)
colorLine=adjustcolor("magenta", alpha = .4)
plot(equis, unname(ngs_magenta[1,]),col="white",
     lty=1,lwd=.4,#(.5+gtm_magenta$moduleMembership[1])^2-min((.5+gtm_magenta$moduleMembership[1])^2)+.1,
     type = "l",
     ylim = c(5,11.8),
     xaxt = "n",#xlab="Sample",
     xlab="",ylab="log2 gene expression",cex.lab=1.8,cex.axis=1.2,las=2)

for(i in 1:nrow(gtm_magenta_cw)){
  if(gtm_magenta_cw$moduleMembership[i] >= .7 && gtm_magenta_cw$absTraitSig[i] >= 0){
    colorLine=adjustcolor("deepskyblue1", alpha = gtm_magenta_cw$moduleMembership[i])
    lineWidth=(gtm_magenta_cw$moduleMembership[i]-.5)*5
    # lineWidth=1.5
    lines(equis, unname(ngs_magenta_cw[i,]),col=colorLine,lty=1,lwd=lineWidth)
  }
}
#NAC007
lineWidth=  (gtm_magenta_cw$moduleMembership[gtm_magenta_cw$gene == "AT1G12260"]-.5)*5
#  lineWidth=1.5
lines(equis, unname(ngs_magenta_tf[1,]),col="firebrick1",lty=1,lwd=lineWidth)
#KNAT7
lineWidth=  (gtm_magenta_cw$moduleMembership[gtm_magenta_cw$gene == "AT1G62990"]-.5)*5
#  lineWidth=1.5  
lines(equis, unname(ngs_magenta_tf[4,]),col="goldenrod2",lty=1,lwd=lineWidth)

par(new = TRUE)
plot(equis,1/shape$Width*1000,col="black", xaxt = "n", yaxt = "n",
     ylab = "", xlab = "",type="l",ylim=c(0.7,1.92),lwd=2,cex.axis=1.3)
axis(side = 4,cex.axis=1.2,las=1)
mtext("1000/Width", side = 4, line = 3.5,cex=1.5)
axis(1,label=colnames(ngs),at=equis,las=2,cex.axis=1.5,las=2)


par(mar=c(0,0,0,0))
plot.new()
legend("top", legend=c("Magenta genes","Cell wall genes","NAC007", "KNAT7","Width"), pch=c(20,1,1,1,NA),
       lty=c(NA,1,1,1,1),lwd=c(NA,1,2,2,4),pt.cex=c(2,2,3,3,2),pt.lwd=c(1,2,3,3,NA),col=c("mediumorchid2","deepskyblue1","firebrick1","goldenrod2","black"),
       bty="n",horiz = TRUE,cex=1.2,
       text.width=c(0.12,0.12,0.12,0.12,0.12),
       inset = c(0, -0.2), x.intersp=0.1,
       xjust=0, yjust=0)


dev.off()

#write.table(gtm_magenta$gene,file="Results/magentaGenes.txt",row.names = FALSE,
 #           col.names = FALSE,quote = FALSE)
