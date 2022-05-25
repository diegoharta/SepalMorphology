###############################################
### Sepal morphology and gene expression
### Author: Diego A. Hartasanchez
### Created: 25/08/2021
### Last edit: 20/04/2022
### - R version used: 4.0.3
### Pipeline step:
### - Generate Figure 5A
###   - GO terms for HVGs and LVGs obtained manually
### - Generarate Figure 5B

## Set working directory
setwd('~/Documents/Academico/Proyectos/SepalGeneNetworks_2020/SepalMorphology_GitHubRepository/')

#Load shape parameters
shapeParams <- read.table("Results/shapeParams_7_pipe.txt",h=T,row.names=1,sep="\t")

#Load plant Data specifies which samples correspond to which plant (D, E and F)
plantData <- read.table('Data/sampleDataType_Sepals_DEF.txt',header=TRUE,row.names=1)

#Define plant colors
plantColors <- c('red2','forestgreen','cornflowerblue','orange')

#read gene expression file
ngs <- read.table(file="Results/normCounts_log2_Sepal_DEF_14085.txt",row.names = 1, h=T)


##############################################################
### CV2 analysis on gene expression ##########################
##############################################################

#set length of genes to test
lengthOfTest <- nrow(ngs)
#lengthOfTest <- 2000

#construct data frame with info for each gene
#info: variance, means, cv2
geneCharac <- data.frame(vars[1:lengthOfTest],means[1:lengthOfTest],cv2[1:lengthOfTest])
colnames(geneCharac) <- c("Variance", "Means","CV2")

#This is figure 5A, with cv2 vs mean showing hvgs and lvgs

#define color palette
cv2Colors <- brewer.pal(n = 8, name = "BrBG")

#Add columng Variability
geneCharac <- geneCharac %>%
  add_column(Variability = "Other")
head(geneCharac)

geneCharacOrderedByCv2 <- geneCharac[order(-geneCharac$CV2),]
head(geneCharacOrderedByCv2)

sizeOfGroups <- 704
geneCharacOrderedByCv2[1:sizeOfGroups,4] = "Top 5%"
geneCharacOrderedByCv2[(sizeOfGroups+1):(sizeOfGroups*2),4] = "Top 5-10%"
geneCharacOrderedByCv2[(sizeOfGroups*2+1):(sizeOfGroups*3),4] = "Top 10-15%"

totNumOfGenes <- nrow(geneCharacOrderedByCv2)
geneCharacOrderedByCv2[(totNumOfGenes-sizeOfGroups+1):totNumOfGenes,4] = "Bottom 5%"
geneCharacOrderedByCv2[(totNumOfGenes-2*sizeOfGroups+1):(totNumOfGenes-sizeOfGroups),4] = "Bottom 5-10%"
geneCharacOrderedByCv2[(totNumOfGenes-3*sizeOfGroups+1):(totNumOfGenes-2*sizeOfGroups),4] = "Bottom 10-15%"

hvg_lvg_plot <- ggplot() +
  geom_point(
    data=geneCharacOrderedByCv2[geneCharacOrderedByCv2$Variability == "Top 5%",],
    aes(x=Means, y=CV2), 
    color = cv2Colors[1],
    size=.1) +
  geom_point(
    data=geneCharacOrderedByCv2[geneCharacOrderedByCv2$Variability == "Top 5-10%",],
    aes(x=Means, y=CV2), 
    color = cv2Colors[3],
    size=.1) +
  geom_point(
    data=geneCharacOrderedByCv2[geneCharacOrderedByCv2$Variability == "Top 10-15%",],
    aes(x=Means, y=CV2), 
    color = cv2Colors[4],
    size=.1) +
  geom_point(
    data=geneCharacOrderedByCv2[geneCharacOrderedByCv2$Variability == "Other",],
    aes(x=Means, y=CV2), 
    color = cv2Colors[5],
    size=.1) +
  geom_point(
    data=geneCharacOrderedByCv2[geneCharacOrderedByCv2$Variability == "Bottom 10-15%",],
    aes(x=Means, y=CV2), 
    color = cv2Colors[6],
    size=.1) +
  geom_point(
    data=geneCharacOrderedByCv2[geneCharacOrderedByCv2$Variability == "Bottom 5-10%",],
    aes(x=Means, y=CV2), 
    color = cv2Colors[7],
    size=.1) +
  geom_point(
    data=geneCharacOrderedByCv2[geneCharacOrderedByCv2$Variability == "Bottom 5%",],
    aes(x=Means, y=CV2), 
    color = cv2Colors[8],
    size=.1) +
  scale_x_continuous(trans = 'log10', limits = c(80,600000)) +
  scale_y_continuous(trans = 'log10') +
  labs(x= expression(paste(mu,sep="")), 
       y = expression(paste(CV^2,sep="")))

xdensity_hvg_lvg <- ggplot() +
  geom_density(
    data=geneCharacOrderedByCv2[geneCharacOrderedByCv2$Variability == "Bottom 5%",],
    aes(x=Means),fill=cv2Colors[8],alpha=.4) + 
  geom_density(
    data=geneCharacOrderedByCv2[geneCharacOrderedByCv2$Variability == "Bottom 5-10%",],
    aes(x=Means),fill=cv2Colors[7],alpha=.4) + 
  geom_density(
    data=geneCharacOrderedByCv2[geneCharacOrderedByCv2$Variability == "Bottom 10-15%",],
    aes(x=Means),fill=cv2Colors[6],alpha=.4) + 
  geom_density(
    data=geneCharacOrderedByCv2[geneCharacOrderedByCv2$Variability == "Other",],
    aes(x=Means),fill=cv2Colors[5],alpha=.4) + 
  geom_density(
    data=geneCharacOrderedByCv2[geneCharacOrderedByCv2$Variability == "Top 10-15%",],
    aes(x=Means),fill=cv2Colors[4],alpha=.4) + 
  geom_density(
    data=geneCharacOrderedByCv2[geneCharacOrderedByCv2$Variability == "Top 5-10%",],
    aes(x=Means),fill=cv2Colors[3],alpha=.4) + 
  geom_density(
    data=geneCharacOrderedByCv2[geneCharacOrderedByCv2$Variability == "Top 5%",],
    aes(x=Means),fill=cv2Colors[2],alpha=.4) + 
  scale_x_continuous(trans = 'log10',limits = c(80,600000)) +
  labs(y = "Density")+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = unit(c(.2,.2,.1,.465), "cm"))
xdensity_hvg_lvg


geneCharacOrderedByCv2$Variability <- factor(geneCharacOrderedByCv2$Variability, levels=c("Top 5%","Top 5-10%","Top 10-15%","Other","Bottom 10-15%","Bottom 5-10%","Bottom 5%"),
                                             labels=c("Top 5%","Top 5-10%","Top 10-15%","Other","Bottom 10-15%","Bottom 5-10%","Bottom 5%"))

hvg_lvg_legend <- ggplot(as.data.frame(geneCharacOrderedByCv2),aes(x=Means, y=CV2, color=Variability)) + 
  #geom_point() + 
  scale_color_manual(labels=c("Top 5%","Top 5-10%","Top 10-15%","Other","Bottom 10-15%","Bottom 5-10%","Bottom 5%"),
                     values = c(cv2Colors[2],cv2Colors[3],cv2Colors[4],cv2Colors[5],
                                cv2Colors[6],cv2Colors[7],cv2Colors[8])) + 
  # scale_color_brewer(palette="Dark2") + theme_minimal()+
  geom_point(shape=20, size=4, alpha=1)+
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  theme(legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10)) 
#  theme(legend.position=c(0,1), legend.justification=c(0,1))

hvg_legend <- cowplot::get_legend(hvg_lvg_legend)

blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  )


pdf("Plots/Figure5A_CV2_means_hvg_lvg_classification.pdf",width=8, height=4)
grid.arrange(xdensity_hvg_lvg, blankPlot, hvg_lvg_plot, hvg_legend, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.7, 4))
dev.off()


####
# From the list of HVGs and LVGs obtained above GO enrichment terms were obtained using PANTHER
# Using PANTHER's hierarchical term classification, a list of high-hierarchical
# terms was obtained and saved in to files with fold change enrichment and false
# discovery rate in the Folder Data/GeneratedData.


## Load libraries
library("RColorBrewer")
#install.packages("pheatmap")
library("pheatmap")


topHierarchy_FC <- as.data.frame(read.table('Data/GeneratedData/FC_HVG_LVG_goTerms_panther.txt',sep='\t',h=T,row.names = 1))
logFC <- log2(topHierarchy_FC) 

topHierarchy_FDR <- as.data.frame(read.table('Data/GeneratedData/FDR_HVG_LVG_goTerms_panther.txt',sep='\t',h=T,row.names = 1))
logFDR <- log10(topHierarchy_FDR) 

# Plot heatmaps for FC and FDR for the 6 variability groups and 14 top GO terms
pdf('Plots/Figure5B_topHierachy_GoFor6_HvgsLvgs_FDR_FC.pdf',h=3,w=6)
layout(matrix(c(1, 1,1,2,2), nr=1, byrow=T))
par(mar=c(4,1,3,15))
FC_t <- t(logFC)
colfunc <- colorRampPalette(c("white","navyblue"))
#colorFC <- brewer.pal(n = 8, name = "Reds")
#colfunc <- colorRampPalette(c("white","#4292C6", "#2171B5", "#084594"))
image(1:nrow(FC_t),1:ncol(FC_t), FC_t, axes=F, ylab="", xlab="",
      # col= colorRampPalette((brewer.pal(n = 7, name ="Blues")))(100)
      col=colfunc(80),
      main="log2 fold change"
)
abline(v=3.5,col="black",lwd=2)
axis(1, at=c(.5,1.5,2.5,3.5,4.5,5.5,6.5), labels=c("0","5","10","15","10","5","0"), 
     las=1,cex.axis=1)
#1:nrow(FC_t), rownames(FC_t), las=1,cex.axis=.7)
axis(4, 1:ncol(FC_t), colnames(FC_t), las=1,cex.axis = 1)
mtext(side = 1, at = 2, line = 2.5, las=1, cex=.7, text = "Top %")
mtext(side = 1, at = 5, line = 2.5, las=1, cex=.7, text = "Bottom %")

FDR_t <- t(-logFDR)
colfunc <- colorRampPalette(c("white","red","firebrick","firebrick4"))

par(mar=c(4,3,3,4))
image(1:nrow(FDR_t),1:ncol(FDR_t), FDR_t#-log10(FDR_t)
      , axes=F, ylab="", xlab="",
      main="-log10 FDR",
      col=colfunc(80)
      # col= colorRampPalette((brewer.pal(n = 7, name ="Reds")))(100)
)
#axis(1, 1:nrow(FDR_t), rownames(FDR_t), las=1,cex.axis=.7)
axis(2, 1:ncol(FDR_t),label=FALSE)

abline(v=3.5,col="black",lwd=2)
axis(1, at=c(.5,1.5,2.5,3.5,4.5,5.5,6.5), labels=c("0","5","10","15","10","5","0"), 
     las=1,cex.axis=1)
mtext(side = 1, at = 2, line = 2.5, las=1, cex=.7, text = "Top %")
mtext(side = 1, at = 5, line = 2.5, las=1, cex=.7, text = "Bottom %")
dev.off()


# Two additional plots to create legends for the above plot
pdf('Plots/Figure5B_FC_legend_topHierachy_GoFor6_HvgsLvgs.pdf',h=4,w=6)
colfunc <- colorRampPalette(c("white","navyblue"))
pheatmap::pheatmap(data.matrix(logFC),cluster_cols = FALSE,
                   color = colfunc(80))
dev.off()

pdf('Plots/Figure5B_FDR_legend_topHierachy_GoFor6_HvgsLvgs.pdf',h=4,w=6)
colfunc <- colorRampPalette(c("white","red","firebrick","firebrick4"))
#log10_fdr <- -log10(FDR_heatMapGo)
pheatmap::pheatmap(data.matrix(-logFDR),cluster_cols = FALSE,color=colfunc(80))
dev.off()





