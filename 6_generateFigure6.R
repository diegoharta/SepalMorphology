###############################################
### Sepal morphology and gene expression
### Author: Diego A. Hartasanchez
### Created: 20/12/2021
### Last edit: 04/04/2022
### - R version used: 4.0.3
### Pipeline:
###   R script to generate Figure 6 
###   CWRG list and CW_GO CV2 compared to all genes
###   enrichment of cw genes in each module 

library(ggpval)
library(tidyverse)
library(viridis)
library(ggpubr)

### Comparison between cell wall organization and biogenesis go term with cell-wall related gene list
setwd('~/Documents/Academico/Proyectos/SepalGeneNetworks_2020/SepalMorphology_GitHubRepository/')

# load data
cellWallGenes_Charline <- read.table('Data/cellWallGenesList_1585.txt',h=F)
ngs <- read.table(file="Results/normCounts_log2_Sepal_DEF_14085.txt",row.names = 1, h=T)
geneNames <- read.table('Data/geneNames.txt',h=F)
allGenes <- rownames(ngs)
CharlineIn14085 <- cellWallGenes_Charline[cellWallGenes_Charline[,1] %in% allGenes,1]
length(CharlineIn14085) #718 of Charlines list in 14085

###calculate CV2 for every gene
ed<-2^ngs
means <- rowMeans(ed)
vars <- apply(ed,1,var)
cv2 <- vars/means^2

##############################################################
### CV2 analysis on gene expression ##########################
##############################################################

#set length of genes to test
lengthOfTest <- nrow(ngs)

#construct geneCharac data frame with info for each gene
#info: variance, means, cv2
geneCharac <- data.frame(vars[1:lengthOfTest],means[1:lengthOfTest],cv2[1:lengthOfTest])
colnames(geneCharac) <- c("Variance", "Means","CV2")
head(geneCharac)

#Assign "All_non_CW" as a Type. 
geneCharac <- geneCharac %>%
  add_column(Type = "All_genes")
head(geneCharac)

geneCharac[rownames(geneCharac) %in% CharlineIn14085,4]="CW_list"
cw_list_genes <- geneCharac[geneCharac$Type == "CW_list",]

geneCharac2 <- data.frame(vars[1:lengthOfTest],means[1:lengthOfTest],cv2[1:lengthOfTest])
colnames(geneCharac2) <- c("Variance", "Means","CV2")
geneCharac2 <- geneCharac2 %>% add_column(Type = "All genes")

geneCharac <- data.frame(vars[1:lengthOfTest],means[1:lengthOfTest],cv2[1:lengthOfTest])
colnames(geneCharac) <- c("Variance", "Means","CV2")
geneCharac <- geneCharac %>%  add_column(Type = "All genes")

geneCharac[rownames(geneCharac) %in% CharlineIn14085,4]="In CWRG list"
cw_list_genes <- geneCharac[geneCharac$Type == "In CWRG list",]

geneCharac2 <- rbind(geneCharac2, cw_list_genes)
dim(geneCharac2)

geneCharac2$Type <- factor(geneCharac2$Type,
                           levels = c('All genes','In CWRG list'),ordered = TRUE)

my_comparisons <- list( c('All genes','In CWRG list'))

#pp+scale_color_manual(values=c("#440154FF","#2A788EFF"))
#pp + scale_fill_brewer(palette="Set1") + theme_classic()
theme_set(theme_bw())

pdf(file='Plots/Figure6B_Violin_cw_list_vs_all.pdf',height=5.2,width=4)
ggplot(geneCharac2, aes(x=Type, y=CV2, fill=Type)) + 
  geom_violin()+
  scale_y_continuous(trans = 'log10')+
  geom_boxplot(width=0.2, color="grey", alpha=0.15) +
  labs(y = expression(paste(CV^2,sep="")))+#, 
  #   x = "Classification according to presence in cell wall related lists") + 
  scale_fill_manual(values=c("#440154FF","#2A788EFF"))+
  # scale_fill_viridis(discrete = TRUE) +
  annotate("text", x=1, y=6.5, label= nrow(geneCharac2[geneCharac2$Type == "All genes",])) +
  annotate("text", x=2, y=6.5, label= nrow(geneCharac2[geneCharac2$Type == "In CWRG list",]))+
  stat_compare_means(label.y =c(1),label.x = 1)
dev.off()


##### Evaluate how many cell-wall related genes are in each module

moduleColors <- read.table("Data/GeneratedData/moduleColors_DEF_27samp_7param_14085gene_50min_.35thre_3_50.txt",h=F)
geneModule <- as.data.frame(cbind(allGenes,moduleColors[,1]))
head(geneModule)
colnames(geneModule) <- c("Gene","Module")
table(geneModule$Module)

dim(geneCharac)
geneCharac$Module = geneModule$Module

tbl <- with(geneCharac, table(Module, Type))
tbl_df <- as.data.frame(tbl)

cwList_perModule <- as.data.frame(tbl)[seq(17,32,1),c(1,3)]
expectedRate <- 718/14085
enrichment_Cw_list <- cwList_perModule$Freq/unname(table(geneModule$Module))
enrich_cw_list <- enrichment_Cw_list/expectedRate
names(enrich_cw_list) <- names(table(geneModule$Module))

enrich_cw_df_2 <- as.data.frame(enrich_cw_list)
colnames(enrich_cw_df_2) <- c("Module","Enrichment")
positions <- enrich_cw_df_2[c(4,11,16,13,14,5,1,10,9,12,8,3,15,2,6,7),1]



pdf(file='Plots/Figure6A_Barplot_cwList_enrichment_perModule_newnew.pdf',height=5.2,width=4)

ggplot(enrich_cw_df_2, aes(y=Enrichment, x=Module)) + 
  geom_bar(position="dodge", stat="identity", fill = "#2A788EFF")+ 
  labs(y = "Enrichment in CWRG list")+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  geom_abline(slope=0, intercept=1,  col = "black",lty=2)+
  scale_x_discrete(limits = positions)

dev.off()
