###############################################
### Sepal morphology and gene expression
### Author: Diego A. Hartasanchez
### Created: 01/03/2022
### Last edit: 15/04/2022
### - R version used: 4.0.3
### Pipeline:
### - Generate Figure 4C network after running GENIE3 setup
### 

setwd('~/Documents/Academico/Proyectos/SepalGeneNetworks_2020/SepalMorphology_GitHubRepository/')

#Read ngs and gtm_magenta
ngs <- read.table("Results/normCounts_log2_Sepal_DEF_14085.txt",h=T)
gtm_magenta <- read.table("Data/GeneratedData/geneTraitModule_DEF_27samp_7param_14085gene_50min_.35thre_3_Width_inModule_magenta.txt",h=T)

#This list of cw genes in magenta is extracted from the list of genes with GO term cell wall organization or biogenesis associated to them
cwGenes_ngs <- read.table("Data/cw_genes_in_14085.txt",h=F,sep=";")
length(cwGenes_ngs)
gtm_magenta_cw <- gtm_magenta[gtm_magenta$gene %in% cwGenes_ngs,]
dim(gtm_magenta_cw)

ngs_magenta_cw <- ngs[rownames(ngs) %in% gtm_magenta_cw$gene,]
dim(ngs_magenta_cw)

expr.matrix <- ngs_magenta_cw
geneNames <- read.table('Data/geneNames.txt',h=F)
geneNames <- as.data.frame(geneNames)
colnames(geneNames)=c("Gene_ID","Gene_Name")
for(i in 1:nrow(geneNames)){
  if(geneNames$Gene_Name[i] == "araport11"){
    geneNames$Gene_Name[i] = geneNames$Gene_ID[i]
  }
}

dim(expr.matrix)
newNames <- geneNames[geneNames$Gene_ID %in% rownames(expr.matrix),2]
rownames(expr.matrix) <- newNames

weightMat <- GENIE3(expr.matrix, regulators = NULL, targets = NULL,treeMethod = "RF", K = "sqrt", nTrees = 1000, nCores = 1,verbose = TRUE)

threshLink = 0.045
e <- rf
e[e<threshLink]=0
e[e>=threshLink]=1
e_filtered <- e
cw_graph<-graph_from_adjacency_matrix(e_filtered,mode="undirected")

pdf(file='Plots/Figure4C_magenta_cw_grn_.045.pdf',height=6,width=12)
ggnet2(cw_graph, size = "degree", label = TRUE, label.size = 4, label.color="black",
       node.color="deepskyblue1", edge.color="grey",edge.size = .5)
dev.off()


