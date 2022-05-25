###############################################
### Sepal morphology and gene expression
### Author: Diego A. Hartasanchez
### Created: 12/10/2021
### Last edit: 04/04/2022
### - R version used: 4.0.3
### Pipeline:
### - WGCNA to be run in a cluster
### Based on WGCNA https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html  

#Set R path
#export PATH=/usr/local/R-3.6.0/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/snap/bin

setwd('~/WGCNA_pipe/')
experiment <- "DEF_27samp_7param_14085gene_50min_.35thre_pipe"

#Read data files
NGS_data <- read.table(file='normCounts_log2_Sepal_DEF_14085.txt',header=TRUE,row.names = 1)
morph_data  <- read.table(file='shapeParams_7_pipe.txt',header = TRUE, sep="\t",row.names = 1)
head(morph_data)
head(NGS_data)
dim(morph_data)
dim(NGS_data)

library(WGCNA);

options(stringsAsFactors = FALSE);
enableWGCNAThreads()

datExpr0 = as.data.frame(t(NGS_data));
rownames(datExpr0) = names(NGS_data);

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
dim(datExpr0)

sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
pdf(file = paste("sampleClustering_",experiment,".pdf",sep=""), width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)    
dev.off()

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)


datSamples <- as.data.frame(morph_data)
collectGarbage();

sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datSamples, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
pdf(file = paste("sampleDendogram_",experiment,"2.pdf",sep=""), width = 12, height = 9);
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datSamples),
                    cex.dendroLabels = 1.8, cex.colorLabels = 1.8,
                    marAll = c(1, 8, 3, 1),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

save(datExpr0, datSamples, file = "MorphDataAll-01-dataInput.RData")


## part 2

options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. 
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "MorphDataAll-01-dataInput.RData");
# The variable lnames contains the names of loaded variables.
lnames

# Choose a set of soft-thresholding powers
powers = c(c(1:15), seq(from = 16, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
# Plot the results:
pdf(file = paste("scaleFreeTopology_",experiment,".pdf",sep=""), width = 20, height = 15);
sizeGrWindow(20, 15)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
#dev.off()
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

softPower = 11;
adjacency = adjacency(datExpr0, power = softPower);
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
pdf(file = paste("geneClustering_",experiment,"_softPower",softPower,".pdf",sep=""), width = 12, height = 9);
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()

# Set minimum modules size:
minModuleSize = 50;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf(file = paste("geneClustering_",experiment,"_modSize",minModuleSize,".pdf",sep=""), width = 8, height = 9);
#par(mfrow = c(1,1));
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

# Plot the result
pdf(file = paste("geneClustering_",experiment,".pdf",sep=""), width = 8, height = 9);
sizeGrWindow(8, 9)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")


#Set threshold to merge modules that are close
MEDissThres = 0.35
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()
# Call an automatic merging function
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


pdf(file = paste("geneDendro_",experiment,"_",minModuleSize,".pdf",sep=""), wi = 9, he = 6)
sizeGrWindow(9, 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
newModuleColors = moduleColors
#Replace colors by standard colors
newColors = c("black","brown","cyan","green","gold","magenta","red","purple",
              "grey","salmon","maroon","blue","orange","navy","violet","khaki");
for(i in 1:length(table(moduleColors))){
  newModuleColors[which(newModuleColors==names(table(moduleColors))[i])]=newColors[i]
}
moduleColors = newModuleColors
colorOrder = c("red", "khaki","purple",
               "cyan","grey","navy","black","salmon","magenta","orange",
               "gold","brown","green","maroon","violet","blue")

# Construct numerical labels corresponding to the colors
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = paste("sepals-",experiment,"_minModuleSize",minModuleSize,"-networkConstruction-stepByStep.RData",sep=""))

write.table(moduleColors, file = paste("moduleColors_",experiment,"_",minModuleSize,".txt",sep=""),
            row.names = FALSE, col.names = FALSE, quote = FALSE);

geneModule <- cbind(rownames(NGS_data),moduleColors)
colnames(geneModule) <- c("Gene_ID","Module")
write.table(geneModule, file = paste("geneModuleColors_",experiment,"_",minModuleSize,".txt",sep=""),
            row.names = FALSE, col.names = TRUE, quote = FALSE,sep="\t");


pdf(file = paste("moduleTraitHeatMap_",experiment,"_",minModuleSize,"_2b.pdf",sep=""), wi = 8, he = 6)
labeledHeatmap(Matrix = (as.matrix(table(moduleColors))),
               yLabels = rownames(as.matrix(table(moduleColors))),
               xLabels = "Size",
               # ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               #    colors=c("#543005", "#8C510A" ,"#BF812D" ,"#DFC27D", "#F6E8C3", "#F5F5F5" ,"#C7EAE5", "#80CDC1", "#35978F","#01665E", "#003C30"),
               # colors=c("#F7FCFD" ,"#E0ECF4", "#BFD3E6","#9EBCDA", "#8C96C6", "#8C6BB1" ,"#88419D" ,"#6E016B"),
               textMatrix = unname(table(moduleColors)),
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(0,3000),
               main = paste("Module sizes"))

### part 3

lnames = load(file = "MorphDataAll-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = paste("sepals-",experiment,"_minModuleSize",minModuleSize,"-networkConstruction-stepByStep.RData",sep=""));
lnames

datExpr=datExpr0
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datSamples, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

datTraits=datSamples
sizeGrWindow(8,6)
pdf(file = paste("moduleTraitHeatMap_",experiment,"_",minModuleSize,"_2b.pdf",sep=""), wi = 8, he = 6)
# Will display correlations and their p-values
pvalAsterisk <- symnum(moduleTraitPvalue, corr = FALSE, na = FALSE, 
                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                       symbols = c("***", "**", "*", ".", " "))
textMatrix = paste(signif(moduleTraitCor, 2), " ",
                   pvalAsterisk, sep = "");
#textMatrix= signif(moduleTraitPvalue, 1);
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors=c("#8C510A" ,"#BF812D" ,"#DFC27D", "#F6E8C3", "#F5F5F5" ,"#C7EAE5", "#80CDC1", "#35978F","#01665E"),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


##################

intModules = names(table(moduleColors))

for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  probes=names(datExpr)
  caseGenes = probes[modGenes]
  fileName = paste(module,"_GenesForGO_",experiment,".txt", sep="");
  write.table(as.data.frame(caseGenes), file = fileName,
              row.names = FALSE, col.names = FALSE,quote = FALSE)
}


probes=names(datExpr)
write.table(as.data.frame(probes),file = "allGenesForGO.txt",row.names=FALSE,
            col.names=FALSE,quote = FALSE)


##########################################

# Define variable weight containing the weight column of datTrait
traitInterest = as.data.frame(datTraits$Length);
names(traitInterest) = "Length"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, traitInterest, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(traitInterest), sep="");
names(GSPvalue) = paste("p.GS.", names(traitInterest), sep="");



for (module in intModules)
{
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  pdf(file = paste("scatterPlot_",experiment,"_",module,"_",names(traitInterest),".pdf",sep=""), wi = 7, he = 7)
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("Gene significance for ",names(traitInterest)),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  dev.off()
  
  
  ###Print to file geneTraitSignificance of trait of interest
  fileName = paste("traitSignificanceList_",experiment,"_",names(traitInterest),".txt", sep="");
  write.table(as.data.frame(geneTraitSignificance), file = fileName,
              row.names = TRUE, col.names = TRUE, quote = FALSE);
  
  ##order genes in a particular module by their geneTraitSignificance of a particular trait
  absGeneTraitSig <- abs(geneTraitSignificance)
  rownames(absGeneTraitSig)<- rownames(geneTraitSignificance)
  absGeneTraitSig$gene <- rownames(absGeneTraitSig)
  geneTraitSig_ord <- absGeneTraitSig[order(-absGeneTraitSig[,1]),] ### here I have genes ordered by their geneTraitSig value
  selectedDatExpr <- datExpr[,moduleGenes]
  signifList <- colnames(selectedDatExpr) ## this is the list of genes in the module
  geneTraitSigOrd_inModule <- geneTraitSig_ord[geneTraitSig_ord$gene %in% signifList,] ##list of genes in module ordered by GenetraitSignificance
  fileName = paste("geneTraitSigOrd_Trait_",experiment,"_",names(traitInterest),"_inModule_",module,".txt", sep="");
  write.table(as.data.frame(geneTraitSigOrd_inModule), file = fileName,
              row.names = FALSE, col.names = TRUE, quote = FALSE);
  
  ##print to file genemembers, trait significance and module membership
  absModuleMembershipSelec <- abs(geneModuleMembership[moduleGenes, column])
  geneTraitModule <- cbind(geneTraitSignificance[moduleGenes,], absModuleMembershipSelec)
  selectedDatExpr <- datExpr[,moduleGenes]
  signifList <- colnames(selectedDatExpr) ## this is the list of genes in the module
  geneTraitModule <- cbind(signifList,geneTraitModule)
  colnames(geneTraitModule) =c("gene","traitSignificance","moduleMembership")
  fileName = paste("geneTraitModule_",experiment,"_",names(traitInterest),"_inModule_",module,".txt", sep="");
  write.table(as.data.frame(geneTraitModule), file = fileName,
              row.names = FALSE, col.names = TRUE, quote = FALSE,sep="\t");
  
}

###############################################


# Define variable weight containing the weight column of datTrait
traitInterest = as.data.frame(datTraits$Width);
names(traitInterest) = "Width"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, traitInterest, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(traitInterest), sep="");
names(GSPvalue) = paste("p.GS.", names(traitInterest), sep="");



for (module in intModules)
{
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  pdf(file = paste("scatterPlot_",experiment,"_",module,"_",names(traitInterest),".pdf",sep=""), wi = 7, he = 7)
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("Gene significance for ",names(traitInterest)),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  dev.off()
  
  
  ###Print to file geneTraitSignificance of trait of interest
  fileName = paste("traitSignificanceList_",experiment,"_",names(traitInterest),".txt", sep="");
  write.table(as.data.frame(geneTraitSignificance), file = fileName,
              row.names = TRUE, col.names = TRUE, quote = FALSE);
  
  ##order genes in a particular module by their geneTraitSignificance of a particular trait
  absGeneTraitSig <- abs(geneTraitSignificance)
  rownames(absGeneTraitSig)<- rownames(geneTraitSignificance)
  absGeneTraitSig$gene <- rownames(absGeneTraitSig)
  geneTraitSig_ord <- absGeneTraitSig[order(-absGeneTraitSig[,1]),] ### here I have genes ordered by their geneTraitSig value
  selectedDatExpr <- datExpr[,moduleGenes]
  signifList <- colnames(selectedDatExpr) ## this is the list of genes in the module
  geneTraitSigOrd_inModule <- geneTraitSig_ord[geneTraitSig_ord$gene %in% signifList,] ##list of genes in module ordered by GenetraitSignificance
  fileName = paste("geneTraitSigOrd_Trait_",experiment,"_",names(traitInterest),"_inModule_",module,".txt", sep="");
  write.table(as.data.frame(geneTraitSigOrd_inModule), file = fileName,
              row.names = FALSE, col.names = TRUE, quote = FALSE);
  
  ##print to file genemembers, trait significance and module membership
  absModuleMembershipSelec <- abs(geneModuleMembership[moduleGenes, column])
  geneTraitModule <- cbind(geneTraitSignificance[moduleGenes,], absModuleMembershipSelec)
  selectedDatExpr <- datExpr[,moduleGenes]
  signifList <- colnames(selectedDatExpr) ## this is the list of genes in the module
  geneTraitModule <- cbind(signifList,geneTraitModule)
  colnames(geneTraitModule) =c("gene","traitSignificance","moduleMembership")
  fileName = paste("geneTraitModule_",experiment,"_",names(traitInterest),"_inModule_",module,".txt", sep="");
  write.table(as.data.frame(geneTraitModule), file = fileName,
              row.names = FALSE, col.names = TRUE, quote = FALSE,sep="\t");
  
}


# Define variable weight containing the weight column of datTrait
traitInterest = as.data.frame(datTraits$AspRatio);
names(traitInterest) = "AspRatio"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, traitInterest, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(traitInterest), sep="");
names(GSPvalue) = paste("p.GS.", names(traitInterest), sep="");



for (module in intModules)
{
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  pdf(file = paste("scatterPlot_",experiment,"_",module,"_",names(traitInterest),".pdf",sep=""), wi = 7, he = 7)
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("Gene significance for ",names(traitInterest)),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  dev.off()
  
  
  ###Print to file geneTraitSignificance of trait of interest
  fileName = paste("traitSignificanceList_",experiment,"_",names(traitInterest),".txt", sep="");
  write.table(as.data.frame(geneTraitSignificance), file = fileName,
              row.names = TRUE, col.names = TRUE, quote = FALSE);
  
  ##order genes in a particular module by their geneTraitSignificance of a particular trait
  absGeneTraitSig <- abs(geneTraitSignificance)
  rownames(absGeneTraitSig)<- rownames(geneTraitSignificance)
  absGeneTraitSig$gene <- rownames(absGeneTraitSig)
  geneTraitSig_ord <- absGeneTraitSig[order(-absGeneTraitSig[,1]),] ### here I have genes ordered by their geneTraitSig value
  selectedDatExpr <- datExpr[,moduleGenes]
  signifList <- colnames(selectedDatExpr) ## this is the list of genes in the module
  geneTraitSigOrd_inModule <- geneTraitSig_ord[geneTraitSig_ord$gene %in% signifList,] ##list of genes in module ordered by GenetraitSignificance
  fileName = paste("geneTraitSigOrd_Trait_",experiment,"_",names(traitInterest),"_inModule_",module,".txt", sep="");
  write.table(as.data.frame(geneTraitSigOrd_inModule), file = fileName,
              row.names = FALSE, col.names = TRUE, quote = FALSE);
  
  ##print to file genemembers, trait significance and module membership
  absModuleMembershipSelec <- abs(geneModuleMembership[moduleGenes, column])
  geneTraitModule <- cbind(geneTraitSignificance[moduleGenes,], absModuleMembershipSelec)
  selectedDatExpr <- datExpr[,moduleGenes]
  signifList <- colnames(selectedDatExpr) ## this is the list of genes in the module
  geneTraitModule <- cbind(signifList,geneTraitModule)
  colnames(geneTraitModule) =c("gene","traitSignificance","moduleMembership")
  fileName = paste("geneTraitModule_",experiment,"_",names(traitInterest),"_inModule_",module,".txt", sep="");
  write.table(as.data.frame(geneTraitModule), file = fileName,
              row.names = FALSE, col.names = TRUE, quote = FALSE,sep="\t");
  
}



# Define variable weight containing the weight column of datTrait
traitInterest = as.data.frame(datTraits$Area);
names(traitInterest) = "Area"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, traitInterest, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(traitInterest), sep="");
names(GSPvalue) = paste("p.GS.", names(traitInterest), sep="");



for (module in intModules)
{
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  pdf(file = paste("scatterPlot_",experiment,"_",module,"_",names(traitInterest),".pdf",sep=""), wi = 7, he = 7)
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("Gene significance for ",names(traitInterest)),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  dev.off()
  
  
  ###Print to file geneTraitSignificance of trait of interest
  fileName = paste("traitSignificanceList_",experiment,"_",names(traitInterest),".txt", sep="");
  write.table(as.data.frame(geneTraitSignificance), file = fileName,
              row.names = TRUE, col.names = TRUE, quote = FALSE);
  
  ##order genes in a particular module by their geneTraitSignificance of a particular trait
  absGeneTraitSig <- abs(geneTraitSignificance)
  rownames(absGeneTraitSig)<- rownames(geneTraitSignificance)
  absGeneTraitSig$gene <- rownames(absGeneTraitSig)
  geneTraitSig_ord <- absGeneTraitSig[order(-absGeneTraitSig[,1]),] ### here I have genes ordered by their geneTraitSig value
  selectedDatExpr <- datExpr[,moduleGenes]
  signifList <- colnames(selectedDatExpr) ## this is the list of genes in the module
  geneTraitSigOrd_inModule <- geneTraitSig_ord[geneTraitSig_ord$gene %in% signifList,] ##list of genes in module ordered by GenetraitSignificance
  fileName = paste("geneTraitSigOrd_Trait_",experiment,"_",names(traitInterest),"_inModule_",module,".txt", sep="");
  write.table(as.data.frame(geneTraitSigOrd_inModule), file = fileName,
              row.names = FALSE, col.names = TRUE, quote = FALSE);
  
  ##print to file genemembers, trait significance and module membership
  absModuleMembershipSelec <- abs(geneModuleMembership[moduleGenes, column])
  geneTraitModule <- cbind(geneTraitSignificance[moduleGenes,], absModuleMembershipSelec)
  selectedDatExpr <- datExpr[,moduleGenes]
  signifList <- colnames(selectedDatExpr) ## this is the list of genes in the module
  geneTraitModule <- cbind(signifList,geneTraitModule)
  colnames(geneTraitModule) =c("gene","traitSignificance","moduleMembership")
  fileName = paste("geneTraitModule_",experiment,"_",names(traitInterest),"_inModule_",module,".txt", sep="");
  write.table(as.data.frame(geneTraitModule), file = fileName,
              row.names = FALSE, col.names = TRUE, quote = FALSE,sep="\t");
  
}


# Define variable weight containing the weight column of datTrait
traitInterest = as.data.frame(datTraits$LongCurv);
names(traitInterest) = "LongCurv"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, traitInterest, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(traitInterest), sep="");
names(GSPvalue) = paste("p.GS.", names(traitInterest), sep="");



for (module in intModules)
{
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  pdf(file = paste("scatterPlot_",experiment,"_",module,"_",names(traitInterest),".pdf",sep=""), wi = 7, he = 7)
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("Gene significance for ",names(traitInterest)),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  dev.off()
  
  
  ###Print to file geneTraitSignificance of trait of interest
  fileName = paste("traitSignificanceList_",experiment,"_",names(traitInterest),".txt", sep="");
  write.table(as.data.frame(geneTraitSignificance), file = fileName,
              row.names = TRUE, col.names = TRUE, quote = FALSE);
  
  ##order genes in a particular module by their geneTraitSignificance of a particular trait
  absGeneTraitSig <- abs(geneTraitSignificance)
  rownames(absGeneTraitSig)<- rownames(geneTraitSignificance)
  absGeneTraitSig$gene <- rownames(absGeneTraitSig)
  geneTraitSig_ord <- absGeneTraitSig[order(-absGeneTraitSig[,1]),] ### here I have genes ordered by their geneTraitSig value
  selectedDatExpr <- datExpr[,moduleGenes]
  signifList <- colnames(selectedDatExpr) ## this is the list of genes in the module
  geneTraitSigOrd_inModule <- geneTraitSig_ord[geneTraitSig_ord$gene %in% signifList,] ##list of genes in module ordered by GenetraitSignificance
  fileName = paste("geneTraitSigOrd_Trait_",experiment,"_",names(traitInterest),"_inModule_",module,".txt", sep="");
  write.table(as.data.frame(geneTraitSigOrd_inModule), file = fileName,
              row.names = FALSE, col.names = TRUE, quote = FALSE);
  
  ##print to file genemembers, trait significance and module membership
  absModuleMembershipSelec <- abs(geneModuleMembership[moduleGenes, column])
  geneTraitModule <- cbind(geneTraitSignificance[moduleGenes,], absModuleMembershipSelec)
  selectedDatExpr <- datExpr[,moduleGenes]
  signifList <- colnames(selectedDatExpr) ## this is the list of genes in the module
  geneTraitModule <- cbind(signifList,geneTraitModule)
  colnames(geneTraitModule) =c("gene","traitSignificance","moduleMembership")
  fileName = paste("geneTraitModule_",experiment,"_",names(traitInterest),"_inModule_",module,".txt", sep="");
  write.table(as.data.frame(geneTraitModule), file = fileName,
              row.names = FALSE, col.names = TRUE, quote = FALSE,sep="\t");
  
}


# Define variable weight containing the weight column of datTrait
traitInterest = as.data.frame(datTraits$TransCurv);
names(traitInterest) = "TransCurv"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, traitInterest, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(traitInterest), sep="");
names(GSPvalue) = paste("p.GS.", names(traitInterest), sep="");



for (module in intModules)
{
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  pdf(file = paste("scatterPlot_",experiment,"_",module,"_",names(traitInterest),".pdf",sep=""), wi = 7, he = 7)
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("Gene significance for ",names(traitInterest)),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  dev.off()
  
  
  ###Print to file geneTraitSignificance of trait of interest
  fileName = paste("traitSignificanceList_",experiment,"_",names(traitInterest),".txt", sep="");
  write.table(as.data.frame(geneTraitSignificance), file = fileName,
              row.names = TRUE, col.names = TRUE, quote = FALSE);
  
  ##order genes in a particular module by their geneTraitSignificance of a particular trait
  absGeneTraitSig <- abs(geneTraitSignificance)
  rownames(absGeneTraitSig)<- rownames(geneTraitSignificance)
  absGeneTraitSig$gene <- rownames(absGeneTraitSig)
  geneTraitSig_ord <- absGeneTraitSig[order(-absGeneTraitSig[,1]),] ### here I have genes ordered by their geneTraitSig value
  selectedDatExpr <- datExpr[,moduleGenes]
  signifList <- colnames(selectedDatExpr) ## this is the list of genes in the module
  geneTraitSigOrd_inModule <- geneTraitSig_ord[geneTraitSig_ord$gene %in% signifList,] ##list of genes in module ordered by GenetraitSignificance
  fileName = paste("geneTraitSigOrd_Trait_",experiment,"_",names(traitInterest),"_inModule_",module,".txt", sep="");
  write.table(as.data.frame(geneTraitSigOrd_inModule), file = fileName,
              row.names = FALSE, col.names = TRUE, quote = FALSE);
  
  ##print to file genemembers, trait significance and module membership
  absModuleMembershipSelec <- abs(geneModuleMembership[moduleGenes, column])
  geneTraitModule <- cbind(geneTraitSignificance[moduleGenes,], absModuleMembershipSelec)
  selectedDatExpr <- datExpr[,moduleGenes]
  signifList <- colnames(selectedDatExpr) ## this is the list of genes in the module
  geneTraitModule <- cbind(signifList,geneTraitModule)
  colnames(geneTraitModule) =c("gene","traitSignificance","moduleMembership")
  fileName = paste("geneTraitModule_",experiment,"_",names(traitInterest),"_inModule_",module,".txt", sep="");
  write.table(as.data.frame(geneTraitModule), file = fileName,
              row.names = FALSE, col.names = TRUE, quote = FALSE,sep="\t");
  
}


# Define variable weight containing the weight column of datTrait
traitInterest = as.data.frame(datTraits$GaussCurv);
names(traitInterest) = "GaussCurv"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, traitInterest, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(traitInterest), sep="");
names(GSPvalue) = paste("p.GS.", names(traitInterest), sep="");



for (module in intModules)
{
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  pdf(file = paste("scatterPlot_",experiment,"_",module,"_",names(traitInterest),".pdf",sep=""), wi = 7, he = 7)
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("Gene significance for ",names(traitInterest)),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  dev.off()
  
  
  ###Print to file geneTraitSignificance of trait of interest
  fileName = paste("traitSignificanceList_",experiment,"_",names(traitInterest),".txt", sep="");
  write.table(as.data.frame(geneTraitSignificance), file = fileName,
              row.names = TRUE, col.names = TRUE, quote = FALSE);
  
  ##order genes in a particular module by their geneTraitSignificance of a particular trait
  absGeneTraitSig <- abs(geneTraitSignificance)
  rownames(absGeneTraitSig)<- rownames(geneTraitSignificance)
  absGeneTraitSig$gene <- rownames(absGeneTraitSig)
  geneTraitSig_ord <- absGeneTraitSig[order(-absGeneTraitSig[,1]),] ### here I have genes ordered by their geneTraitSig value
  selectedDatExpr <- datExpr[,moduleGenes]
  signifList <- colnames(selectedDatExpr) ## this is the list of genes in the module
  geneTraitSigOrd_inModule <- geneTraitSig_ord[geneTraitSig_ord$gene %in% signifList,] ##list of genes in module ordered by GenetraitSignificance
  fileName = paste("geneTraitSigOrd_Trait_",experiment,"_",names(traitInterest),"_inModule_",module,".txt", sep="");
  write.table(as.data.frame(geneTraitSigOrd_inModule), file = fileName,
              row.names = FALSE, col.names = TRUE, quote = FALSE);
  
  ##print to file genemembers, trait significance and module membership
  absModuleMembershipSelec <- abs(geneModuleMembership[moduleGenes, column])
  geneTraitModule <- cbind(geneTraitSignificance[moduleGenes,], absModuleMembershipSelec)
  selectedDatExpr <- datExpr[,moduleGenes]
  signifList <- colnames(selectedDatExpr) ## this is the list of genes in the module
  geneTraitModule <- cbind(signifList,geneTraitModule)
  colnames(geneTraitModule) =c("gene","traitSignificance","moduleMembership")
  fileName = paste("geneTraitModule_",experiment,"_",names(traitInterest),"_inModule_",module,".txt", sep="");
  write.table(as.data.frame(geneTraitModule), file = fileName,
              row.names = FALSE, col.names = TRUE, quote = FALSE,sep="\t");
  
}

