###############################################
### Sepal morphology and gene expression
### Author: Diego A. Hartasanchez
### Created: 25/08/2021
### Last edit: 20/04/2022
### - R version used: 4.0.3
### Pipeline:
### - Data processing: normalization, quality control, filtering
### 

## Set working directory
setwd('~/Documents/Academico/Proyectos/SepalGeneNetworks_2020/SepalMorphology_GitHubRepository/')

###########################################################
### Correlation plot (ggpairs) of morphology measurements #
###########################################################

shapeParams <- read.table("Results/shapeParams_7_pipe.txt",h=T,row.names=1,sep="\t")

#Plant Data specifies which samples correspond to which plant (D, E and F)
plantData <- read.table('Data/sampleDataType_Sepals_DEF.txt',header=TRUE,row.names=1)
dim(plantData)

plantColors <- c('red2','forestgreen','cornflowerblue','orange')

#Add plant type to data frame of morphology measurements
shapeDataNorm <- as.data.frame(shapeParams)
str(shapeDataNorm)
shapeDataNormPlant <- as.data.frame(cbind(shapeDataNorm, plantData$Plant))
colnames(shapeDataNormPlant)[8] <- "Plant"
shapeDataNormPlant$Plant <- as.factor(shapeDataNormPlant$Plant)
str(shapeDataNormPlant)

#Define function to fit line to datapoints
my_fitModel_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    #  geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(formula = y ~ x, method=lm, fill=plantColors[4], color="grey40", ...)
  p
}

#pdf(file = paste("Plots/corrBetweenPairs_morphology_noNorm.pdf",sep=""), wi = 6, he = 6)
#require(ggplot2)
#require(GGally)
ggplot2::theme_set(ggplot2::theme_bw())  ## (this just globally sets the ggplot2 theme to theme_bw) 
g = ggpairs(shapeDataNormPlant,columns = 1:8,
            lower = list(continuous = my_fitModel_fn),
            upper = list(continuous = wrap("cor", size = 2.5)),
            diag = list(continuous = "densityDiag"),
            aes(color = Plant,alpha = 0.9))+
  scale_colour_manual(name="Plant", values=plantColors)+
  scale_fill_manual(name="Plant", values=plantColors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7))+
  theme(axis.text.y = element_text(size=7))
#g
#dev.off()

g_density = ggpairs(shapeDataNormPlant[,-8],columns = 1:7,
                    lower = list(continuous = "blank"),
                    upper = list(continuous = "blank"),
                    diag = list(continuous = "densityDiag"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7))+
  theme(axis.text.y = element_text(size=7))
#g_density

densityPlots <- list()
densityPlots <- lapply(1:7, function(i) getPlot(g_density, i = i, j = i)) 
length(densityPlots)

plots = list()
length(plots)
plots <- c(plots, lapply(1:(g$nrow-1), function(i) getPlot(g, i = i, j = 8)+coord_flip())) #+ coord_flip() 
for (i in 1:7){
  plots <- c(plots, lapply(1:7, function(j) getPlot(g, i = i, j = j)))
}
plots[8]<-densityPlots[1]
plots[16]<-densityPlots[2]
plots[24]<-densityPlots[3]
plots[32]<-densityPlots[4]
plots[40]<-densityPlots[5]
plots[48]<-densityPlots[6]
plots[56]<-densityPlots[7]

#pdf(file = paste("Plots/Figure2A_corrBetweenPairs_morphology_noNorm_boxplot.pdf",sep=""), wi = 6.5, he = 6.5)
morph_pairPlot <- ggmatrix(plots, 
                           nrow=8,
                           ncol=7,
                           xAxisLabels = g$xAxisLabels[1:7], 
                           yAxisLabels = c(g$yAxisLabels[8],g$yAxisLabels[1:7]))
#morph_pairPlot
#dev.off()


#############################################
### PCA for morphological measurements ######
#############################################
# -------------------------------------------------------------------------


project.pca.shape <- prcomp(shapeDataNormPlant[, -8],  scale = TRUE)

pca1 <- fviz_pca_ind(project.pca.shape, habillage=shapeDataNormPlant$Plant,
                     addEllipses=TRUE, ellipse.level=0.90,repel=TRUE,
                     xlab="PC1 (56.5%)", ylab="PC2 (22.3%)",
                     title=NULL,legend="none"
) + 
  theme(legend.position = "none",
        plot.title = element_blank())
#pca1

pca2 <- fviz_pca_biplot(project.pca.shape, repel=TRUE,
                        col.var = "dodgerblue4",#2E9FDF", # Variables color
                        col.ind = "grey70", ##696969",  # Individuals color
                        xlab="PC1 (56.5%)", ylab="PC2 (22.3%)",
                        title=NULL,legend="none"
) + 
  theme(legend.position = "none",
        plot.title = element_blank())
#pca2


####################################################
### Generate Pdf for Figure2_Morphology_PCA_pairs ##
####################################################
# -------------------------------------------------------------------------


pdf(file = paste("Plots/Figure2_morphology_PCA_pairs.pdf",sep=""), wi = 10.5, he = 7)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 3)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(morph_pairPlot, vp = define_region(row = 1:2, col = 1:2))
print(pca1, vp = define_region(row = 1, col = 3))   # Span over two columns
print(pca2, vp = define_region(row = 2, col = 3))
dev.off()


