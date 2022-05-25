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
### Compare CV2 between gene expression data and morphology
###########################################################
# -------------------------------------------------------------------------
#read gene expression file
ngs <- read.table(file="Results/normCounts_log2_Sepal_DEF_14085.txt",row.names = 1, h=T)
dim(ngs)
head(ngs)

#read shape params file
shape <- read.csv(file="Data/morphologyData_DEF_mgx.csv",h=T,row.names = 1)
head(shape)


InvLRC <- 1/shape$LR423 
InvTRC <- 1/shape$TR423 
shapeMod <- shape[,c(1,2,5,7)]
shapeMod <- cbind(shapeMod,InvLRC)
shapeMod <- cbind(shapeMod,InvTRC)
shapeMod <- cbind(shapeMod,shape$GaussianCurvature)
shapePerParam <- t(shapeMod)
rownames(shapePerParam)=c("Length","Width","AspRatio","Area",
                          "LongCurv","TransCurv","GaussCurv")

write.table(t(shapePerParam),file="Results/shapeParams_7_pipe.txt",row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")
#newParams <- read.table("Results/shapeParams_7_pipe.txt",h=T,row.names=1,sep="\t")

###calculate CV2 for every gene
ed<-2^ngs
means <- rowMeans(ed)
vars <- apply(ed,1,var)
cv2 <- vars/means^2

#calculate Cv2 for every shape param
meansShape <- rowMeans(shapePerParam)
varsShape <- apply(shapePerParam,1,var)
cv2Shape <- varsShape/meansShape^2

colorBoxes <- brewer.pal(n = 8, name = "Dark2")
pdf("Plots/Figure3A_CV2_GeneExpAndShape_DensityAndBoxplot.pdf",width=5, height=3.5)

par(mar=c(3.5,3,.5,.5),mgp=c(2,0.65,0),cex=0.9, mfrow=c(1,1))

plot(density(log(cv2)),ylim=c(-.25,.4),xlim=c(-8,5),col="grey60",lwd=4,main="",
     xlab="",yaxt="n",ylab="",bty="n",xaxt="n")
mtext('Density', side=2, line=2, at=0.2,cex=1)
mtext(expression(log(CV^2)),side=1,line=2,at=-3)
lines(density(log(cv2Shape)),col=colorBoxes[6],lwd=4)
axis(2,at=c(0,.2,.4),las=2,cex=.8)
axis(1,at=c(-8,-6,-4,-2,0,2),cex=.8)

boxplot(log(cv2),log(cv2Shape),
        at=c(-.1,-.2),
        horizontal=TRUE,
        col=c("grey60",colorBoxes[6]),
        boxwex=c(.08,.08),add=TRUE,alpha=.5,xaxt="n",yaxt="n",frame.plot=F)
points(log(unname(cv2Shape[1])),-.2,pch=20,cex=2,col=colorBoxes[1])
points(log(unname(cv2Shape[2])),-.175,pch=20,cex=2,col=colorBoxes[2])
points(log(unname(cv2Shape[3])),-.22,pch=20,cex=2,col=colorBoxes[3])
points(log(unname(cv2Shape[4])),-.2,pch=20,cex=2,col=colorBoxes[4])
points(log(unname(cv2Shape[5])),-.19,pch=20,cex=2,col=colorBoxes[5])
points(log(unname(cv2Shape[6])),-.215,pch=20,cex=2,col=colorBoxes[7])
points(log(unname(cv2Shape[7])),-.2,pch=20,cex=2,col=colorBoxes[8])


legend("topright",legend=c("Gene expression","Morphology","Length","Width",
                           "Aspect ratio","Area","Long. curvature","Trans. curvature","Gauss. curvature"),
       col=c("grey60",colorBoxes[6],colorBoxes[1],colorBoxes[2],colorBoxes[3],colorBoxes[4],
             colorBoxes[5],colorBoxes[7],colorBoxes[8]),
       lty=c(1,1,NA,NA,NA,NA,NA,NA,NA),
       lwd=c(4,4,NA,NA,NA,NA,NA,NA,NA),
       pch=c(NA,NA,16,16,16,16,16,16,16),
       cex=.8,
       pt.cex=1.5,
       bty="n"
)

dev.off()

