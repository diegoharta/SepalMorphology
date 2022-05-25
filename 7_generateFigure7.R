###############################################
### Sepal morphology and gene expression
### Author: Diego A. Hartasanchez
### Created: 09/12/2021
### Last edit: 20/04/2022
### - R version used: 4.0.3
### Pipeline:
### - Mutant analysis 
### - Figure 7 and supp figure 5

library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(ggpubr)
library(RColorBrewer)
colorBoxes <- brewer.pal(n = 8, name = "Dark2")
library(ggthemes) # Load
library(yarrr)
color16Pal <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", 
                "palegreen2", "#FDBF6F", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4", "brown") 


setwd('~/Documents/Academico/Proyectos/SepalGeneNetworks_2020/SepalMorphology_GitHubRepository/')

mutantData <- read.table(file="Data/mutantData_DiffCV2.txt",h=T)
mathildeCol0Data <- read.table(file="Data/mathildeCol0Phenotypes_Clean.tab",h=T)

mutantData$Batch
col0DataFiltered <- col0Data[col0Data$batch %in% mutantData$Batch,]
col0_cv2_area <- (col0DataFiltered$Area.sd/col0DataFiltered$Area.mean)^2 
mean(col0_cv2_area)
sd(col0_cv2_area)

fig7a <- ggplot(data = mutantData, aes(Gene_CV2, abs(AreaDiff))) +
  geom_smooth(method=lm,  color="grey",fill="khaki")+
  labs(x = "", y = "Knockout mutant change in Area (abs)")  + 
  scale_x_continuous(trans='log10')+
  scale_color_manual(values=color16Pal)+
  geom_point(aes(color = Gene_Name)) +
  theme_bw()+
  stat_cor(method = "pearson",label.x.npc = "left",label.y.npc = "top")+
  theme(axis.title=element_text(size=13),
        plot.margin=unit(c(.2,0,0,.5), 'cm'))+
  theme(legend.title = element_blank())

lengendFig7 <- get_legend(fig7a)
pdf("Plots/legend_Fig7.pdf",width=8, height=6.5)
plot(lengendFig7)
dev.off()


fig7a <- ggplot(data = mutantData, aes(Gene_CV2, abs(AreaDiff))) +
  #  geom_smooth(method=lm,  color="grey")+#,fill="khaki")+
  labs(x = "", y = "Area change (abs) in mutant")  + 
  scale_x_continuous(trans='log10')+
  #  scale_y_continuous(limits = c(min(testSetMutsClean$AreaDiff)-min(testSetMutsClean$AreaDiff)*.15,max(testSetMutsClean$AreaDiff)+max(testSetMutsClean$AreaDiff)*.15))+
  coord_cartesian(ylim = c(0.005, .145)) +
  scale_color_manual(values=color16Pal)+
  geom_point(aes(color = Gene_Name)) +
  stat_cor(method = "pearson",label.x.npc = "left",label.y.npc = "top")+
  theme(axis.title=element_text(size=13),
        plot.margin=unit(c(.7,0,.05,.6), 'cm'),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x =element_blank(),
        axis.ticks.x=element_blank())

fig7b  <- ggplot(data = mutantData, aes(Gene_meanExp, abs(AreaDiff))) +
  #geom_smooth(method=lm,  color="grey")+#,fill="khaki")+
  labs(x = "", y = "")+#ifference between knockout mutant and WT")  + 
  scale_x_continuous(trans='log10')+
  coord_cartesian(ylim = c(0.005, .145)) +
  #  scale_y_continuous(limits = c(min(testSetMutsClean$AreaDiff)-min(testSetMutsClean$AreaDiff)*.15,max(testSetMutsClean$AreaDiff)+max(testSetMutsClean$AreaDiff)*.15))+
  stat_cor(method = "pearson",label.x.npc = "left",label.y.npc = "top")+
  theme(axis.title=element_text(size=13),
        plot.margin=unit(c(.7,.7,.05,0), 'cm'),
        axis.title.y=element_blank(),
        axis.text.y =element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x =element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_rect(fill = "white", colour =  NA),
        legend.position='none',
        axis.line = element_line(colour = "black"))+
  scale_color_manual(values=color16Pal)+
  geom_point(aes(color = Gene_Name))

fig7c  <- ggplot(data = mutantData, aes(Gene_CV2, Area_CV2)) +
  #geom_smooth(method=lm,  color="grey")+#,fill="khaki")+
  labs(x = expression(paste("Gene ",CV^2," in WT",sep="")), 
       y = expression(paste("Area ",CV^2," in mutant"))) + 
  scale_x_continuous(trans='log10')+
  coord_cartesian(ylim = c(0, .05)) +
  # scale_y_continuous(limits = c(min(testSetMutsClean$Area_CV2)-(max(testSetMutsClean$Area_CV2)-min(testSetMutsClean$Area_CV2))*.1,max(testSetMutsClean$Area_CV2)+(max(testSetMutsClean$Area_CV2)-min(testSetMutsClean$Area_CV2))*.1))+
  stat_cor(method = "pearson",label.x.npc = "left",label.y.npc = "top")+
  theme(axis.title=element_text(size=13),
        plot.margin=unit(c(0,0,.5,.5), 'cm'),
        panel.background = element_rect(fill = "white", colour =  NA),
        legend.position='none',
        axis.line = element_line(colour = "black"))+
  #geom_point(color= colorsInfo[1]) +
  scale_color_manual(values=color16Pal)+
  geom_point(aes(color = Gene_Name))+
  geom_hline(yintercept=mean(col0_cv2_area), linetype="dashed", color = "darkgrey")+
  geom_hline(yintercept=mean(col0_cv2_area)+sd(col0_cv2_area), linetype="dotted", color = "darkgrey")+
  geom_hline(yintercept=mean(col0_cv2_area)-sd(col0_cv2_area), linetype="dotted", color = "darkgrey")

fig7d  <- ggplot(data = mutantData, aes(Gene_meanExp, Area_CV2)) +
  geom_smooth(method=lm,  color="grey",fill="khaki")+
  labs(x = "Gene mean expression in WT", y = "") + 
  scale_x_continuous(trans='log10')+
  #scale_y_continuous(limits = c(min(testSetMutsClean$Area_CV2)-(max(testSetMutsClean$Area_CV2)-min(testSetMutsClean$Area_CV2))*.2,max(testSetMutsClean$Area_CV2)+(max(testSetMutsClean$Area_CV2)-min(testSetMutsClean$Area_CV2))*.2))+
  coord_cartesian(ylim = c(0, .05)) +
  stat_cor(method = "pearson",label.x.npc = "left",label.y.npc = "top")+
  theme(axis.title=element_text(size=13),
        plot.margin=unit(c(0,.7,.63,0), 'cm'),
        axis.title.y=element_blank(),
        axis.text.y =element_blank(),
        axis.ticks.y=element_blank(),
        legend.position='none',
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = "white", colour = NA))+
  #geom_point(color= colorsInfo[1]) +
  scale_color_manual(values=color16Pal)+
  geom_point(aes(color = Gene_Name)) +
  geom_hline(yintercept=mean(col0_cv2_area), linetype="dashed", color = "darkgrey")+
  geom_hline(yintercept=mean(col0_cv2_area)+sd(col0_cv2_area), linetype="dotted", color = "darkgrey")+
  geom_hline(yintercept=mean(col0_cv2_area)-sd(col0_cv2_area), linetype="dotted", color = "darkgrey")

pdf("Plots/Figure7_mutant_Area_4plots.pdf",width=8, height=6)
grid.arrange(fig7a,fig7b,fig7c,fig7d, 
             ncol=2, nrow=2, widths=c(4, 3.5), heights=c(3, 3.5))
dev.off()

#####
##supplementary figure 5

figSupp5a <- ggplot(data = mutantData, aes(Gene_meanExp, Length_CV2)) +
  geom_smooth(method=lm,  color="grey",fill="khaki")+
  labs(x = "Gene mean expression in WT",  
       y = expression(paste("Length ",CV^2," in gene knockout mutant"))) + 
  scale_x_continuous(trans='log10')+
  scale_color_manual(values=color16Pal)+
  geom_point(aes(color = Gene_Name)) +
  stat_cor(method = "pearson",label.x.npc = "left",label.y.npc = "top")+
  ggtitle("Length (16 mutants)") +
  theme(axis.title=element_text(size=13),
        plot.margin=unit(c(.5,.5,.5,.5), 'cm'),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
  )
#figSupp5a

figSupp5b <- ggplot(data = mutantData, aes(Gene_meanExp, Width_CV2)) +
  geom_smooth(method=lm,  color="grey",fill="khaki")+
  labs(x = "Gene mean expression in WT",  
       y = expression(paste("Width ",CV^2," in gene knockout mutant"))) + 
  scale_x_continuous(trans='log10')+
  scale_color_manual(values=color16Pal)+
  geom_point(aes(color = Gene_Name)) +
  stat_cor(method = "pearson",label.x.npc = "left",label.y.npc = "top")+
  ggtitle("Width (16 mutants)") +
  theme(axis.title=element_text(size=13),
        plot.margin=unit(c(.5,.5,.5,.5), 'cm'),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
  )
#figSupp5b

figSupp5c <- ggplot(data = mutantData, aes(Gene_meanExp, AspectRatio_CV2)) +
  geom_smooth(method=lm,  color="grey",fill="khaki")+
  labs(x = "Gene mean expression in WT",  
       y = expression(paste("AspRatio ",CV^2," in gene knockout mutant"))) + 
  scale_x_continuous(trans='log10')+
  scale_color_manual(values=color16Pal)+
  geom_point(aes(color = Gene_Name)) +
  stat_cor(method = "pearson",label.x.npc = "left",label.y.npc = "top")+
  ggtitle("AspRatio (16 mutants)") +
  theme(axis.title=element_text(size=13),
        plot.margin=unit(c(.5,.5,.5,.5), 'cm'),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
  )
#figSupp5c

pdf("Plots/Supp_Figure5_mutant_LengthWidthAspRatio_3plots.pdf",width=12, height=4)
grid.arrange(figSupp5a,figSupp5b,figSupp5c, 
             ncol=3, nrow=1, widths=c(4,4,4), heights=c(4))
dev.off()



