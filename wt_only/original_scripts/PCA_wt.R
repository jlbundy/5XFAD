# R script for analyzing mass spec protomics and/or RNA-seq data
# Joseph L Bundy
# project started 10_9_2013

##################################################################################1
#Block 1 is basic I/O and data structuring 
#generates objects used in later analyses

.libPaths(new = "C:/Program Files/RRO/R-3.2.2/library")

#import libraries
library("stats")
library("rgl")
library("MASS")
#library("DESeq")
library("edgeR")
library("vegan")
library("gplots")
library("DESeq2")
library("tcltk2")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
#library(RUVSeq)
library(ggplot2)
library(matrixStats)
library("gridExtra")



JLB_PCA = function(object, intgroup="condition", ntop=500, returnData=FALSE)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=" : "))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2],PC3=pca$x[,3],PC4=pca$x[,4],PC5=pca$x[,5], group=group, intgroup.df, name=colnames(object))
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:5]
    return(d)
  }
  
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    coord_fixed()
}



#set working directory to input file location
setwd("C:/Users/joseph.bundy/OneDrive/cloud/R/5XFAD/ngs")
#setwd("C:/Users/josep/OneDrive/Cloud/R/5XFAD/ngs")

in_data_RNA <- read.csv("input/5XFAD_counts_unnorm.csv", header = T, stringsAsFactors=FALSE)
raw_data_RNA <- in_data_RNA[-c(1:3)]
raw_data_RNA <- data.frame(raw_data_RNA)
row.names(raw_data_RNA) <- in_data_RNA$Gene.ID

meta= read.csv("input/5XFAD_NGS_meta.csv",header = F)
meta = data.frame(meta)
meta = meta[,2:61]


#parse out information for specific pairwise comparison
raw_data_RNA = cbind(raw_data_RNA[1:10],raw_data_RNA[21:30],raw_data_RNA[41:50])
meta = cbind(meta[1:10],meta[21:30],meta[41:50])

setwd("C:/Users/joseph.bundy/OneDrive/cloud/R/5XFAD/ngs/WT_only/PCA")
#setwd("C:/Users/josep/OneDrive/Cloud/R/5XFAD/ngs/WT_only/PCA")

# filter <- apply(raw_data_RNA, 1, function(x) length(x[x>5])>=5)
# raw_data_RNA <- raw_data_RNA[filter,]
# genes <- rownames(raw_data_RNA)[grep("ENS", rownames(raw_data_RNA))]

#meta 1 for sex by genotype by timepoint group
#meta 2 for sex
#meta 3 for genotype
#meta 4 for time point

group=meta[2,]
group = as.matrix(group)
group = as.character(group)

 
sex=meta[3,]
sex = as.matrix(sex)
sex = as.character(sex)

timepoint=meta[5,]
timepoint = as.matrix(timepoint)
timepoint = as.character(timepoint)

batch=meta[30,]
batch = as.matrix(batch)
batch = as.character(batch)

group_3=meta[7,]
group_3 = as.matrix(group_3)
group_3 = as.character(group_3)

stage=meta[9,]
stage = as.matrix(stage)
stage = as.character(stage)

# set <- newSeqExpressionSet(as.matrix(raw_data_RNA),
#                            phenoData = data.frame(group, row.names=colnames(raw_data_RNA)))
# 
# 
# differences <- seq(1:60)
# differences <- matrix(differences, nrow = 12, byrow = TRUE)
# set_RUVs <- RUVs(set, genes, k=1, differences)
# pData(set_RUVs)
# 
# SV <- as.numeric(pData(set_RUVs)[,2])

rownames(raw_data_RNA) = NULL

#DESeq2 analysis
SampleName=colnames(raw_data_RNA)
SampleTable = data.frame(SampleName = SampleName, sex = sex,  timepoint = timepoint, group = group, stage = stage)

CountData <- as.matrix(raw_data_RNA)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = CountData,
  colData = SampleTable,
  design = ~group)

dds = ddsFullCountTable
dds <- DESeq(dds)

# detach("package:RUVSeq", unload=TRUE)
# detach("package:EDASeq", unload=TRUE)
# detach("package:DESeq", unload=TRUE)


rld <- rlogTransformation(dds, blind=FALSE)
sampleDists <- dist( t( assay(rld) ) )
sampleDists

png(filename="figures/sample_correlations_rld.png", res = 500, unit = "in", width = 10, height = 10)
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$genotype, rld$timepoint,rld$sex, sep="-" )
colnames(sampleDistMatrix) <- paste( rld$genotype, rld$timepoint,rld$sex, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

png(filename="figures/sample_correlations_poisson.png", res = 500, unit = "in", width = 10, height = 10)
poisd <- PoissonDistance(t(counts(dds, normalized=TRUE)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$genotype, rld$timepoint,rld$sex, sep="-" )
colnames(samplePoisDistMatrix) <- paste( rld$genotype, rld$timepoint,rld$sex, sep="-" )
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors)
dev.off()


vst <- varianceStabilizingTransformation(dds, blind = FALSE,fitType = "parametric")
sampleDists <- dist( t( assay(vst) ) )
sampleDists
png(filename="figures/sample_correlations_vst.png", res = 500, unit = "in", width = 10, height = 10)
sampleDistMatrix_vst <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$genotype, rld$timepoint,rld$sex, sep="-" )
colnames(sampleDistMatrix) <- paste( rld$genotype, rld$timepoint,rld$sex, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

#print(plotPCA(rld, intgroup=c("sex","genotype")))
PCA_object <- plotPCA(rld, intgroup = c( "sex","timepoint"), returnData=TRUE, n = 500)
percentVar <- round(100 * attr(PCA_object, "percentVar"))

# PCA_object_1M = subset(PCA_object, PCA_object$timepoint=="1m")
# PCA_object_2M = subset(PCA_object, PCA_object$timepoint=="2m")
# PCA_object_4M = subset(PCA_object, PCA_object$timepoint=="4m")
# PCA_object_M = subset(PCA_object, PCA_object$sex=="m")
# PCA_object_F = subset(PCA_object, PCA_object$sex=="f")
# PCA_object_AD = subset(PCA_object, PCA_object$genotype=="ad")
# PCA_object_WT = subset(PCA_object, PCA_object$genotype=="wt")


pch = c(rep(16,10),rep(17,10),rep(15,10))

scatterplot3d_pch = c(rep(1,10),rep(2,10),rep(0,10))

bg = c(rep(c(rep("#5a5aff",5),rep("#FFB6C1",5)),3))

col = c(rep(c(rep("#5a5aff",5),rep("#FFB6C1",5)),3))


png(filename="figures/PCA_rld_500.png", res = 500, unit = "in", width = 5, height = 5)
ggplot(PCA_object, aes(PC1, PC2, shape=timepoint,fill = factor(col))) + geom_point( size=4, aes(shape = timepoint)) +
  scale_colour_manual(values=c("#FFB6C1", "#5a5aff", "black"))+
  scale_fill_manual(values=c("#FFB6C1", "#5a5aff"))+
  scale_shape_manual(values = c(21,24,22))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()




PCA_object <- JLB_PCA(rld, intgroup = c( "sex","timepoint"), returnData=TRUE, n = 500)
percentVar <- round(100 * attr(PCA_object, "percentVar"))

#write.csv(PCA_object,file = "input/PCA_data.csv")


empty1<- ggplot()+geom_point(aes(1,1), colour="white") +
  geom_text(aes(x = 1, y = 1, label = "PC1"), color="black", size=10, parse = TRUE)+
  theme(          
    plot.background = element_blank(), 
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


empty2<- ggplot()+geom_point(aes(1,1), colour="white") +
  geom_text(aes(x = 1, y = 1, label = "PC2"), color="black", size=10, parse = TRUE)+
  theme(                              
    plot.background = element_blank(), 
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
empty3<- ggplot()+geom_point(aes(1,1), colour="white") +
  geom_text(aes(x = 1, y = 1, label = "PC3"), color="black", size=10, parse = TRUE)+
  theme(                              
    plot.background = element_blank(), 
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

plot1 = ggplot(PCA_object, aes(PC1, PC2, shape=timepoint,fill = factor(col))) + geom_point( size=4, aes(shape = timepoint)) +
  scale_colour_manual(values=c("#FFB6C1", "#5a5aff", "black"))+
  scale_fill_manual(values=c("#FFB6C1", "#5a5aff"), labels = c("females", "males"))+
  scale_shape_manual(values = c(21,24,22))+
  guides(fill=guide_legend(title = "sex", override.aes = list(colour = c("#FFB6C1", "#5a5aff"))))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #theme(legend.position = "none")
  theme(legend.position=c(.5, .5))


plot2 = ggplot(PCA_object, aes(PC1, PC3, shape=timepoint,fill = factor(col))) + geom_point( size=4, aes(shape = timepoint)) +
  scale_colour_manual(values=c("#FFB6C1", "#5a5aff", "black"))+
  scale_fill_manual(values=c("#FFB6C1", "#5a5aff"))+
  scale_shape_manual(values = c(21,24,22))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC3: ",percentVar[3],"% variance")) +
theme(legend.position = "none")

plot3 = ggplot(PCA_object, aes(PC2, PC1, shape=timepoint,fill = factor(col))) + geom_point( size=4, aes(shape = timepoint)) +
  scale_colour_manual(values=c("#FFB6C1", "#5a5aff", "black"))+
  scale_fill_manual(values=c("#FFB6C1", "#5a5aff"))+
  scale_shape_manual(values = c(21,24,22))+
  xlab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylab(paste0("PC1: ",percentVar[1],"% variance")) +
theme(legend.position = "none")

plot4 = ggplot(PCA_object, aes(PC2, PC3, shape=timepoint,fill = factor(col))) + geom_point( size=4, aes(shape = timepoint)) +
  scale_colour_manual(values=c("#FFB6C1", "#5a5aff", "black"))+
  scale_fill_manual(values=c("#FFB6C1", "#5a5aff"))+
  scale_shape_manual(values = c(21,24,22))+
  xlab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylab(paste0("PC3: ",percentVar[3],"% variance")) +
theme(legend.position = "none")

plot5 = ggplot(PCA_object, aes(PC3, PC1, shape=timepoint,fill = factor(col))) + geom_point( size=4, aes(shape = timepoint)) +
  scale_colour_manual(values=c("#FFB6C1", "#5a5aff", "black"))+
  scale_fill_manual(values=c("#FFB6C1", "#5a5aff"))+
  scale_shape_manual(values = c(21,24,22))+
  xlab(paste0("PC3: ",percentVar[3],"% variance")) +
  ylab(paste0("PC1: ",percentVar[1],"% variance")) +
theme(legend.position = "none")

plot6 = ggplot(PCA_object, aes(PC3, PC2, shape=timepoint,fill = factor(col))) + geom_point( size=4, aes(shape = timepoint)) +
  scale_colour_manual(values=c("#FFB6C1", "#5a5aff", "black"))+
  scale_fill_manual(values=c("#FFB6C1", "#5a5aff"))+
  scale_shape_manual(values = c(21,24,22))+
  xlab(paste0("PC3: ",percentVar[3],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
theme(legend.position = "none")

png(filename="figures/PCA_rld_3.png", res = 500, unit = "in", width = 15, height = 15)
figure_01 = grid.arrange(empty1, plot3, plot5, plot1, empty2, plot6, plot2, plot4, empty3,  ncol=3, nrow=3, widths=c(5,5,5), heights=c(5,5,5))
dev.off()


#these two make interactive 3d scatter plots but there's not way to change the symbols
library(rgl)

#plot3d(PCA_object[,1:3], col = col, pch = pch, size = 10)
 
library(car)

#scatter3d(PCA_object$PC1,PCA_object$PC2,PCA_object$PC3, surface = FALSE, shapes = cube3d(), groups = as.factor(group), elipsoid = TRUE)
 
library(scatterplot3d)

png(filename="figures/3D_PCA_1.png", res = 500, unit = "in", width = 5, height = 5)
s3d <- scatterplot3d(PCA_object$PC1,PCA_object$PC2,PCA_object$PC3, angle = -300, 
pch = pch, color = col, bg="grey",
xlab = paste0("PC1: ",percentVar[1],"% variance"), 
ylab = paste0("PC2: ",percentVar[2],"% variance"), 
zlab = paste0("PC3: ",percentVar[3],"% variance"),
tick.marks = TRUE)
legend(s3d$xyz.convert(0, -8, -8), pch = c(21,24,22), yjust=0,
       legend = c("1M", "2M", "4M"), cex = 1)
dev.off()


png(filename="figures/3D_PCA_2.png", res = 500, unit = "in", width = 5, height = 5)
s3d <- scatterplot3d(PCA_object$PC1,PCA_object$PC2,PCA_object$PC3, angle = -200, 
pch = pch, color = col, bg="grey",
xlab = "PC 1", ylab = "PC 2", zlab = "PC 3", tick.marks = TRUE)
legend(s3d$xyz.convert(-15, 2, -2), pch = c(21,24,22), yjust=0,
       legend = c("1M", "2M", "4M"), cex = 1)
dev.off()


points= c(rep("s",10),rep("c",10),rep("b",10))
cols = col = c(rep(c(rep("#FFB6C1",5),rep("#5a5aff",5)),3))
#FFB6C1 #4563e7
library(vrmlgen)
library(misc3d)
 
PCA_matrix <- as.matrix(PCA_object[1:3])

group = factor(group, levels = c("wt_1m_f", "wt_1m_m" ,
                                    "wt_2m_f", "wt_2m_m", 
                                    "wt_4m_f" ,"wt_4m_m"))

 # create VRML ouput file
 cloud3d(PCA_matrix, labels = timepoint,cols = cols,pointstyle = points,filename = "PCA_WT_only.vrml", 
  lab.axis = c("PC 1", "PC 2", "PC 3"), metalabels = group, showlegend = TRUE)
 
 
#  png(filename="figures/all/PCA_rld_3d.png", res = 500, unit = "in", width = 5, height = 5)
#  scatterplot3d(PCA_object[,1:3], color = col, pch = pch)
# dev.off()

# png(filename="figures/all/PCA_1M_rld_1000.png", res = 500, unit = "in", width = 5, height = 5)
# ggplot(PCA_object_1M, aes(PC1, PC2, colour=genotype, shape=sex)) + geom_point(size=4) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance"))
# dev.off()
# 
# png(filename="figures/all/PCA_2M_rld_1000.png", res = 500, unit = "in", width = 5, height = 5)
# ggplot(PCA_object_2M, aes(PC1, PC2, colour=genotype, shape=sex)) + geom_point(size=4) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance"))
# dev.off()
# 
# png(filename="figures/all/PCA_4M_rld_1000.png", res = 500, unit = "in", width = 5, height = 5)
# ggplot(PCA_object_4M, aes(PC1, PC2, colour=genotype, shape=sex)) + geom_point(size=4) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance"))
# dev.off()
# 
# png(filename="figures/all/PCA_AD_rld_1000.png", res = 500, unit = "in", width = 5, height = 5)
# ggplot(PCA_object_AD, aes(PC1, PC2, colour=timepoint, shape=sex)) + geom_point(size=4) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance"))
# dev.off()
# 
# png(filename="figures/all/PCA_WT_rld_1000.png", res = 500, unit = "in", width = 5, height = 5)
# ggplot(PCA_object_WT, aes(PC1, PC2, colour=timepoint, shape=sex)) + geom_point(size=4) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance"))
# dev.off()

 
 pch = c(rep(16,10),rep(17,10),rep(15,10))
 
 bg = c(rep(c(rep("#FF69B4",5),rep("#5a5aff",5)),3))
 
 col = c(rep(c(rep("#FF69B4",5),rep("#5a5aff",5)),3))
 
 
 
 PCA_object <- JLB_PCA(vst, intgroup = c( "sex","timepoint"), returnData=TRUE, n = 500)
 percentVar <- round(100 * attr(PCA_object, "percentVar"))
 
 #write.csv(PCA_object,file = "input/PCA_data.csv")
 
 
 empty1<- ggplot()+geom_point(aes(1,1), colour="white") +
   geom_text(aes(x = 1, y = 1, label = "PC1"), color="black", size=10, parse = TRUE)+
   theme(                              
     plot.background = element_blank(), 
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
 
 
 empty2<- ggplot()+geom_point(aes(1,1), colour="white") +
   geom_text(aes(x = 1, y = 1, label = "PC2"), color="black", size=10, parse = TRUE)+
   theme(                              
     plot.background = element_blank(), 
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
 empty3<- ggplot()+geom_point(aes(1,1), colour="white") +
   geom_text(aes(x = 1, y = 1, label = "PC3"), color="black", size=10, parse = TRUE)+
   theme(                              
     plot.background = element_blank(), 
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
 
 plot1 = ggplot(PCA_object, aes(PC1, PC2, shape=timepoint,fill = factor(col))) + geom_point( size=4, aes(shape = timepoint)) +
   scale_colour_manual(values=c("#FF69B4", "#5a5aff", "black"))+
   scale_fill_manual(values=c("#FF69B4", "#5a5aff"), labels = c("females", "males"))+
   scale_shape_manual(values = c(21,24,22))+
   guides(fill=guide_legend(title = "sex", override.aes = list(colour = c("#FF69B4", "#5a5aff"))))+
   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
   #theme(legend.position = "none")
   theme(legend.position=c(1.25, .5))
 
 
 plot2 = ggplot(PCA_object, aes(PC1, PC3, shape=timepoint,fill = factor(col))) + geom_point( size=4, aes(shape = timepoint)) +
   scale_colour_manual(values=c("#FF69B4", "#5a5aff", "black"))+
   scale_fill_manual(values=c("#FF69B4", "#5a5aff"))+
   scale_shape_manual(values = c(21,24,22))+
   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
   ylab(paste0("PC3: ",percentVar[3],"% variance")) +
   theme(legend.position = "none")
 
 plot3 = ggplot(PCA_object, aes(PC2, PC1, shape=timepoint,fill = factor(col))) + geom_point( size=4, aes(shape = timepoint)) +
   scale_colour_manual(values=c("#FF69B4", "#5a5aff", "black"))+
   scale_fill_manual(values=c("#FF69B4", "#5a5aff"))+
   scale_shape_manual(values = c(21,24,22))+
   xlab(paste0("PC2: ",percentVar[2],"% variance")) +
   ylab(paste0("PC1: ",percentVar[1],"% variance")) +
   theme(legend.position = "none")
 
 plot4 = ggplot(PCA_object, aes(PC2, PC3, shape=timepoint,fill = factor(col))) + geom_point( size=4, aes(shape = timepoint)) +
   scale_colour_manual(values=c("#FF69B4", "#5a5aff", "black"))+
   scale_fill_manual(values=c("#FF69B4", "#5a5aff"))+
   scale_shape_manual(values = c(21,24,22))+
   xlab(paste0("PC2: ",percentVar[2],"% variance")) +
   ylab(paste0("PC3: ",percentVar[3],"% variance")) +
   theme(legend.position = "none")
 
 plot5 = ggplot(PCA_object, aes(PC3, PC1, shape=timepoint,fill = factor(col))) + geom_point( size=4, aes(shape = timepoint)) +
   scale_colour_manual(values=c("#FF69B4", "#5a5aff", "black"))+
   scale_fill_manual(values=c("#FF69B4", "#5a5aff"))+
   scale_shape_manual(values = c(21,24,22))+
   xlab(paste0("PC3: ",percentVar[3],"% variance")) +
   ylab(paste0("PC1: ",percentVar[1],"% variance")) +
   theme(legend.position = "none")
 
 plot6 = ggplot(PCA_object, aes(PC3, PC2, shape=timepoint,fill = factor(col))) + geom_point( size=4, aes(shape = timepoint)) +
   scale_colour_manual(values=c("#FF69B4", "#5a5aff", "black"))+
   scale_fill_manual(values=c("#FF69B4", "#5a5aff"))+
   scale_shape_manual(values = c(21,24,22))+
   xlab(paste0("PC3: ",percentVar[3],"% variance")) +
   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
   theme(legend.position = "none")
 
tiff(filename="figures/PCA_vst_3.tiff", res = 400, unit = "in", width = 10, height = 10)
 figure_01 = grid.arrange(empty1, plot3, plot5, plot1, empty2, plot6, plot2, plot4, empty3,  ncol=3, nrow=3, widths=c(5,5,5), heights=c(5,5,5))
 dev.off()
 
 
 
 
#  png(filename="figures/PCA_rld_500_stage.png", res = 500, unit = "in", width = 5, height = 5)
#  ggplot(PCA_object, aes(PC1, PC2, shape=timepoint,fill = factor(stage))) + geom_point( size=4, aes(shape = timepoint)) +
#    scale_colour_manual(values=c("#FFB6C1", "#5a5aff", "black"))+
#    # scale_fill_manual(values=c("#FFB6C1", "#5a5aff"))+
#    scale_shape_manual(values = c(21,24,22))+
#    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#    ylab(paste0("PC2: ",percentVar[2],"% variance"))
#  dev.off()
#  
#  
 
 
 
 
#  plot1 = ggplot(PCA_object, aes(PC1, PC2, shape=timepoint,fill = factor(stage))) + geom_point( size=4, aes(shape = timepoint)) +
#    scale_colour_manual(values=c("#FFB6C1", "#5a5aff", "black"))+
#    scale_fill_manual(values=brewer.pal(n=5, name = "Set1"))+
#    scale_shape_manual(values = c(21,24,22))+
#    guides(fill=guide_legend(title = "stage", override.aes = list(colour = brewer.pal(n=5, name = "Set1"))))+
#    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#    theme(legend.position = "none")
 
 
 
 
 
 
 
 
 pch = c(rep(21,10),rep(24,10),rep(22,10))
 
 
 
 #these two make interactive 3d scatter plots but there's not way to change the symbols
 library(rgl)
 
 #plot3d(PCA_object[,1:3], col = col, pch = pch, size = 10)
 
 library(car)
 
 scatter3d(PCA_object$PC1,PCA_object$PC2,PCA_object$PC3, surface = FALSE, shapes = cube3d(), ellipsoid = TRUE, groups = as.factor(group), grid = FALSE)
 
 library(scatterplot3d)
 
 png(filename="figures/3D_PCA_vst_1.png", res = 400, unit = "in", width = 10, height = 10)
 s3d <- scatterplot3d(PCA_object$PC1,PCA_object$PC2,PCA_object$PC3, #angle = -200, 
                      pch = pch, color = "black", bg= col,
                      cex.symbols= 3, cex.axis = 2, cex.lab= 2,type="h",
                      xlab = paste0("PC1: ",percentVar[1],"% variance"), 
                      ylab = paste0("PC2: ",percentVar[2],"% variance"), 
                      zlab = paste0("PC3: ",percentVar[3],"% variance"),
  tick.marks = TRUE)
 legend(s3d$xyz.convert(0, -8, 0), pch = c(21,24,22), yjust=0,
        legend = c("1M", "2M", "4M"), cex = 2, pt.cex = 3)
 dev.off()
 
 
 png(filename="figures/3D_PCA_vst_2.png", res = 400, unit = "in", width = 5, height = 5)
 s3d <- scatterplot3d(PCA_object$PC1,PCA_object$PC2,PCA_object$PC3, angle = -200, 
                      pch = pch, color = col, bg="grey",type="h",
                      xlab = "PC 1", ylab = "PC 2", zlab = "PC 3", tick.marks = TRUE)
 legend(s3d$xyz.convert(-15, 2, -2), pch = c(21,24,22), yjust=0,
        legend = c("1M", "2M", "4M"), cex = 2, pt.cex = 3)
 dev.off()
 
 
 points= c(rep("s",20),rep("c",20),rep("b",20))
 cols = col = c(rep(c(rep("#FFB6C1",5),rep("#5a5aff",5),rep("#ff1e3f",5),rep("#4563e7",5)),3))
 #FFB6C1 #4563e7
 library(vrmlgen)
 library(misc3d)
 
 PCA_matrix <- as.matrix(PCA_object[1:3])
 
 group = factor(group, levels = c("wt_1m_f", "wt_1m_m" ,
                                  "wt_2m_f", "wt_2m_m", 
                                  "wt_4m_f" ,"wt_4m_m" ))
 
 # create VRML ouput file
 cloud3d(PCA_matrix, labels = timepoint,cols = cols,pointstyle = points,filename = "PCA_WT_only_vst.vrml", 
         lab.axis = c("PC 1", "PC 2", "PC 3"), metalabels = group, showlegend = TRUE)
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 #estrous cycle PCA with vst
 timepoint= c(timepoint[11:15],timepoint[21:25])
stage= c(stage[11:15],stage[21:25])
 
 #DESeq2 analysis
 SampleName=colnames(data.frame(raw_data_RNA[11:15],raw_data_RNA[21:25]) )
 SampleTable = data.frame(SampleName = SampleName,  timepoint = timepoint, stage = stage)
 
 CountData <- as.matrix(data.frame(raw_data_RNA[11:15],raw_data_RNA[21:25]))
 
 ddsFullCountTable <- DESeqDataSetFromMatrix(
   countData = CountData,
   colData = SampleTable,
   design = ~stage + timepoint)
 
 dds = ddsFullCountTable
 dds <- DESeq(dds)
 
 vst <- varianceStabilizingTransformation(dds, blind = FALSE,fitType = "parametric")
 sampleDists <- dist( t( assay(vst) ) )
 
 PCA_object <- JLB_PCA(vst, intgroup = c( "stage","timepoint"), returnData=TRUE, n = 500)
 percentVar <- round(100 * attr(PCA_object, "percentVar"))
 
 
 
 plot1 = ggplot(PCA_object, aes(PC1, PC2, shape=timepoint,fill = factor(stage))) + geom_point( size=4, aes(shape = timepoint)) +
   scale_colour_manual(values=c("#FFB6C1", "#5a5aff", "black"))+
   scale_fill_manual(values=brewer.pal(n=3, name = "Set1"), labels = c("Estrus", "Diestrus", "Metestrus"))+
   scale_shape_manual(values = c(24,22))+
   guides(fill=guide_legend(title = "stage", override.aes = list(colour = brewer.pal(n=3, name = "Set1"))))+
   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
   #theme(legend.position = "none")
   theme(legend.position=c(1.5, .5))
 
 
 plot2 = ggplot(PCA_object, aes(PC1, PC3, shape=timepoint,fill = factor(stage))) + geom_point( size=4, aes(shape = timepoint)) +
   scale_colour_manual(values=c("#FFB6C1", "#5a5aff", "black"))+
   scale_fill_manual(values=brewer.pal(n=3, name = "Set1"))+
   scale_shape_manual(values = c(24,22))+
   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
   ylab(paste0("PC3: ",percentVar[3],"% variance")) +
   theme(legend.position = "none")
 
 plot3 = ggplot(PCA_object, aes(PC2, PC1, shape=timepoint,fill = factor(stage))) + geom_point( size=4, aes(shape = timepoint)) +
   scale_colour_manual(values=c("#FFB6C1", "#5a5aff", "black"))+
   scale_fill_manual(values=brewer.pal(n=3, name = "Set1"))+
   scale_shape_manual(values = c(24,22))+
   xlab(paste0("PC2: ",percentVar[2],"% variance")) +
   ylab(paste0("PC1: ",percentVar[1],"% variance")) +
   theme(legend.position = "none")
 
 plot4 = ggplot(PCA_object, aes(PC2, PC3, shape=timepoint,fill = factor(stage))) + geom_point( size=4, aes(shape = timepoint)) +
   scale_colour_manual(values=c("#FFB6C1", "#5a5aff", "black"))+
   scale_fill_manual(values=brewer.pal(n=3, name = "Set1"))+
   scale_shape_manual(values = c(24,22))+
   xlab(paste0("PC2: ",percentVar[2],"% variance")) +
   ylab(paste0("PC3: ",percentVar[3],"% variance")) +
   theme(legend.position = "none")
 
 plot5 = ggplot(PCA_object, aes(PC3, PC1, shape=timepoint,fill = factor(stage))) + geom_point( size=4, aes(shape = timepoint)) +
   scale_colour_manual(values=c("#FFB6C1", "#5a5aff", "black"))+
   scale_fill_manual(values=brewer.pal(n=3, name = "Set1"))+
   scale_shape_manual(values = c(24,22))+
   xlab(paste0("PC3: ",percentVar[3],"% variance")) +
   ylab(paste0("PC1: ",percentVar[1],"% variance")) +
   theme(legend.position = "none")
 
 plot6 = ggplot(PCA_object, aes(PC3, PC2, shape=timepoint,fill = factor(stage))) + geom_point( size=4, aes(shape = timepoint)) +
   scale_colour_manual(values=c("#FFB6C1", "#5a5aff", "black"))+
   scale_fill_manual(values=brewer.pal(n=3, name = "Set1"))+
   scale_shape_manual(values = c(24,22))+
   xlab(paste0("PC3: ",percentVar[3],"% variance")) +
   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
   theme(legend.position = "none")
 
png(filename="figures/PCA_vst_3_stage.png", res = 300, unit = "in", width = 8, height = 8)
 figure_01 = grid.arrange(empty1, plot3, plot5, plot1, empty2, plot6, plot2, plot4, empty3,  ncol=3, nrow=3, widths=c(5,5,5), heights=c(5,5,5))
 dev.off()
 
