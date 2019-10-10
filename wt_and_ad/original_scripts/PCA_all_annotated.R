# R script for analyzing mass spec protomics and/or RNA-seq data
# Joseph L Bundy
# project started 10_9_2013

##################################################################################1
#Block 1 is basic I/O and data structuring 
#generates objects used in later analyses

#import libraries

library("DESeq2")
library("ggplot2")
library("gridExtra")
library("scatterplot3d")


#I stole and modified this function from the one in the DESeq2 package - I had to modify it so it would spit out 
#3 PCs instead of only 2

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
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2],PC3=pca$x[,3], group=group, intgroup.df, name=colnames(object))
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:3]
    return(d)
  }
  
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    coord_fixed()
}



#set working directory
setwd("D:/OneDrive/Cloud/R/5XFAD/ngs")

#read in your data - I have all my input files in a subdirectory called "input", so it'll look there
in_data_RNA <- read.csv("input/5XFAD_counts_unnorm.csv", header = T, stringsAsFactors=FALSE)

#this just splits off the first 3 columns metadata so you have a nice clean dataframe - adjust as needed
raw_data_RNA <- in_data_RNA[-c(1:3)]

#converty into a dataframe and give these Gene.ID rownames
raw_data_RNA <- data.frame(raw_data_RNA)
row.names(raw_data_RNA) <- in_data_RNA$Gene.ID

#I stored all the metadata in this file - you can make a similar one if you like and import it the same way
meta= read.csv("input/5XFAD_NGS_meta.csv",header = F)
meta = data.frame(meta)
meta = meta[,2:61]

#change working directory again if you want to do this work somewhere that isn't your main directory with the file inputs
setwd("D:/OneDrive/Cloud/R/5XFAD/ngs/AD_and_WT/PCA")

#I split the groups off from the meta sheet into their own objects - this is probably unecessary but it made the code more clear to me
group=meta[2,]
group = as.matrix(group)
group = as.character(group)
group = group[1:60]
 
sex=meta[3,]
sex = as.matrix(sex)
sex = as.character(sex)

genotype=meta[4,]
genotype = as.matrix(genotype)
genotype = as.character(genotype)

timepoint=meta[5,]
timepoint = as.matrix(timepoint)
timepoint = as.character(timepoint)

batch=meta[30,]
batch = as.matrix(batch)
batch = as.character(batch)

group_3=meta[7,]
group_3 = as.matrix(group_3)
group_3 = as.character(group_3)

rownames(raw_data_RNA) = NULL

#I went ahead and followed through a DESeq2 analysis - this allows you to get normalized 
#read counts, which can then be transformed with other Deseq2 options

#DESeq2 analysis
SampleName=colnames(raw_data_RNA)
SampleTable = data.frame(SampleName = SampleName, sex = sex, genotype = genotype, timepoint = timepoint, group = group)

CountData <- as.matrix(raw_data_RNA)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = CountData,
  colData = SampleTable,
  design = ~group)

dds = ddsFullCountTable
dds <- DESeq(dds)


#there are two transformations built into DESeq2 - we'll make objects of transformed data for both
#the rld and vst transformations - can try PCA on data from one/both of those

#these take a sec to run so be patient

rld <- rlogTransformation(dds, blind=FALSE)
sampleDists <- dist( t( assay(rld) ) )
sampleDists

vst <- varianceStabilizingTransformation(dds, blind = FALSE,fitType = "parametric")
sampleDists <- dist( t( assay(vst) ) )
sampleDists



#here's the code of making the basic 2 dimensional PCA 
PCA_object <- plotPCA(rld, intgroup = c( "genotype", "sex","timepoint"), returnData=TRUE, n = 500)
percentVar <- round(100 * attr(PCA_object, "percentVar"))


#make some objects for plotting the different conditions

#first make pch that are different for ages and for WT and AD (AD have outlines)
pch = c(rep(21,20),rep(24,20),rep(22,20))

bg = c(rep(c(rep("#5a5aff",5),rep("#FF69B4",5)),3))

#make pink and blue outlines for WT to make it look like they don't have one
#and then black outlines for AD samples
outline = c(rep(c(rep("#5a5aff",5),rep("#FF69B4",5),rep("black",10)),3))

#fill (ping for females, blue for males)
col = c(rep(c(rep("#5a5aff",5),rep("#FF69B4",5)),6))

size = c(rep(1,20),rep(2,20),rep(4,20))

#this will go ahead and print out a "simple" version of the pca
png(filename="figures/PCA_rld_500.png", res = 500, unit = "in", width = 5, height = 5)
ggplot(PCA_object, aes(PC1, PC2, shape=timepoint,fill = factor(col), colour = factor(outline))) + geom_point( size=4, aes(shape = timepoint)) +
  scale_colour_manual(values=c("#FF69B4", "#5a5aff", "black"))+
  scale_fill_manual(values=c("#FF69B4", "#5a5aff"))+
  scale_shape_manual(values = c(21,24,22))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()



#make a PCA object using the custom function that will return the top 3 PCAs 
PCA_object <- JLB_PCA(rld, intgroup = c( "genotype", "sex","timepoint"), returnData=TRUE, n = 500)
percentVar <- round(100 * attr(PCA_object, "percentVar"))


#now we'll go ahead and generate a bunch of ggplot objects and names them - we'll start with some "empty" 
#ones that will fill the spaces in the final figure we want to build

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

#These are all basically the same kind of plot but with the axes changed up with different PCs
#represented on different axes

plot1 = ggplot(PCA_object, aes(PC1, PC2, shape=timepoint,fill = factor(col), colour = factor(outline))) + geom_point( size=4, aes(shape = timepoint)) +
  scale_colour_manual(values=c("#FF69B4", "#5a5aff", "black"))+
  scale_fill_manual(values=c("#FF69B4", "#5a5aff"))+
  scale_shape_manual(values = c(21,24,22))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
theme(legend.position = "none")


plot2 = ggplot(PCA_object, aes(PC1, PC3, shape=timepoint,fill = factor(col), colour = factor(outline))) + geom_point( size=4, aes(shape = timepoint)) +
  scale_colour_manual(values=c("#FF69B4", "#5a5aff", "black"))+
  scale_fill_manual(values=c("#FF69B4", "#5a5aff"))+
  scale_shape_manual(values = c(21,24,22))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC3: ",percentVar[3],"% variance")) +
theme(legend.position = "none")

plot3 = ggplot(PCA_object, aes(PC2, PC1, shape=timepoint,fill = factor(col), colour = factor(outline))) + geom_point( size=4, aes(shape = timepoint)) +
  scale_colour_manual(values=c("#FF69B4", "#5a5aff", "black"))+
  scale_fill_manual(values=c("#FF69B4", "#5a5aff"))+
  scale_shape_manual(values = c(21,24,22))+
  xlab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylab(paste0("PC1: ",percentVar[1],"% variance")) +
theme(legend.position = "none")

plot4 = ggplot(PCA_object, aes(PC2, PC3, shape=timepoint,fill = factor(col), colour = factor(outline))) + geom_point( size=4, aes(shape = timepoint)) +
  scale_colour_manual(values=c("#FF69B4", "#5a5aff", "black"))+
  scale_fill_manual(values=c("#FF69B4", "#5a5aff"))+
  scale_shape_manual(values = c(21,24,22))+
  xlab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylab(paste0("PC3: ",percentVar[3],"% variance")) +
theme(legend.position = "none")

plot5 = ggplot(PCA_object, aes(PC3, PC1, shape=timepoint,fill = factor(col), colour = factor(outline))) + geom_point( size=4, aes(shape = timepoint)) +
  scale_colour_manual(values=c("#FF69B4", "#5a5aff", "black"))+
  scale_fill_manual(values=c("#FF69B4", "#5a5aff"))+
  scale_shape_manual(values = c(21,24,22))+
  xlab(paste0("PC3: ",percentVar[3],"% variance")) +
  ylab(paste0("PC1: ",percentVar[1],"% variance")) +
theme(legend.position = "none")

plot6 = ggplot(PCA_object, aes(PC3, PC2, shape=timepoint,fill = factor(col), colour = factor(outline))) + geom_point( size=4, aes(shape = timepoint)) +
  scale_colour_manual(values=c("#FF69B4", "#5a5aff", "black"))+
  scale_fill_manual(values=c("#FF69B4", "#5a5aff"))+
  scale_shape_manual(values = c(21,24,22))+
  xlab(paste0("PC3: ",percentVar[3],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
theme(legend.position = "none")


#now we'll open a figure file and arrange all the objects in the correct order 
png(filename="figures/PCA_rld_3.png", res = 500, unit = "in", width = 7, height = 7)
figure_01 = grid.arrange(empty1, plot3, plot5, plot1, empty2, plot6, plot2, plot4, empty3,  ncol=3, nrow=3, widths=c(5,5,5), heights=c(5,5,5))
dev.off()








#generate 3d scatterplot

#generate new plot characters for scatterplot3d
scatterplot3d_pch = c(rep(21,20),rep(24,20),rep(22,20))

png(filename="figures/3D_PCA_1.png", res = 500, unit = "in", width = 5, height = 5)
s3d <- scatterplot3d(PCA_object$PC1,PCA_object$PC2,PCA_object$PC3, angle = -50, 
                     pch = scatterplot3d_pch, bg = col, color = outline,
                     xlab = "PC 1", ylab = "PC 2", zlab = "PC 3", tick.marks = TRUE)
legend(s3d$xyz.convert(0, -8, -2), pch = c(21,24,22), yjust=0,
       legend = c("1M", "2M", "4M"), cex = 1)
dev.off()