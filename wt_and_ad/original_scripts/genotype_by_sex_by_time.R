# R script for analyzing mass spec protomics and/or RNA-seq data
# Joseph L Bundy
# project started 10_9_2013

##################################################################################1
#Block 1 is basic I/O and data structuring 
#generates objects used in later analyses

#import libraries
library("stats")
library("rgl")
library("MASS")
library("DESeq")
library("edgeR")
library("vegan")
library("gplots")
library("DESeq2")
library("tcltk2")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")

#set working directory to input file location
#setwd("C:/Users/joseph.bundy/OneDrive/cloud/R/5XFAD/ngs")
#setwd("C:/Users/joseph.bundy.MED/OneDrive/cloud/R/5XFAD/ngs")
setwd("C:/Users/josep/OneDrive/cloud/R/5XFAD/ngs")


in_data_RNA <- read.csv("input/5XFAD_counts_unnorm_transgenes_added.csv", header = T, stringsAsFactors=FALSE)
row.names(in_data_RNA) = in_data_RNA[,1]
raw_data_RNA <- in_data_RNA[-c(1:3)]
raw_data_RNA <- data.frame(raw_data_RNA)

meta= read.csv("input/5XFAD_NGS_meta.csv",header = F)
meta = data.frame(meta)
meta = meta[,2:61]

#setwd("C:/Users/joseph.bundy/OneDrive/cloud/R/5XFAD/ngs/AD_and_WT/sex_by_genotype_by_time")
#setwd("C:/Users/joseph.bundy.MED/OneDrive/cloud/R/5XFAD/ngs/AD_and_WT/sex_by_genotype_by_time")
setwd("C:/Users/josep/OneDrive/cloud/R/5XFAD/ngs/AD_and_WT/sex_by_genotype_by_time")
#parse out information for specific pairwise comparison
raw_data_RNA = (raw_data_RNA[1:60])



#meta 1 for sex by genotype by timepoint group
#meta 2 for sex
#meta 3 for genotype
#meta 4 for time point

 group=meta[2,]
 group = as.matrix(group)
 group = as.character(group)
 #group = group[1:60]
 

sex=meta[3,]
sex = as.matrix(sex)
sex = as.character(sex)
sex = sex[1:60]


genotype=meta[4,]
genotype = as.matrix(genotype)
genotype = as.character(genotype)


timepoint=meta[5,]
timepoint = as.matrix(timepoint)
timepoint = as.character(timepoint)


# group_2=meta[6,]
# group_2 = as.matrix(group_2)
# group_2 = as.character(group_2)
# group_2 = group_2[1:60]
# 
group_3=meta[7,]
group_3 = as.matrix(group_3)
group_3 = as.character(group_3)
group_3 = group_3[1:60]

batch=meta[30,]
batch = as.matrix(batch)
batch = as.character(batch)

#DESeq2 analysis
SampleName=colnames(raw_data_RNA)
SampleTable = data.frame(SampleName = SampleName,group = group, sex = sex, genotype = genotype, timepoint = timepoint, group_3 = group_3, batch = batch)
#SampleTable = data.frame(SampleName = SampleName, group = group,group_2 = group_2 ,group_3 = group_3)

CountData <- as.matrix(raw_data_RNA)
#CountData <- matrix(CountData, ncol = ncol(CountData), dimnames = NULL)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = CountData,
  colData = SampleTable,
  design = ~sex*genotype*timepoint)
#  design = ~sex+genotype+timepoint+sex:genotype:timepoint)



#ddsFullCountTable$genotype <- relevel(ddsFullCountTable$genotype, "wt")
# dds$genotype <- ordered(dds$genotype, levels = c("wt","ad)"))


# dds <- estimateSizeFactors(dds)
# dds <- estimateDispersions(dds, fitType = "local")
# dds <- nbinomWaldTest(dds)

dds = ddsFullCountTable
dds <- DESeq(dds, test = "LRT", reduced = ~sex+genotype+timepoint+sex:genotype+sex:timepoint+genotype:timepoint,betaPrior=FALSE)
#dds <- DESeq(dds, test = "LRT", reduced = ~sex+genotype+timepoint,betaPrior=FALSE)


attr(dds,"modelMatrix")

res = results(dds, alpha = 0.05)

summary(res)

# normalizedCounts <- t( t(counts(dds)) / sizeFactors(dds) )
# colnames(normalizedCounts) = group

out_file <- cbind(in_data_RNA[2:3], res)

write.csv(out_file,file = "output/all_interactions.csv")




pch = c(rep(21,20),rep(24,20),rep(22,20))
outline = c(rep(c(rep(c(rep("#5a5aff",5),rep("#FF69B4",5)),6))))
col = c(rep(c(rep("#5a5aff",5),rep("#FF69B4",5)),6))


data <- plotCounts(dds, gene= "ENSG00000080815", 
intgroup=c("timepoint","genotype","sex","group_3"), returnData=TRUE)

gene = "PSEN1"

tiff(filename=paste0("figures/", gene, ".tiff"), res = 500, unit = "in", width = 3.5, height = 2)
plot1 = ggplot(data, aes(x=timepoint, y=count, shape=genotype,fill = factor(col), colour = factor(outline), group = group_3)) + 
  geom_point( size=1.5, aes(shape = genotype)) +
 # geom_jitter(width = 0.2,size=4, aes(shape = genotype))+
  guides(fill=guide_legend(title = "sex", override.aes = list(colour = c("#FF69B4", "#5a5aff"))), shape = guide_legend(title = "age"))+
  scale_colour_manual(values=c("#FF69B4", "#5a5aff"), guide = FALSE)+
  ggtitle(gene)+
  geom_errorbar(stat = "summary", fun.y = "mean", width=0.3,aes(ymax=..y..,ymin=..y..))+
  stat_smooth(se=FALSE,method="loess") +
  scale_fill_manual(values=c("#FF69B4", "#5a5aff"), labels = c("females", "males"))+
  scale_shape_manual(values = c(21,24,22))+
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(legend.position=c(1.25, .5),
      plot.title = element_text(colour = "black", face= "italic"))
plot1
dev.off()



png(filename="figures/dispersion_estimate.png", res = 500, unit = "in", width = 5, height = 5)
plotDispEsts(dds)
dev.off()

plotMA(res,ylim=c(-2,2),main="DESeq2")

