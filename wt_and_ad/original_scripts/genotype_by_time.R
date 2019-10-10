# R script for analyzing mass spec protomics and/or RNA-seq data
# Joseph L Bundy
# project started 10_9_2013

##################################################################################1
#Block 1 is basic I/O and data structuring 
#generates objects used in later analyses

Females = FALSE

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
setwd("C:/Users/josep/OneDrive/Cloud/R/5XFAD/ngs")

in_data_RNA <- read.csv("input/5XFAD_counts_unnorm.csv", header = T, stringsAsFactors=FALSE)
raw_data_RNA <- in_data_RNA[-c(1:3)]
raw_data_RNA <- data.frame(raw_data_RNA)

meta= read.csv("input/5XFAD_NGS_meta.csv",header = F)
meta = data.frame(meta)
meta = meta[,2:61]

setwd("C:/Users/josep/OneDrive/Cloud/R/5XFAD/ngs/AD_and_WT/genotype_by_time_interaction")
#setwd("C:/Users/joseph.bundy/OneDrive/Cloud/R/5XFAD/ngs/AD_and_WT/genotype_by_time_interaction")

if (Females == TRUE){
  
  #females
  raw_data_RNA = cbind(raw_data_RNA[1:5],raw_data_RNA[11:15],raw_data_RNA[21:25],raw_data_RNA[31:35],raw_data_RNA[41:45],raw_data_RNA[51:55])
  meta = cbind(meta[1:5],meta[11:15],meta[21:25],meta[31:35],meta[41:45],meta[51:55])
} else{
  
  #males
  raw_data_RNA = cbind(raw_data_RNA[6:10],raw_data_RNA[16:20],raw_data_RNA[26:30],raw_data_RNA[36:40],raw_data_RNA[46:50],raw_data_RNA[56:60])
  meta = cbind(meta[6:10],meta[16:20],meta[26:30],meta[36:40],meta[46:50],meta[56:60])
}


#filtering
filter <- apply(raw_data_RNA, 1, function(x) length(x[x>5])>=5)
filtered <- raw_data_RNA[filter,]
genes <- rownames(filtered)#[grep("ENS", rownames(filtered))]

in_data_RNA_filtered <- in_data_RNA[filter,]


group=meta[2,]
group = as.matrix(group)
group = as.character(group)

genotype=meta[4,]
genotype = as.matrix(genotype)
genotype = as.character(genotype)
 
stage=meta[9,]
stage = as.matrix(stage)
stage = as.character(stage)


genotype=meta[4,]
genotype = as.matrix(genotype)
genotype = as.character(genotype)


timepoint=meta[5,]
timepoint = as.matrix(timepoint)
timepoint = as.character(timepoint)


#DESeq2 analysis
SampleName=colnames(raw_data_RNA)
SampleTable = data.frame(SampleName = SampleName, genotype = genotype, timepoint = timepoint)


CountData <- as.matrix(filtered)
#CountData <- matrix(CountData, ncol = ncol(CountData), dimnames = NULL)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = CountData,
  colData = SampleTable,
  design = ~genotype+timepoint+genotype:timepoint)

# dds <- estimateSizeFactors(dds)
# dds <- estimateDispersions(dds, fitType = "local")
# dds <- nbinomWaldTest(dds)

dds = ddsFullCountTable
dds <- DESeq(dds, test = "LRT", reduced = ~genotype+timepoint)


alpha = 0.05
res = results(dds, alpha = alpha, independentFiltering = FALSE)

table(res$padj > 0.05)

normalizedCounts <- t( t(counts(dds)) / sizeFactors(dds) )
colnames(normalizedCounts) = group

out_file <- cbind(in_data_RNA_filtered[2:3],normalizedCounts, res)

if (Females == TRUE){
  write.csv(out_file,file = "output/female_genotype_by_time.csv")
} else{
  write.csv(out_file,file = "output/male_genotype_by_time.csv")
}





data <- plotCounts(dds, which.min(res$pvalue),
                   
#data <- plotCounts(dds, 21352,            
intgroup=c("timepoint","genotype"), returnData=TRUE)

png(filename="figures/top_transcript_female_ns.png", res = 500, unit = "in", width = 10, height = 5)
ggplot(data, aes(x=timepoint, y=count, color=genotype, group=genotype)) + 
  geom_point() + stat_smooth(se=FALSE,method="loess") 
dev.off()



png(filename="figures/dispersion_estimate.png", res = 500, unit = "in", width = 5, height = 5)
plotDispEsts(dds)
dev.off()

plotMA(res,ylim=c(-2,2),main="DESeq2")





