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
library("RColorBrewer")
library("xlsx")
library("grid")
library("gridBase")
library("lattice")
library("gridExtra")


options( java.parameters = "-Xmx4g" )
library( "RWeka" )

#set working directory to input file location
#setwd("C:/Users/joseph.bundy.MED/OneDrive/cloud/R/5XFAD/ngs")
#setwd("C:/Users/joseph.bundy/OneDrive/cloud/R/5XFAD/ngs")
setwd("C:/Users/josep/OneDrive/cloud/R/5XFAD/ngs")

in_data_RNA <- read.csv("input/5XFAD_counts_unnorm_transgenes_added.csv", header = T, stringsAsFactors=FALSE)
raw_data_RNA <- in_data_RNA[-c(1:3)]
gene_name <- in_data_RNA[2]

raw_data_RNA <- data.frame(raw_data_RNA)
row.names(raw_data_RNA) = in_data_RNA[,1]

meta= read.csv("input/5XFAD_NGS_meta.csv",header = F)
meta = meta[2:61]

setwd("C:/Users/josep/OneDrive/cloud/R/5XFAD/ngs/AD_and_WT/AD_pairwise")
#setwd("C:/Users/joseph.bundy/OneDrive/cloud/R/5XFAD/ngs/AD_and_WT/AD_pairwise")
#setwd("C:/Users/joseph.bundy.MED/OneDrive/cloud/R/5XFAD/ngs/AD_and_WT/AD_pairwise")


#filtering
# filter <- apply(raw_data_RNA, 1, function(x) length(x[x>5])>=5)
# filtered <- raw_data_RNA[filter,]
# genes <- rownames(filtered)#[grep("ENS", rownames(filtered))]
# 
# in_data_RNA_filtered <- in_data_RNA[filter,]

#RUVseq 

# x <- as.matrix(meta[2,])
# x = as.factor(x)
# x = factor(x, levels = unique(x))

#$$$$$$$$$$$$$$$$$$$$$$$$ If you want to include RUVs
# library("RUVSeq")
# 
# set <- newSeqExpressionSet(as.matrix(filtered),
# phenoData = data.frame(x, row.names=colnames(filtered)))
# 
# differences <- seq(1:15)
# 
# differences <- matrix(differences, ncol = 5, byrow = TRUE)
# set_RUVs <- RUVs(set, genes, k=1, differences, round = TRUE)
# pData(set_RUVs)
# 
# SV = as.numeric(pData(set_RUVs)$W_1)


group=meta[2,]
group = as.matrix(group)
group = as.character(group)

sex=meta[3,]
sex = as.matrix(sex)
sex = as.character(sex)

group_1=meta[6,]
group_1 = as.matrix(group_1)
group_1 = as.character(group_1)

#$$$$$$$$$$$$$$$$$$$$$$$$$$ If you want to include batch effects
# batch=meta[30,]
# batch = as.matrix(batch)
# batch = as.character(batch)


#DESeq2 analysis

SampleName=colnames(raw_data_RNA)
SampleTable = data.frame(SampleName = SampleName, group = group, sex = sex, group_1 = group_1)

CountData <- as.matrix(raw_data_RNA)
#CountData <- matrix(CountData, ncol = ncol(CountData), dimnames = NULL)

ddsFullCountTable <- DESeqDataSetFromMatrix(
countData = CountData,
colData = SampleTable,
design = ~group)
#design = ~group+SV)
#design = ~group+batch)
dds = ddsFullCountTable
dds <- DESeq(dds)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = CountData,
  colData = SampleTable,
  design = ~group_1+sex)
#design = ~group+SV)
#design = ~group+batch)
dds_1 = ddsFullCountTable
dds_1 <- DESeq(dds_1)


alpha = 0.05
independentFiltering = TRUE


normalizedCounts <- t( t(counts(dds)) / sizeFactors(dds) )
normalizedCounts= data.frame(normalizedCounts)

names(normalizedCounts) = SampleName

out_data = data.frame(gene_info, normalizedCounts)
filename = "normalized_counts"
write.csv(out_data,file = paste0("output/",filename, ".csv"))




filename = "wt_v_ad_1m_f"
res = results(dds, contrast=c("group","ad_1m_f","wt_1m_f"), alpha=alpha, independentFiltering =independentFiltering)
summary(res)
out_file <-  cbind(in_data_RNA[2:3],res)
write.csv(out_file,file = paste0("output/",filename, ".csv"))
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
write.xlsx(out_file, file=paste0("output/",filename, ".xlsx"), sheetName="ALL")
gc()
write.xlsx(DE, file=paste0("output/",filename, ".xlsx"), sheetName="DE", append=TRUE)
gc()
write.xlsx(UP, file=paste0("output/",filename, ".xlsx"), sheetName="UP", append=TRUE)
gc()
write.xlsx(DOWN, file=paste0("output/",filename, ".xlsx"), sheetName="DOWN", append=TRUE)

tiff(filename=paste0("figures/", filename,".tiff"), res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-2,2),main=filename)


textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =1.5, label = textup)
text(x = .1, y =-1.5, label = textdown)
text(x = 100000, y =1.9, label = textDE)
dev.off()


filename = "wt_v_ad_1m_m"
res = results(dds, contrast=c("group","ad_1m_m","wt_1m_m"), alpha=alpha, independentFiltering =independentFiltering)
summary(res)
out_file <-  cbind(in_data_RNA[2:3],res)
write.csv(out_file,file = paste0("output/",filename, ".csv"))
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
write.xlsx(out_file, file=paste0("output/",filename, ".xlsx"), sheetName="ALL")
gc()
write.xlsx(DE, file=paste0("output/",filename, ".xlsx"), sheetName="DE", append=TRUE)
gc()
write.xlsx(UP, file=paste0("output/",filename, ".xlsx"), sheetName="UP", append=TRUE)
gc()
write.xlsx(DOWN, file=paste0("output/",filename, ".xlsx"), sheetName="DOWN", append=TRUE)

tiff(filename=paste0("figures/", filename,".tiff"), res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-2,2),main=filename)
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =1.5, label = textup)
text(x = .1, y =-1.5, label = textdown)
text(x = 100000, y =1.9, label = textDE)
dev.off()

filename = "wt_v_ad_2m_f"
res = results(dds, contrast=c("group","ad_2m_f","wt_2m_f"), alpha=alpha, independentFiltering =independentFiltering)
summary(res)
out_file <-  cbind(in_data_RNA[2:3],res)
write.csv(out_file,file = paste0("output/",filename, ".csv"))
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
write.xlsx(out_file, file=paste0("output/",filename, ".xlsx"), sheetName="ALL")
gc()
write.xlsx(DE, file=paste0("output/",filename, ".xlsx"), sheetName="DE", append=TRUE)
gc()
write.xlsx(UP, file=paste0("output/",filename, ".xlsx"), sheetName="UP", append=TRUE)
gc()
write.xlsx(DOWN, file=paste0("output/",filename, ".xlsx"), sheetName="DOWN", append=TRUE)

tiff(filename=paste0("figures/", filename,".tiff"), res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-2,2),main=filename)
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =1.5, label = textup)
text(x = .1, y =-1.5, label = textdown)
text(x = 100000, y =1.9, label = textDE)
dev.off()

filename = "wt_v_ad_2m_m"
res = results(dds, contrast=c("group","ad_2m_m","wt_2m_m"), alpha=alpha, independentFiltering =independentFiltering)
summary(res)
out_file <-  cbind(in_data_RNA[2:3],res)
write.csv(out_file,file = paste0("output/",filename, ".csv"))
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
write.xlsx(out_file, file=paste0("output/",filename, ".xlsx"), sheetName="ALL")
gc()
write.xlsx(DE, file=paste0("output/",filename, ".xlsx"), sheetName="DE", append=TRUE)
gc()
write.xlsx(UP, file=paste0("output/",filename, ".xlsx"), sheetName="UP", append=TRUE)
gc()
write.xlsx(DOWN, file=paste0("output/",filename, ".xlsx"), sheetName="DOWN", append=TRUE)

tiff(filename=paste0("figures/", filename,".tiff"), res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-2,2),main=filename)
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =1.5, label = textup)
text(x = .1, y =-1.5, label = textdown)
text(x = 100000, y =1.9, label = textDE)
dev.off()

filename = "wt_v_ad_4m_f"
res = results(dds, contrast=c("group","ad_4m_f","wt_4m_f"), alpha=alpha, independentFiltering =independentFiltering)
summary(res)
out_file <-  cbind(in_data_RNA[2:3],res)
write.csv(out_file,file = paste0("output/",filename, ".csv"))
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
write.xlsx(out_file, file=paste0("output/",filename, ".xlsx"), sheetName="ALL")
gc()
write.xlsx(DE, file=paste0("output/",filename, ".xlsx"), sheetName="DE", append=TRUE)
gc()
write.xlsx(UP, file=paste0("output/",filename, ".xlsx"), sheetName="UP", append=TRUE)
gc()
write.xlsx(DOWN, file=paste0("output/",filename, ".xlsx"), sheetName="DOWN", append=TRUE)

tiff(filename=paste0("figures/", filename,".tiff"), res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-2,2),main="4 month female WT vs AD")
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =1.5, label = textup)
text(x = .1, y =-1.5, label = textdown)
text(x = 100000, y =1.9, label = textDE)
dev.off()

filename = "wt_v_ad_4m_m"
res = results(dds, contrast=c("group","ad_4m_m","wt_4m_m"), alpha=alpha, independentFiltering =independentFiltering)
summary(res)
out_file <-  cbind(in_data_RNA[2:3],res)
write.csv(out_file,file = paste0("output/",filename, ".csv"))
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
write.xlsx(out_file, file=paste0("output/",filename, ".xlsx"), sheetName="ALL")
gc()
write.xlsx(DE, file=paste0("output/",filename, ".xlsx"), sheetName="DE", append=TRUE)
gc()
write.xlsx(UP, file=paste0("output/",filename, ".xlsx"), sheetName="UP", append=TRUE)
gc()
write.xlsx(DOWN, file=paste0("output/",filename, ".xlsx"), sheetName="DOWN", append=TRUE)

tiff(filename=paste0("figures/", filename,".tiff"), res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-2,2),main="4 month male WT vs AD")
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =1.5, label = textup)
text(x = .1, y =-1.5, label = textdown)
text(x = 100000, y =1.9, label = textDE)
dev.off()



filename = "ad_4m_m_v_f"
res = results(dds, contrast=c("group","ad_4m_f","ad_4m_m"), alpha=alpha, independentFiltering =independentFiltering)
summary(res)
out_file <-  cbind(in_data_RNA[2:3],res)
write.csv(out_file,file = paste0("output/",filename, ".csv"))
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
write.xlsx(out_file, file=paste0("output/",filename, ".xlsx"), sheetName="ALL")
gc()
write.xlsx(DE, file=paste0("output/",filename, ".xlsx"), sheetName="DE", append=TRUE)
gc()
write.xlsx(UP, file=paste0("output/",filename, ".xlsx"), sheetName="UP", append=TRUE)
gc()
write.xlsx(DOWN, file=paste0("output/",filename, ".xlsx"), sheetName="DOWN", append=TRUE)

tiff(filename=paste0("figures/", filename,".tiff"), res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-2,2),main="4 month AD female vs male")
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =1.5, label = textup)
text(x = .1, y =-1.5, label = textdown)
text(x = 100000, y =1.9, label = textDE)
dev.off()

filename = "ad_2m_m_v_f"
res = results(dds, contrast=c("group","ad_2m_f","ad_2m_m"), alpha=alpha, independentFiltering =independentFiltering)
summary(res)
out_file <-  cbind(in_data_RNA[2:3],res)
write.csv(out_file,file = paste0("output/",filename, ".csv"))
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
write.xlsx(out_file, file=paste0("output/",filename, ".xlsx"), sheetName="ALL")
gc()
write.xlsx(DE, file=paste0("output/",filename, ".xlsx"), sheetName="DE", append=TRUE)
gc()
write.xlsx(UP, file=paste0("output/",filename, ".xlsx"), sheetName="UP", append=TRUE)
gc()
write.xlsx(DOWN, file=paste0("output/",filename, ".xlsx"), sheetName="DOWN", append=TRUE)

tiff(filename=paste0("figures/", filename,".tiff"), res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-2,2),main=filename)
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =1.5, label = textup)
text(x = .1, y =-1.5, label = textdown)
text(x = 100000, y =1.9, label = textDE)
dev.off()

filename = "ad_1m_m_v_f"
res = results(dds, contrast=c("group","ad_1m_f","ad_1m_m"), alpha=alpha, independentFiltering =independentFiltering)
summary(res)
out_file <-  cbind(in_data_RNA[2:3],res)
write.csv(out_file,file = paste0("output/",filename, ".csv"))
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
write.xlsx(out_file, file=paste0("output/",filename, ".xlsx"), sheetName="ALL")
gc()
write.xlsx(DE, file=paste0("output/",filename, ".xlsx"), sheetName="DE", append=TRUE)
gc()
write.xlsx(UP, file=paste0("output/",filename, ".xlsx"), sheetName="UP", append=TRUE)
gc()
write.xlsx(DOWN, file=paste0("output/",filename, ".xlsx"), sheetName="DOWN", append=TRUE)

tiff(filename=paste0("figures/", filename,".tiff"), res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-2,2),main=filename)
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =1.5, label = textup)
text(x = .1, y =-1.5, label = textdown)
text(x = 100000, y =1.9, label = textDE)
dev.off()

filename = "ad_1m_v_2m_f"
res = results(dds, contrast=c("group","ad_2m_f","ad_1m_f"), alpha=alpha, independentFiltering =independentFiltering)
summary(res)
out_file <-  cbind(in_data_RNA[2:3],res)
write.csv(out_file,file = paste0("output/",filename, ".csv"))
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
write.xlsx(out_file, file=paste0("output/",filename, ".xlsx"), sheetName="ALL")
gc()
write.xlsx(DE, file=paste0("output/",filename, ".xlsx"), sheetName="DE", append=TRUE)
gc()
write.xlsx(UP, file=paste0("output/",filename, ".xlsx"), sheetName="UP", append=TRUE)
gc()
write.xlsx(DOWN, file=paste0("output/",filename, ".xlsx"), sheetName="DOWN", append=TRUE)

tiff(filename=paste0("figures/", filename,".tiff"), res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-2,2),main=filename)
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =1.5, label = textup)
text(x = .1, y =-1.5, label = textdown)
text(x = 100000, y =1.9, label = textDE)
dev.off()

filename = "ad_2m_v_4m_f"
res = results(dds, contrast=c("group","ad_4m_f","ad_2m_f"), alpha=alpha, independentFiltering =independentFiltering)
summary(res)
out_file <-  cbind(in_data_RNA[2:3],res)
write.csv(out_file,file = paste0("output/",filename, ".csv"))
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
write.xlsx(out_file, file=paste0("output/",filename, ".xlsx"), sheetName="ALL")
gc()
write.xlsx(DE, file=paste0("output/",filename, ".xlsx"), sheetName="DE", append=TRUE)
gc()
write.xlsx(UP, file=paste0("output/",filename, ".xlsx"), sheetName="UP", append=TRUE)
gc()
write.xlsx(DOWN, file=paste0("output/",filename, ".xlsx"), sheetName="DOWN", append=TRUE)

tiff(filename=paste0("figures/", filename,".tiff"), res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-2,2),main=filename)
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =1.5, label = textup)
text(x = .1, y =-1.5, label = textdown)
text(x = 100000, y =1.9, label = textDE)
dev.off()

filename = "ad_1m_v_2m_m"
res = results(dds, contrast=c("group","ad_2m_m","ad_1m_m"), alpha=alpha, independentFiltering =independentFiltering)
summary(res)
out_file <-  cbind(in_data_RNA[2:3],res)
write.csv(out_file,file = paste0("output/",filename, ".csv"))
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
write.xlsx(out_file, file=paste0("output/",filename, ".xlsx"), sheetName="ALL")
gc()
write.xlsx(DE, file=paste0("output/",filename, ".xlsx"), sheetName="DE", append=TRUE)
gc()
write.xlsx(UP, file=paste0("output/",filename, ".xlsx"), sheetName="UP", append=TRUE)
gc()
write.xlsx(DOWN, file=paste0("output/",filename, ".xlsx"), sheetName="DOWN", append=TRUE)

tiff(filename=paste0("figures/", filename,".tiff"), res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-2,2),main=filename)
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =1.5, label = textup)
text(x = .1, y =-1.5, label = textdown)
text(x = 100000, y =1.9, label = textDE)
dev.off()

filename = "ad_2m_v_4m_m"
res = results(dds, contrast=c("group","ad_4m_m","ad_2m_m"), alpha=alpha, independentFiltering =independentFiltering)
summary(res)
out_file <-  cbind(in_data_RNA[2:3],res)
write.csv(out_file,file = paste0("output/",filename, ".csv"))
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
write.xlsx(out_file, file=paste0("output/",filename, ".xlsx"), sheetName="ALL")
gc()
write.xlsx(DE, file=paste0("output/",filename, ".xlsx"), sheetName="DE", append=TRUE)
gc()
write.xlsx(UP, file=paste0("output/",filename, ".xlsx"), sheetName="UP", append=TRUE)
gc()
write.xlsx(DOWN, file=paste0("output/",filename, ".xlsx"), sheetName="DOWN", append=TRUE)

tiff(filename=paste0("figures/", filename,".tiff"), res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-2,2),main=filename)
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =1.5, label = textup)
text(x = .1, y =-1.5, label = textdown)
text(x = 100000, y =1.9, label = textDE)
dev.off()





filename = "wt_v_ad_1m_pooled"
res = results(dds_1, contrast=c("group_1","ad_1m","wt_1m"), alpha=alpha, independentFiltering =independentFiltering)
summary(res)
out_file <-  cbind(in_data_RNA[2:3],res)
write.csv(out_file,file = paste0("output/",filename, ".csv"))
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
write.xlsx(out_file, file=paste0("output/",filename, ".xlsx"), sheetName="ALL")
gc()
write.xlsx(DE, file=paste0("output/",filename, ".xlsx"), sheetName="DE", append=TRUE)
gc()
write.xlsx(UP, file=paste0("output/",filename, ".xlsx"), sheetName="UP", append=TRUE)
gc()
write.xlsx(DOWN, file=paste0("output/",filename, ".xlsx"), sheetName="DOWN", append=TRUE)

tiff(filename=paste0("figures/", filename,".tiff"), res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-2,2),main="1 month WT vs AD", alpha = alpha)
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =1.5, label = textup)
text(x = .1, y =-1.5, label = textdown)
text(x = 100000, y =1.9, label = textDE)
dev.off()

filename = "wt_v_ad_2m_pooled"
res = results(dds_1, contrast=c("group_1","ad_2m","wt_2m"), alpha=alpha, independentFiltering =independentFiltering)
summary(res)
out_file <-  cbind(in_data_RNA[2:3],res)
write.csv(out_file,file = paste0("output/",filename, ".csv"))
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
write.xlsx(out_file, file=paste0("output/",filename, ".xlsx"), sheetName="ALL")
gc()
write.xlsx(DE, file=paste0("output/",filename, ".xlsx"), sheetName="DE", append=TRUE)
gc()
write.xlsx(UP, file=paste0("output/",filename, ".xlsx"), sheetName="UP", append=TRUE)
gc()
write.xlsx(DOWN, file=paste0("output/",filename, ".xlsx"), sheetName="DOWN", append=TRUE)

tiff(filename=paste0("figures/", filename,".tiff"), res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-2,2),main="2 month WT vs AD", alpha = alpha)
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =1.5, label = textup)
text(x = .1, y =-1.5, label = textdown)
text(x = 100000, y =1.9, label = textDE)
dev.off()

filename = "wt_v_ad_4m_pooled"
res = results(dds_1, contrast=c("group_1","ad_4m","wt_4m"), alpha=alpha, independentFiltering =independentFiltering)
summary(res)
out_file <-  cbind(in_data_RNA[2:3],res)
write.csv(out_file,file = paste0("output/",filename, ".csv"))
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
write.xlsx(out_file, file=paste0("output/",filename, ".xlsx"), sheetName="ALL")
gc()
write.xlsx(DE, file=paste0("output/",filename, ".xlsx"), sheetName="DE", append=TRUE)
gc()
write.xlsx(UP, file=paste0("output/",filename, ".xlsx"), sheetName="UP", append=TRUE)
gc()
write.xlsx(DOWN, file=paste0("output/",filename, ".xlsx"), sheetName="DOWN", append=TRUE)

tiff(filename=paste0("figures/", filename,".tiff"), res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-2,2),main="4 month WT vs AD", alpha = alpha)
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =1.5, label = textup)
text(x = .1, y =-1.5, label = textdown)
text(x = 100000, y =1.9, label = textDE)
dev.off()


res.shr <- lfcShrink(dds=dds_1, coef=2, res=res)
res.shr <- lfcShrink(dds=dds_1, contrast=c("group_1","ad_4m","wt_4m"), res=res)


tiff(filename=paste0("figures/", "TOC_figure",".tiff"), unit = "px", width = 300, height = 300)
plotMA(res.shr,ylim=c(-2,2),main="4 month WT vs AD", alpha = alpha)
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =1.5, label = textup)
text(x = .1, y =-1.5, label = textdown)
text(x = 100000, y =1.9, label = textDE)
dev.off()


png(filename="figures/dispersion_estimate.png", res = 500, unit = "in", width = 5, height = 5)
plotDispEsts(dds)
dev.off()

# png(filename="figures/MA_plot.png", res = 500, unit = "in", width = 5, height = 5)
# plotMA(res,ylim=c(-2,2),main="DESeq2")
# dev.off()

data <- plotCounts(dds, gene="ENSG00000142192", 
intgroup=c("group", "group_1","sex", "SampleName"), returnData=TRUE)




#compare fold changes for DE transcripts at 4 months of age with their fold changes at only 2 months of age

wt_vs_ad_4m = data.frame(results(dds_1, contrast=c("group_1","ad_4m","wt_4m"), alpha=alpha, independentFiltering =independentFiltering))
wt_vs_ad_2m = data.frame(results(dds_1, contrast=c("group_1","ad_2m","wt_2m"), alpha=alpha, independentFiltering =independentFiltering))

wt_vs_ad_4m = data.frame(in_data_RNA[2:3],wt_vs_ad_4m )

wt_vs_ad_4m = subset(wt_vs_ad_4m, wt_vs_ad_4m$padj < 0.05)
merged = merge(wt_vs_ad_4m, wt_vs_ad_2m, by = 0)

recalculated_2m_padj = p.adjust(merged$pvalue.y, method = "fdr")
merged = data.frame(merged, recalculated_2m_padj)
length(subset(merged, merged$padj.y<0.05))
length(subset(merged, merged$recalculated_2m_padj<0.05))

write.csv(merged,file = paste0("output/","merged_data_fdr_readjusted_2m", ".csv"))

merged_agree = subset(merged, merged$log2FoldChange.x>0&merged$log2FoldChange.y>0|merged$log2FoldChange.x<0&merged$log2FoldChange.y<0)


model = lm(merged$log2FoldChange.y~merged$log2FoldChange.x)
test = cor.test(merged$log2FoldChange.x, merged$log2FoldChange.y, method = "pearson")

label = paste0("R = ", sprintf(" %.2f",test$estimate))



plot_01 <- ggplot(merged, aes(merged$log2FoldChange.x,merged$log2FoldChange.y )) + 
  geom_point(aes(color=ifelse(merged$log2FoldChange.x>0&merged$log2FoldChange.y>0|merged$log2FoldChange.x<0&merged$log2FoldChange.y<0, "green", 
                              ifelse(merged$log2FoldChange.x!=0&merged$log2FoldChange.y==0, "red", "black"))), size = 1.2) + 
  scale_color_manual(values = c("red", "blue","black"))+
  # geom_smooth(method = "loess", se = FALSE, colour = "green",  lwd = 1.2, span = 0.6) +
  geom_smooth(method = "lm", formula = y~x, se = FALSE, color = "black", fullrange=F, lwd = 1.2, linetype = 2) +
  geom_text(aes(x = 0.0, y = 0.7, label = label), color="black", size=5, parse = FALSE)+
  scale_x_continuous(lim = c(-1,1))+
  scale_y_continuous(lim = c(-1,1))+
  labs(x = "4 month WT vs AD",
       y = "2 month WT vs AD",
       title = element_blank())+
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0),"cm"),
        axis.title.x =  element_text(colour="black"),
        axis.title.y =  element_text(colour="black"), 
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.text=element_text(size=8),
        plot.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(size=.4)
  )


tiff(filename = "figures/2m_4m_FC.tiff", width = 4, height = 4, units = "in", res = 600)
plot_01
dev.off()








