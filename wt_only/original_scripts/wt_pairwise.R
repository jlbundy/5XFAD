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
library("gplots")
library("DESeq2")
library("tcltk2")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("RColorBrewer")
library("xlsx")
#library("openxlsx")
#library("WriteXLS")
#library("dataframes2xls")


#set working directory to input file location
#setwd("C:/Users/joseph.bundy/OneDrive/cloud/R/5XFAD/ngs")
#setwd("C:/Users/joseph.bundy.MED/OneDrive/cloud/R/5XFAD/ngs")
setwd("C:/Users/josep/OneDrive/cloud/R/5XFAD/ngs")

in_data_RNA <- read.csv("input/5XFAD_counts_unnorm.csv", header = T, stringsAsFactors=FALSE)
row.names(in_data_RNA) = in_data_RNA[,1]
raw_data_RNA <- in_data_RNA[-c(1:3)]
raw_data_RNA <- data.frame(raw_data_RNA)


meta= read.csv("input/5XFAD_NGS_meta.csv",header = F)
meta = meta[2:60]

raw_data_RNA = cbind(raw_data_RNA[1:10],raw_data_RNA[21:30],raw_data_RNA[41:50])
meta = cbind(meta[1:10],meta[21:30],meta[41:50])

#setwd("C:/Users/joseph.bundy/OneDrive/cloud/R/5XFAD/ngs/WT_only/pairwise")
#setwd("C:/Users/joseph.bundy.MED/OneDrive/cloud/R/5XFAD/ngs/WT_only/pairwise")
setwd("C:/Users/josep/OneDrive/cloud/R/5XFAD/ngs/WT_only/pairwise")

x <- as.matrix(meta[2,])
x = as.factor(x)
x = factor(x, levels = unique(x))

geneinfo<- in_data_RNA[1:3]

# filter for removing genes with few counts
# filter <- apply(raw_data_RNA, 1, function(x) length(x[x>5])>=5)
# filtered <- raw_data_RNA[filter,]
# geneinfo <- in_data_RNA[filter,]
# geneinfo <- geneinfo[2:3]
# genes <- rownames(filtered)#[grep("ENS", rownames(filtered))]



#RUVseq 

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


timepoint=meta[5,]
timepoint = as.matrix(timepoint)
timepoint = as.character(timepoint)


#$$$$$$$$$$$$$$$$$$$$$$$$$$ If you want to include batch effects
# batch=meta[30,]
# batch = as.matrix(batch)
# batch = as.character(batch)


#DESeq2 analysis

SampleName=colnames(raw_data_RNA)
SampleTable = data.frame(SampleName = SampleName, group = group, sex = sex, timepoint = timepoint)

#CountData <- as.matrix(filtered)
CountData <- as.matrix(raw_data_RNA)

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
  design = ~sex + timepoint)
dds_p = ddsFullCountTable
dds_p <- DESeq(dds_p)


alpha = 0.05
independentFiltering = TRUE

normalizedCounts <- t( t(counts(dds)) / sizeFactors(dds) )
colnames(normalizedCounts) = group

means  = data.frame("female_4m" = rowMeans(normalizedCounts[,21:25]),"male_4m" = rowMeans(normalizedCounts[,26:30]))


res = results(dds, contrast=c("group","wt_4m_f","wt_4m_m"), alpha=alpha, independentFiltering = independentFiltering)
table(res$padj > 0.05)
out_file <-  cbind(geneinfo,means,res)
write.csv(out_file,file = "output/wt_4m_m_v_f.csv")
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
# write.xlsx(out_file, file="output/wt_4m_m_v_f.xlsx", sheetName="ALL")
# write.xlsx(DE, file="output/wt_4m_m_v_f.xlsx", sheetName="DE", append=TRUE)
# write.xlsx(UP, file="output/wt_4m_m_v_f.xlsx", sheetName="UP", append=TRUE)
# write.xlsx(DOWN, file="output/wt_4m_m_v_f.xlsx", sheetName="DOWN", append=TRUE)


tiff(filename="figures/wt_4m_m_v_f.tiff", res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-1,1),main="4 month male vs female", alpha = alpha)
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =0.75, label = textup)
text(x = .1, y =-0.75, label = textdown)
text(x = 100000, y =0.9, label = textDE)
dev.off()

res = results(dds, contrast=c("group","wt_2m_f","wt_2m_m"), alpha=alpha, independentFiltering = independentFiltering)
table(res$padj > 0.05)
out_file <-  cbind(geneinfo,res)
write.csv(out_file,file = "output/wt_2m_m_v_f.csv")
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
# write.xlsx(out_file, file="output/wt_2m_m_v_f.xlsx", sheetName="ALL")
# write.xlsx(DE, file="output/wt_2m_m_v_f.xlsx", sheetName="DE", append=TRUE)
# write.xlsx(UP, file="output/wt_2m_m_v_f.xlsx", sheetName="UP", append=TRUE)
# write.xlsx(DOWN, file="output/wt_2m_m_v_f.xlsx", sheetName="DOWN", append=TRUE)

tiff(filename="figures/wt_2m_m_v_f.tiff", res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-1,1),main="2 month male vs female", alpha = alpha)
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =0.75, label = textup)
text(x = .1, y =-0.75, label = textdown)
text(x = 100000, y =0.9, label = textDE)
dev.off()

res = results(dds, contrast=c("group","wt_1m_f","wt_1m_m"), alpha=alpha, independentFiltering = independentFiltering)
table(res$padj > 0.05)
out_file <-  cbind(geneinfo,res)
write.csv(out_file,file = "output/wt_1m_m_v_f.csv")
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
# write.xlsx(out_file, file="output/wt_1m_m_v_f.xlsx", sheetName="ALL")
# write.xlsx(DE, file="output/wt_1m_m_v_f.xlsx", sheetName="DE", append=TRUE)
# write.xlsx(UP, file="output/wt_1m_m_v_f.xlsx", sheetName="UP", append=TRUE)
# write.xlsx(DOWN, file="output/wt_1m_m_v_f.xlsx", sheetName="DOWN", append=TRUE)

tiff(filename="figures/wt_1m_m_v_f.tiff", res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-1,1),main="1 month male vs female", alpha = alpha)
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =0.75, label = textup)
text(x = .1, y =-0.75, label = textdown)
text(x = 100000, y =0.9, label = textDE)
dev.off()

res = results(dds, contrast=c("group","wt_2m_f","wt_1m_f"), alpha=alpha, independentFiltering = independentFiltering)
table(res$padj > 0.05)
out_file <-  cbind(geneinfo,res)
write.csv(out_file,file = "output/wt_1m_v_2m_f.csv")
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
# write.xlsx(out_file, file="output/wt_1m_v_2m_f.xlsx", sheetName="ALL")
# write.xlsx(DE, file="output/wt_1m_v_2m_f.xlsx", sheetName="DE", append=TRUE)
# write.xlsx(UP, file="output/wt_1m_v_2m_f.xlsx", sheetName="UP", append=TRUE)
# write.xlsx(DOWN, file="output/wt_1m_v_2m_f.xlsx", sheetName="DOWN", append=TRUE)

tiff(filename="figures/wt_1m_v_2m_f.tiff", res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-2,2),main="wt_1m_v_2m_f")
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =1.5, label = textup)
text(x = .1, y =-1.5, label = textdown)
text(x = 100000, y =1.9, label = textDE)
dev.off()

res = results(dds, contrast=c("group","wt_4m_f","wt_2m_f"), alpha=alpha, independentFiltering = independentFiltering)
table(res$padj > 0.05)
out_file <-  cbind(geneinfo,res)
write.csv(out_file,file = "output/wt_2m_v_4m_f.csv")
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
# write.xlsx(out_file, file="output/wt_2m_v_4m_f.xlsx", sheetName="ALL")
# write.xlsx(DE, file="output/wt_2m_v_4m_f.xlsx", sheetName="DE", append=TRUE)
# write.xlsx(UP, file="output/wt_2m_v_4m_f.xlsx", sheetName="UP", append=TRUE)
# write.xlsx(DOWN, file="output/wt_2m_v_4m_f.xlsx", sheetName="DOWN", append=TRUE)

tiff(filename="figures/wt_2m_v_4m_f.tiff", res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-2,2),main="wt_2m_v_4m_f")
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =1.5, label = textup)
text(x = .1, y =-1.5, label = textdown)
text(x = 100000, y =1.9, label = textDE)
dev.off()

res = results(dds, contrast=c("group","wt_2m_m","wt_1m_m"), alpha=alpha, independentFiltering = independentFiltering)
table(res$padj > 0.05)
out_file <-  cbind(geneinfo,res)
write.csv(out_file,file = "output/wt_1m_v_2m_m.csv")
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
# write.xlsx(out_file, file="output/wt_1m_v_2m_m.xlsx", sheetName="ALL")
# write.xlsx(DE, file="output/wt_1m_v_2m_m.xlsx", sheetName="DE", append=TRUE)
# write.xlsx(UP, file="output/wt_1m_v_2m_m.xlsx", sheetName="UP", append=TRUE)
# write.xlsx(DOWN, file="output/wt_1m_v_2m_m.xlsx", sheetName="DOWN", append=TRUE)

tiff(filename="figures/wt_1m_v_2m_m.tiff", res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-2,2),main="wt_1m_v_2m_m")
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =1.5, label = textup)
text(x = .1, y =-1.5, label = textdown)
text(x = 100000, y =1.9, label = textDE)
dev.off()

res = results(dds, contrast=c("group","wt_4m_m","wt_2m_m"), alpha=alpha, independentFiltering = independentFiltering)
table(res$padj > 0.05)
out_file <-  cbind(geneinfo,res)
write.csv(out_file,file = "output/wt_2m_v_4m_m.csv")
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
# write.xlsx(out_file, file="output/wt_2m_v_4m_m.xlsx", sheetName="ALL")
# write.xlsx(DE, file="output/wt_2m_v_4m_m.xlsx", sheetName="DE", append=TRUE)
# write.xlsx(UP, file="output/wt_2m_v_4m_m.xlsx", sheetName="UP", append=TRUE)
# write.xlsx(DOWN, file="output/wt_2m_v_4m_m.xlsx", sheetName="DOWN", append=TRUE)

tiff(filename="figures/wt_2m_v_4m_m.tiff", res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-2,2),main="wt_2m_v_4m_m")
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =1.5, label = textup)
text(x = .1, y =-1.5, label = textdown)
text(x = 100000, y =1.9, label = textDE)
dev.off()

res = results(dds_p, contrast=c("timepoint","2 month","1 month"), alpha=alpha, independentFiltering = independentFiltering)
table(res$padj > 0.05)
out_file <-  cbind(geneinfo,res)
write.csv(out_file,file = "output/wt_1m_v_2m_pooled.csv")
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
# write.xlsx(out_file, file="output/wt_1m_v_2m_pooled.xlsx", sheetName="ALL")
# gc()
# write.xlsx(DE, file="output/wt_1m_v_2m_pooled.xlsx", sheetName="DE", append=TRUE)
# gc()
# write.xlsx(UP, file="output/wt_1m_v_2m_pooled.xlsx", sheetName="UP", append=TRUE)
# gc()
# write.xlsx(DOWN, file="output/wt_1m_v_2m_pooled.xlsx", sheetName="DOWN", append=TRUE)

tiff(filename="figures/wt_1m_v_2m_pooled.tiff", res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-1,1),main="1 month vs 2 month", alpha = alpha)
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =0.75, label = textup)
text(x = .1, y =-0.75, label = textdown)
text(x = 100000, y =0.9, label = textDE)
dev.off()

res = results(dds_p, contrast=c("timepoint","4 month","2 month"), alpha=alpha, independentFiltering = independentFiltering)
table(res$padj > 0.05)
out_file <-  cbind(geneinfo,res)
write.csv(out_file,file = "output/wt_2m_v_4m_pooled.csv")
DE = subset(out_file, out_file$padj < 0.05)
UP = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange >0)
DOWN = subset(out_file, out_file$padj < 0.05& out_file$log2FoldChange <0)
# write.xlsx(out_file, file="output/wt_2m_v_4m_pooled.xlsx", sheetName="ALL")
# gc()
# write.xlsx(DE, file="output/wt_2m_v_4m_pooled.xlsx", sheetName="DE", append=TRUE)
# gc()
# write.xlsx(UP, file="output/wt_2m_v_4m_pooled.xlsx", sheetName="UP", append=TRUE)
# gc()
# write.xlsx(DOWN, file="output/wt_2m_v_4m_pooled.xlsx", sheetName="DOWN", append=TRUE)

tiff(filename="figures/wt_2m_v_4m_pooled.tiff", res = 300, unit = "in", width = 5, height = 5)
plotMA(res,ylim=c(-1,1),main="2 month vs 4 month", alpha = alpha)
textup = paste0("??? ", length(UP$Associated.Gene.Name))
textdown =  paste0("??? ", length(DOWN$Associated.Gene.Name))
textDE =  paste0(length(DE$Associated.Gene.Name), " DE")
text(x = .1, y =0.75, label = textup)
text(x = .1, y =-0.75, label = textdown)
text(x = 100000, y =0.9, label = textDE)
dev.off()




png(filename="figures/dispersion_estimate.png", res = 500, unit = "in", width = 5, height = 5)
plotDispEsts(dds)
dev.off()




genes = c("ENSMUSG00000021270","ENSMUSG00000023944","ENSMUSG00000020048","ENSMUSG00000090877","ENSMUSG00000026864",
             "ENSMUSG00000015656","ENSMUSG00000025980","ENSMUSG00000029657","ENSMUSG00000090197")

#box and whisker blots for pairwise comparisons

res = results(dds, contrast=c("group","wt_4m_f","wt_4m_m"), alpha=alpha, independentFiltering = independentFiltering)
table(res$padj > 0.05)
out_file <-  cbind(geneinfo,res)


#list of heat shock factor genes for plotting

for (i in 1: length(genes))
  
{
  
  index = genes[i]  
  
  gene_name = data.frame(out_file[index,][2])    
  
  filename = paste0("figures/",gene_name, ".tiff")
  
  data <- plotCounts(dds, index, 
                     intgroup=c("timepoint","sex"), returnData=TRUE)
  
  data = subset(data, data$timepoint=="4 month")
  
  tiff(filename=filename, res = 300, unit = "in", width = 2.5, height = 2.5)
  print(ggplot(data, aes(x=sex, y=count)) +ggtitle(paste0(gene_name))+
          geom_boxplot(aes(fill = sex)) + scale_fill_manual(values=c("#FFB6C1","#89cff0"))+ stat_smooth(se=FALSE,method="loess")+
          labs(x = "sex",
               y = "normalized counts")+
          theme(legend.position = "none",
                #plot.margin = unit(c(0,0,0,0),"cm"),
                plot.title = element_text(size = rel(.5)), 
                axis.text.x=element_text(colour="black",size = rel(.5)),
                axis.text.y=element_text(colour="black",size = rel(.5)),
                axis.ticks = element_line(colour="black"),
                plot.background = element_blank(), 
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), 
                panel.border = element_blank(), 
                panel.background = element_blank(),
                axis.line = element_line(size=.4),
                axis.title.x = element_text(size = rel(.5), colour = "black"),
                axis.title.y = element_text(size = rel(.5), angle = 90, colour = "black")))
  
  dev.off()
}













#box and whisker plots of HSPs at 4 months


res = results(dds, contrast=c("group","wt_4m_f","wt_4m_m"), alpha=alpha, independentFiltering = independentFiltering)
table(res$padj > 0.05)
out_file <-  cbind(geneinfo,res)



for (i in 1: length(genes))
  
{
  
  index = genes[i]  
  
  gene_name = data.frame(out_file[index,][1])    
  
  filename = paste0("figures/",gene_name, ".tiff")
  
  data <- plotCounts(dds, index, 
                     intgroup=c("timepoint","sex"), returnData=TRUE)
  
  data = subset(data, data$timepoint=="4 month")
  
  tiff(filename=filename, res = 300, unit = "in", width = 2.5, height = 2.5)
  print(ggplot(data, aes(x=sex, y=count)) +ggtitle(paste0(gene_name))+
          geom_boxplot(aes(fill = sex)) + scale_fill_manual(values=c("#FFB6C1","#89cff0"))+
          labs(x = "sex",
               y = "normalized counts")+
          theme(legend.position = "none",
                #plot.margin = unit(c(0,0,0,0),"cm"),
                plot.title = element_text(size = rel(.5)), 
                axis.text.x=element_text(colour="black",size = rel(.5)),
                axis.text.y=element_text(colour="black",size = rel(.5)),
                axis.ticks = element_line(colour="black"),
                plot.background = element_blank(), 
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), 
                panel.border = element_blank(), 
                panel.background = element_blank(),
                axis.line = element_line(size=.4),
                axis.title.x = element_text(size = rel(.5), colour = "black"),
                axis.title.y = element_text(size = rel(.5), angle = 90, colour = "black")))
  
  dev.off()
}








# data <- plotCounts(dds, gene = "ENSMUSG00000027375", 
#                    intgroup=c("sex"), returnData=TRUE)
# 
# png(filename="figures/top_transcript_male_Zfp945.png", res = 500, unit = "in", width = 10, height = 5)
# ggplot(data, aes(x=timepoint, y=count,  group=1)) + 
#   geom_point() + stat_smooth(se=FALSE,method="loess") 
# dev.off()


