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

#set working directory to input file location\
#setwd("C:/Users/joseph.bundy.MED/OneDrive/cloud/R/5XFAD/ngs")
setwd("C:/Users/joseph.bundy/OneDrive/cloud/R/5XFAD/ngs")
#setwd("C:/Users/josep/OneDrive/cloud/R/5XFAD/ngs")


in_data_RNA <- read.csv("input/5XFAD_counts_unnorm_transgenes_added.csv", header = T, stringsAsFactors=FALSE)
raw_data_RNA <- in_data_RNA[-c(1:3)]
raw_data_RNA <- data.frame(raw_data_RNA)
row.names(raw_data_RNA) = in_data_RNA[,1]

meta= read.csv("input/5XFAD_NGS_meta.csv",header = F)
meta = data.frame(meta)
meta = meta[,2:61]

setwd("C:/Users/joseph.bundy/OneDrive/cloud/R/5XFAD/ngs/AD_and_WT/sex_by_genotype")
#setwd("C:/Users/joseph.bundy.MED/OneDrive/cloud/R/5XFAD/ngs/AD_and_WT/sex_by_genotype")
#setwd("C:/Users/josep/OneDrive/cloud/R/5XFAD/ngs/AD_and_WT/sex_by_genotype_by_time")
#parse out information for specific pairwise comparison
raw_data_RNA = (raw_data_RNA[41:60])
meta = meta[41:60]



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



genotype=meta[4,]
genotype = as.matrix(genotype)
genotype = as.character(genotype)


# timepoint=meta[5,]
# timepoint = as.matrix(timepoint)
# timepoint = as.character(timepoint)


# group_2=meta[6,]
# group_2 = as.matrix(group_2)
# group_2 = as.character(group_2)
# group_2 = group_2[1:60]
# 
group_3=meta[7,]
group_3 = as.matrix(group_3)
group_3 = as.character(group_3)


batch=meta[30,]
batch = as.matrix(batch)
batch = as.character(batch)

#DESeq2 analysis
SampleName=colnames(raw_data_RNA)
SampleTable = data.frame(SampleName = SampleName,group = group, sex = sex, genotype = genotype, group_3 = group_3, batch = batch)
#SampleTable = data.frame(SampleName = SampleName, group = group,group_2 = group_2 ,group_3 = group_3)

CountData <- as.matrix(raw_data_RNA)
#CountData <- matrix(CountData, ncol = ncol(CountData), dimnames = NULL)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = CountData,
  colData = SampleTable,
  design = ~genotype+sex+genotype:sex)




#ddsFullCountTable$genotype <- relevel(ddsFullCountTable$genotype, "wt")
# dds$genotype <- ordered(dds$genotype, levels = c("wt","ad)"))


# dds <- estimateSizeFactors(dds)
# dds <- estimateDispersions(dds, fitType = "local")
# dds <- nbinomWaldTest(dds)

dds = ddsFullCountTable
dds <- DESeq(dds, test = "LRT", reduced = ~genotype+sex)


res = results(dds, alpha = 0.05, independentFiltering = TRUE)
summary(res)


normalizedCounts <- t( t(counts(dds)) / sizeFactors(dds) )
colnames(normalizedCounts) = group
out_file <- cbind(in_data_RNA[2:3], res)
write.csv(out_file,file = "output/4m_sex_by_genotype.csv")


# png(filename="figures/Ubc.png", res = 500, unit = "in", width = 10, height = 5)
# 
# data <- plotCounts(dds, gene = "ENSMUSG00000046805", intgroup=c("genotype","sex"), returnData=TRUE)
# ggplot(data, aes(x=sex, y=count, color=group_3, shape = sex,group=group_3)) + 
#   labs(title = "Ubc")+
#   geom_point() +
#   ylab("mRNA expression")+
#   guides(color = guide_legend("condition") )
# dev.off()


genes = c("ENSG00000142192", "ENSMUSG00000020932", "ENSMUSG00000046805","ENSMUSG00000024621")
for (i in 1: length(genes))
{
  index = genes[i]  
  gene_name = data.frame(out_file[index,][1])    
  label = paste0(" pval= ",round(out_file[index,][7]$pvalue, digits = 5), " padj= " ,round(out_file[index,][8]$padj, digits = 5))
  filename = paste0("figures/",gene_name, ".tiff")
  data <- plotCounts(dds, index, 
                     intgroup=c("genotype","sex"), returnData=TRUE)
  
  tiff(filename=filename, res = 300, unit = "in", width = 10, height = 5)
  print(ggplot(data, aes(x=sex, y=count, color=group_3, shape = sex,group=group_3)) +ggtitle(paste0(gene_name,label))+
          geom_point() + scale_color_manual(values=c("#ff1e3f","#4563e7","#FFB6C1","#89cff0")))
  dev.off()
}



png(filename="figures/dispersion_estimate.png", res = 500, unit = "in", width = 5, height = 5)
plotDispEsts(dds)
dev.off()

plotMA(res,ylim=c(-2,2),main="DESeq2")

