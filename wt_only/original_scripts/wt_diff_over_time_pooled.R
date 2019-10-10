# R script for analyzing mass spec protomics and/or RNA-seq data
# Joseph L Bundy
# project started 10_9_2013

##################################################################################1
#Block 1 is basic I/O and data structuring 
#generates objects used in later analyses

.libPaths(new = "C:/Program Files/RRO/R-3.2.2/library")

#import libraries
# library("stats")
# library("rgl")
# library("MASS")
# library("DESeq")
# library("edgeR")
# library("vegan")
# library("gplots")
library("DESeq2")
library("grid")
# library("tcltk2")
 library("ggplot2")
 library("pheatmap")
# library("RColorBrewer")
# library("PoiClaClu")


#set working directory to input file location
#setwd("C:/Users/joseph.bundy/OneDrive/cloud/R/5XFAD/ngs")
#setwd("C:/Users/joseph.bundy.MED/OneDrive/cloud/R/5XFAD/ngs")
setwd("C:/Users/josep/OneDrive/Cloud/R/5XFAD/ngs")

in_data_RNA <- read.csv("input/5XFAD_counts_unnorm.csv", header = T, stringsAsFactors=FALSE)
row.names(in_data_RNA) = in_data_RNA[,1]
raw_data_RNA <- in_data_RNA[-c(1:3)]
raw_data_RNA <- data.frame(raw_data_RNA)

meta= read.csv("input/5XFAD_NGS_meta.csv",header = F)
meta = data.frame(meta)
meta = meta[,2:61]

raw_data_RNA = cbind(raw_data_RNA[1:10],raw_data_RNA[21:30],raw_data_RNA[41:50])
meta = cbind(meta[1:10],meta[21:30],meta[41:50])

#setwd("C:/Users/joseph.bundy/OneDrive/cloud/R/5XFAD/ngs/WT_only/wt_diff_over_time")
#setwd("C:/Users/joseph.bundy.MED/OneDrive/cloud/R/5XFAD/ngs/WT_only/wt_diff_over_time")
setwd("C:/Users/josep/OneDrive/Cloud/R/5XFAD/ngs/WT_only/wt_diff_over_time")

geneinfo<- in_data_RNA[2:3]

# filter <- apply(raw_data_RNA, 1, function(x) length(x[x>5])>=5)
# filtered <- raw_data_RNA[filter,]
# geneinfo <- in_data_RNA[filter,]
# geneinfo <- geneinfo[2:3]
# genes <- rownames(filtered)#[grep("ENS", rownames(filtered))]


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

stage=meta[9,]
stage = as.matrix(stage)
stage = as.character(stage)

batch=meta[30,]
batch = as.matrix(batch)
batch = as.character(batch)

#DESeq2 analysis
SampleName=colnames(raw_data_RNA)
SampleTable = data.frame(SampleName = SampleName, timepoint = timepoint, sex = sex)
#SampleTable = data.frame(SampleName = SampleName, group = group,group_2 = group_2 ,group_3 = group_3)

#CountData <- as.matrix(filtered)
CountData <- as.matrix(raw_data_RNA)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = CountData,
  colData = SampleTable,

  
#differences over time
design = ~timepoint+sex)
#design = ~timepoint+SV)
#design = ~timepoint+batch)

dds = ddsFullCountTable

mcols(dds) <- DataFrame(mcols(dds), geneinfo)

#sex_over_time
dds <- DESeq(dds, reduced = ~sex, test = "LRT")
#dds <- DESeq(dds, reduced = ~SV, test = "LRT")
#dds <- DESeq(dds, reduced = ~batch, test = "LRT")


res = results(dds, alpha = 0.05, independentFiltering = TRUE)
table(res$padj > 0.05)


# #edgeR analysis (regular)
# edgeR_data <- data.frame(raw_data_RNA)
# edgeR_object <- DGEList(counts = edgeR_data, group = meta)
# edgeR_object <- calcNormFactors(edgeR_object)
# #edgeR_object <- estimateDisp(edgeR_object, span = .15)
#
# edgeR_object <- estimateCommonDisp(edgeR_object, rowsum.filter = 1, verbose = TRUE)
# edgeR_object <-estimateTagwiseDisp(edgeR_object, verbose = TRUE)
# et <- exactTest(edgeR_object)
# er.FDR.list <-(p.adjust(et$table$PValue,method = "BH"))

normalizedCounts <- t( t(counts(dds)) / sizeFactors(dds) )
colnames(normalizedCounts) = group

means  = data.frame(month1 = rowMeans(normalizedCounts[,1:10]),month2 = rowMeans(normalizedCounts[,11:20]),month4 = rowMeans(normalizedCounts[,21:30]))

out_file <-  cbind(geneinfo, means,res)



write.csv(out_file,file = "output/wt_over_time_pooled.csv")



genes = c("ENSMUSG00000026879","ENSMUSG00000038301", "ENSMUSG00000038319","ENSMUSG00000024810", "ENSMUSG00000020914", "ENSMUSG00000032517" )

for (i in 1: length(genes))
  
{
  
  index = genes[i]  
  
  gene_name = data.frame(out_file[index,][1])    
  
  filename = paste0("figures/",gene_name, ".tiff")
  
  data <- plotCounts(dds, index, 
                     intgroup=c("timepoint","sex"), returnData=TRUE)
  
  
  tiff(filename=filename, res = 300, unit = "in", width = 5, height = 1.25)
  print(ggplot(data, color = "black",aes(x=timepoint, y=count, group=1)) +ggtitle(paste0(gene_name))+
          geom_point() + stat_smooth(se=FALSE,method="loess") +
          labs(x = element_blank(),
               y = "normalized counts")+
          theme(legend.position = "none",
                plot.margin = unit(c(0,0,0,0),"cm"),
                plot.title = element_text(size = rel(.8)), 
                panel.grid.minor = element_blank(), 
                panel.border = element_blank(), 
                panel.background = element_blank(),
                panel.grid.major.y = element_line(colour = "black", size = 0.1, linetype = "dashed"),
                axis.text.x=element_blank(),
                axis.text.y=element_text(colour="black",size = rel(.7)),
                axis.ticks = element_blank(), 
                axis.title.x = element_blank(),
                axis.title.y = element_text(size = rel(.7), angle = 90, colour = "black")))
  dev.off()
}




png(filename="figures/dispersion_estimate.png", res = 500, unit = "in", width = 5, height = 5)
plotDispEsts(dds)
dev.off()

plotMA(res,ylim=c(-2,2),main="DESeq2")









ann_colors = list(
  sex = c(male = "#89cff0", female = "#FFB6C1"),
  timepoint = c("1 month" = "#7570B3", "2 month" = "#E7298A", "4 month" = "#66A61E")
)


tiff(filename="figures/change_over_time_heatmap.tiff", res = 300, unit = "in", width =5, height = 2.5)

vst <- varianceStabilizingTransformation(dds, blind = FALSE,fitType = "parametric")
siggenes <- head(order(res$padj),50)
mat <- assay(vst)[ siggenes, ]
names <- mcols(vst)[ siggenes, ]
names <- names[,1]
#names = data.frame(out_file[siggenes,][1])
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vst)[,c("sex","timepoint")])
#pheatmap(scale = "row", fontsize = 4.5,mat, annotation_col=df, cluster_col=FALSE, show_colnames = FALSE, labels_row = names, annotation_colors = ann_colors, col=topo.colors(100))
pheatmap(scale = "row", fontsize = 7,mat, annotation_col=df, cluster_col=FALSE, show_rownames = FALSE,show_colnames = FALSE,  annotation_colors = ann_colors, col=topo.colors(100), border_color = NA)
dev.off()



