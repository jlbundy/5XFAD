# R script for analyzing mass spec protomics and/or RNA-seq data
# Joseph L Bundy
# project started 10_9_2013

#independent filtering left as default

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
# library("tcltk2")
 library("ggplot2")
 library("pheatmap")
library("grid")
# library("RColorBrewer")
# library("PoiClaClu")
library("genefilter")


#set working directory to input file location
setwd("C:/Users/joseph.bundy/OneDrive/cloud/R/5XFAD/ngs")
#setwd("C:/Users/joseph.bundy.MED/OneDrive/cloud/R/5XFAD/ngs")
#setwd("C:/Users/josep/OneDrive/cloud/R/5XFAD/ngs")

in_data_RNA <- read.csv("input/5XFAD_counts_unnorm.csv", header = T, stringsAsFactors=FALSE)
row.names(in_data_RNA) = in_data_RNA[,1]
raw_data_RNA <- in_data_RNA[-c(1:3)]
raw_data_RNA <- data.frame(raw_data_RNA)

meta= read.csv("input/5XFAD_NGS_meta.csv",header = F)
meta = data.frame(meta)
meta = meta[,2:61]

setwd("C:/Users/joseph.bundy/OneDrive/Cloud/R/5XFAD/ngs/WT_only/wt_sex_by_time")
#setwd("C:/Users/joseph.bundy.MED/OneDrive/Cloud/R/5XFAD/ngs/WT_only/wt_sex_by_time")
#setwd("C:/Users/josep/OneDrive/cloud/R/5XFAD/ngs/WT_only/wt_sex_by_time")


raw_data_RNA = cbind(raw_data_RNA[1:10],raw_data_RNA[21:30],raw_data_RNA[41:50])
meta = cbind(meta[1:10],meta[21:30],meta[41:50])

geneinfo <- in_data_RNA[2:3]

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
 #group = group[1:60]
 
 
sex=meta[3,]
sex = as.matrix(sex)
sex = as.character(sex)


timepoint=meta[5,]
timepoint = as.matrix(timepoint)
timepoint = as.character(timepoint)


batch=meta[30,]
batch = as.matrix(batch)
batch = as.character(batch)


#DESeq2 analysis
SampleName=colnames(raw_data_RNA)
SampleTable = data.frame(SampleName = SampleName, sex = sex, timepoint = timepoint, batch = batch)
#SampleTable = data.frame(SampleName = SampleName, group = group,group_2 = group_2 ,group_3 = group_3)

#CountData <- as.matrix(filtered)
CountData <- as.matrix(raw_data_RNA)


ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = CountData,
  colData = SampleTable,

#sex by time interaction
design = ~sex+timepoint+sex:timepoint)
#design = ~sex+timepoint+SV+sex:timepoint)
#design = ~sex+timepoint+batch+sex:timepoint)





dds = ddsFullCountTable


mcols(dds) <- DataFrame(mcols(dds), geneinfo)


#sex_by_time
dds <- DESeq(dds, reduced = ~sex+timepoint, test = "LRT")
#dds <- DESeq(dds, reduced = ~sex+timepoint+SV, test = "LRT")
#dds <- DESeq(dds, reduced = ~sex+timepoint+batch, test = "LRT")



res = results(dds, alpha = 0.05)

# normalizedCounts <- t( t(counts(dds)) / sizeFactors(dds) )
# colnames(normalizedCounts) = group
table(res$padj > 0.05)

out_file <-  cbind(geneinfo,res)
write.csv(out_file,file = "output/sex_by_time.csv")


genes = subset(row.names(out_file),out_file$padj<0.05)
#add heatshock factor genes to list
genes = c(genes, "ENSMUSG00000022556","ENSMUSG00000019878","ENSMUSG00000045802","ENSMUSG00000033249","ENSMUSG00000070345")

#add sex-linked genes to list
genes = c(genes, "ENSMUSG00000068457","ENSMUSG00000037369", "ENSMUSG00000069045", "ENSMUSG00000000787","ENSMUSG00000069049","ENSMUSG00000035150","ENSMUSG00000056673","ENSMUSG00000025332")

for (i in 1: length(genes))
  
{
  
  index = genes[i]  
  
  gene_name = data.frame(out_file[index,][1])    
  
  filename = paste0("figures/",gene_name, ".tiff")
  
  data <- plotCounts(dds, index, 
                     intgroup=c("timepoint","sex"), returnData=TRUE)
  
  
  tiff(filename=filename, res = 300, unit = "in", width = 2.5, height = 1.25)
  print(ggplot(data, aes(x=timepoint, y=count, color=sex, group=sex)) +ggtitle(paste0(gene_name))+
          geom_point(size = 1) + scale_color_manual(values=c("#FF69B4","#3232FF"))+ stat_smooth(se=FALSE,method="loess")+
        labs(x = "age",
             y = "normalized counts")+
          theme(legend.position = "none",
                plot.margin = unit(c(0,0,0,0),"cm"),
                plot.title = element_text(size = rel(.8), face= "italic"), 
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







#change so female and male regressions are on same graph 



SEX="male"

#8 heat-shock proteins
genes = c("ENSMUSG00000029657","ENSMUSG00000090877","ENSMUSG00000091971","ENSMUSG00000026864","ENSMUSG00000004951","ENSMUSG00000032285","ENSMUSG00000015656","ENSMUSG00000036052")

#get HSF1 data for use
HSF1_data <- plotCounts(dds, gene = "ENSMUSG00000022556", 
                        intgroup=c("timepoint","sex"), returnData=TRUE)
colnames(HSF1_data)[1] <- "HSF1"


for (i in 1: length(genes))
  
{
  
  index = genes[i]  
  
  gene_name = data.frame(out_file[index,][1])    
  
  
  
filename = paste0("figures/HSF1_regressions/correlation_HSF1_vs_",gene_name, "_all.tiff")

if (SEX=="female")  {
  filename = paste0("figures/HSF1_regressions/correlation_HSF1_vs_",gene_name, "_female.tiff")
}
if (SEX=="male")  {
  filename = paste0("figures/HSF1_regressions/correlation_HSF1_vs_",gene_name, "_male.tiff")
}  
  

tiff(filename=filename, res = 300, unit = "in", width = 2.5, height = 2.5)
HSP_data <- plotCounts(dds, gene = index, 
intgroup=c("timepoint","sex"), returnData=TRUE)
colnames(HSP_data)[1] <- "HSP"
merged_data <- cbind(HSF1_data, HSP_data)


if (SEX=="female")  {
  merged_data <- subset(merged_data, merged_data$sex =="female")
}
if (SEX=="male")  {
  merged_data <- subset(merged_data, merged_data$sex =="male")
}  


model <- lm(merged_data$HSP~merged_data$HSF1)
test = cor.test(merged_data$HSF1,merged_data$HSP, method = "spearman")
text = paste0("slope == ",sprintf("%.2f",model$coefficients[2]),"~~~R^2 == ", 
                 sprintf("%.2f",summary(model)$r.squared)#, "~~~Rho == ", sprintf("%.2f",test$estimate), "~~~pvalue == ",sprintf("%.4f",summary(model)$coefficients[8], d = 4))
)
if (SEX=="all"){
print(ggplot(merged_data, aes(merged_data$HSF1, merged_data$HSP,merged_data$timepoint)) + 
  geom_point(aes(color=as.factor(merged_data$sex), shape = as.factor(merged_data$timepoint)), size = 2.0) + 
    scale_shape_manual(values = c(16,17,15))+
   scale_color_manual(values = c("#FFB6C1", "#89cff0"))+
   geom_smooth(method = "lm", formula = y~x, se = FALSE, color = "black", fullrange=F, linetype = 2, lwd = 0.9) +
  geom_text(aes(x = 0.83 * max(merged_data$HSF1), y = 1.05 * max(merged_data$HSP), label = text), color="black", size=3.0, parse = TRUE)+


  labs(size = 4,
       x = "Hsf1",
       y = as.character(gene_name),
       title = " ")+
  theme(legend.position = "none",
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        axis.ticks = element_line(colour="black"),
        plot.title = element_blank(), 
        plot.margin = unit(c(.2,.2,.2,.2), units = "lines"),
        axis.title.x = element_text(size = rel(1), colour = "black"),
        axis.title.y = element_text(size = rel(1), angle = 90, colour = "black"),
        axis.text=element_text(size=7),
        plot.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(size=.4)
  ))
}

if (SEX=="female"){
  print(ggplot(merged_data, aes(merged_data$HSF1, merged_data$HSP,merged_data$timepoint)) + 
          geom_point(aes(color=as.factor(merged_data$sex), shape = as.factor(merged_data$timepoint)), size = 2.0) + 
          scale_shape_manual(values = c(16,17,15))+
          scale_color_manual(values = c("#FF69B4"))+
          geom_smooth(method = "lm", formula = y~x, se = FALSE, color = "black", fullrange=F, linetype = 2, lwd = 0.9)+
          geom_text(aes(x = 0.90 * max(merged_data$HSF1), y = 1.05 * max(merged_data$HSP), label = text), color="black", size=3.75, parse = TRUE)+
          labs(size = 4,
               x = "Hsf1",
               y = as.character(gene_name),
               title = " ")+
          theme(legend.position = "none",
                axis.text.x=element_text(size = rel(1.35),colour="black"),
                axis.text.y=element_text(size = rel(1.35),colour="black"),
                axis.ticks = element_line(colour="black"),
                plot.title = element_blank(), 
                plot.margin = unit(c(.2,.2,.2,.2), units = "lines"),
                axis.title.x = element_text(size = rel(1), colour = "black", face= "italic"),
                axis.title.y = element_text(size = rel(1), angle = 90, colour = "black", face= "italic"),
                axis.text=element_text(size=7),
                plot.background = element_blank(), 
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), 
                panel.border = element_blank(), 
                panel.background = element_blank(),
                axis.line = element_line(size=.4)
          ))
}

if (SEX=="male"){
  print(ggplot(merged_data, aes(merged_data$HSF1, merged_data$HSP,merged_data$timepoint)) + 
          geom_point(aes(color=as.factor(merged_data$sex), shape = as.factor(merged_data$timepoint)), size = 2.0) + 
          scale_shape_manual(values = c(16,17,15))+
          scale_color_manual(values = c("#3232FF"))+
          geom_smooth(method = "lm", formula = y~x, se = FALSE, color = "black", fullrange=F, linetype = 2, lwd = 0.9) +
        geom_text(aes(x = 0.82 * max(merged_data$HSF1), y = 1.05 * max(merged_data$HSP), label = text), color="black", size=3.75, parse = TRUE)+
          labs(size = 4,
               x = "Hsf1",
               y = as.character(gene_name),
               title = " ")+
          theme(legend.position = "none",
                axis.text.x=element_text(size = rel(1.35),colour="black"),
                axis.text.y=element_text(size = rel(1.35),colour="black"),
                axis.ticks = element_line(colour="black"),
                plot.title = element_blank(), 
                plot.margin = unit(c(.2,.2,.2,.2), units = "lines"),
                axis.title.x = element_text(size = rel(1), colour = "black", face= "italic"),
                axis.title.y = element_text(size = rel(1), angle = 90, colour = "black", face= "italic"),
                axis.text=element_text(size=7),
                plot.background = element_blank(), 
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), 
                panel.border = element_blank(), 
                panel.background = element_blank(),
                axis.line = element_line(size=.4)
          ))
}



dev.off()
}






SEX = "male"


HSF1_data <- plotCounts(dds, gene = "ENSMUSG00000022556", 
                        intgroup=c("timepoint","sex"), returnData=TRUE)
colnames(HSF1_data)[1] <- "HSF1"

HSBP1_data <- plotCounts(dds, gene = "ENSMUSG00000031839", 
                       intgroup=c("timepoint","sex"), returnData=TRUE)
colnames(HSBP1_data)[1] <- "HSBP1"

merged_data <- cbind(HSF1_data, HSBP1_data)

if (SEX=="female")  {
  merged_data <- subset(merged_data, merged_data$sex =="female")
}
if (SEX=="male")  {
  merged_data <- subset(merged_data, merged_data$sex =="male")
}  


model <- lm(merged_data$HSF1~merged_data$HSBP1)
test = cor.test(merged_data$HSBP1,merged_data$HSF1, method = "spearman")
text = paste0("slope == ",sprintf("%.2f",model$coefficients[2]),"~~~R^2 == ", 
              sprintf("%.2f",summary(model)$r.squared)#, "~~~Rho == ", sprintf("%.2f",test$estimate), "~~~pvalue == ",sprintf("%.4f",summary(model)$coefficients[8], d = 4))
)
plot <- ggplot(merged_data, aes(merged_data$HSBP1, merged_data$HSF1, merged_data$timepoint)) + 
        geom_point(aes(color=as.factor(merged_data$sex), shape = merged_data$timepoint), size = 2.0) + 
         scale_shape_manual(values = c(16,17,15))+ 
          # scale_color_manual(values = c("#FFB6C1", "#89cff0"))+
            scale_color_manual(values = c("#3232FF"))+
           #scale_color_manual(values = c("#FF69B4"))+
        geom_smooth(method = "lm", formula = y~x, se = FALSE, color = "black", fullrange=F, linetype = 2, lwd = 0.9) +
        geom_text(aes(x = 0.83 * max(merged_data$HSBP1), y = 1.05 * max(merged_data$HSF1), label = text), color="black", size=3.75, parse = TRUE)+
        labs(size = 4,
             x = "Hsbp1",
             y = "Hsf1",
             title = " ")+
        theme(legend.position = "none",
              axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black"),
              axis.ticks = element_line(colour="black"),
              plot.title = element_blank(), 
              plot.margin = unit(c(.2,.2,.2,.2), units = "lines"),
              axis.title.x = element_text(size = rel(1), colour = "black", face= "italic"),
              axis.title.y = element_text(size = rel(1), angle = 90, colour = "black", face= "italic"),
              axis.text=element_text(size=7),
              plot.background = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.border = element_blank(), 
              panel.background = element_blank(),
              axis.line = element_line(size=.4))



filename = paste0("figures/HSF1_regressions/correlation_HSF1_vs_HSBP1_all.tiff")

if (SEX=="female")  {
  filename = paste0("figures/HSF1_regressions/correlation_HSF1_vs_HSBP1_female.tiff")
}
if (SEX=="male")  {
  filename = paste0("figures/HSF1_regressions/correlation_HSF1_vs_HSBP1_male.tiff")
}  


tiff(filename=filename, res = 300, unit = "in", width = 2.5, height = 2.5)
plot
dev.off()
































# data <- plotCounts(dds, which.min(res$pvalue), 
# intgroup=c("timepoint","sex"), returnData=TRUE)

data <- plotCounts(dds, gene = "ENSMUSG00000031839", 
intgroup=c("timepoint","sex"), returnData=TRUE)


png(filename="figures/HSBP1.png", res = 500, unit = "in", width = 10, height = 5)

ggplot(data, aes(x=timepoint, y=count, color=sex, group=sex), main = "derp") + 
geom_point() + stat_smooth(se=FALSE,method="loess") +ggtitle(paste0("HSBP1"))+
scale_color_manual(values=c("#FFB6C1", "#89cff0"))+
  
#is this more appropriate to force lower bound to 0?  
scale_y_continuous(expand = c(0,0))+
expand_limits(y = c(0,1.05 * max(data$count)))  
dev.off()




png(filename="figures/dispersion_estimate.png", res = 500, unit = "in", width = 5, height = 5)
plotDispEsts(dds)
dev.off()

plotMA(res,ylim=c(-2,2),main="DESeq2")



# vst <- varianceStabilizingTransformation(dds, blind = FALSE,fitType = "parametric")
# topVarGenes <- head(order(rowVars(assay(vst)),decreasing=TRUE),20)
# mat <- assay(vst)[ topVarGenes, ]
# mat <- mat - rowMeans(mat)
# df <- as.data.frame(colData(vst)[,c("sex","timepoint")])
# pheatmap(mat, annotation_col=df, cluster_col=FALSE)


ann_colors = list(
  sex = c(male = "#89cff0", female = "#FFB6C1"),
  timepoint = c("1 month" = "#7570B3", "2 month" = "#E7298A", "4 month" = "#66A61E")
)


tiff(filename="figures/sex_by_time_heatmap.tiff", res = 300, unit = "in", width =5, height = 2.5)

vst <- varianceStabilizingTransformation(dds, blind = FALSE,fitType = "parametric")
siggenes <- head(order(res$padj),68)
mat <- assay(vst)[ siggenes, ]
names <- mcols(vst)[ siggenes, ]
names <- names[,1]
#names = data.frame(out_file[siggenes,][1])
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vst)[,c("sex","timepoint")])
pheatmap(scale = "row", fontsize = 7,mat, annotation_col=df, cluster_col=FALSE, show_rownames = FALSE,show_colnames = FALSE,  annotation_colors = ann_colors, col=topo.colors(100), border_color = NA)
dev.off()




