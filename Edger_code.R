
#########################Transcriptome-aligned
setwd("C:/Users/thula/Documents/researchR/htseq_illumina_transcriptome/count/")
getwd()
library(parallel)
library(data.table)
library("edgeR")
count_file_dir <- "C:/Users/thula/Documents/researchR/htseq_illumina_transcriptome/count/"
count_file_names <- list.files(count_file_dir)
count_files <- paste0(count_file_dir, count_file_names)
count_files
# read and merge count files, while removing last 5 lines that correspond to summary stats
count_table_list <- mclapply(count_files, function(x) {fread(x)[1:(.N-5)]})
count_table_list
names(count_table_list) <- sub(".text", "", count_file_names)

# check the dimensions of count files
t(sapply(count_table_list, dim))

# check if the count files have the same gene ordering
# compare gene names of the first file with all
gene_names <- sapply(count_table_list, function(x) x$V1)
as.data.frame(apply(gene_names, 2, function(x) all(x == gene_names[,1])))
gene_names

# add file name as a unique column name for gene counts
for(i in names(count_table_list)) {
  names(count_table_list[[i]])[2] <- i
}

# merge read counts and convert to a matrix
count_matrix <- Reduce(merge, count_table_list)
count_matrix_rownames <- count_matrix$V1
count_matrix_colnames = count_matrix$V2
count_matrix$V1 <- NULL
count_matrix <- as.matrix(count_matrix)

head(count_matrix)
count_matrix_rownames = as.character(count_matrix_rownames)
count_matrix=cbind(count_matrix_rownames,count_matrix)
count_matrix[is.na(count_matrix)] = 0
write.table(count_matrix, "count_matrix", sep = '\t', row.names = F, col.names = T, quote = F)
#coldata=data.frame(sample=c('WT1_out','WT2_out','WT3_out','WT4_out','WT5_out'),condition = c('1AtPOL2b','2AtPOL2b','3AtPOL2b','AtPOL2a','Wild_type'), type = c('paired-end','paired-end','paired-end','paired-end','paired-end') )
#coldata=as.matrix(coldata)
#coldata
dim(count_matrix)



#dds <- DESeqDataSetFromMatrix(countData = count_matrix,
#colData = coldata ,
#design = ~condition)
Treatment <- factor(c("N5n.text","N6n.text","N809n.text","N810n.text","WTn.text"))

x=read.delim(file="count_matrix",header = FALSE, sep = "\t", row.names = 1)
x=as.matrix(x)
class(x) <- "numeric"
x=x[-1,]
colnames(x)=c("AtPOL2b1","AtPOL2a","AtPOL2b2","AtPOL2b3","WT")
d <- DGEList(x,group=Treatment)
d=calcNormFactors(d)
d=estimateGLMCommonDisp(d,method="Pearson")
?estimateGLMCommonDisp
d$common.dispersion
d=estimateGLMTagwiseDisp(d)
d$common.dispersion
et=exactTest(d, pair = c(2,5),dispersion=0.3238099^2)
top=topTags(et,adjust.method="BH")

top=as.data.frame(top$table)
write.table(topTags(et,n=nrow(et$table)), "AtPOL2b3_WT_edger_sig.txt", quote = FALSE, sep = ",")
file=read.csv("AtPOL2b3_WT_edger_sig.txt")
file_N6=as.data.frame(file)


sample_N6 = c("atpol2a-1")
T_G=c("Transcriptome")
file_N6 = cbind(file_N6, sample_N6,T_G)
sample_N809 = c("atpol2b-3")
file_N809 = cbind(file_N809,sample_N809, T_G)
sample_N810 = c("atpol2b-1")
file_N810 = cbind(file_N810,sample_N810, T_G)
sample_N5 = c("atpol2b-2")
file_N5 = cbind(file_N5,sample_N5, T_G)

write.table(res2, "table.txt", sep = ',', row.names = F, col.names = T, quote = F)
MA=plotMA(res2)
addmargins( table( res_sig = res2$padj < .05, res2_sig = res2$padj < .05 ) )

#Filter the DEGs by defining conditions as log2fc >= 1 or log2fc <=-1 and FDR <=0.05
sig=subset(file_N6,file_N6$logFC>=1|file_N6$logFC<=-1)
sig=subset(sig,sig$FDR<=0.05)
write.table(sig, "sig_N6.txt", sep = ',',quote = FALSE)


install.packages("ggdraw")
library(reshape)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(cowplot)

#thresholding
threshold_OE <- file_N810$FDR < 0.01 & (file_N810$logFC >=1 | file_N810$logFC<=-1)
length(which(threshold_OE))
file_N5$threshold <- threshold_OE 
file_N809$threshold <- threshold_OE
file_N810$threshold <- threshold_OE 
file_N6$threshold <- threshold_OE 
## Volcano plot


#volcano plots of samples
N5=ggplot(file_N5) +
  geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x =element_blank())+
  coord_cartesian(ylim = c(0, 15), xlim = c(-30, 30))+
  facet_grid(sample_N5~T_G) 

N809=ggplot(file_N809) +
  geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x =element_blank()) +
  coord_cartesian(ylim = c(0, 15), xlim = c(-30, 30))+
  facet_grid(sample_N809~T_G)

N810=ggplot(file_N810) +
  geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x =element_blank()) +
  coord_cartesian(ylim = c(0, 15), xlim = c(-30, 30))+
  facet_grid(sample_N810~T_G)

N6=ggplot(file_N6) +
  geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x =element_blank()) +
  coord_cartesian(ylim = c(0, 15), xlim = c(-30, 30))+
  facet_grid(sample_N6~T_G)

#plotmatrix(N5,N6,N809,N810)
#labeller=labeller(cyl=c())
c=cowplot::plot_grid(N6,N810,N5,N809, ncol = 1)
mtext("log2 fold change",side=1,line=0,outer=TRUE,cex=1.3)
mtext("-log10 adjusted p-value",side=2,line=0,outer=TRUE,cex=1.3,las=0)
?cowplot::plot_grid

annotate_figure(c,
                left = text_grob("-log10 P value", color = "black", rot = 90),
                bottom = text_grob("log2 fold change", color = "black"))

c$facet
library(ggpubr) 



top=subset(file,(file$logFC>=1|file$logFC<=-1))
top=subset(top,top$FDR<=0.05)
write.table(rownames(top), "AtPOL2a_WT_edger_sig_id.txt", quote = FALSE, sep = ",")
?exactTest



########################################Genome aligned
setwd("C:/Users/thula/Documents/researchR/htseq_illumina_transcriptome/count/")
getwd()
library(parallel)
library(data.table)
library("edgeR")
count_file_dir <- "C:/Users/thula/Documents/researchR/htseq_illumina_transcriptome/count/"
count_file_names <- list.files(count_file_dir)
count_files <- paste0(count_file_dir, count_file_names)
count_files
# read and merge count files, while removing last 5 lines that correspond to summary stats
count_table_list <- mclapply(count_files, function(x) {fread(x)[1:(.N-5)]})
count_table_list
names(count_table_list) <- sub(".txt", "", count_file_names)

# check the dimensions of count files
t(sapply(count_table_list, dim))

# check if the count files have the same gene ordering
# compare gene names of the first file with all
gene_names <- sapply(count_table_list, function(x) x$V1)
as.data.frame(apply(gene_names, 2, function(x) all(x == gene_names[,1])))
gene_names

# add file name as a unique column name for gene counts
for(i in names(count_table_list)) {
  names(count_table_list[[i]])[2] <- i
}

# merge read counts and convert to a matrix
count_matrix <- Reduce(merge, count_table_list)
count_matrix_rownames <- count_matrix$V1
count_matrix_colnames = count_matrix$V2
count_matrix$V1 <- NULL
count_matrix <- as.matrix(count_matrix)

head(count_matrix)
count_matrix_rownames = as.character(count_matrix_rownames)
count_matrix=cbind(count_matrix_rownames,count_matrix)
count_matrix[is.na(count_matrix)] = 0
write.table(count_matrix, "count_matrix", sep = '\t', row.names = F, col.names = T, quote = F)
#coldata=data.frame(sample=c('WT1_out','WT2_out','WT3_out','WT4_out','WT5_out'),condition = c('1AtPOL2b','2AtPOL2b','3AtPOL2b','AtPOL2a','Wild_type'), type = c('paired-end','paired-end','paired-end','paired-end','paired-end') )
#coldata=as.matrix(coldata)
#coldata
dim(count_matrix)



#dds <- DESeqDataSetFromMatrix(countData = count_matrix,
#colData = coldata ,
#design = ~condition)
Treatment <- factor(c("N5.txt","N6.txt","N809.txt","N810.txt","WT.txt"))

x=read.delim(file="count_matrix",header = FALSE, sep = "\t", row.names = 1)
x=as.matrix(x)
class(x) <- "numeric"
x=x[-1,]
colnames(x)=c("AtPOL2b1","AtPOL2a","AtPOL2b2","AtPOL2b3","WT")
d <- DGEList(x,group=Treatment)
d=calcNormFactors(d)
d=estimateGLMCommonDisp(d,method="Pearson")
?estimateGLMCommonDisp
d$common.dispersion
d=estimateGLMTagwiseDisp(d)
d$common.dispersion
et=exactTest(d, pair = c(2,5),dispersion=0.3238099^2)
top=topTags(et,adjust.method="BH")

top=as.data.frame(top$table)
write.table(topTags(et,n=nrow(et$table)), "AtPOL2b3_WT_edger_sig.txt", quote = FALSE, sep = ",")
file=read.csv("AtPOL2b3_WT_edger_sig.txt")
file_N6_G=as.data.frame(file)

#Differential expression analysis using nbiomtest
sample_N6 = c("atpol2a-1")
T_G=c("Genome")
file_N6_G = cbind(file_N6_G, sample_N6,T_G)
sample_N809 = c("atpol2b-3")
file_N809_G = cbind(file_N809_G,sample_N809, T_G)
sample_N810 = c("atpol2b-1")
file_N810_G = cbind(file_N810_G,sample_N810, T_G)
sample_N5 = c("atpol2b-2")
file_N5_G = cbind(file_N5_G,sample_N5, T_G)

write.table(res2, "table.txt", sep = ',', row.names = F, col.names = T, quote = F)
MA=plotMA(res2)
addmargins( table( res_sig = res2$padj < .05, res2_sig = res2$padj < .05 ) )

#Filter the DEGs by defining conditions as log2fc >= 1 or log2fc <=-1 and FDR <=0.05
sig=subset(file_N6_G,file_N6_G$logFC>=1|file_N6_G$logFC<=-1)
sig=subset(sig,sig$FDR<=0.05)
write.table(sig, "sig_N6_G.txt", sep = ',',quote = FALSE)


install.packages("ggdraw")
library(reshape)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(cowplot)

#thresholding
threshold_OE <- file_N810_G$FDR < 0.01 & (file_N810_G$logFC >=1 | file_N810_G$logFC<=-1)
length(which(threshold_OE))
file_N5_G$threshold <- threshold_OE 
file_N809_G$threshold <- threshold_OE
file_N810_G$threshold <- threshold_OE 
file_N6_G$threshold <- threshold_OE 

## Volcano plots of samples
N5_G=ggplot(file_N5_G) +
  geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x =element_blank())+
  coord_cartesian(ylim = c(0, 15), xlim = c(-30, 30))+
  facet_grid(sample_N5~T_G) 

N809_G=ggplot(file_N809_G) +
  geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x =element_blank()) +
  coord_cartesian(ylim = c(0, 15), xlim = c(-30, 30))+
  facet_grid(sample_N809~T_G)

N810_G=ggplot(file_N810_G) +
  geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x =element_blank()) +
  coord_cartesian(ylim = c(0, 15), xlim = c(-30, 30))+
  facet_grid(sample_N810~T_G)

N6_G=ggplot(file_N6_G) +
  geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x =element_blank()) +
  coord_cartesian(ylim = c(0, 15), xlim = c(-30, 30))+
  facet_grid(sample_N6~T_G)

#plotmatrix(N5,N6,N809,N810)
#labeller=labeller(cyl=c())
c=cowplot::plot_grid(N6,N6_G,N810,N810_G,N5,N5_G,N809,N809_G, ncol = 2)
mtext("log2 fold change",side=1,line=0,outer=TRUE,cex=1.3)
mtext("-log10 adjusted p-value",side=2,line=0,outer=TRUE,cex=1.3,las=0)
?cowplot::plot_grid

annotate_figure(c,
                left = text_grob("-log10 P value", color = "black", rot = 90),
                bottom = text_grob("log2 fold change", color = "black"))

c$facet
library(ggpubr) 



top=subset(file,(file$logFC>=1|file$logFC<=-1))
top=subset(top,top$FDR<=0.05)
write.table(rownames(top), "AtPOL2a_WT_edger_sig_id.txt", quote = FALSE, sep = ",")
?exactTest
