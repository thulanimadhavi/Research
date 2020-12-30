
library("data.table")
library(ggplot2)
library(parallel)
setwd("C:/Users/thula/Documents/researchR/htseq_illumina_transcriptome/count/")

#Defining directory containg the expression data
directory <- "C:/Users/thula/Documents/researchR/htseq_illumina_transcriptome/count/"
sampleFiles <- grep(".text",list.files(directory),value=TRUE)
#Defining samples
Treatment <- factor(c("N5n.text","N6n.text","N809n.text","N810n.text","WTn.text"))

sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = Treatment)

library( "DESeq" )
#Creating the count dataset (matrix)
cds = newCountDataSetFromHTSeqCount(sampleTable,directory = directory)

#Normalization of data
cds = estimateSizeFactors( cds )
sizeFactors( cds )
cds2 = cds[ ,c( "N809n.text","N810n.text","N5n.text","N6n.text", "WTn.text" ) ]

#Dispersion across samples
cds2 = estimateDispersions( cds2, method="blind", sharingMode="fit-only" )

#Differential expression analysis using nbiomtest
sample_N6 = c("atpol2a-1")
T_G=c("Transcriptome")
res2_N6 = cbind(res2_N6, sample_N6,T_G)
sample_N809 = c("atpol2b-3")
res2_N809 = cbind(res2_N809,sample_N809, T_G)
sample_N810 = c("atpol2b-1")
res2_N810 = cbind(res2_N810,sample_N810, T_G)
sample_N5 = c("atpol2b-2")
res2_N5 = cbind(res2_N5,sample_N5, T_G)
res2_N5 = nbinomTest( cds2, "N5n.text", "WTn.text" )
write.table(res2, "table.txt", sep = ',', row.names = F, col.names = T, quote = F)
MA=plotMA(res2)
addmargins( table( res_sig = res2$padj < .05, res2_sig = res2$padj < .05 ) )

#Filter the DEGs by defining conditions as log2fc >= 2 or log2fc <=-2 and FDR <=0.05
sig=subset(res2_N5,res2_N5$log2FoldChange>=1|res2_N5$log2FoldChange<=-1)
sig=subset(sig,sig$padj<=0.05)
write.table(sig, "sig_N5.txt", sep = ',',quote = FALSE)


install.packages("ggdraw")
library(reshape)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(cowplot)

threshold_OE <- res2_N6$padj < 0.01 & (res2_N6$log2FoldChange>=1|res2_N6$log2FoldChange<=-1)
length(which(threshold_OE))
res2_N5$threshold <- threshold_OE 
res2_N809$threshold <- threshold_OE
res2_N810$threshold <- threshold_OE 
res2_N6$threshold <- threshold_OE 
## Volcano plot

ggplot(res2) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  ggtitle("N501413 vs WT") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  

par(mfrow=c(2,2),oma=c(2,2,0,0))

N5=ggplot(res2_N5) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x =element_blank())+
  coord_cartesian(ylim = c(0, 15), xlim = c(-30, 30))+
  facet_wrap(~T_G) 

N809=ggplot(res2_N809) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x =element_blank()) +
  coord_cartesian(ylim = c(0, 15), xlim = c(-30, 30))+
  facet_wrap(~T_G)

N810=ggplot(res2_N810) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x =element_blank()) +
  coord_cartesian(ylim = c(0, 15), xlim = c(-30, 30))+
  facet_wrap(~T_G)

N6=ggplot(res2_N6) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x =element_blank()) +
  coord_cartesian(ylim = c(0, 15), xlim = c(-30, 30))+
  facet_wrap(~T_G)

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


library("data.table")
library(ggplot2)
library(parallel)
setwd("C:/Users/thula/Documents/researchR/htseq_illumina_genome/count_matrix/")

#Defining directory containg the expression data
directory <- "C:/Users/thula/Documents/researchR/htseq_illumina_transcriptome/count/"
sampleFiles <- grep(".txt",list.files(directory),value=TRUE)
#Defining samples
Treatment <- factor(c("N5.txt","N6.txt","N809.txt","N810.txt","WT.txt"))

sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = Treatment)

library( "DESeq" )
#Creating the count dataset (matrix)
cds = newCountDataSetFromHTSeqCount(sampleTable,directory = directory)

#Normalization of data
cds = estimateSizeFactors( cds )
sizeFactors( cds )
cds2 = cds[ ,c( "N809.txt","N810.txt","N5.txt","N6.txt", "WT.txt" ) ]

#Dispersion across samples
cds2 = estimateDispersions( cds2, method="blind", sharingMode="fit-only" )

#Differential expression analysis using nbiomtest
sample_N6 = c("atpol2a-1")
T_G=c("genome")
res2_N6_G = cbind(res2_N6_G, sample_N6,T_G)
sample_N809 = c("atpol2b-3")
res2_N809_G = cbind(res2_N809_G,sample_N809, T_G)
sample_N810 = c("atpol2b-1")
res2_N810_G = cbind(res2_N810_G,sample_N810, T_G)
sample_N5 = c("atpol2b-2")
res2_N5_G = cbind(res2_N5_G,sample_N5, T_G)
res2_N6_G = nbinomTest( cds2, "N6.txt", "WT.txt" )
write.table(res2, "table.txt", sep = ',', row.names = F, col.names = T, quote = F)
MA=plotMA(res2)
addmargins( table( res_sig = res2$padj < .05, res2_sig = res2$padj < .05 ) )

#Filter the DEGs by defining conditions as log2fc >= 2 or log2fc <=-2 and FDR <=0.05
sig=subset(res2_N6_G,res2_N6_G$log2FoldChange>=1|res2_N6_G$log2FoldChange<=-1)
sig=subset(sig,sig$padj<=0.05)
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

threshold_OE <- res2_N809_G$padj < 0.01 & (res2_N809_G$log2FoldChange>=1|res2_N809_G$log2FoldChange<=-1)
length(which(threshold_OE))
res2_N5_G$threshold <- threshold_OE 
res2_N809_G$threshold <- threshold_OE
res2_N810_G$threshold <- threshold_OE 
res2_N6_G$threshold <- threshold_OE 
## Volcano plot

ggplot(res2) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  ggtitle("N501413 vs WT") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  

par(mfrow=c(2,2),oma=c(2,2,0,0))

N5_G=ggplot(res2_N5_G) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x =element_blank())+
  coord_cartesian(ylim = c(0, 15), xlim = c(-30, 30))+
  facet_grid(sample_N5~T_G) 

N809_G=ggplot(res2_N809_G) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x =element_blank()) +
  coord_cartesian(ylim = c(0, 15), xlim = c(-30, 30))+
  facet_grid(sample_N809~T_G)

N810_G=ggplot(res2_N810_G) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x =element_blank()) +
  coord_cartesian(ylim = c(0, 15), xlim = c(-30, 30))+
  facet_grid(sample_N810~T_G)

N6_G=ggplot(res2_N6_G) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x =element_blank()) +
  coord_cartesian(ylim = c(0, 15), xlim = c(-30, 30))+
  facet_grid(sample_N6~T_G)

#plotmatrix(N5,N6,N809,N810)
#labeller=labeller(cyl=c())
c=cowplot::plot_grid(N6,N6_G,N810,N810_G,N5,N5_G,N809,N809_G,ncol = 2)
mtext("log2 fold change",side=1,line=0,outer=TRUE,cex=1.3)
mtext("-log10 adjusted p-value",side=2,line=0,outer=TRUE,cex=1.3,las=0)
?cowplot::plot_grid

annotate_figure(c,
                left = text_grob("-log10 P value", color = "black", rot = 90),
                bottom = text_grob("log2 fold change", color = "black"))

c$facet
library(ggpubr)          

