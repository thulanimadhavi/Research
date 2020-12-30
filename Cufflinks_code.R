
#####################Transcriptome aligned
#Filter the DEGs by defining conditions as log2fc >= 1 or log2fc <=-1 and FDR <=0.05
setwd("C:/Users/thula/Documents/researchR/Transcriptome/cuffdiff_illumina_transcriptome/output")
file=read.csv("gene_exp.diff.csv")
top_N5=subset(file,(file$sample_1== "N501413" & file$sample_2== "wildtype"))
top_N6=subset(file,(file$sample_1== "N676289" & file$sample_2== "wildtype"))
top_N809=subset(file,(file$sample_1== "N859809" & file$sample_2== "wildtype"))
top_N810=subset(file,(file$sample_1== "N859810" & file$sample_2== "wildtype"))

top_N6=subset(top_N6,(top_N6$log2fold_change>=1|top_N6$log2fold_change<=-1))
top_N6=subset(top_N6,top_N6$q_value<=0.05)
write.table(top_N6, "Sig_N6_T.txt", sep = ',',quote = FALSE)


sample_N6 = c("atpol2a-1")
T_G=c("Transcriptome")
top_N6 = cbind(top_N6, sample_N6,T_G)
sample_N809 = c("atpol2b-3")
top_N809 = cbind(top_N809,sample_N809, T_G)
sample_N810 = c("atpol2b-1")
top_N810 = cbind(top_N810,sample_N810, T_G)
sample_N5 = c("atpol2b-2")
top_N5 = cbind(top_N5,sample_N5, T_G)



install.packages("ggdraw")
library(reshape)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(cowplot)
library(ggpubr)

#mark the threshold
threshold_OE <- (top_N810$log2fold_change>=1|top_N810$log2fold_change<=-1) & top_N810$q_value<=0.01 
length(which(threshold_OE))
top_N5$threshold <- threshold_OE 
top_N809$threshold <- threshold_OE
top_N810$threshold <- threshold_OE 
top_N6$threshold <- threshold_OE 


## Volcano plots of samples
N5=ggplot(top_N5) +
  geom_point(aes(x=log2fold_change, y=-log10(q_value), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x =element_blank())+
  coord_cartesian(ylim = c(0, 3), xlim = c(-15, 15))+
  facet_grid(sample_N5~T_G) 

N809=ggplot(top_N809) +
  geom_point(aes(x=log2fold_change, y=-log10(q_value), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x =element_blank()) +
  coord_cartesian(ylim = c(0, 3), xlim = c(-15, 15))+
  facet_grid(sample_N809~T_G)

N810=ggplot(top_N810) +
  geom_point(aes(x=log2fold_change, y=-log10(q_value), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x =element_blank()) +
  coord_cartesian(ylim = c(0, 3), xlim = c(-15, 15))+
  facet_grid(sample_N810~T_G)

N6=ggplot(top_N6) +
  geom_point(aes(x=log2fold_change, y=-log10(q_value), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x =element_blank()) +
  coord_cartesian(ylim = c(0, 3), xlim = c(-15, 15))+
  facet_grid(sample_N6~T_G)
?facet_wrap


####################Genome aligned
#Filter the DEGs by defining conditions as log2fc >= 1 or log2fc <=-1 and FDR <=0.05
file=read.csv("gene_exp_g.diff.csv")
top_N5_G=subset(file,(file$sample_1== "N501413" & file$sample_2== "wildtype"))
top_N6_G=subset(file,(file$sample_1== "N676289" & file$sample_2== "wildtype"))
top_N809_G=subset(file,(file$sample_1== "N859809" & file$sample_2== "wildtype"))
top_N810_G=subset(file,(file$sample_1== "N859810" & file$sample_2== "wildtype"))

top_N810_G=subset(top_N810_G,(top_N810_G$log2fold_change>=1|top_N810_G$log2fold_change<=-1))
top_N810_G=subset(top_N810_G,top_N810_G$q_value<=0.05)
write.table(top_N810_G, "Sig_N810_G.txt", sep = ',',quote = FALSE)

sample_N6 = c("atpol2a-1")
G_T=c("Genome")
top_N6_G = cbind(top_N6_G, sample_N6,G_T)
sample_N809 = c("atpol2b-3")
top_N809_G = cbind(top_N809_G,sample_N809, G_T)
sample_N810 = c("atpol2b-1")
top_N810_G = cbind(top_N810_G,sample_N810, G_T)
sample_N5 = c("atpol2b-2")
top_N5_G = cbind(top_N5_G,sample_N5, G_T)





install.packages("ggdraw")
library(reshape)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(cowplot)
library(ggpubr)

#mark the threshold
threshold_OE <- (top_N810_G$log2fold_change>=1|top_N810_G$log2fold_change<=-1) & top_N810_G$q_value<=0.01 
length(which(threshold_OE))
top_N5_G$threshold <- threshold_OE 
top_N809_G$threshold <- threshold_OE
top_N810_G$threshold <- threshold_OE 
top_N6_G$threshold <- threshold_OE 

## Volcano plots of samples
N5_G=ggplot(top_N5_G) +
  geom_point(aes(x=log2fold_change, y=-log10(q_value), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x =element_blank())+
  coord_cartesian(ylim = c(0, 3), xlim = c(-15, 15))+
  facet_grid(sample_N5~G_T) 

N809_G=ggplot(top_N809_G) +
  geom_point(aes(x=log2fold_change, y=-log10(q_value), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x =element_blank()) +
  coord_cartesian(ylim = c(0, 3), xlim = c(-15, 15))+
  facet_grid(sample_N809~G_T)

N810_G=ggplot(top_N810_G) +
  geom_point(aes(x=log2fold_change, y=-log10(q_value), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x =element_blank()) +
  coord_cartesian(ylim = c(0, 3), xlim = c(-15, 15))+
  facet_grid(sample_N810~G_T)

N6_G=ggplot(top_N6_G) +
  geom_point(aes(x=log2fold_change, y=-log10(q_value), colour=threshold)) + 
  scale_colour_manual(values = c("black","red")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x =element_blank()) +
  coord_cartesian(ylim = c(0, 3), xlim = c(-15, 15))+
  facet_grid(sample_N6~G_T)

#plotmatrix(N5,N6,N809,N810)
#labeller=labeller(cyl=c())
c=cowplot::plot_grid(N6,N6_G,N810,N810_G,N5,N5_G,N809,N809_G,ncol = 2)

?cowplot::plot_grid

annotate_figure(c,
                left = text_grob("-log10 P value", color = "black", rot = 90),
                bottom = text_grob("log2 fold change", color = "black"))
