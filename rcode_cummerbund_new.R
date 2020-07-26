if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('clusterProfiler')
library(EnhancedVolcano)

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('Gviz')
BiocManager::install('dplyr')

#install.packages("RSQLite")
install.packages("Rtools")
setwd("C:/Users/thula/Desktop/researchR/Transcriptome/cuffdiff_illumina_transcriptome")

library("cummeRbund")
#library("RSQLite")
#Reading output of cuffdiff
cuffdata <- readCufflinks("output") 
cuffdata

#####Genes
#
dis = dispersionPlot(genes(cuffdata))
densityplot=csDensity(genes(cuffdata),replicates=F)
densityplot
#Volcanoplot of all samples
volcanoplot = csVolcanoMatrix(genes(cuffdata))
volcanoplot
#Volcanoplots of comparing ATPOL2b
?csVolcano
volcanoplot1_2 = csVolcano(genes(cuffdata),'N501413','N859809',features=FALSE,alpha=0.05, showSignificant=TRUE)
volcanoplot1_2
volcanoplot1_3 = csVolcano(genes(cuffdata),'N501413','N859810',features=FALSE,alpha=0.05, showSignificant=TRUE)
volcanoplot1_3
volcanoplot2_3 = csVolcano(genes(cuffdata),'N859809','N859810',features=FALSE,alpha=0.05, showSignificant=TRUE)
volcanoplot2_3
#Volcanoplots of ATPOL2b and wildtype
volcanoplot1_5 = csVolcano(genes(cuffdata),'N501413','wildtype',features=FALSE,alpha=0.05, showSignificant=TRUE)
volcanoplot1_5
volcanoplot2_5 = csVolcano(genes(cuffdata),'N859809','wildtype',features=FALSE,alpha=0.05, showSignificant=TRUE)
volcanoplot2_5
volcanoplot3_5 = csVolcano(genes(cuffdata),'N859810','wildtype',features=FALSE,alpha=0.05, showSignificant=TRUE)
volcanoplot3_5
#Volcanoplots comparing ATPOL2a and ATPOL2b types
volcanoplot1_4 = csVolcano(genes(cuffdata),'N501413','N676289',features=FALSE,alpha=0.05, showSignificant=TRUE)
volcanoplot1_4
volcanoplot2_4 = csVolcano(genes(cuffdata),'N859809','N676289',features=FALSE,alpha=0.05, showSignificant=TRUE)
volcanoplot2_4
volcanoplot3_4 = csVolcano(genes(cuffdata),'N859810','N676289',features=FALSE,alpha=0.05, showSignificant=TRUE)
volcanoplot3_4
#Volcano plot of AtPOL2a and wildtype
volcanoplot4_5 = csVolcano(genes(cuffdata),'N676289','wildtype',features=FALSE,alpha=0.05, showSignificant=TRUE)
volcanoplot4_5

#####

scatterplot=csScatterMatrix(genes(cuffdata))
scatterplot
scatterplot1_5=csScatter(genes(cuffdata1_5),'X1ATPOL2b','wild_type',features=FALSE,alpha=0.05, showSignificant=TRUE)
scatterplot1_5
?csScatter
boxplot=csBoxplot(isoforms(cuffdata))
boxplot
dendogram=csDendro(genes(cuffdata))
dendogram
?csDendro

MAplot=MAplot(genes(cuffdata),'N676289','wildtype')
csVolcanoMatrix(isoforms(cuffdata))
?csDendro

?csVolcano
scatterplot1_3=csScatterMatrix(genes(cuffdata1_3),genes(cuffdata1_2),replicates=F)
scatterplot1_3
boxplot1_3=csBoxplot(genes(cuffdata1_3))
boxplot1_3
dendogram1_3=csDendro(genes(cuffdata1_3))
csVolcanoMatrix(isoforms(cuffdata1_3))

#isoforms
csVolcanoMatrix(isoforms(cuffdata1_5))

########
#differentially and non differentially expressed genes between ATPOL2b samples
gene_diff_data1_2<- diffData(genes(cuffdata),'N501413','N859809') 
gene_diff_data1_2
gene_diff_data1_3<- diffData(genes(cuffdata),'N501413','N859810') 
gene_diff_data1_3
gene_diff_data2_3<- diffData(genes(cuffdata),'N859809','N859810') 
gene_diff_data2_3
#differentially and non differentially expressed genes between ATPOL2b samples and wildtype
gene_diff_data1_5<- diffData(genes(cuffdata),'N501413','wildtype') 
gene_diff_data1_5
gene_diff_data2_5<- diffData(genes(cuffdata),'N859809','wildtype') 
gene_diff_data2_5
gene_diff_data3_5<- diffData(genes(cuffdata),'N859810','wildtype') 
gene_diff_data3_5
write.table(gene_diff_data4_5, 'AtPOL2a_WT.txt', sep = ',', row.names = F, col.names = T, quote = F)
#differentially and non differentially expressed genes between ATPOL2b samples and ATPOL2a
gene_diff_data1_4<- diffData(genes(cuffdata),'N501413','N676289') 
gene_diff_data1_4
gene_diff_data2_4<- diffData(genes(cuffdata),'N859809','N676289') 
gene_diff_data2_4
gene_diff_data3_4<- diffData(genes(cuffdata),'N859810','N676289') 
gene_diff_data3_4
#differentially and non differentially expressed genes between ATPOL2a and wildtype
gene_diff_data4_5<- diffData(genes(cuffdata),'N676289','wildtype') 
gene_diff_data4_5

#####
#gene_id's of diffdata of differentially and non differentially expressed genes between ATPOL2b samples
gene_diff12_geneid=gene_diff_data1_2$gene_id
gene_diff13_geneid=gene_diff_data1_3$gene_id
gene_diff23_geneid=gene_diff_data2_3$gene_id
#gene_id's of diffdata of differentially and non differentially expressed genes between ATPOL2b samples and wildtype
gene_diff15_geneid=gene_diff_data1_5$gene_id
gene_diff25_geneid=gene_diff_data2_5$gene_id
gene_diff35_geneid=gene_diff_data3_5$gene_id
#gene_id's of diffdata of differentially and non differentially expressed genes between ATPOL2b samples and ATPOL2a
gene_diff14_geneid=gene_diff_data1_4$gene_id
gene_diff24_geneid=gene_diff_data2_4$gene_id
gene_diff34_geneid=gene_diff_data3_4$gene_id
#gene_id's of diffdata of differentially and non differentially expressed genes between ATPOL2a and wildtype
gene_diff45_geneid=gene_diff_data4_5$gene_id


#########
##Upregulated
#Between ATPOL2b types
#1_2
sig_gene_dataup1_2 <- subset(gene_diff_data1_2,(significant=='yes') & (log2_fold_change>=2)) 
sig_gene_dataup1_2 
nrow(sig_gene_dataup1_2) 
upregsig_genid1_2=sig_gene_dataup1_2$gene_id
#1_3
sig_gene_dataup1_3 <- subset(gene_diff_data1_3,(significant=='yes') && (log2_fold_change>=2)) 
sig_gene_dataup1_3 
nrow(sig_gene_dataup1_3) 
upregsig_genid1_3=sig_gene_dataup1_3$gene_id
#2_3
sig_gene_dataup2_3 <- subset(gene_diff_data2_3,(significant=='yes') & (log2_fold_change>=2)) 
sig_gene_dataup2_3 
nrow(sig_gene_dataup2_3) 
upregsig_genid2_3=sig_gene_dataup2_3$gene_id
#Between ATPOL2b types and wildtype
#1_5
sig_gene_dataup1_5 <- subset(gene_diff_data1_5,(significant=='yes') & (log2_fold_change>=2)) 
sig_gene_dataup1_5 
nrow(sig_gene_dataup1_5) 
upregsig_genid1_5=sig_gene_dataup1_5$gene_id
#2_5
sig_gene_dataup2_5 <- subset(gene_diff_data2_5,(significant=='yes') & (log2_fold_change>=2)) 
sig_gene_dataup2_5 
nrow(sig_gene_dataup2_5) 
upregsig_genid2_5=sig_gene_dataup2_5$gene_id
#3_5
sig_gene_dataup3_5 <- subset(gene_diff_data3_5,(significant=='yes') & (log2_fold_change>=2)) 
sig_gene_dataup3_5 
nrow(sig_gene_dataup3_5) 
upregsig_genid3_5=sig_gene_dataup3_5$gene_id
#Between ATPOL2b types and AtPOL2a
#1_4
sig_gene_dataup1_4 <- subset(gene_diff_data1_4,(significant=='yes') & (log2_fold_change>=2)) 
sig_gene_dataup1_4 
nrow(sig_gene_dataup1_4) 
upregsig_genid1_4=sig_gene_dataup1_4$gene_id
#2_4
sig_gene_dataup2_4 <- subset(gene_diff_data2_4,(significant=='yes') & (log2_fold_change>=2)) 
sig_gene_dataup2_4 
nrow(sig_gene_dataup2_4) 
upregsig_genid2_4=sig_gene_dataup2_4$gene_id
#3_4
sig_gene_dataup3_4 <- subset(gene_diff_data3_4,(significant=='yes') & (log2_fold_change>=2)) 
sig_gene_dataup3_4 
nrow(sig_gene_dataup3_4) 
upregsig_genid3_4=sig_gene_dataup3_4$gene_id
#Between ATPOL2a and wildtype
sig_gene_dataup4_5 <- subset(gene_diff_data4_5,(significant=='yes') & (log2_fold_change>=2)) 
sig_gene_dataup4_5 
nrow(sig_gene_dataup4_5) 
upregsig_genid4_5=sig_gene_dataup4_5$gene_id

####
##Downregulated
#Between ATPOL2b types
#1_2
sig_gene_datadown1_2 <- subset(gene_diff_data1_2,(significant=='yes') & (log2_fold_change<=-2)) 
sig_gene_datadown1_2 
nrow(sig_gene_datadown1_2) 
downregsig_genid1_2=sig_gene_datadown1_2$gene_id
#1_3
sig_gene_datadown1_3 <- subset(gene_diff_data1_3,(significant=='yes') & (log2_fold_change<=-2)) 
sig_gene_datadown1_3 
nrow(sig_gene_datadown1_3) 
downregsig_genid1_3=sig_gene_datadown1_3$gene_id
#2_3
sig_gene_datadown2_3 <- subset(gene_diff_data2_3,(significant=='yes') & (log2_fold_change<=-2)) 
sig_gene_datadown2_3 
nrow(sig_gene_datadown2_3) 
downregsig_genid2_3=sig_gene_datadown2_3$gene_id
#Between ATPOL2b types and wildtype
#1_5
sig_gene_datadown1_5 <- subset(gene_diff_data1_5,(significant=='yes') & (log2_fold_change<=-2)) 
sig_gene_datadown1_5 
nrow(sig_gene_datadown1_5) 
downregsig_genid1_5=sig_gene_datadown1_5$gene_id
#2_5
sig_gene_datadown2_5 <- subset(gene_diff_data2_5,(significant=='yes') & (log2_fold_change<=-2)) 
sig_gene_datadown2_5 
nrow(sig_gene_datadown2_5) 
downregsig_genid2_5=sig_gene_datadown2_5$gene_id
#3_5
sig_gene_datadown3_5 <- subset(gene_diff_data3_5,(significant=='yes') & (log2_fold_change<=-2)) 
sig_gene_datadown3_5 
nrow(sig_gene_datadown3_5) 
downregsig_genid3_5=sig_gene_datadown3_5$gene_id
#Between ATPOL2b types and AtPOL2a
#1_4
sig_gene_datadown1_4 <- subset(gene_diff_data1_4,(significant=='yes') & (log2_fold_change<=-2)) 
sig_gene_datadown1_4 
nrow(sig_gene_datadown1_4) 
downregsig_genid1_4=sig_gene_datadown1_4$gene_id
#2_4
sig_gene_datadown2_4 <- subset(gene_diff_data2_4,(significant=='yes') & (log2_fold_change<=-2)) 
sig_gene_datadown2_4 
nrow(sig_gene_datadown2_4) 
downregsig_genid2_4=sig_gene_datadown2_4$gene_id
#3_4
sig_gene_datadown3_4 <- subset(gene_diff_data3_4,(significant=='yes') & (log2_fold_change<=-2)) 
sig_gene_datadown3_4 
nrow(sig_gene_datadown3_4) 
downregsig_genid3_4=sig_gene_datadown3_4$gene_id
#Between ATPOL2a and wildtype
#4_5
sig_gene_datadown4_5 <- subset(gene_diff_data4_5,(significant=='yes') & (log2_fold_change<=-2)) 
sig_gene_datadown4_5 
nrow(sig_gene_datadown4_5) 
downregsig_genid4_5=sig_gene_datadown4_5$gene_id

###########
#gene_name

gene.features<-annotation(genes(cuffdata))
gene.features
gene.features$gene_short_name
gene.features$gene_id
########
gene_id=gene.features$gene_id
gene_id
?genes
###########
#downregulated significant genes
downgenname1_5=numeric(0)
for (i in 1:length(downregsig_genid1_5)){
  for (k in 1:length(gene_id)){
    if (downregsig_genid1_5[i]==gene_id[k])
      gen_name=print(gene.features$gene_short_name[k])
    downgenname1_5[i]=c(gen_name)
  }
}
downgenname1_5

#upregulated significant genes
upgenname1_5=numeric(0)
for (i in 1:length(upregsig_genid1_5)){
  for (k in 1:length(gene_id)){
    if (upregsig_genid1_5[i]==gene_id[k])
      gen_name=print(gene.features$gene_short_name[k])
    upgenname1_5[i]=c(gen_name)
  }
}
upgenname1_5

#to get differentially expressed data(either upregulated or downregulated)
sig_gene_data4_5 <- subset(gene_diff_data4_5,(significant=='yes')) 
sig_gene_data4_5 
nrow(sig_gene_data4_5)
sig_gengen_id4_5=sig_gene_data4_5$gene_id
write.table(sig_gene_data1_2, 'diff_genes_1_2.txt', sep = '/t', row.names = F, col.names = T, quote = F)

#########
#differentially expressed genes(both up and down)
#AtPOL2b's
sig_gene_data1_2 = rbind(sig_gene_datadown1_2,sig_gene_dataup1_2)
sig_gene_data1_2
sig_gene_data1_3=rbind(sig_gene_datadown1_3,sig_gene_dataup1_3)
sig_gene_data1_3
sig_gene_data2_3=rbind(sig_gene_datadown2_3,sig_gene_dataup2_3)
sig_gene_data2_3
sig_gene_data_AtPOL2bs=rbind(sig_gene_data1_2,sig_gene_data1_3,sig_gene_data2_3)
sig_gene_data_AtPOL2bs

#AtPOL2b_WT
sig_gene_data1_5=rbind(sig_gene_datadown1_5,sig_gene_dataup1_5)
sig_gene_data1_5
sig_gene_data2_5=rbind(sig_gene_datadown2_5,sig_gene_dataup2_5)
sig_gene_data2_5
sig_gene_data3_5=rbind(sig_gene_datadown3_5,sig_gene_dataup3_5)
sig_gene_data3_5
sig_gene_data_AtPOL2b_WT=rbind(sig_gene_data1_5,sig_gene_data2_5,sig_gene_data3_5)

#AtPOL2b_AtPOL2a
sig_gene_data1_4=rbind(sig_gene_datadown1_4,sig_gene_dataup1_4)
sig_gene_data1_4
sig_gene_data2_4=rbind(sig_gene_datadown2_4,sig_gene_dataup2_4)
sig_gene_data2_4
sig_gene_data3_4=rbind(sig_gene_datadown3_4,sig_gene_dataup3_4)
sig_gene_data3_4
sig_gene_data_AtPOL2b_AtPOL2a=rbind(sig_gene_data1_4,sig_gene_data2_4,sig_gene_data3_4)

#AtPOL2a_WT
sig_gene_data_AtPOL2a_WT=rbind(sig_gene_datadown4_5,sig_gene_dataup4_5)
sig_gene_data_AtPOL2a_WT

sig_gene_data1_5 <- subset(gene_diff_data1_5,(significant=='yes')) 
sig_gene_data1_5 
nrow(sig_gene_data1_5) 
sig_genid1_5=sig_gene_dataup1_5$gene_id
gene_diff15_geneid=gene_diff_data1_5$gene_id
#to obtain gene names which are significantly differentially expressed
vecAtPOL2b1_WT=numeric(0)
for (i in 1:length(gene_diff_data5_1_gene_id)){
  for (k in 1:length(gene_id)){
    if (gene_diff_data4_5_gene_id[i]==gene_id[k])
      gen_name=print(gene.features$gene_short_name[k])
    vecAtPOL2b4_WT[i]=c(gen_name)
  }
}
vecAtPOL2b4_WT
write.table(vecAtPOL2b4_WT, 'AtPOL2b4genid.txt', sep = ',', row.names = F, col.names = T, quote = F)

######
#FPKM value expression of genes across samples
a=read.csv("flowering_related.csv")

GR1
RAD51
BRCA1
MRE11
KU70
MSH2
MSH6
BARD1
PARP1
PARP2
AtATM
AtATR
RNR1
RNR2
SRP2
SRP3
RAD17
SYN2
DNA ligase iv
RAD50

RUB1
UBP7

gene_int <- getGene(cuffdata, 'SSEP3') 
SEP3=expressionBarplot(gene_int)
cowplot::plot_grid(AG,AGL1,SEP1,SEP2,SEP3,AGL24,ncol = 3,labels = LETTERS[1:6])

#dispersion plot
disp<-dispersionPlot(genes(cuffdata))
disp

#scv plot
genes.scv<-fpkmSCVPlot(genes(cuffdata))
genes.scv
runInfo(cuffdata)

###############
#genes FPKM values
gene.fpkm<-fpkm(genes(cuffdata))
gene.fpkm
head(gene.fpkm)
#gene_count
gene.counts<-count(genes(cuffdata))
gene.counts
head(gene.counts)
#isoform FPKM
isoform.fpkm<-fpkm(isoforms(cuffdata))
head(isoform.fpkm)
#differential expression
gene.diff<-diffData(genes(cuffdata))
gene.diff
head(gene.diff)
#sample_names
sample.names<-samples(genes(cuffdata))
head(sample.names)
#genefeature_names
gene.featurenames<-featureNames(genes(cuffdata))
gene.featurenames
head(gene.featurenames)

#Gene_FPKM matrix
gene.matrix<-fpkmMatrix(genes(cuffdata))
head(gene.matrix)
#gene_count matrix
gene.count.matrix<-countMatrix(genes(cuffdata))
head(gene.count.matrix)

#########
sig_gene_data1_2 <- subset(gene_diff_data1_2,(significant=='yes')) 
sig_gene_data1_2 
myGeneIds=sig_gene_data1_2$gene_id
Fi=head(myGeneIds)
myGenes<-getGenes(cuffdata,myGeneIds)
myGenes
fpkm(myGenes)
head(fpkm(myGenes))
head(fpkm(isoforms(myGenes)))
h<-csHeatmap(myGenes,cluster='both')
h
?csHeatmap
b<-expressionBarplot(myGenes)
b
s<-csScatter(myGenes,"N501413","N859809",smooth=T)
s
v<-csVolcano(myGenes,"N501413","N859809")
v
gl<-expressionPlot(myGenes)
gl
igb<-expressionBarplot(isoforms(myGenes))
igb


#Venn diagrams - clustering
install.packages('VennDiagram')
library('VennDiagram')
install.packages('gplots')
require(gplots)
require(VennDiagram)
###differentially expressed
d1=venn.diagram(list(AtPOL2b1_AtPOL2b2=vec1_2,AtPOL2b1_AtPOL2b3=vec1_3,AtPOL2b2_AtPOL2b3=vec2_3),filename=NULL, fill=rainbow(3))
d2=venn.diagram(list(AtPOL2b1_wildtype=vec1_5,AtPOL2b2_wildtype=vec2_5,AtPOL2b3_wildtype=vec3_5),filename=NULL, fill=rainbow(3))
d3=venn.diagram(list(AtPOL2b1_AtPOL2a=vec1_4,AtPOL2b2_AtPOL2a=vec2_4,AtPOL2b3_AtPOL2a=vec3_4),filename=NULL, fill=rainbow(3))
grid.newpage()
grid.draw(d3)

###downregulated
d4=venn.diagram(list(AtPOL2b1_wildtype=downgenname1_5,AtPOL2b2_wildtype=downgenname2_5,AtPOL2b3_wildtype=downgenname3_5),filename=NULL, fill=rainbow(3))
d5=venn.diagram(list(AtPOL2b1_AtPOL2b2=downgenname1_2,AtPOL2b1_AtPOL2b3=downgenname1_3,AtPOL2b2_AtPOL2b3=downgenname2_3),filename=NULL, fill=rainbow(3))
d6=venn.diagram(list(AtPOL2b1_AtPOL2a=downgenname1_4,AtPOL2b2_AtPOL2a=downgenname2_4,AtPOL2b3_AtPOL2a=downgenname3_4),filename=NULL, fill=rainbow(3))
grid.newpage()
grid.draw(d6)
#upregulated
d7=venn.diagram(list(AtPOL2b1_wildtype=upgenname1_5,AtPOL2b2_wildtype=upgenname2_5,AtPOL2b3_wildtype=upgenname3_5),filename=NULL, fill=rainbow(3))
d8=venn.diagram(list(AtPOL2b1_AtPOL2b2=upgenname1_2,AtPOL2b1_AtPOL2b3=upgenname1_3,AtPOL2b2_AtPOL2b3=upgenname2_3),filename=NULL, fill=rainbow(3))
d9=venn.diagram(list(AtPOL2b1_AtPOL2a=upgenname1_4,AtPOL2b2_AtPOL2a=upgenname2_4,AtPOL2b3_AtPOL2a=upgenname3_4),filename=NULL, fill=rainbow(3))
grid.newpage()
grid.draw(d9)

#Biological Id TranslatoR
vecx=c(upgenname1_4,upgenname2_4,upgenname3_4,downgenname1_4,downgenname2_4,downgenname3_4)
eg = bitr(vecx, fromType="SYMBOL", toType="TAIR", OrgDb="org.At.tair.db")
eg

install.packages('KEGGREST')
source("http://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")
install.packages("doseplot")
BiocManager::install("pathfindR")
BiocManager::install("doseplot")
BiocManager::install("VennDiagram")
BiocManager::install("KEGGREST")
BiocManager::install("cufflinks")
BiocManager::install("BiocGenerics")
library(pathview)
library(pathfindR)
library(KEGGREST)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(doseplot)
install.packages('KEGGREST')
require(doseplot)

################Functional enrichment-KEGG
#TAIR ids of diff exp genes
#AtPOL2b's
AtPOL2b_diff=read.csv("diff_genes_AtPOL2bs_names_genelist.csv")
AtPOL2b_diff=as.vector(AtPOL2b_diff$Gene_ID)
AtPOL2b_diff
#AtPOL2b's and WT
AtPOL2b_WT_diff=read.csv("diff_genes_AtPOL2b_WT_names.csv")
AtPOL2b_WT_diff=as.vector(AtPOL2b_WT_diff$Gene_ID)
AtPOL2b_WT_diff
#AtPOL2b's and AtPOL2a
AtPOL2b_AtPOL2a_diff=read.csv("diff_genes_AtPOL2b_AtPOL2a_names.csv")
AtPOL2b_AtPOL2a_diff=as.vector(AtPOL2b_AtPOL2a_diff$Gene_ID)
AtPOL2b_AtPOL2a_diff
#AtPOL2a and WT
AtPOL2a_WT_diff=read.csv("diff_genes_AtPOL2a_WT_names.csv")
AtPOL2a_WT_diff=as.vector(AtPOL2a_WT_diff$Gene_ID)
AtPOL2a_WT_diff

#pathview file of all significant genes
pv.out=pathview(gene.data =AtPOL2b_diff , pathway.id = "ath00950", gene.idtype = "KEGG", species = "ath",out.suffix = "kegg.get.all", kegg.native = T)

#vector of gene names of all significant genes
#vector of log fold values of all significant genes
logvec=c(sig_gene_data1_2$log2_fold_change,sig_gene_data1_3$log2_fold_change,sig_gene_data1_4$log2_fold_change,sig_gene_data1_5$log2_fold_change,sig_gene_data2_3$log2_fold_change,sig_gene_data2_4$log2_fold_change,sig_gene_data2_5$log2_fold_change,sig_gene_data3_4$log2_fold_change,sig_gene_data3_5$log2_fold_change,sig_gene_data4_5$log2_fold_change)

#####KEGG
#AtPOL2b
AtPOL2b_KEGG=enrichKEGG(gene         = diff_genes_AtPOL2bs_names_genelist[,3],
                       organism     = 'ath',
                       pvalueCutoff=1, qvalueCutoff=1)
#AtPOL2b_WT
AtPOL2b_WT_KEGG=enrichKEGG(gene         = diff_genes_AtPOL2b_WT_names_genelist[,3],
                           organism     = 'ath',
                           pvalueCutoff=1, qvalueCutoff=1)
#AtPOL2bs_AtPOL2a
AtPOL2b_AtPOL2a_KEGG=enrichKEGG(gene         = diff_genes_AtPOL2b_AtPOL2a_names_genelist[,3],
                                organism     = 'ath',
                                pvalueCutoff=1, qvalueCutoff=1)
#AtPOL2a_WT
AtPOL2a_WT_KEGG=enrichKEGG(gene         = diff_genes_AtPOL2a_WT_names_genelist[,3],
                                organism     = 'ath',
                                pvalueCutoff=1, qvalueCutoff=1)

######used one
geneList=log2FC
names(geneList)=list
geneList
list=c('AT1G02205','AT1G02310','AT1G02340','AT1G04800','AT1G06360','AT1G08310','AT1G09750','AT1G11545','AT1G11850','AT1G12570','AT1G13609','AT1G13650','AT1G17710','AT1G18620','AT1G20190','AT1G21100','AT1G22500','AT1G23110','AT1G26390','AT1G27020','AT1G30700','AT1G31173','AT1G34060','AT1G47395','AT1G49500','AT1G52000','AT1G52040','AT1G52400','AT1G53520','AT1G54010','AT1G54575','AT1G55330','AT1G58225','AT1G58270','AT1G60590','AT1G61800','AT1G62510','AT1G65450','AT1G65690','AT1G66100','AT1G66390','AT1G67810','AT1G68740','AT1G69530','AT1G70640','AT1G70800','AT1G73010','AT1G73260','AT1G73325','AT1G73330','AT1G73480','AT1G74670','AT1G80130','AT2G02010','AT2G02990','AT2G04135','AT2G04460','AT2G05790','AT2G06850','AT2G10940','AT2G14247','AT2G14560','AT2G14610','AT2G20670','AT2G20750','AT2G23130','AT2G26010','AT2G26020','AT2G27120','AT2G28190','AT2G29350','AT2G29460','AT2G30750','AT2G30760','AT2G30766','AT2G30770','AT2G32830','AT2G34210','AT2G34430','AT2G34620','AT2G37130','AT2G38940','AT2G39030','AT2G39800','AT2G40610','AT2G41240','AT2G42380','AT2G42870','AT2G43510','AT2G45130','AT2G45220','AT2G45570','AT2G46880','AT3G01420','AT3G01500','AT3G02030','AT3G02040','AT3G03260','AT3G03820','AT3G04290','AT3G04530','AT3G04720','AT3G05730','AT3G06880','AT3G07350','AT3G08030','AT3G08860','AT3G09922','AT3G13750','AT3G14210','AT3G15450','AT3G15950','AT3G16420','AT3G16460','AT3G16670','AT3G17520','AT3G17790','AT3G18250','AT3G21500','AT3G21720','AT3G22120','AT3G22121','AT3G22142','AT3G22600','AT3G23550','AT3G24460','AT3G25240','AT3G25760','AT3G27690','AT3G28220','AT3G29590','AT3G43110','AT3G44300','AT3G44510','AT3G46900','AT3G47340','AT3G47420','AT3G48850','AT3G55646','AT3G55970','AT3G56970','AT3G56980','AT3G58120','AT3G60140','AT4G01390','AT4G01460','AT4G11650','AT4G12480','AT4G12490','AT4G12730','AT4G13575','AT4G14130','AT4G15210','AT4G15610','AT4G16260','AT4G16540','AT4G18670','AT4G18970','AT4G19810','AT4G22010','AT4G22070','AT4G22470','AT4G22880','AT4G23140','AT4G23600','AT4G23700','AT4G23820','AT4G28780','AT4G29020','AT4G30290','AT4G31290','AT4G33790','AT4G34250','AT4G37990','AT4G38400','AT5G01015','AT5G03210','AT5G03545','AT5G04150','AT5G04950','AT5G05250','AT5G07100','AT5G09440','AT5G13170','AT5G15500','AT5G15780','AT5G17220','AT5G19120','AT5G20150','AT5G20250','AT5G20630','AT5G20790','AT5G22580','AT5G23020','AT5G24770','AT5G24780','AT5G25190','AT5G25460','AT5G25980','AT5G26000','AT5G38780','AT5G40890','AT5G42800','AT5G43350','AT5G43580','AT5G44380','AT5G44420','AT5G44430','AT5G48490','AT5G48900','AT5G49100','AT5G49360','AT5G50740','AT5G52740','AT5G54060','AT5G54190','AT5G56870','AT5G57760','AT5G61160','AT5G62920','AT5G63130','AT5G65080')
label_name=c('CR88','LHCB2.2','PR1','PDF1.2b','CYP71A13','LHB1B1','AT3G06530','PYK10','Cpn60beta2','AHP4','EMB2729','DXPS1','CYP82G1','TSO2','LHCB2.3','AT3G43670','ASN1','AT3G49160','CML41','F3H','NAP57','PSBQ-2','AT4G14090','BAM5','AT4G19130','PARP2','ATPMEPCRB','LDOX','AMY1','FSD1','FIB2','IAA29','KCS16','KCS17','ELI3-2','CHIL','TT7','SAUR22','TT4','GSTF12','TGG2','TGG1','DFR','PDF1.2','PDF1.2c','AT5G48900','UF3GT','LHCB3','PSAN','CER1','PSAD-2','CAT3','AT1G25230','PSAF','HEMA1','TPS04','4CL3','AT1G67120','AT2G14247','AT4G01985','AT1G03495','AT1G13609','AT1G35300','AT1G36460','AT1G39110','AT1G42050','AT1G44060','AT2G14230','AT2G06720','AT2G04770','AT2G06340','AT2G06180','AT2G13000','AT2G06800','AT2G14245','AT2G11780','AT3G30396','AT3G42716','AT3G42719','AT3G42720','AT3G43304','AT4G08091','AT4G03790','AT4G06720','AT4G08050','AT4G07850','AT5G28165','AT5G28464','AT5G32516','AT5G32624','AT5G33389','AT5G45082','AT1G10657','AT3G22142','AT1G53887','AT2G30766','AT1G51402','RNS1','AT2G10940','AT2G18193','EXPB1','MEE3','scpl10','SCA3','PDF1.3','SIP4','AT2G30760','RLP24','AT2G34170','OFP15','AT2G36885','ChlAKR','AT2G38240','FER4','EXPA8','BHLH100','AT2G41800','BZIP34','PAR1','CYCP4;1','AT3G06070','AT3G10020','BSMT1','EMB3120','AT3G15450','AT3G16670','BTS','ELIP1','AT3G25290','AT3G27630','AT5MAT','CYP94B3','CAX3','FLN1','AT3G55240','JRG21','bHLH38','bHLH39','BZIP61','PRN','TT8','AT4G10290','EARLI1','ELIP2','AT4G22960','AT4G23680','AT4G27450','XTH19','SEN1','AT4G36010','AT4G39795','BHLH101','AT5G05250','PDE340','SAG29','ENH1','FER1','AT5G02160','AT5G22580','GMI1','RLP53','RVE2','BXL1','CP1','NEET','AT5G52050','AT5G52390','ATM2','LTP4','LTP3','AT5G60250','MAF5','DAR3','AT1G07500','ELP','AT1G15125','AT1G21140','PIP1C','AT1G03940','AT1G04800','AT1G04540','COR413IM1','LOL1','AT1G47395','AT1G47400','AT1G51400','NAC019','PUB19','MYB90','BBX27','FLN2','AT1G70640','GASA6','AT1G74770','AT1G76800','AT1G77960','J8','AT1G10682','AT1G67105','AT2G15555','AT2G30362','AT3G24982','AT3G51238','AT4G16370')
log2FC=c(3.19714861,-3.902692582,2.927937264,-3.655930866,-5.12664771,4.388346154,5.021176004,-5.245552268,3.731272659,3.588589898,4.53536481,-6.280485094,3.379733981,4.2871801,5.802007872,3.344582908,3.313216677,3.350363473,3.369362995,-10.13847252,-8.716128778,-10.03093095,-7.18089562,-9.751843176,-6.082103178,-7.835157031,2.960528011,3.014880243,-3.467613093,6.317528341,3.29631221,-3.593324967,-6.315684637,-3.533312215,-4.4148577,4.676793141,-3.065954386,3.054510302,-3.412453393,-4.34214583,4.501604508,-4.836964822,4.912749868,3.767907615,2.948864909,-5.599604168,-3.045060057,-9.032309287,3.371378442,-9.0524489,-9.072311235,-10.22973431,-9.167717089,4.998408012,-5.980183205,-9.167717089,-6.571075911,-11.50965567,-4.043941192,5.206822461,-8.765765521,-3.929933841,4.300273722,3.954878235,-3.212854052,-3.057459587,7.005839981,5.842709187,-3.698724999,-3.43521752,-4.918604334,-4.969968387,-3.257537539,-5.580249454,3.583044806,5.681101073,4.435093552,4.067054312,-4.061789354,-5.121583431,4.73295223,3.417930225,-6.841271724,-7.841318759,3.877728317,5.758115394,5.343266371,3.09681952,-3.482948223,-3.383447512,3.646620733,-7.851705487,-3.002576259,-4.725083637,3.443930338,5.153142828,9.199688273,-3.687671937,-4.405466775,-4.603681676,5.291808899,-3.085089105,-3.936341627,-7.300003258,-3.147693273,-5.081726418,-5.77601366,6.378058879,-5.495689624,-9.072311235,-10.23855083,-9.111233898,-10.05108968,-8.611418102,4.032836477,3.905927029,-7.529334892,5.111722395,-4.042819657,-3.315969831,-3.299820245,-3.861269539,-3.013594551,8.047330319,-3.716696303,-5.737345768,-5.941727177,-3.057580768,3.553707343,-5.371319252,-5.787085862,-4.274005928,-4.596418997,-8.716128778,2.9188106,-9.239744503,-4.449558536,-8.664723218,-9.373833459,-5.52605046,-9.751843176,-3.068177405,-3.967677164,-5.933676413,-3.866519321,-3.154216879,-4.651717445,-6.250713305,-5.933870819,-6.243095724,-3.158896264,4.212402041,-3.295955219,3.371564047,-4.415863254,5.464953315,3.893080044,-7.254078907,4.582640776,-3.25912646,-5.498560278,5.30236931,3.411837594,3.014507525,-3.590761755,-3.851942207,-4.420789273,-4.203762197,-3.903760504,-5.856211004,-3.961318359,3.000504082,-4.397285312,5.174049855,4.017269207,-3.375474527,3.792260923,4.887129321,-4.300219849,-9.48173627,-8.789957824,-8.860191936,-7.885916362,-9.0524489,-3.529366758,-5.337611311,5.212037085,7.452238912,-9.239744503,3.416580802,3.613495148,-3.140178283,8.513552707,-5.493560744,-6.321739345,-5.671079709,2.990540817,3.520807801,-10.39779386,-5.759364861,-6.573315391,3.651870849,6.057421897,4.065948949)
log2FC=c(10.90284899,2.46407044,2.658567117,7.499412573,9.906888527,7.798552235,8.273643226,1.992889385,-2.560715303,1.614932051,9.308088615,6.210066748,6.385412546,1.615890145,7.286791167,6.529657566,5.988452779,7.424228356,5.652482203,5.705987496,7.09143615,5.671807686,6.751699274,1.79879344,0.412490489,6.629276387,6.381354208,5.809854973,6.813543066,1.150436753,5.100831366,3.760116963,5.310930833,7.476122937,0.85091356,4.65329641,5.384500173,6.602005354,4.689357174,7.828943725,4.711981562,-0.617178464,3.230595695,4.901851431,2.505692179,-2.125970533,5.150039237,1.022991178,4.682795224,5.14874891,5.918845787,-0.610597225,4.685835356,5.070072149,5.224288002,0.875525967,7.435622624,4.878478826,4.523952063,4.259140326,0.909172287,5.720313642,4.319128814,5.684562362,4.167693372,5.61517792,4.381427134,0.913970808,4.523650684,1.600978365,0.832459272,0.770466133,1.58351529,3.685945875,1.134937789,3.899839922,4.047304093,-1.409447241,0.815538011,1.755291696,3.678700305,0.775047425,5.266171703,1.480071556,3.816455704,3.114716858,1.437098881,0.813288258,1.390600603,3.733785336,-3.379188981,5.124762223,0.721651086,3.926882374,3.481816563,3.010293531,3.658895611,3.299795458,3.279265375,0.641269872,1.411586507,3.056727455,2.819861345,3.323302,1.352083581,3.019718217,3.271274324,3.473028338,-1.203098682,1.319474344,3.125766165,2.832071042,3.179911959,1.323870351,0.059834367,1.496137331,1.164139708,1.226064169,2.789312559,2.730957243,2.71015607,1.228934552,1.524357911,2.932471797,2.611857709,0.356350909,-0.991679159,2.503262821,3.282689848,-1.033464967,3.18697608,1.213643122,0.193092063,2.342131677,2.703532772,1.559795631,0.120222038,2.480990906,2.285344709,-1.111639959,2.598561023,2.277345301,2.193406793,2.89685604,1.125708556,1.140671807,0.145897954,2.341384763,2.733238307,0.04242219,1.063194357,-0.911153631,2.08091362,-0.845389384,-0.744524587,0.984630112,-0.743992144,0.165453514,-0.13015942,-0.693808837,-0.772250487,0.761033613,0.226673206,-0.578047984,3.099691383,0.006873164,0.188046926,1.988689989,-0.570375307,0.229078709,0.686613486,-0.530556329,-0.58181696,-0.4637758,-0.184997445,-0.161360732,-1.139299599,0.113803808,1.071395311,0.132195789,0.252587195,0.133355828,-0.463651626,0.358918885,-0.341321306,-0.149996888,-0.427783387,-0.907417268,0.103845054,0.085338086,0.403458107,0.810375623,-0.291743753,0.120582278,-0.287949511,-0.305103773,0.707818785,-0.211754949,0.06498241,-0.20020379,-0.228468119,-0.206182315,-0.01705561,-0.826300456,-0.005507542,0.190558175,-0.334554508,-0.199332007,-0.110859307,0.006646426,-0.135551108,-0.003651032,0.074017224,0.005675777,-0.034446275,0.001508066,-0.089965155,-0.038721279,0.0118155,-0.006610145,0.325976775,0.028184933)
#genes <- names(geneList)[abs(geneList) > 2]
kk <- enrichKEGG(gene         = labelHM,
                 organism     = 'ath',
                 pvalueCutoff=1, qvalueCutoff=1)
head(kk)
kk$ID
browseKEGG(kk, 'ath00750')
head(summary(kk))
b2=barplot(kk, showCategory=100)
clusterProfiler::dotplot(kk,showCategory=200)
cnet1=cnetplot(kk, foldChange=geneList, showCategory=100,fixed=FALSE)
cnet1
?cnetplot

cowplot::plot_grid(cnet1,cnet2,ncol = 2,labels = c("(A)","(B)"))
?cowplot::plot_grid
p2 <- cnetplot(kk, node_label="gene")
heatplot(kk,showCategory=100)
emapplot(kk,showCategory=100)
upsetplot(kk)
ridgeplot(kk)
gseaplot(kk)
?barplot
####
x2 <- setReadable(kk, OrgDb ='org.At.tair.db', 'TAIR')
head(summary(x2))


geneList
gene.idtype.list
library(topGO)
library(org.At.tair.db)
library(Rgraphviz)
BiocManager::install(c("topGO", "KEGGREST", "org.At.tair.db", "Rgraphviz"))

#Biological Id TranslatoR
??bitr
labelHM=c('AT1G01620','AT1G02205','AT1G03130','AT1G03495','AT1G03940','AT1G04540','AT1G04800','AT1G07500','AT1G10657','AT1G10682','AT1G12090','AT1G13609','AT1G15125','AT1G20620','AT1G21140','AT1G25230','AT1G29395','AT1G31330','AT1G32540','AT1G35300','AT1G36460','AT1G39110','AT1G42050','AT1G44060','AT1G47395','AT1G47400','AT1G51400','AT1G51402','AT1G52890','AT1G53887','AT1G58290','AT1G60190','AT1G61120','AT1G65060','AT1G66390','AT1G67105','AT1G67120','AT1G68190','AT1G69200','AT1G70640','AT1G74670','AT1G74770','AT1G76800','AT1G77960','AT1G80920','AT2G02990','AT2G04030','AT2G04770','AT2G05070','AT2G06180','AT2G06340','AT2G06720','AT2G06800','AT2G10940','AT2G11780','AT2G13000','AT2G14230','AT2G14245','AT2G14247','AT2G14610','AT2G15555','AT2G18193','AT2G20750','AT2G21650','AT2G23000','AT2G24120','AT2G26010','AT2G26020','AT2G30360','AT2G30362','AT2G30760','AT2G30766','AT2G30770','AT2G33020','AT2G34170','AT2G34430','AT2G36050','AT2G36885','AT2G37770','AT2G38240','AT2G40300','AT2G40610','AT2G41240','AT2G41800','AT2G42380','AT2G42870','AT2G44740','AT3G06070','AT3G06530','AT3G09260','AT3G10020','AT3G11480','AT3G13470','AT3G14900','AT3G15450','AT3G16360','AT3G16670','AT3G18290','AT3G20440','AT3G21500','AT3G22142','AT3G22840','AT3G24982','AT3G25180','AT3G25290','AT3G27060','AT3G27630','AT3G27690','AT3G29590','AT3G30396','AT3G42716','AT3G42719','AT3G42720','AT3G43304','AT3G43670','AT3G47340','AT3G48520','AT3G49160','AT3G50770','AT3G51238','AT3G51240','AT3G51860','AT3G54090','AT3G55240','AT3G55970','AT3G56970','AT3G56980','AT3G57150','AT3G58120','AT3G59220','AT4G01985','AT4G02330','AT4G02390','AT4G03790','AT4G05180','AT4G06720','AT4G07850','AT4G08050','AT4G08091','AT4G09820','AT4G10290','AT4G12480','AT4G14090','AT4G14690','AT4G15210','AT4G16370','AT4G19130','AT4G22880','AT4G22960','AT4G23680','AT4G25000','AT4G25100','AT4G25630','AT4G27450','AT4G30290','AT4G32280','AT4G34250','AT4G34510','AT4G35770','AT4G36010','AT4G37990','AT4G39795','AT5G01600','AT5G02160','AT5G04150','AT5G05250','AT5G05270','AT5G07990','AT5G08610','AT5G13170','AT5G13930','AT5G17170','AT5G17220','AT5G18050','AT5G22580','AT5G24280','AT5G25980','AT5G26000','AT5G27060','AT5G28165','AT5G28464','AT5G32516','AT5G32624','AT5G33389','AT5G37260','AT5G42800','AT5G44420','AT5G44430','AT5G45082','AT5G48900','AT5G49360','AT5G49480','AT5G51720','AT5G52050','AT5G52390','AT5G54060','AT5G54270','AT5G54280','AT5G59310','AT5G59320','AT5G60250','AT5G64040','AT5G65080','AT5G66640')
eg = bitr(labelHM, fromType="TAIR", toType="SYMBOL", OrgDb="org.At.tair.db")
eg
cnet$data
names <- keggGet(kk$ID[2])[[1]]$GENE
names
namesodd <-  names[seq(0,length(names),2)]
namestrue <- gsub("\\;.*","",namesodd)
write.csv(namestrue, file = "ath00941",quote = F, row.names = F)


#gtf_file_read_to_be_used_in_statistics
setwd("F:/Illumina_transcriptome_Alignment")
gtf <- rtracklayer::import('genes.gtf')
gtf_df=as.data.frame(gtf)
write.table(gtf_df, 'gtf_df.txt', sep = ',', row.names = F, col.names = T, quote = F)