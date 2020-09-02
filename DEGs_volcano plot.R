source("https://bioconductor.org/biocLite.R")
biocLite("pheatmap")
library("DESeq2")
library(ggplot2)
library(dplyr)
library(pheatmap)

## Input the count data
countData <- read.table("chan2_vs_5", header=TRUE, row.names=1)
dim(countData)
head(countData)

##Import the metadata that inlcudes the treatment assignments and if single or paired reads.
colData <- read.table("chan2_vs_5_sample.txt", header=TRUE, row.names=1)
dim(colData)
head(colData)

# check the order of treatment assignments in colData and countData
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

##e Creat the dataset and define model
?DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)
dds

##########remove genes with no counts
dds <- dds[rowSums(counts(dds)) > 1,]
dds

###########PCA plot
#rlog This function transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, 
#and which normalizes with respect to library size.
rld<-rlog(dds)
plotPCA(rld)
############

#perform DESeq, DESeq will do the following things:
#estimating dispersions
#found already estimated dispersions, replacing these
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
jumpy <- DESeq(dds)
jumpy
resultsNames(jumpy)
#make sure results look okay
res <- results(jumpy)
write.csv(as.data.frame(res),file="res1.csv")
head(res)
summary(res)

#reduce alpha to .05 for analysis of results
?results
res05 <- results(jumpy, alpha =.05)
head(res05)
names(res05)
#find out how many genes are differentially expressed below the maximum adjusted p-value cutoff
#sum(res05$padj < .05, na.rm=TRUE)
sum(res05$padj < .05, na.rm=TRUE)
#order all the genes by padj
#resOrdered <- res05[order(res05$padj),]
resOrdered <- res05[order(res05$padj),]
#store the number of DE genes as a variable
#numDEgenes <- sum(res05$padj < .05, na.rm=TRUE) should use padj
numDEgenes <- sum(res05$padj < .05, na.rm=TRUE)
numDEgenes
#get only the DEgenes ordered by padj
DEgenes <- head(resOrdered, numDEgenes)
head(DEgenes)
#write only these genes to a CSV file
write.csv(as.data.frame(DEgenes), file="starfishtopDEgenes1.csv")


#makes volcano plot of data and copies to pdf file#
?plotMA
plotMA(res05, main='DESeq2', ylim=c(-10,10))

#gets the number of numDEgenes (or top xxx most) differentially expressed genes#
select <- order(rowMeans(counts(jumpy,normalized=T)),decreasing=T)[1:numDEgenes]

head(select)
#transforms data and gets sample information ready for heatmap#
?normTransform
nt <- normTransform(jumpy)
nt
#df <- as.data.frame(colData(jumpy)[,c('treatment','type')])
df <- as.data.frame(colData(jumpy)[,c("condition")])
df
colnames(df) <- c('condition') 
df
trans.counts <- assay(nt)[select,]

head(trans.counts)
?assay
##makes heatmap and prints to file###
rownames(df) <- colnames(trans.counts)
df
pheatmap(trans.counts, cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=df)

?pheatmap
###test first 100 and plot heatmap
select100 <- order(rowMeans(counts(jumpy,normalized=T)),decreasing=T)[1:100]
trans.counts100 <- assay(nt)[select100,]
pheatmap(trans.counts100, cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=df)

###########Principal component plot of the samples, no replicate, nothing
#vsd <- vst(jumpy, blind=FALSE)
#head(assay(vsd), 3)
#plotPCA(vsd, intgroup=c("condition", "type"))
########################

#####volcano.plot#################
alpha <- 0.05 # Threshold on the adjusted p-value
cols <- densCols(res$log2FoldChange, -log10(res$pvalue))

#plot(res$log2FoldChange, -log10(res$padj), col=cols, panel.first=grid(),main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",pch=20, cex=0.6)
plot(res$log2FoldChange, -log10(res$pvalue), col=cols, panel.first=grid(),main="Volcano plot", xlab="log2(fold-change)", ylab="-log10(adjusted p-value)",pch=19, cex=0.6)

abline(v=0)
abline(v=c(-1,1), col="red")
abline(h=-log10(alpha), col="blue")
?abline
#gn.selected <- abs(res$log2FoldChange) > 2.5 & res$padj < alpha
#text(res$log2FoldChange[gn.selected],-log10(res$padj)[gn.selected],lab=rownames(res)[gn.selected ], cex=0.7)



###black and red volcano plot###
res <- read.csv("res1.csv", header = TRUE)
plot(res$log2FoldChange, -log10(res$padj),
     pch=19,cex=0.06,
     col=c("black","red")[res$level],xlab="Log2 Fold change", ylab="-log10(padj)", font.lab=2,cex.lab=1, xlim=c(-8,8), ylim=c(3,100),xaxt="n", yaxt="n")
a<- c(-10,-8,-6,-4,-2,0,2,4,6,8,10)
b <-c(-20,0,20,40,60,80,100)
axis(1,a, font.axis = 2,lwd=2)
axis(2,b, las=1,font.axis=2,lwd=2)
grid(nx=NULL,ny=NULL, col="lightgrey")
title("Channel day 2 VS. Channel day 5", cex = 1, col = "black", font = 3)
