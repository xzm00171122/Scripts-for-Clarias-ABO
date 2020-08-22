library("BiocGenerics")
library("parallel")
library("Biobase")
#if (!requireNamespace("BiocManager", quietly = TRUE)
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("Mfuzz")
library(Mfuzz)
browseVignettes("Mfuzz")
x <- read.table(
  "channel_DEGs_fpkm_avgs.txt",
  row.names = 1,
  sep = "\t",
  header = T)
count <- data.matrix(x)
eset <- new("ExpressionSet",exprs = count)
eset <- filter.std(eset,min.std=0)
eset <- standardise(eset)
c <- 10
m <- mestimate(eset)
cl <- mfuzz(eset, c = c, m = m)
cl$size
cl$cluster[cl$cluster == 2]
# write the cluster data to csv
for (ii in 1:10){
  filename <- paste(ii, ".csv", sep="")
  write.csv(cl$cluster[cl$cluster == ii],filename)
}
cl$membership
# plot pdf
pdf("Mfuzz.plot.pdf")
mfuzz.plot(
  eset,
  cl,
  mfrow=c(2,3),
  new.window= FALSE)

par(mfrow=c(4,3))
pdf("1.pdf")
#tiff("Mfuzz.plot.tiff", width=800, height=403)
mfuzz.plot(
  eset,
  cl,
  mfrow=c(4,3),
  new.window= FALSE)
dev.off()

