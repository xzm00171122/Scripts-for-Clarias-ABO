if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")

BiocManager::install(c("clusterProfiler", "AnnotationHub"))
BiocManager::install("Biobase")
BiocManager::install("IRanges")
BiocManager::install("topGO")
BiocManager::install("Rgraphviz")


library(clusterProfiler)
library(AnnotationHub)
library(topGO)
library(Rgraphviz)


hub <- AnnotationHub()

query(hub, "Cricetulus")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Dr.eg.db", version = "3.8")

library(org.Dr.eg.db)

resistantgene <- read.table("clarias.csv",sep = "\t")
resistantgene <- round(as.matrix(resistantgene))


#molecular function
ego <- enrichGO(gene= resistantgene,
                OrgDb= org.Dr.eg.db,
                keyType = "ENTREZID",
                ont = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05, qvalueCutoff = 0.05
)

enrich_go <- as.data.frame(ego)

write.csv(as.data.frame(ego), file="clarias_molecular function.csv")

barplot(ego, showCategory=10)

dotplot(ego,font.size=10, showCategory=10)

library("enrichplot")
cnetplot(ego)

#### enrichMap(ego)

?plotGOgraph
plotGOgraph(ego)

plotGOgraph(ego, firstSigNodes = 10, useInfo = "all", sigForAll = TRUE,
            useFullNames = TRUE)


#biology process
ego1 <- enrichGO(gene= resistantgene,
                 OrgDb= org.Dr.eg.db,
                 keyType = "ENTREZID",
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.05)

enrich_go1 <- as.data.frame(ego1)
write.csv(as.data.frame(ego1), file="clarias_biological process.csv")

barplot(ego1, showCategory=10)

dotplot(ego1,font.size=10)
enrichMap(ego1)
plotGOgraph(ego1)


#cellular component
ego2 <- enrichGO(gene= resistantgene,
                 OrgDb= org.Dr.eg.db,
                 keyType = "ENTREZID",
                 ont = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.05)

enrich_go2 <- as.data.frame(ego2)
write.csv(as.data.frame(ego2), file="clarias_cellular component.csv")

barplot(ego2, showCategory=10)

dotplot(ego2,font.size=10,showCategory=10)
enrichMap(ego2)
plotGOgraph(ego2)
