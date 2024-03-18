### Perform DE analysis
### reference : https://satijalab.org/seurat/articles/de_vignette.html

# load package
library(DESeq2)
library(Seurat)
library(SeuratData)
library(EnhancedVolcano)
library(dplyr)

# load data
plus <- readRDS("annotation_object")

# Assigning cell type identity to clusters
clusters <- DimPlot(plus, reduction = "umap", 
                    group.by = "seurat_clusters", label = T)
treat <- DimPlot(plus, reduction = "umap", group.by = "orig.ident")

new.cluster.ids <- c("Immature neurons", "Myofibroblasts", "Fibroblasts", "Unknown", 
                     "Fibroblasts", "Glutamatergic neurons",
                     "GABAergic neurons")
names(new.cluster.ids) <- levels(plus)
plus <- RenameIdents(plus, new.cluster.ids)

# find markers for every cluster compared to all remaining cells, report only the positive ones
plus.markers <- FindAllMarkers(plus, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
plus.markers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC) -> marker.gene.list
View(marker.gene.list)
write.csv(marker.gene.list, file = "marker gene list 50.csv")

# plotting the top 5 markers for each cluster
plus.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5

DoHeatmap(plus, features = top5$gene, label = F) 
ggsave("heatmap.png", dpi=600)

# visualizing marker expression
VlnPlot(plus, features = "PTBP1") + NoLegend()
ggsave("ptbp1_violinplot.png", dpi = 600)

# Test for DE features using the DESeq2 package
glu <- FindMarkers(plus, ident.1 = "Glutamatergic neurons", 
                   ident.2 = NULL, only.pos = TRUE,
                   test.use = "DESeq2", max.cells.per.ident = 50)
View(glu)
write.csv(glu,file = "glutamatergic.csv")

gaba <- FindMarkers(plus, ident.1 = "GABAergic neurons", 
                    ident.2 = NULL, only.pos = TRUE,
                    test.use = "DESeq2", max.cells.per.ident = 50)
View(gaba)
write.csv(gaba, file = "gabaergic.csv")

neuro <- FindMarkers(plus, ident.1 = c("GABAergic neurons","Glutamatergic neurons"),
                     ident.2 = NULL, only.pos = TRUE,
                     test.use = "DESeq2", max.cells.per.ident = 50)
View(neuro)
write.csv(neuro, file = "neuro.csv")


myo <- FindMarkers(plus, ident.1 = "Myofibroblasts", 
                   ident.2 = NULL, only.pos = TRUE,
                   test.use = "DESeq2", max.cells.per.ident = 50)
View(myo)
write.csv(myo, file = "Myofibroblasts.csv")

fib <- FindMarkers(plus, ident.1 = "Fibroblasts", 
                   ident.2 = NULL, only.pos = TRUE,
                   test.use = "DESeq2", max.cells.per.ident = 50)
View(fib)
write.csv(fib, file = "Fibroblasts.csv")

imn <- FindMarkers(plus, ident.1 = "Immature neurons", 
                   ident.2 = NULL, only.pos = TRUE,
                   test.use = "DESeq2", max.cells.per.ident = 50)
View(imn)
write.csv(imn, file = "Immature neurons.csv")

unk <- FindMarkers(plus, ident.1 = "Unknown", 
                   ident.2 = NULL, only.pos = TRUE,
                   test.use = "DESeq2", max.cells.per.ident = 50)
View(unk)
write.csv(unk, file = "Unknown.csv")

# volcano plot of DE genes 
neuro.df <- as.data.frame(neuro)
EnhancedVolcano(neuro, x="avg_log2FC", y = "p_val_adj", 
                lab = rownames(neuro), pCutoff = 1e-4, FCcutoff = 1,
                title = NULL, subtitle = NULL)

ggsave("volcanoplot.png", dpi=600)

### Gene Ontology - barplot

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

rownames(neuro[neuro$avg_log2FC > 1.5,])
genes_to_test <- rownames(neuro[neuro$avg_log2FC > 1.5,])

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db",
                       keyType = "SYMBOL", ont = "BP")


as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 10))
ggsave("barplot.png", dpi=600)
