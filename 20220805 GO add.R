library(Seurat)
library(dplyr)
library(tidyverse)
library(DESeq2)


memory.limit(size = 50000)

# load data
plus <- readRDS("20220805-1")
str(plus)
View(plus@meta.data)



# visualize data
clusters <- DimPlot(plus, reduction = "umap", 
                    group.by = "seurat_clusters", label = T)
treat <- DimPlot(plus, reduction = "umap", group.by = "orig.ident")

clusters|treat


# Assigning cell type identity to clusters
new.cluster.ids <- c("Immature neurons", "Myofibroblasts", "Fibroblasts", "Unknown", 
                     "Fibroblasts", "Glutamatergic neurons",
                     "GABAergic neurons")
  names(new.cluster.ids) <- levels(plus)
plus <- RenameIdents(plus, new.cluster.ids)
DimPlot(plus, reduction = "umap", label = TRUE, pt.size = 0.5,
        split.by = "orig.ident") + NoLegend()
DimPlot(plus, reduction = "umap", label = F, pt.size = 0.5,
        split.by = "orig.ident") 
ggsave("20220808 legend add umap split version.png", dpi = 600)

# find markers for every cluster compared to all remaining cells, report only the positive ones

plus.markers <- FindAllMarkers(plus, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
plus.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC) -> marker.gene.list
View(marker.gene.list)
write.csv(marker.gene.list, file = "20220711 marker gene list 20.csv")

# plotting the top 10 markers for each cluster
plus.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

plus.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5

DoHeatmap(plus, features = top10$gene) + NoLegend()
DoHeatmap(plus, features = top5$gene, label = F) 
ggsave("20220817 heatmap top5.png", dpi=600)


# visualizing marker expression
VlnPlot(plus, features = "PTBP1") + NoLegend()
ggsave("20220808 ptbp1 violin plot.png", dpi = 600)
FeaturePlot(plus, features = "PTBP1")

# Various gene diagrams are integrated
Gene1 = c("HMOX1","MMP1","GDF15","HSPB7","LINC00520",
          "CYSTM1","SLC6A15","TP53I11",
          "G0S2","IGFBP5", "PDK4","RGS4")
VlnPlot(plus,features=Gene1, pt.size = 0, stack=T, flip=T) + NoLegend()
ggsave("20220808 gene1 violin plot.png", dpi = 600)
Gene2 = c("SERPINB2","BMP2", "SPON2",
          "SQSTM1","WISP2","CYSTM1","PBX3")
VlnPlot(plus,features=Gene2, pt.size = 0, stack=T, flip=T) + NoLegend()


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


# Gene Ontology - barplot

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

rownames(neuro[neuro$avg_log2FC > 0.5,])
genes_to_test <- rownames(neuro[neuro$avg_log2FC > 0.5,])

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db",
                       keyType = "SYMBOL", ont = "BP")


as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 10))
ggsave("20220805 GO barplot.png", dpi=600, width = 6.43, height = 6, 
       units = c("in"))

# volcano plot of DE genes 

library(EnhancedVolcano)

neuro.df <- as.data.frame(neuro)
EnhancedVolcano(neuro, x="avg_log2FC", y = "p_val_adj", 
                lab = rownames(neuro), pCutoff = 1e-4, FCcutoff = 1,
                title = NULL, subtitle = NULL)

selected = c("CA12", "FTH1", "PSAP", "HMOX1", "CTSK")
EnhancedVolcano(neuro, x="avg_log2FC", y = "p_val_adj", 
                lab = rownames(neuro), pCutoff = 1e-4, FCcutoff = 1,
                title = NULL, subtitle = NULL, selectLab = selected)

ggsave("20220805 volcano plot.png", dpi=600, width = 6.43, height = 6, 
       units = c("in"))


# Trajectory set up
# https://satijalab.org/seurat/articles/conversion_vignette.html
library(SingleCellExperiment)
library(destiny)
library(scater)
library(clusterExperiment)
library(gam)
library(corrplot)
library(ggplot2)
library(ggthemes)
library(remotes)
library(SeuratDisk)
library(SeuratData)
library(patchwork)
library(Matrix)
library(monocle)


# Trajectory analysis
# https://broadinstitute.github.io/2020_scWorkshop/trajectory-analysis.html

plus.sce <- as.SingleCellExperiment(plus)
p1 <- plotExpression(plus.sce, features = "HMOX1", x = "ident") + theme(axis.text.x = element_text(angle = 45,
                                                                                                   hjust = 1))
p2 <- plotPCA(plus.sce, colour_by = "ident")
p1 + p2


pca.sce <- runPCA(plus.sce,ncomponents=50)
pca <- reducedDim(plus.sce, "PCA")
head(pca)
dim(pca)

plus.sce$PC1 <- pca[,1]
plus.sce$PC2 <- pca[,2]
head(colData(plus.sce))
tail(colData(plus.sce))

# Plot PC biplot with cells colored by customclassif. 
# colData(plus.sce) accesses the cell metadata DataFrame object for plus.sce.


library(ggbeeswarm)
ggplot(as.data.frame(colData(plus.sce)), aes(x = PC1, y = PC2, color = customclassif)) + geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("PC1") + ylab("PC2") + ggtitle("PC biplot") + 
  theme(axis.text.x = element_text(size=10, face = "bold"),
        axis.text.y = element_text(size =10, face = "bold"),
        legend.text = element_text(size = 10, face = "bold"), 
        legend.title = element_blank())
ggsave("20220808 PC biplot.png", dpi = 600)


plus.sce$pseudotime_PC1 <- rank(plus.sce$PC1)  
ggplot(as.data.frame(colData(plus.sce)), aes(x = pseudotime_PC1, y = customclassif, 
                                             colour = customclassif)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("PC1") + ylab(NULL) +
  ggtitle("Cells ordered by first principal component") +
  theme(axis.text.x = element_text(size=15, face = "bold"),
        axis.text.y = element_text(size =15, face = "bold")) + NoLegend()

ggsave("20220808 Cells ordered by PC1.png", dpi = 600)


# Diffusion map pseudotime
sdj <- logcounts(plus.sce) 
cellLabels <- plus.sce$customclassif
colnames(sdj) <- cellLabels

sdj=as.data.frame(summary(sdj))

library(pheatmap)
pheatmap(sdj,
         cutree_rows=6,cutree_cols=6,
         cluster_rows = F, cluster_cols = F,
         show_colnames=F, show_rownames=T)

dm <- DiffusionMap(t(sdj))

rownames(pca) <- cellLabels
dm <- DiffusionMap(pca)
View(dm)
dpt <- DPT(dm)
plus.sce$pseudotime_diffusionmap <- rank(eigenvectors(dm)[,1])
ggplot(as.data.frame(colData(plus.sce)),
       aes(x=pseudotime_diffusionmap,
           y=customclassif, colour = customclassif)) +
  geom_quasirandom(groupOnX = F) +
  scale_color_tableau() + theme_classic() +
  xlab("Diffusion component 1 (DC1)") + ylab(NULL) +
  ggtitle("Cells ordered by DC1")  +
  theme(axis.text.x = element_text(size=10, face = "bold"),
        axis.text.y = element_text(size =10, face = "bold")) + NoLegend()
ggsave("20220808 cells ordered by DC1.png", dpi = 600)

library(slingshot)
sce <- slingshot(plus.sce, reducedDim = 'PCA')

# Plot PC1 vs PC2 colored by Slingshot pseudotime.
colors <- rainbow(50, alpha = 1)
plot(reducedDims(sce)$PCA, col = colors[cut(sce$slingPseudotime_1,breaks=50)], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2)

# Plot Slingshot pseudotime vs cell stage. 
ggplot(as.data.frame(colData(plus.sce)), aes(x = sce$slingPseudotime_1, y = customclassif, 
                                             colour = customclassif)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("Slingshot pseudotime") + ylab(NULL)+
  ggtitle("Cells ordered by Slingshot pseudotime") + 
  theme(axis.text.x = element_text(size=10, face = "bold"),
        axis.text.y = element_text(size =10, face = "bold")) + NoLegend()
ggsave("20220808 cells ordered by slingshot pseudotime.png", dpi = 600)

# Cluster cells using the Seurat workflow below.
gcdata <- CreateSeuratObject(counts = counts(plus.sce), project = "slingshot")

gcdata <- NormalizeData(gcdata, normalization.method = "LogNormalize", scale.factor = 10000)
gcdata <- FindVariableFeatures(gcdata, selection.method = "vst", nfeatures = 2000)
gcdata <- ScaleData(object = gcdata, do.center = T, do.scale = F)

gcdata <- RunPCA(gcdata, features = VariableFeatures(gcdata), npcs = 40, ndims.print = 1:5, nfeatures.print = 5)

# Cluster the cells using the first twenty principal components.
gcdata <- FindNeighbors(gcdata, reduction = "pca", dims = 1:20, k.param = 20)

gcdata <- FindClusters(gcdata, resolution = 0.2, algorithm = 1, random.seed = 100)

# Add clustering information from Seurat to the plus.sce object
plus.sce$slingPseudotime_1 <- NULL  # remove old slingshot pseudotime data
colData(plus.sce)$Seurat_clusters <- as.character(Idents(gcdata))  # go from factor to character
head(colData(plus.sce))

# Then run Slingshot using these cluster assignments.
plus.sce <- slingshot(plus.sce, clusterLabels = 'Seurat_clusters', reducedDim = 'PCA')

# Plot PC1 vs PC2 colored by Slingshot pseudotime.
colors <- rainbow(50, alpha = 1)
plot(reducedDims(plus.sce)$PCA, 
     col = colors[cut(plus.sce$slingPseudotime_1,breaks=50)], 
     pch=16, asp = 1, cex.axis=1.5)
lines(SlingshotDataSet(plus.sce), lwd=3)


# Visualize how some of the temporally expressed genes change in time.
plotExpression(plus.sce, "HMOX1", x = "slingPseudotime_1", 
               colour_by = "customclassif", show_violin = FALSE,
               show_smooth = TRUE)
plotExpression(plus.sce, "PTBP1", x = "slingPseudotime_1", 
               colour_by = "customclassif", show_violin = FALSE,
               show_smooth = TRUE)
