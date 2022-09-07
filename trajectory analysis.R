### Converting to/from SingleCellExperiment
### reference : https://satijalab.org/seurat/articles/conversion_vignette.html

# load package
library(scater)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(patchwork)

# load data
plus <- readRDS("annotation_object")

plus.sce <- as.SingleCellExperiment(plus)

### Trajectory analysis
### reference : https://broadinstitute.github.io/2020_scWorkshop/trajectory-analysis.html

# load package
library(SingleCellExperiment)
library(destiny)
library(clusterExperiment)
library(gam)
library(corrplot)
library(ggplot2)
library(ggthemes)
library(remotes)
library(Matrix)


pca.sce <- runPCA(plus.sce,ncomponents=50)
pca <- reducedDim(plus.sce, "PCA")

plus.sce$PC1 <- pca[,1]
plus.sce$PC2 <- pca[,2]

# Plot PC biplot with cells colored by customclassif. 
library(ggbeeswarm)
ggplot(as.data.frame(colData(plus.sce)), aes(x = PC1, y = PC2, color = customclassif)) + geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("PC1") + ylab("PC2") + ggtitle("PC biplot") + 
  theme(axis.text.x = element_text(size=10, face = "bold"),
        axis.text.y = element_text(size =10, face = "bold"),
        legend.text = element_text(size = 10, face = "bold"), 
        legend.title = element_blank())
ggsave("PC_biplot.png", dpi = 600)

# colData(plus.sce) accesses the cell metadata DataFrame object for plus.sce.
plus.sce$pseudotime_PC1 <- rank(plus.sce$PC1)  
ggplot(as.data.frame(colData(plus.sce)), aes(x = pseudotime_PC1, y = customclassif, 
                                             colour = customclassif)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("PC1") + ylab(NULL) +
  ggtitle("Cells ordered by first principal component") +
  theme(axis.text.x = element_text(size=15, face = "bold"),
        axis.text.y = element_text(size =15, face = "bold")) + NoLegend()

ggsave("Cells_ordered_by_PC1.png", dpi = 600)

# Diffusion map pseudotime
sdj <- logcounts(plus.sce) 
cellLabels <- plus.sce$customclassif
colnames(sdj) <- cellLabels

sdj=as.data.frame(summary(sdj))

dm <- DiffusionMap(t(sdj))

rownames(pca) <- cellLabels
dm <- DiffusionMap(pca)
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
ggsave("cells_ordered_by_DC1.png", dpi = 600)

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
  ggtitle("Cells ordered by Slingshot pseudotime") + NoLegend()
ggsave("cells_ordered_by_slingshot_pseudotime.png", dpi = 600)

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
plot(reducedDims(plus.sce)$PCA, col = colors[cut(plus.sce$slingPseudotime_1,breaks=50)], pch=16, asp = 1)
lines(SlingshotDataSet(plus.sce), lwd=2)
