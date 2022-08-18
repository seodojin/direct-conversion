# Reference : <https://satijalab.org/seurat/articles/pbmc3k_tutorial.html>
# Setup the Seurat Object

library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)
# memory.limit()
memory.limit(size = 50000)


counts = read.table("1129_neuron_SampleTag01_hs_empty_RSEC_MolsPerCell.csv", skip = 5, sep = ",", header = TRUE, row.names = 1)
empty = CreateSeuratObject(counts = t(counts), project = "empty")

counts_sh = read.table("1129_neuron_SampleTag02_hs_sh_RSEC_MolsPerCell.csv", skip = 5, sep = ",", header = TRUE, row.names = 1)
ptbp1 = CreateSeuratObject(counts = t(counts_sh), project = "ptbp1")

plus <- merge(empty, ptbp1, add.cell.ids = c("empty","ptbp1"), project = "both")



plus[["percent.mt"]] = PercentageFeatureSet(plus, pattern = "^MT.")

head(plus@meta.data)

VlnPlot(plus, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

# visualize feature-feature relationships
plot1 <- FeatureScatter(plus, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(plus, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1 + plot2

plus <- subset(plus, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 20)

# Normalizing the data
plus <- NormalizeData(plus, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features (feature selection)
plus <- FindVariableFeatures(plus, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(plus), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(plus)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# scaling the data

plus <- ScaleData(plus)

# linear dimensional reduction
plus <- RunPCA(plus, features = VariableFeatures(object = plus))
VizDimLoadings(plus, dims = 1:2, reduction = "pca")
DimPlot(plus, reduction = "pca")
DimHeatmap(plus, dims = 1, cells = 500, balanced = T)

# Determine the ¡®dimensionality¡¯ of the dataset
plus <- JackStraw(plus, num.replicate = 100)
plus <- ScoreJackStraw(plus, dims = 1:20)
JackStrawPlot(plus, dims = 15)
ElbowPlot(plus)

# cluster the cells
plus <- FindNeighbors(plus, dims = 1:10)
plus <- FindClusters(plus, resolution = 0.2)
head(Idents(plus), 5)

# Run non-linear dimensional reduction (UMAP)
plus <- RunUMAP(plus, dims = 1:10)
DimPlot(plus, reduction = "umap")
DimPlot(plus, reduction = "pca")
DimPlot(plus, split.by = "orig.ident", reduction = "umap")

saveRDS(plus, "20220805")


# Reference : <https://github.com/IanevskiAleksandr/sc-type/>
# set up Assigning cell type identity to clusters
# automatically assign cell types using ScType
library(openxlsx)
library(tidyverse)
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

# load gene set preparation function
source("D:/seodojin/Rworld/20220509 sctype function-1.R")
# load cell type annotation function
source("D:/seodojin/Rworld/20220509 sctype function-2.R")

# DB file
db_ = "D:/seodojin/Rworld/ScTypeDB_full.xlsx";
tissue = "Brain" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = plus[["RNA"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(plus@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(plus@meta.data[plus@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(plus@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

# overlay the identified cell types on UMAP plot
plus@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  plus@meta.data$customclassif[plus@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(plus, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  

saveRDS(plus, "20220805-1")


# visualize data
clusters <- DimPlot(plus, reduction = "umap", 
                    group.by = "customclassif", label = T)
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

# reference : <https://satijalab.org/seurat/articles/de_vignette.html>
# Perform DE analysis using alternative tests
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)


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
DoHeatmap(plus, features = top10$gene) + NoLegend()
DoHeatmap(plus, features = top10$gene, label = F) 

# visualizing marker expression
VlnPlot(plus, features = "PTBP1") + NoLegend()
FeaturePlot(plus, features = "PTBP1")


FeaturePlot(plus, features = c("HMOX1","MMP1","GDF15","SLC6A15","LINC00520",
                               "SERPINB2","BMP2", "SPON2","SQSTM1"))
# Various gene diagrams are integrated
Gene1 = c("HMOX1","MMP1","GDF15","HSPB7","LINC00520",
          "CYSTM1","SLC6A15","TP53I11",
          "G0S2","IGFBP5", "PDK4","RGS4")
VlnPlot(plus,features=Gene1, pt.size = 0, stack=T, flip=T) + NoLegend()

Gene2 = c("SERPINB2","BMP2", "SPON2",
          "SQSTM1","WISP2","CYSTM1","PBX3")
VlnPlot(plus,features=Gene2, pt.size = 0, stack=T, flip=T) + NoLegend()


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

class(plus.sce)
structure(plus.sce)

table(plus.sce$customclassif)

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
  xlab("PC1") + ylab("PC2") + ggtitle("PC biplot")

plus.sce$pseudotime_PC1 <- rank(plus.sce$PC1)  
ggplot(as.data.frame(colData(plus.sce)), aes(x = pseudotime_PC1, y = customclassif, 
                                             colour = customclassif)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("PC1") + ylab("Timepoint") +
  ggtitle("Cells ordered by first principal component")

ggsave("Cells ordered by PC1.png", dpi = 600)

# Diffusion map pseudotime
sdj <- logcounts(plus.sce) 
cellLabels <- plus.sce$customclassif
colnames(sdj) <- cellLabels

sdj=as.data.frame(summary(sdj))
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
  xlab("Diffusion component 1 (DC1)") + ylab("Timepoint") +
  ggtitle("Cells ordered by DC1") + NoLegend()
ggsave("timepoint_DC1.png", dpi = 600)



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
ggsave("cells ordered by slingshot pseudotime.png", dpi = 600)


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

# Plot correlation between different pseudotime measures.
corrplot.mixed(cor(df_pseudotime, use = "na.or.complete"), 
               order = "hclust", tl.col = "black",
               main = "Correlation matrix for pseudotime results",
               mar = c(0, 0, 3.1, 0))

# Visualize how some of the temporally expressed genes change in time.
plotExpression(plus.sce, "HMOX1", x = "slingPseudotime_1", 
               colour_by = "customclassif", show_violin = FALSE,
               show_smooth = TRUE)
plotExpression(plus.sce, "PTBP1", x = "slingPseudotime_1", 
               colour_by = "customclassif", show_violin = FALSE,
               show_smooth = TRUE)
