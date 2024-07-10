### store metadata for a package
### reference : https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
## input file : 1129_neuron_SampleTag01_hs_empty_RSEC_MolsPerCell.csv, 1129_neuron_SampleTag02_hs_sh_RSEC_MolsPerCell.csv
## output file : seurat_object.rds                
# load package
library(dplyr)
library(Seurat)
library(patchwork)

load("data/preprocessing.RData")
ls()

# Read scRNA-seq data: control (empty vector) 
# counts = read.table("data/1129_neuron_SampleTag01_hs_empty_RSEC_MolsPerCell.csv", skip = 5, sep = ",", header = TRUE, row.names = 1)
# shCtrl = CreateSeuratObject(counts = t(counts), project = "shCtrl")

# Read scRNA-seq data: experimental group (sh-PTBP1) 
# counts_sh = read.table("data/1129_neuron_SampleTag02_hs_sh_RSEC_MolsPerCell.csv", skip = 5, sep = ",", header = TRUE, row.names = 1)
# shPTBP1 = CreateSeuratObject(counts = t(counts_sh), project = "shPTBP1")

plus <- merge(shCtrl, shPTBP1, add.cell.ids = c("shCtrl","shPTBP1"), project = "both")

# QC and selecting cells
plus[["percent.mt"]] = PercentageFeatureSet(plus, pattern = "^MT.")

# head(plus@meta.data)
# tail(plus@meta.data)
# 
# VlnPlot(plus, 
#         features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
#         ncol = 3)

# filter cells : feature counts >1000 & <4000,  < 20% mitochondrial counts
plus <- subset(plus, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 20)

# Normalizing the data
plus_normalized <- NormalizeData(plus, normalization.method = "LogNormalize", scale.factor = 10000)

# find variable features 
plus2 <- FindVariableFeatures(plus_normalized, selection.method = "vst", nfeatures = 2000)

# scaling the data
plus_scaled <- ScaleData(plus2)

# linear dimensional reduction
plus_reduced <- RunPCA(plus_scaled, features = VariableFeatures(object = plus_scaled))
# DimPlot(plus_reduced, reduction = "pca")

# cluster the cells
plus_reduced <- FindNeighbors(plus_reduced, dims = 1:10)
plus_reduced <- FindClusters(plus_reduced, resolution = 0.2)

# Run non-linear dimensional reduction (UMAP)
plus_umap <- RunUMAP(plus_reduced, dims = 1:10)
# DimPlot(plus_umap, reduction = "umap")
DimPlot(plus_umap, split.by = "orig.ident", reduction = "umap")

# save RDS
saveRDS(plus_umap, "data/seurat_object")
