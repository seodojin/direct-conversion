### store metadata for a package
### reference : https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
## input file : 1129_neuron_SampleTag01_hs_empty_RSEC_MolsPerCell.csv, 1129_neuron_SampleTag02_hs_sh_RSEC_MolsPerCell.csv
## output file : seurat_object.rds                
# load package
library(dplyr)
library(Seurat)
library(patchwork)

counts = read.table("1129_neuron_SampleTag01_hs_empty_RSEC_MolsPerCell.csv", skip = 5, sep = ",", header = TRUE, row.names = 1)
shCtrl = CreateSeuratObject(counts = t(counts), project = "shCtrl")

counts_sh = read.table("1129_neuron_SampleTag02_hs_sh_RSEC_MolsPerCell.csv", skip = 5, sep = ",", header = TRUE, row.names = 1)
shPTBP1 = CreateSeuratObject(counts = t(counts_sh), project = "shPTBP1")

plus <- merge(shCtrl, shPTBP1, add.cell.ids = c("shCtrl","shPTBP1"), project = "both")

# QC and selecting cells
plus[["percent.mt"]] = PercentageFeatureSet(plus, pattern = "^MT.")

head(plus@meta.data)
tail(plus@meta.data)

VlnPlot(plus, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

# filter cell : feature counts <4000 or >1000,  >20% mitochondrial counts
plus <- subset(plus, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 20)

# Normalizing the data
plus <- NormalizeData(plus, normalization.method = "LogNormalize", scale.factor = 10000)

# scaling the data
plus <- ScaleData(plus)

# linear dimensional reduction
plus <- RunPCA(plus, features = VariableFeatures(object = plus))
DimPlot(plus, reduction = "pca")

# cluster the cells
plus <- FindNeighbors(plus, dims = 1:10)
plus <- FindClusters(plus, resolution = 0.2)

# Run non-linear dimensional reduction (UMAP)
plus <- RunUMAP(plus, dims = 1:10)
DimPlot(plus, reduction = "umap")
DimPlot(plus, split.by = "orig.ident", reduction = "umap")

# save RDS
saveRDS(plus, "seurat_object")
