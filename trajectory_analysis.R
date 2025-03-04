# Load required libraries
library(SingleCellExperiment)
library(Matrix)
library(tradeSeq)
library(Seurat)
library(slingshot)

# Load Seurat object
plus <- readRDS("seurat_object2")

# Select shPTBP1 cells
shPTBP1_cells <- rownames(plus@meta.data[plus@meta.data$orig.ident == "shPTBP1", ])

# Extract counts, normalized data, and metadata
counts <- plus[["RNA"]]$counts.shPTBP1[, shPTBP1_cells]
logcounts <- plus[["RNA"]]$data.shPTBP1[, shPTBP1_cells]
metadata <- plus@meta.data[shPTBP1_cells, ]

# Extract dimensional reductions
pca <- Embeddings(plus, "pca")[shPTBP1_cells, ]
umap <- Embeddings(plus, "umap")[shPTBP1_cells, ]

# Create SCE object
sce <- SingleCellExperiment(
  assays = list(counts = counts, logcounts = logcounts),
  colData = metadata,
  reducedDims = list(PCA = pca, UMAP = umap)
)

# Merge neuron types into "Neurons"
sce$updated_customclassif <- ifelse(
  sce$customclassif %in% c("GABAergic neurons", "Glutamatergic neurons"),
  "Neurons",
  sce$customclassif
)

# Gene filtering: keep genes expressed in at least 3 cells with >=10 counts
counts_sparse <- as(assay(sce, "counts"), "sparseMatrix")
gene_filter <- rowSums(counts_sparse >= 10) >= 3
sce <- sce[gene_filter, ]

# Cell filtering: keep cells with at least 500 expressed genes
cell_filter <- colSums(assay(sce, "counts") > 0) >= 500
sce <- sce[, cell_filter]

# Run Slingshot on filtered data
sds <- slingshot(
  sce,
  clusterLabels = 'updated_customclassif',
  reducedDim = 'UMAP',
  start.clus = "Fibroblasts",
  end.clus = c("Myofibroblasts", "Immature neurons", "Neurons"),
  allow.breaks = TRUE,
  extend = 'n',
  approx_points = 300
)

# Extract pseudotime and cell weights
pseudotime <- slingPseudotime(sds)
weights <- slingCurveWeights(sds)

# Keep cells with valid pseudotime
valid_cells <- rowSums(is.na(pseudotime)) == 0
pseudotime <- pseudotime[valid_cells, ]
weights <- weights[valid_cells, ]
sce_final <- sce[, valid_cells]

# Normalize pseudotime (0-1 range)
normalize_pseudotime <- function(pt) {
  apply(pt, 2, function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  })
}
pseudotime <- normalize_pseudotime(pseudotime)

# Normalize cell weights to sum to 1
weights <- weights / rowSums(weights)

# Add pseudotime and weights to colData
colData(sce_final) <- cbind(
  colData(sce_final),
  as.data.frame(pseudotime),
  as.data.frame(weights)
)

# Fit tradeSeq model
sce_tradeSeq <- fitGAM(
  counts = assay(sce_final, "counts"),
  pseudotime = as.matrix(pseudotime),
  cellWeights = as.matrix(weights),
  nknots = 5
)

