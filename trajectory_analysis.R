# Load necessary libraries
library(SingleCellExperiment)
library(slingshot)
library(scater)
library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)
library(tradeSeq)

# Convert Seurat object to SingleCellExperiment object
plus <- readRDS("seurat_object2")

# Extract and merge counts and normalized data
counts <- cbind(plus[["RNA"]]$counts.shCtrl, plus[["RNA"]]$counts.shPTBP1)
data <- cbind(plus[["RNA"]]$data.shCtrl, plus[["RNA"]]$data.shPTBP1)

# Create SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(counts = counts, logcounts = data),
  colData = plus@meta.data,
  reducedDims = list(PCA = Embeddings(plus, "pca"), UMAP = Embeddings(plus, "umap"))
)

# Add active identifiers
sce$cell_type <- plus@active.ident

# Update cluster labels
update_cluster_labels <- function(label) {
  if (label %in% c("GABAergic neurons", "Glutamatergic neurons")) {
    return("Neurons")
  } else {
    return(label)
  }
}
sce$updated_customclassif <- sapply(sce$customclassif, update_cluster_labels)

# Run Slingshot with updated labels
sds <- slingshot(sce, clusterLabels = 'updated_customclassif', reducedDim = 'UMAP',
                 start.clus = "Fibroblasts", end.clus = c("Neurons", "Myofibroblasts", "Immature neurons"))

# Define the color for each lineage
lineage_colors <- c("Immature neurons" = "#440154", 
                    "Myofibroblasts" = "#21908c", 
                    "Neurons" = "#fde725")

# Plot trajectories with ggplot2
umap_coords <- reducedDims(sce)$UMAP
plot_data <- data.frame(UMAP1 = umap_coords[,1], UMAP2 = umap_coords[,2], CellType = sce$updated_customclassif)
cell_colors <- c("Immature neurons" = "#00bef3", 
                 "Myofibroblasts" = "#ff8c8c", 
                 "Fibroblasts" = "#19c3a3", 
                 "Neurons" = "#d4a600")

p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(values = cell_colors) +
  theme_minimal() +
  labs(title = "Cell Types with Trajectories") +
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)))

# Add trajectory lines with the correct colors
for (i in seq_along(slingCurves(sds))) {
  curve_data <- slingCurves(sds)[[i]]$s[slingCurves(sds)[[i]]$ord, ]
  end_cluster <- slingLineages(sds)[[i]][length(slingLineages(sds)[[i]])]
  end_cluster_cells <- which(sce$updated_customclassif == end_cluster)
  end_point <- colMeans(umap_coords[end_cluster_cells,])
  distances <- sqrt(rowSums((curve_data - matrix(end_point, nrow = nrow(curve_data), ncol = 2, byrow = TRUE))^2))
  closest_point <- which.min(distances)
  p <- p + geom_path(data = data.frame(UMAP1 = curve_data[1:closest_point,1], UMAP2 = curve_data[1:closest_point,2]),
                     aes(x = UMAP1, y = UMAP2), color = lineage_colors[end_cluster], linewidth = 1, alpha = 0.7)
}

ggsave("20240806_trajectory_umap_plot.png", plot = p, width = 7, height = 6, dpi = 1000)

# Pseudotime violin plot
set.seed(123)
pseudotime <- as.data.frame(slingPseudotime(sds))
pseudotime$Cell <- rownames(pseudotime)

# Extract common cells
common_cells <- intersect(colnames(sce), rownames(pseudotime))

# Subset sce and pseudotime to include only common cells
sce_subset <- sce[, common_cells]
pseudotime_subset <- pseudotime[common_cells, ]

# Ensure the cluster information is aligned
plot_data <- data.frame(Cell = common_cells, 
                        Cluster = sce_subset$updated_customclassif, 
                        Pseudotime.Lineage1 = pseudotime_subset$Lineage1, 
                        Pseudotime.Lineage2 = pseudotime_subset$Lineage2, 
                        Pseudotime.Lineage3 = pseudotime_subset$Lineage3)

plot_data <- plot_data %>%
  mutate(Lineage = case_when(
    Cluster == "Myofibroblasts" ~ "Lineage1",
    Cluster == "Immature neurons" ~ "Lineage2",
    Cluster == "Neurons" ~ "Lineage3",
    Cluster == "Fibroblasts" ~ sample(c("Lineage1", "Lineage2", "Lineage3"), n(), replace = TRUE)
  )) %>%
  mutate(Pseudotime = case_when(
    Lineage == "Lineage1" ~ Pseudotime.Lineage1,
    Lineage == "Lineage2" ~ Pseudotime.Lineage2,
    Lineage == "Lineage3" ~ Pseudotime.Lineage3
  )) %>%
  select(Cell, Cluster, Lineage, Pseudotime) %>%
  filter(!is.na(Pseudotime))

# Create violin plot
ggplot(plot_data, aes(x = Pseudotime, y = Cluster, fill = Cluster)) +
  geom_violin(scale = "width", adjust = 1.5) +
  scale_fill_manual(values = cell_colors) +
  labs(title = "Pseudotime Distribution by Cluster", x = "Pseudotime", y = "Cluster") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

# Normalize pseudotime within each cluster
plot_data_normalized <- plot_data %>%
  group_by(Cluster) %>%
  mutate(Normalized_Pseudotime = (Pseudotime - min(Pseudotime)) / (max(Pseudotime) - min(Pseudotime))) %>%
  ungroup()

# Create violin plot with normalized pseudotime
ggplot(plot_data_normalized, aes(x = Normalized_Pseudotime, y = Cluster, fill = Cluster)) +
  geom_violin(scale = "width", adjust = 1.5) +
  scale_fill_manual(values = cell_colors) +
  labs(title = "Normalized Pseudotime Distribution by Cluster", x = "Normalized Pseudotime", y = "Cluster") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

ggsave("20240808_normalize_pseudotime_violin_plot.png", width = 7, height = 6, dpi = 1000)

# Filter genes by counts and cells
min_counts <- 10
min_cells <- 3
counts_sparse <- as(counts, "sparseMatrix")
filtered_counts <- counts_sparse[rowSums(counts_sparse >= min_counts) >= min_cells, ]

# Create SingleCellExperiment object for filtered data
sce_filtered <- SingleCellExperiment(
  assays = list(counts = filtered_counts, logcounts = logcounts(sce)[rownames(filtered_counts), ]),
  colData = colData(sce),
  reducedDims = reducedDims(sce)
)

# Run Slingshot with updated labels on filtered data
sds_filtered <- slingshot(sce_filtered, clusterLabels = 'updated_customclassif', reducedDim = 'UMAP',
                          start.clus = "Fibroblasts", end.clus = c("Neurons", "Myofibroblasts", "Immature neurons"))

# Ensure sds_filtered is a SlingshotDataSet
sds_filtered <- as.SlingshotDataSet(sds_filtered)

# Fit GAM models using the filtered dataset
sce_fitted <- fitGAM(counts = assay(sce_filtered, "counts"), sds = sds_filtered, nknots = 6, verbose = TRUE, parallel = FALSE)

# Save fitted GAM results with proper Slingshot data
metadata(sce_fitted)$slingshot <- sds_filtered
saveRDS(sce_fitted, file = "20240808_fitGAM_results_full.rds")
