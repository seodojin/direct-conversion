# Load necessary libraries
library(SingleCellExperiment)
library(slingshot)
library(scater)
library(Seurat)
library(Rcpp)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(data.table)
library(Matrix)
library(tradeSeq)
library(splines)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

# Convert Seurat object to SingleCellExperiment object
plus <- readRDS("seurat_object2")
str(plus, max.level = 2)

# Extract and merge counts data
counts_shCtrl <- plus[["RNA"]]$counts.shCtrl
counts_shPTBP1 <- plus[["RNA"]]$counts.shPTBP1
counts <- cbind(counts_shCtrl, counts_shPTBP1)

# Extract and merge normalized data
data_shCtrl <- plus[["RNA"]]$data.shCtrl
data_shPTBP1 <- plus[["RNA"]]$data.shPTBP1
data <- cbind(data_shCtrl, data_shPTBP1)

# Create SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(counts = counts, logcounts = data),
  colData = plus@meta.data,
  reducedDims = list(PCA = Embeddings(plus, "pca"), UMAP = Embeddings(plus, "umap")))

# Add active identifiers
sce$cell_type <- plus@active.ident

# Scaling and PCA
sce <- logNormCounts(sce) # Log normalization (skip if already done)
sce <- scater::runPCA(sce, exprs_values = "logcounts") 

# Define cluster label update function
update_cluster_labels <- function(label) {
  if (label %in% c("GABAergic neurons", "Glutamatergic neurons")) {
    return("Neurons")
  } else {
    return(label)
  }
}

# Update cluster labels in SingleCellExperiment object
sce$updated_customclassif <- sapply(sce$customclassif, update_cluster_labels)
unique_clusters <- unique(sce$updated_customclassif)
print(unique_clusters)

# Run Slingshot with updated labels
sds <- slingshot(sce, clusterLabels = 'updated_customclassif', reducedDim = 'UMAP',
                 start.clus = "Fibroblasts", end.clus = c("Neurons", "Myofibroblasts", "Immature neurons"))
print(sds)
print(class(sds))
print(slingLineages(sds))

# Define the color for each lineage based on the correct matching
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

print(p)

ggsave("20240806_trajectory_umap_plot.png", plot = p, width = 7, height = 6, dpi = 1000)
print(slingLineages(sds))

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

plot_data$Lineage <- NA
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
  ))

# plot_data_long 생성
plot_data_long <- plot_data %>%
  select(Cell, Cluster, Lineage, Pseudotime)

# Pseudotime에 따른 클러스터 분포 시각화
ggplot(plot_data_long, aes(x = Pseudotime, y = Cluster, fill = Cluster)) +
  geom_violin(scale = "width", adjust = 1.5) +
#  facet_wrap(~ Lineage, scales = "free_y") +
  labs(title = "Pseudotime vs Cluster", x = "Pseudotime", y = "Cluster") +
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 45, hjust = 1))

######################### 

# Merge pseudotime data with primary trajectory information
plot_data <- cbind(plot_data, pseudotime)

# Reshape data for plotting
plot_data_long <- plot_data %>%
  pivot_longer(cols = starts_with("Lineage"), 
               names_to = "Trajectory", 
               values_to = "Pseudotime") %>%

plot_data_long_clean <- plot_data %>%
  select(Cell, Cluster, Lineage, Pseudotime) %>%
  filter(!is.na(Pseudotime))

# Generate table for checking
print(table(plot_data_long_clean$Cluster, plot_data_long_clean$Lineage))

# Define cell colors
cell_colors <- c("Immature neurons" = "#00bef3", 
                 "Myofibroblasts" = "#ff8c8c", 
                 "Fibroblasts" = "#19c3a3", 
                 "Neurons" = "#d4a600")
names(cell_colors) <- unique(plot_data_long_clean$Cluster)

# Create violin plot
ggplot(plot_data_long_clean, aes(x = Pseudotime, y = Cluster, fill = Cluster)) +
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
    axis.title.y = element_text(size = 14))

ggsave("20240806_pseudotime_violin_plot.png", width = 7, height = 6, dpi = 1000)

# Save the SingleCellExperiment object
saveRDS(sds, file = "20240730_slingshot_sds_results.rds")
sds <- readRDS("20240730_slingshot_sds_results.rds")
print(sds)
print(class(sds))
print(slingLineages(sds))

# Additional steps for sampled cells and GAM fitting
set.seed(123)
# 각 클러스터의 10%를 선택하여 샘플링
sampled_cells <- as.data.table(plus@meta.data)[, {
  num_cells <- .N  # 해당 클러스터의 총 세포 수
  sample_size <- ceiling(num_cells * 0.1)  # 10%에 해당하는 세포 수 계산
  .(cell = .I[sample(.N, sample_size)])  # 해당 클러스터에서 샘플링된 세포의 인덱스
}, by = seurat_clusters]$cell

# 샘플링된 세포를 사용하여 Seurat 객체를 서브셋
plus_sampled <- subset(plus, cells = sampled_cells)

# 샘플링 결과 확인
table(plus_sampled@meta.data$seurat_clusters)

counts_shCtrl <- plus_sampled[["RNA"]]$counts.shCtrl
counts_shPTBP1 <- plus_sampled[["RNA"]]$counts.shPTBP1
counts_sampled <- cbind(counts_shCtrl, counts_shPTBP1)
counts_sparse_sampled <- as(counts_sampled, "sparseMatrix")

min_counts <- 10
min_cells <- 3
filtered_counts <- counts_sparse_sampled[rowSums(counts_sparse_sampled >= min_counts) >= min_cells, ]

umap_coords <- Embeddings(plus_sampled, reduction = "umap")
cell_metadata <- plus_sampled@meta.data
sce_sampled <- SingleCellExperiment(assays = list(counts = filtered_counts), colData = cell_metadata, reducedDims = list(UMAP = umap_coords))

sce_sampled$updated_customclassif <- sapply(sce_sampled$customclassif, update_cluster_labels)

set.seed(123)
sds_sampled <- slingshot(sce_sampled, clusterLabels = 'updated_customclassif', reducedDim = 'UMAP',
                         start.clus = "Fibroblasts", end.clus = c("Neurons", "Myofibroblasts", "Immature neurons"))
print(slingLineages(sds_sampled))

sds <- SlingshotDataSet(sds_sampled)

sce_fitted <- fitGAM(counts = assay(sce_sampled, "counts"), sds = sds, nknots = 6, verbose = TRUE, parallel = FALSE)
print(str(sce_fitted))
print(names(sce_fitted))
metadata(sce_fitted)$slingshot <- sds
print(metadata(sce_fitted)$slingshot)

n_lineages <- length(slingLineages(metadata(sce_fitted)$slingshot))
print(paste("Number of lineages:", n_lineages))

saveRDS(list(sce_fitted = sce_fitted, slingshot_data = metadata(sce_fitted)$slingshot), 
        file = "20240802_fitGAM_results010.rds")

pseudotime <- slingPseudotime(sds_sampled)
print(summary(pseudotime))

