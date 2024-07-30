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
  reducedDims = list(PCA = Embeddings(plus, "pca"), UMAP = Embeddings(plus, "umap"))
)

# Add active identifiers
sce$cell_type <- plus@active.ident

# Scaling and PCA
sce <- logNormCounts(sce) # Log normalization (skip if already done)
sce <- scater::runPCA(sce, exprs_values = "logcounts") # Run PCA

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
cell_colors <- c("Immature neurons" = "#e41a1c", 
                 "Myofibroblasts" = "#ff7f00", 
                 "Fibroblasts" = "#d147a3", 
                 "Neurons" = "#3772b8")

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

ggsave("20240730_trajectory_umap_plot.png", plot = p, width = 7, height = 6, dpi = 1000)
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

plot_data_long_clean <- plot_data %>%
  select(Cell, Cluster, Lineage, Pseudotime) %>%
  filter(!is.na(Pseudotime))

# Generate table for checking
print(table(plot_data_long_clean$Cluster, plot_data_long_clean$Lineage))

# Define cell colors
cell_colors <- c("Immature neurons" = "#e41a1c", 
                 "Myofibroblasts" = "#ff7f00", 
                 "Fibroblasts" = "#d147a3", 
                 "Neurons" = "#3772b8")
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
    axis.title.y = element_text(size = 14)
  )

ggsave("20240730_pseudotime_violin_plot.png", width = 7, height = 6, dpi = 1000)

# Save the SingleCellExperiment object
saveRDS(sds, file = "20240730_slingshot_sds_results.rds")
sds <- readRDS("20240730_slingshot_sds_results.rds")
print(sds)
print(class(sds))
print(slingLineages(sds))

# Additional steps for sampled cells and GAM fitting
set.seed(123)
sampled_cells <- as.data.table(plus@meta.data)[, .(cell = .I[sample(.N, min(.N, 50))]), by = seurat_clusters]$cell
plus_sampled <- subset(plus, cells = sampled_cells)

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

sce_fitted <- fitGAM(counts = assay(sce_sampled, "counts"), sds = sds, nknots = 5, verbose = TRUE, parallel = FALSE)
print(str(sce_fitted))
print(names(sce_fitted))
metadata(sce_fitted)$slingshot <- sds
print(metadata(sce_fitted)$slingshot)

n_lineages <- length(slingLineages(metadata(sce_fitted)$slingshot))
print(paste("Number of lineages:", n_lineages))

saveRDS(list(sce_fitted = sce_fitted, slingshot_data = metadata(sce_fitted)$slingshot), file = "20240730_fitGAM_results50.rds")

# Perform patternTest
pattern_res <- patternTest(sce_fitted)

# Extract results for curve 3
curve3_pattern <- pattern_res[, 3]

# Assign gene names to curve3_pattern
names(curve3_pattern) <- rownames(pattern_res)

# Print the top 10 significant genes for curve 3 by p-value
print("Top 10 significant genes for curve 3 by p-value:")
print(head(sort(curve3_pattern), 10))

# Select significant genes with p-value < 0.05 for curve 3
significant_genes_pvalue <- names(curve3_pattern)[curve3_pattern < 0.05]
print(paste("Number of significant genes for curve 3 (p-value < 0.05):", length(significant_genes_pvalue)))

# Calculate FDR and select significant genes with FDR < 0.05 for curve 3
fdr_values <- p.adjust(curve3_pattern, method = "fdr")
significant_genes_fdr_0.05 <- names(curve3_pattern)[fdr_values < 0.05]
print(paste("Number of significant genes for curve 3 (FDR < 0.05):", length(significant_genes_fdr_0.05)))

# Select significant genes for GO analysis (using FDR < 0.05)
significant_genes <- significant_genes_fdr_0.05

# Perform GO analysis
go_results <- enrichGO(gene = significant_genes,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",  # Biological Process
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05)

# Check the results
go_results_df <- as.data.frame(go_results)
print(head(go_results_df, 10))

# Calculate term similarity
go_results_sim <- pairwise_termsim(go_results)

# Create and save emapplot
p <- emapplot(go_results_sim, showCategory = 20, color = "p.adjust",
              layout.params = list(layout = "kk"))
p
ggsave("20240726_emapplot_high_res.png", plot = p, width = 12, height = 10, dpi = 1000, units = "in")

# Save the results to a dataframe
write.csv(go_results_df, file = "20240726_GO_analysis_results_curve3_FDR0.05.csv", row.names = FALSE)

# Perform associationTest
# Select significant genes with FDR < 0.05
wald_stats <- associationTest(sce_fitted, lineages = TRUE)

# Remove NA values
wald_stats <- wald_stats[!is.na(wald_stats$waldStat), ]

# Select top 10 genes by p-value
top_genes <- rownames(wald_stats[order(wald_stats$pvalue, decreasing = FALSE), ])[1:10]

# Check the selected genes
print(top_genes)

# Plot the expression patterns of the top 10 genes along pseudotime
top_genes <- rownames(wald_stats[order(wald_stats$pvalue, decreasing = FALSE), ])[1:10]

# List to store plots
plots <- list()

# Run plotSmoothers for each gene
for (gene in top_genes) {
  p <- plotSmoothers(sce_fitted, counts = assay(sce_sampled, "counts"), gene = gene, lwd = 1.5) +
    ggtitle(paste("Expression of", gene, "along Pseudotime"))
  plots[[gene]] <- p
}

# Check the plots
plots[[1]]  # Plot for the first gene

# Perform earlyDETest
early_de_results <- earlyDETest(sce_fitted)

# Extract significant genes with FDR < 0.05
significant_genes_early_de <- rownames(early_de_results)[p.adjust(early_de_results$pvalue, method = "fdr") < 0.05]
print(paste("Number of significant genes from earlyDETest (FDR < 0.05):", length(significant_genes_early_de)))

# Calculate FDR from p-values
early_de_results$FDR <- p.adjust(early_de_results$pvalue, method = "fdr")

# Select top significant genes based on FDR
top_genes <- rownames(early_de_results[order(early_de_results$FDR), ])[1:20]

# Check the selected top genes
print(top_genes)

# Extract expression data for top genes
expression_data <- assay(sce_sampled, "counts")[top_genes, ]

# Log-transform the data
log_expression_data <- log1p(expression_data)

# Extract cluster information from the metadata
cluster_info <- colData(sce_sampled)$updated_customclassif

# Combine log-transformed expression data with cluster information
log_expression_df <- as.data.frame(t(log_expression_data))
log_expression_df$Cluster <- cluster_info

# Calculate mean expression per cluster
cluster_expression <- log_expression_df %>%
  group_by(Cluster) %>%
  summarise(across(everything(), mean))

# Convert to matrix and set row names
cluster_expression_matrix <- as.matrix(cluster_expression[,-1])
rownames(cluster_expression_matrix) <- cluster_expression$Cluster

# Save the heatmap to a file
png(filename = "20240730_Top_20_Significant_Genes_Heatmap_by_Cluster.png", width = 8, height = 6, units = "in", res = 600)
pheatmap(cluster_expression_matrix, scale = "row", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", clustering_method = "complete", 
         show_rownames = TRUE, show_colnames = TRUE, 
         main = "Top 20 Significant Genes Heatmap by Cluster")
dev.off()
