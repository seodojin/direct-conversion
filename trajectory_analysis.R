# Load required libraries
library(SingleCellExperiment)
library(slingshot)
library(scater)

# Convert Seurat object to SingleCellExperiment object

str(plus, max.level = 2)

# counts 데이터 추출 및 병합
counts_shCtrl <- plus[["RNA"]]$counts.shCtrl
counts_shPTBP1 <- plus[["RNA"]]$counts.shPTBP1
counts <- cbind(counts_shCtrl, counts_shPTBP1)

# data (normalized) 데이터 추출 및 병합
data_shCtrl <- plus[["RNA"]]$data.shCtrl
data_shPTBP1 <- plus[["RNA"]]$data.shPTBP1
data <- cbind(data_shCtrl, data_shPTBP1)

# SingleCellExperiment 객체 생성
sce <- SingleCellExperiment(
  assays = list(
    counts = counts,
    logcounts = data
  ),
  colData = plus@meta.data,
  reducedDims = list(
    PCA = Embeddings(plus, "pca"),
    UMAP = Embeddings(plus, "umap")
  )
)

# 활성 식별자 추가
sce$cell_type <- plus@active.ident

# scaling
sce <- logNormCounts(sce)  # 로그 정규화 (이미 수행되었다면 건너뛰기)
sce <- scater::runPCA(sce, exprs_values = "logcounts")  # PCA 실행

# Run Slingshot with specified start and end clusters
# 클러스터 레이블 업데이트 함수
update_cluster_labels <- function(label) {
  if (label %in% c("GABAergic neurons", "Glutamatergic neurons")) {
    return("Neurons")
  } else {
    return(label)
  }
}

# SingleCellExperiment 객체의 클러스터 레이블 업데이트
sce$updated_customclassif <- sapply(sce$customclassif, update_cluster_labels)

# 업데이트된 고유 클러스터 확인
unique_clusters <- unique(sce$updated_customclassif)
print(unique_clusters)

# Slingshot 실행 (업데이트된 레이블 사용)
sds <- slingshot(sce, clusterLabels = 'updated_customclassif', reducedDim = 'UMAP',
                 start.clus = "Fibroblasts", 
                 end.clus = c("Neurons", "Myofibroblasts", "Immature neurons"))

# 궤적 그리기
# UMAP 좌표 추출
umap_coords <- reducedDims(sce)$UMAP

# 데이터 프레임 생성
plot_data <- data.frame(
  UMAP1 = umap_coords[,1],
  UMAP2 = umap_coords[,2],
  CellType = sce$updated_customclassif
)

# 색상 정의
cell_colors <- c("Neurons" = "blue", "Immature neurons" = "yellow", 
                 "Myofibroblasts" = "green", "Fibroblasts" = "red")

# 궤적 색상 정의
trajectory_colors <- c("#FF9999", "#66B2FF", "#c0d84d")

# ggplot으로 그리기
p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(values = cell_colors) +
  theme_minimal() +
  labs(title = "Cell Types with Trajectories")

# 궤적 추가
for (i in seq_along(slingCurves(sds))) {
  curve_data <- slingCurves(sds)[[i]]$s[slingCurves(sds)[[i]]$ord, ]
  
  # 종착점 클러스터 찾기
  end_cluster <- slingLineages(sds)[[i]][length(slingLineages(sds)[[i]])]
  
  # 종착점 클러스터의 중심점 찾기
  end_cluster_cells <- which(sce$updated_customclassif == end_cluster)
  end_point <- colMeans(umap_coords[end_cluster_cells,])
  
  # 종착점까지의 거리 계산
  distances <- sqrt(rowSums((curve_data - matrix(end_point, nrow = nrow(curve_data), ncol = 2, byrow = TRUE))^2))
  
  # 종착점에 가장 가까운 점 찾기
  closest_point <- which.min(distances)
  
  # 종착점까지의 궤적만 그리기
  p <- p + geom_path(data = data.frame(UMAP1 = curve_data[1:closest_point,1], 
                                       UMAP2 = curve_data[1:closest_point,2]),
                     aes(x = UMAP1, y = UMAP2), 
                     color = trajectory_colors[i], 
                     linewidth = 1, 
                     alpha = 0.7)
}

p

ggsave("20240712_trajectory_umap_plot.png", plot = p, width = 5, height = 4, dpi = 1000)


########## 20240708 pseudotime vln plot
# load packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Pseudotime 값 추출
pseudotime <- slingPseudotime(sds)

# 데이터 프레임 생성
plot_data <- data.frame(
  Cell = rownames(pseudotime),
  Trajectory1 = pseudotime[,1],
  Trajectory2 = pseudotime[,2],
  Trajectory3 = pseudotime[,3],
  Cluster = sce$updated_customclassif
)

# 클러스터 색상 정의
cluster_colors <- c("Neurons" = "blue", 
                    "Immature neurons" = "yellow", 
                    "Myofibroblasts" = "green", 
                    "Fibroblasts" = "red")

# Long 형식으로 변환
plot_data_long <- plot_data %>%
  pivot_longer(cols = c(Trajectory1, Trajectory2, Trajectory3), 
               names_to = "Trajectory", 
               values_to = "Pseudotime")

# NA 값 제거
plot_data_long <- plot_data_long %>% filter(!is.na(Pseudotime))

# Checking presence of clusters in trajectories
for (i in 1:3) {
  traj_cells <- plot_data_long %>% filter(Trajectory == paste0("Trajectory", i))
  print(unique(traj_cells$Cluster))
}

# Check for non-finite values in Pseudotime
non_finite_rows <- plot_data_long %>% filter(!is.finite(Pseudotime))
print(non_finite_rows)

# Normalize pseudotime values within each trajectory
plot_data_long <- plot_data_long %>%
  group_by(Trajectory) %>%
  mutate(Pseudotime = scales::rescale(Pseudotime))

# Filter out rows with non-finite pseudotime values
plot_data_long_filtered <- plot_data_long %>% filter(is.finite(Pseudotime))

# Identify groups with fewer than two datapoints
group_counts <- plot_data_long_filtered %>%
  group_by(Trajectory, Cluster) %>%
  summarise(n = n(), .groups = 'drop')

sparse_groups <- group_counts %>% filter(n < 2)

# Add dummy points for sparse groups
dummy_points <- sparse_groups %>%
  rowwise() %>%
  mutate(Pseudotime = mean(plot_data_long_filtered$Pseudotime, na.rm = TRUE)) %>%
  select(Trajectory, Cluster, Pseudotime)

# Combine original data with dummy points
plot_data_long_adjusted <- bind_rows(plot_data_long_filtered, dummy_points)

# Violin plot with adjusted data
vp <- ggplot(plot_data_long_adjusted, aes(x = Pseudotime, y = Cluster, fill = Cluster)) +
  geom_violin(scale = "width") +
  facet_wrap(~Trajectory, scales = "free_x", drop = FALSE) +
  theme_minimal() +
  labs(title = "Pseudotime Distribution by Cluster and Trajectory",
       x = "Pseudotime",
       y = "Cluster") +
  scale_fill_manual(values = cluster_colors) +  # 사용자 정의 색상 적용
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

print(vp)

ggsave("20240712_pseudotime_vlnplot_by_trajectory.png", plot = vp, width = 9, height = 6, dpi = 1000)

################

library(data.table)
library(Matrix)
library(tradeSeq)
library(splines)

# Step 1: Create SlingshotDataSet object
sds <- SlingshotDataSet(sds)
print(sds)
print(class(sds))

# Step 2: Sample 10 cells from each cluster
set.seed(123)
sampled_cells <- as.data.table(plus@meta.data)[, .(cell = .I[sample(.N, min(.N, 10))]), by = seurat_clusters]
sampled_cells <- sampled_cells$cell

# Step 3: Create a new Seurat object with the sampled cells
plus_sampled <- subset(plus, cells = sampled_cells)

# Step 4: Convert counts to a sparse matrix
counts_shCtrl <- plus_sampled[["RNA"]]$counts.shCtrl
counts_shPTBP1 <- plus_sampled[["RNA"]]$counts.shPTBP1
counts_sampled <- cbind(counts_shCtrl, counts_shPTBP1)
counts_sparse_sampled <- as(counts_sampled, "sparseMatrix")

# Step 5: Create SingleCellExperiment object
umap_coords <- Embeddings(plus_sampled, reduction = "umap")
cell_metadata <- plus_sampled@meta.data
sce_sampled <- SingleCellExperiment(assays = list(counts = counts_sparse_sampled),
                                    colData = cell_metadata,
                                    reducedDims = list(UMAP = umap_coords))

# Step 6: Update cluster labels
update_cluster_labels <- function(label) {
  if (label %in% c("GABAergic neurons", "Glutamatergic neurons")) {
    return("Neurons")
  } else {
    return(label)
  }
}

sce_sampled$updated_customclassif <- sapply(sce_sampled$customclassif, update_cluster_labels)

# Step 7: Run Slingshot with specified start and end clusters
set.seed(123)
sds_sampled <- slingshot(sce_sampled, clusterLabels = 'updated_customclassif', reducedDim = 'UMAP',
                         start.clus = "Fibroblasts", 
                         end.clus = c("Neurons", "Myofibroblasts", "Immature neurons"))
print(slingLineages(sds_sampled))

# Step 8: Convert to SlingshotDataSet
sds <- SlingshotDataSet(sds_sampled)

# Step 9: Fit a generalized additive model (GAM) to the counts
sce_fitted <- fitGAM(counts = assay(sce_sampled, "counts"), 
                     sds = sds,
                     nknots = 5,
                     verbose = TRUE,
                     parallel = FALSE)

# Step 10: Verify the results
print(str(sce_fitted))
print(names(sce_fitted))

metadata(sce_fitted)$slingshot <- sds

print(metadata(sce_fitted)$slingshot)

n_lineages <- length(slingLineages(metadata(sce_fitted)$slingshot))

print(paste("Number of lineages:", n_lineages))

saveRDS(sce_fitted, file = "20240712_fitGAM_results.rds")

