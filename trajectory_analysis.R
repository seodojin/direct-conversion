library(SingleCellExperiment)
library(slingshot)
library(scater)
library(Seurat)
library(Rcpp)

sessionInfo()

# Convert Seurat object to SingleCellExperiment object

plus = readRDS("seurat_object2")
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
library(RColorBrewer)
cell_colors <- brewer.pal(4, "Paired")

# 궤적 색상 정의
trajectory_colors <- c("red", "green", "blue")

# ggplot으로 그리기
p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(values = cell_colors) +
  theme_minimal() +
  labs(title = "Cell Types with Trajectories") +
  theme(
    legend.text = element_text(size = 14),  # 레전드 텍스트 크기
    legend.title = element_text(size = 14)  # 레전드 제목 크기
  )

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

ggsave("20240719_trajectory_umap_plot.png", plot = p, width = 7, height = 6, dpi = 1000)


########## 20240708 pseudotime vln plot
# load packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Assign each cell to the primary trajectory based on the lowest pseudotime value
pseudotime <- as.data.frame(slingPseudotime(sds))
plot_data <- data.frame(
  Cell = rownames(pseudotime),
  Cluster = sce$updated_customclassif
)

# Get the primary trajectory for each cell
plot_data$PrimaryTrajectory <- apply(pseudotime, 1, function(x) {
  traj <- which.min(x)
  if (is.infinite(x[traj])) {
    return(NA)
  } else {
    return(paste0("Trajectory", traj))
  }
})

# Merge pseudotime data with primary trajectory information
plot_data <- cbind(plot_data, pseudotime)

# Reshape data for plotting
plot_data_long <- plot_data %>%
  pivot_longer(cols = starts_with("Lineage"), 
               names_to = "Trajectory", 
               values_to = "Pseudotime") %>%
  filter(!is.na(Pseudotime))

# Trajectory 열의 값을 "Trajectory1", "Trajectory2", "Trajectory3"로 변경
plot_data_long$Trajectory <- paste0("Trajectory", substr(plot_data_long$Trajectory, 8, 8))

# Filter to keep only the primary trajectory for each cell
plot_data_long <- plot_data_long %>%
  filter(Trajectory == PrimaryTrajectory)

# 결과 확인
print(summary(plot_data_long$Pseudotime))
print(summary(plot_data_long$Trajectory))
print(head(plot_data_long))

# Normalize pseudotime values within each trajectory
plot_data_long_scaled <- plot_data_long %>%
  group_by(Trajectory) %>%
  mutate(Pseudotime = scales::rescale(Pseudotime)) %>%
  ungroup()

# Define cluster colors
cluster_colors <- brewer.pal(length(unique(plot_data_long_scaled$Cluster)), "Paired")



# 1. 각 궤적의 end cluster 지정
trajectory_endpoints <- data.frame(
  Trajectory = c("Trajectory1", "Trajectory2", "Trajectory3"),
  End_Cluster = c("Myofibroblasts", "Immature neurons", "Neurons")
)

# 2. 데이터 필터링
plot_data_filtered <- plot_data_long_scaled %>%
  left_join(trajectory_endpoints, by = "Trajectory") %>%
  filter(Cluster == "Fibroblasts" | Cluster == End_Cluster)

# 3. 궤적별 바이올린 플롯
vp <- ggplot(plot_data_filtered, aes(x = Pseudotime, y = Cluster, fill = Cluster)) +
  geom_violin(scale = "width") +
  facet_wrap(~Trajectory, scales = "free_x", drop = TRUE) +
  theme_minimal() +
  labs(title = "Pseudotime Distribution by Cluster and Trajectory",
       x = "Pseudotime",
       y = "Cluster") +
  scale_fill_manual(values = cluster_colors) +
  theme(
    legend.position = "none",  # 레전드 제거
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )
vp
# 4. 전체 유사시간에 따른 클러스터 분포 시각화
overall_plot <- ggplot(plot_data_long_scaled, aes(x = Pseudotime, y = Cluster, color = Cluster)) +
  geom_jitter(alpha = 0.5, height = 0.2) +
  theme_minimal() +
  labs(title = "Overall Pseudotime Distribution by Cluster",
       x = "Pseudotime",
       y = "Cluster") +
  scale_color_manual(values = cluster_colors) +
  theme(
    legend.position = "none",  # 레전드 제거
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )
overall_plot
# 플롯 저장
ggsave("20240719_pseudotime_distribution_by_cluster_and_trajectory.png", 
       plot = vp, 
       width = 10, height = 6, dpi = 1000)

ggsave("20240719_overall_pseudotime_distribution_by_cluster.png", 
       plot = overall_plot, 
       width = 10, height = 6, dpi = 1000)

# Save the SingleCellExperiment object
saveRDS(sce, file = "20240719_sce")


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

