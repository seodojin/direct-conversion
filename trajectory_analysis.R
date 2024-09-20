# 필요한 라이브러리 로드
library(SingleCellExperiment)
library(slingshot)
library(Matrix)
library(tradeSeq)
library(Seurat)
library(clusterExperiment)
library(scater)
library(dplyr)
library(mice)
library(tidyr)
library(ggplot2)
library(ggridges)

# Seurat 객체 불러오기 및 SingleCellExperiment 변환
plus <- readRDS("seurat_object2")

# Seurat 객체에서 meta.data에서 shPTBP1에 해당하는 세포만 선택
shPTBP1_cells <- rownames(plus@meta.data[plus@meta.data$orig.ident == "shPTBP1", ])

# 선택된 shPTBP1 세포들에 대한 RNA 카운트 및 정규화 데이터 추출
counts_shPTBP1 <- plus[["RNA"]]$counts.shPTBP1[, shPTBP1_cells]
data_shPTBP1 <- plus[["RNA"]]$data.shPTBP1[, shPTBP1_cells]

# 필터링된 메타데이터
meta_data_shPTBP1 <- plus@meta.data[shPTBP1_cells, ]

# 필터링된 PCA 및 UMAP 차원 축소 결과
pca_shPTBP1 <- Embeddings(plus, "pca")[shPTBP1_cells, ]
umap_shPTBP1 <- Embeddings(plus, "umap")[shPTBP1_cells, ]

# 필터링된 데이터를 사용하여 SingleCellExperiment 객체 생성
sce_shPTBP1 <- SingleCellExperiment(
  assays = list(counts = counts_shPTBP1, logcounts = data_shPTBP1),
  colData = meta_data_shPTBP1,
  reducedDims = list(PCA = pca_shPTBP1, UMAP = umap_shPTBP1)
)

# 필터링된 SCE 객체 확인
sce_shPTBP1

sce <- sce_shPTBP1

# 신경세포(GABAergic 및 Glutamatergic neurons)를 "Neurons"로 통합하고 나머지는 기존 라벨 유지
sce$updated_customclassif <- sapply(sce$customclassif, function(label) {
  if (label %in% c("GABAergic neurons", "Glutamatergic neurons")) {
    return("Neurons")
  } else {
    return(label)
  }
})

# Slingshot을 사용하여 세포 궤적 분석 수행, 섬유아세포를 시작점으로 설정
sds <- slingshot(sce, clusterLabels = 'updated_customclassif', reducedDim = 'UMAP',
                         start.clus = "Fibroblasts", end.clus = c("Myofibroblasts", "Immature neurons", "Neurons"))

# 유사시간(pseudotime)과 각 리니지에 대한 세포 가중치(cell weights)를 추출
pseudotime <- as.data.frame(slingPseudotime(sds))
cell_weights <- slingCurveWeights(sds)

# Fibroblasts는 무작위로 리니지를 할당하고, 나머지 세포는 가중치가 가장 높은 리니지에 할당
pseudotime$MainLineage <- sapply(seq_len(nrow(pseudotime)), function(i) {
  if (sce$updated_customclassif[i] == "Fibroblasts") {
    return(sample(c("Lineage1", "Lineage2", "Lineage3"), 1))  # Fibroblasts는 무작위로 리니지 할당
  } else {
    lineages <- which(!is.na(pseudotime[i, ]))
    cell_weights_row <- cell_weights[i, lineages]
    main_lineage <- lineages[which.max(cell_weights_row)]
    return(paste0("Lineage", main_lineage))
  }
})

# Pseudotime 정규화 함수 정의
normalize_pseudotime <- function(pseudotime) {
  apply(pseudotime, 2, function(x) {
    valid_idx <- !is.na(x)
    x[valid_idx] <- (x[valid_idx] - min(x[valid_idx])) / (max(x[valid_idx]) - min(x[valid_idx]))
    return(x)
  })
}

# 각 리니지의 유사시간을 0과 1 사이로 정규화
filtered_pseudotime <- pseudotime[, c("Lineage1", "Lineage2", "Lineage3")]
normalized_pseudotime <- normalize_pseudotime(filtered_pseudotime)

# 정규화된 유사시간에 해당하는 세포 이름을 추출
valid_cells_pseudotime <- rownames(normalized_pseudotime)

# Slingshot 결과와 일치하는 세포만 필터링하여 새로운 SCE 객체 생성
sce_filtered_subset <- sce[, valid_cells_pseudotime]
filtered_weights <- cell_weights[valid_cells_pseudotime, , drop = FALSE]
filtered_counts <- assay(sce_filtered_subset, "counts")

# 필터링 후 각 데이터셋의 차원을 확인
cat("Pseudotime dimensions:", dim(normalized_pseudotime), "\n")
cat("CellWeights dimensions:", dim(filtered_weights), "\n")
cat("Counts dimensions:", dim(filtered_counts), "\n")

# 유전자 필터링 기준 설정: 최소 발현 카운트 및 최소 세포 수
min_counts <- 10  # 최소 발현 카운트
min_cells <- 3   # 최소 세포 수

# 희소 행렬로 변환 후 일정 수 이상의 세포에서 일정 수 이상 발현된 유전자만 유지
counts_sparse <- as(filtered_counts, "sparseMatrix")
filtered_genes <- rowSums(counts_sparse >= min_counts) >= min_cells
filtered_counts <- filtered_counts[filtered_genes, ]

# 각 세포에서 발현된 유전자 수를 계산한 후 500개 이상의 유전자가 발현된 세포만 선택
num_genes_per_cell <- colSums(filtered_counts > 0)
cat("500개 이상의 유전자가 발현된 세포 수:", sum(num_genes_per_cell >= 500), "\n")

# 500개 이상의 유전자가 발현된 세포만 유지
valid_cells <- colSums(filtered_counts > 0) >= 500
filtered_counts <- filtered_counts[, valid_cells]
sce_filtered_subset <- sce_filtered_subset[, valid_cells]

# 필터링 후 유전자 및 세포 수 확인
cat("Filtered gene counts dimensions:", dim(filtered_counts), "\n")
cat("Filtered cell counts dimensions:", dim(filtered_counts)[2], "\n")

# 필터링된 유전자와 세포에 맞게 logcounts와 메타데이터도 업데이트
filtered_logcounts <- assay(sce_filtered_subset, "logcounts")[rownames(filtered_counts), ]
sce_filtered <- SingleCellExperiment(
  assays = list(counts = filtered_counts, logcounts = filtered_logcounts),
  colData = colData(sce_filtered_subset),
  reducedDims = reducedDims(sce_filtered_subset)
)

# 필터링된 데이터를 사용하여 Slingshot을 다시 실행하여 새로운 궤적 계산
sds_filtered <- slingshot(
  sce_filtered,
  clusterLabels = 'updated_customclassif',
  reducedDim = 'UMAP',
  start.clus = "Fibroblasts",
  end.clus = c("Myofibroblasts", "Immature neurons", "Neurons"),  # 종점 설정
  approx_points = 100  # 궤적을 매끄럽게 하기 위한 유사시간 경로 설정
)

# Slingshot 재실행 결과에서 유사시간과 세포 가중치를 다시 추출
pseudotime_filtered <- as.data.frame(slingPseudotime(sds_filtered))
cell_weights_filtered <- slingCurveWeights(sds_filtered)

# 각 세포를 가장 높은 가중치를 가진 리니지에 할당하고, 리니지별 세포 수 계산
pseudotime_filtered$MainLineage <- apply(cell_weights_filtered, 1, function(x) {
  if (all(is.na(x))) {
    return(NA)
  } else {
    return(paste0("Lineage", which.max(x)))
  }
})

# 리니지별 세포 수 출력
lineage_counts_filtered <- table(pseudotime_filtered$MainLineage)
cat("리니지별 세포 수 (Slingshot 재실행 후):\n")
print(lineage_counts_filtered)

# Slingshot에서 추가 설정을 사용해 궤적 분석을 다시 실행
sds_new <- slingshot(sce_filtered, clusterLabels = 'updated_customclassif', reducedDim = 'UMAP',
                     start.clus = "Fibroblasts",
                     end.clus = c("Myofibroblasts", "Immature neurons", "Neurons"),
                     allow.breaks = TRUE,  # 궤적에서 단절 허용
                     extend = 'n',         # 궤적 확장하지 않음
                     approx_points = 300)   # 유사시간 경로 매끄럽게 설정

# 재실행된 Slingshot 결과에서 유사시간 및 가중치 추출
new_pseudotime <- slingPseudotime(sds_new)
new_weights <- slingCurveWeights(sds_new)

print("새로운 weights의 NA 개수:")
print(colSums(is.na(new_weights)))

# NA가 아닌 세포에 대해 가중치 기반으로 리니지를 할당하는 함수 정의
assign_lineage <- function(pseudotime, weights) {
  lineage <- rep(NA, nrow(pseudotime))
  for (i in 1:nrow(pseudotime)) {
    valid_lineages <- which(!is.na(pseudotime[i,]))
    if (length(valid_lineages) > 0) {
      lineage[i] <- valid_lineages[which.max(weights[i, valid_lineages])]
    }
  }
  return(lineage)
}

# 리니지 할당
lineage_assignment <- assign_lineage(new_pseudotime, new_weights)

# 리니지 할당 결과 출력
print("리니지 할당 결과:")
print(table(lineage_assignment, useNA = "ifany"))

# 유효한 세포만 선택 (NA가 없는 세포만)
valid_cells <- !is.na(lineage_assignment)
final_pseudotime <- new_pseudotime[valid_cells, ]
final_weights <- new_weights[valid_cells, ]
final_counts <- assay(sce_filtered, "counts")[, valid_cells]

# 최종 데이터 차원 확인
print("최종 데이터 차원:")
print(dim(final_pseudotime))
print(dim(final_weights))
print(dim(final_counts))

###########################################################

# 0과 1 사이로 유사시간을 정규화
normalize_pseudotime <- function(pseudotime) {
  apply(pseudotime, 2, function(x) {
    valid_idx <- !is.na(x)  # NA가 아닌 값들만 정규화
    if (sum(valid_idx) > 0) {
      x[valid_idx] <- (x[valid_idx] - min(x[valid_idx])) / (max(x[valid_idx]) - min(x[valid_idx]))
    }
    return(x)
  })
}

# 유사시간 정규화 적용
normalized_pseudotime <- normalize_pseudotime(final_pseudotime)  # balanced_pseudotime 대신 final_pseudotime 사용

# 정규화된 유사시간 확인
head(normalized_pseudotime)

# 결측값을 처리하기 위해 'mice' 패키지의 다중 대체(Multiple Imputation) 방법 사용
imputed_data <- mice(normalized_pseudotime, method = 'pmm', m = 5, maxit = 50, seed = 123)

# 대체된 데이터셋 중 첫 번째를 선택하여 사용
imputed_pseudotime <- complete(imputed_data, action = 1)

# 대체된 pseudotime 확인
head(imputed_pseudotime)

#######################################################

# 유효한 세포의 인덱스에 해당하는 클러스터 정보 추출
valid_clusters <- sce_filtered$updated_customclassif[valid_cells]

# 각 세포에서 가장 높은 가중치를 가진 리니지에 해당하는 Pseudotime만 남김
pseudotime_long <- data.frame(
  Pseudotime = imputed_pseudotime[cbind(1:nrow(imputed_pseudotime), lineage_assignment)],  # 가장 높은 가중치의 유사시간만 선택
  Lineage = factor(lineage_assignment),  # 각 세포가 할당된 리니지
  CellWeights = final_weights[cbind(1:nrow(final_weights), lineage_assignment)],  # 각 세포의 리니지별 가중치
  Cluster = valid_clusters  # 각 세포가 속한 클러스터 정보
)


# 밀도 기반 플롯 생성 (4x1 배치)

cluster_colors <- c("Immature neurons" = "#00bef3", 
                    "Myofibroblasts" = "#ff8c8c", 
                    "Fibroblasts" = "#19c3a3", 
                    "Neurons" = "#d4a600")

density_plot <- ggplot(pseudotime_long, aes(x = Pseudotime, fill = Cluster, weight = CellWeights)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = cluster_colors) +  
  facet_grid(Cluster ~ ., scales = "free_y") +  # 각 클러스터를 행으로 설정하여 4x1 배치로 변경
  labs(
    title = "Pseudotime Progression Along Differentiation Trajectories",
    x = "Pseudotime",
    y = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# 플롯 확인
print(density_plot)

# 그림을 600dpi로 저장
ggsave("20240920_pseudotime_vs_cluster_weights.png", plot = density_plot, dpi = 600, width = 8, height = 12)

##########################################################

# 대체된 pseudotime을 사용하여 fitGAM 실행
sce_fitted <- fitGAM(
  counts = final_counts,  # balanced_counts 대신 final_counts 사용
  pseudotime = imputed_pseudotime,
  cellWeights = final_weights,  # balanced_weights 대신 final_weights 사용
  nknots = 5,    # 궤적에 대한 매끄러움을 조절하는 매개변수
  verbose = TRUE,
  parallel = FALSE
)

######################################################
# 각 세포를 가장 높은 가중치를 가진 리니지에 할당
lineage_assignment <- apply(final_weights, 1, which.max)  
# 리니지별 세포 수 계산
lineage_counts <- table(lineage_assignment)
# 결과 출력
print("리니지별 세포 수:")
print(lineage_counts)

# fitGAM 결과 저장
saveRDS(sce_fitted, file = "20240910_post_fitGAM.rds")
# Slingshot 객체 저장
saveRDS(sds_filtered, file = "slingshot_obj.rds")
# pseudotime 저장
saveRDS(imputed_pseudotime, file = "pseudotime.rds")
# cell weights 저장
saveRDS(final_weights, file = "cell_weights.rds")  
# counts 객체 저장
saveRDS(final_counts, file = "final_counts.rds")  
