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

# Seurat 객체 불러오기 및 SingleCellExperiment 변환
plus <- readRDS("seurat_object2")

# Seurat 객체에서 RNA 카운트 및 정규화 데이터를 추출하여 두 조건을 결합
counts <- cbind(plus[["RNA"]]$counts.shCtrl, plus[["RNA"]]$counts.shPTBP1)
data <- cbind(plus[["RNA"]]$data.shCtrl, plus[["RNA"]]$data.shPTBP1)

# Seurat 객체에서 추출한 카운트, 정규화 데이터, 메타데이터 및 차원 축소 결과를 사용하여 SCE 객체 생성
sce <- SingleCellExperiment(
  assays = list(counts = counts, logcounts = data),
  colData = plus@meta.data,
  reducedDims = list(PCA = Embeddings(plus, "pca"), UMAP = Embeddings(plus, "umap"))
)

# 신경세포(GABAergic 및 Glutamatergic neurons)를 "Neurons"로 통합하고 나머지는 기존 라벨 유지
sce$updated_customclassif <- sapply(sce$customclassif, function(label) {
  if (label %in% c("GABAergic neurons", "Glutamatergic neurons")) {
    return("Neurons")
  } else {
    return(label)
  }
})

# Fibroblasts 세포 추출 및 무작위로 100개 샘플링
# (분석의 초점을 신경세포 및 다른 세포들의 궤적과 발달 경로에 맞추기 위한 것)
fibroblast_cells <- which(sce$updated_customclassif == "Fibroblasts")
non_fibroblast_cells <- which(sce$updated_customclassif != "Fibroblasts")
set.seed(123)  # 샘플링의 재현성을 위해 시드 설정
sampled_fibroblasts <- sample(fibroblast_cells, size = 100)

# 샘플링된 섬유아세포와 나머지 세포를 병합하여 최종 데이터셋 생성
sce_reduced <- sce[, c(sampled_fibroblasts, non_fibroblast_cells)]

# Slingshot을 사용하여 세포 궤적 분석 수행, 섬유아세포를 시작점으로 설정
sds_reduced <- slingshot(sce_reduced, clusterLabels = 'updated_customclassif', reducedDim = 'UMAP',
                         start.clus = "Fibroblasts", end.clus = c("Myofibroblasts", "Immature neurons", "Neurons"))

# 유사시간(pseudotime)과 각 리니지에 대한 세포 가중치(cell weights)를 추출
pseudotime_reduced <- as.data.frame(slingPseudotime(sds_reduced))
cell_weights_reduced <- slingCurveWeights(sds_reduced)

# Fibroblasts는 무작위로 리니지를 할당하고, 나머지 세포는 가중치가 가장 높은 리니지에 할당
pseudotime_reduced$MainLineage <- sapply(seq_len(nrow(pseudotime_reduced)), function(i) {
  if (sce_reduced$updated_customclassif[i] == "Fibroblasts") {
    return(sample(c("Lineage1", "Lineage2", "Lineage3"), 1))  # Fibroblasts는 무작위로 리니지 할당
  } else {
    lineages <- which(!is.na(pseudotime_reduced[i, ]))
    cell_weights_row <- cell_weights_reduced[i, lineages]
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
filtered_pseudotime <- pseudotime_reduced[, c("Lineage1", "Lineage2", "Lineage3")]
normalized_pseudotime <- normalize_pseudotime(filtered_pseudotime)

# 정규화된 유사시간에 해당하는 세포 이름을 추출
valid_cells_pseudotime <- rownames(normalized_pseudotime)

# Slingshot 결과와 일치하는 세포만 필터링하여 새로운 SCE 객체 생성
sce_filtered_subset <- sce_reduced[, valid_cells_pseudotime]
filtered_weights <- cell_weights_reduced[valid_cells_pseudotime, , drop = FALSE]
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

# NA 값 확인
print("새로운 pseudotime의 NA 개수:")
print(colSums(is.na(new_pseudotime)))

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

######################################################

# 리니지별 최소 세포 수를 기반으로 각 리니지에서 균형 잡힌 샘플링 진행
min_cells <- min(table(lineage_assignment[valid_cells]))
set.seed(42)  # 샘플링의 재현성을 위해 시드 설정
balanced_cells <- unlist(lapply(unique(lineage_assignment[valid_cells]), function(i) {
  sample(which(lineage_assignment[valid_cells] == i), min_cells)
}))

# 균형 잡힌 데이터셋 생성
balanced_pseudotime <- final_pseudotime[balanced_cells, ]
balanced_weights <- final_weights[balanced_cells, ]
balanced_counts <- final_counts[, balanced_cells]

# 균형 잡힌 데이터셋의 리니지별 세포 수 출력
print("균형 잡힌 데이터셋의 리니지별 세포 수:")
print(table(lineage_assignment[valid_cells][balanced_cells]))

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
normalized_pseudotime <- normalize_pseudotime(balanced_pseudotime)

# 정규화된 유사시간 확인
head(normalized_pseudotime)

# 결측값을 처리하기 위해 'mice' 패키지의 다중 대체(Multiple Imputation) 방법 사용
imputed_data <- mice(normalized_pseudotime, method = 'pmm', m = 5, maxit = 50, seed = 123)

# 대체된 데이터셋 중 첫 번째를 선택하여 사용
imputed_pseudotime <- complete(imputed_data, action = 1)

# 대체된 pseudotime 확인
head(imputed_pseudotime)

# 대체된 pseudotime을 사용하여 fitGAM 실행
sce_fitted <- fitGAM(
  counts = balanced_counts,
  pseudotime = imputed_pseudotime,
  cellWeights = balanced_weights,
  nknots = 5,    # 궤적에 대한 매끄러움을 조절하는 매개변수
  verbose = TRUE,
  parallel = FALSE
)


######################################################
# 각 세포를 가장 높은 가중치를 가진 리니지에 할당
lineage_assignment <- apply(balanced_weights, 1, which.max)

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
saveRDS(balanced_weights, file = "cell_weights.rds")
