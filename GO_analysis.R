# Load necessary libraries
library(SingleCellExperiment)
library(slingshot)
library(tradeSeq)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(cluster)
library(RColorBrewer)

# Load the saved object
loaded_results <- readRDS("20240808_fitGAM_results_full.rds")

# Access the fitted SingleCellExperiment object
sce_fitted <- loaded_results

# Extract Slingshot data from the metadata
sds <- metadata(sce_fitted)$slingshot

# Access pseudotime data
pseudotime <- slingPseudotime(sds)
lineages <- slingLineages(sds)

# Print summary of pseudotime
print(summary(pseudotime))

# Extract pseudotime data
pseudotime_matrix <- slingPseudotime(sds)


# 리니지 1과 2 간의 패턴 비교
selected_lineages <- c("Lineage1", "Lineage2")
pseudotime_subset <- pseudotime_matrix[, selected_lineages]

# NA가 아닌 셀들만 선택
valid_cells <- rownames(pseudotime_subset)[apply(!is.na(pseudotime_subset), 1, all)]

# 새로운 SingleCellExperiment 객체 생성 (유효한 셀만 포함)
sce_fitted_subset <- sce_fitted[, valid_cells]

# 의사시간 데이터를 sce_fitted_subset에 맞게 줄임
pseudotime_subset <- pseudotime_subset[valid_cells, ]

# slingshot 메타데이터 업데이트
metadata(sce_fitted_subset)$slingshot <- list(pseudotime = pseudotime_subset)

# 패턴 테스트 수행
pattern_results_subset <- patternTest(sce_fitted_subset)

# p-value 조정 (FDR)
pattern_results_subset$FDR <- p.adjust(pattern_results_subset$pvalue, method = "fdr")

# 유의한 유전자 선택 (FDR < 0.005)
significant_genes_subset <- pattern_results_subset[pattern_results_subset$FDR < 0.005, ]
significant_genes_names_subset <- rownames(significant_genes_subset)

cat("Number of significant DE genes between Lineage 1 and Lineage 2:", length(significant_genes_names_subset), "\n")


# 리니지 1과 2 간의 유의미한 유전자 목록
significant_genes_1_vs_2 <- significant_genes_names_subset

# 리니지 1과 2 간의 유의미한 유전자 목록을 CSV 파일로 저장
write.csv(significant_genes_1_vs_2, 
          file = "Lineage1_vs_Lineage2_significant_genes.csv", 
          row.names = FALSE, quote = FALSE)


# 리니지 2와 3 간의 패턴 비교
selected_lineages <- c("Lineage2", "Lineage3")
pseudotime_subset <- pseudotime_matrix[, selected_lineages]

# NA가 아닌 셀들만 선택
valid_cells <- rownames(pseudotime_subset)[apply(!is.na(pseudotime_subset), 1, all)]

# 새로운 SingleCellExperiment 객체 생성 (유효한 셀만 포함)
sce_fitted_subset <- sce_fitted[, valid_cells]

# 의사시간 데이터를 sce_fitted_subset에 맞게 줄임
pseudotime_subset <- pseudotime_subset[valid_cells, ]

# slingshot 메타데이터 업데이트
metadata(sce_fitted_subset)$slingshot <- list(pseudotime = pseudotime_subset)

# 패턴 테스트 수행
pattern_results_subset <- patternTest(sce_fitted_subset)

# p-value 조정 (FDR)
pattern_results_subset$FDR <- p.adjust(pattern_results_subset$pvalue, method = "fdr")

# 유의한 유전자 선택 (FDR < 0.005)
significant_genes_subset <- pattern_results_subset[pattern_results_subset$FDR < 0.005, ]
significant_genes_names_subset <- rownames(significant_genes_subset)

cat("Number of significant DE genes between Lineage 2 and Lineage 3:", length(significant_genes_names_subset), "\n")

significant_genes_2_vs_3 <- significant_genes_names_subset

write.csv(significant_genes_2_vs_3, 
          file = "Lineage2_vs_Lineage3_significant_genes.csv", 
          row.names = FALSE, quote = FALSE)



# 리니지 2와 3의 유의미한 유전자들이 리니지 1과 2에도 포함되는지 확인
common_genes <- intersect(significant_genes_2_vs_3, significant_genes_1_vs_2)

# 교집합의 수 출력
cat("Number of common genes between Lineage 2 vs 3 and Lineage 1 vs 2:", length(common_genes), "\n")

# 공통 유전자 목록 출력
print(common_genes)



# 리니지 1과 2의 유의미한 유전자 중에서 리니지 2와 3에서는 유의미하지 않은 유전자
unique_genes_1_vs_2 <- setdiff(significant_genes_1_vs_2, common_genes)

# 리니지 2와 3의 유의미한 유전자 중에서 리니지 1과 2에서는 유의미하지 않은 유전자
unique_genes_2_vs_3 <- setdiff(significant_genes_2_vs_3, common_genes)

# 결과 출력
cat("Number of unique genes in Lineage 1 vs Lineage 2:", length(unique_genes_1_vs_2), "\n")
cat("Number of unique genes in Lineage 2 vs Lineage 3:", length(unique_genes_2_vs_3), "\n")

# 고유 유전자 목록 출력
print("Unique genes in Lineage 1 vs Lineage 2:")
print(unique_genes_1_vs_2)

print("Unique genes in Lineage 2 vs Lineage 3:")
print(unique_genes_2_vs_3)

