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


loaded_data <- readRDS("20240802_fitGAM_results010.rds")
sce_fitted <- loaded_data$sce_fitted
metadata(sce_fitted)$slingshot <- loaded_data$slingshot_data

# slingshot 관련 정보 확인
print(metadata(sce_fitted)$slingshot)

# lineage 정보 확인
print(slingLineages(metadata(sce_fitted)$slingshot))

# lineage 수 확인
n_lineages <- length(slingLineages(metadata(sce_fitted)$slingshot))
print(paste("Number of lineages:", n_lineages))



# 객체의 메타데이터 확인
print(metadata(sce_fitted))

# slingshot 데이터 확인
print(metadata(sce_fitted)$slingshot)

# assay 데이터 구조 확인
print(assayNames(sce_fitted))
print(dim(assay(sce_fitted, "counts")))

# 세포의 그룹 정보 확인
print(colData(sce_fitted)$lineage)




slingshot_data <- metadata(sce_fitted)$slingshot
print(slingLineages(slingshot_data))





# SlingshotDataSet 객체 확인 
sds <- metadata(sce_fitted)$slingshot

# 의사 시간 추출
pseudotime_matrix <- slingPseudotime(sds)

# 각 세포에 가장 높은 의사 시간을 가진 궤적을 할당
cell_lineages <- apply(pseudotime_matrix, 1, function(x) {
  if (all(is.na(x))) {
    return(NA)
  } else {
    return(names(which.max(x)))
  }
})

# 세포별 궤적 정보 추가
colData(sce_fitted)$lineage <- cell_lineages

# 의사시간 정보 추가 (가장 높은 값을 가진 궤적의 의사시간 사용)
colData(sce_fitted)$pseudotime <- apply(pseudotime_matrix, 1, max, na.rm = TRUE)

# 세포별 궤적 정보 확인
print(table(colData(sce_fitted)$lineage))

# 의사시간 정보 확인
summary(colData(sce_fitted)$pseudotime)






#### patternTest 수행
pattern_results <- patternTest(sce_fitted)

# p-value 조정 (FDR)
pattern_results$FDR <- p.adjust(pattern_results$pvalue, method = "fdr")

# 유의미한 유전자 선택 (FDR < 0.005)
significant_genes <- pattern_results[pattern_results$FDR < 0.005, ]

# 유의미한 유전자 이름 추출
significant_genes_names <- rownames(significant_genes)

# 유의미한 유전자 수 확인
cat("Significant genes across all lineages:", length(significant_genes_names), "\n")

# 적합값 추출 및 표준화
fitted_values <- predictSmooth(sce_fitted, gene = significant_genes_names, nPoints = 100)

# 유전자별로 yhat 표준화
standardized_fitted <- fitted_values %>%
  group_by(gene) %>%
  mutate(standardized_yhat = scale(yhat, center = TRUE, scale = TRUE)) %>%
  ungroup()

# 시각화를 위한 데이터 준비
plot_data <- standardized_fitted %>%
  select(gene, time, lineage, standardized_yhat) %>%
  rename(Pseudotime = time, Expression = standardized_yhat)

# 유전자 발현 패턴 시각화
ggplot(plot_data, aes(x = Pseudotime, y = Expression, color = factor(lineage))) +
  geom_line(aes(group = gene), alpha = 0.3) +
  facet_wrap(~lineage, scales = "free_y") +
  theme_minimal() +
  labs(title = "Gene Expression Patterns Across Lineages",
       x = "Pseudotime",
       y = "Standardized Expression") +
  scale_color_manual(values = c("1" = "#21908c", "2" = "#440154", "3" = "#fde725")) +
  theme(legend.position = "none")


# 각 궤적별로 클러스터링 수행
perform_clustering <- function(data, n_clusters = 5) {
  # 데이터 준비
  wide_data <- data %>%
    pivot_wider(names_from = Pseudotime, values_from = Expression) %>%
    select(-lineage)
  
  # 클러스터링
  set.seed(123)  # 재현성을 위해
  kmeans_result <- kmeans(wide_data[,-1], centers = n_clusters)
  
  # 결과 반환
  wide_data$cluster <- kmeans_result$cluster
  return(wide_data %>% select(gene, cluster))
}

# 각 궤적별로 클러스터링 수행
lineage1_clusters <- plot_data %>% filter(lineage == 1) %>% perform_clustering()
lineage2_clusters <- plot_data %>% filter(lineage == 2) %>% perform_clustering()
lineage3_clusters <- plot_data %>% filter(lineage == 3) %>% perform_clustering()

# 클러스터 정보를 원본 데이터에 추가
plot_data_with_clusters <- plot_data %>%
  left_join(lineage1_clusters, by = "gene") %>%
  rename(cluster1 = cluster) %>%
  left_join(lineage2_clusters, by = "gene") %>%
  rename(cluster2 = cluster) %>%
  left_join(lineage3_clusters, by = "gene") %>%
  rename(cluster3 = cluster) %>%
  mutate(cluster = case_when(
    lineage == 1 ~ cluster1,
    lineage == 2 ~ cluster2,
    lineage == 3 ~ cluster3
  ))

# 색상 팔레트 생성
n_colors <- 5  # 클러스터 수
color_palette <- colorRampPalette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))(n_colors)

# 유전자 발현 패턴 시각화
ggplot(plot_data_with_clusters, aes(x = Pseudotime, y = Expression, color = factor(cluster))) +
  geom_line(aes(group = gene), alpha = 0.5) +
  facet_wrap(~lineage, scales = "free_y", ncol = 3) +
  theme_minimal() +
  labs(title = "Gene Expression Patterns Across Lineages",
       x = "Pseudotime",
       y = "Standardized Expression",
       color = "Cluster") +
  scale_color_manual(values = color_palette) +
  theme(legend.position = "bottom")

# 그래프를 파일로 저장
ggsave("20240806_gene_expression_patterns_clustered.png", width = 7.5, height = 4, dpi = 1000)













### 1. earlyDETest 수행
early_de_results <- earlyDETest(sce_fitted)

# 2. p-value 조정 (FDR 계산)
early_de_results <- early_de_results %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  mutate(fdr = p.adjust(pvalue, method = "fdr"))

# 3. 유의미한 유전자 선택 (FDR < 0.005)
significant_genes <- early_de_results %>%
  dplyr::filter(fdr < 0.005) %>%
  dplyr::arrange(fdr) %>%
  dplyr::slice_head(n = 20)

# 4. 발현 데이터 추출
expression_data <- counts(sce_fitted)[significant_genes$gene, ]

# 5. 데이터 정규화 (log2 변환 및 스케일링)
expression_data_normalized <- t(scale(t(log2(expression_data + 1))))

# 6. 히트맵 생성 (색 범위 조정)
color_palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

# 데이터의 95% 범위 계산
quantiles <- quantile(expression_data_normalized, probs = c(0.025, 0.975))
color_breaks <- seq(quantiles[1], quantiles[2], length.out = 101)

# 히트맵 생성 및 저장
png("20240806_heatmap_high_res.png", width = 10, height = 8, units = "in", res = 600)

pheatmap(
  expression_data_normalized,
  color = color_palette,
  breaks = color_breaks,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  main = "Early Differential Expression Heatmap",
  fontsize_row = 12,
  annotation_col = data.frame(
    Lineage = factor(colData(sce_fitted)$lineage)
  ),
  annotation_colors = list(
    Lineage = c("Lineage1" = "#21908c", "Lineage2" = "#440154", "Lineage3" = "#fde725")
  ),
  annotation_legend = TRUE,
  legend_breaks = seq(-3, 3, by = 1),
  legend_labels = seq(-3, 3, by = 1),
  labels_row = significant_genes$gene  
)

dev.off()
