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


loaded_data <- readRDS("data/20240802_fitGAM_results010.rds")
sce_fitted <- loaded_data$sce_fitted
metadata(sce_fitted)$slingshot <- loaded_data$slingshot_data

# slingshot 관련 정보 확인
print(metadata(sce_fitted)$slingshot)

# lineage 정보 확인
print(slingLineages(metadata(sce_fitted)$slingshot))

# 객체의 메타데이터 확인
print(metadata(sce_fitted))

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



# 각 궤적별로 클러스터 유전자 목록 추출
extract_cluster_genes <- function(data, lineage_num) {
  data %>%
    filter(lineage == lineage_num) %>%
    group_by(cluster) %>%
    summarise(genes = list(unique(gene))) %>%
    mutate(cluster_name = paste0("lineage", lineage_num, "_cluster", cluster))
}

lineage1_clusters <- extract_cluster_genes(plot_data_with_clusters, 1)
lineage2_clusters <- extract_cluster_genes(plot_data_with_clusters, 2)
lineage3_clusters <- extract_cluster_genes(plot_data_with_clusters, 3)

# 모든 클러스터 정보 합치기
all_clusters <- bind_rows(lineage1_clusters, lineage2_clusters, lineage3_clusters)




# GO 분석 함수 (이전과 동일)
perform_go_analysis <- function(gene_list, cluster_name) {
  ego <- enrichGO(gene = gene_list,
                  OrgDb = org.Hs.eg.db,
                  keyType = "SYMBOL",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2)
  
  return(list(result = ego, name = cluster_name))
}

# 모든 클러스터에 대해 GO 분석 수행
go_results <- lapply(1:nrow(all_clusters), function(i) {
  perform_go_analysis(all_clusters$genes[[i]], all_clusters$cluster_name[i])
})

# emapplot을 사용한 결과 시각화
for (result in go_results) {
  # 결과가 비어있지 않은 경우에만 그래프 생성
  if (nrow(result$result) > 0) {
    # 엣지 계산
    ego_sim <- pairwise_termsim(result$result)
    
    # emapplot 생성
    p <- emapplot(ego_sim, showCategory = 20)  # 상위 30개 카테고리 표시
    
    # 그래프 출력
    print(p + ggtitle(paste("Enrichment Map -", result$name)))
    
    # 그래프 저장
    ggsave(
      filename = paste0("emapplot_", result$name, ".png"),
      plot = p,
      width = 12,  # 인치 단위, 필요에 따라 조정
      height = 8, # 인치 단위, 필요에 따라 조정
      dpi = 600,
      bg = "white" # 배경색 설정
    )
  } else {
    cat("No significant GO terms found for", result$name, "\n")
  }
}




# 특정 클러스터의 시간에 따른 발현 변화 시각화
plot_time_course <- function(data, lineage_num, cluster_num) {
  data %>%
    filter(lineage == lineage_num, cluster == cluster_num) %>%
    ggplot(aes(x = Pseudotime, y = Expression, group = gene, color = gene)) +
    geom_line() +
    theme_minimal() +
    labs(title = paste("Expression over time - Lineage", lineage_num, "Cluster", cluster_num),
         x = "Pseudotime", y = "Expression") +
    theme(legend.position = "none")
}

# 예: 궤적 3의 클러스터 5 시각화
print(plot_time_course(plot_data_with_clusters, 3, 5))
print(plot_time_course(plot_data_with_clusters, 3, 3))
print(plot_time_course(plot_data_with_clusters, 1, 5))
print(plot_time_course(plot_data_with_clusters, 1, 2))
print(plot_time_course(plot_data_with_clusters, 2, 4))



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
