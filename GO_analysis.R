# 필요한 라이브러리 로드
library(SingleCellExperiment)
library(slingshot)
library(tradeSeq)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# 데이터 로드 및 Slingshot 데이터 추출
sce_fitted <- readRDS("20240808_fitGAM_results_full.rds")
sds <- metadata(sce_fitted)$slingshot
pseudotime <- slingPseudotime(sds)

### 리니지 1과 2에 대한 분석 ###

# 리니지 1과 2의 셀 선택 및 Fitted values 예측
cells_lineage_1_2 <- which(!is.na(pseudotime[, "Lineage1"]) & !is.na(pseudotime[, "Lineage2"]))
sce_fitted_1_2 <- sce_fitted[, cells_lineage_1_2]
pattern_results_1_2 <- patternTest(sce_fitted_1_2)
pattern_results_1_2$FDR <- p.adjust(pattern_results_1_2$pvalue, method = "fdr")
significant_genes_1_2 <- pattern_results_1_2 %>%
  filter(FDR < 0.05)
cat("리니지 1과 2에서 유의미한 유전자의 수:", nrow(significant_genes_1_2), "\n")
fitted_values_1_2 <- predictSmooth(sce_fitted_1_2, gene = intersect(rownames(significant_genes_1_2), rownames(sce_fitted)), nPoints = 100)

# 데이터 정리 및 유사시간 기반 클러스터링 준비
plot_data_1_2 <- fitted_values_1_2 %>%
  filter(lineage %in% c(1, 2)) %>%
  group_by(gene) %>%
  mutate(Expression = scale(yhat, center = TRUE, scale = TRUE)[,1]) %>%
  ungroup() %>%
  dplyr::select(gene, Pseudotime = time, lineage, Expression)

# 유사시간 기반 클러스터링 함수 정의 및 적용
perform_clustering_pseudotime <- function(data, n_clusters = 5) {
  set.seed(123)
  wide_data <- data %>%
    group_by(gene, Pseudotime) %>%
    summarize(Expression = mean(Expression, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = Pseudotime, values_from = Expression) %>%
    mutate(across(-gene, ~ifelse(is.na(.) | is.nan(.) | is.infinite(.), mean(., na.rm = TRUE), .)))
  kmeans_result <- kmeans(wide_data[,-1], centers = n_clusters)
  tibble(gene = wide_data$gene, cluster = kmeans_result$cluster)
}

clusters_pseudotime_1_2 <- map(1:2, ~plot_data_1_2 %>%
                                 filter(lineage == .x) %>%
                                 perform_clustering_pseudotime()) %>%
  set_names(paste0("cluster", 1:2))

# 클러스터링 결과 병합
plot_data_with_clusters_1_2 <- plot_data_1_2 %>%
  left_join(clusters_pseudotime_1_2$cluster1, by = "gene") %>%
  left_join(clusters_pseudotime_1_2$cluster2, by = "gene") %>%
  mutate(cluster = case_when(
    lineage == 1 ~ cluster.x,
    lineage == 2 ~ cluster.y
  )) %>%
  dplyr::select(-cluster.x, -cluster.y)

# 리니지 1과 2에서 서로 다른 클러스터에 속하는 유전자 식별
different_clusters_1_2 <- plot_data_with_clusters_1_2 %>%
  group_by(gene, lineage) %>%
  summarize(cluster = first(cluster), .groups = 'drop') %>%
  pivot_wider(names_from = lineage, values_from = cluster, names_prefix = "Lineage_") %>%
  filter(Lineage_1 != Lineage_2)

cat("리니지 1과 2에서 서로 다른 클러스터에 속하는 유전자의 수:", nrow(different_clusters_1_2), "\n")

# GO 분석 수행 및 시각화
entrez_ids_1_2 <- bitr(different_clusters_1_2$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
go_results_1_2 <- enrichGO(gene = entrez_ids_1_2$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", 
                           pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
emap_plot_1_2 <- emapplot(pairwise_termsim(go_results_1_2), showCategory = 20, layout.params = list(layout = "kk")) +
  ggtitle("Enrichment Map of GO Enrichment Analysis (BP) for Lineage 1 and 2")
ggsave("20240902_Enrichment_Map_Lineage1_vs_2.png", plot = emap_plot_1_2, width = 10, height = 8, dpi = 600)

# 유사시간을 0에서 1로 정규화 (리니지 1과 2만 사용)
plot_data_1_2 <- plot_data_with_clusters_1_2 %>%
  dplyr::filter(lineage %in% c(1, 2)) %>%
  dplyr::group_by(lineage) %>%
  dplyr::mutate(Pseudotime = (Pseudotime - min(Pseudotime)) / (max(Pseudotime) - min(Pseudotime))) %>%
  ungroup()

# 클러스터 별로 유전자 그룹화 (예: 클러스터 1~5)
cluster1_genes <- plot_data_with_clusters_1_2 %>% dplyr::filter(cluster == 1) %>% dplyr::pull(gene)
cluster2_genes <- plot_data_with_clusters_1_2 %>% dplyr::filter(cluster == 2) %>% dplyr::pull(gene)
cluster3_genes <- plot_data_with_clusters_1_2 %>% dplyr::filter(cluster == 3) %>% dplyr::pull(gene)
cluster4_genes <- plot_data_with_clusters_1_2 %>% dplyr::filter(cluster == 4) %>% dplyr::pull(gene)
cluster5_genes <- plot_data_with_clusters_1_2 %>% dplyr::filter(cluster == 5) %>% dplyr::pull(gene)

# 시각화 및 저장 함수 정의
plot_and_save_cluster <- function(cluster_genes, cluster_number) {
  plot_data_cluster <- plot_data_1_2 %>%
    dplyr::filter(gene %in% cluster_genes)
  
  p <- ggplot(plot_data_cluster, aes(x = Pseudotime, y = Expression, color = factor(lineage))) +
    geom_line(aes(group = interaction(gene, lineage)), alpha = 0.7) +
    scale_color_manual(values = c("1" = "#21908c", "2" = "#440154")) +
    theme_minimal() +
    labs(title = paste("Cluster", cluster_number, "Expression Patterns for Lineage 1 and Lineage 2"),
         x = "Pseudotime",
         y = "Standardized log(count + 1)",
         color = "Lineage") +
    theme(legend.position = "bottom")
  
  # 그래프 저장
  ggsave(filename = paste0("Cluster", cluster_number, "_Lineage1_vs_Lineage2_normalized.png"), 
         plot = p, width = 7, height = 5, dpi = 600)
}

# 각 클러스터별 시각화 및 저장
plot_and_save_cluster(cluster1_genes, 1)
plot_and_save_cluster(cluster2_genes, 2)
plot_and_save_cluster(cluster3_genes, 3)
plot_and_save_cluster(cluster4_genes, 4)
plot_and_save_cluster(cluster5_genes, 5)


### 리니지 2와 3에 대한 분석 ###

# 리니지 2와 3의 셀 선택 및 Fitted values 예측
cells_lineage_2_3 <- which(!is.na(pseudotime[, "Lineage2"]) & !is.na(pseudotime[, "Lineage3"]))
sce_fitted_2_3 <- sce_fitted[, cells_lineage_2_3]
pattern_results_2_3 <- patternTest(sce_fitted_2_3)
pattern_results_2_3$FDR <- p.adjust(pattern_results_2_3$pvalue, method = "fdr")
significant_genes_2_3 <- pattern_results_2_3 %>%
  filter(FDR < 0.05)
cat("리니지 2와 3에서 유의미한 유전자의 수:", nrow(significant_genes_2_3), "\n")

fitted_values_2_3 <- predictSmooth(sce_fitted_2_3, gene = intersect(rownames(significant_genes_2_3), rownames(sce_fitted)), nPoints = 100)

# 데이터 정리 및 유사시간 기반 클러스터링 준비
plot_data_2_3 <- fitted_values_2_3 %>%
  filter(lineage %in% c(2, 3)) %>%
  group_by(gene) %>%
  mutate(Expression = scale(yhat, center = TRUE, scale = TRUE)[,1]) %>%
  ungroup() %>%
  dplyr::select(gene, Pseudotime = time, lineage, Expression)

# 클러스터링 수행 및 결과 병합
clusters_pseudotime_2_3 <- map(2:3, ~plot_data_2_3 %>%
                                 filter(lineage == .x) %>%
                                 perform_clustering_pseudotime()) %>%
  set_names(paste0("cluster", 2:3))
plot_data_with_clusters_2_3 <- plot_data_2_3 %>%
  left_join(clusters_pseudotime_2_3$cluster2, by = "gene") %>%
  left_join(clusters_pseudotime_2_3$cluster3, by = "gene") %>%
  mutate(cluster = case_when(
    lineage == 2 ~ cluster.x,
    lineage == 3 ~ cluster.y
  )) %>%
  dplyr::select(-cluster.x, -cluster.y)

# 리니지 2와 3에서 서로 다른 클러스터에 속하는 유전자 식별
different_clusters_2_3 <- plot_data_with_clusters_2_3 %>%
  group_by(gene, lineage) %>%
  summarize(cluster = first(cluster), .groups = 'drop') %>%
  pivot_wider(names_from = lineage, values_from = cluster, names_prefix = "Lineage_") %>%
  filter(Lineage_2 != Lineage_3)
cat("리니지 2와 3에서 서로 다른 클러스터에 속하는 유전자의 수:", nrow(different_clusters_2_3), "\n")

# GO 분석 수행 및 시각화
entrez_ids_2_3 <- bitr(different_clusters_2_3$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
go_results_2_3 <- enrichGO(gene = entrez_ids_2_3$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", 
                           pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
emap_plot_2_3 <- emapplot(pairwise_termsim(go_results_2_3), showCategory = 20, layout.params = list(layout = "kk")) +
  ggtitle("Enrichment Map of GO Enrichment Analysis (BP) for Lineage 2 and 3")
ggsave("20240902_Enrichment_Map_Lineage2_vs_3.png", plot = emap_plot_2_3, width = 10, height = 8, dpi = 600)

# 유사시간을 0에서 1로 정규화 (리니지 2와 3만 사용)
plot_data_2_3 <- plot_data_with_clusters_2_3 %>%
  dplyr::filter(lineage %in% c(2, 3)) %>%
  dplyr::group_by(lineage) %>%
  dplyr::mutate(Pseudotime = (Pseudotime - min(Pseudotime)) / (max(Pseudotime) - min(Pseudotime))) %>%
  ungroup()

# 클러스터 별로 유전자 그룹화 (예: 클러스터 1~5)
cluster1_genes <- plot_data_with_clusters_2_3 %>% dplyr::filter(cluster == 1) %>% dplyr::pull(gene)
cluster2_genes <- plot_data_with_clusters_2_3 %>% dplyr::filter(cluster == 2) %>% dplyr::pull(gene)
cluster3_genes <- plot_data_with_clusters_2_3 %>% dplyr::filter(cluster == 3) %>% dplyr::pull(gene)
cluster4_genes <- plot_data_with_clusters_2_3 %>% dplyr::filter(cluster == 4) %>% dplyr::pull(gene)
cluster5_genes <- plot_data_with_clusters_2_3 %>% dplyr::filter(cluster == 5) %>% dplyr::pull(gene)

# 각 클러스터별 시각화 및 저장
plot_and_save_cluster <- function(cluster_genes, cluster_number) {
  plot_data_cluster <- plot_data_2_3 %>%
    dplyr::filter(gene %in% cluster_genes)
  
  p <- ggplot(plot_data_cluster, aes(x = Pseudotime, y = Expression, color = factor(lineage))) +
    geom_line(aes(group = interaction(gene, lineage)), alpha = 0.7) +
    scale_color_manual(values = c("2" = "#440154", "3" = "#fde725")) +
    theme_minimal() +
    labs(title = paste("Cluster", cluster_number, "Expression Patterns for Lineage 2 and Lineage 3"),
         x = "Pseudotime",
         y = "Standardized log(count + 1)",
         color = "Lineage") +
    theme(legend.position = "bottom")
  
  # 그래프 저장
  ggsave(filename = paste0("Cluster", cluster_number, "_Lineage2_vs_Lineage3_normalized.png"), 
         plot = p, width = 7, height = 5, dpi = 600)
}

# 각 클러스터별 시각화 및 저장
plot_and_save_cluster(cluster1_genes, 1)
plot_and_save_cluster(cluster2_genes, 2)
plot_and_save_cluster(cluster3_genes, 3)
plot_and_save_cluster(cluster4_genes, 4)
plot_and_save_cluster(cluster5_genes, 5)
