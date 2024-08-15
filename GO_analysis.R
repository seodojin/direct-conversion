# 필요한 라이브러리 로드
library(SingleCellExperiment)
library(slingshot)
library(tradeSeq)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dorothea)
library(viper)
library(corrplot)

# 데이터 로드 및 전처리
loaded_results <- readRDS("20240808_fitGAM_results_full.rds")
sce_fitted <- loaded_results
sds <- metadata(sce_fitted)$slingshot
pseudotime <- slingPseudotime(sds)
lineages <- slingLineages(sds)

# 리니지 간 패턴 비교 함수
compare_lineages <- function(lineage1, lineage2) {
  selected_lineages <- c(lineage1, lineage2)
  pseudotime_subset <- pseudotime[, selected_lineages]
  valid_cells <- rownames(pseudotime_subset)[apply(!is.na(pseudotime_subset), 1, all)]
  sce_fitted_subset <- sce_fitted[, valid_cells]
  pseudotime_subset <- pseudotime_subset[valid_cells, ]
  metadata(sce_fitted_subset)$slingshot <- list(pseudotime = pseudotime_subset)
  
  pattern_results <- patternTest(sce_fitted_subset)
  pattern_results$FDR <- p.adjust(pattern_results$pvalue, method = "fdr")
  significant_genes <- rownames(pattern_results[pattern_results$FDR < 0.005, ])
  
  return(significant_genes)
}

# 리니지 1 vs 2 및 2 vs 3 비교
significant_genes_1_vs_2 <- compare_lineages("Lineage1", "Lineage2")
significant_genes_2_vs_3 <- compare_lineages("Lineage2", "Lineage3")

# GO 분석
perform_go_analysis <- function(gene_list) {
  ego <- enrichGO(gene = gene_list,
                  OrgDb = org.Hs.eg.db,
                  keyType = "SYMBOL",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
  return(ego)
}

ego_1_vs_2 <- perform_go_analysis(significant_genes_1_vs_2)
ego_2_vs_3 <- perform_go_analysis(significant_genes_2_vs_3)

# Barplot for ego_1_vs_2
barplot(ego_1_vs_2, showCategory=20)  # 상위 20개 용어 보여주기

# Barplot for ego_2_vs_3
barplot(ego_2_vs_3, showCategory=20)  # 상위 20개 용어 보여주기



# VIPER 분석 함수
perform_viper_analysis <- function(gene_list) {
  common_genes <- intersect(rownames(sce_fitted), gene_list)
  expr_data <- assay(sce_fitted, "counts")[common_genes, ]
  
  dorothea_regulon_human <- dorothea_hs %>% dplyr::filter(confidence %in% c("A", "B", "C"))
  regulon <- dorothea::df2regulon(dorothea_regulon_human)
  
  common_genes <- intersect(rownames(expr_data), unique(unlist(lapply(regulon, function(x) names(x$tfmode)))))
  expr_data_filtered <- expr_data[common_genes, ]
  
  tf_activities <- viper::viper(expr_data_filtered, regulon,
                                method = "none", minsize = 4,
                                eset.filter = FALSE, verbose = TRUE)
  
  return(list(tf_activities = tf_activities, regulon = regulon))
}

# Lineage 1 vs 2 VIPER 분석
viper_results_1_vs_2 <- perform_viper_analysis(significant_genes_1_vs_2)
tf_activities_1_vs_2 <- viper_results_1_vs_2$tf_activities
regulon_1_vs_2 <- viper_results_1_vs_2$regulon

# Lineage 2 vs 3 VIPER 분석
viper_results_2_vs_3 <- perform_viper_analysis(significant_genes_2_vs_3)
tf_activities_2_vs_3 <- viper_results_2_vs_3$tf_activities
regulon_2_vs_3 <- viper_results_2_vs_3$regulon

# 결과 시각화 함수 수정
visualize_tf_activities <- function(tf_activities, title, sample_size = 1000) {
  set.seed(123)
  sample_cols <- sample(1:ncol(tf_activities), min(sample_size, ncol(tf_activities)))
  tf_activities_sampled <- tf_activities[, sample_cols]
  
  top_tfs <- sort(rowMeans(tf_activities_sampled), decreasing = TRUE)
  top_20_tfs <- names(head(top_tfs, 20))
  
  tf_activities_mean <- rowMeans(tf_activities_sampled)
  df <- data.frame(TF = names(tf_activities_mean), Activity = tf_activities_mean)
  
  p1 <- ggplot(df, aes(x = reorder(TF, -Activity), y = Activity)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = paste("Mean TF Activities -", title), x = "Transcription Factors", y = "Mean Activity")
  ggsave(paste0("mean_tf_activities_", gsub(" ", "_", title), ".png"), p1, width = 9, height = 6, dpi = 600)
  
  png(paste0("heatmap_top20_tf_activities_", gsub(" ", "_", title), ".png"), width = 6, height = 5, units = "in", res = 600)
  pheatmap(tf_activities_sampled[top_20_tfs,], show_colnames = FALSE,
           main = paste("Heatmap of Top 20 TF Activities -", title))
  dev.off()
  
  cor_matrix <- cor(t(tf_activities_sampled[top_20_tfs,]))
  
  # 여백을 늘려 제목이 잘리지 않게 수정
  png(paste0("correlation_top20_tf_activities_", gsub(" ", "_", title), ".png"), width = 6, height = 6, units = "in", res = 600)
  par(mar = c(4, 4, 6, 2))  # 상단 여백을 6으로 설정하여 제목 공간을 늘림
  corrplot(cor_matrix, method = "color", type = "upper",
           order = "hclust", tl.col = "black", tl.srt = 45)
  
  # mtext로 제목 추가
  mtext(paste("Correlation of Top 20 TF Activities -", title), side = 3, line = 3, cex = 1.2)
  
  dev.off()
  
  gc()  # 메모리 정리
  return(top_20_tfs)
}

# Lineage 1 vs 2 결과 시각화
top_20_tfs_1_vs_2 <- visualize_tf_activities(tf_activities_1_vs_2, "Lineage 1 vs 2")

# Lineage 2 vs 3 결과 시각화
top_20_tfs_2_vs_3 <- visualize_tf_activities(tf_activities_2_vs_3, "Lineage 2 vs 3")

# 상위 활성 전사 인자들의 타겟 유전자 분석
analyze_top_tf_targets <- function(top_tfs, regulon, title) {
  for (tf in top_tfs[1:5]) {  # 상위 5개 전사 인자에 대해 분석
    targets <- names(regulon[[tf]]$tfmode)
    cat("Top targets for", tf, "in", title, ":\n")
    print(head(targets, 10))
    cat("\n")
  }
}

# Lineage 1 vs 2의 상위 전사 인자 타겟 분석
analyze_top_tf_targets(top_20_tfs_1_vs_2, regulon_1_vs_2, "Lineage 1 vs 2")

# Lineage 2 vs 3의 상위 전사 인자 타겟 분석
analyze_top_tf_targets(top_20_tfs_2_vs_3, regulon_2_vs_3, "Lineage 2 vs 3")

