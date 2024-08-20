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
library(igraph)
library(ggraph)

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

# VIPER 분석 실행
viper_results_1_vs_2 <- perform_viper_analysis(significant_genes_1_vs_2)
viper_results_2_vs_3 <- perform_viper_analysis(significant_genes_2_vs_3)

# VIPER Network Visualization
# 결과 시각화 함수
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
  
  png(paste0("correlation_top20_tf_activities_", gsub(" ", "_", title), ".png"), width = 6, height = 6, units = "in", res = 600)
  par(mar = c(4, 4, 6, 2))
  corrplot(cor_matrix, method = "color", type = "upper",
           order = "hclust", tl.col = "black", tl.srt = 45)
  mtext(paste("Correlation of Top 20 TF Activities -", title), side = 3, line = 3, cex = 1.2)
  dev.off()
  
  gc()
  return(top_20_tfs)
}

# 시각화 실행
top_20_tfs_1_vs_2 <- visualize_tf_activities(viper_results_1_vs_2$tf_activities, "Lineage 1 vs 2")
top_20_tfs_2_vs_3 <- visualize_tf_activities(viper_results_2_vs_3$tf_activities, "Lineage 2 vs 3")

# 상위 전사 인자 타겟 분석
analyze_top_tf_targets <- function(top_tfs, regulon, title) {
  for (tf in top_tfs[1:5]) {
    targets <- names(regulon[[tf]]$tfmode)
    cat("Top targets for", tf, "in", title, ":\n")
    print(head(targets, 10))
    cat("\n")
  }
}

# 타겟 분석 실행
analyze_top_tf_targets(top_20_tfs_1_vs_2, viper_results_1_vs_2$regulon, "Lineage 1 vs 2")
analyze_top_tf_targets(top_20_tfs_2_vs_3, viper_results_2_vs_3$regulon, "Lineage 2 vs 3")

# 네트워크 시각화 함수
visualize_tf_target_network <- function(top_tfs, regulon, title, top_n = 3, max_targets_per_tf = 10, file_name) {
  edges <- data.frame()
  
  for (tf in top_tfs[1:top_n]) {
    targets <- names(regulon[[tf]]$tfmode)
    targets <- head(targets, max_targets_per_tf)
    df <- data.frame(from = rep(tf, length(targets)), to = targets)
    edges <- rbind(edges, df)
  }
  
  graph <- graph_from_data_frame(edges)
  
  plot <- ggraph(graph, layout = 'kk') +
    geom_edge_link(aes(edge_alpha = 0.8), show.legend = FALSE) +
    geom_node_point(color = "blue", size = 5) +
    geom_node_text(aes(label = name), vjust = 1, hjust = 1) +
    ggtitle(paste("TF-Target Network -", title)) +
    theme_void()
  
  print(plot)
  ggsave(filename = file_name, plot = plot, width = 10, height = 8, dpi = 300)
}

# 네트워크 시각화 저장
visualize_tf_target_network(top_20_tfs_1_vs_2, viper_results_1_vs_2$regulon, "Lineage 1 vs 2", top_n = 5, max_targets_per_tf = 5, file_name = "lineage_1_vs_2_network.png")
visualize_tf_target_network(top_20_tfs_2_vs_3, viper_results_2_vs_3$regulon, "Lineage 2 vs 3", top_n = 5, max_targets_per_tf = 5, file_name = "lineage_2_vs_3_network.png")




# TRRUST 네트워크 시각화 함수
visualize_trrust_network <- function(tf_genes, title, file_name) {
  graph <- graph_from_data_frame(tf_genes, directed = FALSE)
  
  plot <- ggraph(graph, layout = "kk") +
    geom_edge_link(aes(edge_alpha = 0.8), color = "grey") +
    geom_node_point(size = 5, color = "skyblue") +
    geom_node_text(aes(label = name), vjust = 1.5, hjust = 0.5, size = 4) +
    theme_void() +
    labs(title = title) +
    guides(edge_alpha = "none")
  
  print(plot)
  ggsave(filename = file_name, plot = plot, width = 8, height = 6, dpi = 600)
}

# Lineage 1 vs 2 데이터 및 시각화
tf_genes12 <- data.frame(
  Key_TF = c("RB1", "RB1", "RB1", "MYB", "MYB", "MYB", "RBL1", "RBL1", "NF1", "NF1", "ETS2", "ETS2", "TP53", "TP53", "TP53", "NR3C1", "NR3C1", "HNF4A", "HNF4A", "WT1", "WT1", "AR", "AR", "E2F1", "E2F1", "SP1", "SP1", "SP1", "RELA", "RELA", "NFKB1", "NFKB1"),
  Gene = c("MYC", "ATF2", "TFDP1", "KLF1", "MYC", "SP3", "TFDP1", "MYC", "SP3", "MYC", "MYC", "ERG", "MYC", "SMAD3", "E2F7", "SRF", "ATF2", "TBP", "MYC", "SMAD3", "MYC", "MYC", "TBP", "ETV4", "MYC", "ATF2", "SP3", "MYC", "MYC", "TBP", "TBP", "MYC")
)

visualize_trrust_network(tf_genes12, "Network of Key TFs and Overlapped Genes - Lineage 1 vs 2", "20240820_lineage1vs2_trrust_network.png")

# Lineage 2 vs 3 데이터 및 시각화
tf_genes23 <- data.frame(
  Key_TF = c("RB1", "RB1", "RB1", "MYB", "MYB", "MYB", "CEBPE", "CEBPE", "MYBL2", "MYBL2", "NF1", "NF1", "HDAC2", "HDAC2", "ETS2", "ETS2", "TP53", "TP53", "TP53", "NR3C1", "NR3C1", "HNF4A", "HNF4A", "CEBPA", "CEBPA", "FOS", "FOS", "GATA1", "GATA1", "WT1", "WT1", "ESR1", "ESR1", "AR", "AR", "JUN", "JUN", "SP1", "SP1", "SP1", "RELA", "RELA", "NFKB1", "NFKB1"),
  Gene = c("MYC", "ATF2", "SP1", "KLF1", "MYC", "SP3", "SPI1", "MYC", "SP1", "MYC", "SP3", "MYC", "MYC", "SP1", "MYC", "ERG", "MYC", "SMAD3", "E2F7", "SRF", "ATF2", "TBP", "MYC", "SPI1", "MYC", "SPI1", "MYC", "SPI1", "KLF1", "SMAD3", "MYC", "MYC", "SP1", "MYC", "TBP", "SPI1", "MYC", "ATF2", "SP3", "MYC", "MYC", "TBP", "TBP", "MYC")
)

visualize_trrust_network(tf_genes23, "Network of Key TFs and Overlapped Genes - Lineage 2 vs 3", "20240820_lineage2vs3_trrust_network.png")



##### trrust GO
# 유전자 목록 준비 (예시로 tf_genes12 사용)
genes12 <- unique(c(tf_genes12$Key_TF, tf_genes12$Gene))

# 유전자 이름을 Entrez ID로 변환
entrez_ids12 <- mapIds(org.Hs.eg.db, keys = genes12, keytype = "SYMBOL", column = "ENTREZID")
entrez_ids12 <- entrez_ids12[!is.na(entrez_ids12)]

# GO 분석 수행
go_results12 <- enrichGO(gene = entrez_ids12,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",  # 모든 GO 범주 (BP, CC, MF)
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

# 결과 확인
head(go_results12)

# 결과 시각화
# 상위 10개 GO 용어의 막대 그래프
barplot(go_results12, showCategory = 15)
ggsave("TRRUST_lineage1vs2_GO_barplot.png", width = 12, height = 12, dpi=600)

# GO 네트워크 그래프
# 용어 간 유사성 계산
library(enrichplot)
go_results_sim12 <- pairwise_termsim(go_results12)

# GO 네트워크 그래프
go_network12 <- enrichplot::emapplot(go_results_sim12, showCategory = 20)
ggsave("TRRUST_lineage1vs2_GO_network.png", plot = go_network12, width = 15, height = 12.5)


# 결과를 CSV 파일로 저장
write.csv(as.data.frame(go_results12), "TRRUST_lineage1vs2_GO_analysis_results.csv", row.names = FALSE)



# 유전자 목록 준비 (예시로 tf_genes12 사용)
genes23 <- unique(c(tf_genes23$Key_TF, tf_genes23$Gene))

# 유전자 이름을 Entrez ID로 변환
entrez_ids23 <- mapIds(org.Hs.eg.db, keys = genes23, keytype = "SYMBOL", column = "ENTREZID")
entrez_ids23 <- entrez_ids23[!is.na(entrez_ids23)]

# GO 분석 수행
go_results23 <- enrichGO(gene = entrez_ids23,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",  # 모든 GO 범주 (BP, CC, MF)
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

# 결과 확인
head(go_results23)

# 결과 시각화
# 상위 20개 GO 용어의 막대 그래프
barplot(go_results23, showCategory = 20)
ggsave("TRRUST_lineage2vs3_GO_barplot.png", width = 12, height = 12, dpi=600)

# GO 네트워크 그래프

go_results_sim23 <- pairwise_termsim(go_results23)

# GO 네트워크 그래프
go_network23 <- enrichplot::emapplot(go_results_sim23, showCategory = 20)
ggsave("TRRUST_lineage2vs3_GO_network.png", plot = go_network23, width = 15, height = 12.5)


# 결과를 CSV 파일로 저장
write.csv(as.data.frame(go_results23), "TRRUST_lineage2vs3_GO_analysis_results.csv", row.names = FALSE)
