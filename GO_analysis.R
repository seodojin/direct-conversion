# 저장된 객체 불러오기
sce_fitted <- readRDS("20241114_sce_fitted.rds")

# associationTest 실행하여 유의미한 유전자 식별
assoc_test_results <- associationTest(sce_fitted, lineages = 1)  # lineage1을 대상으로 설정

# 유의미한 유전자들 필터링 (예: p값 기준)
alpha <- 0.05  # 유의 수준 설정
sig_genes <- rownames(assoc_test_results)[assoc_test_results$pvalue < alpha]

cat("p value 기준으로 선택된 유의미한 유전자 수:", length(sig_genes), "\n")

# 'sce_fitted'와 'sig_genes' 사이의 공통 유전자 확인
common_genes <- intersect(rownames(sce_fitted), sig_genes)

# 공통 유전자들만 추출
counts_data <- assay(sce_fitted, "counts")[common_genes, ]
significant_counts <- log1p(counts_data)  # log 변환

cat("유의미한 공통 유전자 수:", length(common_genes), "\n")

# 유의미한 공통 유전자에 대한 counts 데이터 추출 및 로그 변환
all_counts_data <- assay(sce_fitted, "counts")[common_genes, ]
all_significant_counts <- log1p(all_counts_data)  # log1p는 log(x + 1)

# 발현 값 분포 확인
summary(as.vector(all_significant_counts))

# 색상 범위 제한 설정 (-2에서 2 사이로 설정, 필요시 조정 가능)
breaks <- seq(-2, 2, length.out = 101)

lineage1_clusters <- sce_filtered$updated_customclassif[colnames(sce_filtered) %in% lineage1_cells]
print(table(lineage1_clusters))

lineage1_clusters <- sce_filtered$updated_customclassif[match(lineage1_cells, colnames(sce_filtered))]
print(table(lineage1_clusters))

lineage1_clusters <- colData(sce_filtered)$updated_customclassif[match(lineage1_cells, colnames(sce_filtered))]
print(table(lineage1_clusters))

# Neurons와 Fibroblasts 클러스터에 속하는 세포만 선택
selected_clusters <- c("Neurons", "Fibroblasts")
selected_cells <- lineage1_cells[lineage1_clusters %in% selected_clusters]

# 선택된 세포에 대한 발현 데이터 추출
selected_counts <- all_significant_counts[, selected_cells]

# 가장 변동성이 큰 상위 100개 유전자 선택
gene_variances <- apply(selected_counts, 1, var)
top_variable_genes <- names(sort(gene_variances, decreasing = TRUE))[1:100]

# 선택된 유전자와 세포로 데이터 준비
plot_data <- selected_counts[top_variable_genes, ]

# 클러스터 정보 준비
cluster_info <- lineage1_clusters[match(selected_cells, lineage1_cells)]

# 주석 데이터프레임 생성
annotation_col <- data.frame(
  Cluster = cluster_info,
  row.names = colnames(plot_data)
)

# 클러스터별 색상 지정
cluster_colors <- c("Fibroblasts" = "yellow", "Neurons" = "green")

# 히트맵 생성 및 저장
png("20241113_lineage1_neurons_fibroblasts_heatmap.png", width = 7.5, height = 6, units = "in", res = 600)
pheatmap(
  plot_data,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Heatmap of Top Variable Genes in Lineage 1 (Neurons and Fibroblasts)",
  breaks = breaks,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  annotation_col = annotation_col,
  annotation_colors = list(Cluster = cluster_colors)
)
dev.off()

##########################

neuron_cells <- colnames(plot_data)[cluster_info == "Neurons"]

neuron_expression <- rowMeans(plot_data[, neuron_cells])

non_neuron_cells <- colnames(plot_data)[cluster_info != "Neurons"]
non_neuron_expression <- rowMeans(plot_data[, non_neuron_cells])

expression_difference <- neuron_expression - non_neuron_expression

top_neuron_genes <- names(sort(expression_difference, decreasing = TRUE)[1:20])  # 상위 20개 선택

# 파일 저장을 위한 png 디바이스 열기
png("20241114_neuron_specific_genes_heatmap.png", width = 6, height = 5, units = "in", res = 600)

# 히트맵 생성
pheatmap(plot_data[top_neuron_genes, ],
         scale = "row",
         annotation_col = annotation_col,
         main = "Top 20 Neuron-specific Genes",
         show_rownames = TRUE,  # 유전자 이름 표시
         show_colnames = FALSE,  # 세포 이름 숨기기
         fontsize_row = 8,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         annotation_colors = list(Cluster = cluster_colors))

# 디바이스 닫기
dev.off()

###################################################

# 1. earlyDETest 실행
early_de_results <- earlyDETest(sce_fitted)

# 2. 유의미한 유전자 추출 (예: p-value < 0.05)
significant_genes <- rownames(early_de_results)[early_de_results$pvalue < 0.05]

# 3. dorothea를 사용하여 전사인자 식별
data(dorothea_hs, package = "dorothea") # 인간 데이터 사용. 마우스의 경우 dorothea_mm 사용
regulons <- dorothea_hs %>% 
  dplyr::filter(confidence %in% c("A", "B", "C"))

# 유의미한 유전자 중 전사인자 식별
tf_in_sig_genes <- intersect(significant_genes, unique(regulons$tf))

# 결과 출력
print("유의미한 전사인자:")
print(tf_in_sig_genes)

# 전사인자 수 출력
print(paste("유의미한 전사인자 수:", length(tf_in_sig_genes)))

# 전사인자 목록
tf_list <- c("ATF4", "CEBPB", "CEBPD", "EGR1", "FOSL1", "HIF1A", "JUN", "JUND", 
             "KLF6", "NFIC", "NR2F2", "PBX3", "RBPJ", "STAT1")

# counts 데이터 추출 (SingleCellExperiment 객체에서)
counts_data <- counts(sce_fitted)

# 각 전사인자에 대해 개별적으로 그래프 생성 및 저장
for (tf in tf_list) {
  tryCatch({
    p <- plotSmoothers(sce_fitted, counts = counts_data, gene = tf, border = TRUE) +
      ggtitle(paste("Expression of", tf)) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    # 그래프 저장
    ggsave(paste0("TF_expression_", tf, ".png"), plot = p, width = 4, height = 3, dpi = 600)
    
    print(paste("Saved plot for", tf))
  }, error = function(e) {
    print(paste("Error plotting", tf, ":", e$message))
  })
}


################################# 이하 차등발현분석
library(limma)

# 디자인 매트릭스 생성
design <- model.matrix(~ imputed_pseudotime$Lineage1 + 
                         imputed_pseudotime$Lineage2 + 
                         imputed_pseudotime$Lineage3)

# limma를 사용한 차등 발현 분석
fit <- lmFit(tf_expression, design)
fit <- eBayes(fit)

# 각 리니지에 대한 결과 추출
lineage1_results <- topTable(fit, coef = "imputed_pseudotime$Lineage1", number = Inf)
lineage2_results <- topTable(fit, coef = "imputed_pseudotime$Lineage2", number = Inf)
lineage3_results <- topTable(fit, coef = "imputed_pseudotime$Lineage3", number = Inf)

# 결과 병합
de_results <- data.frame(
  TF = rownames(lineage1_results),
  Lineage1_logFC = lineage1_results$logFC,
  Lineage2_logFC = lineage2_results$logFC,
  Lineage3_logFC = lineage3_results$logFC
)

# 결과 시각화
de_results_long <- pivot_longer(de_results, cols = c("Lineage1_logFC", "Lineage2_logFC", "Lineage3_logFC"), 
                                names_to = "Lineage", values_to = "logFC")
de_results_long$Lineage <- gsub("_logFC", "", de_results_long$Lineage)

# 리니지 이름 변경
lineage_names <- c("Lineage1" = "Neurons", 
                   "Lineage2" = "Myofibroblasts", 
                   "Lineage3" = "Immature neurons")

# 결과 데이터프레임의 리니지 이름 변경
de_results_long$Lineage <- lineage_names[de_results_long$Lineage]

# 수정된 차등 발현 분석 그래프 생성
de_plot <- ggplot(de_results_long, aes(x = TF, y = logFC, fill = Lineage)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Differential Expression of TFs across Lineages", 
       x = "Transcription Factor", y = "log2 Fold Change") +
  scale_fill_manual(values = c("Neurons" = "blue", 
                               "Myofibroblasts" = "green", 
                               "Immature neurons" = "red"))

# 그래프를 고해상도로 저장
ggsave("differential_expression_hires.png", de_plot, width = 9, height = 6, dpi = 600)
