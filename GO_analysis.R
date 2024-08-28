# Load necessary libraries
library(SingleCellExperiment)
library(slingshot)
library(tradeSeq)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(dorothea)
library(viper)
library(corrplot)
library(igraph)
library(ggraph)
library(ggplot2)
library(pheatmap)
library(VennDiagram)

# 데이터 로드 및 전처리
loaded_results <- readRDS("20240808_fitGAM_results_full.rds")
sce_fitted <- loaded_results
sds <- metadata(sce_fitted)$slingshot
pseudotime <- slingPseudotime(sds)
lineages <- slingLineages(sds)

# Lineage 정보 확인
table(sce_fitted$Lineage)

# slingshot 정보 확인
sds <- metadata(sce_fitted)$slingshot
str(sds)

# 전체 데이터에서 유전자 발현에 대한 association test 수행
asso_test <- associationTest(sce_fitted)

# 유의미한 유전자 식별 (p-value를 FDR 방식으로 조정)
signif_genes <- rownames(asso_test)[which(p.adjust(asso_test$pvalue, method = "fdr") < 0.05)]


# slingshot에서 pseudotime 정보 추출
pseudotime <- slingPseudotime(sds)

# 각 리니지별로 NA가 아닌 pseudotime을 가진 셀들 추출
lin1_cells <- which(!is.na(pseudotime[, 1]))  # 리니지 1에 속하는 셀
lin2_cells <- which(!is.na(pseudotime[, 2]))  # 리니지 2에 속하는 셀
lin3_cells <- which(!is.na(pseudotime[, 3]))  # 리니지 3에 속하는 셀

# 리니지별 유의미한 유전자 추출 (signif_genes는 associationTest 결과에서 나온 유전자들)
lin1_genes <- signif_genes[lin1_cells]
lin2_genes <- signif_genes[lin2_cells]
lin3_genes <- signif_genes[lin3_cells]

# 리니지별 유의미한 유전자 출력
head(lin1_genes)
head(lin2_genes)
head(lin3_genes)


# ENTREZID로 변환 (SYMBOL에서)
lin1_genes_entrez <- bitr(lin1_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
lin2_genes_entrez <- bitr(lin2_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
lin3_genes_entrez <- bitr(lin3_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO Enrichment 분석 (Biological Process)
ego_lin1 <- enrichGO(gene = lin1_genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "fdr")
ego_lin2 <- enrichGO(gene = lin2_genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "fdr")
ego_lin3 <- enrichGO(gene = lin3_genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "fdr")

# 결과 시각화 (예: barplot)
barplot(ego_lin1, showCategory = 10, title = "Lineage 1 GO Enrichment")
barplot(ego_lin2, showCategory = 10, title = "Lineage 2 GO Enrichment")
barplot(ego_lin3, showCategory = 10, title = "Lineage 3 GO Enrichment")



# KEGG Pathway Enrichment 분석
kegg_lin1 <- enrichKEGG(gene = lin1_genes_entrez$ENTREZID, organism = 'hsa', pAdjustMethod = "fdr")
kegg_lin2 <- enrichKEGG(gene = lin2_genes_entrez$ENTREZID, organism = 'hsa', pAdjustMethod = "fdr")
kegg_lin3 <- enrichKEGG(gene = lin3_genes_entrez$ENTREZID, organism = 'hsa', pAdjustMethod = "fdr")

# 결과 시각화
barplot(kegg_lin1, showCategory = 10, title = "Lineage 1 KEGG Pathway Enrichment")
barplot(kegg_lin2, showCategory = 10, title = "Lineage 2 KEGG Pathway Enrichment")
barplot(kegg_lin3, showCategory = 10, title = "Lineage 3 KEGG Pathway Enrichment")




# GO Enrichment barplot을 파일로 저장 (Lineage 1)
p1 <- barplot(ego_lin1, showCategory = 10, title = "Lineage 1 GO Enrichment")
ggsave(filename = "20240827_Lineage1_GO_Enrichment.png", plot = p1, dpi = 600, width = 8, height = 6)

# GO Enrichment barplot을 파일로 저장 (Lineage 2)
p2 <- barplot(ego_lin2, showCategory = 10, title = "Lineage 2 GO Enrichment")
ggsave(filename = "20240827_Lineage2_GO_Enrichment.png", plot = p2, dpi = 600, width = 8, height = 6)

# GO Enrichment barplot을 파일로 저장 (Lineage 3)
p3 <- barplot(ego_lin3, showCategory = 10, title = "Lineage 3 GO Enrichment")
ggsave(filename = "20240827_Lineage3_GO_Enrichment.png", plot = p3, dpi = 600, width = 8, height = 6)

# KEGG Pathway barplot을 파일로 저장 (Lineage 1)
p4 <- dotplot(kegg_lin1, showCategory = 10, title = "Lineage 1 KEGG Pathway Enrichment")
ggsave(filename = "20240827_Lineage1_KEGG_Enrichment.png", plot = p4, dpi = 600, width = 8, height = 6)

# KEGG Pathway barplot을 파일로 저장 (Lineage 2)
p5 <- dotplot(kegg_lin2, showCategory = 10, title = "Lineage 2 KEGG Pathway Enrichment")
ggsave(filename = "20240827_Lineage2_KEGG_Enrichment.png", plot = p5, dpi = 600, width = 8, height = 6)

# KEGG Pathway barplot을 파일로 저장 (Lineage 3)
p6 <- dotplot(kegg_lin3, showCategory = 10, title = "Lineage 3 KEGG Pathway Enrichment")
ggsave(filename = "20240827_Lineage3_KEGG_Enrichment.png", plot = p6, dpi = 600, width = 8, height = 6)






# Venn Diagram 그리기
# NA 값을 제거한 후 Venn Diagram 그리기
venn_data <- list(
  Lineage1 = lin1_genes[!is.na(lin1_genes)],  # NA 제거
  Lineage2 = lin2_genes[!is.na(lin2_genes)],  # NA 제거
  Lineage3 = lin3_genes[!is.na(lin3_genes)]   # NA 제거
)

# Venn Diagram 그리기
venn.plot <- venn.diagram(
  venn_data,
  filename = "20240827_venn_diagram.png", 
  fill = c("#21908c", "#440154", "#fde725"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  cat.col = c("#21908c", "#440154", "#fde725"),
  cat.pos = c(-30, 30, 0),    
  cat.dist = c(0.05, 0.05, 0.05) 
)

# Venn Diagram 출력
grid::grid.draw(venn.plot)



# Dorothea regulon 데이터 로드 및 필터링 (human, high confidence)
regulon_viper <- dorothea::dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B", "C"))  # 신뢰성 높은 regulon만 필터링

# 리니지별로 유의미한 유전자들에 해당하는 전사인자 추출 (리니지 1)
tf_lin1 <- regulon_viper %>%
  dplyr::filter(target %in% lin1_genes) %>%
  dplyr::pull(tf) %>%
  unique()

# 리니지별로 유의미한 유전자들에 해당하는 전사인자 추출 (리니지 2)
tf_lin2 <- regulon_viper %>%
  dplyr::filter(target %in% lin2_genes) %>%
  dplyr::pull(tf) %>%
  unique()

# 리니지별로 유의미한 유전자들에 해당하는 전사인자 추출 (리니지 3)
tf_lin3 <- regulon_viper %>%
  dplyr::filter(target %in% lin3_genes) %>%
  dplyr::pull(tf) %>%
  unique()

# 결과 출력
cat("Lineage 1 TFs: ", tf_lin1, "\n")
cat("Lineage 2 TFs: ", tf_lin2, "\n")
cat("Lineage 3 TFs: ", tf_lin3, "\n")

# 리니지별 전사인자 목록을 파일로 저장
write.csv(tf_lin1, file = "TF_list_Lineage1.csv", row.names = FALSE)
write.csv(tf_lin2, file = "TF_list_Lineage2.csv", row.names = FALSE)
write.csv(tf_lin3, file = "TF_list_Lineage3.csv", row.names = FALSE)


# 교집합
# 리니지 1과 2에서만 겹치는 전사인자
common_tfs_12 <- intersect(tf_lin1, tf_lin2)
cat("Lineage 1과 2에서 겹치는 전사인자: ", common_tfs_12, "\n")

# 리니지 1과 3에서만 겹치는 전사인자
common_tfs_13 <- intersect(tf_lin1, tf_lin3)
cat("Lineage 1과 3에서 겹치는 전사인자: ", common_tfs_13, "\n")

# 리니지 2와 3에서만 겹치는 전사인자
common_tfs_23 <- intersect(tf_lin2, tf_lin3)
cat("Lineage 2와 3에서 겹치는 전사인자: ", common_tfs_23, "\n")

# 차집합
# 리니지 1에만 있는 전사인자 (리니지 2와 3에 속하지 않음)
unique_tfs_lin1 <- setdiff(tf_lin1, union(tf_lin2, tf_lin3))
cat("Lineage 1에만 있는 고유 전사인자: ", unique_tfs_lin1, "\n")

# 리니지 2에만 있는 전사인자 (리니지 1과 3에 속하지 않음)
unique_tfs_lin2 <- setdiff(tf_lin2, union(tf_lin1, tf_lin3))
cat("Lineage 2에만 있는 고유 전사인자: ", unique_tfs_lin2, "\n")

# 리니지 3에만 있는 전사인자 (리니지 1과 2에 속하지 않음)
unique_tfs_lin3 <- setdiff(tf_lin3, union(tf_lin1, tf_lin2))
cat("Lineage 3에만 있는 고유 전사인자: ", unique_tfs_lin3, "\n")






# TRRUST에서 얻은 key regulator 데이터를 로드
key_regulator_data3 <- read.csv("20240827_lineage3.tsv", sep = "\t")

# 네트워크 그래프 생성
# 데이터 프레임에서 key TF와 관련된 타겟 유전자 연결 정보 추출
edges3 <- key_regulator_data3 %>%
  mutate(genes = strsplit(List.of.overlapped.genes, ",")) %>%
  unnest(genes) %>%
  dplyr::select(Key.TF, genes)

# 5개 이상의 타겟을 가진 TF만 필터링하여 시각화
filtered_edges3 <- edges3 %>%
  group_by(Key.TF) %>%
  filter(n() > 5)  # 5개 이상의 타겟을 가진 TF

# 필터링된 네트워크 생성
filtered_net3 <- graph_from_data_frame(filtered_edges3, directed = TRUE)

# 각 노드의 Degree 계산
node_degree3 <- degree(filtered_net3)

# Degree 값에 따라 필터링 (예: Degree가 2 이상인 노드만 시각화)
high_degree_nodes3 <- V(filtered_net3)[node_degree3 >= 2]

# 필터링된 서브 네트워크 생성 (Degree 2 이상인 노드만 포함)
filtered_net_degree3 <- induced_subgraph(filtered_net3, high_degree_nodes3)

# Degree에 따라 노드 색상 설정: Degree가 높을수록 빨간색, 낮으면 파란색
V(filtered_net_degree3)$color <- ifelse(degree(filtered_net_degree3) > 3, "High Degree", "Low Degree")

# Degree에 따라 노드 크기 조정
V(filtered_net_degree3)$size <- degree(filtered_net_degree3) * 3

# 네트워크 시각화 (Degree 기반 노드 색상 및 크기 조정)
ggraph(filtered_net_degree3, layout = "kk") +  # Kamada-Kawai 레이아웃 사용
  geom_edge_link(aes(edge_alpha = 0.3), edge_width = 0.5, show.legend = FALSE) +  # Edge Alpha 레전드 제거
  geom_node_point(aes(size = degree(filtered_net_degree3), color = V(filtered_net_degree3)$color)) +  # 노드 크기와 색상 설정
  geom_node_text(aes(label = name), size = 3, vjust = 1.5, repel = TRUE) +  # 레이블 설정
  scale_color_manual(values = c("High Degree" = "red", "Low Degree" = "lightblue")) +  # 색상 범례 추가
  theme_void() +
  ggtitle("Filtered Key Regulator Network (Lineage 3, Degree-based)") +
  # 범례 추가: 노드 크기와 색상
  guides(size = guide_legend("Node Degree"), color = guide_legend("Node Color"))

# 네트워크 그래프를 파일로 저장 (PNG 예시)
ggsave("20240827_asso_lin3_tfkr.png", width = 10, height = 8, dpi = 600)






# TRRUST에서 얻은 key regulator 데이터를 로드
key_regulator_data2 <- read.csv("20240827_lineage2.tsv", sep = "\t")

# 네트워크 그래프 생성
# 데이터 프레임에서 key TF와 관련된 타겟 유전자 연결 정보 추출
edges2 <- key_regulator_data2 %>%
  mutate(genes = strsplit(List.of.overlapped.genes, ",")) %>%
  unnest(genes) %>%
  dplyr::select(Key.TF, genes)

# 5개 이상의 타겟을 가진 TF만 필터링하여 시각화
filtered_edges2 <- edges2 %>%
  group_by(Key.TF) %>%
  filter(n() > 5)  # 5개 이상의 타겟을 가진 TF

# 필터링된 네트워크 생성
filtered_net2 <- graph_from_data_frame(filtered_edges2, directed = TRUE)

# 각 노드의 Degree 계산
node_degree2 <- degree(filtered_net2)

# Degree 값에 따라 필터링 (예: Degree가 2 이상인 노드만 시각화)
high_degree_nodes2 <- V(filtered_net2)[node_degree2 >= 2]

# 필터링된 서브 네트워크 생성 (Degree 2 이상인 노드만 포함)
filtered_net_degree2 <- induced_subgraph(filtered_net2, high_degree_nodes2)

# Degree에 따라 노드 색상 설정: Degree가 높을수록 빨간색, 낮으면 파란색
V(filtered_net_degree2)$color <- ifelse(degree(filtered_net_degree2) > 3, "High Degree", "Low Degree")

# Degree에 따라 노드 크기 조정
V(filtered_net_degree2)$size <- degree(filtered_net_degree2) * 3

# 네트워크 시각화 (Degree 기반 노드 색상 및 크기 조정)
ggraph(filtered_net_degree2, layout = "kk") +  # Kamada-Kawai 레이아웃 사용
  geom_edge_link(aes(edge_alpha = 0.3), edge_width = 0.5, show.legend = FALSE) +  # Edge Alpha 레전드 제거
  geom_node_point(aes(size = degree(filtered_net_degree2), color = V(filtered_net_degree2)$color)) +  # 노드 크기와 색상 설정
  geom_node_text(aes(label = name), size = 3, vjust = 1.5, repel = TRUE) +  # 레이블 설정
  scale_color_manual(values = c("High Degree" = "red", "Low Degree" = "lightblue")) +  # 색상 범례 추가
  theme_void() +
  ggtitle("Filtered Key Regulator Network (Lineage 2, Degree-based)") +
  # 범례 추가: 노드 크기와 색상
  guides(size = guide_legend("Node Degree"), color = guide_legend("Node Color"))

# 네트워크 그래프를 파일로 저장 (PNG 예시)
ggsave("20240827_asso_lin2_tfkr.png", width = 10, height = 8, dpi = 600)






# TRRUST에서 얻은 key regulator 데이터를 로드
key_regulator_data1 <- read.csv("20240827_lineage1.tsv", sep = "\t")

# 네트워크 그래프 생성
# 데이터 프레임에서 key TF와 관련된 타겟 유전자 연결 정보 추출
edges1 <- key_regulator_data1 %>%
  mutate(genes = strsplit(List.of.overlapped.genes, ",")) %>%
  unnest(genes) %>%
  dplyr::select(Key.TF, genes)

# 5개 이상의 타겟을 가진 TF만 필터링하여 시각화
filtered_edges1 <- edges1 %>%
  group_by(Key.TF) %>%
  filter(n() > 5)  # 5개 이상의 타겟을 가진 TF

# 필터링된 네트워크 생성
filtered_net1 <- graph_from_data_frame(filtered_edges1, directed = TRUE)

# 각 노드의 Degree 계산
node_degree1 <- degree(filtered_net1)

# Degree 값에 따라 필터링 (예: Degree가 2 이상인 노드만 시각화)
high_degree_nodes1 <- V(filtered_net1)[node_degree1 >= 2]

# 필터링된 서브 네트워크 생성 (Degree 2 이상인 노드만 포함)
filtered_net_degree1 <- induced_subgraph(filtered_net1, high_degree_nodes1)

# Degree에 따라 노드 색상 설정: Degree가 높을수록 빨간색, 낮으면 파란색
V(filtered_net_degree1)$color <- ifelse(degree(filtered_net_degree1) > 3, "High Degree", "Low Degree")

# Degree에 따라 노드 크기 조정
V(filtered_net_degree1)$size <- degree(filtered_net_degree1) * 3

# 네트워크 시각화 (Degree 기반 노드 색상 및 크기 조정)
ggraph(filtered_net_degree1, layout = "kk") +  # Kamada-Kawai 레이아웃 사용
  geom_edge_link(aes(edge_alpha = 0.3), edge_width = 0.5, show.legend = FALSE) +  # Edge Alpha 레전드 제거
  geom_node_point(aes(size = degree(filtered_net_degree1), color = V(filtered_net_degree1)$color)) +  # 노드 크기와 색상 설정
  geom_node_text(aes(label = name), size = 3, vjust = 1.5, repel = TRUE) +  # 레이블 설정
  scale_color_manual(values = c("High Degree" = "red", "Low Degree" = "lightblue")) +  # 색상 범례 추가
  theme_void() +
  ggtitle("Filtered Key Regulator Network (Lineage 1, Degree-based)") +
  # 범례 추가: 노드 크기와 색상
  guides(size = guide_legend("Node Degree"), color = guide_legend("Node Color"))

# 네트워크 그래프를 파일로 저장 (PNG 예시)
ggsave("20240827_asso_lin1_tfkr.png", width = 10, height = 8, dpi = 600)





# 함수: 리니지별 허브 노드 식별 및 GO 분석
analyze_hub_nodes_and_go <- function(key_regulator_data, lineage_name) {
  
  # 네트워크 그래프 생성
  edges <- key_regulator_data %>%
    mutate(genes = strsplit(List.of.overlapped.genes, ",")) %>%
    unnest(genes) %>%
    dplyr::select(Key.TF, genes)
  
  # 5개 이상의 타겟을 가진 TF만 필터링하여 시각화
  filtered_edges <- edges %>%
    group_by(Key.TF) %>%
    filter(n() > 5)
  
  # 필터링된 네트워크 생성
  filtered_net <- graph_from_data_frame(filtered_edges, directed = TRUE)
  
  # 각 노드의 Degree 계산
  node_degree <- degree(filtered_net)
  
  # Degree 값에 따라 필터링 (예: Degree가 2 이상인 노드만 시각화)
  high_degree_nodes <- V(filtered_net)[node_degree >= 2]
  
  # 필터링된 서브 네트워크 생성 (Degree 2 이상인 노드만 포함)
  filtered_net_degree <- induced_subgraph(filtered_net, high_degree_nodes)
  
  # 허브 노드 식별 (상위 5% Degree 노드)
  hub_nodes <- V(filtered_net_degree)[degree(filtered_net_degree) > quantile(degree(filtered_net_degree), 0.95)]
  
  # 허브 노드의 이름을 추출
  hub_node_names <- names(hub_nodes)
  
  # 허브 노드에 대한 GO 분석 수행
  ego <- enrichGO(gene = hub_node_names, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "fdr")
  
  # GO 분석 결과 시각화
  barplot(ego, showCategory = 10, title = paste("GO Enrichment for Hub Nodes (", lineage_name, ")", sep = ""))
  
  # 허브 노드와 GO 결과 반환
  return(list(hub_nodes = hub_node_names, go_results = ego))
}

# 리니지 1 분석
key_regulator_data1 <- read.csv("20240827_lineage1.tsv", sep = "\t")
result1 <- analyze_hub_nodes_and_go(key_regulator_data1, "Lineage 1")

# 리니지 2 분석
key_regulator_data2 <- read.csv("20240827_lineage2.tsv", sep = "\t")
result2 <- analyze_hub_nodes_and_go(key_regulator_data2, "Lineage 2")

# 리니지 3 분석
key_regulator_data3 <- read.csv("20240827_lineage3.tsv", sep = "\t")
result3 <- analyze_hub_nodes_and_go(key_regulator_data3, "Lineage 3")


# 리니지 1에 대한 GO 분석 결과 시각화
barplot(result1$go_results, showCategory = 10, title = "GO Enrichment for Hub Nodes (Lineage 1)")

# PNG 파일로 저장
ggsave("20240827_asso_lineage1_GO.png", width = 10, height = 8, dpi = 600)

# 리니지 2에 대한 GO 분석 결과 시각화
barplot(result2$go_results, showCategory = 10, title = "GO Enrichment for Hub Nodes (Lineage 2)")

# PNG 파일로 저장
ggsave("20240827_asso_lineage2_GO.png", width = 10, height = 8, dpi = 600)

# 리니지 3에 대한 GO 분석 결과 시각화
barplot(result3$go_results, showCategory = 10, title = "GO Enrichment for Hub Nodes (Lineage 3)")

# PNG 파일로 저장
ggsave("20240827_asso_lineage3_GO.png", width = 10, height = 8, dpi = 600)
