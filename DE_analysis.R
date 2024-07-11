### Perform DE analysis
### reference : https://satijalab.org/seurat/articles/de_vignette.html
## input file : annotation_object.rds
## output file : heatmap.png, neuro.csv, volcanoplot.png
# load package
library(DESeq2)
library(Seurat)
library(SeuratData)
library(EnhancedVolcano)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# load data
plus <- readRDS("seurat_object")

new.cluster.ids <- c("Immature neurons", "Myofibroblasts", "Fibroblasts", "Unknown", 
                     "Fibroblasts", "Glutamatergic neurons",
                     "GABAergic neurons")
names(new.cluster.ids) <- levels(plus)
plus <- RenameIdents(plus, new.cluster.ids)

# find markers for every cluster compared to all remaining cells
plus.markers <- FindAllMarkers(plus, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
plus.markers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC) -> marker.gene.list
View(marker.gene.list)
write.csv(marker.gene.list, file = "marker gene list 50.csv")

# top 5 markers for each cluster
plus.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5

# Extract the data you need from Seurat objects
data <- GetAssayData(plus, assay = "RNA", slot = "data")

genes <- top5$gene

# Extract data matrix for genes of interest
mat <- data[genes, ]

# Convert a sparse matrix to a regular matrix
mat <- as.matrix(mat)

# data scaling
mat <- t(scale(t(log1p(mat))))

# Enable row dendrogram (gene) reordering
row_dend_reorder <- TRUE 
  
cluster_anno<- plus@meta.data$customclassif
quantile(mat, c(0.05, 0.95))
  
col_fun = circlize::colorRamp2(c(-1, 0, 2), c("blue", "white", "red"))

my_colors <- c("Glutamatergic neurons" = "#619cff",
               "GABAergic neurons" = "#f564e3",
               "Immature neurons" = "#f8766d",
               "Myofibroblasts" = "#b79f00",
               "Fibroblasts" = "#00ba38",
               "Unknown" = "#00BFC4")

# HeatmapAnnotation
top_annotation <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = my_colors[levels(factor(cluster_anno))])), 
                                    show_legend = TRUE)

cluster_anno <- factor(cluster_anno, levels = names(my_colors))
top_annotation <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = my_colors[levels(cluster_anno)])))

# heatmap
png("heatmap.png", width = 6.4, height = 4.3, units = "in", res = 1000)

Heatmap(mat, 
        name = "Expression", 
        top_annotation = top_annotation,  # 수정된 top_annotation 사용
        column_split = factor(cluster_anno, levels = names(my_colors)),
        row_dend_reorder = row_dend_reorder, 
        cluster_rows = TRUE, 
        cluster_columns = F,
        show_row_names = TRUE, 
        show_column_names = F,
        cluster_column_slices = F,
        row_names_gp = gpar(fontsize = 8),
        column_title_rot = 30,
        column_title_gp = gpar(fontsize=10),
        col = col_fun,
        show_heatmap_legend = TRUE,
        use_raster = TRUE,
        raster_device = c("png"),
        raster_quality = 10)

dev.off()


# visualizing marker expression
VlnPlot(plus, features = "PTBP1") + 
  NoLegend() + 
  theme(axis.title.x = element_blank(), plot.title = element_blank()) + 
  ylab("PTBP1 Expression Level")
ggsave("20240321ptbp1violinplot.png", dpi = 1000)

# Test for DE features using the DESeq2 package
glu <- FindMarkers(plus, ident.1 = "Glutamatergic neurons", 
                   ident.2 = NULL, only.pos = TRUE,
                   test.use = "DESeq2", max.cells.per.ident = 50)
View(glu)
write.csv(glu,file = "glutamatergic.csv")

gaba <- FindMarkers(plus, ident.1 = "GABAergic neurons", 
                    ident.2 = NULL, only.pos = TRUE,
                    test.use = "DESeq2", max.cells.per.ident = 50)
View(gaba)
write.csv(gaba, file = "gabaergic.csv")

neuro <- FindMarkers(plus, ident.1 = c("GABAergic neurons","Glutamatergic neurons"),
                     ident.2 = NULL, only.pos = TRUE,
                     test.use = "DESeq2", max.cells.per.ident = 50)
View(neuro)
write.csv(neuro, file = "neuro.csv")


myo <- FindMarkers(plus, ident.1 = "Myofibroblasts", 
                   ident.2 = NULL, only.pos = TRUE,
                   test.use = "DESeq2", max.cells.per.ident = 50)
View(myo)
write.csv(myo, file = "Myofibroblasts.csv")

fib <- FindMarkers(plus, ident.1 = "Fibroblasts", 
                   ident.2 = NULL, only.pos = TRUE,
                   test.use = "DESeq2", max.cells.per.ident = 50)
View(fib)
write.csv(fib, file = "Fibroblasts.csv")

imn <- FindMarkers(plus, ident.1 = "Immature neurons", 
                   ident.2 = NULL, only.pos = TRUE,
                   test.use = "DESeq2", max.cells.per.ident = 50)
View(imn)
write.csv(imn, file = "Immature neurons.csv")

unk <- FindMarkers(plus, ident.1 = "Unknown", 
                   ident.2 = NULL, only.pos = TRUE,
                   test.use = "DESeq2", max.cells.per.ident = 50)
View(unk)
write.csv(unk, file = "Unknown.csv")

# volcano plot of DE genes 

neuro.df <- as.data.frame(neuro)
EnhancedVolcano(neuro, x="avg_log2FC", y = "p_val_adj", 
                lab = rownames(neuro), pCutoff = 1e-4, FCcutoff = 1,
                title = NULL, subtitle = NULL)
ggsave("volcanoplot.png", dpi=1000, width = 6.43, height = 6, 
       units = c("in"))


# load packages
library(SingleCellExperiment)
library(slingshot)
library(Seurat)
library(ggplot2)
library(tradeSeq)
library(Matrix)
library(openxlsx)
library(tidyverse)
library(scales)
library(BiocParallel)
library(data.table)
library(patchwork)
library(uwot)
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

# load data
plus <- readRDS("seurat_object")

# Print the current cluster levels
current.cluster.levels <- levels(plus)
print(current.cluster.levels)


# load gene set preparation function
source("D:/seodojin/Rworld/20220509 sctype function-1.R")

# load cell type annotation function
source("D:/seodojin/Rworld/20220509 sctype function-2.R")

# DB file
db_ = "D:/seodojin/Rworld/ScTypeDB_full6.xlsx";
tissue = "Brain" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)





# scale.data 가져오기 (layer 사용)
scale_data <- GetAssayData(plus, layer = "scale.data")

# sctype_score 계산
es.max <- sctype_score(scRNAseqData = scale_data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# 클러스터별로 병합
cL_results <- do.call("rbind", lapply(unique(plus@meta.data$seurat_clusters), function(cl) {
  es.max.cl <- sort(rowSums(es.max[, rownames(plus@meta.data[plus@meta.data$seurat_clusters == cl, ])]), decreasing = TRUE)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(plus@meta.data$seurat_clusters == cl)), 10)
}))

# 가장 높은 점수의 세포 타입 선택
sctype_scores <- cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# 신뢰도가 낮은 클러스터를 "Unknown"으로 설정
sctype_scores$type[as.numeric(sctype_scores$scores) < sctype_scores$ncells / 4] <- "Unknown"
print(sctype_scores[, 1:3])

# UMAP 플롯에 식별된 세포 타입 오버레이
plus@meta.data$customclassif <- ""
for (j in unique(sctype_scores$cluster)) {
  cl_type <- sctype_scores[sctype_scores$cluster == j,]
  plus@meta.data$customclassif[plus@meta.data$seurat_clusters == j] <- as.character(cl_type$type[1])
}




####### 4th july
# Filter out "Unknown" cells
plus@meta.data$customclassif[plus@meta.data$customclassif == ""] <- "Unknown"
plus_filtered <- subset(plus, subset = customclassif != "Unknown")

# Create UMAP plot without "Unknown" cells
DimPlot(plus_filtered, reduction = "umap", group.by = "customclassif", label = TRUE, repel = TRUE) +
  ggtitle("Cell Types (Excluding Unknown)") +
  theme(plot.title = element_text(hjust = 0.5))








# Define custom colors with transparency
cc <- c("shCtrl" = alpha("red", 0.2), "shPTBP1" = alpha("blue", 0.2))
treat <- DimPlot(plus_filtered, 
                 reduction = "umap", 
                 group.by = "orig.ident", 
                 cols = cc)
treat

# 고해상도로 플롯 저장
ggsave("20240705_unknown_filtered_umap_2colors.png", plot = treat, dpi = 1000)

# 클러스터에 세포 타입 이름 할당
new.cluster.ids <- c("Immature\nneurons", "Myofibroblasts", "Fibroblasts", "Unknown", 
                     "Fibroblasts", "Neurons", "Neurons")
names(new.cluster.ids) <- current.cluster.levels
plus <- RenameIdents(plus_filtered, new.cluster.ids)

# 클러스터 색상을 지정
custom_colors <- c("Immature\nneurons" = "yellow", "Myofibroblasts" = "green", 
                   "Fibroblasts" = "red", 
                   "Neurons" = "blue")

# UMAP 플롯 생성
DimPlot(plus, reduction = "umap", label = TRUE, pt.size = 0.5, 
        split.by = "orig.ident", cols = custom_colors, 
        label.size = 4,     # Adjust the label size
        alpha = 0.6) + NoLegend()

# 플롯 저장
ggsave("20240706_unknown_filtered_umap_split_version.png", dpi = 1000)

# 각 clster의 대표 유전자를 feature plot으로... 
FeaturePlot(plus, features = c("PTGIR","SULF1","CALD1"))

ggsave("20240702_4type_gene_featureplot.png", dpi = 1000)






##### cell composition plot
# Load required libraries
library(ggplot2)
library(dplyr)

# Extract metadata
metadata <- plus@meta.data

# Add cell type identity to metadata
metadata$cell_type <- Idents(plus)

# Create a data frame for plotting
cell_type_counts <- metadata %>%
  group_by(orig.ident, cell_type) %>%
  summarize(count = n()) %>%
  ungroup()

# Define the colors for each cell type
cell_type_colors <- c("Immature\nneurons" = "yellow", 
                      "Myofibroblasts" = "green", 
                      "Fibroblasts" = "red", 
                      "Neurons" = "blue")

# Plot the composition of cell types
ggplot(cell_type_counts, aes(x = orig.ident, y = count, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = cell_type_colors) +  # Apply the custom colors
  labs(x = NULL, y = NULL, fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 14, face = "bold"), 
        axis.text.y = element_text(size = 14, face = "bold"))

# Save the plot
ggsave("20240705_cell_type_composition_barplot.png", dpi = 1000)


# Convert Seurat object to SingleCellExperiment object

str(plus, max.level = 2)

library(SingleCellExperiment)

# counts 데이터 추출 및 병합
counts_shCtrl <- plus[["RNA"]]$counts.shCtrl
counts_shPTBP1 <- plus[["RNA"]]$counts.shPTBP1
counts <- cbind(counts_shCtrl, counts_shPTBP1)

# data (normalized) 데이터 추출 및 병합
data_shCtrl <- plus[["RNA"]]$data.shCtrl
data_shPTBP1 <- plus[["RNA"]]$data.shPTBP1
data <- cbind(data_shCtrl, data_shPTBP1)

# SingleCellExperiment 객체 생성
library(SingleCellExperiment)

sce <- SingleCellExperiment(
  assays = list(
    counts = counts,
    logcounts = data
  ),
  colData = plus@meta.data,
  reducedDims = list(
    PCA = Embeddings(plus, "pca"),
    UMAP = Embeddings(plus, "umap")
  )
)

# 활성 식별자 추가
sce$cell_type <- plus@active.ident

# scaling
library(scater)

sce <- logNormCounts(sce)  # 로그 정규화 (이미 수행되었다면 건너뛰기)
sce <- scater::runPCA(sce, exprs_values = "logcounts")  # PCA 실행


####### 5th July



# Run Slingshot with specified start and end clusters


# 클러스터 레이블 업데이트 함수
update_cluster_labels <- function(label) {
  if (label %in% c("GABAergic neurons", "Glutamatergic neurons")) {
    return("Neurons")
  } else {
    return(label)
  }
}

# SingleCellExperiment 객체의 클러스터 레이블 업데이트
sce$updated_customclassif <- sapply(sce$customclassif, update_cluster_labels)

# 업데이트된 고유 클러스터 확인
unique_clusters <- unique(sce$updated_customclassif)
print(unique_clusters)

# Slingshot 실행 (업데이트된 레이블 사용)
sds <- slingshot(sce, clusterLabels = 'updated_customclassif', reducedDim = 'UMAP',
                 start.clus = "Fibroblasts", 
                 end.clus = c("Neurons", "Myofibroblasts"))

# 궤적 그리기
# UMAP 좌표 추출
umap_coords <- reducedDims(sce)$UMAP

# 데이터 프레임 생성
plot_data <- data.frame(
  UMAP1 = umap_coords[,1],
  UMAP2 = umap_coords[,2],
  CellType = sce$updated_customclassif
)

# 색상 정의
cell_colors <- c("Neurons" = "blue", "Immature neurons" = "yellow", 
                 "Myofibroblasts" = "green", "Fibroblasts" = "red")

# 궤적 색상 정의
trajectory_colors <- c("#FF9999", "#66B2FF")

# ggplot으로 그리기
p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(values = cell_colors) +
  theme_minimal() +
  labs(title = "Cell Types with Trajectories")

# 궤적 추가
for (i in seq_along(slingCurves(sds))) {
  curve_data <- slingCurves(sds)[[i]]$s[slingCurves(sds)[[i]]$ord, ]
  
  # 종착점 클러스터 찾기
  end_cluster <- slingLineages(sds)[[i]][length(slingLineages(sds)[[i]])]
  
  # 종착점 클러스터의 중심점 찾기
  end_cluster_cells <- which(sce$updated_customclassif == end_cluster)
  end_point <- colMeans(umap_coords[end_cluster_cells,])
  
  # 종착점까지의 거리 계산
  distances <- sqrt(rowSums((curve_data - matrix(end_point, nrow = nrow(curve_data), ncol = 2, byrow = TRUE))^2))
  
  # 종착점에 가장 가까운 점 찾기
  closest_point <- which.min(distances)
  
  # 종착점까지의 궤적만 그리기
  p <- p + geom_path(data = data.frame(UMAP1 = curve_data[1:closest_point,1], 
                                       UMAP2 = curve_data[1:closest_point,2]),
                     aes(x = UMAP1, y = UMAP2), 
                     color = trajectory_colors[i], 
                     linewidth = 1, 
                     alpha = 0.7)
}

p

# 플롯 저장
ggsave("20240706_updated_trajectory_umap_plot_with_endpoints.png", plot = p, width = 10, height = 8, dpi = 1000)


########## 20240708 pseudotime vln plot
str(pseudotime)

library(ggplot2)
library(dplyr)
library(tidyr)

# Pseudotime 값 추출
pseudotime <- slingPseudotime(sds)

# 데이터 프레임 생성
plot_data <- data.frame(
  Cell = rownames(pseudotime),
  Trajectory1 = pseudotime[,1],
  Trajectory2 = pseudotime[,2],
  Cluster = sce$updated_customclassif
)

# Long 형식으로 변환
plot_data_long <- plot_data %>%
  pivot_longer(cols = c(Trajectory1, Trajectory2), 
               names_to = "Trajectory", 
               values_to = "Pseudotime")

# NA 값 제거
plot_data_long <- plot_data_long %>% filter(!is.na(Pseudotime))

# 클러스터 색상 정의
cluster_colors <- c("Neurons" = "blue", 
                    "Immature neurons" = "yellow", 
                    "Myofibroblasts" = "green", 
                    "Fibroblasts" = "red")

# Violin plot 생성
vp <- ggplot(plot_data_long, aes(x = Pseudotime, y = Cluster, fill = Cluster)) +
  geom_violin(scale = "width") +
  facet_wrap(~Trajectory, scales = "free_x") +
  theme_minimal() +
  labs(title = "Pseudotime Distribution by Cluster and Trajectory",
       x = "Pseudotime",
       y = "Cluster") +
  scale_fill_manual(values = cluster_colors) +  # 여기서 사용자 정의 색상을 적용
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

vp

# 플롯 저장
ggsave("20240708_pseudotime_violin_plot_by_trajectory_custom_colors.png", plot = vp, width = 12, height = 8, dpi = 600)

# SlingshotDataSet 객체 생성
sds <- SlingshotDataSet(sds)

# sds 객체 확인
print(sds)
print(class(sds))



# 1. 각 클러스터에서 10개 셀 샘플링
sampled_cells <- as.data.table(plus@meta.data)[, .(cell = .I[sample(.N, min(.N, 10))]), by = seurat_clusters]
sampled_cells <- sampled_cells$cell

# 2. 샘플링된 셀만 사용하여 새로운 Seurat 객체 생성
plus_sampled <- subset(plus, cells = sampled_cells)

# 3. 희소 행렬로 변환
counts_shCtrl <- plus_sampled[["RNA"]]$counts.shCtrl
counts_shPTBP1 <- plus_sampled[["RNA"]]$counts.shPTBP1
counts_sampled <- cbind(counts_shCtrl, counts_shPTBP1)

counts_sparse_sampled <- as(counts_sampled, "sparseMatrix")






library(SingleCellExperiment)

# counts 데이터 추출 (이전에 생성한 counts_sampled 사용)
counts_matrix <- counts_sampled

# UMAP 좌표 추출
umap_coords <- Embeddings(plus_sampled, reduction = "umap")

# 메타데이터 추출
cell_metadata <- plus_sampled@meta.data

# SingleCellExperiment 객체 생성
sce_sampled <- SingleCellExperiment(assays = list(counts = counts_matrix),
                                    colData = cell_metadata,
                                    reducedDims = list(UMAP = umap_coords))

########### 8th july
# 클러스터 레이블 수정 함수
update_cluster_labels <- function(label) {
  if (label %in% c("GABAergic neurons", "Glutamatergic neurons")) {
    return("Neurons")
  } else {
    return(label)
  }
}

# SingleCellExperiment 객체의 클러스터 레이블 업데이트
sce_sampled$updated_customclassif <- sapply(sce_sampled$customclassif, update_cluster_labels)

# Slingshot 재실행
set.seed(123)
sds_sampled <- slingshot(sce_sampled, clusterLabels = 'updated_customclassif', reducedDim = 'UMAP')

# 업데이트된 lineages 확인
print(slingLineages(sds_sampled))

set.seed(123)
sds_sampled <- slingshot(sce_sampled, clusterLabels = 'updated_customclassif', reducedDim = 'UMAP',
                         start.clus = "Fibroblasts", 
                         end.clus = c("Neurons", "Myofibroblasts"))

print(slingLineages(sds_sampled))

class(sds_sampled)
print(class(sds))

print(dim(assay(sds_sampled, "counts")))
print(dim(slingCurveWeights(sds)))

print(ncol(sds_sampled))
print(ncol(slingCurveWeights(sds)))

# Slingshot 재실행
sds_new <- slingshot(sds_sampled, clusterLabels = 'updated_customclassif', reducedDim = 'UMAP',
                     start.clus = "Fibroblasts", 
                     end.clus = c("Neurons", "Myofibroblasts"))

# SlingshotDataSet 객체 생성
sds <- SlingshotDataSet(sds_new)

# 차원 다시 확인
print(dim(assay(sds_sampled, "counts")))
print(dim(slingCurveWeights(sds)))


sce_fitted <- fitGAM(counts = assay(sds_sampled, "counts"), 
                     sds = sds,
                     nknots = 5,
                     verbose = TRUE,
                     parallel = FALSE)


# fitGAM 실행 후
saveRDS(sce_fitted, file = "fitGAM_results.rds")

# 새로운 R 세션에서
sce_fitted <- readRDS("fitGAM_results.rds")

# sce_fitted 객체의 구조 확인
str(sce_fitted)

# 맞춰진 GAM 모델의 수 확인
print(nrow(sce_fitted))

pattern_results <- patternTest(sce_fitted)
head(pattern_results)

sorted_results <- pattern_results[order(pattern_results$pvalue), ]
sorted_results$adj_pvalue <- p.adjust(sorted_results$pvalue, method = "BH")
significant_genes <- sorted_results[sorted_results$adj_pvalue < 0.05, ]



assayNames(sce_fitted)

# 상위 10개 유전자 선택
top_genes <- rownames(significant_genes)[1:10]

# 각 유전자에 대한 플롯 생성
for (gene in top_genes) {
  p <- plotSmoothers(sce_fitted, gene = gene, counts = assay(sce_fitted, "counts")) +
    ggtitle(paste("Gene:", gene)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(p)
  
  ggsave(paste0("gene_", gene, "_plot.png"), p, width = 8, height = 6)
}





###### 9th july

# earlyDETest 실행
early_results <- earlyDETest(sce_fitted)

# 상위 10개 유전자 선택
top_10_genes <- rownames(early_results)[order(early_results$pvalue)][1:10]

# counts 데이터 추출
counts_data <- assay(sce_fitted, "counts")

# 각 유전자에 대한 플롯 생성
for (gene in top_10_genes) {
  p <- plotSmoothers(sce_fitted, gene = gene, counts = counts_data) +
    ggtitle(paste("Early DE Gene:", gene)) +
    theme_minimal()
  
  print(p)
  
  # 필요하다면 각 플롯을 파일로 저장
  ggsave(paste0("early_DE_gene_", gene, "_plot.png"), p, width = 8, height = 6)
}

# diffEndTest 실행
end_results <- diffEndTest(sce_fitted)

# 상위 10개 유전자 선택
top_10_genes <- rownames(end_results)[order(end_results$pvalue)][1:10]

# counts 데이터 추출
counts_data <- assay(sce_fitted, "counts")

# 각 유전자에 대한 플롯 생성
for (gene in top_10_genes) {
  p <- plotSmoothers(sce_fitted, gene = gene, counts = counts_data) +
    ggtitle(paste("End DE Gene:", gene)) +
    theme_minimal()
  
  print(p)
  
  # 필요하다면 각 플롯을 파일로 저장
  ggsave(paste0("end_DE_gene_", gene, "_plot.png"), p, width = 8, height = 6)
}


# 각 테스트 결과에서 유의미한 유전자 선택 (예: p-value < 0.05)
pattern_genes <- rownames(pattern_results)[pattern_results$pvalue < 0.05]
early_genes <- rownames(early_results)[early_results$pvalue < 0.05]
end_genes <- rownames(end_results)[end_results$pvalue < 0.05]

# 공통 유전자 찾기
common_genes <- Reduce(intersect, list(pattern_genes, early_genes, end_genes))

# 결과 출력
print(paste("공통 유전자 수:", length(common_genes)))
print("공통 유전자 목록:")
print(common_genes)

# 공통 유전자 목록을 파일로 저장
write.csv(common_genes, "common_genes.csv", row.names = FALSE)


library(clusterProfiler)
library(org.Hs.eg.db)  # 인간 유전자 주석 데이터베이스


# GO 분석 수행 (유전자 심볼 사용)
go_result <- enrichGO(gene = common_genes, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "SYMBOL",  # 유전자 심볼을 사용한다고 명시
                      ont = "BP",  # Biological Process
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05)

# 결과 시각화
dotplot(go_result)

# GO 분석 결과 시각화
p <- dotplot(go_result)

# 고해상도로 저장 (DPI 1000)
ggsave("GO_analysis_dotplot.png", plot = p, width = 12, height = 10, units = "in", dpi = 1000)


library(STRINGdb)

# STRINGdb 객체 생성
string_db <- STRINGdb$new(version="11", species=9606)

# 유전자 매핑
mapped_genes <- string_db$map(data.frame(gene = common_genes), "gene", removeUnmappedRows = TRUE)

# 상위 30개 유전자 선택
top_genes <- head(mapped_genes$STRING_id, 20)

# 네트워크 플롯 생성
string_db$plot_network(top_genes, required_score = 700)

# 고해상도 TIFF 파일 생성 시작
tiff("20240709_20_network_plot_high_res.tiff", width = 8, height = 8, units = "in", res = 600, compression = "lzw")

# 네트워크 플롯 생성
string_db$plot_network(top_genes)

# 장치 종료 및 파일 저장
dev.off()

