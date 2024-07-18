### set up Assigning cell type identity to clusters
### automatically assign cell types using ScType
### reference : https://github.com/IanevskiAleksandr/sc-type/
## input file : seurat_object.rds, ScTypeDB_full.xlsx
# load packages
library(openxlsx)
library(tidyverse)

lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

# load data
plus <- readRDS("data/seurat_object")

# load gene set preparation function
source("data/20220509 sctype function-1.R")

# load cell type annotation function
source("data/20220509 sctype function-2.R")

# DB file
db_ = "data/ScTypeDB_full.xlsx";
tissue = "Brain" 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get scale.data
scale_data <- GetAssayData(plus, layer = "scale.data")

# calculate sctype_score 
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

# We can also overlay the identified cell types on UMAP plot
plus@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  plus@meta.data$customclassif[plus@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(plus, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  

# Filter out "Unknown" cells
plus@meta.data$customclassif[plus@meta.data$customclassif == ""] <- "Unknown"
plus_filtered <- subset(plus, subset = customclassif != "Unknown")

# Create UMAP plot without "Unknown" cells
DimPlot(plus_filtered, reduction = "umap", group.by = "customclassif", label = TRUE, repel = TRUE) +
  ggtitle("Cell Types (Excluding Unknown)") +
  theme(plot.title = element_text(hjust = 0.5))

# Define custom colors with transparency
cc <- c("shCtrl" = alpha("red", 0.2), "shPTBP1" = alpha("blue", 0.2))
plot1 <- DimPlot(plus_filtered, 
                 reduction = "umap", 
                 group.by = "orig.ident", 
                 cols = cc)
plot1

# 고해상도로 플롯 저장
ggsave("20240711_unknown_filtered_umap_2colors.png", plot = plot1, dpi = 1000)

######### ---------------------------- 11th july

# Print the current cluster levels
current.cluster.levels <- levels(plus)
print(current.cluster.levels)

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

ggsave("20240712_unknown_filtered_umap_split_version.png", dpi = 1000)
