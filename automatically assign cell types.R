### set up Assigning cell type identity to clusters
### automatically assign cell types using ScType
### reference : https://github.com/IanevskiAleksandr/sc-type/

# load packages
library(openxlsx)
library(tidyverse)
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

# load data
plus <- readRDS("20220831")

# load gene set preparation function
source("D:/seodojin/Rworld/20220509 sctype function-1.R")

# load cell type annotation function
source("D:/seodojin/Rworld/20220509 sctype function-2.R")

# DB file
db_ = "D:/seodojin/Rworld/ScTypeDB_full.xlsx";
tissue = "Brain" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = plus[["RNA"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(plus@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(plus@meta.data[plus@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(plus@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

# overlay the identified cell types on UMAP plot
plus@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  plus@meta.data$customclassif[plus@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(plus, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  

saveRDS(plus, "annotation 20220831")