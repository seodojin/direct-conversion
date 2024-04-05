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
plus <- readRDS("annotation_object")

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
       units = c("in")
