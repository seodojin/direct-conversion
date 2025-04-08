# Load required packages

library(tradeSeq)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(pheatmap)
library(dorothea)
library(enrichplot)

# Fit tradeSeq model

run_fitGAM <- function(counts, pseudotime, cell_weights, nknots = 5) {
  fitGAM(
    counts = counts,
    pseudotime = as.matrix(pseudotime),
    cellWeights = as.matrix(cell_weights),
    nknots = nknots
  )
}

# Run patternTest or earlyDETest

get_significant_genes <- function(sce, test = c("pattern", "early"), fdr_cutoff = 0.05) {
  test <- match.arg(test)
  res <- if (test == "pattern") patternTest(sce) else earlyDETest(sce)
  res$fdr <- p.adjust(res$pvalue, method = "BH")
  rownames(res)[res$fdr < fdr_cutoff]
}

# Draw heatmap

draw_heatmap <- function(expr_mat, annot, cluster_by = "Cluster", filename, show_rownames = FALSE) {
  annot <- annot[colnames(expr_mat), , drop = FALSE]
  order_idx <- order(annot[[cluster_by]])
  expr_mat <- expr_mat[, order_idx]
  annot <- annot[order_idx, , drop = FALSE]
  
  ann_colors <- list(
    Cluster = c("Fibroblasts" = "purple", "Immature neurons" = "orange", 
                "Myofibroblasts" = "yellow", "Neurons" = "cyan"),
    Lineage = c("Lineage1" = "blue", "Lineage2" = "green", "Lineage3" = "red")
  )
  
  pheatmap(expr_mat,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_colnames = FALSE,
           show_rownames = show_rownames,
           annotation_col = annot[cluster_by, drop = FALSE],
           annotation_colors = ann_colors[cluster_by],
           scale = "row",
           filename = filename)
}

# GO / KEGG enrichment analysis

run_enrichment <- function(genes, ont = "BP", orgdb = org.Hs.eg.db) {
  entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)
  
  ego <- enrichGO(gene = entrez$ENTREZID, OrgDb = orgdb, ont = ont, readable = TRUE)
  ekegg <- enrichKEGG(gene = entrez$ENTREZID, organism = "hsa")
  
  list(GO = ego, KEGG = ekegg)
}

# Identify transcription factors and plot expression

plot_tf_smoothers <- function(sce, gene_list, lineages = NULL, output_prefix = "TF") {
  counts_data <- counts(sce)
  for (tf in gene_list) {
    tryCatch({
      p <- plotSmoothers(sce, counts = counts_data, gene = tf,
                         lineages = lineages, border = TRUE) +
        ggtitle(tf) + theme_minimal()
      ggsave(paste0(output_prefix, "_", tf, ".png"), p, width = 4, height = 3, dpi = 600)
    }, error = function(e) {
      message("Error plotting ", tf)
    })
  }
}

# Run full analysis

run_analysis <- function() {
  sce <- run_fitGAM(counts = assay(sce_filtered_final, "counts"),
                    pseudotime = pseudotime_df_final,
                    cell_weights = cell_weights_df_final)
  
  sig_genes_pattern <- get_significant_genes(sce, "pattern")
  sig_genes_early <- get_significant_genes(sce, "early")
  
  expr <- log2(assay(sce)[sig_genes_pattern, ] + 1)
  expr_scaled <- t(scale(t(expr)))
  
  common_cells <- intersect(colnames(expr_scaled), rownames(combined_data))
  expr_scaled <- expr_scaled[, common_cells]
  metadata <- combined_data[common_cells, , drop = FALSE]
  
  draw_heatmap(expr_scaled, metadata, cluster_by = "Cluster", filename = "heatmap_patternTest.png")
  
  enrich_res <- run_enrichment(sig_genes_early)
  ggsave("GO_dotplot.png", dotplot(enrich_res$GO, showCategory = 10), width = 8, height = 6, dpi = 600)
  ggsave("KEGG_dotplot.png", dotplot(enrich_res$KEGG, showCategory = 10), width = 8, height = 6, dpi = 600)
  
  data(dorothea_hs, package = "dorothea")
  regulons <- dorothea_hs %>% filter(confidence %in% c("A", "B", "C"))
  tfs <- intersect(sig_genes_early, unique(regulons$tf))
  cat("Number of significant transcription factors:", length(tfs), "\n")
  
  plot_tf_smoothers(sce, tfs, lineages = c(1, 3), output_prefix = "TF_lineage1_3")
}

# Run it
run_analysis()
