#!/usr/bin/env Rscript
# 02a_global_qc.R
# Global QC plots: PCA, sample distance heatmap (모든 샘플 대상)
# RNA-Seq_DE_GO_analysis/src/analysis/02a_generate_qc_plots.R 패턴 계승
#
# Usage:
#   Rscript src/analysis/02a_global_qc.R <config> <output_dir>

suppressPackageStartupMessages({
  library(yaml)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tibble)
})

source("src/utils/load_data.R")

`%||%` <- function(a, b) if (!is.null(a)) a else b

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 02a_global_qc.R <config> <output_dir>")
config_file <- args[1]
output_dir  <- args[2]

config <- yaml::read_yaml(config_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("\n=== Global QC Plots ===\n")

# ── Load data ──────────────────────────────────────────────────
metadata  <- load_metadata(config$metadata_path)
data_list <- load_featurecounts(config$count_matrix_path, metadata)
count_matrix <- data_list$count_matrix
group_var    <- config$da_analysis$group_variable

min_count    <- config$da_analysis$min_peak_count %||% 10
min_sample_n <- config$da_analysis$min_sample_n   %||% 2
keep_peaks   <- rowSums(count_matrix >= min_count) >= min_sample_n
count_matrix <- count_matrix[keep_peaks, , drop = FALSE]

# DESeq2 object for normalization
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = metadata[colnames(count_matrix), , drop = FALSE],
  design    = ~ 1
)
dds <- estimateSizeFactors(dds)
vst_mat <- tryCatch(assay(vst(dds, blind = TRUE)),
                    error = function(e) assay(normTransform(dds)))

# Top variable peaks for PCA / heatmap
top_n   <- config$qc_plots$top_peaks_n %||% 500
rv      <- rowVars(vst_mat)
top_idx <- order(rv, decreasing = TRUE)[seq_len(min(top_n, length(rv)))]
vst_top <- vst_mat[top_idx, ]

# ── PCA ───────────────────────────────────────────────────────
pca_res  <- prcomp(t(vst_top), scale. = FALSE)
pca_df   <- as.data.frame(pca_res$x[, 1:2]) %>%
  rownames_to_column("sample_id") %>%
  left_join(rownames_to_column(metadata, "sample_id"), by = "sample_id")
var_exp  <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2,
                              color = .data[[group_var]],
                              label = sample_id)) +
  geom_point(size = 4, alpha = 0.85) +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
  labs(
    title    = "PCA — All Samples (top variable peaks)",
    x        = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y        = sprintf("PC2 (%.1f%%)", var_exp[2]),
    color    = group_var
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right")

ggsave(file.path(output_dir, "pca_all_samples.png"), p_pca,
       width = 7, height = 5, dpi = 150)
message("PCA 저장 완료")

# ── Sample distance heatmap ────────────────────────────────────
dist_method <- config$qc_plots$sample_distance_method %||% "euclidean"
samp_dists  <- as.matrix(dist(t(vst_top), method = dist_method))

ann_col <- data.frame(
  group = metadata[colnames(samp_dists), group_var, drop = TRUE],
  row.names = colnames(samp_dists)
)
n_groups <- length(unique(ann_col$group))
pal      <- setNames(brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)],
                     unique(ann_col$group))

png(file.path(output_dir, "sample_distance_heatmap.png"),
    width = 800, height = 700, res = 120)
pheatmap(samp_dists,
         clustering_distance_rows = dist_method,
         clustering_distance_cols = dist_method,
         annotation_col = ann_col,
         annotation_colors = list(group = pal),
         main = "Sample Distance Heatmap",
         fontsize = 11)
dev.off()
message("Sample distance heatmap 저장 완료")

# ── Top variable peaks heatmap ─────────────────────────────────
top_heatmap_n <- min(100, nrow(vst_top))
heat_mat <- vst_top[seq_len(top_heatmap_n), ]
heat_mat <- t(scale(t(heat_mat)))

png(file.path(output_dir, "top_peaks_heatmap.png"),
    width = 800, height = 900, res = 120)
pheatmap(heat_mat,
         show_rownames   = FALSE,
         annotation_col  = ann_col,
         annotation_colors = list(group = pal),
         main     = sprintf("Top %d Variable Peaks (z-score)", top_heatmap_n),
         fontsize = 11,
         cluster_rows = config$qc_plots$heatmap_cluster_rows %||% TRUE,
         cluster_cols = config$qc_plots$heatmap_cluster_cols %||% TRUE)
dev.off()
message("Top peaks heatmap 저장 완료")

message("\nGlobal QC 완료. 출력: ", output_dir)
