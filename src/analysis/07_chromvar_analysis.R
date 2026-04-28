#!/usr/bin/env Rscript
# 07_chromvar_analysis.R
# TF Activity Analysis using chromVAR + JASPAR2020
# Global analysis: deviation scores computed across ALL samples once,
# then differential TF test run for each pairwise comparison in config.
#
# Usage:
#   Rscript src/analysis/07_chromvar_analysis.R <config> <output_dir>

suppressPackageStartupMessages({
  library(yaml)
  library(chromVAR)
  library(motifmatchr)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(BiocParallel)
  library(JASPAR2020)
  library(TFBSTools)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

source("src/utils/load_data.R")

`%||%` <- function(a, b) if (!is.null(a)) a else b

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 07_chromvar_analysis.R <config> <output_dir>")
}
config_file <- args[1]
output_dir  <- args[2]

config <- yaml::read_yaml(config_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "differential_tf"), recursive = TRUE, showWarnings = FALSE)

message("\n=== chromVAR TF Activity Analysis (Global) ===\n")

# ── Config ─────────────────────────────────────────────────────
chromvar_cfg      <- config$chromvar %||% list()
species           <- config$species
db_cfg            <- config$databases[[species]]
bsgenome_pkg      <- db_cfg$bsgenome
top_tfs_n         <- chromvar_cfg$top_tfs_n         %||% 30
jaspar_collection <- chromvar_cfg$jaspar_collection %||% "CORE"
n_workers         <- chromvar_cfg$n_workers         %||% 1
padj_cutoff       <- config$da_analysis$padj_cutoff  %||% 0.05
group_var         <- config$da_analysis$group_variable
pairs             <- config$da_analysis$pairwise_comparisons

if (.Platform$OS.type == "unix" && n_workers > 1) {
  register(MulticoreParam(n_workers))
} else {
  register(SerialParam())
}

# ── Load BSgenome ──────────────────────────────────────────────
message("BSgenome 로딩: ", bsgenome_pkg)
library(bsgenome_pkg, character.only = TRUE)
bsgenome <- get(bsgenome_pkg)

# ── Load ALL samples ───────────────────────────────────────────
metadata     <- load_metadata(config$metadata_path)
data_list    <- load_featurecounts(config$count_matrix_path, metadata)
count_matrix <- data_list$count_matrix
peak_info    <- data_list$peak_info

message(sprintf("전체 샘플: %d, peaks: %d", ncol(count_matrix), nrow(count_matrix)))

# Peak filtering (전체 샘플 기준)
min_count    <- config$da_analysis$min_peak_count %||% 10
min_sample_n <- config$da_analysis$min_sample_n   %||% 2
keep_peaks   <- rowSums(count_matrix >= min_count) >= min_sample_n
count_matrix <- count_matrix[keep_peaks, , drop = FALSE]
peak_info    <- peak_info[keep_peaks, , drop = FALSE]
message(sprintf("필터링 후 peaks: %d", nrow(count_matrix)))

# ── Build SummarizedExperiment (전체 샘플) ─────────────────────
chr_names <- peak_info$Chr
if (!any(grepl("^chr", chr_names))) {
  chr_names <- paste0("chr", chr_names)
  message("  'chr' prefix 추가 (UCSC 형식)")
}
peaks_gr <- GRanges(
  seqnames = chr_names,
  ranges   = IRanges(start = as.integer(peak_info$Start),
                     end   = as.integer(peak_info$End))
)
se <- SummarizedExperiment(
  assays    = list(counts = count_matrix),
  rowRanges = peaks_gr,
  colData   = DataFrame(metadata[colnames(count_matrix), , drop = FALSE])
)

# ── GC bias ────────────────────────────────────────────────────
message("GC bias 계산 중...")
se <- addGCBias(se, genome = bsgenome)

# ── JASPAR2020 motifs ──────────────────────────────────────────
message("JASPAR2020 motif 로딩: collection = ", jaspar_collection)
motifs <- getMatrixSet(
  JASPAR2020,
  opts = list(collection   = jaspar_collection,
              tax_group    = "vertebrates",
              all_versions = FALSE)
)
message(sprintf("  %d motifs 로드됨", length(motifs)))

# ── Motif-peak matching ────────────────────────────────────────
message("Motif-peak 매칭 중...")
motif_ix <- matchMotifs(motifs, se, genome = bsgenome)

# ── Background peaks ───────────────────────────────────────────
message("Background peaks 계산 중...")
bg <- getBackgroundPeaks(se, niterations = 200)

# ── Compute deviations (전체 샘플, 1회) ───────────────────────
message("Chromatin deviations 계산 중...")
dev   <- computeDeviations(object = se, annotations = motif_ix,
                           background_peaks = bg)
z_mat <- deviationScores(dev)   # motifs × all samples

# ── TF variability (전체 샘플 기준) ───────────────────────────
var_res <- computeVariability(dev)
var_df  <- as.data.frame(var_res) %>%
  rownames_to_column("motif") %>%
  arrange(desc(variability))

write.csv(var_df, file.path(output_dir, "tf_variability.csv"), row.names = FALSE)
message(sprintf("TF variability 저장: %d motifs", nrow(var_df)))

# Variability barplot (전체 샘플)
top_var <- var_df %>% slice_head(n = top_tfs_n)
p_var <- ggplot(top_var, aes(x = reorder(motif, variability), y = variability)) +
  geom_col(fill = "#5B8DB8", alpha = 0.85) +
  coord_flip() +
  labs(title = sprintf("Top %d Variable TFs (All Samples)", top_tfs_n),
       x = NULL, y = "Variability") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(output_dir, "tf_variability_plot.png"),
       p_var, width = 8, height = max(4, top_tfs_n * 0.32), dpi = 150)

# ── Deviation heatmap (전체 샘플 × top TFs) ───────────────────
top_motifs <- head(var_df$motif, top_tfs_n)
heat_mat   <- z_mat[top_motifs, , drop = FALSE]

conditions  <- metadata[colnames(heat_mat), group_var, drop = TRUE]
cond_levels <- unique(conditions)
n_grp  <- length(cond_levels)
pal    <- setNames(
  colorRampPalette(brewer.pal(min(n_grp, 9), "Set1"))(n_grp),
  cond_levels
)
ann_col <- data.frame(condition = conditions, row.names = colnames(heat_mat))

png(file.path(output_dir, "tf_deviation_heatmap.png"),
    width = max(900, ncol(heat_mat) * 40 + 300),
    height = max(600, top_tfs_n * 22), res = 120)
pheatmap(heat_mat,
         annotation_col    = ann_col,
         annotation_colors = list(condition = pal),
         scale             = "row",
         show_rownames     = TRUE,
         show_colnames     = TRUE,
         fontsize_row      = 8,
         fontsize          = 10,
         main = sprintf("TF Activity Z-score — All Samples (top %d TFs)", top_tfs_n))
dev.off()
message("TF deviation heatmap 저장 완료")

# ── Pairwise differential TF activity ─────────────────────────
message("\n=== Pairwise Differential TF Analysis ===")
for (pair in pairs) {
  compare_group <- pair[[1]]
  base_group    <- pair[[2]]
  pair_name     <- sprintf("%s_vs_%s", compare_group, base_group)
  message(sprintf("\n  [%s]", pair_name))

  compare_idx <- which(metadata[colnames(z_mat), group_var] == compare_group)
  base_idx    <- which(metadata[colnames(z_mat), group_var] == base_group)

  if (length(compare_idx) == 0 || length(base_idx) == 0) {
    message(sprintf("    경고: 그룹 없음 (%s 또는 %s), 건너뜀", compare_group, base_group))
    next
  }

  diff_list <- lapply(rownames(z_mat), function(motif) {
    g1 <- z_mat[motif, compare_idx]
    g2 <- z_mat[motif, base_idx]
    wt <- tryCatch(
      wilcox.test(g1, g2, exact = FALSE),
      error = function(e) list(p.value = NA_real_)
    )
    data.frame(
      motif        = motif,
      mean_compare = mean(g1, na.rm = TRUE),
      mean_base    = mean(g2, na.rm = TRUE),
      delta        = mean(g1, na.rm = TRUE) - mean(g2, na.rm = TRUE),
      p_value      = wt$p.value,
      stringsAsFactors = FALSE
    )
  })
  diff_df      <- bind_rows(diff_list)
  diff_df$padj <- p.adjust(diff_df$p_value, method = "BH")
  diff_df      <- diff_df %>% arrange(padj)

  out_csv <- file.path(output_dir, "differential_tf", sprintf("%s_diff_tf.csv", pair_name))
  write.csv(diff_df, out_csv, row.names = FALSE)
  n_sig <- sum(!is.na(diff_df$padj) & diff_df$padj < padj_cutoff)
  message(sprintf("    저장: %s (significant: %d)", basename(out_csv), n_sig))

  # Top differential TF barplot
  top_diff <- diff_df %>%
    filter(!is.na(padj)) %>%
    slice_head(n = min(20, nrow(.)))
  if (nrow(top_diff) > 0) {
    top_diff$direction <- ifelse(top_diff$delta > 0, "Up", "Down")
    p_diff <- ggplot(top_diff,
                     aes(x = reorder(motif, delta), y = delta,
                         fill = direction)) +
      geom_col(alpha = 0.85) +
      scale_fill_manual(values = c("Up" = "#E57373", "Down" = "#64B5F6")) +
      coord_flip() +
      labs(title = sprintf("Differential TF Activity: %s", pair_name),
           x = NULL, y = "Delta Z-score (compare - base)", fill = NULL) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
    out_png <- file.path(output_dir, "differential_tf",
                         sprintf("%s_diff_tf_plot.png", pair_name))
    ggsave(out_png, p_diff, width = 8, height = max(4, nrow(top_diff) * 0.35), dpi = 150)
  }
}

message(sprintf("\nchromVAR 분석 완료. 출력: %s", output_dir))
