#!/usr/bin/env Rscript
# 01_run_pairwise_da.R
# Differential Accessibility Analysis (DESeq2)
# RNA-Seq_DE_GO_analysis/src/analysis/01b_run_pairwise_de.R 패턴 계승
#
# Usage:
#   Rscript src/analysis/01_run_pairwise_da.R <config> <compare> <base> <output_dir>

suppressPackageStartupMessages({
  library(yaml)
  library(DESeq2)
  library(apeglm)
  library(ashr)
  library(dplyr)
  library(tibble)
})

source("src/utils/load_data.R")

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Arguments ──────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript 01_run_pairwise_da.R <config> <compare_group> <base_group> <output_dir>")
}
config_file   <- args[1]
compare_group <- args[2]
base_group    <- args[3]
output_dir    <- args[4]

config <- yaml::read_yaml(config_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

message(sprintf("\n========================================"))
message(sprintf("DA Analysis: %s vs %s", compare_group, base_group))
message(sprintf("========================================\n"))

# ── Load data ──────────────────────────────────────────────────
metadata  <- load_metadata(config$metadata_path)
data_list <- load_featurecounts(config$count_matrix_path, metadata)
count_matrix <- data_list$count_matrix

group_var <- config$da_analysis$group_variable

# 비교군/기준군 샘플만 추출
keep_samples <- rownames(metadata)[metadata[[group_var]] %in% c(compare_group, base_group)]
if (length(keep_samples) == 0) {
  stop("지정한 그룹이 metadata에 없습니다: ", compare_group, ", ", base_group,
       "\n사용 가능한 그룹: ", paste(unique(metadata[[group_var]]), collapse = ", "))
}

count_matrix <- count_matrix[, keep_samples, drop = FALSE]
meta_sub     <- metadata[keep_samples, , drop = FALSE]
meta_sub[[group_var]] <- factor(meta_sub[[group_var]], levels = c(base_group, compare_group))

message(sprintf("샘플 수 - %s: %d, %s: %d",
  compare_group, sum(meta_sub[[group_var]] == compare_group),
  base_group,    sum(meta_sub[[group_var]] == base_group)))

# ── Peak filtering ─────────────────────────────────────────────
min_count    <- config$da_analysis$min_peak_count %||% 10
min_sample_n <- config$da_analysis$min_sample_n   %||% 2

keep_peaks <- rowSums(count_matrix >= min_count) >= min_sample_n
count_matrix <- count_matrix[keep_peaks, , drop = FALSE]
message(sprintf("Peak 필터링: %d / %d 통과 (min count=%d, min samples=%d)",
                sum(keep_peaks), length(keep_peaks), min_count, min_sample_n))

# ── DESeq2 ─────────────────────────────────────────────────────
design_formula <- as.formula(config$da_analysis$design_formula)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = meta_sub,
  design    = design_formula
)
dds <- DESeq(dds, quiet = FALSE)

# ── Results ────────────────────────────────────────────────────
padj_cutoff <- config$da_analysis$padj_cutoff  %||% 0.05
lfc_cutoff  <- config$da_analysis$log2fc_cutoff %||% 1.0

# 기본 results
res_name <- resultsNames(dds)
compare_coef <- res_name[grepl(paste0("_", compare_group, "$"), res_name)]

res <- results(dds,
               contrast = c(group_var, compare_group, base_group),
               alpha    = padj_cutoff)

# LFC shrinkage (apeglm): coef 이름이 있을 때만
res_shrunk <- tryCatch({
  if (length(compare_coef) == 1) {
    lfcShrink(dds, coef = compare_coef, type = "apeglm", res = res, quiet = TRUE)
  } else {
    message("apeglm shrinkage 건너뜀 (coef 이름 불명확). ashr 사용.")
    lfcShrink(dds, contrast = c(group_var, compare_group, base_group),
              type = "ashr", res = res, quiet = TRUE)
  }
}, error = function(e) {
  message("LFC shrinkage 실패, shrinkage 없이 진행: ", conditionMessage(e))
  res
})

res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("peak_id") %>%
  arrange(padj)

# ── Save results ───────────────────────────────────────────────
write.csv(res_df, file.path(output_dir, "final_da_results.csv"), row.names = FALSE)
message("DA 결과 저장: ", file.path(output_dir, "final_da_results.csv"))

# Normalized counts (VST)
vst_mat <- tryCatch(
  assay(vst(dds, blind = FALSE)),
  error = function(e) assay(normTransform(dds))
)
norm_df <- as.data.frame(vst_mat) %>% rownames_to_column("peak_id")
write.csv(norm_df, file.path(output_dir, "normalized_counts.csv"), row.names = FALSE)

# Config copy (재현성)
file.copy(config_file, file.path(output_dir, "config_used.yml"), overwrite = TRUE)

# ── Summary ────────────────────────────────────────────────────
sig <- res_df %>% filter(!is.na(padj), padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)
message(sprintf("\n[Summary] %s vs %s", compare_group, base_group))
message(sprintf("  분석된 peaks:   %d", nrow(res_df)))
message(sprintf("  Significant:    %d  (padj < %.2f, |log2FC| > %.1f)",
                nrow(sig), padj_cutoff, lfc_cutoff))
message(sprintf("  More open (%s): %d", compare_group, sum(sig$log2FoldChange > 0)))
message(sprintf("  Less open (%s): %d", compare_group, sum(sig$log2FoldChange < 0)))
message("완료.")
