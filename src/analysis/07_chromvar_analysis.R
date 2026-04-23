#!/usr/bin/env Rscript
# 07_chromvar_analysis.R
# TF Activity Analysis using chromVAR + JASPAR2020
# ATAC-seq 고유 분석: peak 전체에서 TF binding motif 활성도(deviation score) 계산
#
# 요구사항: bioconductor-chromvar, motifmatchr, BSgenome, JASPAR2020, TFBSTools
#   conda install -c bioconda bioconductor-chromvar bioconductor-motifmatchr \
#     bioconductor-jaspar2020 bioconductor-tfbstools \
#     bioconductor-bsgenome.hsapiens.ucsc.hg38 bioconductor-bsgenome.mmusculus.ucsc.mm10
#
# Usage:
#   Rscript src/analysis/07_chromvar_analysis.R <config> <compare> <base> <output_dir>

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
if (length(args) < 4) {
  stop("Usage: Rscript 07_chromvar_analysis.R <config> <compare_group> <base_group> <output_dir>")
}
config_file   <- args[1]
compare_group <- args[2]
base_group    <- args[3]
output_dir    <- args[4]

config <- yaml::read_yaml(config_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("\n=== chromVAR TF Activity Analysis ===")
message(sprintf("Comparison: %s vs %s\n", compare_group, base_group))

# ── Config ─────────────────────────────────────────────────────
chromvar_cfg      <- config$chromvar %||% list()
species           <- config$species
db_cfg            <- config$databases[[species]]
bsgenome_pkg      <- db_cfg$bsgenome
top_tfs_n         <- chromvar_cfg$top_tfs_n         %||% 30
jaspar_collection <- chromvar_cfg$jaspar_collection %||% "CORE"
n_workers         <- chromvar_cfg$n_workers         %||% 1
padj_cutoff       <- config$da_analysis$padj_cutoff  %||% 0.05

if (.Platform$OS.type == "unix" && n_workers > 1) {
  register(MulticoreParam(n_workers))
} else {
  register(SerialParam())
}

# ── Load BSgenome ──────────────────────────────────────────────
message("BSgenome 로딩: ", bsgenome_pkg)
library(bsgenome_pkg, character.only = TRUE)
bsgenome <- get(bsgenome_pkg)

# ── Load data ──────────────────────────────────────────────────
metadata     <- load_metadata(config$metadata_path)
data_list    <- load_featurecounts(config$count_matrix_path, metadata)
count_matrix <- data_list$count_matrix
peak_info    <- data_list$peak_info

group_var    <- config$da_analysis$group_variable
keep_samples <- rownames(metadata)[metadata[[group_var]] %in% c(compare_group, base_group)]
if (length(keep_samples) == 0) {
  stop("지정한 그룹이 metadata에 없습니다: ", compare_group, ", ", base_group,
       "\n사용 가능한 그룹: ", paste(unique(metadata[[group_var]]), collapse = ", "))
}

count_matrix <- count_matrix[, keep_samples, drop = FALSE]
meta_sub     <- metadata[keep_samples, , drop = FALSE]
meta_sub[[group_var]] <- factor(meta_sub[[group_var]], levels = c(base_group, compare_group))

# ── Peak filtering ─────────────────────────────────────────────
min_count    <- config$da_analysis$min_peak_count %||% 10
min_sample_n <- config$da_analysis$min_sample_n   %||% 2
keep_peaks   <- rowSums(count_matrix >= min_count) >= min_sample_n
count_matrix <- count_matrix[keep_peaks, , drop = FALSE]
peak_info    <- peak_info[keep_peaks, , drop = FALSE]
message(sprintf("Peaks: %d, Samples: %d", nrow(count_matrix), ncol(count_matrix)))

# ── Build SummarizedExperiment ─────────────────────────────────
peaks_gr <- GRanges(
  seqnames = peak_info$Chr,
  ranges   = IRanges(start = as.integer(peak_info$Start),
                     end   = as.integer(peak_info$End))
)

se <- SummarizedExperiment(
  assays    = list(counts = count_matrix),
  rowRanges = peaks_gr,
  colData   = DataFrame(meta_sub)
)

# ── GC bias ────────────────────────────────────────────────────
message("GC bias 계산 중...")
se <- addGCBias(se, genome = bsgenome)

# ── JASPAR2020 motifs ──────────────────────────────────────────
message("JASPAR2020 motif 로딩: collection = ", jaspar_collection)
motifs <- getMatrixSet(
  JASPAR2020,
  opts = list(collection  = jaspar_collection,
              tax_group   = "vertebrates",
              all_versions = FALSE)
)
message(sprintf("  %d motifs 로드됨", length(motifs)))

# ── Motif-peak matching ────────────────────────────────────────
message("Motif-peak 매칭 중...")
motif_ix <- matchMotifs(motifs, se, genome = bsgenome)

# ── Background peaks ───────────────────────────────────────────
message("Background peaks 계산 중...")
bg <- getBackgroundPeaks(se, niterations = 200)

# ── Compute deviations ─────────────────────────────────────────
message("Chromatin deviations 계산 중...")
dev   <- computeDeviations(object = se, annotations = motif_ix,
                           background_peaks = bg)
z_mat <- deviationScores(dev)   # motifs × samples (z-score)

# ── TF variability ─────────────────────────────────────────────
var_res <- computeVariability(dev)
var_df  <- as.data.frame(var_res) %>%
  rownames_to_column("motif") %>%
  arrange(desc(variability))

write.csv(var_df, file.path(output_dir, "tf_variability.csv"), row.names = FALSE)
message(sprintf("TF variability 저장: %d motifs", nrow(var_df)))

# Variability barplot
top_var <- var_df %>% slice_head(n = top_tfs_n)
p_var <- ggplot(top_var, aes(x = reorder(motif, variability), y = variability)) +
  geom_col(fill = "#5B8DB8", alpha = 0.85) +
  coord_flip() +
  labs(
    title = sprintf("Top %d Variable TFs (%s vs %s)", top_tfs_n, compare_group, base_group),
    x = NULL, y = "Variability"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "tf_variability_plot.png"),
       p_var, width = 8, height = max(4, top_tfs_n * 0.32), dpi = 150)
message("Variability plot 저장 완료")

# ── Differential TF activity (Wilcoxon per motif) ──────────────
message("Differential TF activity 계산 중...")
compare_idx <- which(meta_sub[[group_var]] == compare_group)
base_idx    <- which(meta_sub[[group_var]] == base_group)

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
diff_df       <- bind_rows(diff_list)
diff_df$padj  <- p.adjust(diff_df$p_value, method = "BH")
diff_df       <- diff_df %>% arrange(padj)

write.csv(diff_df, file.path(output_dir, "differential_tf_activity.csv"),
          row.names = FALSE)
n_sig <- sum(!is.na(diff_df$padj) & diff_df$padj < padj_cutoff)
message(sprintf("Differential TF activity 저장: %d significant (padj < %.2f)",
                n_sig, padj_cutoff))

# ── Deviation heatmap (top variable TFs) ──────────────────────
top_motifs <- head(var_df$motif, top_tfs_n)
heat_mat   <- z_mat[top_motifs, , drop = FALSE]

ann_col <- data.frame(group = meta_sub[[group_var]],
                      row.names = colnames(heat_mat))
n_grp <- nlevels(meta_sub[[group_var]])
pal   <- setNames(brewer.pal(max(3, n_grp), "Set2")[seq_len(n_grp)],
                  levels(meta_sub[[group_var]]))

png(file.path(output_dir, "tf_deviation_heatmap.png"),
    width = 900, height = max(600, top_tfs_n * 22), res = 120)
pheatmap(heat_mat,
         annotation_col    = ann_col,
         annotation_colors = list(group = pal),
         scale             = "row",
         show_rownames     = TRUE,
         show_colnames     = TRUE,
         fontsize_row      = 8,
         fontsize          = 10,
         main = sprintf("TF Activity (z-score) — %s vs %s", compare_group, base_group))
dev.off()
message("TF deviation heatmap 저장 완료")

message(sprintf("\nchromVAR 분석 완료. 출력: %s", output_dir))
