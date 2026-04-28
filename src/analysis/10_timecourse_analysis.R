#!/usr/bin/env Rscript
# 10_timecourse_analysis.R
# Temporal trajectory analysis for time-series ATAC-seq data
# Union of all DA peaks → Z-score across timepoints → K-means clustering
# Trend plots + temporal heatmap + GO enrichment per cluster
#
# Usage:
#   Rscript src/analysis/10_timecourse_analysis.R <config> <output_dir>

suppressPackageStartupMessages({
  library(yaml)
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(factoextra)
  library(clusterProfiler)
})

source("src/utils/load_data.R")

`%||%` <- function(a, b) if (!is.null(a)) a else b

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 10_timecourse_analysis.R <config> <output_dir>")
config_file <- args[1]
output_dir  <- args[2]

config <- yaml::read_yaml(config_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "go_enrichment"), recursive = TRUE, showWarnings = FALSE)

message("\n=== Time-course Temporal Analysis ===\n")

# ── Config ─────────────────────────────────────────────────────
tc_cfg       <- config$timecourse %||% list()
padj_cutoff  <- config$da_analysis$padj_cutoff   %||% 0.05
lfc_cutoff   <- config$da_analysis$log2fc_cutoff %||% 1.0
time_order   <- tc_cfg$time_order   %||% stop("timecourse.time_order 설정 필요")
time_hours   <- tc_cfg$time_hours   %||% seq_along(time_order)
n_clusters   <- tc_cfg$n_clusters   # null = auto
cluster_meth <- tc_cfg$cluster_method %||% "kmeans"
go_ontology  <- tc_cfg$go_ontology  %||% "BP"
min_peaks_go <- tc_cfg$min_peaks_for_go %||% 50
group_var    <- config$da_analysis$group_variable
pairs        <- config$da_analysis$pairwise_comparisons
output_base  <- config$output_dir

species      <- config$species
db_cfg       <- config$databases[[species]]
org_db_name  <- db_cfg$organism_db
library(org_db_name, character.only = TRUE)
org_db       <- get(org_db_name)

pair_names   <- sapply(pairs, function(p) sprintf("%s_vs_%s", p[[1]], p[[2]]))

# ── Load DA results → union peaks ─────────────────────────────
message("DA results 로딩 및 union peaks 추출...")
sig_peaks <- character(0)
for (pair in pair_names) {
  path <- file.path(output_base, "pairwise", pair, "final_da_results.csv")
  if (!file.exists(path)) { message(sprintf("  경고: 파일 없음 %s", path)); next }
  df <- load_da_results(path)
  sig <- df %>%
    filter(!is.na(padj), padj < padj_cutoff, abs(log2FoldChange) >= lfc_cutoff) %>%
    pull(peak_id)
  sig_peaks <- union(sig_peaks, sig)
  message(sprintf("  [%s] %d DA peaks (union: %d)", pair, length(sig), length(sig_peaks)))
}
message(sprintf("Union DA peaks: %d", length(sig_peaks)))
if (length(sig_peaks) < 10) stop("Union DA peaks가 너무 적습니다 (< 10).")

# ── Load all samples raw counts → VST (전체 샘플 기준) ─────────
# 각 pairwise normalized_counts는 해당 비교군 6샘플만 포함하므로
# 전체 샘플 VST는 원본 featureCounts에서 직접 계산
message("\nAll-sample VST 정규화 계산...")
data_list    <- load_featurecounts(config$count_matrix_path,
                                   load_metadata(config$metadata_path))
raw_counts   <- data_list$count_matrix

min_count    <- config$da_analysis$min_peak_count %||% 10
min_sample_n <- config$da_analysis$min_sample_n   %||% 2
keep_peaks   <- rowSums(raw_counts >= min_count) >= min_sample_n
raw_counts   <- raw_counts[keep_peaks, , drop = FALSE]
message(sprintf("  필터링 후 peaks: %d × %d samples", nrow(raw_counts), ncol(raw_counts)))

metadata_all <- load_metadata(config$metadata_path)
col_data     <- data.frame(
  condition = factor(metadata_all[colnames(raw_counts), group_var]),
  row.names = colnames(raw_counts)
)
dds_all   <- DESeqDataSetFromMatrix(countData = raw_counts,
                                    colData   = col_data,
                                    design    = ~ condition)
vst_all   <- assay(vst(dds_all, blind = TRUE))
norm_counts <- as.data.frame(vst_all)
message(sprintf("  VST 완료: %d peaks × %d samples", nrow(norm_counts), ncol(norm_counts)))

# union peaks 서브셋
available <- intersect(sig_peaks, rownames(norm_counts))
if (length(available) < length(sig_peaks)) {
  message(sprintf("  경고: %d peaks가 count matrix에 없어 제외 (저카운트 필터)",
                  length(sig_peaks) - length(available)))
}
norm_sub <- norm_counts[available, , drop = FALSE]

# ── Condition별 평균 (timepoint order) ────────────────────────
message("\nCondition별 평균 계산...")
metadata <- load_metadata(config$metadata_path)
cond_means <- sapply(time_order, function(cond) {
  samps <- rownames(metadata)[metadata[[group_var]] == cond]
  samps <- intersect(samps, colnames(norm_sub))
  if (length(samps) == 0) {
    message(sprintf("  경고: condition '%s' 샘플 없음", cond))
    return(rep(NA_real_, nrow(norm_sub)))
  }
  rowMeans(norm_sub[, samps, drop = FALSE], na.rm = TRUE)
})
colnames(cond_means) <- time_order
cond_means <- as.data.frame(cond_means)

# Remove peaks with any NA timepoint
na_peaks <- rowSums(is.na(cond_means)) > 0
if (any(na_peaks)) {
  message(sprintf("  NA 포함 peaks 제외: %d", sum(na_peaks)))
  cond_means <- cond_means[!na_peaks, , drop = FALSE]
}

# ── Z-score normalization (peak 단위) ──────────────────────────
message("Z-score 정규화...")
z_mat <- t(scale(t(cond_means)))
# peaks with sd=0 → NaN, 제외
valid <- complete.cases(z_mat)
if (!all(valid)) {
  message(sprintf("  분산=0 peaks 제외: %d", sum(!valid)))
  z_mat     <- z_mat[valid, , drop = FALSE]
  cond_means <- cond_means[valid, , drop = FALSE]
}
message(sprintf("  분석 가능 peaks: %d", nrow(z_mat)))

# ── Optimal k (elbow + silhouette) ────────────────────────────
if (is.null(n_clusters)) {
  message("\nOptimal k 계산 중 (elbow + silhouette, k=3~10)...")
  set.seed(42)
  k_range <- 3:min(10, nrow(z_mat) - 1)

  # Elbow (WSS)
  wss <- sapply(k_range, function(k) {
    km <- kmeans(z_mat, centers = k, nstart = 25, iter.max = 100)
    km$tot.withinss
  })
  wss_df <- data.frame(k = k_range, wss = wss)
  p_elbow <- ggplot(wss_df, aes(k, wss)) +
    geom_line(color = "#5B8DB8") + geom_point(size = 3, color = "#5B8DB8") +
    scale_x_continuous(breaks = k_range) +
    labs(title = "Elbow Method (WSS)", x = "Number of clusters (k)", y = "Total WSS") +
    theme_bw(base_size = 12)
  ggsave(file.path(output_dir, "elbow_plot.png"), p_elbow, width = 7, height = 5, dpi = 150)

  # Silhouette (factoextra)
  sil_res <- fviz_nbclust(z_mat, kmeans, method = "silhouette",
                           k.max = max(k_range), nstart = 25)
  ggsave(file.path(output_dir, "silhouette_plot.png"), sil_res, width = 7, height = 5, dpi = 150)

  # 최적 k: silhouette 최대값
  sil_data <- sil_res$data
  n_clusters <- as.integer(sil_data$clusters[which.max(sil_data$y)])
  message(sprintf("  → 최적 k = %d (silhouette)", n_clusters))
} else {
  message(sprintf("  지정된 k = %d", n_clusters))
}

# ── K-means clustering ─────────────────────────────────────────
message(sprintf("\nK-means clustering (k=%d)...", n_clusters))
set.seed(42)
if (cluster_meth == "hierarchical") {
  hc         <- hclust(dist(z_mat), method = "ward.D2")
  cluster_id <- cutree(hc, k = n_clusters)
} else {
  km         <- kmeans(z_mat, centers = n_clusters, nstart = 50, iter.max = 200)
  cluster_id <- km$cluster
}

cluster_df <- data.frame(
  peak_id    = rownames(z_mat),
  cluster_id = cluster_id,
  as.data.frame(z_mat),
  stringsAsFactors = FALSE
)
write.csv(cluster_df, file.path(output_dir, "da_peaks_temporal_clusters.csv"),
          row.names = FALSE)

cluster_summary <- cluster_df %>%
  group_by(cluster_id) %>%
  summarise(n_peaks = n(), .groups = "drop") %>%
  arrange(cluster_id)
write.csv(cluster_summary, file.path(output_dir, "cluster_summary.csv"), row.names = FALSE)
message("Cluster assignments 저장 완료")
print(cluster_summary)

# ── Temporal heatmap ───────────────────────────────────────────
message("\nTemporal heatmap 생성 중...")
row_order  <- order(cluster_id)
heat_mat   <- z_mat[row_order, , drop = FALSE]
ann_row    <- data.frame(
  Cluster = factor(paste0("C", cluster_id[row_order])),
  row.names = rownames(heat_mat)
)
cluster_colors <- setNames(
  colorRampPalette(brewer.pal(max(3, min(n_clusters, 9)), "Set1"))(n_clusters),
  paste0("C", seq_len(n_clusters))
)

png(file.path(output_dir, "temporal_heatmap.png"),
    width = max(900, length(time_order) * 80 + 300),
    height = min(4000, max(600, nrow(heat_mat) * 2 + 200)),
    res = 120)
pheatmap(
  heat_mat,
  cluster_rows   = FALSE,
  cluster_cols   = FALSE,
  annotation_row = ann_row,
  annotation_colors = list(Cluster = cluster_colors),
  color          = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
  show_rownames  = FALSE,
  show_colnames  = TRUE,
  fontsize        = 11,
  main = sprintf("Temporal DA Peak Heatmap (k=%d clusters, %d peaks)",
                 n_clusters, nrow(heat_mat))
)
dev.off()
message("Temporal heatmap 저장 완료")

# ── Trend plots per cluster ────────────────────────────────────
message("\nTrend plots 생성 중...")
# SE 계산을 위해 per-sample z-score도 필요 → 원 counts에서 샘플별로
norm_sub_ordered <- norm_sub[rownames(z_mat), , drop = FALSE]

for (k in seq_len(n_clusters)) {
  k_peaks <- rownames(z_mat)[cluster_id == k]
  if (length(k_peaks) == 0) next

  # 각 condition × sample의 z-score 값 (peak 단위로 정규화된 값)
  k_z    <- z_mat[k_peaks, , drop = FALSE]

  trend_df <- lapply(seq_along(time_order), function(i) {
    cond   <- time_order[i]
    vals   <- k_z[, i]
    data.frame(
      condition = cond,
      hours     = time_hours[i],
      mean_z    = mean(vals, na.rm = TRUE),
      se_z      = sd(vals, na.rm = TRUE) / sqrt(sum(!is.na(vals))),
      stringsAsFactors = FALSE
    )
  })
  trend_df <- bind_rows(trend_df)

  p_trend <- ggplot(trend_df, aes(x = hours, y = mean_z)) +
    geom_ribbon(aes(ymin = mean_z - se_z, ymax = mean_z + se_z),
                fill = cluster_colors[paste0("C", k)], alpha = 0.25) +
    geom_line(color = cluster_colors[paste0("C", k)], linewidth = 1.2) +
    geom_point(color = cluster_colors[paste0("C", k)], size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_x_continuous(breaks = time_hours, labels = time_order) +
    labs(
      title    = sprintf("Cluster %d (n=%d peaks)", k, length(k_peaks)),
      subtitle = "Mean ± SE of Z-score across timepoints",
      x        = "Timepoint",
      y        = "Mean Z-score"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(output_dir, sprintf("trend_cluster_%d.png", k)),
         p_trend, width = 7, height = 4.5, dpi = 150)
}
message(sprintf("  %d개 trend plots 저장 완료", n_clusters))

# ── GO enrichment per cluster ──────────────────────────────────
message("\n=== GO Enrichment per Cluster ===")
peak_info_path <- config$peak_info_path %||% NULL
peak_info_df   <- NULL
if (!is.null(peak_info_path) && file.exists(peak_info_path)) {
  peak_info_df <- read.csv(peak_info_path, stringsAsFactors = FALSE)
  message(sprintf("  peak_info 로드: %d peaks", nrow(peak_info_df)))
}

for (k in seq_len(n_clusters)) {
  k_peaks <- rownames(z_mat)[cluster_id == k]
  message(sprintf("\n  [Cluster %d, n=%d]", k, length(k_peaks)))
  if (length(k_peaks) < min_peaks_go) {
    message(sprintf("    peak 수 부족 (%d < %d), GO 건너뜀", length(k_peaks), min_peaks_go))
    next
  }
  if (is.null(peak_info_df)) {
    message("    peak_info 없어 GO 불가")
    next
  }
  genes <- peak_info_df %>%
    filter(peak_id %in% k_peaks, !is.na(gene_id), gene_id != "") %>%
    pull(gene_id) %>% unique()
  entrez <- tryCatch(
    bitr(genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org_db)$ENTREZID,
    error = function(e) character(0)
  )
  if (length(entrez) < 5) {
    message(sprintf("    Entrez gene 수 부족 (%d), 건너뜀", length(entrez)))
    next
  }
  res <- tryCatch(
    enrichGO(gene          = entrez,
             OrgDb         = org_db,
             ont           = go_ontology,
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05,
             qvalueCutoff  = 0.25,
             readable      = TRUE),
    error = function(e) NULL
  )
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    message("    → 유의한 GO term 없음")
    next
  }
  n_terms <- nrow(as.data.frame(res))
  message(sprintf("    → %d GO terms", n_terms))
  out_csv <- file.path(output_dir, "go_enrichment",
                       sprintf("cluster_%d_%s.csv", k, go_ontology))
  write.csv(as.data.frame(res), out_csv, row.names = FALSE)
}

message(sprintf("\nTime-course 분석 완료. 출력: %s", output_dir))
