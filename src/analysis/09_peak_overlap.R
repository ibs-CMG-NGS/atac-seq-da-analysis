#!/usr/bin/env Rscript
# 09_peak_overlap.R
# Multi-comparison DA peak overlap analysis
# UpSet plot (all comparisons) + Venn diagrams (subsets) + GO enrichment per overlap group
#
# Usage:
#   Rscript src/analysis/09_peak_overlap.R <config> <output_dir>

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ComplexUpset)
  library(VennDiagram)
  library(grid)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

source("src/utils/load_data.R")

`%||%` <- function(a, b) if (!is.null(a)) a else b

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 09_peak_overlap.R <config> <output_dir>")
config_file <- args[1]
output_dir  <- args[2]

config <- yaml::read_yaml(config_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "go_enrichment"), recursive = TRUE, showWarnings = FALSE)

message("\n=== Multi-Comparison Peak Overlap Analysis ===\n")

# ── Config ─────────────────────────────────────────────────────
overlap_cfg    <- config$peak_overlap %||% list()
padj_cutoff    <- config$da_analysis$padj_cutoff   %||% 0.05
lfc_cutoff     <- config$da_analysis$log2fc_cutoff %||% 1.0
min_peaks_go   <- overlap_cfg$min_peaks_for_go     %||% 20
go_ontology    <- overlap_cfg$go_ontology          %||% "BP"
venn_subsets   <- overlap_cfg$venn_subsets         %||% list()
pairs          <- config$da_analysis$pairwise_comparisons
output_base    <- config$output_dir

species        <- config$species
db_cfg         <- config$databases[[species]]
org_db_name    <- db_cfg$organism_db
kegg_code      <- db_cfg$kegg_code
library(org_db_name, character.only = TRUE)
org_db         <- get(org_db_name)

# ── Load DA results ────────────────────────────────────────────
message("DA results 로딩...")
pair_names <- sapply(pairs, function(p) sprintf("%s_vs_%s", p[[1]], p[[2]]))
da_list <- lapply(pair_names, function(pair) {
  path <- file.path(output_base, "pairwise", pair, "final_da_results.csv")
  if (!file.exists(path)) stop(sprintf("파일 없음: %s", path))
  df <- load_da_results(path)
  df$pair <- pair
  df
})
names(da_list) <- pair_names
message(sprintf("  %d개 비교 로드됨", length(da_list)))

# ── Peak membership matrix ─────────────────────────────────────
message("Peak membership matrix 생성 중...")
all_peaks <- unique(unlist(lapply(da_list, function(df) df$peak_id)))

sig_sets <- lapply(da_list, function(df) {
  df %>%
    filter(!is.na(padj), padj < padj_cutoff, abs(log2FoldChange) >= lfc_cutoff) %>%
    pull(peak_id)
})

membership <- data.frame(peak_id = all_peaks, stringsAsFactors = FALSE)
for (pair in pair_names) {
  membership[[pair]] <- membership$peak_id %in% sig_sets[[pair]]
}

write.csv(membership, file.path(output_dir, "peak_membership_matrix.csv"),
          row.names = FALSE)
message(sprintf("  총 peaks: %d (DA인 것: %d)",
                nrow(membership),
                sum(rowSums(membership[, -1]) > 0)))

# ── Overlap summary ────────────────────────────────────────────
overlap_df <- membership %>%
  mutate(pattern = apply(.[, pair_names, drop = FALSE], 1,
                         function(r) paste(pair_names[r], collapse = "|"))) %>%
  filter(pattern != "") %>%
  group_by(pattern) %>%
  summarise(n_peaks = n(), .groups = "drop") %>%
  arrange(desc(n_peaks))

write.csv(overlap_df, file.path(output_dir, "overlap_summary.csv"), row.names = FALSE)
message("Overlap summary 저장 완료")
print(overlap_df)

# ── UpSet plot ─────────────────────────────────────────────────
message("\nUpSet plot 생성 중...")
upset_data <- membership %>% dplyr::select(all_of(pair_names))
# ComplexUpset requires data frame with logical columns
upset_data[] <- lapply(upset_data, as.logical)

p_upset <- upset(
  upset_data,
  intersect  = pair_names,
  name       = "Comparison",
  width_ratio = 0.2,
  set_sizes  = upset_set_size(
    geom = geom_bar(fill = "#5B8DB8", width = 0.6)
  ) + ylab("Set size"),
  base_annotations = list(
    "Intersection size" = intersection_size(
      counts = TRUE,
      text   = list(size = 3)
    )
  )
) +
  labs(title = "DA Peak Overlap Across Comparisons") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(output_dir, "upset_plot.png"),
       p_upset, width = max(8, length(pair_names) * 1.5), height = 7, dpi = 150)
message("UpSet plot 저장 완료")

# ── Venn diagrams (subsets) ────────────────────────────────────
if (length(venn_subsets) > 0) {
  message("\nVenn diagram 생성 중...")
  for (subset in venn_subsets) {
    subset_pairs <- unlist(subset)
    if (!all(subset_pairs %in% pair_names)) {
      message(sprintf("  경고: subset에 없는 비교명 포함 (%s), 건너뜀",
                      paste(subset_pairs[!subset_pairs %in% pair_names], collapse = ", ")))
      next
    }
    # 3개 이하만 VennDiagram으로 그릴 수 있음
    n <- length(subset_pairs)
    if (n > 4) {
      message(sprintf("  경고: Venn은 최대 4-way, %d-way 건너뜀", n))
      next
    }
    sets   <- lapply(subset_pairs, function(p) sig_sets[[p]])
    labels <- gsub("_vs_", "\nvs\n", subset_pairs)
    safe_name <- gsub("[^a-zA-Z0-9_]", "_", paste(subset_pairs, collapse = "_"))
    out_venn  <- file.path(output_dir, sprintf("venn_%s.png", safe_name))

    cols <- c("#E57373", "#64B5F6", "#81C784", "#FFB74D")[seq_len(n)]
    venn.diagram(
      x          = sets,
      category.names = labels,
      filename   = out_venn,
      output     = TRUE,
      imagetype  = "png",
      height     = 1800, width = 1800, resolution = 150,
      col        = "white",
      fill       = cols,
      alpha      = 0.5,
      cex        = 1.2,
      fontface   = "bold",
      cat.cex    = 1.0,
      cat.fontface = "bold",
      margin     = 0.1,
      main       = sprintf("%d-way Venn Diagram", n),
      main.cex   = 1.3,
      main.fontface = "bold"
    )
    message(sprintf("  Venn (%d-way) 저장: %s", n, basename(out_venn)))
  }
}

# ── GO enrichment per overlap group ───────────────────────────
message("\n=== GO Enrichment per Overlap Group ===")
peak_info_path <- config$peak_info_path %||% NULL
peak_info_df   <- NULL
if (!is.null(peak_info_path) && file.exists(peak_info_path)) {
  peak_info_df <- read.csv(peak_info_path, stringsAsFactors = FALSE)
  message(sprintf("  peak_info 로드: %d peaks", nrow(peak_info_df)))
} else {
  # count_matrix 옆에서 자동 탐색
  candidates <- Sys.glob("/home/ngs/data/ygkim/2026/count_matrix/*_peak_info.csv")
  for (p in candidates) {
    if (grepl(basename(config$output_dir), p) ||
        grepl("ben|pym", p)) {
      peak_info_df <- read.csv(p, stringsAsFactors = FALSE)
      message(sprintf("  peak_info 자동 탐색: %s", p))
      break
    }
  }
}

run_go_for_group <- function(peak_ids, group_label) {
  if (length(peak_ids) < min_peaks_go) {
    message(sprintf("    [%s] peak 수 부족 (%d < %d), 건너뜀",
                    group_label, length(peak_ids), min_peaks_go))
    return(NULL)
  }
  # peak_id → gene mapping
  if (!is.null(peak_info_df)) {
    genes <- peak_info_df %>%
      filter(peak_id %in% peak_ids, !is.na(gene_id), gene_id != "") %>%
      pull(gene_id) %>% unique()
  } else {
    message(sprintf("    [%s] peak_info 없어 GO 불가", group_label))
    return(NULL)
  }
  if (length(genes) < 5) {
    message(sprintf("    [%s] gene 수 부족 (%d), 건너뜀", group_label, length(genes)))
    return(NULL)
  }
  # Ensembl → Entrez 변환
  entrez <- tryCatch(
    bitr(genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org_db)$ENTREZID,
    error = function(e) character(0)
  )
  if (length(entrez) < 5) {
    message(sprintf("    [%s] Entrez 변환 후 gene 수 부족 (%d), 건너뜀",
                    group_label, length(entrez)))
    return(NULL)
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
    message(sprintf("    [%s] → 유의한 GO term 없음", group_label))
    return(NULL)
  }
  message(sprintf("    [%s] → %d GO terms", group_label, nrow(as.data.frame(res))))
  as.data.frame(res)
}

# 각 overlap 패턴별 GO
for (i in seq_len(nrow(overlap_df))) {
  pattern   <- overlap_df$pattern[i]
  n_peaks   <- overlap_df$n_peaks[i]
  peak_ids  <- membership %>%
    filter(apply(.[, pair_names, drop = FALSE], 1,
                 function(r) paste(pair_names[r], collapse = "|")) == pattern) %>%
    pull(peak_id)

  safe_pat  <- gsub("[^a-zA-Z0-9]", "_", pattern)
  safe_pat  <- substr(safe_pat, 1, 80)
  go_res    <- run_go_for_group(peak_ids, safe_pat)
  if (!is.null(go_res)) {
    out_csv <- file.path(output_dir, "go_enrichment",
                         sprintf("group_%s_%s.csv", safe_pat, go_ontology))
    write.csv(go_res, out_csv, row.names = FALSE)
  }
}

message(sprintf("\nPeak overlap 분석 완료. 출력: %s", output_dir))
