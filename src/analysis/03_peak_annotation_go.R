#!/usr/bin/env Rscript
# 03_peak_annotation_go.R
# Peak Annotation (ChIPseeker) + GO/KEGG Enrichment (clusterProfiler)
# ATAC-seq 고유 스텝: DA peaks → 가장 가까운 유전자 → GO/KEGG
#
# Usage:
#   Rscript src/analysis/03_peak_annotation_go.R <config> <compare> <base> <output_dir>

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(tibble)
  library(ChIPseeker)
  library(GenomicRanges)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
})

source("src/utils/load_data.R")

`%||%` <- function(a, b) if (!is.null(a)) a else b

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: Rscript 03_peak_annotation_go.R <config> <compare> <base> <output_dir>")
config_file   <- args[1]
compare_group <- args[2]
base_group    <- args[3]
output_dir    <- args[4]

config <- yaml::read_yaml(config_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("\n=== Peak Annotation + GO/KEGG Enrichment ===")
message(sprintf("Comparison: %s vs %s\n", compare_group, base_group))

# ── Load database by species ───────────────────────────────────
species  <- config$species
db_cfg   <- config$databases[[species]]

org_db_name <- db_cfg$organism_db
kegg_code   <- db_cfg$kegg_code
txdb_name   <- db_cfg$txdb

message("Loading TxDb: ", txdb_name)
library(txdb_name, character.only = TRUE)
txdb <- get(txdb_name)
library(org_db_name, character.only = TRUE)
org_db <- get(org_db_name)

# ── Load DA results ────────────────────────────────────────────
padj_cutoff <- config$da_analysis$padj_cutoff   %||% 0.05
lfc_cutoff  <- config$da_analysis$log2fc_cutoff  %||% 1.0

res_df <- load_da_results(file.path(output_dir, "final_da_results.csv"))
da_sets <- classify_da_peaks(res_df, padj_cutoff, lfc_cutoff)

message(sprintf("DA peaks — total: %d, up: %d, down: %d",
                length(da_sets$total), length(da_sets$up), length(da_sets$down)))

# ── Parse peak IDs → GRanges ───────────────────────────────────
# peak_id 형식: "chr1:100000-100500" 또는 "consensus_peak_1"
# nf-core featureCounts의 Geneid는 "chr:start-end" 형식
parse_peaks_to_gr <- function(peak_ids) {
  # chr:start-end 형식 파싱
  parts <- strsplit(peak_ids, "[:-]")
  valid <- sapply(parts, length) == 3
  if (!all(valid)) {
    message(sprintf("  %d개 peak ID를 파싱할 수 없어 제외됨", sum(!valid)))
  }
  parts <- parts[valid]
  GRanges(
    seqnames = sapply(parts, `[[`, 1),
    ranges   = IRanges(
      start = as.integer(sapply(parts, `[[`, 2)),
      end   = as.integer(sapply(parts, `[[`, 3))
    ),
    peak_id  = peak_ids[valid]
  )
}

# ── Peak → Gene mapping via ChIPseeker ────────────────────────
tss_dist <- config$enrichment$tss_distance %||% 10000

annotate_peaks_to_genes <- function(peak_ids, txdb, tss_dist) {
  if (length(peak_ids) == 0) return(character(0))

  gr <- tryCatch(parse_peaks_to_gr(peak_ids), error = function(e) {
    message("  peak_id 파싱 오류: ", conditionMessage(e))
    return(GRanges())
  })
  if (length(gr) == 0) return(character(0))

  anno <- suppressMessages(
    annotatePeak(gr, tssRegion = c(-tss_dist, tss_dist),
                 TxDb = txdb, verbose = FALSE)
  )
  entrez_ids <- unique(na.omit(anno@anno$geneId))
  entrez_ids
}

message("Peak annotation 수행 중...")
gene_sets_entrez <- lapply(da_sets, annotate_peaks_to_genes,
                            txdb = txdb, tss_dist = tss_dist)

for (nm in names(gene_sets_entrez)) {
  message(sprintf("  [%s] → %d genes", nm, length(gene_sets_entrez[[nm]])))
}

# ── Peak Genomic Distribution ──────────────────────────────────
FEATURE_ORDER  <- c("Promoter", "5' UTR", "3' UTR", "Exon",
                    "Intron", "Downstream", "Distal Intergenic")
FEATURE_COLORS <- c(
  "Promoter"          = "#E57373",
  "5' UTR"            = "#FFB74D",
  "3' UTR"            = "#FFF176",
  "Exon"              = "#81C784",
  "Intron"            = "#64B5F6",
  "Downstream"        = "#CE93D8",
  "Distal Intergenic" = "#B0BEC5"
)

plot_peak_distribution <- function(peak_ids, label, txdb, tss_dist, out_dir) {
  if (length(peak_ids) == 0) {
    message(sprintf("  [%s] peaks 없음 — distribution 건너뜀", label))
    return(invisible(NULL))
  }

  gr <- tryCatch(parse_peaks_to_gr(peak_ids), error = function(e) GRanges())
  if (length(gr) == 0) return(invisible(NULL))

  anno <- suppressMessages(
    annotatePeak(gr, tssRegion = c(-tss_dist, tss_dist),
                 TxDb = txdb, verbose = FALSE)
  )

  write.csv(as.data.frame(anno@anno),
            file.path(out_dir, sprintf("peak_annotation_%s.csv", label)),
            row.names = FALSE)

  feat_cat <- dplyr::case_when(
    grepl("Promoter",   anno@anno$annotation) ~ "Promoter",
    grepl("5' UTR",     anno@anno$annotation) ~ "5' UTR",
    grepl("3' UTR",     anno@anno$annotation) ~ "3' UTR",
    grepl("Exon",       anno@anno$annotation) ~ "Exon",
    grepl("Intron",     anno@anno$annotation) ~ "Intron",
    grepl("Downstream", anno@anno$annotation) ~ "Downstream",
    TRUE                                      ~ "Distal Intergenic"
  )

  dist_df <- as.data.frame(table(feat_cat)) %>%
    dplyr::rename(Feature = feat_cat, Count = Freq) %>%
    dplyr::mutate(
      Pct     = Count / sum(Count) * 100,
      Feature = factor(Feature, levels = rev(FEATURE_ORDER))
    )

  p <- ggplot(dist_df, aes(x = Feature, y = Pct, fill = Feature)) +
    geom_col(alpha = 0.9, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", Pct, Count)),
              hjust = -0.05, size = 3) +
    coord_flip() +
    scale_fill_manual(values = FEATURE_COLORS, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(
      title = sprintf("Genomic Distribution — %s peaks (%s vs %s)",
                      label, compare_group, base_group),
      x = NULL, y = "Peaks (%)"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 11))

  ggsave(file.path(out_dir, sprintf("peak_distribution_%s.png", label)),
         p, width = 7, height = 4, dpi = 150)
  message(sprintf("  [%s] distribution plot 저장 (%d peaks)", label, length(peak_ids)))
}

message("\nPeak genomic distribution 시각화 중...")
for (gs_name in names(da_sets)) {
  plot_peak_distribution(da_sets[[gs_name]], gs_name, txdb, tss_dist, output_dir)
}

# ── GO Enrichment ─────────────────────────────────────────────
gene_lists  <- config$enrichment$gene_lists  %||% c("total", "up", "down")
go_ontos    <- config$enrichment$go_ontologies %||% c("BP", "CC", "MF")
pval_cut    <- config$enrichment$pvalue_cutoff %||% 0.05
qval_cut    <- config$enrichment$qvalue_cutoff %||% 0.25
min_gs      <- config$enrichment$min_gs_size   %||% 5
max_gs      <- config$enrichment$max_gs_size   %||% 500
min_count   <- config$enrichment$min_gene_count %||% 2
plot_top_n  <- config$enrichment$plot_top_n    %||% 15

run_go <- function(entrez_ids, ont, org_db,
                   pval_cut, qval_cut, min_gs, max_gs) {
  if (length(entrez_ids) < 5) {
    message(sprintf("    유전자 수 부족 (%d), GO 건너뜀", length(entrez_ids)))
    return(NULL)
  }
  tryCatch(
    enrichGO(gene         = entrez_ids,
             OrgDb        = org_db,
             ont          = ont,
             keyType      = "ENTREZID",
             pAdjustMethod = "BH",
             pvalueCutoff = pval_cut,
             qvalueCutoff = qval_cut,
             minGSSize    = min_gs,
             maxGSSize    = max_gs,
             readable     = TRUE),
    error = function(e) {
      message("    GO 오류: ", conditionMessage(e))
      NULL
    }
  )
}

message("\nGO Enrichment 실행 중...")
for (gs_name in gene_lists) {
  entrez <- gene_sets_entrez[[gs_name]]
  for (ont in go_ontos) {
    message(sprintf("  [%s / %s]", gs_name, ont))
    go_res <- run_go(entrez, ont, org_db, pval_cut, qval_cut, min_gs, max_gs)

    out_csv <- file.path(output_dir, sprintf("go_enrichment_%s_%s.csv", gs_name, ont))
    if (!is.null(go_res) && nrow(go_res@result) > 0) {
      res_filt <- go_res@result %>% filter(Count >= min_count)
      write.csv(res_filt, out_csv, row.names = FALSE)
      message(sprintf("    → %d terms 저장", nrow(res_filt)))

      # Dotplot
      n_terms <- sum(res_filt$p.adjust < pval_cut)
      if (n_terms > 0) {
        png(file.path(output_dir, sprintf("go_dotplot_%s_%s.png", gs_name, ont)),
            width = 900, height = 600, res = 120)
        print(dotplot(go_res, showCategory = plot_top_n,
                      title = sprintf("GO %s — %s peaks (%s vs %s)", ont, gs_name,
                                      compare_group, base_group)))
        dev.off()
      }
    } else {
      write.csv(data.frame(), out_csv, row.names = FALSE)
      message("    → 결과 없음")
    }
  }
}

# ── KEGG Enrichment ───────────────────────────────────────────
message("\nKEGG Enrichment 실행 중...")
for (gs_name in gene_lists) {
  entrez <- gene_sets_entrez[[gs_name]]
  message(sprintf("  [%s / KEGG]", gs_name))

  kegg_res <- tryCatch(
    enrichKEGG(gene         = entrez,
               organism     = kegg_code,
               pAdjustMethod = "BH",
               pvalueCutoff = pval_cut,
               qvalueCutoff = qval_cut,
               minGSSize    = min_gs,
               maxGSSize    = max_gs),
    error = function(e) {
      message("    KEGG 오류: ", conditionMessage(e))
      NULL
    }
  )

  out_csv <- file.path(output_dir, sprintf("kegg_enrichment_%s.csv", gs_name))
  if (!is.null(kegg_res) && nrow(kegg_res@result) > 0) {
    res_filt <- kegg_res@result %>% filter(Count >= min_count)
    write.csv(res_filt, out_csv, row.names = FALSE)
    message(sprintf("    → %d pathways 저장", nrow(res_filt)))

    n_terms <- sum(res_filt$p.adjust < pval_cut)
    if (n_terms > 0) {
      png(file.path(output_dir, sprintf("kegg_dotplot_%s.png", gs_name)),
          width = 900, height = 600, res = 120)
      print(dotplot(kegg_res, showCategory = plot_top_n,
                    title = sprintf("KEGG — %s peaks (%s vs %s)", gs_name,
                                    compare_group, base_group)))
      dev.off()
    }
  } else {
    write.csv(data.frame(), out_csv, row.names = FALSE)
    message("    → 결과 없음")
  }
}

message("\nPeak annotation + GO/KEGG 완료. 출력: ", output_dir)
