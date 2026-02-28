#!/usr/bin/env Rscript
# 05_motif_analysis.R
# ATAC-seq Motif Enrichment Analysis (HOMER)
# ATAC-seq 고유 분석 — RNA-seq에는 없는 스텝
#
# 요구사항: conda 환경에 homer 설치 필요
#   conda install -c bioconda homer
#   perl $(which configureHomer.pl) -install <genome>  (예: mm10, hg38)
#
# Usage:
#   Rscript src/analysis/05_motif_analysis.R <config> <compare> <base> <output_dir>

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(tibble)
  library(GenomicRanges)
  library(rtracklayer)
})

source("src/utils/load_data.R")

`%||%` <- function(a, b) if (!is.null(a)) a else b

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: Rscript 05_motif_analysis.R <config> <compare> <base> <output_dir>")
config_file   <- args[1]
compare_group <- args[2]
base_group    <- args[3]
output_dir    <- args[4]

config <- yaml::read_yaml(config_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("\n=== Motif Enrichment Analysis (HOMER) ===")
message(sprintf("Comparison: %s vs %s\n", compare_group, base_group))

motif_cfg    <- config$motif %||% list()
homer_genome <- motif_cfg$homer_genome %||% "mm10"
homer_opts   <- motif_cfg$homer_options %||% "-size 200 -mask"
gene_lists   <- motif_cfg$gene_lists   %||% c("up", "down")

padj_cutoff  <- config$da_analysis$padj_cutoff   %||% 0.05
lfc_cutoff   <- config$da_analysis$log2fc_cutoff  %||% 1.0

# ── HOMER 설치 확인 ────────────────────────────────────────────
homer_path <- Sys.which("findMotifsGenome.pl")
if (nchar(homer_path) == 0) {
  stop("HOMER를 찾을 수 없습니다. conda 환경에 homer가 설치되어 있는지 확인하세요.\n",
       "  conda install -c bioconda homer\n",
       "  perl $(which configureHomer.pl) -install ", homer_genome)
}
message("HOMER 경로: ", homer_path)

# ── Load DA results ────────────────────────────────────────────
da_dir  <- dirname(output_dir)   # output/pairwise/{pair}
res_df  <- load_da_results(file.path(da_dir, "final_da_results.csv"))
da_sets <- classify_da_peaks(res_df, padj_cutoff, lfc_cutoff)

# consensus peaks BED (background)
bg_bed <- config$consensus_peaks_path
if (!file.exists(bg_bed)) {
  stop("consensus_peaks_path 파일을 찾을 수 없습니다: ", bg_bed,
       "\nconfig.yml의 consensus_peaks_path를 확인하세요.")
}

# ── Peak ID → BED 변환 ─────────────────────────────────────────
peaks_to_bed <- function(peak_ids, out_file) {
  if (length(peak_ids) == 0) return(FALSE)

  # "chr:start-end" 파싱
  parts <- strsplit(peak_ids, "[:-]")
  valid <- sapply(parts, length) == 3
  peak_ids <- peak_ids[valid]
  parts    <- parts[valid]

  if (length(parts) == 0) {
    message("  BED 변환 가능한 peak이 없습니다.")
    return(FALSE)
  }

  bed_df <- data.frame(
    chr   = sapply(parts, `[[`, 1),
    start = as.integer(sapply(parts, `[[`, 2)),
    end   = as.integer(sapply(parts, `[[`, 3)),
    name  = peak_ids,
    score = 0,
    strand = "."
  )
  write.table(bed_df, out_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  message(sprintf("  BED 파일 저장: %s (%d peaks)", out_file, nrow(bed_df)))
  TRUE
}

# ── HOMER 실행 ─────────────────────────────────────────────────
run_homer <- function(target_bed, output_dir, genome, opts) {
  cmd <- sprintf(
    "findMotifsGenome.pl %s %s %s %s -bg %s -p 4",
    target_bed, genome, output_dir, opts, bg_bed
  )
  message("  실행: ", cmd)
  ret <- system(cmd)
  if (ret != 0) warning(sprintf("HOMER 종료 코드: %d", ret))
  ret
}

# ── 각 gene_list별 실행 ────────────────────────────────────────
for (gs_name in gene_lists) {
  peaks <- da_sets[[gs_name]]
  if (length(peaks) == 0) {
    message(sprintf("\n[%s] DA peaks 없음 — 건너뜀", gs_name))
    next
  }

  message(sprintf("\n[%s] %d peaks → HOMER 실행", gs_name, length(peaks)))

  gs_dir  <- file.path(output_dir, gs_name)
  bed_file <- file.path(output_dir, sprintf("%s_peaks.bed", gs_name))

  dir.create(gs_dir, recursive = TRUE, showWarnings = FALSE)

  ok <- peaks_to_bed(peaks, bed_file)
  if (!ok) next

  run_homer(bed_file, gs_dir, homer_genome, homer_opts)
}

# ── 결과 요약 ─────────────────────────────────────────────────
message("\n=== HOMER 결과 요약 ===")
for (gs_name in gene_lists) {
  gs_dir     <- file.path(output_dir, gs_name)
  known_html <- file.path(gs_dir, "knownResults.html")
  denovo_html<- file.path(gs_dir, "homerResults.html")

  if (file.exists(known_html)) {
    message(sprintf("[%s] Known motifs: %s", gs_name, known_html))
  }
  if (file.exists(denovo_html)) {
    message(sprintf("[%s] De novo motifs: %s", gs_name, denovo_html))
  }

  # Top known motifs 요약 CSV
  known_txt <- file.path(gs_dir, "knownResults.txt")
  if (file.exists(known_txt)) {
    motif_df <- tryCatch(
      read.table(known_txt, header = TRUE, sep = "\t", fill = TRUE,
                 check.names = FALSE, quote = ""),
      error = function(e) NULL
    )
    if (!is.null(motif_df) && nrow(motif_df) > 0) {
      top_motifs <- motif_df %>% slice_head(n = 20)
      write.csv(top_motifs,
                file.path(output_dir, sprintf("top_known_motifs_%s.csv", gs_name)),
                row.names = FALSE)
      message(sprintf("  Top 20 known motifs 저장: top_known_motifs_%s.csv", gs_name))
    }
  }
}

message("\nMotif analysis 완료. 출력: ", output_dir)
