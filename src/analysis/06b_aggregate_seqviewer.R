#!/usr/bin/env Rscript
# 파일 경로: src/analysis/06b_aggregate_seqviewer.R
# 목적: 모든 pairwise staging JSON을 하나의 CMG-SeqViewer import 패키지로 집계
# RNA-Seq_DE_GO_analysis/06b_aggregate_seqviewer.R 과 구조/출력 포맷 통일
#
# Usage:
#   Rscript 06b_aggregate_seqviewer.R <config_path> <output_dir>
#
# Output:
#   {output_dir}/seqviewer/seqviewer_manifest.json   ← SeqViewer import용
#   {output_dir}/seqviewer/seqviewer_import.zip       ← parquet + manifest 묶음

suppressPackageStartupMessages({
  library(yaml)
  library(jsonlite)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# ─────────────────────────────────────────────────────────────
# 인자 파싱
# ─────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript 06b_aggregate_seqviewer.R <config_path> <output_dir>")
}
config_path <- args[1]
output_dir  <- args[2]

config <- yaml.load_file(config_path)

seqviewer_dir <- file.path(output_dir, "seqviewer")
staging_dir   <- file.path(seqviewer_dir, "staging")
datasets_dir  <- file.path(seqviewer_dir, "datasets")
manifest_path <- file.path(seqviewer_dir, "seqviewer_manifest.json")
zip_path      <- file.path(seqviewer_dir, "seqviewer_import.zip")

cat("\n=== SeqViewer Aggregate Export ===\n")
cat(sprintf("Output dir: %s\n\n", output_dir))

# ─────────────────────────────────────────────────────────────
# staging JSON 수집
# ─────────────────────────────────────────────────────────────
staging_files <- list.files(staging_dir, pattern = "_entries\\.json$",
                            full.names = TRUE)

if (length(staging_files) == 0) {
  cat("[WARN] No staging JSON files found. Run 06_export_seqviewer.R first.\n")
  quit(status = 0)
}

cat(sprintf("Found %d staging file(s):\n", length(staging_files)))
for (f in staging_files) cat(sprintf("  %s\n", basename(f)))

all_entries <- lapply(staging_files, function(f) {
  tryCatch(
    fromJSON(f, simplifyVector = FALSE),
    error = function(e) { cat(sprintf("[WARN] Failed to read %s: %s\n", f, e$message)); list() }
  )
})
all_entries <- unlist(all_entries, recursive = FALSE)
cat(sprintf("\nTotal entries: %d\n", length(all_entries)))

# ─────────────────────────────────────────────────────────────
# 프로젝트 레벨 매니페스트 생성 (RNA-seq 버전과 동일 스키마)
# ─────────────────────────────────────────────────────────────
project_id <- basename(config$output_dir %||% output_dir)
pairs_raw  <- config$da_analysis$pairwise_comparisons %||% list()
pairs      <- vapply(pairs_raw, function(p) paste0(p[[1]], "_vs_", p[[2]]), character(1))

# 타입별 집계
da_entries  <- Filter(function(e) e$dataset_type == "differential_accessibility", all_entries)
go_entries  <- Filter(function(e) e$dataset_type == "go_analysis", all_entries)

total_peaks <- sum(vapply(da_entries, function(e) e$peak_count %||% 0L, integer(1)))
total_sig   <- sum(vapply(da_entries, function(e) e$significant_peaks %||% 0L, integer(1)))

manifest <- list(
  # ── 프로젝트 메타 (RNA-seq 버전과 동일 필드) ──
  project_id        = project_id,
  pipeline_type     = "atac-seq",          # RNA-seq: "rna-seq"
  assay_type        = "atac-seq",
  pipeline_version  = "1.0.0",
  generated_date    = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  species           = config$species %||% "unknown",
  comparisons       = as.list(pairs),
  # ── 집계 통계 ──
  summary = list(
    total_datasets       = length(all_entries),
    da_datasets          = length(da_entries),   # RNA-seq: de_datasets
    go_datasets          = length(go_entries),
    total_peaks_analyzed = total_peaks,           # RNA-seq: total_genes_analyzed
    total_significant    = total_sig,             # RNA-seq: total_significant_genes
    padj_cutoff          = config$da_analysis$padj_cutoff   %||% 0.05,
    log2fc_cutoff        = config$da_analysis$log2fc_cutoff %||% 1.0
  ),
  # ── 데이터셋 목록 ──
  datasets = all_entries,
  # ── 파일 경로 ──
  directory_structure = list(
    seqviewer_dir = seqviewer_dir,
    datasets_dir  = datasets_dir,
    staging_dir   = staging_dir
  )
)

write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
cat(sprintf("\nManifest saved: %s\n", manifest_path))

# ─────────────────────────────────────────────────────────────
# ZIP 패키지 생성 (manifest + 모든 parquet)
# ─────────────────────────────────────────────────────────────
parquet_files <- list.files(datasets_dir, pattern = "\\.parquet$", full.names = TRUE)

if (length(parquet_files) > 0) {
  files_to_zip <- c(manifest_path, parquet_files)
  # zip은 현재 디렉터리 기준 상대경로로 저장
  old_wd <- getwd()
  setwd(seqviewer_dir)

  zip_files_rel <- c(
    "seqviewer_manifest.json",
    file.path("datasets", basename(parquet_files))
  )

  tryCatch({
    zip(zipfile = zip_path, files = zip_files_rel)
    cat(sprintf("ZIP package saved: %s\n", zip_path))
    cat(sprintf("  Contents: 1 manifest + %d parquet file(s)\n", length(parquet_files)))
  }, error = function(e) {
    cat(sprintf("[WARN] ZIP creation failed: %s\n", e$message))
    cat("  Parquet files and manifest are available individually.\n")
  })

  setwd(old_wd)
} else {
  cat("[WARN] No parquet files found for zipping.\n")
}

# ─────────────────────────────────────────────────────────────
# 요약 출력
# ─────────────────────────────────────────────────────────────
cat(sprintf("\n=== Export Summary ===\n"))
cat(sprintf("Project:              %s\n",   project_id))
cat(sprintf("Pipeline type:        atac-seq\n"))
cat(sprintf("Comparisons:          %s\n",   paste(pairs, collapse = ", ")))
cat(sprintf("Total datasets:       %d\n",   length(all_entries)))
cat(sprintf("  DA datasets:        %d\n",   length(da_entries)))
cat(sprintf("  GO+KEGG datasets:   %d\n",   length(go_entries)))
cat(sprintf("Total peaks analyzed: %s\n",   format(total_peaks, big.mark = ",")))
cat(sprintf("Total significant:    %s\n",   format(total_sig,   big.mark = ",")))
cat(sprintf("\nImport files:\n"))
cat(sprintf("  Manifest: %s\n",  manifest_path))
if (file.exists(zip_path))
  cat(sprintf("  ZIP:      %s\n",  zip_path))
cat("\nDone.\n")
