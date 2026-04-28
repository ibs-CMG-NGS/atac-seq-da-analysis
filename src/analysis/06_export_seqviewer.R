#!/usr/bin/env Rscript
# 파일 경로: src/analysis/06_export_seqviewer.R
# cmg-seqviewer용 parquet + staging metadata JSON 생성 (per pair)
# RNA-Seq_DE_GO_analysis/06_export_seqviewer.R 과 스키마 통일
#
# Usage:
#   Rscript 06_export_seqviewer.R <config_path> <compare_group> <base_group> <pair_output_dir>
#
# Output (per pair):
#   {output_dir}/seqviewer/datasets/{pair}_da.parquet
#   {output_dir}/seqviewer/datasets/{pair}_go_kegg.parquet   (if enrichment exists)
#   {output_dir}/seqviewer/staging/{pair}_entries.json

suppressPackageStartupMessages({
  library(yaml)
  library(arrow)
  library(jsonlite)
  library(dplyr)
  library(openxlsx)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# ─────────────────────────────────────────────────────────────
# 1. 인자 파싱
# ─────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript 06_export_seqviewer.R <config_path> <compare_group> <base_group> <pair_output_dir>")
}
config_path     <- args[1]
compare_group   <- args[2]
base_group      <- args[3]
pair_output_dir <- args[4]   # e.g. output/2026-pym-mouse-atac/pairwise/GABA_vs_Veh

config <- yaml.load_file(config_path)

pair_label    <- paste0(compare_group, "_vs_", base_group)
output_dir    <- config$output_dir
seqviewer_dir <- file.path(output_dir, "seqviewer")
datasets_dir  <- file.path(seqviewer_dir, "datasets")
staging_dir   <- file.path(seqviewer_dir, "staging")

dir.create(datasets_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(staging_dir,  recursive = TRUE, showWarnings = FALSE)

cat(sprintf("\n=== SeqViewer Export: %s ===\n", pair_label))

# ─────────────────────────────────────────────────────────────
# 헬퍼 함수 (RNA-seq 06_export_seqviewer.R 동일)
# ─────────────────────────────────────────────────────────────
make_alias_slug <- function(alias, max_len = 40) {
  slug <- gsub("[^\\w\uac00-\ud7a3]+", "_", alias, perl = TRUE)
  slug <- gsub("^_|_$", "", slug)
  substr(slug, 1, max_len)
}

new_uuid <- function() {
  hex <- paste0(sample(c(0:9, letters[1:6]), 32, replace = TRUE), collapse = "")
  paste(
    substr(hex,  1,  8),
    substr(hex,  9, 12),
    paste0("4", substr(hex, 14, 16)),
    paste0(sample(c("8", "9", "a", "b"), 1), substr(hex, 18, 20)),
    substr(hex, 21, 32),
    sep = "-"
  )
}

write_parquet_dataset <- function(df, alias, datasets_dir) {
  uid      <- new_uuid()
  slug     <- make_alias_slug(alias)
  filename <- paste0(slug, ".parquet")
  # 재실행 시 중복 제거
  old_files <- list.files(datasets_dir, pattern = paste0("^", slug, "\\.parquet$"), full.names = TRUE)
  if (length(old_files) > 0) file.remove(old_files)
  write_parquet(df, file.path(datasets_dir, filename))
  list(uid = uid, filename = filename)
}

# ─────────────────────────────────────────────────────────────
# 2. DA 결과 → parquet
# ─────────────────────────────────────────────────────────────
# 입력: final_da_summary.xlsx (DA_Results 시트) 또는 final_da_results.csv
# RNA-seq의 DE_Results 시트와 동일한 컬럼 규격으로 표준화:
#   base_mean, log2fc, lfcse, adj_pvalue, peak_id, nearest_gene, annotation

da_xlsx <- file.path(pair_output_dir, "final_da_result.xlsx")
da_csv  <- file.path(pair_output_dir, "final_da_results.csv")

if (file.exists(da_xlsx)) {
  da_raw <- tryCatch(
    read.xlsx(da_xlsx, sheet = "DA_Results", check.names = FALSE),
    error = function(e) NULL
  )
} else {
  da_raw <- NULL
}

# xlsx가 없거나 비어있으면 CSV fallback
if (is.null(da_raw) || nrow(da_raw) == 0) {
  if (!file.exists(da_csv)) stop(paste("DA results not found:", da_csv))
  da_raw <- read.csv(da_csv, check.names = FALSE)
}

# ── 컬럼 표준화 (RNA-seq 규격과 통일) ──────────────────────
# RNA-seq:  gene_id, symbol, baseMean, log2FoldChange, lfcSE, padj, ...sample counts
# ATAC-seq: peak_id, nearest_gene, annotation, baseMean, log2FoldChange, lfcSE, padj

da_std <- da_raw %>%
  rename(any_of(c(
    base_mean    = "baseMean",
    log2fc       = "log2FoldChange",
    lfcse        = "lfcSE",
    adj_pvalue   = "padj",
    nearest_gene = "SYMBOL",         # ChIPseeker 출력
    nearest_gene = "geneSymbol",
    annotation   = "annotation",
    peak_chr     = "Chr",
    peak_start   = "Start",
    peak_end     = "End"
  )))

# peak_id 컬럼 보장
if (!"peak_id" %in% colnames(da_std) && "Geneid" %in% colnames(da_std)) {
  da_std <- da_std %>% rename(peak_id = Geneid)
}

# nearest_gene 없으면 peak_id 사용
if (!"nearest_gene" %in% colnames(da_std)) {
  da_std$nearest_gene <- da_std$peak_id %||% ""
}

# baseMean = 0 또는 NA 제거 (통계값 전부 NA)
if ("base_mean" %in% colnames(da_std)) {
  n_before <- nrow(da_std)
  da_std   <- da_std %>% filter(!is.na(base_mean) & base_mean > 0)
  cat(sprintf("Filtered %d zero-accessibility peaks; %d peaks retained\n",
              n_before - nrow(da_std), nrow(da_std)))
}

organism  <- config$species %||% ""
condition <- paste(compare_group, "vs", base_group)
da_alias  <- paste(condition, "DA")

da_info   <- write_parquet_dataset(da_std, da_alias, datasets_dir)

padj_cut <- config$da_analysis$padj_cutoff   %||% 0.05
lfc_cut  <- config$da_analysis$log2fc_cutoff %||% 1.0

sig_count <- sum(
  !is.na(da_std$adj_pvalue) & !is.na(da_std$log2fc) &
  da_std$adj_pvalue < padj_cut & abs(da_std$log2fc) >= lfc_cut,
  na.rm = TRUE
)

# ── staging entry (RNA-seq 스키마 + ATAC 전용 필드) ──────────
da_entry <- list(
  dataset_id           = da_info$uid,
  alias                = da_alias,
  original_filename    = da_info$filename,
  dataset_type         = "differential_accessibility",   # RNA-seq: "differential_expression"
  assay_type           = "atac-seq",                     # ATAC 전용
  experiment_condition = condition,
  organism             = organism,
  cell_type            = "",
  tissue               = "",
  timepoint            = "",
  row_count            = nrow(da_std),
  peak_count           = nrow(da_std),                   # RNA-seq: gene_count
  significant_peaks    = sig_count,                      # RNA-seq: significant_genes
  padj_cutoff          = padj_cut,
  log2fc_cutoff        = lfc_cut,
  import_date          = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  file_path            = da_info$filename,
  notes                = paste0(config$da_analysis$method %||% "DESeq2"),
  tags                 = as.list(c("DA", "ATAC", compare_group, base_group))
)
cat(sprintf("DA parquet saved: %s (%d peaks, %d significant)\n",
            da_info$filename, nrow(da_std), sig_count))

# ─────────────────────────────────────────────────────────────
# 3. GO/KEGG 결과 → parquet (RNA-seq 스키마와 완전 동일)
# ─────────────────────────────────────────────────────────────
go_entries <- list()

go_xlsx <- file.path(pair_output_dir, "final_go_result.xlsx")
if (file.exists(go_xlsx)) {
  sheets <- tryCatch(getSheetNames(go_xlsx), error = function(e) character(0))
  # No_Results 시트(DA peak 없는 경우 placeholder) 제외
  data_sheets <- sheets[sheets != "No_Results"]

  all_go <- lapply(data_sheets, function(sh) {
    df <- tryCatch(
      read.xlsx(go_xlsx, sheet = sh, check.names = FALSE),
      error = function(e) NULL
    )
    if (is.null(df) || nrow(df) == 0) return(NULL)

    # RNA-seq GO 시트와 동일한 컬럼 표준화
    df <- df %>%
      rename(any_of(c(
        term_id      = "GO ID",
        term_id      = "KEGG ID",
        description  = "GO Term",
        description  = "KEGG Pathway",
        gene_set     = "Gene Set",
        gene_ratio   = "Gene Ratio",
        bg_ratio     = "Background Ratio",
        pvalue       = "P-value",
        fdr          = "Adjusted P-value",
        qvalue       = "Q-value",
        gene_count   = "Gene Count",
        gene_symbols = "Gene Symbols"
      )))
    if (!"Ontology" %in% colnames(df)) df$Ontology <- "KEGG"
    df <- rename(df, any_of(c(ontology = "Ontology")))
    df
  })

  all_go <- Filter(Negate(is.null), all_go)

  if (length(all_go) > 0) {
    go_combined <- bind_rows(all_go)
    go_alias    <- paste(condition, "GO+KEGG")
    go_info     <- write_parquet_dataset(go_combined, go_alias, datasets_dir)

    go_entries[[1]] <- list(
      dataset_id           = go_info$uid,
      alias                = go_alias,
      original_filename    = go_info$filename,
      dataset_type         = "go_analysis",
      assay_type           = "atac-seq",
      experiment_condition = condition,
      organism             = organism,
      cell_type            = "",
      tissue               = "",
      timepoint            = "",
      row_count            = nrow(go_combined),
      peak_count           = 0L,
      significant_peaks    = 0L,
      import_date          = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      file_path            = go_info$filename,
      notes                = "clusterProfiler enrichGO + enrichKEGG (peak-associated genes)",
      tags                 = as.list(c("GO", "KEGG", "ATAC", compare_group, base_group))
    )
    cat(sprintf("GO+KEGG parquet saved: %s (%d terms)\n",
                go_info$filename, nrow(go_combined)))
  }
}

# ─────────────────────────────────────────────────────────────
# 4. staging JSON 저장 (RNA-seq와 동일 구조)
# ─────────────────────────────────────────────────────────────
all_entries  <- c(list(da_entry), go_entries)
staging_path <- file.path(staging_dir, paste0(pair_label, "_entries.json"))
write_json(all_entries, staging_path, pretty = TRUE, auto_unbox = TRUE)
cat(sprintf("Staging JSON saved: %s (%d entries)\n",
            staging_path, length(all_entries)))

# 완료 플래그
flag_path <- file.path(pair_output_dir, ".seqviewer_export_done.flag")
file.create(flag_path)
cat("seqviewer export done.\n")
