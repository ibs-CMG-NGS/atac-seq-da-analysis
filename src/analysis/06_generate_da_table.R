#!/usr/bin/env Rscript
# 06_generate_da_table.R
# Publication-ready Excel output — two separate workbooks:
#   final_da_result.xlsx  : DA_Results  | Normalized_Counts | Significant_Peaks
#   final_go_result.xlsx  : UP_BP/CC/MF | DOWN_BP/CC/MF | TOTAL_BP/CC/MF | KEGG_UP/DOWN/TOTAL
#
# Usage:
#   Rscript src/analysis/06_generate_da_table.R <config> <output_dir>

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(openxlsx)
})

source("src/utils/load_data.R")

`%||%` <- function(a, b) if (!is.null(a)) a else b

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 06_generate_da_table.R <config> <output_dir>")
config_file <- args[1]
output_dir  <- args[2]

config <- yaml::read_yaml(config_file)

if (!isTRUE(config$export$export_to_excel %||% TRUE)) {
  message("export_to_excel=false — 건너뜀")
  quit(status = 0)
}

padj_cutoff <- config$da_analysis$padj_cutoff   %||% 0.05
lfc_cutoff  <- config$da_analysis$log2fc_cutoff  %||% 1.0
gene_lists  <- config$enrichment$gene_lists  %||% c("total", "up", "down")
go_ontos    <- config$enrichment$go_ontologies %||% c("BP", "CC", "MF")

safe_read_csv <- function(f) {
  if (!file.exists(f) || file.size(f) == 0) return(NULL)
  df <- tryCatch(read.csv(f), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df
}

# ── Shared styles ──────────────────────────────────────────────
make_header_style <- function() {
  createStyle(fontColour = "#FFFFFF", fgFill = "#2E4057",
              halign = "CENTER", textDecoration = "bold",
              border = "Bottom", borderColour = "#FFFFFF")
}
up_style   <- createStyle(fgFill = "#FFE0D9")
down_style <- createStyle(fgFill = "#D9E8FF")
sig_style  <- createStyle(fgFill = "#E8F5E9")

write_sheet <- function(wb, sheet_name, df, table_style = "TableStyleMedium9",
                        row_colors = NULL) {
  addWorksheet(wb, sheet_name)
  writeDataTable(wb, sheet_name, df, tableStyle = table_style)
  addStyle(wb, sheet_name, make_header_style(),
           rows = 1, cols = seq_len(ncol(df)))
  if (!is.null(row_colors)) {
    for (rc in row_colors) {
      addStyle(wb, sheet_name, rc$style, rows = rc$rows,
               cols = seq_len(ncol(df)), gridExpand = TRUE, stack = TRUE)
    }
  }
  setColWidths(wb, sheet_name, cols = seq_len(ncol(df)), widths = "auto")
}

message("\n=== Generating Excel Files ===\n")

# ── Load peak_info for annotation join ────────────────────────
# peak_info.csv: peak_id, chr, start, end, gene_name, gene_id,
#                annotation, distance_to_tss  (from nf-core/atacseq)
peak_info <- NULL
peak_info_path <- config$peak_info_path %||% NULL
if (!is.null(peak_info_path) && file.exists(peak_info_path)) {
  raw_info <- read.csv(peak_info_path, stringsAsFactors = FALSE)
  # 필요한 컬럼만 선택 (컬럼명이 다를 수 있으므로 유연하게)
  keep_cols <- intersect(
    c("peak_id", "chr", "start", "end", "gene_name", "gene_id",
      "annotation", "distance_to_tss"),
    colnames(raw_info)
  )
  peak_info <- raw_info[, keep_cols, drop = FALSE]
  message(sprintf("peak_info 로드: %d peaks, 컬럼: %s",
                  nrow(peak_info), paste(keep_cols, collapse = ", ")))
} else {
  message("peak_info_path 없음 — peak 좌표/유전자 정보 없이 진행")
}

# DA 결과에 peak_info join 후 컬럼 순서 정렬
# 컬럼 순서: peak_id | 좌표 | 유전자 정보 | DESeq2 통계 | direction
annotate_da <- function(df) {
  if (is.null(peak_info) || nrow(peak_info) == 0) return(df)
  df <- left_join(df, peak_info, by = "peak_id")

  # 컬럼 순서: peak_id, 좌표, 유전자, annotation, distance, DESeq2, direction
  front_cols <- intersect(
    c("peak_id", "chr", "start", "end",
      "gene_name", "gene_id", "annotation", "distance_to_tss"),
    colnames(df)
  )
  rest_cols <- setdiff(colnames(df), front_cols)
  df[, c(front_cols, rest_cols), drop = FALSE]
}

row_color_list <- function(df) {
  up_rows   <- which(df$direction == "up")   + 1
  down_rows <- which(df$direction == "down") + 1
  Filter(function(x) length(x$rows) > 0, list(
    list(style = up_style,   rows = up_rows),
    list(style = down_style, rows = down_rows)
  ))
}

# ══════════════════════════════════════════════════════════════
# File 1: final_da_result.xlsx
# ══════════════════════════════════════════════════════════════
wb_da <- createWorkbook()

# ── Sheet 1: DA_Results (all peaks) ───────────────────────────
da_file <- file.path(output_dir, "final_da_results.csv")
da <- NULL
if (file.exists(da_file)) {
  da <- read.csv(da_file) %>%
    mutate(direction = case_when(
      !is.na(padj) & padj < padj_cutoff & log2FoldChange >  lfc_cutoff ~ "up",
      !is.na(padj) & padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "down",
      TRUE ~ "ns"
    )) %>%
    annotate_da()

  rc <- row_color_list(da)
  write_sheet(wb_da, "DA_Results", da,
              table_style = "TableStyleMedium9",
              row_colors  = if (length(rc) > 0) rc else NULL)
  message(sprintf("Sheet 'DA_Results': %d peaks (%d up, %d down, %d ns)",
                  nrow(da),
                  sum(da$direction == "up"),
                  sum(da$direction == "down"),
                  sum(da$direction == "ns")))
}

# ── Sheet 2: Normalized_Counts (VST) ──────────────────────────
nc_file <- file.path(output_dir, "normalized_counts.csv")
if (file.exists(nc_file)) {
  nc <- read.csv(nc_file)
  # peak_id 뒤에 좌표 + 유전자 정보 삽입
  if (!is.null(peak_info)) {
    nc <- left_join(nc, peak_info[, intersect(
      c("peak_id", "chr", "start", "end", "gene_name", "gene_id"),
      colnames(peak_info))], by = "peak_id")
    front <- intersect(c("peak_id", "chr", "start", "end", "gene_name", "gene_id"),
                       colnames(nc))
    nc <- nc[, c(front, setdiff(colnames(nc), front))]
  }
  write_sheet(wb_da, "Normalized_Counts", nc, table_style = "TableStyleMedium2")
  message(sprintf("Sheet 'Normalized_Counts': %d peaks × %d samples",
                  nrow(nc), ncol(nc) - 1 - length(intersect(
                    c("chr","start","end","gene_name","gene_id"), colnames(nc)))))
}

# ── Sheet 3: Significant_Peaks (DA only) ──────────────────────
if (!is.null(da)) {
  sig <- da %>% filter(direction != "ns") %>% arrange(padj, desc(abs(log2FoldChange)))

  rc <- row_color_list(sig)
  write_sheet(wb_da, "Significant_Peaks", sig,
              table_style = "TableStyleMedium4",
              row_colors  = if (length(rc) > 0) rc else NULL)
  message(sprintf("Sheet 'Significant_Peaks': %d peaks", nrow(sig)))
}

out_da <- file.path(output_dir, "final_da_result.xlsx")
saveWorkbook(wb_da, out_da, overwrite = TRUE)
message("\nSaved: ", out_da)

# ══════════════════════════════════════════════════════════════
# File 2: final_go_result.xlsx
# ══════════════════════════════════════════════════════════════
wb_go <- createWorkbook()

# Tab order: UP first, then DOWN, then TOTAL (matches RNA-seq convention)
gs_order  <- c("up", "down", "total")
gs_labels <- c(up = "UP", down = "DOWN", total = "TOTAL")

go_table_style <- c(
  up    = "TableStyleMedium3",
  down  = "TableStyleMedium2",
  total = "TableStyleMedium9"
)

n_go_sheets <- 0L

for (gs_name in gs_order) {
  for (ont in go_ontos) {
    csv_file <- file.path(output_dir, sprintf("go_enrichment_%s_%s.csv", gs_name, ont))
    df <- safe_read_csv(csv_file)
    if (is.null(df)) next

    sheet_name <- sprintf("%s_%s", gs_labels[gs_name], ont)
    write_sheet(wb_go, sheet_name, df, table_style = go_table_style[gs_name])
    message(sprintf("Sheet '%s': %d terms", sheet_name, nrow(df)))
    n_go_sheets <- n_go_sheets + 1L
  }
}

# ── KEGG sheets ────────────────────────────────────────────────
kegg_table_style <- c(
  up    = "TableStyleLight9",
  down  = "TableStyleLight8",
  total = "TableStyleLight6"
)

for (gs_name in gs_order) {
  csv_file <- file.path(output_dir, sprintf("kegg_enrichment_%s.csv", gs_name))
  df <- safe_read_csv(csv_file)
  if (is.null(df)) next

  sheet_name <- sprintf("KEGG_%s", gs_labels[gs_name])
  write_sheet(wb_go, sheet_name, df, table_style = kegg_table_style[gs_name])
  message(sprintf("Sheet '%s': %d pathways", sheet_name, nrow(df)))
  n_go_sheets <- n_go_sheets + 1L
}

out_go <- file.path(output_dir, "final_go_result.xlsx")
if (n_go_sheets > 0L) {
  saveWorkbook(wb_go, out_go, overwrite = TRUE)
  message("\nSaved: ", out_go)
} else {
  # 빈 workbook — GO 결과 없음 (DA peak 없는 경우 등)
  addWorksheet(wb_go, "No_Results")
  df_empty <- data.frame(Note = "No significant DA peaks — GO/KEGG enrichment not performed.")
  writeData(wb_go, "No_Results", df_empty)
  saveWorkbook(wb_go, out_go, overwrite = TRUE)
  message("\nSaved (no enrichment results): ", out_go)
}
