#!/usr/bin/env Rscript
# 06_generate_da_table.R
# Publication-ready Excel summary (DA results + GO/KEGG 통합)
# RNA-Seq_DE_GO_analysis/src/analysis/05_generate_go_table.R 패턴 계승
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

message("\n=== Generating Summary Excel ===\n")

# ── Styles ────────────────────────────────────────────────────
wb <- createWorkbook()
header_style <- createStyle(fontColour = "#FFFFFF", fgFill = "#2E4057",
                             halign = "CENTER", textDecoration = "bold",
                             border = "Bottom", borderColour = "#FFFFFF")
up_style   <- createStyle(fgFill = "#FFE0D9")
down_style <- createStyle(fgFill = "#D9E8FF")

# ── Sheet 1: DA Results ───────────────────────────────────────
da_file <- file.path(output_dir, "final_da_results.csv")
if (file.exists(da_file)) {
  da <- read.csv(da_file) %>%
    mutate(direction = case_when(
      !is.na(padj) & padj < padj_cutoff & log2FoldChange >  lfc_cutoff ~ "up",
      !is.na(padj) & padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "down",
      TRUE ~ "ns"
    ))

  addWorksheet(wb, "DA_Results")
  writeDataTable(wb, "DA_Results", da, tableStyle = "TableStyleMedium9")
  addStyle(wb, "DA_Results", header_style, rows = 1, cols = seq_len(ncol(da)))

  # 색상: up = 연빨강, down = 연파랑
  up_rows   <- which(da$direction == "up")   + 1
  down_rows <- which(da$direction == "down") + 1
  if (length(up_rows) > 0)
    addStyle(wb, "DA_Results", up_style,   rows = up_rows,   cols = seq_len(ncol(da)), gridExpand = TRUE, stack = TRUE)
  if (length(down_rows) > 0)
    addStyle(wb, "DA_Results", down_style, rows = down_rows, cols = seq_len(ncol(da)), gridExpand = TRUE, stack = TRUE)

  setColWidths(wb, "DA_Results", cols = seq_len(ncol(da)), widths = "auto")
  message(sprintf("Sheet 'DA_Results': %d peaks (%d up, %d down)",
                  nrow(da), sum(da$direction == "up"), sum(da$direction == "down")))
}

# ── Sheets: GO per ontology × geneset ────────────────────────
for (ont in go_ontos) {
  for (gs_name in gene_lists) {
    csv_file <- file.path(output_dir, sprintf("go_enrichment_%s_%s.csv", gs_name, ont))
    if (!file.exists(csv_file)) next
    df <- read.csv(csv_file)
    if (nrow(df) == 0) next

    sheet_name <- sprintf("GO_%s_%s", ont, gs_name)
    addWorksheet(wb, sheet_name)
    writeDataTable(wb, sheet_name, df, tableStyle = "TableStyleMedium2")
    addStyle(wb, sheet_name, header_style, rows = 1, cols = seq_len(ncol(df)))
    setColWidths(wb, sheet_name, cols = seq_len(ncol(df)), widths = "auto")
    message(sprintf("Sheet '%s': %d terms", sheet_name, nrow(df)))
  }
}

# ── Sheets: KEGG per geneset ──────────────────────────────────
for (gs_name in gene_lists) {
  csv_file <- file.path(output_dir, sprintf("kegg_enrichment_%s.csv", gs_name))
  if (!file.exists(csv_file)) next
  df <- read.csv(csv_file)
  if (nrow(df) == 0) next

  sheet_name <- sprintf("KEGG_%s", gs_name)
  addWorksheet(wb, sheet_name)
  writeDataTable(wb, sheet_name, df, tableStyle = "TableStyleMedium3")
  addStyle(wb, sheet_name, header_style, rows = 1, cols = seq_len(ncol(df)))
  setColWidths(wb, sheet_name, cols = seq_len(ncol(df)), widths = "auto")
  message(sprintf("Sheet '%s': %d pathways", sheet_name, nrow(df)))
}

# ── Save ──────────────────────────────────────────────────────
out_xlsx <- file.path(output_dir, "final_da_summary.xlsx")
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message("\nExcel 저장 완료: ", out_xlsx)
