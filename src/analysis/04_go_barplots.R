#!/usr/bin/env Rscript
# 04_go_barplots.R
# GO/KEGG barplot 생성 (03_peak_annotation_go.R 출력 기반)
# RNA-Seq_DE_GO_analysis/src/analysis/04_generate_go_plots.R 패턴 계승
#
# Usage:
#   Rscript src/analysis/04_go_barplots.R <config> <output_dir>

suppressPackageStartupMessages({
  library(yaml)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 04_go_barplots.R <config> <output_dir>")
config_file <- args[1]
output_dir  <- args[2]

config <- yaml::read_yaml(config_file)

gene_lists  <- config$enrichment$gene_lists  %||% c("total", "up", "down")
namespaces  <- config$go_barplot$namespaces   %||% c("BP", "CC", "MF")
top_n       <- config$go_barplot$top_n        %||% 10
ns_colors   <- config$go_barplot$colors       %||%
  list(BP = "#E57373", CC = "#64B5F6", MF = "#81C784")
pval_cut    <- config$enrichment$pvalue_cutoff %||% 0.05

# 빈 CSV(결과 없음)를 안전하게 읽는 헬퍼
safe_read_csv <- function(f) {
  if (!file.exists(f) || file.size(f) == 0) return(NULL)
  df <- tryCatch(read.csv(f), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df
}

message("\n=== GO Barplots ===\n")

make_barplot <- function(df, title, fill_col = "#E57373", top_n = 10) {
  df_plot <- df %>%
    filter(p.adjust < pval_cut) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(Description = factor(Description,
                                levels = rev(Description)))
  if (nrow(df_plot) == 0) return(NULL)

  ggplot(df_plot, aes(x = Count, y = Description)) +
    geom_col(fill = fill_col, alpha = 0.85) +
    geom_text(aes(label = sprintf("p=%.2e", p.adjust)),
              hjust = -0.1, size = 3.2) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(title = title, x = "Gene Count", y = NULL) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12))
}

for (gs_name in gene_lists) {
  # 각 Ontology별 barplot
  for (ont in namespaces) {
    csv_file <- file.path(output_dir, sprintf("go_enrichment_%s_%s.csv", gs_name, ont))
    df <- safe_read_csv(csv_file)
    if (is.null(df)) next

    p <- make_barplot(df,
                      title   = sprintf("GO %s — %s peaks (top %d)", ont, gs_name, top_n),
                      fill_col = ns_colors[[ont]] %||% "steelblue",
                      top_n   = top_n)
    if (!is.null(p)) {
      ggsave(file.path(output_dir, sprintf("go_barplot_%s_%s.png", gs_name, ont)),
             p, width = 8, height = 5, dpi = 150)
    }
  }

  # 통합 barplot (BP + CC + MF)
  combined <- list()
  for (ont in namespaces) {
    csv_file <- file.path(output_dir, sprintf("go_enrichment_%s_%s.csv", gs_name, ont))
    df <- safe_read_csv(csv_file)
    if (is.null(df)) next
    df$ontology <- ont
    combined[[ont]] <- df %>% filter(p.adjust < pval_cut) %>% arrange(p.adjust) %>% slice_head(n = 5)
  }

  if (length(combined) > 0) {
    all_df <- bind_rows(combined) %>%
      arrange(p.adjust) %>%
      mutate(Description = factor(Description, levels = rev(unique(Description))))

    p_combined <- ggplot(all_df, aes(x = Count, y = Description, fill = ontology)) +
      geom_col(alpha = 0.85) +
      geom_text(aes(label = sprintf("p=%.2e", p.adjust)),
                hjust = -0.1, size = 3) +
      scale_fill_manual(values = unlist(ns_colors)) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
      labs(title = sprintf("GO Enrichment — %s peaks (top 5 per ontology)", gs_name),
           x = "Gene Count", y = NULL, fill = "Ontology") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size = 12))

    ggsave(file.path(output_dir, sprintf("go_barplot_%s_combined.png", gs_name)),
           p_combined, width = 9, height = max(4, nrow(all_df) * 0.35), dpi = 150)
    message(sprintf("[%s] Combined barplot 저장", gs_name))
  }
}

message("\nGO barplots 완료. 출력: ", output_dir)
