#!/usr/bin/env Rscript
# 02c_pairwise_plots.R
# Pairwise QC plots: volcano, MA plot
# RNA-Seq_DE_GO_analysis/src/analysis/02_generate_plots.R 패턴 계승
#
# Usage:
#   Rscript src/analysis/02c_pairwise_plots.R <config> <compare> <base> <output_dir>

suppressPackageStartupMessages({
  library(yaml)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

source("src/utils/load_data.R")

`%||%` <- function(a, b) if (!is.null(a)) a else b

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: Rscript 02c_pairwise_plots.R <config> <compare> <base> <output_dir>")
config_file   <- args[1]
compare_group <- args[2]
base_group    <- args[3]
output_dir    <- args[4]

config <- yaml::read_yaml(config_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

padj_cutoff <- config$da_analysis$padj_cutoff   %||% 0.05
lfc_cutoff  <- config$da_analysis$log2fc_cutoff  %||% 1.0

# Plot aesthetics
aes_cfg  <- config$plot_aesthetics$volcano %||% list()
up_col   <- aes_cfg$up_color    %||% "#FF5733"
down_col <- aes_cfg$down_color  %||% "#3375FF"
base_col <- aes_cfg$base_color  %||% "grey70"
pt_size  <- aes_cfg$point_size  %||% 2.0
top_n    <- aes_cfg$label_top_n %||% 10
fs       <- aes_cfg$base_font_size %||% 14

# ── Load DA results ────────────────────────────────────────────
res_df <- load_da_results(file.path(output_dir, "final_da_results.csv"))

res_df <- res_df %>%
  mutate(
    neg_log10_padj = -log10(pmax(padj, 1e-300)),
    direction = case_when(
      !is.na(padj) & padj < padj_cutoff & log2FoldChange >  lfc_cutoff ~ "up",
      !is.na(padj) & padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "down",
      TRUE ~ "ns"
    )
  )

n_up   <- sum(res_df$direction == "up",   na.rm = TRUE)
n_down <- sum(res_df$direction == "down",  na.rm = TRUE)

# ── Volcano plot ───────────────────────────────────────────────
label_df <- res_df %>%
  filter(direction != "ns") %>%
  arrange(padj) %>%
  slice_head(n = top_n)

p_vol <- ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj, color = direction)) +
  geom_point(size = pt_size, alpha = 0.6) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(padj_cutoff),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  scale_color_manual(
    values = c(up = up_col, down = down_col, ns = base_col),
    labels = c(
      up   = sprintf("More open (%s)  n=%d", compare_group, n_up),
      down = sprintf("Less open (%s)  n=%d", compare_group, n_down),
      ns   = "Not significant"
    )
  ) +
  { if (top_n > 0 && nrow(label_df) > 0)
      geom_text_repel(data = label_df, aes(label = peak_id),
                      size = 3, max.overlaps = 20, show.legend = FALSE)
  } +
  labs(
    title    = sprintf("Volcano Plot: %s vs %s", compare_group, base_group),
    x        = expression(log[2]~"Fold Change"),
    y        = expression(-log[10]~"(adjusted p-value)"),
    color    = NULL
  ) +
  theme_bw(base_size = fs) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "volcano_plot.png"), p_vol,
       width = 8, height = 6, dpi = 150)
message("Volcano plot 저장: ", file.path(output_dir, "volcano_plot.png"))

# ── MA plot ────────────────────────────────────────────────────
ylim_val <- config$qc_plots$ma_plot_ylim %||% c(-5, 5)

p_ma <- ggplot(res_df %>% filter(!is.na(log2FoldChange)),
               aes(x = log2(baseMean + 1), y = log2FoldChange, color = direction)) +
  geom_point(size = pt_size * 0.8, alpha = 0.5) +
  geom_hline(yintercept = c(-lfc_cutoff, 0, lfc_cutoff),
             linetype   = c("dashed", "solid", "dashed"),
             color      = c("grey40", "black", "grey40"),
             linewidth  = 0.5) +
  scale_color_manual(values = c(up = up_col, down = down_col, ns = base_col),
                     guide = "none") +
  coord_cartesian(ylim = ylim_val) +
  labs(
    title = sprintf("MA Plot: %s vs %s", compare_group, base_group),
    x     = expression(log[2]~"Mean Accessibility"),
    y     = expression(log[2]~"Fold Change")
  ) +
  theme_bw(base_size = fs)

ggsave(file.path(output_dir, "ma_plot.png"), p_ma,
       width = 7, height = 5, dpi = 150)
message("MA plot 저장: ", file.path(output_dir, "ma_plot.png"))

message("\nPairwise plots 완료.")
