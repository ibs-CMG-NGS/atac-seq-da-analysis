#!/usr/bin/env Rscript
# 파일 경로: src/analysis/02d_generate_global_qc_report.R
# 목적: Global QC 플롯들을 모아 HTML 리포트 생성 (ATAC-seq 특화)
# RNA-Seq_DE_GO_analysis/02c_generate_global_qc_report.R 과 구조/CSS 통일
#
# Usage:
#   Rscript 02d_generate_global_qc_report.R \
#     --config config.yml \
#     --qc_plots_dir output/PROJECT/global_qc \
#     --output_file output/PROJECT/global_qc/global_qc_report.html

suppressPackageStartupMessages({
  library(yaml)
  library(optparse)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

option_list <- list(
  make_option(c("-c", "--config"), type = "character",
              help = "Path to config YAML file", metavar = "FILE"),
  make_option(c("-q", "--qc_plots_dir"), type = "character",
              help = "Directory containing QC plots", metavar = "DIR"),
  make_option(c("-o", "--output_file"), type = "character",
              help = "Output HTML file path", metavar = "FILE")
)
opt_parser <- OptionParser(option_list = option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$qc_plots_dir) || is.null(opt$output_file)) {
  print_help(opt_parser)
  stop("--qc_plots_dir and --output_file are required.", call. = FALSE)
}

config <- yaml.load_file(opt$config)
cat("\n=== ATAC-Seq Global QC Report Generation ===\n")
cat(paste0("QC plots directory: ", opt$qc_plots_dir, "\n"))
cat(paste0("Output HTML file:   ", opt$output_file, "\n\n"))

# ─────────────────────────────────────────────────────────────
# 메타데이터 및 카운트 로딩
# ─────────────────────────────────────────────────────────────
group_var <- config$da_analysis$group_variable %||% "condition"

meta_data <- tryCatch(
  read.csv(config$metadata_path, row.names = 1),
  error = function(e) { cat("[WARN] metadata not found\n"); data.frame() }
)

n_samples <- if (nrow(meta_data) > 0) nrow(meta_data) else "N/A"
groups    <- if (nrow(meta_data) > 0 && group_var %in% colnames(meta_data))
               unique(meta_data[[group_var]]) else character(0)

# featureCounts peak 수 추정
n_peaks <- tryCatch({
  as.integer(system(paste("wc -l <", config$count_matrix_path), intern = TRUE)) - 2L
}, error = function(e) "N/A")

# ─────────────────────────────────────────────────────────────
# 플롯 정의 (ATAC-seq QC 플롯 목록)
# ─────────────────────────────────────────────────────────────
plot_files <- c(
  "pca_plot.png",
  "pca_scree_plot.png",
  "sample_distance_heatmap.png",
  "count_distribution_boxplot.png",
  "dispersion_plot.png",
  "frip_distribution.png",
  "tss_enrichment.png",
  "peak_annotation_pie.png"
)

plot_titles <- c(
  "PCA — Sample Clustering",
  "PCA Scree Plot",
  "Sample Distance Heatmap",
  "Count Distribution Across Samples",
  "Dispersion Estimates",
  "FRiP Distribution",
  "TSS Enrichment Score",
  "Peak Genomic Annotation"
)

plot_descriptions <- c(
  "Principal Component Analysis showing the first two principal components of VST-normalized peak counts. Samples from the same condition should cluster together. Outliers may indicate sample swaps or technical failures.",
  "Variance explained by each principal component. The first PC(s) should capture the treatment effect if the signal is strong. A flat scree plot may indicate low signal or high sample variability.",
  "Hierarchical clustering of samples based on Euclidean distance of VST-normalized counts. Samples within the same condition should show shorter distances. Outliers appear as samples that do not cluster with their group.",
  "Distribution of normalized (VST) peak accessibility across all samples. Samples should show similar distributions after normalization. Large shifts may indicate library size or quality differences.",
  "Peak-wise dispersion estimates from DESeq2. Points should follow the fitted mean-dispersion trend. Genes with unusually high dispersion may reflect technical noise or genuine biological variability.",
  "Fraction of Reads in Peaks (FRiP) per sample. FRiP ≥ 0.20 is recommended for bulk ATAC-seq (ENCODE guideline). Low FRiP may indicate poor nucleosome positioning or library quality issues.",
  "TSS enrichment score per sample. A score ≥ 5.0 indicates good signal-to-noise. This metric reflects nucleosome-free regions around transcription start sites and is a key quality indicator for ATAC-seq.",
  "Genomic distribution of consensus peaks annotated by ChIPseeker. A high proportion of promoter/TSS peaks is expected in bulk ATAC-seq. Elevated intergenic fractions may indicate broad accessibility or technical artifacts."
)

# ─────────────────────────────────────────────────────────────
# HTML 생성 (RNA-seq 버전과 동일 CSS/구조)
# ─────────────────────────────────────────────────────────────
project_id <- basename(config$output_dir %||% "PROJECT")
da_method  <- config$da_analysis$method %||% "DESeq2"

html <- paste0('<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Global QC Report — ', project_id, '</title>
    <style>
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            line-height: 1.6;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }
        .header {
            background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }
        .header h1 { margin: 0 0 10px 0; font-size: 2.5em; }
        .header p  { margin: 5px 0; font-size: 1.1em; opacity: 0.9; }
        .summary {
            background: white;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .summary h2 {
            margin-top: 0;
            color: #11998e;
            border-bottom: 2px solid #11998e;
            padding-bottom: 10px;
        }
        .summary-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }
        .summary-item {
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            border-left: 4px solid #11998e;
        }
        .summary-item .label { font-size: 0.9em; color: #6c757d; margin-bottom: 5px; }
        .summary-item .value { font-size: 1.5em; font-weight: bold; color: #212529; }
        .plot-section {
            background: white;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .plot-section h2 {
            margin-top: 0;
            color: #11998e;
            border-bottom: 2px solid #11998e;
            padding-bottom: 10px;
        }
        .plot-description {
            background: #f8f9fa;
            padding: 15px;
            border-radius: 5px;
            margin: 15px 0;
            font-size: 0.95em;
            color: #495057;
            border-left: 4px solid #38ef7d;
        }
        .plot-container { text-align: center; margin: 20px 0; }
        .plot-container img {
            max-width: 100%;
            height: auto;
            border-radius: 5px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }
        .plot-missing {
            background: #fff3cd;
            color: #856404;
            padding: 20px;
            border-radius: 5px;
            text-align: center;
            border: 1px dashed #ffc107;
        }
        .footer {
            text-align: center;
            margin-top: 40px;
            padding: 20px;
            color: #6c757d;
            font-size: 0.9em;
        }
        table { width: 100%; border-collapse: collapse; margin: 15px 0; }
        th, td { padding: 12px; text-align: left; border-bottom: 1px solid #dee2e6; }
        th { background-color: #11998e; color: white; font-weight: 600; }
        tr:hover { background-color: #f8f9fa; }
        .badge-pass { background: #28a745; color: white; padding: 2px 8px; border-radius: 4px; font-size: 0.85em; }
        .badge-warn { background: #ffc107; color: #212529; padding: 2px 8px; border-radius: 4px; font-size: 0.85em; }
    </style>
</head>
<body>
<div class="header">
    <h1>&#x1F52C; ATAC-Seq Global QC Report</h1>
    <p>Differential Accessibility Analysis Pipeline</p>
    <p>Project: ', project_id, '</p>
    <p>Generated: ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '</p>
</div>

<div class="summary">
    <h2>&#x1F4CA; Dataset Summary</h2>
    <div class="summary-grid">
        <div class="summary-item">
            <div class="label">Total Samples</div>
            <div class="value">', n_samples, '</div>
        </div>
        <div class="summary-item">
            <div class="label">Consensus Peaks</div>
            <div class="value">', format(as.integer(n_peaks), big.mark=","), '</div>
        </div>
        <div class="summary-item">
            <div class="label">Conditions</div>
            <div class="value">', length(groups), '</div>
        </div>
        <div class="summary-item">
            <div class="label">DA Method</div>
            <div class="value">', da_method, '</div>
        </div>
    </div>
')

# 샘플 그룹 테이블
if (nrow(meta_data) > 0 && length(groups) > 0) {
  html <- paste0(html, '
    <h3 style="margin-top: 25px; color: #495057;">Sample Groups</h3>
    <table>
        <thead><tr><th>Condition</th><th>N Samples</th><th>Sample IDs</th></tr></thead>
        <tbody>')
  for (grp in groups) {
    samps <- rownames(meta_data)[meta_data[[group_var]] == grp]
    html <- paste0(html, '<tr><td><strong>', grp, '</strong></td><td>',
                   length(samps), '</td><td>', paste(samps, collapse = ", "), '</td></tr>')
  }
  html <- paste0(html, '</tbody></table>')
}

html <- paste0(html, '</div>\n')

# ── QC 가이드라인 (ATAC-seq 특화) ────────────────────────────
html <- paste0(html, '
<div class="summary">
    <h2>&#x2139; ATAC-Seq QC Guidelines</h2>
    <table>
        <thead><tr><th>Metric</th><th>Recommended Threshold</th><th>Reference</th></tr></thead>
        <tbody>
            <tr><td>FRiP (Fraction of Reads in Peaks)</td><td>≥ 0.20</td><td>ENCODE ATAC-seq Standards</td></tr>
            <tr><td>TSS Enrichment Score</td><td>≥ 5.0</td><td>ENCODE ATAC-seq Standards</td></tr>
            <tr><td>Unique aligned reads</td><td>≥ 25M (bulk)</td><td>Buenrostro et al. 2015</td></tr>
            <tr><td>Mitochondrial read fraction</td><td>≤ 0.20</td><td>ENCODE ATAC-seq Standards</td></tr>
            <tr><td>Non-redundant fraction (NRF)</td><td>≥ 0.70</td><td>ENCODE ATAC-seq Standards</td></tr>
        </tbody>
    </table>
</div>
')

# ── 플롯 섹션 ─────────────────────────────────────────────────
for (i in seq_along(plot_files)) {
  plot_path <- file.path(opt$qc_plots_dir, plot_files[i])
  img_block <- if (file.exists(plot_path)) {
    paste0('<div class="plot-container"><img src="', plot_files[i],
           '" alt="', plot_titles[i], '"></div>')
  } else {
    paste0('<div class="plot-missing">&#x26A0; Plot not generated: <code>',
           plot_files[i], '</code><br>',
           'This plot may not apply to the current dataset or was skipped.</div>')
  }

  html <- paste0(html, '
<div class="plot-section">
    <h2>', i, '. ', plot_titles[i], '</h2>
    <div class="plot-description">
        <strong>&#x1F4A1; Interpretation:</strong> ', plot_descriptions[i], '
    </div>
    ', img_block, '
</div>
')
}

# ── Footer ────────────────────────────────────────────────────
html <- paste0(html, '
<div class="footer">
    <p>Generated by ATAC-Seq DA Analysis Pipeline</p>
    <p>Configuration: ', opt$config, '</p>
</div>
</body>
</html>')

# ── Write ─────────────────────────────────────────────────────
if (!dir.exists(dirname(opt$output_file)))
  dir.create(dirname(opt$output_file), recursive = TRUE)

writeLines(html, opt$output_file)
cat(paste0("\n=== Global QC Report Generation Complete ===\n"))
cat(paste0("HTML report saved to: ", opt$output_file, "\n"))
cat(paste0("Open in browser: file://", normalizePath(opt$output_file), "\n\n"))
