#!/usr/bin/env Rscript
# 파일 경로: src/analysis/07_generate_summary_report.R
# 목적: DA 파이프라인 완료 후 모든 pairwise 비교군 결과를 통합 HTML 요약 리포트로 생성
# RNA-Seq_DE_GO_analysis/07_generate_summary_report.R 과 구조/CSS 통일
#
# 사용법:
#   Rscript 07_generate_summary_report.R \
#     --config configs/config_PROJECT.yml \
#     --output-dir output/PROJECT \
#     --output output/PROJECT/summary_report.html

# --- 1. Setup ---
suppressPackageStartupMessages({
  library(yaml)
  library(optparse)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

option_list <- list(
  make_option(c("-c", "--config"), type = "character",
              help = "Path to the config YAML file", metavar = "FILE"),
  make_option(c("-d", "--output-dir"), type = "character",
              help = "Pipeline output directory (e.g. output/PROJECT)", metavar = "DIR"),
  make_option(c("-o", "--output"), type = "character",
              help = "Output HTML file path", metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$config) || is.null(opt$`output-dir`) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("--config, --output-dir and --output are required.", call. = FALSE)
}

config      <- yaml.load_file(opt$config)
output_dir  <- normalizePath(opt$`output-dir`, mustWork = TRUE)
output_file <- opt$output

cat("\n=== DA Summary Report Generation ===\n")
cat("Config:     ", opt$config,   "\n")
cat("Output dir: ", output_dir,   "\n")
cat("HTML:       ", output_file,  "\n\n")

# --- 2. Derive parameters ---
get_pairs <- function(cfg) {
  pairs <- character(0)
  for (p in cfg$da_analysis$pairwise_comparisons) {
    pairs <- c(pairs, paste0(p[[1]], "_vs_", p[[2]]))
  }
  pairs
}

PAIRS      <- get_pairs(config)
padj_cut   <- config$da_analysis$padj_cutoff   %||% 0.05
lfc_cut    <- config$da_analysis$log2fc_cutoff %||% 1.0
project_id <- basename(config$output_dir %||% output_dir)
method     <- config$da_analysis$method        %||% "DESeq2"
species    <- config$species                   %||% "unknown"
group_var  <- config$da_analysis$group_variable %||% "condition"
lfc_shrink <- config$da_analysis$lfc_shrinkage  %||% "apeglm"
tss_dist   <- config$enrichment$tss_distance    %||% 10000
go_ontos   <- config$enrichment$go_ontologies   %||% c("BP", "CC", "MF")
kegg_code  <- config$databases[[species]]$kegg_code %||% NULL

cat(sprintf("Project: %s | Pairs: %d | padj < %.2f | |log2FC| > %.1f\n\n",
            project_id, length(PAIRS), padj_cut, lfc_cut))

# --- Load metadata ---
condition_samples <- list()
n_total_samples   <- NA_integer_

meta_path <- config$metadata_path %||% NULL
if (!is.null(meta_path) && file.exists(meta_path)) {
  tryCatch({
    meta <- read.csv(meta_path, stringsAsFactors = FALSE)
    id_col <- intersect(c("sample_id", "sample", "SampleID", "ID"), colnames(meta))[1]
    if (is.na(id_col)) id_col <- colnames(meta)[1]
    if (group_var %in% colnames(meta)) {
      condition_samples <- split(meta[[id_col]], meta[[group_var]])
      n_total_samples   <- nrow(meta)
      cat(sprintf("  Metadata loaded: %d samples, %d conditions\n\n",
                  n_total_samples, length(condition_samples)))
    }
  }, error = function(e) cat(sprintf("  [WARN] Could not read metadata: %s\n\n", e$message)))
}

# --- 3. Collect per-pair results ---
pair_data <- list()

for (pair in PAIRS) {
  pair_dir <- file.path(output_dir, "pairwise", pair)
  da_file  <- file.path(pair_dir, "final_da_results.csv")

  if (!file.exists(da_file)) {
    cat(sprintf("  [WARN] %s: final_da_results.csv not found, skipping.\n", pair))
    next
  }
  da <- read.csv(da_file, check.names = FALSE)

  # column detection
  padj_col <- intersect(c("padj", "FDR", "adj.P.Val"), colnames(da))[1]
  lfc_col  <- intersect(c("log2FoldChange", "logFC"), colnames(da))[1]
  peak_col <- intersect(c("peak_id", "PeakID", "interval", "Interval"), colnames(da))[1]
  if (is.na(peak_col)) peak_col <- colnames(da)[1]

  sig  <- da[!is.na(da[[padj_col]]) & da[[padj_col]] < padj_cut &
               abs(da[[lfc_col]]) > lfc_cut, ]
  up   <- sig[sig[[lfc_col]] > 0, ]
  down <- sig[sig[[lfc_col]] < 0, ]

  # Top peaks sorted by |log2FC|
  top_up   <- head(up[order(-abs(up[[lfc_col]])), ],   10)
  top_down <- head(down[order(-abs(down[[lfc_col]])), ], 10)

  # GO BP enrichment (up/down)
  read_go <- function(geneset) {
    f <- file.path(pair_dir, paste0("go_enrichment_", geneset, "_BP.csv"))
    if (!file.exists(f)) return(data.frame())
    tryCatch(read.csv(f, check.names = FALSE), error = function(e) data.frame())
  }
  go_up   <- read_go("up")
  go_down <- read_go("down")

  # Volcano plot (relative from output_dir)
  volcano_rel <- file.path("pairwise", pair, "volcano_plot.png")
  ma_rel      <- file.path("pairwise", pair, "ma_plot.png")

  pair_data[[pair]] <- list(
    pair        = pair,
    n_up        = nrow(up),
    n_down      = nrow(down),
    n_total     = nrow(sig),
    n_peaks     = nrow(da),
    top_up      = top_up,
    top_down    = top_down,
    go_up       = go_up,
    go_down     = go_down,
    peak_col    = peak_col,
    lfc_col     = lfc_col,
    padj_col    = padj_col,
    volcano_rel = volcano_rel,
    ma_rel      = ma_rel
  )

  cat(sprintf("  %s: %d peaks total | %d DA (up: %d / down: %d)\n",
              pair, nrow(da), nrow(sig), nrow(up), nrow(down)))
}

cat("\n")

# --- 4. HTML helpers ---

# Inline CSS — green theme (ATAC), purple theme 유지 (header 그라디언트만 녹색 계열)
css <- '
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            line-height: 1.6;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }
        .header {
            background: linear-gradient(135deg, #2d6a4f 0%, #52b788 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }
        .header h1 { margin: 0 0 10px 0; font-size: 2.2em; }
        .header p  { margin: 4px 0; font-size: 1.05em; opacity: 0.9; }
        .section {
            background: white;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 28px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .section h2 {
            margin-top: 0;
            color: #2d6a4f;
            border-bottom: 2px solid #52b788;
            padding-bottom: 10px;
        }
        .summary-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 16px;
            margin-top: 20px;
        }
        .summary-item {
            background: #f8f9fa;
            padding: 14px;
            border-radius: 8px;
            border-left: 4px solid #52b788;
        }
        .summary-item .label { font-size: 0.88em; color: #6c757d; margin-bottom: 4px; }
        .summary-item .value { font-size: 1.5em; font-weight: bold; color: #212529; }
        table { width: 100%; border-collapse: collapse; margin: 12px 0; font-size: 0.92em; }
        th, td { padding: 10px 12px; text-align: left; border-bottom: 1px solid #dee2e6; }
        th { background-color: #2d6a4f; color: white; font-weight: 600; }
        tr:hover { background-color: #f0faf4; }
        .up   { color: #c0392b; font-weight: 600; }
        .down { color: #2980b9; font-weight: 600; }
        .badge {
            display: inline-block;
            padding: 3px 10px;
            border-radius: 12px;
            font-size: 0.82em;
            font-weight: 600;
        }
        .badge-up   { background: #fdecea; color: #c0392b; }
        .badge-down { background: #eaf4fb; color: #2980b9; }
        .badge-none { background: #f0f0f0; color: #888; }
        details {
            background: white;
            border-radius: 10px;
            margin-bottom: 16px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.08);
            overflow: hidden;
        }
        details > summary {
            padding: 18px 25px;
            cursor: pointer;
            font-size: 1.1em;
            font-weight: 600;
            color: #495057;
            background: #f8f9fa;
            border-radius: 10px;
            list-style: none;
            display: flex;
            align-items: center;
            gap: 12px;
        }
        details > summary::-webkit-details-marker { display: none; }
        details > summary::before {
            content: "\\25B6";
            font-size: 0.75em;
            color: #52b788;
            transition: transform 0.2s;
        }
        details[open] > summary::before { transform: rotate(90deg); }
        details[open] > summary { border-radius: 10px 10px 0 0; }
        .detail-body { padding: 20px 25px; }
        .two-col { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }
        .col-block h4 { margin: 0 0 8px 0; color: #495057; }
        .plot-container { text-align: center; margin: 16px 0; }
        .plot-container img {
            max-width: 80%;
            height: auto;
            border-radius: 6px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }
        .plots-row {
            display: flex;
            gap: 16px;
            justify-content: center;
            margin: 16px 0;
        }
        .plots-row .plot-container { flex: 1; }
        .plots-row .plot-container img { max-width: 100%; }
        .note {
            background: #f0faf4;
            padding: 12px 16px;
            border-left: 4px solid #52b788;
            border-radius: 4px;
            font-size: 0.9em;
            color: #495057;
            margin: 12px 0;
        }
        .note-warn {
            background: #fff8e1;
            padding: 12px 16px;
            border-left: 4px solid #ffc107;
            border-radius: 4px;
            font-size: 0.9em;
            color: #555;
            margin: 12px 0;
        }
        .footer {
            text-align: center;
            margin-top: 40px;
            padding: 20px;
            color: #6c757d;
            font-size: 0.88em;
        }
        @media (max-width: 768px) {
            .two-col { grid-template-columns: 1fr; }
            .summary-grid { grid-template-columns: 1fr 1fr; }
            .plots-row { flex-direction: column; }
        }
'

# Peak table (up/down)
peak_table <- function(df, peak_col, lfc_col, padj_col, direction) {
  if (nrow(df) == 0) {
    return('<p style="color:#888;font-size:0.9em;">No significant DA peaks.</p>')
  }
  colour_class <- if (direction == "up") "up" else "down"
  rows <- apply(df, 1, function(r) {
    lfc_val  <- as.numeric(r[lfc_col])
    padj_val <- as.numeric(r[padj_col])
    peak_val <- r[peak_col]
    sprintf('<tr><td style="font-size:0.88em;">%s</td><td class="%s">%+.3f</td><td>%.2e</td></tr>',
            peak_val, colour_class, lfc_val, padj_val)
  })
  paste0('<table><thead><tr><th>Peak ID</th><th>log2FC</th><th>padj</th></tr></thead><tbody>',
         paste(rows, collapse = "\n"), '</tbody></table>')
}

# GO table
go_table <- function(go_df, top_n = 5) {
  if (nrow(go_df) == 0) {
    return('<p style="color:#888;font-size:0.9em;">No enriched terms.</p>')
  }
  desc_col  <- intersect(c("Description", "description", "Term"), colnames(go_df))[1]
  padj_col2 <- intersect(c("p.adjust", "qvalue", "FDR"), colnames(go_df))[1]
  count_col <- intersect(c("Count", "count"), colnames(go_df))[1]
  if (is.na(desc_col) || is.na(padj_col2)) {
    return('<p style="color:#888;font-size:0.9em;">(column names not recognised)</p>')
  }
  rows <- apply(head(go_df, top_n), 1, function(r) {
    cnt <- if (!is.na(count_col)) r[count_col] else "—"
    sprintf('<tr><td>%s</td><td>%s</td><td>%.3g</td></tr>',
            r[desc_col], cnt, as.numeric(r[padj_col2]))
  })
  paste0('<table><thead><tr><th>GO Term (BP)</th><th>Count</th><th>p.adjust</th></tr></thead><tbody>',
         paste(rows, collapse = "\n"), '</tbody></table>')
}

# Differential TF table (chromVAR)
diff_tf_table <- function(df, top_n = 10, padj_cut = 0.05) {
  if (nrow(df) == 0) {
    return('<p style="color:#888;font-size:0.9em;">No differential TF data.</p>')
  }
  padj_col <- intersect(c("padj", "p_value"), colnames(df))[1]
  df_sorted <- df[order(abs(df$delta), decreasing = TRUE), ]
  df_top    <- head(df_sorted, top_n)
  n_sig     <- if (!is.na(padj_col)) sum(!is.na(df[[padj_col]]) & df[[padj_col]] < padj_cut) else NA_integer_
  note <- if (!is.na(n_sig) && n_sig == 0) {
    sprintf('<p class="note" style="margin-bottom:8px;">No significant TFs (padj &lt; %.2f). Showing top %d by |&Delta; z-score|.</p>', padj_cut, top_n)
  } else if (!is.na(n_sig)) {
    sprintf('<p style="font-size:0.88em;color:#555;margin-bottom:6px;">Significant TFs: <strong>%d</strong> (padj &lt; %.2f)</p>', n_sig, padj_cut)
  } else ""
  rows <- apply(df_top, 1, function(r) {
    delta_val <- as.numeric(r["delta"])
    padj_val  <- if (!is.na(padj_col)) as.numeric(r[padj_col]) else NA_real_
    dir_class <- if (!is.na(delta_val) && delta_val > 0) "up" else "down"
    padj_str  <- if (!is.na(padj_val)) sprintf("%.3g", padj_val) else "&mdash;"
    sprintf('<tr><td style="font-size:0.88em;">%s</td><td class="%s">%+.3f</td><td>%s</td></tr>',
            r["motif"], dir_class, delta_val, padj_str)
  })
  paste0(note,
         '<table><thead><tr><th>TF Motif</th><th>&Delta; z-score</th><th>padj</th></tr></thead><tbody>',
         paste(rows, collapse = "\n"), '</tbody></table>')
}

# --- 5. Build HTML ---

total_pairs <- length(pair_data)
total_da    <- sum(sapply(pair_data, function(x) x$n_total))

html <- paste0('<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>DA Summary Report - ', project_id, '</title>
    <style>', css, '</style>
</head>
<body>

<div class="header">
    <h1>&#127807; ATAC-Seq DA Analysis Summary</h1>
    <p>Project: <strong>', project_id, '</strong></p>
    <p>Species: ', species,
    ' &nbsp;|&nbsp; Method: ', method,
    ' &nbsp;|&nbsp; padj &lt; ', padj_cut,
    ' &nbsp;|&nbsp; |log2FC| &gt; ', lfc_cut, '</p>
    <p>Generated: ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '</p>
</div>

<div class="section">
    <h2>&#128203; Analysis Overview</h2>
    <div class="summary-grid">
        <div class="summary-item">
            <div class="label">Comparisons</div>
            <div class="value">', total_pairs, '</div>
        </div>
        <div class="summary-item">
            <div class="label">Total Samples</div>
            <div class="value">', ifelse(is.na(n_total_samples), "&mdash;", n_total_samples), '</div>
        </div>
        <div class="summary-item">
            <div class="label">Total DA Peaks (all pairs)</div>
            <div class="value">', format(total_da, big.mark = ","), '</div>
        </div>
        <div class="summary-item">
            <div class="label">padj cutoff</div>
            <div class="value">', padj_cut, '</div>
        </div>
        <div class="summary-item">
            <div class="label">|log2FC| cutoff</div>
            <div class="value">', lfc_cut, '</div>
        </div>
    </div>',
  # Sample breakdown table
  if (length(condition_samples) > 0) {
    sample_rows <- paste(sapply(names(condition_samples), function(cond) {
      sids <- condition_samples[[cond]]
      sprintf('<tr><td><strong>%s</strong></td><td>%d</td><td style="font-size:0.88em;color:#555;">%s</td></tr>',
              cond, length(sids), paste(sids, collapse = ", "))
    }), collapse = "\n")
    paste0('
    <h3 style="color:#495057;margin-top:24px;margin-bottom:8px;">Sample Groups</h3>
    <table>
        <thead><tr><th>Condition</th><th>N</th><th>Sample IDs</th></tr></thead>
        <tbody>', sample_rows, '</tbody>
    </table>')
  } else "",
'
</div>
')

# 5b. Results overview table
overview_rows <- sapply(pair_data, function(d) {
  go_up_term <- if (nrow(d$go_up) > 0) {
    dc <- intersect(c("Description","description","Term"), colnames(d$go_up))[1]
    if (!is.na(dc)) d$go_up[[dc]][1] else "&mdash;"
  } else "&mdash;"
  go_down_term <- if (nrow(d$go_down) > 0) {
    dc <- intersect(c("Description","description","Term"), colnames(d$go_down))[1]
    if (!is.na(dc)) d$go_down[[dc]][1] else "&mdash;"
  } else "&mdash;"

  da_badge <- if (d$n_total == 0) {
    '<span class="badge badge-none">0 DA peaks</span>'
  } else {
    sprintf('<strong>%s</strong>', format(d$n_total, big.mark = ","))
  }

  sprintf(
    '<tr>
      <td><strong>%s</strong></td>
      <td>%s</td>
      <td class="up">%s</td>
      <td class="down">%s</td>
      <td style="font-size:0.88em;color:#555;">%s</td>
      <td style="font-size:0.88em;color:#555;">%s</td>
    </tr>',
    d$pair,
    da_badge,
    format(d$n_up,   big.mark = ","),
    format(d$n_down, big.mark = ","),
    go_up_term, go_down_term
  )
})

html <- paste0(html, '
<div class="section">
    <h2>&#128200; Results Overview</h2>
    <table>
        <thead>
            <tr>
                <th>Comparison</th>
                <th>DA Peaks</th>
                <th>Open (up)</th>
                <th>Closed (down)</th>
                <th>Top GO UP (BP)</th>
                <th>Top GO DOWN (BP)</th>
            </tr>
        </thead>
        <tbody>
            ', paste(overview_rows, collapse = "\n"), '
        </tbody>
    </table>
    <div class="note">
        Cutoffs: padj &lt; ', padj_cut, ' &amp; |log2FC| &gt; ', lfc_cut, '.
        GO enrichment on nearest-gene sets; top-ranked Biological Process term shown.
    </div>
    ', if (total_da == 0) {
        '<div class="note-warn">
            &#9888;&#65039; No DA peaks were detected across all comparisons under the current thresholds.
            This may reflect limited chromatin accessibility changes under the tested conditions.
        </div>'
    } else "", '
</div>
')

# 5c. Global PCA plot (if available)
pca_abs <- file.path(output_dir, "global_qc", "pca_all_samples.png")
if (file.exists(pca_abs)) {
  html <- paste0(html, '
<div class="section">
    <h2>&#127760; Global PCA</h2>
    <div class="plot-container">
        <img src="global_qc/pca_all_samples.png" alt="Global PCA">
    </div>
</div>
')
}

# 5d. chromVAR TF Activity (optional)
run_chromvar <- isTRUE(config$chromvar$run_chromvar)
chromvar_dir <- file.path(output_dir, "chromvar")
if (run_chromvar && dir.exists(chromvar_dir)) {
  cat("  Adding chromVAR section...\n")

  var_plot_rel  <- "chromvar/tf_variability_plot.png"
  heat_plot_rel <- "chromvar/tf_deviation_heatmap.png"
  var_plot_tag  <- if (file.exists(file.path(output_dir, var_plot_rel)))
    sprintf('<div class="plot-container"><img src="%s" alt="TF Variability"><p style="font-size:0.85em;color:#888;">TF Variability (all samples)</p></div>', var_plot_rel) else ""
  heat_plot_tag <- if (file.exists(file.path(output_dir, heat_plot_rel)))
    sprintf('<div class="plot-container"><img src="%s" alt="TF Deviation Heatmap"><p style="font-size:0.85em;color:#888;">TF Deviation Z-score Heatmap</p></div>', heat_plot_rel) else ""

  diff_tf_html <- ""
  for (pair in PAIRS) {
    csv_path <- file.path(output_dir, "chromvar", "differential_tf",
                          paste0(pair, "_diff_tf.csv"))
    png_rel  <- file.path("chromvar", "differential_tf",
                          paste0(pair, "_diff_tf_plot.png"))
    if (!file.exists(csv_path)) next
    diff_df  <- tryCatch(read.csv(csv_path, check.names = FALSE),
                         error = function(e) data.frame())
    plot_tag <- if (file.exists(file.path(output_dir, png_rel)))
      sprintf('<div class="plot-container"><img src="%s" alt="Diff TF %s" style="max-width:70%%;"><p style="font-size:0.85em;color:#888;">Top differential TFs</p></div>', png_rel, pair) else ""
    diff_tf_html <- paste0(diff_tf_html, sprintf('
    <details>
        <summary>%s</summary>
        <div class="detail-body">
            %s
            <h4 style="margin-top:14px;">Top Differential TFs (by |&Delta; z-score|)</h4>
            %s
        </div>
    </details>', pair, plot_tag, diff_tf_table(diff_df, top_n = 10, padj_cut = padj_cut)))
  }

  html <- paste0(html, '
<div class="section">
    <h2>&#127981; chromVAR TF Activity</h2>
    <p style="color:#666;font-size:0.93em;">
        TF activity deviation scores computed across all samples using JASPAR2020 motifs.
        Differential TF activity tested per comparison (Wilcoxon rank-sum test).
    </p>
    <div class="plots-row">', var_plot_tag, heat_plot_tag, '</div>',
    if (nchar(diff_tf_html) > 0) paste0('
    <h3 style="color:#495057;margin-top:24px;margin-bottom:10px;">Differential TF Activity per Comparison</h3>
    <p style="color:#666;font-size:0.9em;margin-bottom:12px;">Click a comparison to expand. &Delta; z-score = compare &minus; base.</p>',
    diff_tf_html) else "", '
</div>
')
}

# 5e. Peak Overlap (optional)
run_overlap <- isTRUE(config$peak_overlap$run_overlap)
mc_dir <- file.path(output_dir, "multi_comparison")
if (run_overlap && dir.exists(mc_dir)) {
  cat("  Adding Peak Overlap section...\n")

  overlap_csv <- file.path(mc_dir, "overlap_summary.csv")
  upset_rel   <- "multi_comparison/upset_plot.png"
  upset_tag   <- if (file.exists(file.path(output_dir, upset_rel)))
    sprintf('<div class="plot-container"><img src="%s" alt="UpSet Plot"><p style="font-size:0.85em;color:#888;">UpSet plot — all comparisons</p></div>', upset_rel) else ""

  # Venn diagrams — all venn_*.png (exclude .log files)
  venn_files <- list.files(mc_dir, pattern = "^venn_.*\\.png$", full.names = FALSE)
  venn_tags  <- paste(sapply(venn_files, function(f) {
    rel <- file.path("multi_comparison", f)
    label <- sub("^venn_", "", sub("\\.png$", "", f))
    label <- gsub("_vs_CONTROL|_vs_Veh", "", label)
    label <- gsub("_", " / ", label)
    sprintf('<div class="plot-container"><img src="%s" alt="Venn %s"><p style="font-size:0.85em;color:#888;">%s</p></div>', rel, f, label)
  }), collapse = "\n")

  overlap_table_html <- if (file.exists(overlap_csv)) {
    ov <- tryCatch(read.csv(overlap_csv, check.names = FALSE), error = function(e) data.frame())
    if (nrow(ov) == 0) {
      '<div class="note-warn">&#9888;&#65039; No DA peaks detected in any comparison &mdash; no overlap patterns to display.</div>'
    } else {
      rows <- apply(head(ov, 20), 1, function(r) {
        groups <- strsplit(as.character(r["pattern"]), "\\|")[[1]]
        badges <- paste(sapply(groups, function(g)
          sprintf('<span class="badge badge-up" style="margin-right:4px;font-size:0.78em;">%s</span>', g)),
          collapse = "")
        sprintf('<tr><td>%s</td><td style="text-align:right;font-weight:600;">%s</td></tr>',
                badges, format(as.integer(r["n_peaks"]), big.mark = ","))
      })
      note_more <- if (nrow(ov) > 20) sprintf('<p style="font-size:0.85em;color:#888;">Showing top 20 of %d patterns.</p>', nrow(ov)) else ""
      paste0('<table><thead><tr><th>Overlap Pattern</th><th style="text-align:right;">Peaks</th></tr></thead><tbody>',
             paste(rows, collapse = "\n"), '</tbody></table>', note_more)
    }
  } else '<p style="color:#888;">overlap_summary.csv not found.</p>'

  html <- paste0(html, '
<div class="section">
    <h2>&#128101; Multi-Comparison Peak Overlap</h2>
    <p style="color:#666;font-size:0.93em;">
        Overlap of significant DA peaks (padj &lt; ', padj_cut, ', |log2FC| &gt; ', lfc_cut, ') across all pairwise comparisons.
    </p>
    <h3 style="color:#495057;margin-top:16px;margin-bottom:8px;">Overlap Summary</h3>
    ', overlap_table_html, '
    <div class="plots-row" style="margin-top:20px;">', upset_tag, '</div>',
    if (nchar(venn_tags) > 0) paste0('
    <h3 style="color:#495057;margin-top:20px;margin-bottom:8px;">Venn Diagrams</h3>
    <div class="plots-row">', venn_tags, '</div>') else "", '
</div>
')
}

# 5f. Time-course (optional)
run_timecourse <- isTRUE(config$timecourse$run_timecourse)
tc_dir <- file.path(output_dir, "timecourse")
if (run_timecourse && dir.exists(tc_dir)) {
  cat("  Adding Time-course section...\n")

  cluster_csv <- file.path(tc_dir, "cluster_summary.csv")
  tc_cluster_table <- if (file.exists(cluster_csv)) {
    tc_clust <- tryCatch(read.csv(cluster_csv, check.names = FALSE), error = function(e) data.frame())
    if (nrow(tc_clust) > 0) {
      rows <- apply(tc_clust, 1, function(r)
        sprintf('<tr><td>Cluster %s</td><td style="text-align:right;font-weight:600;">%s</td></tr>',
                r["cluster_id"], format(as.integer(r["n_peaks"]), big.mark = ",")))
      paste0('<table><thead><tr><th>Cluster</th><th style="text-align:right;">Peaks</th></tr></thead><tbody>',
             paste(rows, collapse = "\n"), '</tbody></table>')
    } else ""
  } else ""

  # Diagnostic plots
  elbow_rel <- "timecourse/elbow_plot.png"
  sil_rel   <- "timecourse/silhouette_plot.png"
  heat_rel  <- "timecourse/temporal_heatmap.png"
  diag_tags <- paste(c(
    if (file.exists(file.path(output_dir, elbow_rel)))
      sprintf('<div class="plot-container"><img src="%s" alt="Elbow"><p style="font-size:0.85em;color:#888;">Elbow (WSS)</p></div>', elbow_rel) else NULL,
    if (file.exists(file.path(output_dir, sil_rel)))
      sprintf('<div class="plot-container"><img src="%s" alt="Silhouette"><p style="font-size:0.85em;color:#888;">Silhouette</p></div>', sil_rel) else NULL
  ), collapse = "\n")
  heat_tag <- if (file.exists(file.path(output_dir, heat_rel)))
    sprintf('<div class="plot-container"><img src="%s" alt="Temporal Heatmap"><p style="font-size:0.85em;color:#888;">Temporal DA Peak Heatmap (Z-score)</p></div>', heat_rel) else ""

  # Trend plots
  trend_files <- sort(list.files(tc_dir, pattern = "^trend_cluster_[0-9]+\\.png$", full.names = FALSE))
  trend_tags  <- paste(sapply(trend_files, function(f) {
    rel <- file.path("timecourse", f)
    k   <- sub("trend_cluster_", "", sub("\\.png$", "", f))
    sprintf('<div class="plot-container"><img src="%s" alt="Trend Cluster %s"><p style="font-size:0.85em;color:#888;">Cluster %s</p></div>', rel, k, k)
  }), collapse = "\n")

  # GO per cluster (top 5 terms each)
  go_cluster_html <- ""
  go_files <- sort(list.files(file.path(tc_dir, "go_enrichment"),
                               pattern = "^cluster_[0-9]+_.*\\.csv$",
                               full.names = FALSE))
  if (length(go_files) > 0) {
    go_blocks <- paste(sapply(go_files, function(f) {
      k     <- sub("cluster_", "", strsplit(f, "_")[[1]][2])
      fpath <- file.path(tc_dir, "go_enrichment", f)
      go_df <- tryCatch(read.csv(fpath, check.names = FALSE), error = function(e) data.frame())
      sprintf('<div class="col-block"><h4>Cluster %s GO (BP)</h4>%s</div>', k, go_table(go_df, top_n = 5))
    }), collapse = "\n")
    go_cluster_html <- paste0('<div class="two-col" style="margin-top:16px;">', go_blocks, '</div>')
  }

  time_order <- config$timecourse$time_order %||% NULL
  time_info  <- if (!is.null(time_order))
    sprintf('<div class="note">Timepoint order: %s</div>', paste(time_order, collapse = " &rarr; ")) else ""

  html <- paste0(html, '
<div class="section">
    <h2>&#9200; Time-course Temporal Analysis</h2>
    <p style="color:#666;font-size:0.93em;">
        Union of DA peaks from all comparisons, Z-score normalized per peak across timepoints, then clustered.
    </p>
    ', time_info, '
    <div class="two-col" style="margin-top:16px;">
        <div class="col-block">
            <h3 style="color:#495057;margin-bottom:8px;">Cluster Summary</h3>
            ', tc_cluster_table, '
        </div>
        <div class="col-block">',
    if (nchar(diag_tags) > 0) paste0('<h3 style="color:#495057;margin-bottom:8px;">Optimal k Selection</h3><div class="plots-row">', diag_tags, '</div>') else '', '
        </div>
    </div>
    <h3 style="color:#495057;margin-top:20px;margin-bottom:8px;">Temporal Heatmap</h3>
    ', heat_tag, '
    <h3 style="color:#495057;margin-top:20px;margin-bottom:8px;">Cluster Trend Plots</h3>
    <div class="plots-row">', trend_tags, '</div>',
    go_cluster_html, '
</div>
')
}

# 5g. Per-comparison details
go_onto_str  <- paste(go_ontos, collapse = "/")
kegg_str     <- if (!is.null(kegg_code)) sprintf(" &amp; KEGG (<code>%s</code>)", kegg_code) else ""
species_disp <- switch(tolower(species),
  mouse = "Mouse (<em>Mus musculus</em>)",
  human = "Human (<em>Homo sapiens</em>)",
  species)

html <- paste0(html, sprintf('
<div class="section">
    <h2>&#128300; Per-Comparison Analysis Details</h2>
    <div class="note" style="margin-bottom:18px;">
        <strong>분석 파이프라인 개요</strong>
        <ol style="margin:10px 0 4px 18px;padding:0;line-height:1.9;">
            <li>
                <strong>Differential Accessibility (DA) — %s</strong><br>
                <span style="color:#555;">Raw featureCounts를 입력으로 DESeq2 size-factor 정규화 후 Wald test를 수행하였습니다.
                log2FC는 <code>%s</code> shrinkage로 보정하였으며, 유의 기준은 padj &lt; %s, |log2FC| &gt; %s 입니다.</span>
            </li>
            <li>
                <strong>Peak Annotation — ChIPseeker</strong><br>
                <span style="color:#555;">각 DA peak을 %s 게놈 TxDb에 mapping하여 TSS ±%s bp 이내를 Promoter로 정의하고
                가장 가까운 유전자의 symbol, Ensembl ID, 게놈 영역(Promoter/Intron/Exon/Intergenic 등), TSS까지의 거리를 부여하였습니다.</span>
            </li>
            <li>
                <strong>GO/KEGG Enrichment — clusterProfiler</strong><br>
                <span style="color:#555;">DA peak에 annotation된 유전자를 대상으로 GO (%s)%s enrichment를 수행하였습니다 (BH correction, pAdj &lt; %s).
                Opened(up) / Closed(down) / Total DA peak 유래 유전자 세 그룹에 대해 각각 분석하였습니다.</span>
            </li>
        </ol>
    </div>
    <p style="color:#666;font-size:0.93em;">Click a comparison to expand results.</p>
',
  method, lfc_shrink,
  format(padj_cut, nsmall = 2), format(lfc_cut, nsmall = 1),
  species_disp, format(tss_dist, big.mark = ","),
  go_onto_str, kegg_str, format(padj_cut, nsmall = 2)
))

for (d in pair_data) {
  parts           <- strsplit(d$pair, "_vs_")[[1]]
  compare_label   <- parts[1]
  base_label      <- parts[2]

  compare_samples <- condition_samples[[compare_label]]
  base_samples    <- condition_samples[[base_label]]
  sample_info_html <- if (!is.null(compare_samples) && !is.null(base_samples)) {
    sprintf('
            <div class="note" style="margin-bottom:14px;">
                <strong>%s</strong> (n=%d): %s<br>
                <strong>%s</strong> (n=%d): %s
            </div>',
      compare_label, length(compare_samples), paste(compare_samples, collapse = ", "),
      base_label,    length(base_samples),    paste(base_samples,    collapse = ", "))
  } else ""

  # Check plot availability
  volcano_exists <- file.exists(file.path(output_dir, d$volcano_rel))
  ma_exists      <- file.exists(file.path(output_dir, d$ma_rel))

  plots_html <- if (volcano_exists || ma_exists) {
    v_tag <- if (volcano_exists) sprintf('<div class="plot-container"><img src="%s" alt="Volcano %s"><p style="font-size:0.85em;color:#888;">Volcano plot</p></div>', d$volcano_rel, d$pair) else ""
    m_tag <- if (ma_exists)      sprintf('<div class="plot-container"><img src="%s" alt="MA %s"><p style="font-size:0.85em;color:#888;">MA plot</p></div>',      d$ma_rel,      d$pair) else ""
    paste0('<div class="plots-row">', v_tag, m_tag, '</div>')
  } else ""

  no_da_note <- if (d$n_total == 0) {
    '<div class="note-warn">No DA peaks detected under current thresholds.</div>'
  } else ""

  pct_da <- if (d$n_peaks > 0) sprintf("%.1f%%", d$n_total / d$n_peaks * 100) else "—"
  n_go_up   <- nrow(d$go_up)
  n_go_down <- nrow(d$go_down)
  go_up_note   <- if (n_go_up   > 0) sprintf("%d BP terms", n_go_up)   else "no significant terms"
  go_down_note <- if (n_go_down > 0) sprintf("%d BP terms", n_go_down) else "no significant terms"

  analysis_summary_html <- sprintf('
        <table style="width:100%%;border-collapse:collapse;font-size:0.88em;margin-bottom:16px;">
            <thead>
                <tr style="background:#f0f4f0;">
                    <th style="padding:6px 12px;text-align:left;border-bottom:2px solid #52b788;width:22%%;">항목</th>
                    <th style="padding:6px 12px;text-align:left;border-bottom:2px solid #52b788;width:39%%;">방법</th>
                    <th style="padding:6px 12px;text-align:left;border-bottom:2px solid #52b788;width:39%%;">결과</th>
                </tr>
            </thead>
            <tbody>
                <tr style="background:#fff;">
                    <td style="padding:6px 12px;border-bottom:1px solid #e0e0e0;font-weight:600;">DA 분석</td>
                    <td style="padding:6px 12px;border-bottom:1px solid #e0e0e0;">%s · LFC shrinkage: <code>%s</code><br>padj &lt; %s, |log2FC| &gt; %s</td>
                    <td style="padding:6px 12px;border-bottom:1px solid #e0e0e0;">%s peaks 검사 &rarr; <strong>%s DA</strong> (%s)<br>
                        <span style="color:#c0392b;">&#8593; opened: %s</span> &nbsp;
                        <span style="color:#2980b9;">&#8595; closed: %s</span></td>
                </tr>
                <tr style="background:#fafafa;">
                    <td style="padding:6px 12px;border-bottom:1px solid #e0e0e0;font-weight:600;">Peak Annotation</td>
                    <td style="padding:6px 12px;border-bottom:1px solid #e0e0e0;">ChIPseeker · TSS ±%s bp<br>%s</td>
                    <td style="padding:6px 12px;border-bottom:1px solid #e0e0e0;">Nearest gene symbol, Ensembl ID,<br>genomic region, distance to TSS</td>
                </tr>
                <tr style="background:#fff;">
                    <td style="padding:6px 12px;font-weight:600;">GO/KEGG</td>
                    <td style="padding:6px 12px;">clusterProfiler · %s%s<br>BH correction, pAdj &lt; %s</td>
                    <td style="padding:6px 12px;">Opened peaks: %s<br>Closed peaks: %s</td>
                </tr>
            </tbody>
        </table>',
    method, lfc_shrink,
    format(padj_cut, nsmall = 2), format(lfc_cut, nsmall = 1),
    format(d$n_peaks, big.mark = ","),
    format(d$n_total, big.mark = ","), pct_da,
    format(d$n_up,   big.mark = ","),
    format(d$n_down, big.mark = ","),
    format(tss_dist, big.mark = ","), species_disp,
    go_onto_str, kegg_str, format(padj_cut, nsmall = 2),
    go_up_note, go_down_note
  )

  html <- paste0(html, sprintf('
    <details>
        <summary>
            %s
            &nbsp;<span class="badge badge-up">&#8593; %s open</span>
            <span class="badge badge-down">&#8595; %s closed</span>
        </summary>
        <div class="detail-body">

            %s
            %s
            %s
            %s

            <div class="two-col">
                <div class="col-block">
                    <h4>&#128308; Top Opened Peaks (up)</h4>
                    %s
                </div>
                <div class="col-block">
                    <h4>&#128309; Top Closed Peaks (down)</h4>
                    %s
                </div>
            </div>

            <div class="two-col" style="margin-top:16px;">
                <div class="col-block">
                    <h4>GO Biological Process &mdash; Opened peaks</h4>
                    %s
                </div>
                <div class="col-block">
                    <h4>GO Biological Process &mdash; Closed peaks</h4>
                    %s
                </div>
            </div>

        </div>
    </details>
',
    d$pair,
    format(d$n_up,   big.mark = ","),
    format(d$n_down, big.mark = ","),
    sample_info_html,
    analysis_summary_html,
    no_da_note,
    plots_html,
    peak_table(d$top_up,   d$peak_col, d$lfc_col, d$padj_col, "up"),
    peak_table(d$top_down, d$peak_col, d$lfc_col, d$padj_col, "down"),
    go_table(d$go_up),
    go_table(d$go_down)
  ))
}

html <- paste0(html, '\n</div>\n')  # close .section

# 5e. Footer
html <- paste0(html, '
<div class="footer">
    <p>Generated by ATAC-Seq DA Analysis Pipeline</p>
    <p>Configuration: ', opt$config, '</p>
</div>

</body>
</html>
')

# --- 6. Write output ---
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
writeLines(html, output_file)

cat("=== Summary Report Generation Complete ===\n")
cat(paste0("HTML report saved to: ", output_file, "\n"))
cat(paste0("View in browser: file://", normalizePath(output_file), "\n\n"))
