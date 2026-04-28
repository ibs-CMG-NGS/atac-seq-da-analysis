#!/usr/bin/env Rscript
# 파일 경로: src/analysis/08_generate_methods_section.R
# 목적: DA 분석에 사용된 모든 도구, 파라미터, 버전 정보를 담은 Methods 섹션 Markdown 생성
# RNA-Seq_DE_GO_analysis/08_generate_methods_section.R 과 구조/형식 통일
#
# Usage:
#   Rscript src/analysis/08_generate_methods_section.R \
#     --config  configs/config_PROJECT.yaml \
#     --output-dir output/PROJECT \
#     --output  output/PROJECT/methods_section.md

suppressPackageStartupMessages({
  library(yaml)
  library(optparse)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !identical(a, "")) a else b

# ─────────────────────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────────────────────
option_list <- list(
  make_option(c("-c", "--config"),     type = "character",
              help = "Path to config YAML file", metavar = "FILE"),
  make_option(c("-d", "--output-dir"), type = "character",
              help = "Pipeline output directory", metavar = "DIR"),
  make_option(c("-o", "--output"),     type = "character",
              help = "Output Markdown file path", metavar = "FILE")
)
opt_parser <- OptionParser(option_list = option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$config) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("--config and --output are required.", call. = FALSE)
}

cat("\n=== ATAC-Seq DA Methods Section Generator ===\n")
cat(paste0("Config: ", opt$config, "\n"))
cat(paste0("Output: ", opt$output, "\n\n"))

# ─────────────────────────────────────────────────────────────
# Load config
# ─────────────────────────────────────────────────────────────
cfg <- yaml.load_file(opt$config)

project_id   <- basename(cfg$output_dir %||% opt$`output-dir` %||% "PROJECT")
project_name <- cfg$project_name %||% project_id
species      <- cfg$species      %||% "unknown"

da_cfg   <- cfg$da_analysis %||% list()
enr_cfg  <- cfg$enrichment  %||% list()
mot_cfg  <- cfg$motif       %||% list()
chr_cfg  <- cfg$chromvar    %||% list()
adv      <- da_cfg$advanced_options %||% list()

da_method    <- da_cfg$method           %||% "DESeq2"
design_f     <- da_cfg$design_formula   %||% "~ condition"
group_var    <- da_cfg$group_variable   %||% "condition"
padj_cut     <- da_cfg$padj_cutoff      %||% 0.05
lfc_cut      <- da_cfg$log2fc_cutoff    %||% 1.0
min_count    <- da_cfg$min_peak_count   %||% 10
min_sample_n <- da_cfg$min_sample_n     %||% 2

run_motif   <- mot_cfg$run_motif    %||% FALSE
run_chromvar<- chr_cfg$run_chromvar %||% FALSE

pairs_raw <- da_cfg$pairwise_comparisons %||% list()
pairs     <- vapply(pairs_raw, function(p) paste0(p[[1]], "_vs_", p[[2]]), character(1))

go_ontologies <- paste(enr_cfg$go_ontologies %||% c("BP","CC","MF"), collapse = ", ")
gene_lists    <- paste(enr_cfg$gene_lists    %||% c("total","up","down"), collapse = ", ")
go_pval       <- enr_cfg$pvalue_cutoff %||% 0.05
go_qval       <- enr_cfg$qvalue_cutoff %||% 0.25
min_gs        <- enr_cfg$min_gs_size   %||% 5
max_gs        <- enr_cfg$max_gs_size   %||% 500
plot_top_n    <- enr_cfg$plot_top_n    %||% 15

species_lower <- tolower(species)
db_cfg        <- (cfg$databases %||% list())[[species_lower]] %||% list()
organism_db   <- db_cfg$organism_db %||% if (species_lower == "mouse") "org.Mm.eg.db" else "org.Hs.eg.db"
kegg_code     <- db_cfg$kegg_code   %||% if (species_lower == "mouse") "mmu" else "hsa"
txdb          <- db_cfg$txdb        %||% if (species_lower == "mouse") "TxDb.Mmusculus.UCSC.mm10.knownGene" else "TxDb.Hsapiens.UCSC.hg38.knownGene"
bsgenome      <- db_cfg$bsgenome    %||% if (species_lower == "mouse") "BSgenome.Mmusculus.UCSC.mm10" else "BSgenome.Hsapiens.UCSC.hg38"

# ─────────────────────────────────────────────────────────────
# Sample / condition info from metadata
# ─────────────────────────────────────────────────────────────
n_samples  <- NA_integer_
n_peaks    <- NA_integer_
conditions <- list()

if (!is.null(cfg$count_matrix_path) && file.exists(cfg$count_matrix_path)) {
  tryCatch({
    raw <- read.table(cfg$count_matrix_path, header = TRUE, sep = "\t",
                      skip = 1, nrows = 1, check.names = FALSE)
    n_peaks   <- as.integer(system(paste("wc -l <", cfg$count_matrix_path), intern = TRUE)) - 2L
    n_samples <- ncol(raw) - 6L   # featureCounts 고정 6 메타 컬럼
    cat(sprintf("[INFO] Peaks: ~%d, Samples: %d\n", n_peaks, n_samples))
  }, error = function(e) cat(paste0("[WARN] count_matrix: ", e$message, "\n")))
}

if (!is.null(cfg$metadata_path) && file.exists(cfg$metadata_path)) {
  tryCatch({
    meta <- read.csv(cfg$metadata_path, row.names = 1)
    if (group_var %in% colnames(meta))
      conditions <- split(rownames(meta), meta[[group_var]])
    if (is.na(n_samples)) n_samples <- nrow(meta)
    cat(sprintf("[INFO] Conditions: %s\n", paste(names(conditions), collapse = ", ")))
  }, error = function(e) cat(paste0("[WARN] metadata: ", e$message, "\n")))
}

# ─────────────────────────────────────────────────────────────
# Package versions
# ─────────────────────────────────────────────────────────────
pkg_ver <- function(pkg) tryCatch(as.character(packageVersion(pkg)), error = function(e) "N/A")

r_ver          <- paste0(R.version$major, ".", R.version$minor)
deseq2_ver     <- pkg_ver("DESeq2")
clusterp_ver   <- pkg_ver("clusterProfiler")
annodbi_ver    <- pkg_ver("AnnotationDbi")
orgdb_ver      <- pkg_ver(organism_db)
txdb_ver       <- pkg_ver(sub(".*\\.", "", txdb))   # 마지막 패키지명 추출
chipseekerv    <- pkg_ver("ChIPseeker")
chromvar_ver   <- pkg_ver("chromVAR")
motifmatchr_v  <- pkg_ver("motifmatchr")
enrichplot_ver <- pkg_ver("enrichplot")
ggplot2_ver    <- pkg_ver("ggplot2")
snakemake_ver  <- tryCatch(
  trimws(system2("snakemake", "--version", stdout = TRUE, stderr = FALSE)[1]),
  error = function(e) "N/A"
)

# ─────────────────────────────────────────────────────────────
# Markdown builder (RNA-seq 버전과 동일 헬퍼)
# ─────────────────────────────────────────────────────────────
lines <- character(0)

h     <- function(level, text) lines <<- c(lines, paste0(strrep("#", level), " ", text), "")
p     <- function(text = "")   lines <<- c(lines, text)
blank <- function()             lines <<- c(lines, "")
li    <- function(text)         lines <<- c(lines, paste0("- ", text))

md_table <- function(headers, rows) {
  ncol   <- length(headers)
  widths <- nchar(headers)
  for (row in rows)
    for (i in seq_along(row))
      if (i <= ncol) widths[i] <- max(widths[i], nchar(as.character(row[[i]])), 3)

  fmt_row <- function(cells) {
    padded <- vapply(seq_len(ncol), function(i)
      formatC(as.character(cells[[i]]), width = -widths[i], flag = "-"), character(1))
    paste0("| ", paste(padded, collapse = " | "), " |")
  }
  sep <- paste0("| ", paste(strrep("-", widths), collapse = " | "), " |")
  lines <<- c(lines, fmt_row(headers), sep)
  for (row in rows) lines <<- c(lines, fmt_row(row))
  lines <<- c(lines, "")
}

# ─────────────────────────────────────────────────────────────
# Build document
# ─────────────────────────────────────────────────────────────
today       <- format(Sys.Date(), "%Y-%m-%d")
species_cap <- paste0(toupper(substring(species, 1, 1)), substring(species, 2))

if (length(conditions) > 0) {
  cond_str <- paste(
    vapply(names(conditions), function(g)
      paste0(g, " (n=", length(conditions[[g]]), ")"), character(1)),
    collapse = ", "
  )
} else {
  cond_str <- if (!is.na(n_samples)) paste0("n=", n_samples) else "unknown"
}

# ── Title ──────────────────────────────────────────────────────
h(1, "ATAC-Seq Differential Accessibility & Enrichment Analysis Methods")
p(paste0("**Project:** ", project_name, " (`", project_id, "`)  "))
p(paste0("**Generated:** ", today, "  "))
p(paste0("**Species:** *", species_cap, "*  "))
if (!is.na(n_samples))
  p(paste0("**Samples:** ", n_samples, " (", cond_str, ")  "))
if (!is.na(n_peaks))
  p(paste0("**Consensus peaks (pre-filter):** ", format(n_peaks, big.mark = ","), "  "))
blank()
p("> **Note to researcher:** This document describes all tools, versions, and parameter")
p("> settings used in the analysis. Please review each section and retain only the")
p("> information relevant to your manuscript. Suggested citation keys are in the References section.")
blank()

# ── 1. Input Data ──────────────────────────────────────────────
h(2, "1. Input Data")
p(paste0(
  "ATAC-seq reads were processed using the **nf-core/atacseq** pipeline. ",
  "Consensus peak regions across all samples were identified, and read counts per peak ",
  "were quantified using **featureCounts** as input for downstream differential accessibility analysis."
))
blank()
if (length(conditions) > 0) {
  cond_rows <- lapply(names(conditions), function(g)
    list(g, length(conditions[[g]]), paste(conditions[[g]], collapse = ", ")))
  md_table(c("Condition", "N samples", "Sample IDs"), cond_rows)
}
if (length(pairs) > 0) {
  p("**Pairwise comparisons:**")
  blank()
  pair_rows <- lapply(pairs, function(p_str) {
    parts <- strsplit(p_str, "_vs_")[[1]]
    list(parts[1], parts[2], p_str)
  })
  md_table(c("Treatment", "Control", "Comparison ID"), pair_rows)
}

# ── 2. Differential Accessibility Analysis ─────────────────────
h(2, "2. Differential Accessibility Analysis")
p(paste0(
  "Differential chromatin accessibility analysis was performed using **DESeq2** (v", deseq2_ver,
  ") [Love 2014], applying the same negative binomial framework used for RNA-seq differential expression. ",
  "Raw featureCounts read counts per consensus peak were used as input. ",
  "Library size normalization used the median-of-ratios method. ",
  "Peak-wise dispersion estimates were empirically Bayes-shrunk toward a fitted trend. ",
  "Log2 fold changes (accessibility in treatment vs. control) were shrunk using the **ashr** estimator ",
  "to reduce noise for low-count peaks. ",
  "Statistical significance was assessed using the Wald test."
))
blank()
p(paste0(
  "Low-count peaks were removed prior to fitting: peaks with fewer than **",
  min_count, "** reads in at least **", min_sample_n, "** samples were excluded."
))
blank()
p("**Analysis parameters:**")
blank()
md_table(
  c("Parameter", "Value", "Description"),
  list(
    list("DA method",                  da_method,    "Statistical framework"),
    list("Design formula",             paste0("`", design_f, "`"), "Model formula"),
    list("Group variable",             group_var,    "Metadata column for grouping"),
    list("Pre-filter (min counts)",    min_count,    "Minimum reads per peak"),
    list("Pre-filter (min samples)",   min_sample_n, "Minimum samples passing count threshold"),
    list("LFC shrinkage",              "ashr",       "Log2FC shrinkage estimator")
  )
)
p("**Significance cutoffs:**")
blank()
md_table(
  c("Cutoff", "Value", "Description"),
  list(
    list("Adjusted p-value (padj)", padj_cut, "Benjamini-Hochberg FDR-corrected Wald test p-value"),
    list("log2 fold change",        lfc_cut,  "Minimum absolute log2(treatment / control) accessibility")
  )
)

# ── 3. Peak Annotation ─────────────────────────────────────────
h(2, "3. Peak Annotation")
p(paste0(
  "Consensus peaks were annotated to genomic features using **ChIPseeker** (v", chipseekerv,
  ") [Yu 2015] with the **", txdb, "** transcript database. ",
  "Each peak was assigned to the nearest gene and classified by genomic context ",
  "(promoter, 5'UTR, exon, intron, 3'UTR, distal intergenic). ",
  "Promoter regions were defined as ±3 kb from the transcription start site (TSS). ",
  "Gene symbols were retrieved via **", organism_db, "** (v", orgdb_ver, ") through **AnnotationDbi** (v", annodbi_ver, ")."
))
blank()

# ── 4. Functional Enrichment Analysis ──────────────────────────
h(2, "4. Functional Enrichment Analysis")
p(paste0(
  "Over-representation analysis (ORA) of peak-associated genes was performed using ",
  "**clusterProfiler** (v", clusterp_ver, ") [Wu 2021]. ",
  "GO enrichment was tested across ontologies: Biological Process (BP), ",
  "Cellular Component (CC), and Molecular Function (MF). ",
  "KEGG pathway enrichment was additionally performed (KEGG organism code: `", kegg_code, "`). ",
  "Analysis was conducted separately for all DA peaks, peaks with increased accessibility ",
  "(opened), and peaks with decreased accessibility (closed)."
))
blank()
md_table(
  c("Parameter", "Value", "Description"),
  list(
    list("GO ontologies",       go_ontologies, "GO sub-ontologies tested"),
    list("Peak set directions", gene_lists,    "All DA, opened, and closed peaks"),
    list("Organism DB",         organism_db,   "Bioconductor annotation database"),
    list("KEGG organism code",  kegg_code,     "KEGG species identifier"),
    list("p-value cutoff",      go_pval,       "Maximum raw p-value for terms to report"),
    list("q-value cutoff",      go_qval,       "Maximum FDR-adjusted p-value (q-value)"),
    list("Min gene set size",   min_gs,        "Minimum annotated genes per GO/KEGG term"),
    list("Max gene set size",   max_gs,        "Maximum annotated genes per GO/KEGG term"),
    list("Plot top N",          plot_top_n,    "Top terms shown in dot/bar plots")
  )
)

# ── 5. Motif Analysis (optional) ───────────────────────────────
if (isTRUE(run_motif)) {
  h(2, "5. Motif Enrichment Analysis")
  homer_opts <- mot_cfg$homer_options %||% "-size 200 -mask"
  p(paste0(
    "De novo and known transcription factor motif enrichment was performed using **HOMER** ",
    "[Heinz 2010]. DA peak sequences (", homer_opts, ") were analyzed separately for ",
    "opened and closed peak sets relative to the full consensus peak set as background. ",
    "Both de novo motif discovery and known motif enrichment were performed."
  ))
  blank()
}

# ── 6. TF Activity (chromVAR, optional) ────────────────────────
if (isTRUE(run_chromvar)) {
  h(2, paste0(if (isTRUE(run_motif)) "6" else "5", ". TF Activity Inference (chromVAR)"))
  p(paste0(
    "Transcription factor (TF) activity was inferred across samples using **chromVAR** (v", chromvar_ver,
    ") [Schep 2017] with **motifmatchr** (v", motifmatchr_v, ") for motif-to-peak mapping. ",
    "TF motif positions were obtained from the **JASPAR 2022** database (vertebrates, non-redundant set). ",
    "The **", bsgenome, "** reference genome sequence was used for GC-bias correction. ",
    "Per-sample TF deviation scores (z-scores) were computed and used to identify ",
    "differentially active TFs between conditions using the Wilcoxon rank-sum test."
  ))
  blank()
}

# ── Software Versions ──────────────────────────────────────────
sec_num <- sum(c(TRUE, TRUE, TRUE, TRUE, isTRUE(run_motif), isTRUE(run_chromvar))) + 1
h(2, paste0(sec_num, ". Software Versions"))

sw_rows <- list(
  list("R",              r_ver,        "Statistical computing environment"),
  list("DESeq2",         deseq2_ver,   "Differential accessibility testing"),
  list("ChIPseeker",     chipseekerv,  "Peak annotation"),
  list("clusterProfiler",clusterp_ver, "GO and KEGG enrichment analysis"),
  list("AnnotationDbi",  annodbi_ver,  "Bioconductor annotation interface"),
  list(organism_db,      orgdb_ver,    paste0("Gene annotation (", species_cap, ")")),
  list(txdb,             "auto",       paste0("Transcript database (", species_cap, ")")),
  list("enrichplot",     enrichplot_ver,"Enrichment result visualization"),
  list("ggplot2",        ggplot2_ver,  "Data visualization"),
  list("Snakemake",      snakemake_ver,"Workflow management")
)
if (isTRUE(run_chromvar)) {
  sw_rows <- c(sw_rows, list(
    list("chromVAR",     chromvar_ver,   "TF activity inference"),
    list("motifmatchr",  motifmatchr_v,  "Motif-to-peak mapping"),
    list(bsgenome,       "auto",         paste0("Reference genome sequence (", species_cap, ")"))
  ))
}

md_table(c("Tool / Package", "Version", "Purpose"), sw_rows)
p(paste0(
  "All R packages were managed via conda (environment: `atac-seq-da-analysis`). ",
  "The full environment specification is available at `environment.yml` in the pipeline root directory."
))
blank()

# ── References ─────────────────────────────────────────────────
h(2, paste0(sec_num + 1, ". References"))
li(paste0("Love MI, et al. (2014). Moderated estimation of fold change and dispersion for ",
          "RNA-seq data with DESeq2. *Genome Biology*, 15:550."))
li(paste0("Wu T, et al. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. ",
          "*Innovation*, 2(3):100141."))
li(paste0("Yu G, et al. (2015). ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization. ",
          "*Bioinformatics*, 31(14):2382–2383."))
li(paste0("Zhu A, et al. (2019). Heavy-tailed prior distributions for sequence count data. ",
          "*Bioinformatics*, 35(12):2084–2092."))
if (isTRUE(run_motif))
  li(paste0("Heinz S, et al. (2010). Simple combinations of lineage-determining transcription factors prime cis-regulatory elements. ",
            "*Molecular Cell*, 38(4):576–589."))
if (isTRUE(run_chromvar))
  li(paste0("Schep AN, et al. (2017). chromVAR: inferring transcription-factor-associated accessibility from single-cell epigenomic data. ",
            "*Nature Methods*, 14:975–978."))
li(paste0("Pagès H, et al. (2023). AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor. ",
          "R package version ", annodbi_ver, "."))
li(paste0("Köster J & Rahmann S. (2012). Snakemake — a scalable bioinformatics workflow engine. ",
          "*Bioinformatics*, 28(19):2520–2522."))
blank()

# ── Appendix: Full Configuration ───────────────────────────────
h(2, "Appendix: Full Configuration Parameters")
p(paste0(
  "The following parameters were specified in the project configuration file ",
  "(`", basename(opt$config), "`). This appendix provides a complete record for reproducibility."
))
blank()

flatten_yaml <- function(x, prefix = "") {
  rows <- list()
  for (k in names(x)) {
    key <- if (nchar(prefix) > 0) paste0(prefix, ".", k) else k
    v   <- x[[k]]
    if (is.list(v) && !is.null(names(v))) {
      rows <- c(rows, flatten_yaml(v, key))
    } else {
      display <- paste(as.character(v), collapse = ", ")
      if (nchar(display) > 90) display <- paste0(substr(display, 1, 87), "...")
      rows <- c(rows, list(list(paste0("`", key, "`"), paste0("`", display, "`"))))
    }
  }
  rows
}

md_table(c("Parameter", "Value"), flatten_yaml(cfg))

# ─────────────────────────────────────────────────────────────
# Write output
# ─────────────────────────────────────────────────────────────
output_path <- opt$output
if (!dir.exists(dirname(output_path))) dir.create(dirname(output_path), recursive = TRUE)
writeLines(lines, output_path)
cat(sprintf("[INFO] Methods section written to: %s\n", output_path))
