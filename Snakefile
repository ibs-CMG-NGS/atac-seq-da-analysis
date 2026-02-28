"""
ATAC-Seq Differential Accessibility (DA) Analysis Pipeline
Snakemake workflow - mirrors RNA-Seq_DE_GO_analysis structure

Input:  featureCounts.txt (from nf-core/atacseq consensus peaks)
Output: DA results, peak annotation, GO/KEGG enrichment, motif analysis

Usage:
    snakemake --configfile configs/config_MyProject.yml --cores 8
    snakemake --configfile configs/config_MyProject.yml --cores 8 --dryrun
"""

import os
from pathlib import Path

configfile: "configs/template/config.yml"

# ──────────────────────────────────────────────
# Derived parameters
# ──────────────────────────────────────────────
PAIRS       = [f"{c[0]}_vs_{c[1]}" for c in config["da_analysis"]["pairwise_comparisons"]]
GENE_LISTS  = config["enrichment"]["gene_lists"]          # ["total", "up", "down"]
GO_ONTOS    = config["enrichment"]["go_ontologies"]       # ["BP", "CC", "MF"]
OUT         = config["output_dir"]
RUN_MOTIF   = config.get("motif", {}).get("run_motif", False)
RUN_CHROMVAR= config.get("chromvar", {}).get("run_chromvar", False)

# ──────────────────────────────────────────────
# Target rule
# ──────────────────────────────────────────────
rule all:
    input:
        # Global QC
        expand("{out}/global_qc/.global_qc_done.flag", out=OUT),
        # Per-comparison outputs
        expand("{out}/pairwise/{pair}/final_da_results.csv",   out=OUT, pair=PAIRS),
        expand("{out}/pairwise/{pair}/volcano_plot.png",       out=OUT, pair=PAIRS),
        expand("{out}/pairwise/{pair}/.enrichment_done.flag",  out=OUT, pair=PAIRS),
        expand("{out}/pairwise/{pair}/.go_barplots_done.flag", out=OUT, pair=PAIRS),
        expand("{out}/pairwise/{pair}/final_da_summary.xlsx",  out=OUT, pair=PAIRS),
        # Motif (optional)
        expand("{out}/pairwise/{pair}/.motif_done.flag", out=OUT, pair=PAIRS)
            if RUN_MOTIF else [],


# ──────────────────────────────────────────────
# Rule 0: Global QC (all samples)
# ──────────────────────────────────────────────
rule global_qc_plots:
    input:
        counts   = config["count_matrix_path"],
        metadata = config["metadata_path"],
    output:
        flag = touch("{out}/global_qc/.global_qc_done.flag"),
    params:
        out_dir = "{out}/global_qc",
        config  = lambda wildcards: workflow.configfiles[0],
    log:
        "{out}/logs/global_qc.log"
    shell:
        "Rscript src/analysis/02a_global_qc.R "
        "{params.config} {params.out_dir} > {log} 2>&1"


# ──────────────────────────────────────────────
# Rule 1: Differential Accessibility (per comparison)
# ──────────────────────────────────────────────
rule run_pairwise_da:
    input:
        counts   = config["count_matrix_path"],
        metadata = config["metadata_path"],
    output:
        results = "{out}/pairwise/{pair}/final_da_results.csv",
        counts  = "{out}/pairwise/{pair}/normalized_counts.csv",
    params:
        compare = lambda wc: wc.pair.split("_vs_")[0],
        base    = lambda wc: wc.pair.split("_vs_")[1],
        out_dir = "{out}/pairwise/{pair}",
        config  = lambda wildcards: workflow.configfiles[0],
    log:
        "{out}/logs/{pair}_da.log"
    shell:
        "Rscript src/analysis/01_run_pairwise_da.R "
        "{params.config} {params.compare} {params.base} {params.out_dir} "
        "> {log} 2>&1"


# ──────────────────────────────────────────────
# Rule 2a: Pairwise QC plots (volcano, MA)
# ──────────────────────────────────────────────
rule pairwise_plots:
    input:
        results = "{out}/pairwise/{pair}/final_da_results.csv",
    output:
        volcano = "{out}/pairwise/{pair}/volcano_plot.png",
        flag    = touch("{out}/pairwise/{pair}/.pairwise_qc_done.flag"),
    params:
        compare = lambda wc: wc.pair.split("_vs_")[0],
        base    = lambda wc: wc.pair.split("_vs_")[1],
        out_dir = "{out}/pairwise/{pair}",
        config  = lambda wildcards: workflow.configfiles[0],
    log:
        "{out}/logs/{pair}_plots.log"
    shell:
        "Rscript src/analysis/02c_pairwise_plots.R "
        "{params.config} {params.compare} {params.base} {params.out_dir} "
        "> {log} 2>&1"


# ──────────────────────────────────────────────
# Rule 3: Peak annotation + GO/KEGG enrichment
# ──────────────────────────────────────────────
rule peak_annotation_go:
    input:
        results = "{out}/pairwise/{pair}/final_da_results.csv",
        peaks   = config["consensus_peaks_path"],
    output:
        go  = expand(
            "{{out}}/pairwise/{{pair}}/go_enrichment_{geneset}_{onto}.csv",
            geneset=GENE_LISTS, onto=GO_ONTOS
        ),
        kegg = expand(
            "{{out}}/pairwise/{{pair}}/kegg_enrichment_{geneset}.csv",
            geneset=GENE_LISTS
        ),
        flag = touch("{out}/pairwise/{pair}/.enrichment_done.flag"),
    params:
        compare = lambda wc: wc.pair.split("_vs_")[0],
        base    = lambda wc: wc.pair.split("_vs_")[1],
        out_dir = "{out}/pairwise/{pair}",
        config  = lambda wildcards: workflow.configfiles[0],
    log:
        "{out}/logs/{pair}_enrichment.log"
    shell:
        "Rscript src/analysis/03_peak_annotation_go.R "
        "{params.config} {params.compare} {params.base} {params.out_dir} "
        "> {log} 2>&1"


# ──────────────────────────────────────────────
# Rule 4: GO barplots
# ──────────────────────────────────────────────
rule go_barplots:
    input:
        flag = "{out}/pairwise/{pair}/.enrichment_done.flag",
    output:
        flag = touch("{out}/pairwise/{pair}/.go_barplots_done.flag"),
    params:
        out_dir = "{out}/pairwise/{pair}",
        config  = lambda wildcards: workflow.configfiles[0],
    log:
        "{out}/logs/{pair}_barplots.log"
    shell:
        "Rscript src/analysis/04_go_barplots.R "
        "{params.config} {params.out_dir} > {log} 2>&1"


# ──────────────────────────────────────────────
# Rule 5: Motif enrichment (ATAC-specific, optional)
# ──────────────────────────────────────────────
rule motif_enrichment:
    input:
        results = "{out}/pairwise/{pair}/final_da_results.csv",
        peaks   = config["consensus_peaks_path"],
    output:
        flag = touch("{out}/pairwise/{pair}/.motif_done.flag"),
    params:
        compare = lambda wc: wc.pair.split("_vs_")[0],
        base    = lambda wc: wc.pair.split("_vs_")[1],
        out_dir = "{out}/pairwise/{pair}/motif",
        config  = lambda wildcards: workflow.configfiles[0],
    log:
        "{out}/logs/{pair}_motif.log"
    shell:
        "Rscript src/analysis/05_motif_analysis.R "
        "{params.config} {params.compare} {params.base} {params.out_dir} "
        "> {log} 2>&1"


# ──────────────────────────────────────────────
# Rule 6: Summary Excel table
# ──────────────────────────────────────────────
rule generate_da_summary_table:
    input:
        results = "{out}/pairwise/{pair}/final_da_results.csv",
        flag    = "{out}/pairwise/{pair}/.enrichment_done.flag",
    output:
        xlsx = "{out}/pairwise/{pair}/final_da_summary.xlsx",
    params:
        out_dir = "{out}/pairwise/{pair}",
        config  = lambda wildcards: workflow.configfiles[0],
    log:
        "{out}/logs/{pair}_summary_table.log"
    shell:
        "Rscript src/analysis/06_generate_da_table.R "
        "{params.config} {params.out_dir} > {log} 2>&1"
