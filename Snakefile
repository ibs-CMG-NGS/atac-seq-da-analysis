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
RUN_MOTIF      = config.get("motif",        {}).get("run_motif",      False)
RUN_CHROMVAR   = config.get("chromvar",     {}).get("run_chromvar",   False)
RUN_SEQVIEWER  = config.get("seqviewer",    {}).get("run_export",     True)
RUN_OVERLAP    = config.get("peak_overlap", {}).get("run_overlap",    False)
RUN_TIMECOURSE = config.get("timecourse",   {}).get("run_timecourse", False)

# ──────────────────────────────────────────────
# Target rule
# ──────────────────────────────────────────────
rule all:
    input:
        # Global QC
        expand("{out}/global_qc/.global_qc_done.flag", out=OUT),
        expand("{out}/global_qc/global_qc_report.html", out=OUT),
        # Per-comparison outputs
        expand("{out}/pairwise/{pair}/final_da_results.csv",   out=OUT, pair=PAIRS),
        expand("{out}/pairwise/{pair}/volcano_plot.png",       out=OUT, pair=PAIRS),
        expand("{out}/pairwise/{pair}/.enrichment_done.flag",  out=OUT, pair=PAIRS),
        expand("{out}/pairwise/{pair}/.go_barplots_done.flag", out=OUT, pair=PAIRS),
        expand("{out}/pairwise/{pair}/final_da_result.xlsx",   out=OUT, pair=PAIRS),
        expand("{out}/pairwise/{pair}/final_go_result.xlsx",  out=OUT, pair=PAIRS),
        # Motif (optional)
        expand("{out}/pairwise/{pair}/.motif_done.flag", out=OUT, pair=PAIRS)
            if RUN_MOTIF else [],
        # chromVAR (optional, global)
        expand("{out}/chromvar/.chromvar_done.flag", out=OUT)
            if RUN_CHROMVAR else [],
        # Peak overlap (optional, multi-comparison)
        expand("{out}/multi_comparison/.overlap_done.flag", out=OUT)
            if RUN_OVERLAP else [],
        # Time-course (optional)
        expand("{out}/timecourse/.timecourse_done.flag", out=OUT)
            if RUN_TIMECOURSE else [],
        # SeqViewer export (optional)
        expand("{out}/pairwise/{pair}/.seqviewer_export_done.flag", out=OUT, pair=PAIRS)
            if RUN_SEQVIEWER else [],
        # SeqViewer aggregate manifest (optional)
        expand("{out}/seqviewer/seqviewer_manifest.json", out=OUT)
            if RUN_SEQVIEWER else [],
        # Methods section
        expand("{out}/methods_section.md", out=OUT),
        # Summary report
        expand("{out}/summary_report.html", out=OUT),


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
# Rule 7: chromVAR TF activity (optional, global — runs once for all samples)
# ──────────────────────────────────────────────
rule chromvar_analysis:
    input:
        counts   = config["count_matrix_path"],
        metadata = config["metadata_path"],
        peaks    = config["consensus_peaks_path"],
    output:
        flag = touch("{out}/chromvar/.chromvar_done.flag"),
    params:
        out_dir = "{out}/chromvar",
        config  = lambda wildcards: workflow.configfiles[0],
    log:
        "{out}/logs/chromvar.log"
    shell:
        "Rscript src/analysis/07_chromvar_analysis.R "
        "{params.config} {params.out_dir} "
        "> {log} 2>&1"


# ──────────────────────────────────────────────
# Rule 13: Multi-comparison peak overlap (optional)
# ──────────────────────────────────────────────
rule multi_comparison_overlap:
    input:
        expand("{{out}}/pairwise/{pair}/final_da_results.csv", pair=PAIRS),
    output:
        flag = touch("{out}/multi_comparison/.overlap_done.flag"),
    params:
        out_dir = "{out}/multi_comparison",
        config  = lambda wildcards: workflow.configfiles[0],
    log:
        "{out}/logs/multi_comparison_overlap.log"
    shell:
        "Rscript src/analysis/09_peak_overlap.R "
        "{params.config} {params.out_dir} "
        "> {log} 2>&1"


# ──────────────────────────────────────────────
# Rule 14: Time-course temporal analysis (optional)
# ──────────────────────────────────────────────
rule timecourse_analysis:
    input:
        da_results  = expand("{{out}}/pairwise/{pair}/final_da_results.csv", pair=PAIRS),
        norm_counts = expand("{{out}}/pairwise/{pair}/normalized_counts.csv", pair=PAIRS),
    output:
        flag = touch("{out}/timecourse/.timecourse_done.flag"),
    params:
        out_dir = "{out}/timecourse",
        config  = lambda wildcards: workflow.configfiles[0],
    log:
        "{out}/logs/timecourse.log"
    shell:
        "Rscript src/analysis/10_timecourse_analysis.R "
        "{params.config} {params.out_dir} "
        "> {log} 2>&1"


# ──────────────────────────────────────────────
# Rule 6: Summary Excel table
# ──────────────────────────────────────────────
rule generate_da_summary_table:
    input:
        results = "{out}/pairwise/{pair}/final_da_results.csv",
        flag    = "{out}/pairwise/{pair}/.enrichment_done.flag",
    output:
        da_xlsx = "{out}/pairwise/{pair}/final_da_result.xlsx",
        go_xlsx = "{out}/pairwise/{pair}/final_go_result.xlsx",
    params:
        out_dir = "{out}/pairwise/{pair}",
        config  = lambda wildcards: workflow.configfiles[0],
    log:
        "{out}/logs/{pair}_summary_table.log"
    shell:
        "Rscript src/analysis/06_generate_da_table.R "
        "{params.config} {params.out_dir} > {log} 2>&1"


# ──────────────────────────────────────────────
# Rule 8: SeqViewer parquet export (optional)
# ──────────────────────────────────────────────
rule export_seqviewer:
    input:
        xlsx    = "{out}/pairwise/{pair}/final_da_result.xlsx",
        enr_flag= "{out}/pairwise/{pair}/.enrichment_done.flag",
    output:
        flag = touch("{out}/pairwise/{pair}/.seqviewer_export_done.flag"),
    params:
        compare = lambda wc: wc.pair.split("_vs_")[0],
        base    = lambda wc: wc.pair.split("_vs_")[1],
        out_dir = "{out}/pairwise/{pair}",
        config  = lambda wildcards: workflow.configfiles[0],
    log:
        "{out}/logs/{pair}_seqviewer.log"
    shell:
        "Rscript src/analysis/06_export_seqviewer.R "
        "{params.config} {params.compare} {params.base} {params.out_dir} "
        "> {log} 2>&1"


# ──────────────────────────────────────────────
# Rule 9: Global QC HTML report
# ──────────────────────────────────────────────
rule global_qc_report:
    input:
        flag = "{out}/global_qc/.global_qc_done.flag",
    output:
        html = "{out}/global_qc/global_qc_report.html",
    params:
        qc_dir = "{out}/global_qc",
        config = lambda wildcards: workflow.configfiles[0],
    log:
        "{out}/logs/global_qc_report.log"
    shell:
        "Rscript src/analysis/02d_generate_global_qc_report.R "
        "--config {params.config} "
        "--qc_plots_dir {params.qc_dir} "
        "--output_file {output.html} "
        "> {log} 2>&1"


# ──────────────────────────────────────────────
# Rule 10: Methods section Markdown
# ──────────────────────────────────────────────
rule generate_methods_section:
    input:
        # 모든 pairwise DA 완료 후 실행
        expand("{{out}}/pairwise/{pair}/final_da_result.xlsx", pair=PAIRS),
    output:
        md = "{out}/methods_section.md",
    params:
        out_dir = "{out}",
        config  = lambda wildcards: workflow.configfiles[0],
    log:
        "{out}/logs/methods_section.log"
    shell:
        "Rscript src/analysis/08_generate_methods_section.R "
        "--config {params.config} "
        "--output-dir {params.out_dir} "
        "--output {output.md} "
        "> {log} 2>&1"


# ──────────────────────────────────────────────
# Rule 11: SeqViewer aggregate manifest + ZIP
# ──────────────────────────────────────────────
rule aggregate_seqviewer:
    input:
        flags = expand(
            "{{out}}/pairwise/{pair}/.seqviewer_export_done.flag",
            pair=PAIRS
        ),
    output:
        manifest = "{out}/seqviewer/seqviewer_manifest.json",
    params:
        out_dir = "{out}",
        config  = lambda wildcards: workflow.configfiles[0],
    log:
        "{out}/logs/seqviewer_aggregate.log"
    shell:
        "Rscript src/analysis/06b_aggregate_seqviewer.R "
        "{params.config} {params.out_dir} "
        "> {log} 2>&1"


# ──────────────────────────────────────────────
# Rule 12: DA Summary HTML report
# ──────────────────────────────────────────────
rule generate_summary_report:
    input:
        # DA results + volcano plots for all pairs
        expand("{{out}}/pairwise/{pair}/final_da_results.csv", pair=PAIRS),
        expand("{{out}}/pairwise/{pair}/volcano_plot.png",     pair=PAIRS),
        expand("{{out}}/pairwise/{pair}/.enrichment_done.flag", pair=PAIRS),
        # methods section must exist first (ensures full pipeline finished)
        "{out}/methods_section.md",
        # Optional analyses — included when enabled so report reflects final state
        expand("{out}/chromvar/.chromvar_done.flag", out=OUT)
            if RUN_CHROMVAR else [],
        expand("{out}/multi_comparison/.overlap_done.flag", out=OUT)
            if RUN_OVERLAP else [],
        expand("{out}/timecourse/.timecourse_done.flag", out=OUT)
            if RUN_TIMECOURSE else [],
    output:
        html = "{out}/summary_report.html",
    params:
        out_dir = "{out}",
        config  = lambda wildcards: workflow.configfiles[0],
    log:
        "{out}/logs/summary_report.log"
    shell:
        "Rscript src/analysis/07_generate_summary_report.R "
        "--config {params.config} "
        "--output-dir {params.out_dir} "
        "--output {output.html} "
        "> {log} 2>&1"
