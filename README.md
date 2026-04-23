# ATAC-Seq DA Analysis Pipeline

ATAC-seq Differential Accessibility (DA) 3차 분석 파이프라인.
nf-core/atacseq (2차 분석) 출력을 입력으로 받아 DA 분석 → 피크 어노테이션 → GO/KEGG → 모티프 분석 → TF activity 분석까지 수행합니다.

RNA-Seq_DE_GO_analysis 파이프라인과 동일한 Snakemake 구조를 채택합니다.

---

## 분석 흐름

```
nf-core/atacseq 출력
  └── featureCounts.txt (consensus peaks × samples)
  └── consensus_peaks.bed

          ↓

01_run_pairwise_da.R      DESeq2 Differential Accessibility (LFC shrinkage 포함)
02a_global_qc.R           PCA, Sample distance heatmap, Top peaks heatmap (전체 샘플)
02c_pairwise_plots.R      Volcano, MA plot (비교별)
03_peak_annotation_go.R   ChIPseeker → Genomic distribution + GO/KEGG (clusterProfiler)
04_go_barplots.R          GO barplots (ontology별 + combined)
05_motif_analysis.R       HOMER motif enrichment — known & de novo (선택)
06_generate_da_table.R    Excel 통합 요약
07_chromvar_analysis.R    chromVAR TF activity deviation score (선택)

          ↓

output/{pair}/
  ├── final_da_results.csv
  ├── normalized_counts.csv
  ├── volcano_plot.png / ma_plot.png
  ├── peak_distribution_{total|up|down}.png
  ├── peak_annotation_{total|up|down}.csv
  ├── go_enrichment_*.csv / kegg_enrichment_*.csv
  ├── motif/up/ , motif/down/
  ├── chromvar/
  └── final_da_summary.xlsx
```

---

## 시작하기

### 1. 환경 설정

```bash
conda env create -f environment.yml
conda activate atac-seq-da-analysis

# HOMER genome 설치 (run_motif: true 사용 시, 최초 1회)
perl $(which configureHomer.pl) -install mm10
```

### 2. 데이터 준비

nf-core/atacseq 출력에서 다음 파일을 `data/raw/`에 복사하거나 경로를 config에 지정합니다:

| 파일 | nf-core/atacseq 경로 |
|---|---|
| `consensus_peaks.featureCounts.txt` | `results/.../consensus/*.featureCounts.txt` |
| `consensus_peaks.bed` | `results/.../consensus/*.boolean.bed` |

메타데이터 CSV (`data/raw/metadata.csv`) 형식:

```csv
sample_id,condition,batch
CONTROL_REP1,CONTROL,1
H2O2_100uM_REP1,100uM,1
H2O2_200uM_REP1,200uM,1
```

> **주의**: `sample_id`가 featureCounts.txt의 BAM 컬럼명과 일치해야 합니다.
> BAM 이름 패턴 `SAMPLE_ID.mLb.clN.sorted.bam` → `SAMPLE_ID` 자동 추출.

### 3. Config 작성

```bash
cp configs/template/config.yml configs/config_MyProject.yml
# config_MyProject.yml 편집
```

주요 설정 항목:

```yaml
species: "mouse"                         # human / mouse
count_matrix_path: "data/raw/..."        # featureCounts.txt 경로
consensus_peaks_path: "data/raw/..."     # consensus peaks BED
metadata_path: "data/raw/metadata.csv"

da_analysis:
  pairwise_comparisons:
    - ["100uM", "CONTROL"]

motif:
  run_motif: true          # HOMER 실행 여부
  homer_genome: "mm10"

chromvar:
  run_chromvar: true       # chromVAR 실행 여부
  jaspar_collection: "CORE"
  top_tfs_n: 30
  n_workers: 4             # 병렬 처리 worker 수 (Linux/WSL)
```

### 4. 실행

```bash
# Dry-run (실행 계획 확인)
snakemake --configfile configs/config_MyProject.yml --cores 8 --dryrun

# 실행
snakemake --configfile configs/config_MyProject.yml --cores 8

# 특정 비교만 실행
snakemake output/pairwise/100uM_vs_CONTROL/final_da_results.csv \
  --configfile configs/config_MyProject.yml --cores 4
```

---

## RNA-seq 파이프라인과의 차이점

| 항목 | RNA-Seq_DE_GO_analysis | ATAC-Seq_DA_analysis |
|---|---|---|
| 입력 count | gene counts CSV | featureCounts.txt (peak×sample) |
| 분석 단위 | 유전자 | genomic peak region |
| DE/DA 도구 | DESeq2 / edgeR / limma | DESeq2 |
| GO/KEGG | gene list 직접 사용 | peak → ChIPseeker → gene → GO |
| Genomic distribution | — | Promoter/Exon/Intron 등 분포 시각화 |
| Motif enrichment | — | HOMER known & de novo |
| TF activity | — | chromVAR deviation score (JASPAR2020) |
| 결과 파일 | final_de_results.csv | final_da_results.csv |
| 요약 Excel | final_go_results.xlsx | final_da_summary.xlsx |

---

## 출력 디렉토리 구조

```
output/
├── global_qc/
│   ├── pca_all_samples.png
│   ├── sample_distance_heatmap.png
│   └── top_peaks_heatmap.png
└── pairwise/
    └── 100uM_vs_CONTROL/
        ├── final_da_results.csv            # DESeq2 결과 (전체 peaks)
        ├── normalized_counts.csv           # VST 정규화 counts
        ├── config_used.yml                 # 재현성을 위한 config 복사
        ├── volcano_plot.png
        ├── ma_plot.png
        ├── peak_distribution_total.png     # 전체 DA peaks genomic 분포
        ├── peak_distribution_up.png        # More-open peaks 분포
        ├── peak_distribution_down.png      # Less-open peaks 분포
        ├── peak_annotation_total.csv       # ChIPseeker 상세 어노테이션
        ├── peak_annotation_up.csv
        ├── peak_annotation_down.csv
        ├── go_enrichment_up_BP.csv         # GO (up peaks, Biological Process)
        ├── go_dotplot_up_BP.png
        ├── kegg_enrichment_up.csv
        ├── go_barplot_up_combined.png
        ├── motif/                          # HOMER 결과 (run_motif: true 시)
        │   ├── up/knownResults.html
        │   ├── up/homerResults.html
        │   ├── down/knownResults.html
        │   └── top_known_motifs_up.csv
        ├── chromvar/                       # chromVAR 결과 (run_chromvar: true 시)
        │   ├── tf_variability.csv          # 전체 motif variability 점수
        │   ├── tf_variability_plot.png     # 상위 가변 TF barplot
        │   ├── tf_deviation_heatmap.png    # TF z-score heatmap
        │   └── differential_tf_activity.csv  # 그룹 간 Wilcoxon 검정 결과
        └── final_da_summary.xlsx           # 통합 Excel 보고서
```

---

## NGS Agent 연동

nf-core/atacseq 파이프라인의 `manifest.json`에서 입력 경로를 자동 파악할 수 있습니다:

```python
import json

manifest = json.load(open("path/to/sample/atac-seq/metadata/manifest.json"))
count_matrix = manifest["final_outputs"]["counts"]["featurecounts"]
consensus_bed = manifest["final_outputs"]["peaks"]["consensus_bed"]
```

---

## 참고

- [ChIPseeker](https://bioconductor.org/packages/ChIPseeker/)
- [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/)
- [HOMER](http://homer.ucsd.edu/homer/ngs/peakMotifs.html)
- [chromVAR](https://bioconductor.org/packages/chromVAR/)
- [JASPAR2020](https://bioconductor.org/packages/JASPAR2020/)
- 관련 파이프라인: `../atac-seq-pipeline/` (2차), `../RNA-Seq_DE_GO_analysis/` (RNA-seq 3차)
