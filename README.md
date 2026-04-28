# ATAC-Seq DA Analysis Pipeline

ATAC-seq Differential Accessibility (DA) 3차 분석 파이프라인.  
nf-core/atacseq (2차 분석) 출력을 입력으로 받아 DA 분석 → 피크 어노테이션 → GO/KEGG → 모티프 분석 → TF activity → Peak overlap → Time-course 분석까지 수행합니다.

RNA-Seq_DE_GO_analysis 파이프라인과 동일한 Snakemake 구조를 채택합니다.

---

## 분석 흐름

```
nf-core/atacseq 출력
  └── featureCounts.txt (consensus peaks × samples)
  └── consensus_peaks.bed

          ↓

[전체 샘플 공통]
02a_global_qc.R           PCA, Sample distance heatmap, Top peaks heatmap

[비교군 별 (pairwise)]
01_run_pairwise_da.R      DESeq2 Differential Accessibility (LFC shrinkage 포함)
02c_pairwise_plots.R      Volcano, MA plot
03_peak_annotation_go.R   ChIPseeker → Genomic distribution + GO/KEGG (clusterProfiler)
04_go_barplots.R          GO barplots (ontology별 + combined)
05_motif_analysis.R       HOMER motif enrichment — known & de novo  [선택: run_motif]
06_generate_da_table.R    Excel 통합 요약

[전체 샘플 공통 — 선택 분석]
07_chromvar_analysis.R    chromVAR TF activity (전체 샘플 1회 계산 + 비교군별 Wilcoxon)
                          [선택: run_chromvar]
09_peak_overlap.R         Multi-comparison DA peak overlap: UpSet + Venn + GO
                          [선택: run_overlap]
10_timecourse_analysis.R  시계열 temporal clustering: Z-score → K-means → trend plot + GO
                          [선택: run_timecourse, 시계열 디자인 전용]

[리포트]
08_generate_methods_section.R   Methods 섹션 Markdown 자동 생성
07_generate_summary_report.R    DA 분석 통합 HTML 보고서
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

nf-core/atacseq 출력에서 다음 파일을 `data/raw/`에 복사하거나 경로를 config에 직접 지정합니다:

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

피크 어노테이션 정보 CSV (`peak_info.csv`) — 선택, GO enrichment에 필요:

```csv
peak_id,Chr,Start,End,gene_id,gene_name
Interval_1,chr1,3095447,3096447,ENSMUSG00000051951,Xkr4
```

> nf-core/atacseq의 featureCounts.txt 헤더 4개 컬럼(Chr, Start, End, Strand)에서 자동 생성 가능.

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
peak_info_path: "data/raw/peak_info.csv" # peak → gene 매핑 (GO enrichment에 사용)

da_analysis:
  pairwise_comparisons:
    - ["100uM", "CONTROL"]
    - ["200uM", "CONTROL"]

motif:
  run_motif: true          # HOMER 실행 여부
  homer_genome: "mm10"

chromvar:
  run_chromvar: true       # chromVAR TF activity 실행 여부 (전체 샘플 기준 1회)
  jaspar_collection: "CORE"
  top_tfs_n: 30
  n_workers: 4             # 병렬 처리 worker 수 (Linux/WSL)

peak_overlap:
  run_overlap: true        # Multi-comparison peak overlap 실행 여부
  min_peaks_for_go: 20     # GO 분석을 위한 최소 peak 수
  go_ontology: "BP"
  venn_subsets:            # Venn diagram 조합 (최대 4-way)
    - ["100uM_vs_CONTROL", "200uM_vs_CONTROL"]

timecourse:
  run_timecourse: true     # 시계열 분석 여부 (시계열 디자인 전용)
  time_order: ["CONTROL", "3H", "6H", "1D"]   # 시간순 condition 목록
  time_hours:  [0, 3, 6, 24]                   # 대응 시간(h)
  n_clusters: null         # null → elbow+silhouette 자동 선택, 숫자 → 지정
  cluster_method: "kmeans" # "kmeans" 또는 "hierarchical"
  go_ontology: "BP"
  min_peaks_for_go: 50
```

### 4. 실행

```bash
# Dry-run (실행 계획 확인)
snakemake --configfile configs/config_MyProject.yml --cores 8 --dryrun

# 실행
snakemake --configfile configs/config_MyProject.yml --cores 8

# 특정 분석만 실행
snakemake output/pairwise/100uM_vs_CONTROL/final_da_results.csv \
  --configfile configs/config_MyProject.yml --cores 4
```

---

## 선택 분석 상세

### chromVAR (TF Activity)

`run_chromvar: true` 설정 시 전체 샘플에 대해 한 번 실행됩니다.

- **입력**: raw featureCounts + consensus peaks BED + BSgenome
- **방법**: JASPAR2020 motif matching → GC bias 보정 → deviation score (z-score per motif per sample)
- **출력**:
  - `tf_variability.csv` — 전체 모티프 variability 점수 (전체 샘플 기준)
  - `tf_deviation_heatmap.png` — 상위 TF × 전체 샘플 z-score heatmap
  - `differential_tf/{pair}_diff_tf.csv` — 비교군 별 Wilcoxon 검정 결과 (delta z-score, padj)
  - `differential_tf/{pair}_diff_tf_plot.png` — 상위 differential TF barplot

### Peak Overlap (Multi-comparison)

`run_overlap: true` 설정 시 모든 pairwise 비교 완료 후 실행됩니다.

- **입력**: 각 비교군의 `final_da_results.csv`
- **방법**: DA 유의 peak (padj, LFC cutoff) 추출 → 비교군 간 교집합/합집합 분석
- **출력**:
  - `peak_membership_matrix.csv` — peak × comparison TRUE/FALSE 행렬
  - `overlap_summary.csv` — 교집합 패턴별 peak 수 요약
  - `upset_plot.png` — 전체 비교군 UpSet plot
  - `venn_*.png` — config의 venn_subsets 조합 Venn diagram (2~4-way)
  - `go_enrichment/group_*_BP.csv` — 교집합 패턴별 GO 분석 (min_peaks_for_go 이상)

### Time-course (시계열 분석)

`run_timecourse: true` 설정 시 실행. **시계열 디자인 전용**이며 `pym` 등 독립 처리군 디자인에는 `run_timecourse: false`로 설정합니다.

- **입력**: 모든 비교군의 DA peak union + 원본 featureCounts (전체 샘플 VST)
- **방법**: union DA peaks → 전체 샘플 VST 정규화 → condition별 평균 → Z-score → K-means/hierarchical clustering
- **K 자동 선택**: `n_clusters: null` 시 WSS elbow + silhouette 방법으로 최적 k 탐색 (k=3~10)
- **출력**:
  - `da_peaks_temporal_clusters.csv` — peak별 cluster 할당 + Z-score
  - `cluster_summary.csv` — cluster별 peak 수
  - `elbow_plot.png` / `silhouette_plot.png` — k 선택 진단 플롯 (n_clusters: null 시)
  - `temporal_heatmap.png` — 전체 peaks Z-score heatmap (cluster 순서 정렬)
  - `trend_cluster_{k}.png` — cluster별 시간에 따른 mean ± SE Z-score
  - `go_enrichment/cluster_{k}_BP.csv` — cluster별 GO 분석 (min_peaks_for_go 이상)

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
| Peak overlap | — | UpSet + Venn + GO (다중 비교군) |
| 시계열 분석 | — | Z-score K-means + trend plot + GO |
| 결과 파일 | final_de_results.csv | final_da_results.csv |
| 요약 Excel | final_go_results.xlsx | final_da_summary.xlsx |

---

## 출력 디렉토리 구조

```
output/{project}/
├── global_qc/
│   ├── .global_qc_done.flag
│   ├── pca_all_samples.png
│   ├── sample_distance_heatmap.png
│   ├── top_peaks_heatmap.png
│   └── global_qc_report.html
│
├── pairwise/
│   └── {compare}_vs_{base}/
│       ├── final_da_results.csv            # DESeq2 결과 (전체 peaks)
│       ├── normalized_counts.csv           # VST 정규화 counts (해당 비교군 샘플)
│       ├── volcano_plot.png
│       ├── ma_plot.png
│       ├── peak_distribution_{total|up|down}.png
│       ├── peak_annotation_{total|up|down}.csv
│       ├── go_enrichment_{geneset}_{onto}.csv
│       ├── kegg_enrichment_{geneset}.csv
│       ├── go_barplot_{onto}_combined.png
│       ├── motif/                          # [run_motif: true]
│       │   ├── up/knownResults.html
│       │   ├── up/homerResults.html
│       │   └── top_known_motifs_up.csv
│       └── final_da_summary.xlsx
│
├── chromvar/                               # [run_chromvar: true] 전체 샘플 기준 1회
│   ├── tf_variability.csv
│   ├── tf_variability_plot.png
│   ├── tf_deviation_heatmap.png
│   └── differential_tf/
│       ├── {pair}_diff_tf.csv
│       └── {pair}_diff_tf_plot.png
│
├── multi_comparison/                       # [run_overlap: true]
│   ├── peak_membership_matrix.csv
│   ├── overlap_summary.csv
│   ├── upset_plot.png
│   ├── venn_{subset}.png
│   └── go_enrichment/
│       └── group_{pattern}_BP.csv
│
├── timecourse/                             # [run_timecourse: true] 시계열 전용
│   ├── da_peaks_temporal_clusters.csv
│   ├── cluster_summary.csv
│   ├── elbow_plot.png
│   ├── silhouette_plot.png
│   ├── temporal_heatmap.png
│   ├── trend_cluster_{k}.png
│   └── go_enrichment/
│       └── cluster_{k}_BP.csv
│
├── methods_section.md                      # 분석 방법 자동 생성 (논문용)
└── summary_report.html                     # 전체 분석 통합 HTML 보고서
```

---

## 참고

- [DESeq2](https://bioconductor.org/packages/DESeq2/)
- [ChIPseeker](https://bioconductor.org/packages/ChIPseeker/)
- [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/)
- [HOMER](http://homer.ucsd.edu/homer/ngs/peakMotifs.html)
- [chromVAR](https://bioconductor.org/packages/chromVAR/)
- [JASPAR2020](https://bioconductor.org/packages/JASPAR2020/)
- [ComplexUpset](https://krassowski.github.io/complex-upset/)
- [factoextra](https://rpkgs.datanovia.com/factoextra/)
- 관련 파이프라인: `../atac-seq-pipeline/` (2차), `../RNA-Seq_DE_GO_analysis/` (RNA-seq 3차)
