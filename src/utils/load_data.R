#!/usr/bin/env Rscript
# src/utils/load_data.R
# 데이터 로딩 유틸리티 - RNA-Seq_DE_GO_analysis 패턴 계승
# ATAC-seq featureCounts 형식 특화

library(dplyr)
library(tibble)

#' featureCounts.txt 형식의 count matrix 로딩
#'
#' nf-core/atacseq 출력: 첫 줄은 주석(#), 두 번째 줄이 헤더
#' 컬럼 구조: Geneid, Chr, Start, End, Strand, Length, [sample BAM paths...]
#'
#' @param count_path featureCounts.txt 파일 경로
#' @param metadata  metadata data.frame (행=sample_id)
#' @return list(count_matrix, peak_info)
load_featurecounts <- function(count_path, metadata = NULL) {
  message("Loading featureCounts: ", count_path)

  raw <- read.table(count_path, header = TRUE, sep = "\t",
                    skip = 1, check.names = FALSE, comment.char = "#")

  # 메타 컬럼 (featureCounts 고정 6개)
  meta_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  peak_info <- raw[, meta_cols, drop = FALSE]

  # Count matrix: 7번째 컬럼부터
  count_matrix <- as.matrix(raw[, 7:ncol(raw), drop = FALSE])
  rownames(count_matrix) <- raw$Geneid

  # 컬럼명 정리: BAM 파일 경로 → 샘플 ID
  colnames(count_matrix) <- clean_sample_names(colnames(count_matrix))

  # metadata가 있으면 컬럼 순서 맞추기
  if (!is.null(metadata)) {
    shared <- intersect(colnames(count_matrix), rownames(metadata))
    if (length(shared) == 0) {
      stop("Count matrix 컬럼명과 metadata rowname이 하나도 일치하지 않습니다.\n",
           "Count matrix 샘플: ", paste(head(colnames(count_matrix), 5), collapse=", "), "\n",
           "Metadata 샘플: ",     paste(head(rownames(metadata), 5), collapse=", "))
    }
    count_matrix <- count_matrix[, shared, drop = FALSE]
    message(sprintf("매칭된 샘플: %d개", length(shared)))
  }

  message(sprintf("Peaks: %d, Samples: %d", nrow(count_matrix), ncol(count_matrix)))
  list(count_matrix = count_matrix, peak_info = peak_info)
}


#' BAM 파일 경로에서 샘플 ID 추출
#'
#' nf-core/atacseq BAM 이름 패턴:
#'   /path/to/SAMPLE_ID.mLb.clN.sorted.bam
#'   /path/to/SAMPLE_ID.mRp.clN.sorted.bam
#'
#' @param names character vector of BAM file paths / names
#' @return cleaned sample ID vector
clean_sample_names <- function(names) {
  names <- basename(names)
  # nf-core suffix 제거
  names <- sub("\\.mLb\\.clN\\.sorted\\.bam$", "", names)
  names <- sub("\\.mRp\\.clN\\.sorted\\.bam$", "", names)
  names <- sub("\\.sorted\\.bam$", "", names)
  names <- sub("\\.bam$", "", names)
  names
}


#' 메타데이터 CSV 로딩
#'
#' @param metadata_path CSV 파일 경로 (첫 번째 컬럼이 sample_id)
#' @return data.frame (행=sample_id)
load_metadata <- function(metadata_path) {
  message("Loading metadata: ", metadata_path)
  meta <- read.csv(metadata_path, check.names = FALSE)

  # 첫 번째 컬럼을 row name으로
  rownames(meta) <- meta[[1]]
  meta <- meta[, -1, drop = FALSE]
  meta
}


#' DA 결과 CSV 로딩
#'
#' @param results_path final_da_results.csv 경로
#' @return data.frame
load_da_results <- function(results_path) {
  message("Loading DA results: ", results_path)
  res <- read.csv(results_path, check.names = FALSE)
  res
}


#' DA peak 분류: up / down / total
#'
#' @param da_results  load_da_results() 결과
#' @param padj_cutoff adjusted p-value 임계값
#' @param lfc_cutoff  |log2FoldChange| 임계값
#' @return list(total, up, down) - peak_id 벡터
classify_da_peaks <- function(da_results, padj_cutoff = 0.05, lfc_cutoff = 1.0) {
  sig <- da_results %>%
    filter(!is.na(padj), padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)

  list(
    total = sig$peak_id,
    up    = sig$peak_id[sig$log2FoldChange > 0],
    down  = sig$peak_id[sig$log2FoldChange < 0]
  )
}
