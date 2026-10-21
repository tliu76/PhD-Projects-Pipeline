#!/bin/bash

mkdir -p prs_score_files

# 使用grep更简单直接
exposures=("Baso.z" "Eosin.z" "LMR" "Lymph.z" "Mono.z" "Neutro.z" "NLR" "Platelet.z" "PLR")

for exposure_name in "${exposures[@]}"; do
    echo "Processing $exposure_name..."
    
    # 方法1: 使用grep和awk组合
    grep -w "$exposure_name" MR_Blood_Cell_Instruments_Kachuri_PMID34469753.txt | \
    awk '$13=="TRUE" {print $1, $4, $5, $8, $10}' \
    > prs_score_files/${exposure_name}_score.txt
    
    # 输出统计信息
    n_snps=$(wc -l < prs_score_files/${exposure_name}_score.txt)
    echo "  Extracted $n_snps SNPs for $exposure_name"
done

echo "All extractions completed!"