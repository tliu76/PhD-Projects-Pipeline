#!/usr/bin/env Rscript

# 1. 接收命令行参数 (从 Bash 传进来的 phenotype)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Error: No phenotype argument provided. Usage: Rscript script.R <phenotype_name>")
}

pheno_name <- args[1]
message("Processing phenotype: ", pheno_name)

.libPaths("/home1/tanxinli/R_lib")

suppressPackageStartupMessages({
  library(SeqArray)
  library(SNPRelate)
  library(GENESIS)
  library(Biobase)
  library(data.table)
  library(dplyr)
})

# 2. 设置基础路径
base_dir <- "/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/result/null_model"

# 3. 动态构建该表型的具体路径
# 假设文件夹名称就是表型名称 (例如 .../result0115/n.all.del_INT)
assoc_dir <- file.path(base_dir, pheno_name)

# 检查文件夹是否存在
if (!dir.exists(assoc_dir)) {
  stop("Directory does not exist: ", assoc_dir)
}

# 4. 构建文件列表
files <- file.path(assoc_dir, paste0("assoc_chr", 1:22, ".rds"))
exists <- file.exists(files)

if (!all(exists)) {
  warning("Missing files for ", pheno_name, ": ", paste(files[!exists], collapse = ", "))
  # 如果你希望缺少文件时直接报错停止，取消下面这行的注释:
  # stop("Missing files detected.")
}

# 5. 读取并合并
# 只读取存在的文件
if (sum(exists) > 0) {
    assoc_list <- lapply(files[exists], readRDS)
    assoc_all <- bind_rows(assoc_list)
    
    # 6. 保存结果
    # 结果保存在该表型的子文件夹下
    out_file <- file.path(assoc_dir, "assoc_all_chr1_22.rds")
    saveRDS(assoc_all, out_file)
    
    message("Successfully combined ", sum(exists), " files.")
    message("Result saved to: ", out_file)
} else {
    message("No files found to combine for ", pheno_name)
}