.libPaths("/home1/tanxinli/R_lib")
suppressPackageStartupMessages({
  library(SeqArray)
  library(SNPRelate)
  library(GENESIS)
  library(Biobase)
  library(data.table)
  library(SeqVarTools)
})


args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 6.assoc.R <chr> <pheno_name>")
}

chr   <- as.integer(args[1])
pheno <- args[2]  # 读取第二个参数作为 phenotype 名字
message("=== Start assoc for Chromosome: ", chr, " | Phenotype: ", pheno, " ===")

# --- 2. 动态设置路径 (根据你的目录结构) ---
# 基础目录 (所有结果的父目录)
base_dir <- "/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/result/null_model"

# 输出目录: base_dir 下面以 phenotype 命名的子文件夹
out_dir <- file.path(base_dir, pheno)

# Null Model 文件路径: 假设文件名格式为 "null_model_变量名.RDS"
nullmod_file <- file.path(base_dir, paste0("null_model_", pheno, ".RDS"))

# GDS 和 Annotation 文件 (这些通常是不变的)
gds_file   <- "/project2/desmith_128/tanxin/scratch1_new//mp2prt_gwas/filtered/split/MP2PRT.filtered.merged.gds"
annot_file <- "/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/annot_sorted.rds"

# 检查 Null Model 是否存在，不存在就报错停止
if (!file.exists(nullmod_file)) {
  stop("Null model file not found: ", nullmod_file)
}

# --- 3. 读取数据 ---
gds <- seqOpen(gds_file)
nullmod <- readRDS(nullmod_file)
seqData <- SeqVarData(gds)

# 更新 annot
annot <- readRDS(annot_file)
sampleData(seqData) <- annot


## 按 chr 过滤
seqResetFilter(seqData, verbose = FALSE)
chrom <- seqGetData(seqData, "chromosome")
sel <- chrom == chr
if (!any(sel)) {
  seqClose(seqData) # 记得关闭
  stop("No variants found on chromosome ", chr)
}
seqSetFilter(seqData, variant.sel = sel, verbose = TRUE)

## 设定 iterator
iterator <- SeqVarBlockIterator(seqData, verbose = TRUE)


## 跑 single-variant test 这里的 nullmod 已经包含了对应的 phenotype 信息
assoc_chr <- assocTestSingle(iterator, nullmod, verbose=FALSE,
                         BPPARAM=BiocParallel::SerialParam())

## 加上 chr 信息（如果没有的话）
if (!"chr" %in% names(assoc_chr)) {
  assoc_chr$chr <- chr
}

## 输出结果 自动创建对应 phenotype 的文件夹

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(out_dir, paste0("assoc_chr", chr, ".rds"))
saveRDS(assoc_chr, out_file)

message("Saved results to: ", out_file)

## 关掉 GDS
seqClose(seqData)
message("=== Done chr ", chr, " ===")

