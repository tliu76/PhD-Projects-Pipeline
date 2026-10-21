#!/usr/bin/env Rscript


# ===============================================================
# 4.run_pcrelate_full.R
# One-shot PC-Relate on merged GDS (no per-chromosome split)
# Inputs:
#   - MP2PRT.filtered.merged.gds
#   - ldpruned_allchr.RDS
#   - pcair_PCs_top10.csv
#   - pcair_unrels.txt
# Outputs:
#   - pcrelate_obj.RDS
#   - pcrelate_kinship_mat.RDS
# ===============================================================

.libPaths(c("/home1/tanxinli/R_lib", .libPaths()))
suppressPackageStartupMessages({
  library(SeqArray)
  library(SeqVarTools)
  library(GENESIS)
  library(Biobase)
  library(data.table)
  library(BiocParallel)
})


## -----------------------------------------------------------
## 1. 路径
## -----------------------------------------------------------
gds_file    <- "/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/split/MP2PRT.filtered.merged.gds"

pcair_dir   <- "/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/result/pcair"
pcs_file    <- file.path(pcair_dir, "pcair_PCs_top10.csv")
unrels_file <- file.path(pcair_dir, "pcair_unrels.txt")
pruned_file <- file.path(pcair_dir, "ldpruned_allchr.RDS")

out_dir     <- "/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/result/pcrelate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(gds_file))
stopifnot(file.exists(pcs_file))
stopifnot(file.exists(unrels_file))
stopifnot(file.exists(pruned_file))

## -----------------------------------------------------------
## 2. 打开 GDS & 读取样本
## -----------------------------------------------------------
message("Opening merged GDS: ", gds_file)
gds <- seqOpen(gds_file)
sample.id <- seqGetData(gds, "sample.id")
message("N samples in merged GDS: ", length(sample.id))

## -----------------------------------------------------------
## 3. 读取 PC-AiR PCs & unrels，并按 GDS 顺序对齐
## -----------------------------------------------------------
message("Reading PCs from: ", pcs_file)
pcs <- fread(pcs_file)

# 按 GDS 样本顺序重排
pcs <- pcs[match(sample.id, pcs$sample.id), ]
if (any(is.na(pcs$sample.id))) {
  stop("Some GDS sample.id not found in pcair_PCs_top10.csv")
}

# 只用 PC1 和 PC2 来做 PC-Relate（推荐做法）
pcMat <- as.matrix(pcs[, .(PC1, PC2)])
rownames(pcMat) <- pcs$sample.id

message("PC matrix dim: ", paste(dim(pcMat), collapse = " x "))

# unrels
unrels <- scan(unrels_file, what = "character")
unrels <- intersect(unrels, sample.id)
message("Unrelated training set (PC-AiR): ", length(unrels), " samples")

if (length(unrels) == 0) {
  stop("No unrelated samples found in pcair_unrels.txt")
}

## -----------------------------------------------------------
## 4. 读取 LD-pruned SNP（全基因组）
## -----------------------------------------------------------
pruned_all <- readRDS(pruned_file)

vid <- seqGetData(gds, "variant.id")
pruned <- intersect(vid, pruned_all)

message("# LD-pruned variants present in merged GDS: ", length(pruned))
if (length(pruned) == 0) stop("No LD-pruned variants found in merged GDS!")

## -----------------------------------------------------------
## 5. 构建 SeqVarData + 过滤到 pruned SNPs
## -----------------------------------------------------------
message("Building SeqVarData and applying variant filter...")

# 用 PCs 作为 sampleData（后面 fitNullModel 也很方便）
annot <- AnnotatedDataFrame(pcs)  # pcs 里有 sample.id + PC1...PC10
seqData <- SeqVarData(gds, sampleData = annot)

seqSetFilter(seqData, variant.id = pruned, verbose = TRUE)

## -----------------------------------------------------------
## 6. 设置并行（按 block 迭代 SNP）
## -----------------------------------------------------------
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
message("Using ", ncores, " cores for PC-Relate (BiocParallel::MulticoreParam)")
bp <- MulticoreParam(ncores)

iterator <- SeqVarBlockIterator(seqData, variantBlock = 30000, verbose = FALSE)  #Increase block size to 30k (or even 50k) → fewer blocks → less overhead → faster runtime.

## -----------------------------------------------------------
## 7. 跑 PC-Relate（一次性全基因组）
## -----------------------------------------------------------
message("Running PC-Relate on merged GDS...")

pcrel <- pcrelate(
  iterator,
  pcs          = pcMat,          # 只用 PC1、PC2 调整族群结构
  training.set = unrels,         # PC-AiR 得到的 unrelated 样本
  BPPARAM      = bp
)

pcrel_out <- file.path(out_dir, "pcrelate_obj.RDS")
saveRDS(pcrel, pcrel_out)
message("Saved PC-Relate object to: ", pcrel_out)

## -----------------------------------------------------------
## 8. 转换为 kinship matrix（给 fitNullModel 用）
## -----------------------------------------------------------
message("Converting PC-Relate output to kinship matrix (pcrelateToMatrix)...")

kinMat <- pcrelateToMatrix(pcrel)

kin_out <- file.path(out_dir, "pcrelate_kinship_mat.RDS")
saveRDS(kinMat, kin_out)

message("Kinship matrix dim: ", paste(dim(kinMat), collapse = " x "))
message("Saved kinship matrix to: ", kin_out)

## -----------------------------------------------------------
## 9. 收尾
## -----------------------------------------------------------
seqClose(seqData)


message("==== DONE: PC-Relate full-genome run completed ====")
message("Outputs:")
message("  - ", pcrel_out)
message("  - ", kin_out)