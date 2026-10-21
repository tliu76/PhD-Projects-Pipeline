#!/usr/bin/env Rscript
.libPaths("/home1/tanxinli/R_lib")
install.packages("BiocManager")
BiocManager::install(c("GENESIS", "SeqVarTools", "SeqArray", "SNPRelate"))


.libPaths("/home1/tanxinli/R_lib")
suppressPackageStartupMessages({
  library(GENESIS)
})

pcre_dir <- "/scratch1/tanxinli/mp2prt_gwas/filtered/split/split_by_chr"

message("Loading per-chromosome PC-Relate results...")

files <- file.path(pcre_dir, paste0("pcrelate_chr", 1:22, ".RDS"))
pcrel_list <- lapply(files, readRDS)

message("Merging PC-Relate...")
merged <- GENESIS::pcrelateMerge(pcrel_list)

out <- file.path(pcre_dir, "pcrelate_ALLCHR.RDS")
saveRDS(merged, out)

message("Done! Saved merged file to: ", out)


obj <- readRDS("/scratch1/tanxinli/mp2prt_gwas/filtered/split/split_by_chr/pcrelate_chr1.RDS")
str(obj$snpsUsed)


merged <- readRDS("/scratch1/tanxinli/mp2prt_gwas/filtered/split/split_by_chr/pcrelate_ALLCHR.RDS")
str(merged)
