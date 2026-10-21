#!/usr/bin/env Rscript

# ============================
# 2.ld_prune.R
# Per-chromosome LD pruning
# ============================

.libPaths("/home1/tanxinli/R_lib")

suppressPackageStartupMessages({
    if (requireNamespace("GenomeInfoDb", quietly = TRUE)) {
        library(GenomeInfoDb)
    }
  library(SeqArray)
  library(SNPRelate)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript 2.ld_prune.R <gds_dir> <out_dir> <chr>")
}



gds_dir  <- args[1]   # directory containing MP2PRT.filtered.chr*.gds
out_dir  <- args[2]   # directory to save ldpruned_chr*.RDS
chr      <- as.integer(args[3])

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# file path of this chr GDS
gds_file <- file.path(gds_dir, paste0("MP2PRT.filtered.chr", chr, ".gds"))

message("==== LD-pruning chr ", chr, " ====")
message("Input GDS: ", gds_file)

if (!file.exists(gds_file)) {
    stop("GDS file not found: ", gds_file)
}

# -----------------------
# 1. Open chromosome GDS
# -----------------------
gds <- seqOpen(gds_file)
on.exit(try(seqClose(gds), silent = TRUE), add = TRUE)

vid <- seqGetData(gds, "variant.id")
sample.id <- seqGetData(gds, "sample.id")
if (length(vid) == 0) {
    message("No variants found on chr", chr)
    quit(save = "no")
}

message("Variants on chr", chr, ": ", length(vid))

# -----------------------
# 2. LD pruning (chr-specific)
# -----------------------
set.seed(100) 
pruned_chr <- snpgdsLDpruning(
    gds,
    sample.id = sample.id,
    method        = "corr",  ## "composite"
    ld.threshold  = sqrt(0.1),  ## previous 0.1
    slide.max.bp  = 10e6,   ## previous 500000
    num.thread    = 1   # array job handles parallelism
)

pruned_vec <- unlist(pruned_chr, use.names = FALSE)

# -----------------------
# 3. Save output
# -----------------------
out_file <- file.path(out_dir, paste0("ldpruned_chr", chr, ".RDS"))
saveRDS(pruned_vec, out_file)

message("==== DONE chr ", chr, " | LD-pruned SNPs: ", length(pruned_vec), " ====")

