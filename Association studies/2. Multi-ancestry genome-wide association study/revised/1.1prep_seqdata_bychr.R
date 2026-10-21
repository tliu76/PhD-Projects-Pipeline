#!/usr/bin/env Rscript

.libPaths("/home1/tanxinli/R_lib")
suppressPackageStartupMessages({
  library(gdsfmt)
  if (requireNamespace("GenomeInfoDb", quietly = TRUE)) {
    library(GenomeInfoDb)
  }
  suppressWarnings(library(SeqArray))
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript prep_seqdata_chr.R covar_csv gds_path out_dir chr")
}

covar_csv <- args[1]
gds_path  <- args[2]
out_dir   <- args[3]
chr       <- as.integer(args[4])

message("========== START: prep_seqdata for chr ", chr, " ==========")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(covar_csv))
stopifnot(file.exists(gds_path))

# -----------------------
# Load covariates
# -----------------------
covar <- fread(covar_csv)

# -----------------------
# Open GDS
# -----------------------
gds <- seqOpen(gds_path)
on.exit(try(seqClose(gds), silent = TRUE), add = TRUE)
all_ids <- seqGetData(gds, "sample.id")
all_chr <- seqGetData(gds, "chromosome")
all_varid <- seqGetData(gds, "variant.id")

# -----------------------
# Covariate alignment
# -----------------------
covar_full <- covar[match(all_ids, covar$sample.id), ]
covar_full$sample.id[is.na(covar_full$sample.id)] <- all_ids[is.na(covar_full$sample.id)]
stopifnot(identical(covar_full$sample.id, all_ids))

# -----------------------
# Extract specific chromosome
# -----------------------
chr_idx <- which(all_chr == chr | all_chr == paste0("chr", chr))
message("Chr ", chr, " variants: ", length(chr_idx))

if (length(chr_idx) == 0) {
  stop("No variants found for chr ", chr, " in ", gds_path)
}

seqSetFilter(gds, sample.id = all_ids,
             variant.id = all_varid[chr_idx], verbose = TRUE)

# -----------------------
# Export
# -----------------------
chr_out <- file.path(out_dir, paste0("MP2PRT.filtered.chr", chr, ".gds"))
if (file.exists(chr_out)) file.remove(chr_out)

message("Exporting -> ", chr_out)
seqExport(gds, chr_out)

