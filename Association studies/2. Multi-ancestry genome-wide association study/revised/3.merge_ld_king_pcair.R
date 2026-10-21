#!/usr/bin/env Rscript

# ===============================================================
# 3.run_pcair_merge.R
# Merge LD-pruned SNPs + KING + PC-AiR
# Outputs: PCs, pcair_obj, unrels (needed for PC-Relate)
# ===============================================================

.libPaths("/home1/tanxinli/R_lib")
suppressPackageStartupMessages({
  if (requireNamespace("GenomeInfoDb", quietly = TRUE)) {
    library(GenomeInfoDb)
  }
  library(SeqArray)
  library(SNPRelate)
  library(GENESIS)
  library(Biobase)
  library(SeqVarTools)
  library(data.table)
})

message("==== START: PC-AiR merge step ====")

# ---------------------------------------------------------------
# 1. Paths
# ---------------------------------------------------------------
gds_merged <- "/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/split/MP2PRT.filtered.merged.gds"
pruned_dir <- "/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/result/pcair/pruned_by_chr"
out_dir    <- "/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/result/pcair"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------
# 2. Merge LD-pruned variant lists
# ---------------------------------------------------------------
message("Merging LD-pruned SNPs from chr1–22...")

pruned_all <- c()
for (chr in 1:22) {
  f <- file.path(pruned_dir, paste0("ldpruned_chr", chr, ".RDS"))
  if (!file.exists(f)) next
  v <- readRDS(f)
  pruned_all <- c(pruned_all, v)
}
pruned_all <- unique(pruned_all)

message("Total LD-pruned SNP count: ", length(pruned_all))
saveRDS(pruned_all, file.path(out_dir, "ldpruned_allchr.RDS"))

# ---------------------------------------------------------------
# 3. KING robust kinship
# ---------------------------------------------------------------
message("Running KING robust...")

gds <- seqOpen(gds_merged)
on.exit(try(seqClose(gds), silent = TRUE), add = TRUE)
sample.id <- seqGetData(gds, "sample.id")

king <- snpgdsIBDKING(
  gds,
  snp.id     = pruned_all,
  num.thread = 16
)

kingMat <- king$kinship
dimnames(kingMat) <- list(king$sample.id, king$sample.id)
saveRDS(kingMat, file.path(out_dir, "kingMat.RDS"))

# ---------------------------------------------------------------
# 4. PC-AiR
# ---------------------------------------------------------------
message("Running PC-AiR...")

pca <- pcair(
  gds,
  kinobj         = kingMat,
  divobj         = kingMat,
  kin.thresh     = 2^(-9/2),
  div.thresh     = -2^(-9/2),
  sample.include = sample.id,
  snp.include    = pruned_all,
  verbose        = TRUE
)

seqClose(gds)
saveRDS(pca, file.path(out_dir, "pcair_obj.RDS"))

# ---------------------------------------------------------------
# 5. Output principal components
# ---------------------------------------------------------------
pcs_df <- data.frame(pca$vectors)
colnames(pcs_df) <- paste0("PC", seq_len(ncol(pcs_df)))
pcs_df$sample.id <- pca$sample.id

pcs10 <- pcs_df[, c("sample.id", paste0("PC", 1:10))]
fwrite(pcs10, file.path(out_dir, "pcair_PCs_top10.csv"))



# ---------------------------------------------------------------
# 6. Output unrels (REQUIRED for PC-Relate)
# ---------------------------------------------------------------
unrels <- pca$unrels
rels   <- pca$rels

if (is.null(unrels)) stop("ERROR: PC-AiR did not produce unrels. Something is wrong.")

write.table(unrels,
            file.path(out_dir, "pcair_unrels.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(rels,
            file.path(out_dir, "pcair_rels.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

message("Unrelated samples (unrels): ", length(unrels))
message("Related samples (rels):     ", length(rels))

message("==== DONE: PC-AiR completed ====")
message("Outputs:")
message(" - ldpruned_allchr.RDS")
message(" - kingMat.RDS")
message(" - pcair_obj.RDS")
message(" - pcair_PCs_top10.csv")
message(" - pcair_unrels.txt  (for PC-Relate)")
message(" - pcair_rels.txt")