#!/usr/bin/env Rscript

.libPaths("/home1/tanxinli/R_lib")
suppressPackageStartupMessages({
  library(SeqArray)
  library(SNPRelate)
  library(GENESIS)
  library(Biobase)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(qqman)
  library(GenomicRanges)
  library(ChIPseeker)
  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

# ============================================
# 1. Read trait name
# ============================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("No phenotype argument provided.")
trait <- args[1]

cat(">>> Running plots for phenotype:", trait, "\n")

# ============================================
# 1.5 Phenotype labels (human-readable)
# ============================================
### NEW
pheno_labels <- c(
"n.all_INT"   ="Total number of SVs",                 
"n.nonig_INT" ="Total number of non-Ig/TCR SVs" ,          
"n.ig_INT"  ="Total number of Ig/TCR SVs" ,                  

"n.total.least.rag_INT"="Number of total RAG-mediated SVs (≥1 breakpoint)",
"n.ig.least.rag_INT"="Number of Ig/TCR RAG-mediated SVs (≥1 breakpoint)",
"n.nonig.least.rag_INT"="Number of non-Ig/TCR RAG-mediated SVs (≥1 breakpoint)",
"prop_total_least_rag.sv_INT"="Proportion of total RAG-mediated SVs (≥1 breakpoint)",
"prop_ig_least_rag.sv_INT"="Proportion of Ig/TCR RAG-mediated SVs (≥1 breakpoint)",
"prop_nonig_least_rag.sv_INT"="Proportion of non-Ig/TCR RAG-mediated SVs (≥1 breakpoint)",

"n_raghit_only_one_bkp_INT"="Number of SVs with RSS at only one breakpoint",
"n.ig.rag_only_INT"="Number of Ig/TCR SVs with RSS at only one breakpoint",
"n.nonig.rag_only_INT"="Number of non-Ig/TCR SVs with RSS at only one breakpoint",
"prop_raghit_only_one_bkp_INT"="Proportion of SVs with RSS at only one breakpoint",
"prop_ig_rag_only_INT"="Proportion of Ig/TCR SVs with RSS at only one breakpoint",
"prop_nonig_rag_only_INT"="Proportion of non-Ig/TCR SVs with RSS at only one breakpoint",

"n.all.del_INT"="Number of total deletions",
"n.total.least.rag.del_INT"="Number of total RAG-mediated deletions (≥1 breakpoint)",
"n.ig.least.rag.del_INT"="Number of Ig/TCR RAG-mediated deletions (≥1 breakpoint)",
"n.nonig.least.rag.del_INT"="Number of non-Ig/TCR RAG-mediated deletions (≥1 breakpoint)",
"prop_total_least_rag.del_INT"="Proportion of total RAG-mediated deletions (≥1 breakpoint)",
"prop_ig_least_rag.del_INT"="Proportion of Ig/TCR RAG-mediated deletions (≥1 breakpoint)",
"prop_nonig_least_rag.del_INT"="Proportion of non-Ig/TCR RAG-mediated deletions (≥1 breakpoint)",

"n_raghit_only_del_INT"="Number of total RAG-mediated deletions with RSS at only one breakpoint",
"n.ig.rag_only.del_INT"="Number of Ig/TCR RAG-mediated deletions with RSS at only one breakpoint",
"n.nonig.rag_only.del_INT"="Number of non-Ig/TCR RAG-mediated deletions with RSS at only one breakpoint",
"prop_raghit_only_del_INT"="Proportion of total RAG-mediated deletions with RSS at only one breakpoint",
"prop_ig_rag_only_del_INT"="Proportion of Ig/TCR RAG-mediated deletions with RSS at only one breakpoint",
"prop_nonig_rag_only_del_INT"="Proportion of non-Ig/TCR RAG-mediated deletions with RSS at only one breakpoint",

"n.all.dup_INT"="Number of total duplications",
"n.total.least.rag.dup_INT"="Number of total RAG-mediated duplications (≥1 breakpoint)",
"n.ig.least.rag.dup_INT"="Number of Ig/TCR RAG-mediated duplications (≥1 breakpoint)",
"n.nonig.least.rag.dup_INT"="Number of non-Ig/TCR RAG-mediated duplications (≥1 breakpoint)",
"prop_total_least_rag.dup_INT"="Proportion of total RAG-mediated duplications (≥1 breakpoint)",
"prop_ig_least_rag.dup_INT"="Proportion of Ig/TCR RAG-mediated duplications (≥1 breakpoint)",
"prop_nonig_least_rag.dup_INT"="Proportion of non-Ig/TCR RAG-mediated duplications (≥1 breakpoint)",

"n_raghit_only_dup_INT"="Number of total RAG-mediated duplications with RSS at only one breakpoint",
"n.ig.rag_only.dup_INT"="Number of Ig/TCR RAG-mediated duplications with RSS at only one breakpoint",
"n.nonig.rag_only.dup_INT"="Number of non-Ig/TCR RAG-mediated duplications with RSS at only one breakpoint",
"prop_raghit_only_dup_INT"="Proportion of total RAG-mediated duplications with RSS at only one breakpoint",
"prop_ig_rag_only_dup_INT"="Proportion of Ig/TCR RAG-mediated duplications with RSS at only one breakpoint",
"prop_nonig_rag_only_dup_INT"="Proportion of non-Ig/TCR RAG-mediated duplications with RSS at only one breakpoint",

"n.all.inv_INT"="Number of total inversions",
"n.total.least.rag.inv_INT"="Number of total RAG-mediated inversions (≥1 breakpoint)",
"n.ig.least.rag.inv_INT"="Number of Ig/TCR RAG-mediated inversions (≥1 breakpoint)",
"n.nonig.least.rag.inv_INT"="Number of non-Ig/TCR RAG-mediated inversions (≥1 breakpoint)",
"prop_total_least_rag.inv_INT"="Proportion of total RAG-mediated inversions (≥1 breakpoint)",
"prop_ig_least_rag.inv_INT"="Proportion of Ig/TCR RAG-mediated inversions (≥1 breakpoint)",
"prop_nonig_least_rag.inv_INT"="Proportion of non-Ig/TCR RAG-mediated inversions (≥1 breakpoint)",

"n_raghit_only_inv_INT"="Number of total RAG-mediated inversions with RSS at only one breakpoint",
"n.ig.rag_only.inv_INT"="Number of Ig/TCR RAG-mediated inversions with RSS at only one breakpoint",
"n.nonig.rag_only.inv_INT"="Number of non-Ig/TCR RAG-mediated inversions with RSS at only one breakpoint",
"prop_raghit_only_inv_INT"="Proportion of total RAG-mediated inversions with RSS at only one breakpoint",
"prop_ig_rag_only_inv_INT"="Proportion of Ig/TCR RAG-mediated inversions with RSS at only one breakpoint",
"prop_nonig_rag_only_inv_INT"="Proportion of non-Ig/TCR RAG-mediated inversions with RSS at only one breakpoint",

"n.all.tra_INT"="Number of total translocations",
"n.total.least.rag.tra_INT"="Number of total RAG-mediated translocations (≥1 breakpoint)",
"n.ig.least.rag.tra_INT"="Number of Ig/TCR RAG-mediated translocations (≥1 breakpoint)",
"n.nonig.least.rag.tra_INT"="Number of non-Ig/TCR RAG-mediated translocations (≥1 breakpoint)",
"prop_total_least_rag.tra_INT"="Proportion of total RAG-mediated translocations (≥1 breakpoint)",
"prop_ig_least_rag.tra_INT"="Proportion of Ig/TCR RAG-mediated translocations (≥1 breakpoint)",
"prop_nonig_least_rag.tra_INT"="Proportion of non-Ig/TCR RAG-mediated translocations (≥1 breakpoint)",

"n_raghit_only_tra_INT"="Number of total RAG-mediated translocations with RSS at only one breakpoint",
"n.ig.rag_only.tra_INT"="Number of Ig/TCR RAG-mediated translocations with RSS at only one breakpoint",
"n.nonig.rag_only.tra_INT"="Number of non-Ig/TCR RAG-mediated translocations with RSS at only one breakpoint",
"prop_raghit_only_tra_INT"="Proportion of total RAG-mediated translocations with RSS at only one breakpoint",
"prop_ig_rag_only_tra_INT"="Proportion of Ig/TCR RAG-mediated translocations with RSS at only one breakpoint",
"prop_nonig_rag_only_tra_INT"="Proportion of non-Ig/TCR RAG-mediated translocations with RSS at only one breakpoint"
)

plot_label <- pheno_labels[trait]

if (is.na(plot_label)) {
  stop("No label found for phenotype: ", trait)
}

# ============================================
# 2. Paths
# ============================================
base_dir <- "/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/result/null_model"
assoc_file <- file.path(base_dir, trait, "assoc_all_chr1_22.rds")
outdir <- file.path(base_dir, trait, "plots2")

if (!file.exists(assoc_file)) {
  stop("Assoc RDS file does not exist: ", assoc_file)
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================================
# 3. Load assoc results
# ============================================
assoc <- readRDS(assoc_file) %>%
  mutate(
    chr = as.numeric(as.character(chr)),
    pos = as.numeric(pos),
    Score.pval = as.numeric(Score.pval)
  ) %>%
  filter(!is.na(chr), !is.na(pos), !is.na(Score.pval))
  
# 确保列类型正确
assoc$chr <- as.numeric(as.character(assoc$chr))
assoc$pos <- as.numeric(assoc$pos)
assoc$Score.pval <- as.numeric(assoc$Score.pval)
# qqman::manhattan() requires an SNP column
assoc$SNP <- paste0("chr", assoc$chr, ":", assoc$pos)


# ============================================
# 4. QQ Plot (Custom, with lambda GC)
# ============================================
qqPlot_custom <- function(pval, trait_label="") {
  pval <- pval[!is.na(pval)]
  n <- length(pval)

  x <- 1:n
  chisq <- qchisq(1 - pval, df = 1)
  lambda_gc <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
  lambda_text <- paste0("lambda GC = ", round(lambda_gc, 3))

  dat <- data.frame(
    obs   = sort(pval),
    exp   = x / n,
    upper = qbeta(0.025, x, rev(x)),
    lower = qbeta(0.975, x, rev(x))
  )

  ggplot(dat, aes(-log10(exp), -log10(obs))) +
    geom_line(aes(-log10(exp), -log10(upper)), color = "gray") +
    geom_line(aes(-log10(exp), -log10(lower)), color = "gray") +
    geom_point(size = 1, alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = lambda_text,
      hjust = -0.2, vjust = 2,
      size = 5, fontface = "bold"
    ) +
    xlab(expression(paste(-log[10], "(expected P)"))) +
    ylab(expression(paste(-log[10], "(observed P)"))) +
    ggtitle(trait_label) +   # <<< human-readable label
    theme_bw() +
    theme(panel.grid = element_blank())
}

cat("Generating QQ plot...\n")
qq_fig <- qqPlot_custom(assoc$Score.pval, trait_label = plot_label)

qq_file <- file.path(outdir, paste0("qqplot_", trait, ".png"))
ggsave(qq_file, qq_fig, width = 50, height = 50, dpi = 600)


# ============================================
# 5. SNP → nearest gene (TxDb)
# ============================================
snp_gr <- GRanges(
  seqnames = paste0("chr", assoc$chr),
  ranges   = IRanges(start = assoc$pos, end = assoc$pos)
)

genes_gr <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
nearest_idx <- nearest(snp_gr, genes_gr)

assoc$gene <- mapIds(
  org.Hs.eg.db,
  keys = genes_gr$gene_id[nearest_idx],
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

# ============================================
# 6. Gene-level summary (top 5 genes)
# ============================================
gene_summary <- assoc %>%
  filter(!is.na(gene)) %>%
  group_by(gene) %>%
  summarise(
    min_p = min(Score.pval, na.rm = TRUE),
    chr   = chr[which.min(Score.pval)],
    pos   = pos[which.min(Score.pval)],
    .groups = "drop"
  )

top5_genes <- gene_summary %>%
  arrange(min_p) %>%
  slice_head(n = 5)

cat("Top 5 genes:\n")
print(top5_genes)

# ============================================
# 5. Manhattan Plot with Gene Annotation
# ============================================
cat("Generating Manhattan plot...\n")

man_file <- file.path(outdir, paste0("manhattan_", trait, ".png"))

suggestive_line <- -log10(1e-5)
genomewide_line <- -log10(5e-8)

png(filename = man_file, width = 2000, height = 1000, res = 150)

manhattan(
  assoc,
  chr = "chr",
  bp  = "pos",
  p   = "Score.pval",
  snp = "gene",           
  col = c("#4C72B0", "#E15759"),
  suggestiveline = suggestive_line,
  genomewideline = genomewide_line,
  annotatePval = 5e-8,     
  annotateTop  = FALSE,
  ylim = c(0, 10),         
  main = plot_label,       
  cex = 0.8
)

dev.off()
cat("Saved Manhattan plot:", man_file, "\n")

cat(">>> Finished phenotype:", trait, "\n")
