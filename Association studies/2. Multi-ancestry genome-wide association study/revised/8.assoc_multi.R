#!/usr/bin/env Rscript

.libPaths("/home1/tanxinli/R_lib")
# make sure glpk (igraph dependency) is found
Sys.setenv(
  LD_LIBRARY_PATH = paste(
    "/home1/tanxinli/R_lib",
    "/home1/tanxinli/R_lib/lib64",
    Sys.getenv("LD_LIBRARY_PATH"),
    sep = ":"
  )
)
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



args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]
cat("Running for trait:", trait, "\n")

# ============================================
# 1. Read trait name
# ============================================
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("No phenotype argument provided.")
trait <- args[1]

cat(">>> Running plots for phenotype:", trait, "\n")

# ============================================
# 2. Build file paths
# ============================================
base_dir <- "/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/result/null_model"


assoc_file <- file.path(base_dir, trait, "assoc_all_chr1_22.rds")
outdir     <- file.path(base_dir, trait, "plots")

if (!file.exists(assoc_file)) {
  stop("Assoc RDS file does not exist: ", assoc_file)
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cat("Output directory:", outdir, "\n")

# ============================================
# 3. Load assoc results
# ============================================
cat("Loading data...\n")
assoc <- readRDS(assoc_file)

# 确保列类型正确
assoc$chr <- as.numeric(as.character(assoc$chr))
assoc$pos <- as.numeric(assoc$pos)
assoc$Score.pval <- as.numeric(assoc$Score.pval)

# 移除 NA 值以防画图报错
assoc <- assoc %>% filter(!is.na(Score.pval) & !is.na(chr) & !is.na(pos))


if("gene_name" %in% colnames(assoc)){
  label_col <- "gene_name"
  cat("Using 'gene_name' for annotation.\n")
} else {
  label_col <- "variant.id"
  cat("No gene column found. Using 'variant.id' for annotation.\n")
}

# ============================================
# 4. QQ Plot (Your Custom Function)
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
  
  p <- ggplot(dat, aes(-log10(exp), -log10(obs))) +
    geom_line(aes(-log10(exp), -log10(upper)), color="gray") +
    geom_line(aes(-log10(exp), -log10(lower)), color="gray") +
    geom_point(size=1, alpha=0.5) +
    geom_abline(intercept=0, slope=1, color="red") +
    annotate("text", x = -Inf, y = Inf, label = lambda_text,
             hjust = -0.2, vjust = 2, size = 5, fontface="bold") +
    xlab(expression(paste(-log[10], "(expected P)"))) +
    ylab(expression(paste(-log[10], "(observed P)"))) +
    ggtitle(paste0("QQ plot: ", trait_label)) +
    theme_bw() +
    theme(panel.grid = element_blank())
    
  return(p)
}

cat("Generating QQ plot...\n")
qq_fig <- qqPlot_custom(assoc$Score.pval, trait_label = trait)
qq_file <- file.path(outdir, paste0("qqplot_", trait, ".png"))
ggsave(qq_file, qq_fig, width = 6, height = 6, dpi = 300)

# ============================================
# 5. Manhattan Plot with Annotation
# ============================================
cat("Generating Manhattan plot...\n")

man_file <- file.path(outdir, paste0("manhattan_", trait, ".png"))

# 设置阈值
  
suggestive_line <- -log10(1e-5)
genomewide_line <- -log10(5e-8)

png(filename = man_file, width = 2000, height = 1000, res = 150)

# 调用 manhattan

manhattan(
  assoc,
  chr = "chr", 
  bp = "pos", 
  p = "Score.pval", 
  snp = label_col, # 这里决定了图上显示 Gene 名字还是 rsID
  col = c("#4C72B0", "#E15759"),
  suggestiveline = suggestive_line,
  genomewideline = genomewide_line,
  annotatePval = 5e-8,  # 只有 P < 5e-8 的点才会被标记名字
  annotateTop = FALSE,   # FALSE = 标记所有过线的; TRUE = 每个染色体只标记最高的
  main = paste("Manhattan Plot:", trait),
  cex = 0.8
)

dev.off()
cat("Saved Manhattan plot:", man_file, "\n")

# ============================================
# 6. Collect plots to a common folder
# ============================================
collect_dir <- file.path(base_dir, "all_plots")
dir.create(collect_dir, showWarnings = FALSE, recursive = TRUE)

file.copy(qq_file, file.path(collect_dir, basename(qq_file)), overwrite = TRUE)
file.copy(man_file, file.path(collect_dir, basename(man_file)), overwrite = TRUE)
cat("Copied plots to:", collect_dir, "\n")

cat(">>> Finished phenotype:", trait, "\n")