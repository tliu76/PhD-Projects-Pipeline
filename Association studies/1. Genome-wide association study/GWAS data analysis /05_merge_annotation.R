################################################################################
## 05_merge_annotation.R
##
## Merge the ANNOVAR annotation back onto the meta-analysis summary
##           statistics, report genomic inflation (lambda_GC), and list
##           genome-wide significant and suggestive variants.
##
## Input   : 1. output/meta_fixed_sumstats.txt                     (step 3)
##           2. output/meta_fixed_annotated_hg38.hg38_multianno.csv (step 4)
##
## Output  : 1. output/meta_fixed_sumstats_annotated.txt  - all variants + annotation
##           2. output/meta_fixed_significant.txt         - P < 5e-8
##           3. output/meta_fixed_suggestive.txt          - P < 1e-5
##           4. output/meta_fixed_lead_variants.txt       - one lead variant per locus
################################################################################
.libPaths("/home1/tanxinli/R_lib")
library(data.table)
library(tidyverse)
library(dplyr)


#### 1. Input and output paths ------------------------------------------------

pipe_dir <- "/project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline"
out_dir  <- file.path(pipe_dir, "output")


sumstats_file <- file.path(out_dir, "meta_fixed_sumstats.txt")
annovar_file  <- file.path(out_dir, "meta_fixed_annotated_hg38.hg38_multianno.csv")

annotated_file   <- file.path(out_dir, "meta_fixed_sumstats_annotated.txt")
significant_file <- file.path(out_dir, "meta_fixed_significant.txt")
suggestive_file  <- file.path(out_dir, "meta_fixed_suggestive.txt")
lead_file        <- file.path(out_dir, "meta_fixed_lead_variants.txt")
loci_file       <- file.path(out_dir, "meta_fixed_significant_loci.tsv")


# Thresholds
p_gws       <- 5e-8      # genome-wide significance
p_suggest   <- 1e-5      # suggestive significance
min_studies <- 2         # minimum number of contributing studies
locus_gap_bp <- 500000   # a gap larger than this between significant variants starts a new locus

#### 2. Read the summary statistics -------------------------------------------

cat("Reading", sumstats_file, "\n")
meta <- fread(sumstats_file)
cat("  Variants:", nrow(meta), "\n")
# Variants: 9867245 

#### 3. Read and clean the ANNOVAR output --------------------------------------

cat("Reading", annovar_file, "\n")
anno <- fread(annovar_file, sep = ",", header = FALSE, skip = 1,
              colClasses = "character", fill = TRUE)
# Take the column names from the real header line; the last data column holds
# all the -otherinfo fields together.
header <- strsplit(readLines(annovar_file, n = 1), ",")[[1]]
setnames(anno, c(header[1:(ncol(anno) - 1)], "Otherinfo"))
anno <- anno[Chr != "CHR"]
anno[, c("ID", "rsID") := tstrsplit(Otherinfo, "\t", fixed = TRUE)]
head(anno[, .(ID, rsID, Gene.refGene)])
anno$Otherinfo <- NULL

fwrite(anno, file.path(out_dir, "anno.csv"))
# Keep only the annotation columns we need
anno <- anno[, .(ID,rsID,
                 FUNC        = Func.refGene,
                 GENE        = Gene.refGene,
                 GENE_DETAIL = GeneDetail.refGene,
                 EXONIC      = ExonicFunc.refGene,
                 AACHANGE    = AAChange.refGene,
                 CYTOBAND    = cytoBand)]

#### 4. Merge annotation onto the summary statistics ---------------------------
meta <- meta %>% rename("ID"="SNPID")
meta <- merge(meta, anno, by = c("ID","rsID"), all.x = TRUE, sort = FALSE)


cat("  With ANNOVAR record:", format(meta[!is.na(FUNC), .N], big.mark = ","),
    "| with rsID:", format(meta[!is.na(rsID), .N], big.mark = ","), "\n")
if (meta[is.na(FUNC), .N] > 0) {
  cat("  Examples without annotation:\n")
  print(head(meta[is.na(FUNC), .(ID, CHR, POS)]))
}

# Total: 9,867,245    With ANNOVAR record: 9,768,380 | with rsID: 9,656,249 

# Sort by chromosome (numerically) and position; key columns first
meta[, CHR_NUM := suppressWarnings(as.integer(CHR))]
setorder(meta, CHR_NUM, POS, na.last = TRUE)
meta[, CHR_NUM := NULL]
setcolorder(meta, c("CHR", "POS", "ID", "rsID", "REF", "ALT", "EA", "NEA"))

fwrite(meta, annotated_file, sep = "\t", quote = FALSE, na = "NA")
cat("Wrote", annotated_file, "\n")
# Wrote /project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline/output/meta_fixed_sumstats_annotated.txt

#### 5. Genomic inflation factor (lambda_GC) -----------------------------------
## lambda = median observed chi-square / median expected chi-square (0.4549)

lambda_gc <- function(p) {
  p <- p[!is.na(p) & p > 0 & p <= 1]
  median(qchisq(p, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1)
}

cat(sprintf("\nlambda_GC, all variants              : %.4f\n", lambda_gc(meta$P)))
cat(sprintf("lambda_GC, variants in >= %d studies : %.4f\n",
                min_studies, lambda_gc(meta[N_STUDIES >= min_studies, P])))

# lambda_GC,variants in >= 2 studies             : 1.1299
# lambda_GC,variants in >= 3 studies             :  1.1335
# lambda_GC,variants in >= 4 studies             : 1.1396
# lambda_GC,variants in >= 5 studies             : 1.1432
# lambda_GC,variants in >= 6 studies             :  1.1474
# lambda_GC,variants in >= 7 studies             : 1.1557



#### 6. Suggestive variants ---------------------------------------------------

report_cols <- c("CHR", "POS", "ID", "rsID", "EA", "NEA", "EAF", "BETA", "SE", "P",
                 "N", "N_STUDIES", "Direction", "HetISq", "HetPVal",
                 "FUNC", "GENE", "CYTOBAND")

suggestive <- meta[N_STUDIES >= min_studies & P < p_suggest, ..report_cols]
suggestive[, OR     := exp(BETA)]
suggestive[, OR_L95 := exp(BETA - 1.96 * SE)]
suggestive[, OR_U95 := exp(BETA + 1.96 * SE)]
setorder(suggestive, P)

fwrite(suggestive, suggestive_file, sep = "\t", quote = FALSE, na = "NA")
cat("\nSuggestive variants (P <", p_suggest, ", >=", min_studies, "studies):",
    nrow(suggestive), "\n")
cat("Wrote", suggestive_file, "\n")

# Suggestive variants (P < 1e-05 , >= 2 studies): 2419 
#Wrote /project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline/output/meta_fixed_suggestive.txt 


#### 7. Genome-wide significant loci ------------------------------------------
## Significant variants are sorted by position; a new locus starts on a new
## chromosome or when the gap to the previous significant variant exceeds
## locus_gap_bp. Each locus is summarised by its most significant variant.

# First gene name, e.g. "FAM138A,OR4F5" -> "FAM138A"; used as the plot label
first_gene <- function(g) {
  g <- sub("[,;].*$", "", g)
  fifelse(is.na(g) | g == "", NA_character_, g)
}

sig <- meta[N_STUDIES >= min_studies & !is.na(P) & P < p_gws]

if (nrow(sig) == 0) {
  cat("\nNo variants with P <", p_gws, "in >=", min_studies, "studies\n")
  # remove an old loci file so step 6 does not plot stale labels
  if (file.exists(loci_file)) file.remove(loci_file)
} else {
  setorder(sig, CHR, POS)
  sig[, LOCUS := cumsum(c(TRUE, CHR[-1] != CHR[-.N] | diff(POS) > locus_gap_bp))]
  
  loci <- sig[order(P), {
    top   <- .SD[1]                                   # most significant variant
    genes <- first_gene(GENE)
    genes <- genes[!is.na(genes)]
    label <- if (length(genes) > 0) genes[1]
    else if (!is.na(top$rsID)) top$rsID
    else top$ID
    .(CHR       = top$CHR,
      LEAD_POS  = top$POS,
      LEAD_ID   = top$ID,
      LEAD_rsID = top$rsID,
      START     = min(POS),
      END       = max(POS),
      N_SIG     = .N,
      EA        = top$EA,
      NEA       = top$NEA,
      EAF       = top$EAF,
      BETA      = top$BETA,
      SE        = top$SE,
      OR        = exp(top$BETA),
      OR_L95    = exp(top$BETA - 1.96 * top$SE),
      OR_U95    = exp(top$BETA + 1.96 * top$SE),
      P         = top$P,
      N         = top$N,
      N_STUDIES = top$N_STUDIES,
      Direction = top$Direction,
      HetISq    = top$HetISq,
      HetPVal   = top$HetPVal,
      LEAD_GENE = top$GENE,
      LEAD_FUNC = top$FUNC,
      CYTOBAND  = top$CYTOBAND,
      LOCUS_GENES = paste(unique(unlist(strsplit(GENE[!is.na(GENE)], "[,;]"))), collapse = ","),
      LABEL     = label)
  }, by = LOCUS]
  
  setorder(loci, P)
  fwrite(loci, loci_file, sep = "\t", quote = FALSE, na = "NA")
  
  cat("\nGenome-wide significant loci (P <", p_gws, ", >=", min_studies, "studies):",
      nrow(loci), "\n")
  print(loci[, .(CHR, LEAD_POS, LEAD_rsID, P, OR, N_STUDIES, Direction, LABEL)], row.names = FALSE)
  cat("Wrote", loci_file, "\n")
}

