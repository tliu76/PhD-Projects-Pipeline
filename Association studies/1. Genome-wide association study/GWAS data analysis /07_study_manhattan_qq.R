## Step 7 - Manhattan and QQ plots for each individual study (the 7 METAL inputs
## in 02_metal.txt). 
##
## No variant filter beyond valid P and CHR 1-22 (each study is plotted as it
## enters METAL). lambda_1000 rescales lambda to 1000 cases / 1000 controls.


.libPaths("/home1/tanxinli/R_lib")
suppressMessages({
  library(ggplot2); library(dplyr); library(data.table); library(ggrepel)
})
setDTthreads(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
options(scipen = 999)

args   <- commandArgs(trailingOnly = TRUE)
OUTDIR <- if (length(args) > 0) args[1] else
  "/project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline/output"
PLOTDIR <- file.path(OUTDIR, "study_plots")
dir.create(PLOTDIR, showWarnings = FALSE, recursive = TRUE)
SIG <- 5e-8
SUG <- 1e-5
LOCUS_GAP_BP <- 500000    # same as locus_gap_bp in step 5
ANNO_FILE <- file.path(OUTDIR, "meta_fixed_sumstats_annotated.txt")

# per-study palette (different from the meta-analysis plots)
COL_ODD  <- "#00695C"   # dark teal
COL_EVEN <- "#80CBC4"   # light teal
COL_SIG  <- "#C2185B"   # magenta-red for GWS hits / line
COL_QQ   <- "#00695C"
COL_DIAG <- "grey30"

# same files and order as 02_metal.txt; case/control counts as in
# iterations/scripts/08_study_lambda_diagnostics.R
IN <- "/project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline/input"
studies <- rbindlist(list(
  list("CCRLP_Kaiser", "/project2/desmith_128/datasets/CCRLP/CCRLP_Kaiser_TOPMed_imputation/LAT/binary/remove_outlier/ck_0708.otherallele.txt", 1930L, 8481L),
  list("CCLS", "/project2/desmith_128/datasets/CCLS/LAT/Imputed/binary/merged/remove_outlier/ccls_0708.otherallele.txt", 602L, 510L),
  list("Texas_SOL", "/project2/desmith_128/datasets/Texas_SOL/Texas/Results/binary/remove_dup/merged2/remove_outlier/texas_sol_0708.otherallele.txt", 499L, 2622L),
  list("St_Jude_MESA", file.path(IN, "St_Jude_MESA_Rsq_filtered.otherallele.txt"), 782L, 1146L),
  list("Guatemala", "/project2/desmith_128/tanxin/scratch1_new/GUAT/revise_ref_impute/merge/guat_0708.otherallele.txt", 38L, 49L),
  list("CHLA_CALSEC_MEC", file.path(IN, "CHLA_CALSEC_MEC_20PCs.otherallele.txt"), 742L, 5007L),
  list("Mexican", file.path(IN, "Mexican_20PCs.otherallele.txt"), 807L, 5566L)
))
setnames(studies, c("STUDY", "GWAS", "NCASE", "NCTRL"))

idx <- if (length(args) > 1) as.integer(args[2]) else seq_len(nrow(studies))
stopifnot(all(idx %in% seq_len(nrow(studies))))

lambda_gc   <- function(p) median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
lambda_1000 <- function(l, ncase, nctrl) 1 + (l - 1) * (1 / ncase + 1 / nctrl) / (1 / 1000 + 1 / 1000)

# first gene name, e.g. "FAM138A,OR4F5" -> "FAM138A" (as in step 5)
first_gene <- function(g) {
  g <- sub("[,;].*$", "", g)
  fifelse(is.na(g) | g == "" | g == ".", NA_character_, g)
}

# annotation (ID -> rsID, GENE) from step 5, read once for all studies
cat("Reading annotation from", ANNO_FILE, "\n")
anno <- fread(ANNO_FILE, select = c("ID", "rsID", "GENE"), data.table = TRUE)
setkey(anno, ID)
anno <- unique(anno, by = "ID")

# GWS loci of one study (d: CHR, POS, ID, P), labelled by the lead variant's gene
study_loci <- function(d) {
  sig <- d[P < SIG, .(CHR, POS, ID, P)]
  if (nrow(sig) == 0) return(NULL)
  setorder(sig, CHR, POS)
  sig[, LOCUS := cumsum(c(TRUE, CHR[-1] != CHR[-.N] | diff(POS) > LOCUS_GAP_BP))]
  sig <- anno[sig, on = "ID"]
  loci <- sig[order(P), {
    top <- .SD[1]
    g   <- first_gene(top$GENE)
    .(CHR = top$CHR, LEAD_POS = top$POS, LEAD_ID = top$ID, LEAD_rsID = top$rsID,
      START = min(POS), END = max(POS), N_SIG = .N, P = top$P, LEAD_GENE = top$GENE,
      LABEL = if (!is.na(g)) g else if (!is.na(top$rsID)) top$rsID else top$ID)
  }, by = LOCUS]
  setorder(loci, P)
  loci
}

plot_study <- function(s) {
  cat(sprintf("\n==== %s ====\n%s\n", s$STUDY, s$GWAS))
  d <- fread(s$GWAS, select = c("#CHROM", "POS", "ID", "P"), data.table = TRUE)
  setnames(d, "#CHROM", "CHR")
  d[, CHR := suppressWarnings(as.integer(sub("^chr", "", CHR)))]
  d <- d[!is.na(P) & P > 0 & P <= 1 & !is.na(CHR) & CHR %in% 1:22 & !is.na(POS)]
  cat("Variants analysed:", format(nrow(d), big.mark = ","), "\n")
  d[, logP := -log10(P)]

  lam  <- lambda_gc(d$P)
  lam1 <- lambda_1000(lam, s$NCASE, s$NCTRL)
  n_sig <- sum(d$P < SIG); n_sug <- sum(d$P < SUG)
  cat(sprintf("lambda_GC = %.4f | lambda_1000 = %.4f | P < 5e-8: %d | P < 1e-5: %d\n",
              lam, lam1, n_sig, n_sug))

  # cumulative bp position for the x-axis
  chr_len <- d[, .(chr_len = max(POS)), by = CHR][order(CHR)]
  chr_len[, chr_start := c(0, cumsum(as.numeric(chr_len))[-.N])]
  d <- chr_len[, .(CHR, chr_start)][d, on = "CHR"]
  d[, BP_cum := POS + chr_start]
  axis_df <- d[, .(center = mean(BP_cum)), by = CHR][order(CHR)]

  loci <- study_loci(d)
  out_l <- file.path(PLOTDIR, sprintf("loci_%s.tsv", s$STUDY))
  if (!is.null(loci)) {
    loci <- chr_len[, .(CHR, chr_start)][loci, on = "CHR"]
    loci[, `:=`(BP_cum = LEAD_POS + chr_start, logP = -log10(P))]
    cat("GWS loci:", nrow(loci), "\n")
    print(loci[, .(CHR, LEAD_POS, LEAD_rsID, P, N_SIG, LABEL)], row.names = FALSE)
    fwrite(loci[, !c("chr_start", "BP_cum", "logP")], out_l, sep = "\t", quote = FALSE, na = "NA")
  } else {
    cat("GWS loci: 0\n")
    if (file.exists(out_l)) file.remove(out_l)
  }

  d[, base_col := ifelse(CHR %% 2 == 0, COL_EVEN, COL_ODD)]
  sig_line <- -log10(SIG); sug_line <- -log10(SUG)

  # thin for plotting speed: all P < 0.01 plus a random 10% of the rest
  set.seed(1)
  d_plot <- d[P < 0.01 | runif(.N) < 0.10]

  subtitle <- sprintf("%s cases / %s controls | %s variants | lambda = %.3f (lambda_1000 = %.3f)",
                      format(s$NCASE, big.mark = ","), format(s$NCTRL, big.mark = ","),
                      format(nrow(d), big.mark = ","), lam, lam1)

  p <- ggplot(d_plot, aes(x = BP_cum, y = logP)) +
    geom_point(aes(color = base_col), size = 1.1, alpha = 0.7) +
    geom_point(data = d_plot[logP >= sig_line], color = COL_SIG, size = 1.8, alpha = 0.9) +
    scale_color_identity() +
    geom_hline(yintercept = sig_line, linetype = "dashed", linewidth = 0.6, color = COL_SIG) +
    geom_hline(yintercept = sug_line, linetype = "dotted", linewidth = 0.5) +
    scale_x_continuous(breaks = axis_df$center, labels = axis_df$CHR, expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0.02, 0.02), limits = c(0, max(sig_line + 1, max(d$logP)))) +
    labs(x = "Chromosome", y = expression(-log[10](italic(P))),
         title = sprintf("Manhattan Plot - %s", s$STUDY), subtitle = subtitle) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none",
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          axis.text.x = element_text(size = 10))
  if (!is.null(loci)) {
    p <- p + geom_text_repel(
      data = loci, aes(x = BP_cum, y = logP, label = LABEL),
      inherit.aes = FALSE, size = 3.2, fontface = "italic",
      max.overlaps = Inf, box.padding = 0.6, point.padding = 0.3,
      segment.size = 0.4, min.segment.length = 0, seed = 42)
  }
  out_m <- file.path(PLOTDIR, sprintf("manhattan_%s.png", s$STUDY))
  ggsave(out_m, p, width = 14, height = 7, dpi = 300)
  cat("Saved", out_m, "\n")

  # QQ plot (all variants, thinned in the dense null region)
  obs <- sort(d$logP, decreasing = TRUE)
  qq  <- data.frame(exp = -log10(ppoints(length(obs))), obs = obs)
  qq  <- qq[qq$exp > 2 | seq_len(nrow(qq)) %% 50 == 0, ]
  q <- ggplot(qq, aes(exp, obs)) +
    geom_abline(slope = 1, intercept = 0, colour = COL_DIAG, linetype = "dashed") +
    geom_point(size = 0.8, colour = COL_QQ) +
    labs(x = expression(Expected ~ -log[10](italic(P))), y = expression(Observed ~ -log[10](italic(P))),
         title = sprintf("QQ Plot - %s", s$STUDY),
         subtitle = sprintf("lambda = %.3f; lambda_1000 = %.3f", lam, lam1)) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11))
  out_q <- file.path(PLOTDIR, sprintf("qq_%s.png", s$STUDY))
  ggsave(out_q, q, width = 6, height = 6, dpi = 300)
  cat("Saved", out_q, "\n")

  data.table(STUDY = s$STUDY, NCASE = s$NCASE, NCTRL = s$NCTRL, N_VARIANTS = nrow(d),
             LAMBDA = round(lam, 4), LAMBDA_1000 = round(lam1, 4),
             N_P_5e8 = n_sig, N_P_1e5 = n_sug, N_LOCI = if (is.null(loci)) 0L else nrow(loci))
}

for (i in idx) {
  res <- plot_study(studies[i])
  # one summary row per study, so array tasks don't overwrite each other
  fwrite(res, file.path(PLOTDIR, sprintf("lambda_%s.tsv", studies$STUDY[i])), sep = "\t")
  gc()
}
