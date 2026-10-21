## Step 6 - Manhattan and QQ plots.
## Adapted from ../../meta-analysis1020/manhattan/regen_manhattan_meta_fixed.R:
## one label per genome-wide-significant locus, taken from
## meta_fixed_significant_loci.tsv (gene, else rsID, else variant ID).


.libPaths("/home1/tanxinli/R_lib")
suppressMessages({
  library(ggplot2); library(dplyr); library(data.table); library(ggrepel)
})
setDTthreads(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
options(scipen = 999)
 
args   <- commandArgs(trailingOnly = TRUE)
OUTDIR <- if (length(args) > 0) args[1] else
  "/project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline/output"
OUT    <- file.path(OUTDIR, "manhattan_meta_fixed")
SIG    <- 5e-8
SUG    <- 1e-5
MIN_STUDIES <- 2          # keep in line with min_studies in step 5
TITLE  <- sprintf("Manhattan Plot - META_FIXED (7 studies; variants in >= %d studies)", MIN_STUDIES)
 
d <- fread(file.path(OUTDIR, "meta_fixed_sumstats_annotated.txt"),
           select = c("CHR", "POS", "ID", "EAF", "P", "N_STUDIES"), data.table = FALSE)
d$CHR <- suppressWarnings(as.integer(d$CHR))
d <- d %>% filter(N_STUDIES >= MIN_STUDIES,
                  !is.na(P), P > 0, P <= 1, !is.na(CHR), CHR %in% 1:22, !is.na(POS))
cat("Variants analysed:", format(nrow(d), big.mark = ","), "\n")
d$logP <- -log10(d$P)
 
# lambda from ALL filtered variants (before thinning the plot data below)
lambda_gc <- function(p) median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
lambda <- lambda_gc(d$P)
cat(sprintf("lambda_GC = %.4f\n", lambda))
 
# cumulative bp position for the x-axis
chr_len <- d %>% group_by(CHR) %>% summarize(chr_len = max(POS), .groups = "drop") %>% arrange(CHR)
chr_len$chr_start <- c(0, cumsum(as.numeric(chr_len$chr_len))[-nrow(chr_len)])
d <- d %>% left_join(chr_len[, c("CHR", "chr_start")], by = "CHR") %>% mutate(BP_cum = POS + chr_start)
axis_df <- d %>% group_by(CHR) %>% summarize(center = mean(BP_cum), .groups = "drop")
 
d$base_col <- ifelse(d$CHR %% 2 == 0, "skyblue", "grey60")
sig_line <- -log10(SIG); sug_line <- -log10(SUG)
d$is_sig <- d$logP >= sig_line
 
# For plotting speed: keep every variant with P < 0.01 and a random 10% of the
# rest. The plot looks the same; lambda above already used all variants.
set.seed(1)
d_plot <- d %>% filter(P < 0.01 | runif(n()) < 0.10)
cat("Variants plotted:", format(nrow(d_plot), big.mark = ","), "\n")
 
loci_file <- file.path(OUTDIR, "meta_fixed_significant_loci.tsv")
labels_df <- NULL
if (file.exists(loci_file)) {
  labels_df <- fread(loci_file, data.table = FALSE) %>%
    mutate(CHR = as.integer(CHR)) %>%
    filter(CHR %in% 1:22) %>%
    left_join(chr_len[, c("CHR", "chr_start")], by = "CHR") %>%
    mutate(BP_cum = LEAD_POS + chr_start, logP = -log10(P))
}
 
p <- ggplot(d_plot, aes(x = BP_cum, y = logP)) +
  geom_point(aes(color = base_col), size = 1.1, alpha = 0.7) +
  geom_point(data = d_plot %>% filter(is_sig), color = "orange", size = 1.8, alpha = 0.9) +
  scale_color_identity() +
  geom_hline(yintercept = sig_line, linetype = "dashed", linewidth = 0.6, color = "orange") +
  geom_hline(yintercept = sug_line, linetype = "dotted", linewidth = 0.5) +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$CHR, expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  labs(x = "Chromosome", y = expression(-log[10](italic(P))), title = TITLE) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(size = 10))
 
if (!is.null(labels_df) && nrow(labels_df) > 0) {
  p <- p + geom_text_repel(
    data = labels_df, aes(x = BP_cum, y = logP, label = LABEL),
    inherit.aes = FALSE, size = 3.2, fontface = "italic",
    max.overlaps = Inf, box.padding = 0.6, point.padding = 0.3,
    segment.size = 0.4, min.segment.length = 0, seed = 42)
}
 
ggsave(paste0(OUT, ".png"),  p, width = 14, height = 7, dpi = 300)
ggsave(paste0(OUT, ".tiff"), p, width = 14, height = 7, dpi = 300, compression = "lzw")
cat("Saved", paste0(OUT, ".png"), "and .tiff\n")
 
# ------------------------------------------------------------------ #
# QQ plot (uses all filtered variants, not the thinned plot data)
# ------------------------------------------------------------------ #
obs <- sort(d$logP, decreasing = TRUE)
qq  <- data.frame(exp = -log10(ppoints(length(obs))), obs = obs)
# thin the dense null region for plotting speed
qq  <- qq[qq$exp > 2 | seq_len(nrow(qq)) %% 50 == 0, ]
 
q <- ggplot(qq, aes(exp, obs)) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  geom_point(size = 0.8) +
  labs(x = expression(Expected ~ -log[10](italic(P))), y = expression(Observed ~ -log[10](italic(P))),
       title = sprintf("QQ Plot - META_FIXED (lambda = %.3f)", lambda)) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(OUTDIR, "qq_meta_fixed.png"), q, width = 6, height = 6, dpi = 300)
cat("Saved", file.path(OUTDIR, "qq_meta_fixed.png"), "\n")
# ------------------------------------------------------------------ #
# Inflation checks: lambda without MHC / known loci, and by MAF bin
# ------------------------------------------------------------------ #
# extended MHC, hg38 chr6:25-34 Mb
in_mhc <- d$CHR == 6 & d$POS >= 25e6 & d$POS <= 34e6
# GWS loci from step 5, padded by 500 kb on each side
in_loci <- rep(FALSE, nrow(d))
if (!is.null(labels_df) && nrow(labels_df) > 0) {
  for (i in seq_len(nrow(labels_df))) {
    in_loci <- in_loci | (d$CHR == labels_df$CHR[i] &
                          d$POS >= labels_df$START[i] - 5e5 & d$POS <= labels_df$END[i] + 5e5)
  }
}
d$MAF <- pmin(d$EAF, 1 - d$EAF)
d$MAF_BIN <- cut(d$MAF, breaks = c(-Inf, 0.01, 0.05, 0.20, Inf), right = FALSE,
                 labels = c("MAF < 1%", "1% <= MAF < 5%", "5% <= MAF < 20%", "MAF >= 20%"))

lam_tab <- rbind(
  data.frame(SET = "All variants",                   N = nrow(d),             LAMBDA = lambda),
  data.frame(SET = "Excluding MHC (chr6:25-34Mb)",   N = sum(!in_mhc),        LAMBDA = lambda_gc(d$P[!in_mhc])),
  data.frame(SET = "Excluding GWS loci (+/-500kb)",  N = sum(!in_loci),       LAMBDA = lambda_gc(d$P[!in_loci])),
  data.frame(SET = "Excluding MHC and GWS loci",     N = sum(!in_mhc & !in_loci),
             LAMBDA = lambda_gc(d$P[!in_mhc & !in_loci])),
  do.call(rbind, lapply(levels(d$MAF_BIN), function(b) {
    p <- d$P[d$MAF_BIN %in% b]
    data.frame(SET = b, N = length(p), LAMBDA = if (length(p)) lambda_gc(p) else NA_real_)
  })))
lam_tab$LAMBDA <- round(lam_tab$LAMBDA, 4)
cat("\nGenomic inflation (lambda_GC):\n"); print(lam_tab, row.names = FALSE)
fwrite(lam_tab, file.path(OUTDIR, "lambda_meta_fixed.tsv"), sep = "\t")

# QQ plot stratified by MAF bin (expected P computed within each bin)
qq_maf <- d %>% filter(!is.na(MAF_BIN)) %>% group_by(MAF_BIN) %>%
  arrange(desc(logP), .by_group = TRUE) %>%
  mutate(exp = -log10(ppoints(n())), obs = logP) %>%
  filter(exp > 2 | row_number() %% 50 == 0) %>% ungroup()
bin_lab <- lam_tab %>% filter(SET %in% levels(d$MAF_BIN))
bin_lab <- setNames(sprintf("%s (n = %s, lambda = %.3f)", bin_lab$SET,
                            format(bin_lab$N, big.mark = ","), bin_lab$LAMBDA), bin_lab$SET)

q2 <- ggplot(qq_maf, aes(exp, obs, colour = MAF_BIN)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey30", linetype = "dashed") +
  geom_point(size = 0.8) +
  scale_colour_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"),
                      labels = bin_lab, name = NULL) +
  labs(x = expression(Expected ~ -log[10](italic(P))), y = expression(Observed ~ -log[10](italic(P))),
       title = "QQ Plot by MAF - META_FIXED",
       subtitle = sprintf("lambda: all = %.3f; excl. MHC + GWS loci = %.3f",
                          lambda, lam_tab$LAMBDA[4])) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.text = element_text(size = 9), legend.background = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave(file.path(OUTDIR, "qq_meta_fixed_by_maf.png"), q2, width = 7, height = 7, dpi = 300)
cat("Saved", file.path(OUTDIR, "qq_meta_fixed_by_maf.png"), "\n")
