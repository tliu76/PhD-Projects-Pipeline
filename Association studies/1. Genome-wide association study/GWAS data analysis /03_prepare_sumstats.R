## Step 3 - Format the METAL output, add dbSNP151 rsIDs, write the ANNOVAR input.
##
## Usage: Rscript 03_prepare_sumstats.R <metal.tbl> <dbsnp.txt> <out_dir>
.libPaths("/home1/tanxinli/R_lib")
suppressMessages(library(data.table))



#### 1. Input and output paths ------------------------------------------------

pipe_dir   <- "/project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline"
dbsnp_file <- "/project2/desmith_488/tanxinli/gwaslab/1kg_dbsnp151_hg38_auto.txt"


out_dir        <- file.path(pipe_dir, "output")
sumstats_file  <- file.path(out_dir, "meta_fixed_sumstats.txt")
avinput_file   <- file.path(out_dir, "meta_fixed.avinput")

#### 2. Read and clean the METAL output ---------------------------------------

meta <- fread("/project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline/output/meta-analysis_nogc_20260928_20PCs1.tbl")

# Rename METAL columns to standard names
setnames(meta,
         old = c("MarkerName", "Allele1", "Allele2", "Freq1",  "FreqSE" , "MinFreq" , "MaxFreq",    "Effect", "StdErr", "P-value", "Direction", "HetISq", "HetChiSq",
                 "HetDf", "HetPVal",  "TotalSampleSize"),
         new = c("SNPID",         "EA",      "NEA",     "EAF",   "FreqSE" , "MinFreq" , "MaxFreq",     "BETA",   "SE",     "P", "Direction", "HetISq", "HetChiSq",
                 "HetDf", "HetPVal",      "N"),
         skip_absent = TRUE)

### 11094446 observations ###

# METAL writes alleles in lower case
meta[, EA  := toupper(EA)]
meta[, NEA := toupper(NEA)]

head(meta)
# SNPID     EA    NEA    EAF FreqSE MinFreq MaxFreq    BETA     SE        P Direction HetISq HetChiSq HetDf HetPVal     N
# <char> <char> <char>  <num>  <num>   <num>   <num>   <num>  <num>    <num>    <char>  <num>    <num> <int>   <num> <int>
#   1:    chr11:12541586:A:G      A      G 0.4082 0.0085  0.2948  0.4186 -0.0278 0.0239 0.244000   -+---+-    0.0    1.158     6  0.9789 28748
# 2:     chr4:70553471:A:G      A      G 0.2729 0.0057  0.2467  0.2799  0.0025 0.0266 0.925500   +-+--+-   21.5    7.642     6  0.2655 28724
# 3:     chr3:97729697:A:G      A      G 0.0147 0.0020  0.0109  0.0181 -0.1005 0.1217 0.408900   --+??+-    0.0    1.176     4  0.8821 26763
# 4:     chr1:83758028:C:T      T      C 0.1629 0.0123  0.1410  0.1711 -0.0050 0.0438 0.910100   +-?+???   33.2    2.995     2  0.2237 13451
# 5: chr6:159697350:T:TACC      T   TACC 0.1722 0.0039  0.1649  0.1746  0.1272 0.0426 0.002814   ++?+???    0.0    0.262     2  0.8774 13451
# 6:     chr7:47120610:G:C      C      G 0.0656 0.0092  0.0492  0.0716  0.0179 0.0683 0.793100   -+?-???   23.4    2.610     2  0.2712 13451

# Marker IDs look like chr1:730090:G:C  ->  CHR, POS, REF, ALT
# 3 Texas_SOL-only markers carry two IDs joined by ";" (e.g. chr3:25134959:G:A;chr3:25134959:G:C); drop them
cat("  Dropping markers with a multi-allelic ';' ID:", meta[grepl(";", SNPID, fixed = TRUE), .N], "\n")
meta <- meta[!grepl(";", SNPID, fixed = TRUE)]
meta[, c("CHR", "POS", "REF", "ALT") := tstrsplit(SNPID, ":", fixed = TRUE)]
meta[, CHR := sub("^chr", "", CHR)]
meta[, POS := as.integer(POS)]


meta[, N_STUDIES := nchar(gsub("?", "", Direction, fixed = TRUE))]
meta[, table(N_STUDIES == HetDf + 1)]

############ check the number of variants in each study #########
study_tab <- meta[, .N, keyby = N_STUDIES]
study_tab[, Percent := round(100 * N / sum(N), 1)]
study_tab[, Cum_at_least := rev(cumsum(rev(N)))]   # 至少在 k 个研究中出现的位点数
study_tab
# Key: <N_STUDIES>
#   N_STUDIES       N Percent Cum_at_least
# <int>   <int>   <num>        <int>
#   1:         1 1227201    11.1     11094446
# 2:         2  507423     4.6      9867245
# 3:         3  826900     7.5      9359822
# 4:         4  288765     2.6      8532922
# 5:         5  404035     3.6      8244157
# 6:         6  935780     8.4      7840122
# 7:         7 6904342    62.2      6904342

library(ggplot2)
ggplot(study_tab, aes(x = factor(N_STUDIES), y = N)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = paste0(format(N, big.mark = ","), "\n(", Percent, "%)")),
            vjust = -0.3, size = 3) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Number of contributing studies", y = "Number of variants") +
  theme_classic()

### 每个研究贡献多少位点 ###
study_names <- c("CCRLP_Kaiser","CCLS","Texas_SOL","St_Jude_MESA","GUAT","CHLA_CALSEC_MEC","Mexican")  # 换成你 PROCESS 顺序里的研究名，第 7 个是 Mexican

presence <- meta[, lapply(1:7, function(i) substr(Direction, i, i) != "?")]
setnames(presence, study_names)

colSums(presence)   # 每个研究有多少位点

library(UpSetR)
upset(as.data.frame(presence * 1L), nsets = 7, nintersects = 20,
      order.by = "freq", mainbar.y.label = "Number of variants")
#################################
### Filter out variants in only one studies  N_STUDIES=1
library(dplyr)
library(tidyverse)
meta <- meta %>% filter(N_STUDIES>1)   # remove 1227201 variants:  11094446-1227201 = 9867245

#################################

# Sanity check: the effect/other alleles should be the REF/ALT pair in the ID
n_mismatch <- meta[!((EA == REF & NEA == ALT) | (EA == ALT & NEA == REF)), .N]
cat("  Markers whose EA/NEA do not match the ID's REF/ALT:", n_mismatch, "\n")
# Markers whose EA/NEA do not match the ID's REF/ALT: 0 

cat("  Number of contributing studies per marker:\n")
print(meta[, .N, keyby = N_STUDIES])
# Key: <N_STUDIES>
#   N_STUDIES       N
# <int>   <int>
#   1:         2  507423
# 2:         3  826900
# 3:         4  288765
# 4:         5  404035
# 5:         6  935780
# 6:         7 6904342


#### 3. Add rsIDs from dbSNP151 -----------------------------------------------
dbsnp_file <- "/project2/desmith_488/tanxinli/gwaslab/1kg_dbsnp151_hg38_auto.txt"
# SNPID   rsID    CHR     POS     NEA     EA
# 1:10390:CCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA:C  rs1557426845    1       10390   CCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA    C

### dbsnp  66,989,071 snps  5 variables
dbsnp <- fread(dbsnp_file,
               select     = c("rsID", "CHR", "POS", "NEA", "EA"),
               colClasses = list(character = c("rsID", "CHR", "NEA", "EA"),
                                 integer   = "POS"))


# Keep only dbSNP positions that appear in the meta-analysis (saves memory)
dbsnp <- dbsnp[unique(meta[, .(CHR, POS)]), on = c("CHR", "POS"), nomatch = NULL]

### 9,780,377 dbsnps ###


# Match on CHR, POS and the two alleles regardless of which one is REF.
# Sorting each allele pair (A_lo, A_hi) makes "G/C" and "C/G" identical.
dbsnp[, A_lo := pmin(NEA, EA)]
dbsnp[, A_hi := pmax(NEA, EA)]
meta[,  A_lo := pmin(REF, ALT)]
meta[,  A_hi := pmax(REF, ALT)]


# A few variants have more than one rsID; keep the first one
n_multi <- dbsnp[, .N, by = .(CHR, POS, A_lo, A_hi)][N > 1, .N]
dbsnp   <- unique(dbsnp, by = c("CHR", "POS", "A_lo", "A_hi"))
cat("  Variants with more than one rsID (first kept):", n_multi, "\n")
# Variants with more than one rsID (first kept): 0 

meta[dbsnp, rsID := i.rsID, on = c("CHR", "POS", "A_lo", "A_hi")]
meta[, c("A_lo", "A_hi") := NULL]

n_rsid <- meta[!is.na(rsID), .N]
cat(sprintf("  Markers with an rsID: %d (%.1f%%)\n", n_rsid, 100 * n_rsid / nrow(meta)))
### Of 9867245 variants,  Markers with an rsID: 9656249 (97.9%) ### 

#### 4. Write the formatted summary statistics --------------------------------

# Sort by chromosome (numerically) and position
meta[, CHR_NUM := suppressWarnings(as.integer(CHR))]
setorder(meta, CHR_NUM, POS, na.last = TRUE)
meta[, CHR_NUM := NULL]

# Put the key columns first
first_cols <- c("SNPID","CHR", "POS",  "rsID", "REF", "ALT", "EA", "NEA",  "EAF" ,"FreqSE", 
"MinFreq", "MaxFreq",    "BETA"  ,   "SE"  , "P", "Direction", "HetISq", "HetChiSq", "HetDf",  "HetPVal", "N", "N_STUDIES")
setcolorder(meta, first_cols)

fwrite(meta, sumstats_file, sep = "\t", quote = FALSE, na = "NA")
cat("Wrote", sumstats_file, "\n")
# Wrote /project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline/output/meta_fixed_sumstats.txt 



#### 5. Write the ANNOVAR input -----------------------------------------------
##
## ANNOVAR needs: CHR, START, END, REF, ALT (+ extra columns kept via -otherinfo).
## VCF-style indels share a leading base (e.g. REF=AT, ALT=A); ANNOVAR style
## drops it, as convert2annovar.pl does:
##   deletion    AT -> A   becomes  T -> -   and starts one base later
##   insertion   A  -> AT  becomes  - -> T   and keeps the same position
## The marker ID goes in column 6 so step 5 can merge the annotation back on ID.

av <- meta[, .(CHR, START = POS, REF, ALT, SNPID,rsID)] 
### 9867245 snps ###

is_indel <- nchar(av$REF) + nchar(av$ALT) > 2 &
  substr(av$REF, 1, 1) == substr(av$ALT, 1, 1)

av[is_indel, REF := substring(REF, 2)]
av[is_indel, ALT := substring(ALT, 2)]
av[REF == "", REF := "-"]
av[ALT == "", ALT := "-"]

# Deletions (and complex substitutions) start one base after the anchor base
av[is_indel & REF != "-", START := START + 1L]

# END = START for insertions, otherwise START + length(REF) - 1
av[, END := ifelse(REF == "-", START, START + nchar(REF) - 1L)]

# head(av)
# CHR  START    REF    ALT              SNPID        rsID    END
# <char>  <int> <char> <char>             <char>      <char>  <int>
#   1:      1  11172   CTTG      - chr1:11171:CCTTG:C        <NA>  11175
# 2:      1  60351      A      G     chr1:60351:A:G  rs62637817  60351
# 3:      1  64931      G      A     chr1:64931:G:A  rs62639104  64931
# 4:      1  80346      C      G     chr1:80346:C:G        <NA>  80346
# 5:      1 598941      G      A    chr1:598941:G:A rs536226227 598941
# 6:      1 710225      T      A    chr1:710225:T:A rs185127847 710225

av <- av %>% rename(ID = SNPID)
av <- av %>% select(c("CHR", "START", "END", "REF", "ALT", "ID", "rsID"))

fwrite(av, avinput_file,
       sep = "\t", quote = FALSE, col.names = FALSE)

if (file.size(avinput_file) == 0) stop("ANNOVAR input is empty: ", avinput_file)
cat("Wrote", avinput_file, "\n")
# Wrote /project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline/output/meta_fixed.avinput 
