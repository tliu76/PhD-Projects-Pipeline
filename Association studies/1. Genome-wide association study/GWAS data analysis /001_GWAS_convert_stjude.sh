#!/bin/bash
#SBATCH --job-name=meta0924_s1_stjude
#SBATCH --account=desmith_128
#SBATCH --partition=main
#SBATCH --time=1:00:00
#SBATCH -c 2
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline/logs/s1_stjude_%j.out
#SBATCH --error=/project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline/logs/s1_stjude_%j.err

## Step 1 - Convert the St_Jude_MESA R2-filtered summary stats (Li lab, 2024) to METAL input.
## The source file has OR / LOG(OR)_SE only, so:
##   BETA         = log(OR)
##   SE           = LOG(OR)_SE
##   OTHER_ALLELE = the one of REF/ALT that is not A1
##   A1_FREQ      = A1 frequency from the earlier stjude_0708 run (matched on ID, flipped
##                  if that run's A1 is the other allele); variants without a match are dropped
##   OBS_CT       = $STJUDE_N (not in the source file; default 1928)
## Usage: sbatch 001_GWAS_convert_stjude.sh   (override N with --export=ALL,STJUDE_N=<n>)
set -euo pipefail

PIPE=/project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline
SRC=/project2/lishaobo_902/publicDataShare/forTanxin/LAT_GWAS_SumStats/filtered/St_Jude_MESA_Rsq_filtered_finalized_GWAS_summary_stats.tsv
FREQ=/project2/desmith_128/datasets/St_Jude_MESA_SNP_6/Imputed/Results/binary/remove_dup2/merged_0325/remove_outlier/stjude_0708.otherallele.txt
OUT=$PIPE/input/St_Jude_MESA_Rsq_filtered.otherallele.txt
STJUDE_N=${STJUDE_N:-1928}   # = OBS_CT in stjude_0708; new/old SE ratio ~0.98, same samples

module purge
export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR" "$PIPE/logs"

echo "==== Step 1: St Jude input ($(date)) ===="
# source: chr pos ID REF ALT A1 OR LOG(OR)_SE L95 U95 Z_STAT P
# freq:   #CHROM POS ID REF ALT A1 ... A1_FREQ($9) ... OTHER_ALLELE($16)
awk -F'\t' -v OFS='\t' -v n="$STJUDE_N" '
  NR == FNR { if (FNR > 1) { a1[$3] = $6; oa[$3] = $16; f[$3] = $9 }; next }
  FNR == 1  { print "#CHROM","POS","ID","REF","ALT","A1","A1_FREQ","OBS_CT","BETA","SE","P","OTHER_ALLELE"; next }
  {
    tot++
    other = ($6 == $5) ? $4 : $5
    if (!($3 in f) || f[$3] == "NA") { nofreq++; next }
    if (a1[$3] == $6 && oa[$3] == other)      freq = f[$3]
    else if (a1[$3] == other && oa[$3] == $6) { freq = 1 - f[$3]; flip++ }
    else { mism++; next }
    if ($7 <= 0 || $8 == "NA" || $12 == "NA") { bad++; next }
    print $1, $2, $3, $4, $5, $6, freq, n, log($7), $8, $12, other
    kept++
  }
  END {
    printf "Source variants: %d\nKept: %d (A1 flipped vs stjude_0708: %d)\n", tot, kept, flip > "/dev/stderr"
    printf "Dropped - no frequency in stjude_0708: %d | allele mismatch: %d | missing OR/SE/P: %d\n", nofreq, mism, bad > "/dev/stderr"
  }' "$FREQ" "$SRC" > "$OUT"

[ -s "$OUT" ] || { echo "ERROR: $OUT not written" >&2; exit 1; }
head -3 "$OUT"
echo "Step 1 finished $(date)"
