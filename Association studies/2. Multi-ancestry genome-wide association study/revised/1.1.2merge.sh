#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/merge_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/merge_%A_%a.errcd
#SBATCH --job-name=merge
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --partition=epyc-64
#SBATCH --account=desmith_128

cd /scratch1/tanxinli/mp2prt_gwas/filtered/split/split_by_chr

# Load R (modify version as needed on your HPC)
module load  ver/2506  gcc/14.3.0
module load r/4.5.0

echo "Starting GDS merge at: $(date)"

Rscript - << 'EOF'
# Set custom library path
.libPaths("/home1/tanxinli/R_lib")

suppressPackageStartupMessages({
  library(GENESIS)
  library(GWASTools)
  library(SNPRelate)
  library(SeqArray)
  library(SeqVarTools)
  library(MASS)
  library(AER)
  library(dplyr)
  library(ggplot2)
  library(tidyverse)
  library(data.table)
  library(qqman)
})


# List chromosome GDS files
gds_list <- paste0("MP2PRT.filtered.chr", 1:22, ".gds")

# Merge into final GDS
seqMerge(
  gds_list,
  "/scratch1/tanxinli/mp2prt_gwas/filtered/split/MP2PRT.filtered.merged.gds",
  verbose = TRUE
)

EOF
