#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/prep_seqdata_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/prep_seqdata_%A_%a.errcd
#SBATCH --job-name=prep_seqdata
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --partition=epyc-64
#SBATCH --account=desmith_128
#SBATCH --array=1-22

cd /project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/split
# Load R (modify version as needed on your HPC)
module load  ver/2506  gcc/14.3.0
module load r/4.5.0

COVAR=/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/covar.csv
GDS=/project2/desmith_128/tanxin/scratch1_new//mp2prt_gwas/filtered/split/MP2PRT.filtered.merged.gds
OUTDIR=/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/split/split_by_chr

CHR=${SLURM_ARRAY_TASK_ID}

mkdir -p logs
mkdir -p ${OUTDIR}

echo "=== Running chr ${CHR} ==="

Rscript /project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/revised/1.1prep_seqdata_bychr.R   ${COVAR}   ${GDS}   ${OUTDIR}  ${CHR}

echo "=== Done chr ${CHR} ==="
