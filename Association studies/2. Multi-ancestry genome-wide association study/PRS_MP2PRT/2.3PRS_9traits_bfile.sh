#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/prs_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/prs_%A_%a.errcd
#SBATCH --job-name=prs_9
#SBATCH -c 15
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --partition=main
#SBATCH --account=desmith_128
#SBATCH --array=1-9

# PLINK批量计算PRS
cd /scratch1/tanxinli/mp2prt_gwas/PRS/
BFILE="/scratch1/tanxinli/mp2prt_gwas/GWAS_eq/MP2PRT-QC-FINAL"  
OUTDIR="/scratch1/tanxinli/mp2prt_gwas/PRS/prs_results/bfile_result"
SCORE_DIR="/scratch1/tanxinli/mp2prt_gwas/PRS/prs_score_files"

mkdir -p $OUTDIR

exposure=(`ls /scratch1/tanxinli/mp2prt_gwas/PRS/prs_score_files/*score.txt  | sed -n ${SLURM_ARRAY_TASK_ID}p`)
exp=$(basename ${exposure} _score.txt)


# 基础PRS计算（使用所有SNP）
plink --bfile $BFILE \
       --allow-no-sex \
       --score ${exposure} 1 2 4  sum \
       --out ${OUTDIR}/${exp}_prs
