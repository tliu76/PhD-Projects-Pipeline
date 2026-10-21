#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/prs_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/prs_%A_%a.errcd
#SBATCH --job-name=prs_all
#SBATCH -c 15
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --partition=epyc-64
#SBATCH --account=desmith_128
#SBATCH --array=1

# PLINK批量计算PRS
cd /scratch1/tanxinli/mp2prt_gwas/PRS/
PFILE="/scratch1/tanxinli/mp2prt_gwas/filtered/MP2PRT.impute"  # 你的基因型数据前缀
OUTDIR="/scratch1/tanxinli/mp2prt_gwas/PRS/prs_ALL_results"


score="/scratch1/tanxinli/mp2prt_gwas/PRS/ALL_19_weights.txt"


# 基础PRS计算（使用所有SNP）
plink2 --pfile $PFILE \
       --score ${score} 1 2 3 header list-variants  \
       --out ${OUTDIR}/ALL_prs



           
