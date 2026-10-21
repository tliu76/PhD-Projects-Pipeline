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
PFILE="/scratch1/tanxinli/mp2prt_gwas/filtered/MP2PRT.impute"  # 你的基因型数据前缀

PFILE2="/scratch1/tanxinli/mp2prt_gwas/MP2PRT.dosage"  # before post-imputation filtering


OUTDIR="/scratch1/tanxinli/mp2prt_gwas/PRS/prs_results/pfile_before_postimp"
SCORE_DIR="/scratch1/tanxinli/mp2prt_gwas/PRS/prs_score_files"

mkdir -p $OUTDIR

exposure=(`ls /scratch1/tanxinli/mp2prt_gwas/PRS/prs_score_files/*score.txt  | sed -n ${SLURM_ARRAY_TASK_ID}p`)
exp=$(basename ${exposure} _score.txt)


# 基础PRS计算（使用所有SNP）
plink2 --pfile $PFILE2 \
       --score ${exposure} 1 2 4 header list-variants  \
       --out ${OUTDIR}/${exp}_prs



           




#### 数字含义：`1 2 4`
#- **1** = SNP ID在第1列
#- **2** = 效应等位基因在第2列
#- **4** = Beta系数（权重）在第4列

#### `sum` 的含义
#**加法模型** - 这是关键！
#- `sum`：简单相加模型（最常用）
#  - 计算公式：PRS = β₁×G₁ + β₂×G₂ + ... + βₙ×Gₙ
#  - G 是等位基因剂量（0, 1, 2）
  

 
 