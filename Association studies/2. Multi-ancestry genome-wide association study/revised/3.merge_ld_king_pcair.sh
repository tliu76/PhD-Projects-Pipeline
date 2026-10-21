#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/pcair_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/pcair_%A_%a.errcd
#SBATCH --job-name=pcair
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --partition=epyc-64
#SBATCH --account=desmith_128
#SBATCH --array=1       

cd /project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/revised/
module load  ver/2506  gcc/14.3.0
module load r/4.5.0    

Rscript /project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/revised/3.merge_ld_king_pcair.R