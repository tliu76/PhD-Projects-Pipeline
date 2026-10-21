#!/bin/bash
#SBATCH --job-name=multipanel
#SBATCH --mail-type=END
#SBATCH --account=desmith_128
#SBATCH --partition=epyc-64
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --array=1
#SBATCH --output=/scratch1/tanxinli/output/plot_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/plot_%A_%a.err

module load ver/2506 gcc/14.3.0
module load r/4.5.0


Rscript /scratch1/tanxinli/mp2prt_gwas/filtered/code/revised/9.multipanel.R