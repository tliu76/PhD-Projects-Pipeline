#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/vep_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/vep_%A_%a.errcd
#SBATCH --job-name=vep
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=40G
#SBATCH --partition=main

#source activate /project/desmith_488/kerenxu/miniconda3/envs/maf-summary

module purge
eval "$(conda shell.bash hook)"
module load gcc/11.3.0  openblas/0.3.20
conda activate vep
module load gcc/11.3.0  openblas/0.3.20
vep_install -a cf -s homo_sapiens -y GRCh38 -c /project/desmith_488/tanxinli/vep/GRCh38/vep --CONVERT 
