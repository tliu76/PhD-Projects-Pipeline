#!/bin/bash
#SBATCH --time=00:48:00
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch1/tanxinli/output/METAL_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/METAL_%A_%a.errcd
#SBATCH --job-name=metal_ng
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --partition=main
#SBATCH --array=1
#SBATCH --account=desmith_128

module load gcc/11.3.0 openblas/0.3.20
module load openjdk



/home1/tanxinli/METAL/build/bin/metal /project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline/02_metal.txt

