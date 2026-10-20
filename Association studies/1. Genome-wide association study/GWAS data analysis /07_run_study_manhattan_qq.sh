#!/bin/bash
#SBATCH --job-name=meta0924_s7_study_plots
#SBATCH --account=desmith_128
#SBATCH --partition=main
#SBATCH --time=2:00:00
#SBATCH -c 4
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --array=1-7
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline/logs/s7_study_plots_%A_%a.out
#SBATCH --error=/project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline/logs/s7_study_plots_%A_%a.err

## Step 7 - per-study Manhattan + QQ plots (07_study_manhattan_qq.R).
## One array task per study (1-7, order of 02_metal.txt). Output: output/study_plots/
set -euo pipefail

PIPE=/project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline
OUTDIR=$PIPE/output

module purge
module load r/4.5.3
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR"
cd "$PIPE"
mkdir -p logs


Rscript 07_study_manhattan_qq.R "$OUTDIR" "$SLURM_ARRAY_TASK_ID"


