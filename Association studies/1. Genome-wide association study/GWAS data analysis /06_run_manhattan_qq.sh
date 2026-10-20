#!/bin/bash
#SBATCH --job-name=meta0924_s6_plots
#SBATCH --account=desmith_128
#SBATCH --partition=main
#SBATCH --time=4:00:00
#SBATCH -c 4
#SBATCH --ntasks=1
#SBATCH --mem=48G
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline/logs/s6_plots_%j.out
#SBATCH --error=/project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline/logs/s6_plots_%j.err

## Step 6 - Manhattan + QQ plots (06_manhattan_qq.R).
## Split out of run_pipeline.sh. Needs step 5 output.
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

### step 6: Manhattan + QQ plots
# meta_fixed_significant_loci.tsv is optional (only written when there are GWS hits)
Rscript 06_manhattan_qq.R "$OUTDIR"

