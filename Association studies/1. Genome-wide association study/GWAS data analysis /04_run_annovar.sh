#!/bin/bash
#SBATCH --job-name=meta0924_s4_annovar
#SBATCH --account=desmith_128
#SBATCH --partition=main
#SBATCH --time=8:00:00
#SBATCH -c 8
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline/logs/s4_annovar_%j.out
#SBATCH --error=/project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline/logs/s4_annovar_%j.err

## Step 4 - ANNOVAR refGene/cytoBand annotation (variant ID carried via -otherinfo).
## Next: 05_run_merge_annotation.sh
set -euo pipefail

PIPE=/project2/desmith_128/tanxin/scratch1_new/replication/meta-analysis2026_0924/revised_pipeline
OUTDIR=$PIPE/output
ANNOVAR=/project2/desmith_488/tanxinli/tools/annovar/new/annovar
HUMANDB=/project2/desmith_488/tanxinli/tools/annovar/humandb

module purge
export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR"
mkdir -p "$PIPE/logs"

cd "$OUTDIR"
perl "$ANNOVAR/table_annovar.pl" meta_fixed.avinput "$HUMANDB" \
    -buildver hg38 \
    -out meta_fixed_annotated_hg38 \
    -remove \
    -protocol refGene,cytoBand \
    -operation g,r \
    -nastring . \
    -otherinfo \
    -csvout \
    -polish \
    -thread "${SLURM_CPUS_PER_TASK:-4}"

[ -s meta_fixed_annotated_hg38.hg38_multianno.csv ] || { echo "ERROR: ANNOVAR output not written" >&2; exit 1; }
echo "Step 4 finished $(date)"
