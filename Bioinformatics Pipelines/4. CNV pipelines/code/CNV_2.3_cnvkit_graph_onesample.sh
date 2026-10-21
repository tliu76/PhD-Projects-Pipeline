#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/1CNVkit_plot_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/1CNVkit_plot_%A_%a.errcd
#SBATCH --job-name=1CNVkit_plot
#SBATCH -c 15
#SBATCH --ntasks=1
#SBATCH --mem=59G
#SBATCH --partition=wiemels
#SBATCH --account=wiemels_260
#SBATCH --array=1

# load R
module load gcc/8.3.0 openblas/0.3.8 r/4.0.0

module load bcftools
srun hostname

export PATH="/project/desmith_488/tanxinli/SINGULARITY_CACHEDIR/:$PATH"

cd /project/desmith_128/tanxin/scratch2/tumoroutput/batch2/cnvkit.out_batch2

gene_info=/project/desmith_488/tanxinli/refs/GATK_resource_bundle/somatic_cnv/refFlat.txt

## input:
cnr_file="906484_T.aligned.duplicates_marked.recalibrated.cnr"
base=$(basename ${cnr_file} .aligned.duplicates_marked.recalibrated.cnr)
cns_file="$base.aligned.duplicates_marked.recalibrated.cns"
call_cns_file="$base.aligned.duplicates_marked.recalibrated.call.cns"
bintest_cns_file="$base.aligned.duplicates_marked.recalibrated.bintest.cns"

#scatter change y axis

singularity exec --bind /project/desmith_488,/project/desmith_128/  /project/desmith_488/tanxinli/SINGULARITY_CACHEDIR/cnvkit_latest.sif cnvkit.py scatter  -s ${base}.aligned.duplicates_marked.recalibrated.cn{s,r} \
 --y-min=-5 --y-max=2.5  --output ../cnv_graph/$base.scatter0715.pdf

#diagram

singularity exec --bind /project/desmith_488,/project/desmith_128/  /project/desmith_488/tanxinli/SINGULARITY_CACHEDIR/cnvkit_latest.sif cnvkit.py diagram  -s ${base}.aligned.duplicates_marked.recalibrated.cn{s,r} \
 --output ../cnv_graph/$base.diagram0715.pdf

#heatmap
singularity exec --bind /project/desmith_488,/project/desmith_128/  /project/desmith_488/tanxinli/SINGULARITY_CACHEDIR/cnvkit_latest.sif cnvkit.py heatmap *.recalibrated.cns --output ../cnv_graph/cnvkit_batch2.heatmap.pdf


