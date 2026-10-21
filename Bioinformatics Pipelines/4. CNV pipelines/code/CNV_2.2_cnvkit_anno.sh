#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/CNVkit_anno_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/CNVkit_anno_%A_%a.errcd
#SBATCH --job-name=CNVkit_anno
#SBATCH -c 15
#SBATCH --ntasks=1
#SBATCH --mem=59G
#SBATCH --partition=wiemels
#SBATCH --array=1-6
#SBATCH --account=wiemels_260

srun hostname

export PATH="/project/desmith_488/tanxinli/SINGULARITY_CACHEDIR/:$PATH"

cd /project/desmith_128/tanxin/scratch2/tumoroutput/batch2/cnvkit.out_batch2  ##anti-target target.bed anti-coverage.cnn 
gene_info=/project/desmith_488/kerenxu/refs/GATK_resource_bundle/somatic_cnv/refFlat.txt

## input:
#cnr_file='906484_T.aligned.duplicates_marked.recalibrated.cnr'
cnr_file=(`ls *.aligned.duplicates_marked.recalibrated.cnr | sed -n ${SLURM_ARRAY_TASK_ID}p`)
base=$(basename ${cnr_file} .aligned.duplicates_marked.recalibrated.cnr)
cns_file="$base.aligned.duplicates_marked.recalibrated.cns"
call_cns_file="$base.aligned.duplicates_marked.recalibrated.call.cns"
bintest_cns_file="$base.aligned.duplicates_marked.recalibrated.bintest.cns"
#Update gene names (the ‘gene’ column) in CNVkit .cnn/.cnr files, using gene annotations from another UCSC RefFlat, BED, or GFF file (e.g. refFlat.txt). This may be useful if you notice at the end of an analysis that vendor-annotated targets are not the desired gene names, and want to change the labeling without repeating the analysis with an updated target BED file.

singularity exec --bind  /project/desmith_488/,/project/desmith_128/ /project/desmith_488/kerenxu/SINGULARITY_CACHEDIR/cnvkit_latest.sif cnv_annotate.py  $gene_info  $cnr_file --output ../cnvkit.out.anno/$cnr_file

singularity exec --bind  /project/desmith_488/,/project/desmith_128/ /project/desmith_488/kerenxu/SINGULARITY_CACHEDIR/cnvkit_latest.sif cnv_annotate.py  $gene_info  $cns_file --output ../cnvkit.out.anno/$cns_file

singularity exec --bind  /project/desmith_488/,/project/desmith_128/  /project/desmith_488/tanxinli/SINGULARITY_CACHEDIR/cnvkit_latest.sif cnv_annotate.py  $gene_info  $call_cns_file --output ../cnvkit.out.anno/$call_cns_file

singularity exec --bind  /project/desmith_488/,/project/desmith_128/ /project/desmith_488/tanxinli/SINGULARITY_CACHEDIR/cnvkit_latest.sif cnv_annotate.py  $gene_info  $bintest_cns_file --output ../cnvkit.out.anno/$bintest_cns_file

