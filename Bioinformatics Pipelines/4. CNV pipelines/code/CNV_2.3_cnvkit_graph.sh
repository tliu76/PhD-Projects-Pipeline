#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/CNVkit_plot_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/CNVkit_plot_%A_%a.errcd
#SBATCH --job-name=CNVkit_plot
#SBATCH -c 15
#SBATCH --ntasks=1
#SBATCH --mem=59G
#SBATCH --partition=main
#SBATCH --array=1-5

# load R
module load gcc/8.3.0 openblas/0.3.8 r/4.0.0
module load bcftools
srun hostname

export PATH="/project/desmith_488/tanxinli/SINGULARITY_CACHEDIR/:$PATH"

#cd /project/desmith_128/tanxin/scratch2/tumoroutput/batch2/cnvkit.out.anno
cd /project/desmith_128/tanxin/scratch2/tumoroutput/batch2/cnvkit.out_batch2

gene_info=/project/desmith_488/kerenxu/refs/GATK_resource_bundle/somatic_cnv/refFlat.txt

## input:
cnr_file=(`ls *.aligned.duplicates_marked.recalibrated.cnr | sed -n ${SLURM_ARRAY_TASK_ID}p`)
base=$(basename ${cnr_file} .aligned.duplicates_marked.recalibrated.cnr)
cns_file="$base.aligned.duplicates_marked.recalibrated.cns"
call_cns_file="$base.aligned.duplicates_marked.recalibrated.call.cns"
bintest_cns_file="$base.aligned.duplicates_marked.recalibrated.bintest.cns"

#modifying y axis https://github.com/etal/cnvkit/issues/332

# scatter with BAF
singularity exec --bind /project/desmith_488,/project/desmith_128  /project/desmith_488/kerenxu/SINGULARITY_CACHEDIR/cnvkit_latest.sif cnvkit.py scatter  -s "$base.aligned.duplicates_marked.recalibrated.cn{s,r}" \
-v /project/desmith_128/tanxin/scratch2/tumoroutput/batch2/mutect2_out/gather/sample${base}.vcf.gz -i $tumor_name -n $normal_name --y-min=-5 --y-max=2.5  --output ../cnv_graph/$base.baf.scatter.pdf


singularity exec --bind /project/desmith_488,/project/desmith_128  /project/desmith_488/kerenxu/SINGULARITY_CACHEDIR/cnvkit_latest.sif cnvkit.py scatter  -s 808693_T.aligned.duplicates_marked.recalibrated.cn{s,r} \
-v /project/desmith_128/tanxin/scratch2/tumoroutput/batch2/mutect2_out/gather/sample896578_T.vcf.gz -i 896578_T -n L0012_01RC_G --y-min=-2.5 --y-max=2.5  --output ../cnv_graph/$base.baf.scatter.pdf



# diagram 
singularity exec --bind /project/desmith_488,/project/desmith_128  /project/desmith_488/kerenxu/SINGULARITY_CACHEDIR/cnvkit_latest.sif cnvkit.py diagram -s Sample.cns Sample.cnr -output ../cnv_graph/cnvkit.diagram.pdf

# heatmap
singularity exec --bind /project/desmith_488,/project/desmith_128  /project/desmith_488/kerenxu/SINGULARITY_CACHEDIR/cnvkit_latest.sif cnvkit.py heatmap *.recalibrated.cns --output ../cnv_graph/cnvkit.heatmap.pdf

singularity exec --bind /project/desmith_488,/project/desmith_128  /project/desmith_488/kerenxu/SINGULARITY_CACHEDIR/cnvkit_latest.sif cnvkit.py heatmap *.recalibrated.call.cns --output ../cnv_graph/cnvkit.call.heatmap.pdf

## check arguements

singularity exec --bind /project/desmith_488,/project/desmith_128  /project/desmith_488/kerenxu/SINGULARITY_CACHEDIR/cnvkit_latest.sif cnvkit.py heatmap -h



