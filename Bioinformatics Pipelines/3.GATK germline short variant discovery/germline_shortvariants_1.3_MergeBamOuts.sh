#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/MergeBamOuts_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/MergeBamOuts_%A_%a.errcd
#SBATCH --job-name=MergeBamOuts
#SBATCH -c 15
#SBATCH --ntasks=1
#SBATCH --mem=70G
#SBATCH --partition=main
#SBATCH --array=1-6

export PATH="/project/desmith_488/kerenxu/tools/gatk-4.1.8.1/:$PATH"

srun hostname

cd /scratch1/tanxinli/tumoroutput/bqsr.bam/HaplotypeCaller/0000
ref=/project/desmith_488/kerenxu/refs/GATK_resource_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta

## input:
input_bam=(`ls *.bamout.bam | sed -n ${SLURM_ARRAY_TASK_ID}p`)
base=$(basename ${input_bam} .bamout.bam)
prefix="sample"
sample_name=${base#"$prefix"}

## output:
unsorted_out_bam="$base.HaplotypeCaller.unsorted.out.bam"
output_bam="$base.HaplotypeCaller.bam"

## change directory:
cd /scratch2/kerenxu/HaplotypeCaller.out/

## merge:
srun gatk --java-options "-Xms55G -Djava.io.tmpdir=/scratch/kerenxu/tmp/MergeBamOuts" \
    GatherBamFiles \
    -I /scratch1/tanxinli/tumoroutput/bqsr.bam/HaplotypeCaller/0000/$input_bam \
    -I /scratch1/tanxinli/tumoroutput/bqsr.bam/HaplotypeCaller/0001/$input_bam \
    -I /scratch1/tanxinli/tumoroutput/bqsr.bam/HaplotypeCaller/0002/$input_bam \
    -I /scratch1/tanxinli/tumoroutput/bqsr.bam/HaplotypeCaller/0003/$input_bam \
    -I /scratch1/tanxinli/tumoroutput/bqsr.bam/HaplotypeCaller/0004/$input_bam \
    -O /scratch1/tanxinli/tumoroutput/bqsr.bam/HaplotypeCaller/$unsorted_out_bam \
    -R $ref \
    --TMP_DIR /scratch/kerenxu/tmp/MergeBamOuts

# We must sort because adjacent scatters may have overlapping (padded) assembly regions, hence
# overlapping bamouts

srun gatk --java-options "-Xms55G -Djava.io.tmpdir=/scratch/kerenxu/tmp/MergeBamOuts" \
    SortSam \
    -I /scratch2/kerenxu/HaplotypeCaller.out/$unsorted_out_bam \
    -O /scratch2/kerenxu/HaplotypeCaller.out/$output_bam \
    --SORT_ORDER coordinate \
    -VALIDATION_STRINGENCY LENIENT \
    --TMP_DIR /scratch/kerenxu/tmp/MergeBamOuts

srun gatk --java-options "-Xms55G -Djava.io.tmpdir=/scratch/kerenxu/tmp/MergeBamOuts" \
    BuildBamIndex \
    -I /scratch2/kerenxu/HaplotypeCaller.out/$output_bam \
    -VALIDATION_STRINGENCY LENIENT \
    --TMP_DIR /scratch/kerenxu/tmp/MergeBamOuts