#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/MergeVCFs_Haplotype_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/MergeVCFs_Haplotype_%A_%a.errcd
#SBATCH --job-name=MergeVCFs
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --partition=epyc-64
#SBATCH --array=1-6

module load gcc
module load openjdk
module load gatk/4.5.0.0


cd /scratch1/tanxinli/tumoroutput/bqsr.bam/HaplotypeCaller/0000/

## input:
input_vcf=(`ls *.g.vcf.gz | sed -n ${SLURM_ARRAY_TASK_ID}p`)
base=$(basename ${input_vcf} .g.vcf.gz)
prefix="sample"
sample_name=${base#"$prefix"}

## output:
output_vcf="$input_vcf"

## change directory:
cd /scratch1/tanxinli/tumoroutput/bqsr.bam/HaplotypeCaller

## merge:
srun gatk --java-options "-Xms100G -Djava.io.tmpdir=/scratch1/tanxinli/tmp/MergeVCFsGermline" \
    MergeVcfs \
    -I ./0000/$input_vcf \
    -I ./0001/$input_vcf \
    -I ./0002/$input_vcf \
    -I ./0003/$input_vcf \
    -I ./0004/$input_vcf \
    -O ./$output_vcf \
    --TMP_DIR /scratch1/tanxinli/tmp/MergeVCFsGermline