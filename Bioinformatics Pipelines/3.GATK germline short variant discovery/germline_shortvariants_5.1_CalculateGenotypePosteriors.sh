#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/CalculateGenotypePosteriors_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/CalculateGenotypePosteriors_%A_%a.errcd
#SBATCH --job-name=CalculateGenotypePosteriors
#SBATCH -c 15
#SBATCH --ntasks=1
#SBATCH --mem=59G
#SBATCH --partition=debug

# export PATH="/project/desmith_488/kerenxu/tools/gatk-4.2.2.0/:$PATH"

module load gcc/11.3.0
module load openjdk
module load gatk/4.2.6.1
module load samtools
module load htslib

cd /project/desmith_128/recordstudy/batch3/germline/HaplotypeCaller.out

ref=/project/desmith_488/tanxinli/refs/GATK_resource_bundle/Homo_sapiens_assembly38.fasta
population_vcf=/project/desmith_488/kerenxu/refs/GATK_resource_bundle/1000G_phase3_v4_20130502.sites.hg38.vcf
input_vcf=record_batch3.filtered.vcf.gz
output_vcf=record_batch3_filtered.CalculateGenotypePosteriors.vcf.gz

gatk --java-options "-Xms59G" \
    CalculateGenotypePosteriors \
    -R $ref \
    -V $input_vcf \
    -supporting $population_vcf \
    -O ./after_filter/$output_vcf
