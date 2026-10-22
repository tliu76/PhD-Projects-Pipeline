#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/filter_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/filter_%A_%a.errcd
#SBATCH --job-name=filter
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --partition=epyc-64

#export PATH="/project/desmith_488/kerenxu/tools/htslib/htslib-1.10.2/:$PATH"
#export PATH="/project/desmith_488/kerenxu/tools/bcftools/bcftools-1.10.2/:$PATH"

module load gcc/11.3.0 openblas/0.3.20
module load openjdk
module load gatk/4.2.6.1
module load samtools
module load htslib
module load bcftools

cd /project/desmith_128/recordstudy/batch3/germline/HaplotypeCaller.out/after_filter/

input_vcf=anno.hg38_cadd_wgsnv_topmed.vcf.gz
output_vcf=anno.hg38_cadd_wgsnv_topmed.passed.vcf.gz

# first step to remove variants that did not pass VQSR
bcftools view $input_vcf -Ou -i 'FILTER=="PASS"' -o $output_vcf -O z --threads 20