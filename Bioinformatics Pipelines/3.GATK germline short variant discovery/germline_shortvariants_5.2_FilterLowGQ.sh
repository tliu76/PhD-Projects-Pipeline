#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/VariantFiltration_filterlowGQ_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/VariantFiltration_filterlowGQ_%A_%a.errcd
#SBATCH --job-name=VariantFiltration
#SBATCH -c 15
#SBATCH --ntasks=1
#SBATCH --mem=59G
#SBATCH --partition=debug

#export PATH="/project/desmith_488/kerenxu/tools/gatk-4.1.8.1/:$PATH"
module load gcc/11.3.0
module load openjdk
module load gatk/4.2.6.1
module load samtools
module load htslib
cd /project/desmith_128/recordstudy/batch3/germline/HaplotypeCaller.out/after_filter

ref=/project/desmith_488/kerenxu/refs/GATK_resource_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
input_vcf=record_batch3_filtered.CalculateGenotypePosteriors.vcf.gz
output_vcf=record_batch3_filtered.CalculateGenotypePosteriors.GQ_filtered.vcf.gz
output_vcf_denovo=record_batch3_filtered.CalculateGenotypePosteriors.GQ_filtered.denovo_anno.vcf.gz

srun gatk --java-options "-Xms59G" \
    VariantFiltration \
    -R $ref \
    -V /project/desmith_128/recordstudy/batch3/germline/HaplotypeCaller.out/after_filter/$input_vcf \
    --genotype-filter-expression "GQ < 20" \
    --genotype-filter-name "lowGQ" \
    -O /project/desmith_128/recordstudy/batch3/germline/HaplotypeCaller.out/after_filter/$output_vcf

# find de novo variants would be better to have ped file
# Rarity – Allele frequency across all samples sequenced is low
srun gatk --java-options "-Xms59G" \
    VariantAnnotator \
    -R $ref \
    -V /project/desmith_128/recordstudy/batch3/germline/HaplotypeCaller.out/after_filter/$output_vcf \
    -A PossibleDeNovo \
    -O /project/desmith_128/recordstudy/batch3/germline/HaplotypeCaller.out/after_filter/$output_vcf_denovo

