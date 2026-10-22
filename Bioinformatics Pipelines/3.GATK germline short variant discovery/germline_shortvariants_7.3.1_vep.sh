#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/vep_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/vep_%A_%a.errcd
#SBATCH --job-name=vep
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=40G
#SBATCH --partition=main

#source activate /project/desmith_488/kerenxu/miniconda3/envs/maf-summary

module purge
eval "$(conda shell.bash hook)"
module load gcc/11.3.0  openblas/0.3.20
conda activate vep

cd /project/desmith_128/recordstudy/batch3/germline/HaplotypeCaller.out/after_filter

## input:
input_vcf=anno.hg38_cadd_wgsnv_topmed.passed.filtered.string2float.gnomad_1kg_filtered.vcf.gz
output_vcf=anno.hg38_cadd_wgsnv_topmed.passed.filtered.string2float.gnomad_1kg_filtered.vep.vcf

ref=/project/desmith_488/kerenxu/refs/GATK_resource_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta

# run vcf2maf with VEP annotation 
vep --species homo_sapiens --assembly GRCh38 --offline --dir /project/desmith_488/tanxinli/vep/GRCh38/vep/tmp \
--fasta $ref --input_file $input_vcf --output_file $output_vcf --everything --fork 20 --cache --vcf

