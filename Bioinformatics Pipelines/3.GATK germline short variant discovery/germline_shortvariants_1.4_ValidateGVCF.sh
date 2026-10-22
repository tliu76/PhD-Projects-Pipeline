#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/ValidateGVCFs_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/ValidateGVCFs_%A_%a.errcd
#SBATCH --job-name=ValidateGVCFs
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=110G
#SBATCH --partition=epyc-64
#SBATCH --array=1-6

module load gcc
module load openjdk
module load gatk/4.5.0.0

ref=/project/desmith_488/kerenxu/refs/GATK_resource_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
interval_list=/project/desmith_488/tanxinli/refs/GATK_resource_bundle/wgs_calling_regions.hg38.interval_list
dbsnp_vcf=/project/desmith_488/tanxinli/refs/GATK_resource_bundle/Homo_sapiens_assembly38.dbsnp138.vcf

cd /scratch1/tanxinli/tumoroutput/bqsr.bam/HaplotypeCaller

## input:
input_vcf=(`ls *.g.vcf.gz | sed -n ${SLURM_ARRAY_TASK_ID}p`)
base=$(basename ${input_vcf} .g.vcf.gz)
prefix="sample"
sample_name=${base#"$prefix"}

gatk --java-options "-Xms100G -Djava.io.tmpdir=/scratch1/tanxinli/tmp/ValidateVariants" \
      ValidateVariants \
      -V $input_vcf \
      -R $ref \
      -L $interval_list \
      -gvcf \
      --validation-type-to-exclude ALLELES \
      --dbsnp $dbsnp_vcf \
      --tmp-dir /scratch1/tanxinli/tmp/ValidateVariants
