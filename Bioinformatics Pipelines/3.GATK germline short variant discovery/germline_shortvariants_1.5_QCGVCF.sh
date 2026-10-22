#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/CollectVariantCallingMetrics_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/CollectVariantCallingMetrics_%A_%a.errcd
#SBATCH --job-name=CollectVariantCallingMetrics
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=160G
#SBATCH --partition=epyc-64
#SBATCH --array=1-6

#export PATH="/project/desmith_488/kerenxu/tools/gatk-4.2.2.0/:$PATH"
#srun hostname

module load gcc
module load openjdk
module load gatk/4.5.0.0


##evaluation interval list: https://42basepairs.com/browse/s3/gatk-test-data/intervals?file=wgs_evaluation_regions.hg38.interval_list

ref_dict=/project/desmith_488/tanxinli/refs/GATK_resource_bundle/Homo_sapiens_assembly38.dict 
#evaluation_interval_list=/project/desmith_488/kerenxu/refs/GATK_resource_bundle/wgs_evaluation_regions.hg38.interval_list
evaluation_interval_list=/project/desmith_488/tanxinli/refs/GATK_resource_bundle/wgs_evaluation_regions.hg38.interval_list
dbsnp_vcf=/project/desmith_488/tanxinli/refs/GATK_resource_bundle/Homo_sapiens_assembly38.dbsnp138.vcf

cd /scratch1/tanxinli/tumoroutput/bqsr.bam/HaplotypeCaller

## input:
input_vcf=(`ls *.g.vcf.gz | sed -n ${SLURM_ARRAY_TASK_ID}p`)
base=$(basename ${input_vcf} .g.vcf.gz)
prefix="sample"
sample_name=${base#"$prefix"}

# Collect variant calling metrics from GVCF output
gatk --java-options "-Xms120G -Djava.io.tmpdir=/scratch1/tanxinli/tmp/CollectVariantCallingMetrics" \
      CollectVariantCallingMetrics \
      --INPUT $input_vcf \
      --OUTPUT $base \
      --DBSNP $dbsnp_vcf \
      --SEQUENCE_DICTIONARY $ref_dict \
      --TARGET_INTERVALS $evaluation_interval_list \
      --GVCF_INPUT true \
      --TMP_DIR /scratch1/tanxinli/tmp/CollectVariantCallingMetrics
      

      