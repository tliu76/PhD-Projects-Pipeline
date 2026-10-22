#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/FinalGatherVcf_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/FinalGatherVcf_%A_%a.errcd
#SBATCH --job-name=FinalGatherVcf
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --partition=epyc-64

# export PATH="/project/desmith_488/tanxinli/tools/gatk-4.2.2.0/:$PATH"
module load gcc/11.3.0
module load openjdk
module load gatk/4.2.6.1
module load samtools
module load htslib

cd /project/desmith_128/recordstudy/batch3/germline/HaplotypeCaller.out

output_vcf_name=record_batch3.filtered.vcf.gz

# ls *.filtered.vcf.gz  | sort -V -o  inputs_filtered.list # remmeber to sort the list by chrs.

# --ignore-safety-checks makes a big performance difference so we include it in our invocation.
# This argument disables expensive checks that the file headers contain the same set of
# genotyped samples and that files are in order by position of first record.

gatk --java-options "-Xms80G" \
    GatherVcfsCloud \
    --ignore-safety-checks \
    --gather-type BLOCK \
    --input  inputs_filtered.list  \
    --output $output_vcf_name

tabix $output_vcf_name

ref_dict=/project/desmith_488/tanxinli/refs/GATK_resource_bundle/Homo_sapiens_assembly38.dict  
evaluation_interval_list=/project/desmith_488/tanxinli/refs/GATK_resource_bundle/wgs_evaluation_regions.hg38.interval_list
dbsnp_vcf=/project/desmith_488/tanxinli/refs/GATK_resource_bundle/Homo_sapiens_assembly38.dbsnp138.vcf

# Collect variant calling metrics from GVCF output
gatk --java-options "-Xms80G -Djava.io.tmpdir=/scratch1/tanxinli/tmp/FindGather " \
      CollectVariantCallingMetrics \
      --INPUT $output_vcf_name \
      --OUTPUT record_batch3 \
      --DBSNP $dbsnp_vcf \
      --SEQUENCE_DICTIONARY $ref_dict \
      --TARGET_INTERVALS $evaluation_interval_list \
      --THREAD_COUNT 8 \
      --TMP_DIR /scratch1/tanxinli/tmp/FindGather