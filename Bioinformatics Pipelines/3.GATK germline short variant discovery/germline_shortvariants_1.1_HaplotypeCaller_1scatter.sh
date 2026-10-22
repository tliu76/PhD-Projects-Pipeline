#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/HaplotypeCaller_scatter_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/HaplotypeCaller_scatter_%A_%a.errcd
#SBATCH --job-name=1scatter
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=160G
#SBATCH --partition=epyc-64
#SBATCH --array=1-6


module load gcc
module load openjdk
module load gatk/4.5.0.0

#batch2

cd /scratch1/tanxinli/tumoroutput/bqsr.bam

ref=/project/desmith_488/kerenxu/refs/GATK_resource_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
interval_list=/project/desmith_488/tanxinli/refs/GATK_resource_bundle/interval_files_germline/0001-scattered.interval_list

## input:
normal_bam=(`cat normal.id | sed -n ${SLURM_ARRAY_TASK_ID}p`)
base=$(basename ${normal_bam} .aligned.duplicates_marked.recalibrated.bam)
normal_name=${base}


## output: 
output_vcf="sample$normal_name.g.vcf.gz"

# HaplotypeCaller per-sample in GVCF mode
# needs to set the mem to MAX-1GB
gatk --java-options "-Xms145G -Djava.io.tmpdir=/scratch1/tanxinli/tmp/haplotypecaller" \
      HaplotypeCaller \
      -R $ref \
      -I $normal_bam \
      -L $interval_list \
      -O /scratch1/tanxinli/tumoroutput/bqsr.bam/HaplotypeCaller/0001/$output_vcf \
      -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
      -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
      -ERC GVCF \
      --tmp-dir /scratch1/tanxinli/tmp/haplotypecaller

## for splitting intervals

#    java -Xms1g -jar /usr/gitc/picard.jar \
#      IntervalListTools \
#      SCATTER_COUNT=~{scatter_count} \
#      SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
#      UNIQUE=true \
#      SORT=true \
#     BREAK_BANDS_AT_MULTIPLES_OF=~{break_bands_at_multiples_of} \
#      INPUT=~{interval_list} \
#      OUTPUT=out
# command <<<
#    gatk --java-options "-Xms3000m -Xmx3250m" SplitIntervals \
#      -L ~{interval_list} -O  scatterDir -scatter ~{scatter_count} -R ~{ref_fasta} \
#      -mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW --interval-merging-rule OVERLAPPING_ONLY
# >>>
