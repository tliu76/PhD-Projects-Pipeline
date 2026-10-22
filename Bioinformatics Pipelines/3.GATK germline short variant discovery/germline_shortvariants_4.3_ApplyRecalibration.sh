#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch1/tanxinli/output/ApplyVQSR_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/ApplyVQSR_%A_%a.errcd
#SBATCH --job-name=ApplyVQSR
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --partition=epyc-64
#SBATCH --array=1-25

module load gcc/11.3.0
module load openjdk
module load gatk/4.2.6.1
module load samtools
module load htslib

cd /project/desmith_128/recordstudy/batch3/germline/HaplotypeCaller.out

# input 
input_vcf=(`ls *.scattered.vcf.gz | sed -n ${SLURM_ARRAY_TASK_ID}p`)
base=$(basename ${input_vcf} .scattered.vcf.gz)
indels_recalibration=smoking.38cohort.indels.recal
indels_tranches=smoking.38cohort.indels.tranches
indel_filter_level=99.0
snps_recalibration=smoking.38cohort.snps.recal
snps_tranches=smoking.38cohort.snps.tranches
snp_filter_level=99.7

# output
tmp_indel_recalibrated_vcf=$base.tmp.indel.recalibrated.vcf
recalibrated_vcf_filename=$base.filtered.vcf.gz

srun gatk --java-options "-Xms80G -Djava.io.tmpdir=/scratch1/tanxinli/tmp/ApplyRecalibration " \
      ApplyVQSR \
      -O $tmp_indel_recalibrated_vcf \
      -V $input_vcf \
      --recal-file $indels_recalibration \
      --tranches-file $indels_tranches \
      --truth-sensitivity-filter-level $indel_filter_level \
      --create-output-variant-index true \
      -mode INDEL \
      --tmp-dir /scratch1/tanxinli/tmp/ApplyRecalibration

srun gatk --java-options "-Xms80G -Djava.io.tmpdir=/scratch1/tanxinli/tmp/ApplyRecalibration " \
      ApplyVQSR \
      -O $recalibrated_vcf_filename \
      -V $tmp_indel_recalibrated_vcf \
      --recal-file $snps_recalibration \
      --tranches-file $snps_tranches \
      --truth-sensitivity-filter-level $snp_filter_level \
      --create-output-variant-index true \
      -mode SNP \
      --tmp-dir /scratch1/tanxinli/tmp/ApplyRecalibration