#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/ImportGVCFs_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/ImportGVCFs_%A_%a.errcd
#SBATCH --job-name=ImportGVCFs
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=5G
#SBATCH --partition=main
#SBATCH --array=1-25


module load gcc
module load openjdk
module load gatk/4.5.0.0
module load samtools
# make sample_name-map file

ref="/project/desmith_488/tanxinli/refs/GATK_resource_bundle/Homo_sapiens_assembly38.fasta"
sample_name_map="/scratch1/tanxinli/tumoroutput/bqsr.bam/HaplotypeCaller/cohort.sample_map"
batch_size=6 # 6samples

interval=(` cat /scratch1/tanxinli/tumoroutput/bqsr.bam/HaplotypeCaller/chrom_num.text | sed -n ${SLURM_ARRAY_TASK_ID}p`)
workspace_dir_name="workspace_dir_name_$interval"

# change work dir
cd /scratch1/tanxinli/tumoroutput/bqsr.bam/HaplotypeCaller

# caution: remove workspace_dir_name if there is any before the new run rm -rf ${workspace_dir_name}
# The memory setting here is very important and must be several GB lower
# than the total memory allocated to the VM because this tool uses
# a significant amount of non-heap memory for native libraries.
# Also, testing has shown that the multithreaded reader initialization
# does not scale well beyond 5 threads, so don't increase beyond that.

gatk --java-options "-Xms4G -Djava.io.tmpdir=/scratch1/tanxinli/tmp/GenomicsDBImport" \
    GenomicsDBImport \
    --genomicsdb-workspace-path $workspace_dir_name \
    --batch-size $batch_size \
    -L $interval \
    --sample-name-map $sample_name_map \
    --reader-threads 5 \
    --merge-input-intervals \
    --tmp-dir /scratch1/tanxinli/tmp/GenomicsDBImport

# tar -cf $workspace_dir_name.tar $workspace_dir_name
