#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch1/tanxinli/output/1.1smoove_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/1.1smoove_%A_%a.errcd
#SBATCH --job-name=smoove
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --partition=main
#SBATCH --array=1-6
#SBATCH --account=desmith_488

module purge 
module load gcc/11.3.0 openblas/0.3.20
module load openjdk
#module load bcftools/1.14
#module load htslib
#srun hostname

#export PATH="/project/desmith_488/kerenxu/SINGULARITY_CACHEDIR/:$PATH"
export PATH="/home1/tanxinli/:$PATH"

export SINGULARITY_CACHEDIR=/home1/tanxinli/.singularity
export SINGULARITY_BIND=/home1/tanxinli,/scratch1/tanxinli/,/project/desmith_488/kerenxu/refs/SV_resource/,/project/desmith_488/kerenxu/refs/GATK_resource_bundle/,/project/desmith_128/
cd /project/desmith_128/recordstudy/batch3/bqsr.bam

normal_dir="/project/desmith_128/recordstudy/batch3/bqsr.bam"
ref=/project/desmith_488/kerenxu/refs/GATK_resource_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
exclude_bed=/project/desmith_488/kerenxu/refs/SV_resource/exclude.cnvnator_100bp.GRCh38.20170403.bed
threads=20


## input:
#normal_bam=(`ls L00*.aligned.duplicates_marked.recalibrated.bam | sed -n ${SLURM_ARRAY_TASK_ID}p`)
base=(`awk '{print $1}' batch3.id | sed -n ${SLURM_ARRAY_TASK_ID}p`)
normal_bam="${base}.aligned.duplicates_marked.recalibrated.bam"

sample_name=${base}
cd /project/desmith_128/recordstudy/batch3/bqsr.bam
input_name="${sample_name}"
input_file="batch3.id"
# Use awk to find the matching name from the second column
matching_name=$(awk -v input="$input_name" '$1 == input {print $2}' "$input_file")
tumor_bam="${matching_name}.aligned.duplicates_marked.recalibrated.bam"
tumor_dir="/project/desmith_128/recordstudy/batch3/bqsr.bam"

# joint call matched tumor normal per sample 

singularity exec  /home1/tanxinli/smoove_latest.sif \
  smoove call \
  --outdir /scratch1/tanxinli/tumoroutput/batch3sv/smoove.out/ \
  --name sample$sample_name \
 --exclude $exclude_bed \
  --fasta $ref \
  -p $threads \
  --genotype $tumor_dir/$tumor_bam $normal_dir/$normal_bam 
