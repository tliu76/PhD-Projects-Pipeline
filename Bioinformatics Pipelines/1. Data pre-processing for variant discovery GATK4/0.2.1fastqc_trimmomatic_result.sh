#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/fastqc_trimmomatic_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/fastqc_trimmomatic_%A_%a.errcd
#SBATCH --job-name=fastqc
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=246G
#SBATCH --partition epyc-64
#SBATCH --array=1-15

cd /project/raw_data/

fastq1=(`find . -name "*.paired_1.fq" | sed -n ${SLURM_ARRAY_TASK_ID}p`)
DIR=$(dirname ${fastq1})
base=$(basename ${fastq1} _1.fq)
fastq2="${DIR}/$(basename ${fastq1} _1.fq)_2.fq"

fastqc --noextract $fastq1 $fastq2 --outdir ./ -t 20
