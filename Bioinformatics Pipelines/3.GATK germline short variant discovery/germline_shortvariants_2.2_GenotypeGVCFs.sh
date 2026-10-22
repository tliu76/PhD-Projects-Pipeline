#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/GenotypeGVCFs_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/GenotypeGVCFs_%A_%a.errcd
#SBATCH --job-name=GenotypeGVCFs
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --partition=epyc-64
#SBATCH --array=1-25

module load gcc/11.3.0
module load openjdk
module load gatk/4.2.6.1
module load samtools

ref=/project/desmith_488/tanxinli/refs/GATK_resource_bundle/Homo_sapiens_assembly38.fasta
batch_size=6 # 6 samples
dbsnp_vcf=/project/desmith_488/tanxinli/refs/GATK_resource_bundle/Homo_sapiens_assembly38.dbsnp138.vcf

interval=(` cat /project/desmith_128/recordstudy/batch3/germline/HaplotypeCaller.out/chrom_num.text | sed -n ${SLURM_ARRAY_TASK_ID}p`)

base=$interval
workspace_dir_name="workspace_dir_name_$base"


#workspace_dir_name="workspace_dir_name_$base"
WORKSPACE=$workspace_dir_name

# output 
output_vcf_filename="$base.scattered.vcf.gz"

# change working dir
cd /project/desmith_128/recordstudy/batch3/germline/HaplotypeCaller.out

# tar -xf $workspace_dir_name.tar

gatk --java-options "-Xms80G -Djava.io.tmpdir=/scratch1/tanxinli/tmp/GenotypeGVCFs" \
      GenotypeGVCFs \
      -R $ref \
      -O $output_vcf_filename \
      -D $dbsnp_vcf \
      -G StandardAnnotation -G AS_StandardAnnotation \
      --only-output-calls-starting-in-intervals \
      -V gendb://$WORKSPACE \
      -L $interval \
      --merge-input-intervals \
      --tmp-dir /scratch1/tanxinli/tmp/GenotypeGVCFs