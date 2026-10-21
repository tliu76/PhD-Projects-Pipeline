#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/CNVkit_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/CNVkit_%A_%a.errcd
#SBATCH --job-name=CNVkit
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --partition=wiemels
#SBATCH --account=wiemels_260

module load gcc/11.3.0


#build reference from normal samples and infer copy number ratios
srun hostname

export PATH="/project/desmith_488/tanxinli/SINGULARITY_CACHEDIR/:$PATH"

cd /project/desmith_128/recordstudy/batch3/bqsr.bam

ref=/project/desmith_488/tanxinli/refs/GATK_resource_bundle/Homo_sapiens_assembly38.fasta
exclude_bed=/project/desmith_488/kerenxu/refs/SV_resource/exclude.cnvnator_100bp.GRCh38.20170403.bed

threads=20
tumor_bam=$(cat tumorid)
normal_bam=$(cat normalid)

gene_info=/project/desmith_488/kerenxu/refs/SV_resource/Homo_sapiens.GRCh38.104.gff3
access_bed=/project/desmith_488/kerenxu/refs/GATK_resource_bundle/somatic_cnv/cnvkit_access-excludes.hg38.bed
gene_info=/project/desmith_488/kerenxu/refs/GATK_resource_bundle/somatic_cnv/refFlat.txt

singularity exec --bind  /project/desmith_488/,/project/desmith_128/,/scratch1/tanxinli/ /project/desmith_488/tanxinli/SINGULARITY_CACHEDIR/cnvkit_latest.sif   cnvkit.py  batch $tumor_bam \
 --normal $normal_bam \
 --method wgs \
 -f $ref \
 --annotate $gene_info \
 --access $access_bed \
 --output-reference cnvkit_reference.cnn \
 -d /project/desmith_128/recordstudy/batch3/cnvkit.out_batch3/ \
 --diagram \
 --scatter  \
 -p $threads



