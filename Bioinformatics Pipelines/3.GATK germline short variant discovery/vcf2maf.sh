#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/4.vcf2maf_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/4.vcf2maf_%A_%a.errcd
#SBATCH --job-name=vcf2maf
#SBATCH -c 15
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --partition=epyc-64
#SBATCH --account=desmith_128
#SBATCH --array=1

module purge
module load gcc/11.3.0  openblas/0.3.20
module load samtools
eval "$(conda shell.bash hook)"

conda activate /home1/tanxinli/.conda/envs/vep

ref=/project/desmith_488/kerenxu/refs/GATK_resource_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta

cd /project/desmith_128/recordstudy/batch3/germline/HaplotypeCaller.out/after_filter
input_vep="anno.hg38_cadd_wgsnv_topmed.passed.filtered.string2float.gnomad_1kg_filtered.vep.vcf"
base=$(basename ${input_vep} .vep.vcf)
output_maf=${base}.vep.maf

perl /home1/tanxinli/.conda/envs/vep/bin/vcf2maf.pl --input-vcf  $input_vep  --output-maf $output_maf \
--ref-fasta "$ref"  \
--vep-path "/home1/tanxinli/.conda/envs/vep/bin" --vep-data "/project/desmith_488/tanxinli/vep/GRCh38/vep" \
--species homo_sapiens --vcf-normal-id=6020463_A_G L0019_01RC_G L0018_01RC_G L0022_01RC_G L0023_01RC_G L0004_01R_G  --vep-forks 20  --vep-overwrite --ncbi-build GRCh38

conda deactivate
