#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch1/tanxinli/output/sv_tmr_only_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/sv_tmr_only_%A_%a.errcd
#SBATCH --job-name=sv_tmr_only
#SBATCH -c 20
#SBATCH --ntasks 1
#SBATCH --mem=150G
#SBATCH --partition=epyc-64
#SBATCH --array=1-4
#SBATCH --account=desmith_488

module purge
eval "$(conda shell.bash hook)"
conda activate severus_env

# Single sample SV calling (Tumor-only)
cd /project/desmith_128/recordstudy/cedars/haplo_bam

phased_tumor_bam=(`ls *.haplotagged.bam   | sed -n ${SLURM_ARRAY_TASK_ID}p`)
base=$(basename ${phased_tumor_bam} .haplotagged.bam)
phased_vcf="/project/desmith_128/recordstudy/cedars/vcf/${base}.wf_sv.vcf.gz"

/project/desmith_488/tanxinli/Severus/severus.py --target-bam $phased_tumor_bam \
    --out-dir /project/desmith_128/recordstudy/cedars/SV/${base} \
    -t 16 \
    --phasing-vcf $phased_vcf \
    --vntr-bed /project/desmith_488/tanxinli/Severus/vntrs/human_GRCh38_no_alt_analysis_set.trf.bed \
    --PON /project/desmith_488/tanxinli/Severus/pon/PoN_1000G_hg38.tsv.gz
    
    
 #We are currently providing two PoN files for GRCh38 and chm13 which are generated using 1000G data.
#cite paper https://www.biorxiv.org/content/10.1101/2024.04.18.590093v1



# Haplotagged (phased) alignment input is highly recommended but not required. See below for the detailed instructions on how to prepare haplotagged alignments. If using haplotagged bam, the matching phased VCF file should be provided as --phasing-vcf option.