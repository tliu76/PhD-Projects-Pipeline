########### DO NOT RUN ###################
########### DO NOT RUN ###################
########### DO NOT RUN ###################
########### DO NOT RUN ###################
########### DO NOT RUN ###################

#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch/kerenxu/output/CharGer_%A_%a.out
#SBATCH --error=/scratch/kerenxu/error/CharGer_%A_%a.errcd
#SBATCH --job-name=CharGer
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=246G
#SBATCH --partition epyc-64

# install  conda create --name CharGer python=2.7

source activate /project/desmith_488/kerenxu/miniconda3/envs/CharGer

cd /scratch/kerenxu/wgs.smoking.out/HaplotypeCaller.out

input_vcf=/scratch/kerenxu/wgs.smoking.out/HaplotypeCaller.out/anno.hg38_cadd_wgsnv_topmed.passed.filtered.string2float.gnomad_1kg_filtered.vep.vcf
output_tsv=/scratch/kerenxu/wgs.smoking.out/HaplotypeCaller.out/anno.hg38_cadd_wgsnv_topmed.passed.filtered.string2float.gnomad_1kg_filtered.vep.charger.tsv
clinvar_data=/project/desmith_488/kerenxu/refs/charger_ref/clinvar_alleles.single.b38.tsv.gz

# run charger
/project/desmith_488/kerenxu/miniconda3/envs/CharGer/bin/CharGer-0.5.4/bin/charger -f $input_vcf -o $output_tsv -l --mac-clinvar-tsv $clinvar_data