#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/logistic_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/logistic_%A_%a.errcd
#SBATCH --job-name=chla_calsec_mec_logistic
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --partition=main
#SBATCH --array=1-2
#SBATCH --account=desmith_128

module load gcc

input_dir="/project2/desmith_128/new_tx/chla_calsec_mec/2vcf/filter/merge/"

cd $input_dir

# covar_10pca_chla_calsec_mec.txt
	
	

### 20 PCs ###
plink2  --pfile /project2/desmith_128/new_tx/chla_calsec_mec/2vcf/filter/merge/merge_remove_outlier  \
    --covar /project2/desmith_128/new_tx/chla_calsec_mec/2vcf/filter/merge/covar_${SLURM_ARRAY_TASK_ID}0pca_chla_calsec_mec.txt \
    --threads 200 \
    --logistic hide-covar cols=chrom,pos,ref,alt1,nobs,a1countcc,a1freq,a1freqcc,beta,se,p \
    --maf 0.01 \
    --out /project2/desmith_128/new_tx/chla_calsec_mec/2vcf/filter/merge/CHLA_CALSEC_MEC_${SLURM_ARRAY_TASK_ID}0PCs_Sex_logistic


