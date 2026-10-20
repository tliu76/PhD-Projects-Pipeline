#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/logistic_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/logistic_%A_%a.errcd
#SBATCH --job-name=mex_logistic
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --partition=main
#SBATCH --array=1-2
#SBATCH --account=desmith_128

module load gcc

input_dir="/project2/desmith_128/datasets/For_Tanxin/MIGI_1_2_MXB_Post_QC"

cd $input_dir

### 20 PCs ###
plink2 --pfile 	final_allchrmerg  \
    --covar /project2/desmith_128/datasets/For_Tanxin/MIGI_1_2_MXB_Post_QC/covar_${SLURM_ARRAY_TASK_ID}0PCs_sex.txt \
    --threads 200 \
    --logistic hide-covar cols=chrom,pos,ref,alt1,nobs,a1countcc,a1freq,a1freqcc,beta,se,p \
    --maf 0.01 \
    --out Mexican_${SLURM_ARRAY_TASK_ID}0PCs_Sex_logistic

### 6373 samples: 807 cases  5566 controls ### 8099303 variants loaded from final_allchrmerg.pvar.
  