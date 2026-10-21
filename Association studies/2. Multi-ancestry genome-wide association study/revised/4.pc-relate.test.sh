#!/bin/bash
#SBATCH --job-name=test
#SBATCH --output=/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/revised/log/pcrel_%j.out
#SBATCH --error=/scratch1/tanxinli/error/pcrel_%j.err
#SBATCH --time=48:00:00             
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --ntasks=1        
#SBATCH --mem=30G                 
#SBATCH --partition=epyc-64      
#SBATCH --account=desmith_128

# Load R (modify version as needed on your HPC)
module load  ver/2506  gcc/14.3.0
module load r/4.5.0


Rscript /project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/revised/4.PCrelate.R