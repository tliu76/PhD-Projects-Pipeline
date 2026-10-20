### Tutorial: https://ibg.colorado.edu/cdrom2019/medland/imputation/Imputation.pdf
### https://github.com/marcoralab/postImpute/tree/master/workflow/scripts
###########################################################################################################################################
cd /scratch1/tanxinli/GUAT/revise_ref_impute
module load gcc/11.3.0

### 1. download data from TOPMed

#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch1/tanxinli/output/download_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/download_%A_%a.errcd
#SBATCH --job-name=download
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --partition=epyc-64


curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1690342/9b2d268d6b513ababb33ed5b04ee314abb48836f669ba14471e17ad0b5d2777d | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1690346/f93225cb2690c0702b078a749e0ee3b16a9152b79711e53aac343833368d17ee | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1690348/30dbb164b16b337f00fae5504b81a723231f5b7ae1d5f4b1b6c56e3addd5a921 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1690349/6c20191bbc4fb48da10192ad64b187d8c610f1d77cd1cda2afcec9c6a5ca04b2 | bash

###########################################################################################################################################

### 2. unzip 
#check md5
md5sum -c results.md5

module load p7zip
cd /scratch1/tanxinli/GUAT/revise_ref_impute
key='nkz-sWP9wQU9Cw'
for i in {1..22}; do 7z e chr_$i.zip -p${key}; done


###########################################################################################################################################
### 3. Post-imputation QC steps:
#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/guat_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/guat_%A_%a.errcd
#SBATCH --job-name=guat
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --partition=epyc-64
#SBATCH --array=1-22
#SBATCH --account=desmith_128


###Post-imputation QC steps:
Imputation R^2 greater than 0.3; Exclude Variants with overall Minor Allele Frequency (MAF) less than 0.01

### post-imputation QC steps:
module load gcc   openblas

### module load bcftools 

cd /scratch1/tanxinli/GUAT/revise_ref_impute

plink2 --vcf  /scratch1/tanxinli/GUAT/revise_ref_impute/chr"${SLURM_ARRAY_TASK_ID}".id.dose.vcf.gz  dosage=HDS \
     --max-alleles 2 \
     --min-alleles 2 \
     --maf 0.01 \
     --make-pgen \
     --out /scratch1/tanxinli/GUAT/revise_ref_impute/filtered/chr"${SLURM_ARRAY_TASK_ID}".imputed \
     --threads 200

cd /scratch1/tanxinli/GUAT/revise_ref_impute/filtered

plink2 --pfile chr"${SLURM_ARRAY_TASK_ID}".imputed   --exclude-if-info "R2<=0.3"  --maf 0.01 --make-pgen --out chr"${SLURM_ARRAY_TASK_ID}".imputed_r

### Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)

###########################################################################################################################################
### 4. Add phenotype information:

plink2 --pfile chr"${SLURM_ARRAY_TASK_ID}".imputed_r  \
      --new-id-max-allele-len 100 \
      --set-all-var-ids chr@:#:\$r:\$a \
      --pheno pheno.txt \
      --make-pgen  \
      --out chr"${SLURM_ARRAY_TASK_ID}".imputed_pheno

###########################################################################################################################################
### 5. Add phenotype information:
# create pgen_list.txt

input_dir="/scratch1/tanxinli/chla"

# mkdir merge
seq 2 22 | sed 's/^/chr/' | sed 's/$/.imputed_pheno/'  > pgen_list.txt

plink2  --pfile $input_dir/chr1.imputed_pheno \
        --pmerge-list pgen_list.txt  \
        --out $input_dir/merge/merge

###########################################################################################################################################
### 6. PCA:
plink2 --pfile merge \
    --pca 20 \
    --out merge_20pca \
    --threads 200
###########################################################################################################################################
### 7. Logistic regression:
module load gcc

plink2 --pfile merge  \
    --covar chla_pca.txt \
    --logistic hide-covar cols=chrom,pos,ref,alt1,nobs,a1countcc,a1freq,a1freqcc,beta,se,p \
    --maf 0.01 \
    --out chla_adjust_pc_sex
