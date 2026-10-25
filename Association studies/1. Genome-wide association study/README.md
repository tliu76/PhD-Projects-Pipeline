# GWAS meta-analysis

Short description of this workflow: input data, tools, and output.

## QC and imputation
1. [Pre-imputation QC steps](1.Pre-imputation.sh)
2. [Post-imputation QC steps](2.Post-imputation.sh)

## Association testing 
step 6 and 7 [2.Post-imputation.sh]

## GWAS meta-analysis downstream analysis

## Step 0: Run GWAS per cohort
- [00 CHLA/CALSEC/MEC GWAS](./GWAS%20data%20analysis%20/00_CHLA_CALSEC_MEC_GWAS_rerun.sh)
- [00 Mexican GWAS](./GWAS%20data%20analysis%20/00_Mexican_GWAS_rerun.sh)
- [001 Convert St. Jude GWAS results](./GWAS%20data%20analysis%20/001_GWAS_convert_stjude.sh)

## Step 1: Prepare covariates
- [01 Mexican covariates](./GWAS%20data%20analysis%20/01.Mexican_covariate.R)

## Step 2: Meta-analysis (METAL)
- [02 METAL configuration](02_metal.txt)
- [02 Run METAL](./GWAS%20data%20analysis%20/02_run_metal.sh)

## Step 3: Prepare summary statistics
- [03 Prepare sumstats](./GWAS%20data%20analysis%20/03_prepare_sumstats.R)

## Step 4: Annotate variants
- [04 Run ANNOVAR](./GWAS%20data%20analysis%20/04_run_annovar.sh)
- [05 Merge annotation](./GWAS%20data%20analysis%20/05_merge_annotation.R)

## Step 5: Visualization
- [06 Manhattan & QQ plots (meta-analysis)](./GWAS%20data%20analysis%20/06_manhattan_qq.R) | [run script](./GWAS%20data%20analysis%20/06_run_manhattan_qq.sh)
- [07 Manhattan & QQ plots (per study)](./GWAS%20data%20analysis%20/07_study_manhattan_qq.R) | [run script](./GWAS%20data%20analysis%20/07_run_study_manhattan_qq.sh)





## Conditional testing


## Polygenic risk scores
