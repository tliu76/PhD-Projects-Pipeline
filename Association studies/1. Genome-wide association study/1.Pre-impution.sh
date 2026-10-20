# Guatamalan study

module load   legacy/CentOS7  gcc/11.3.0  openblas/0.3.21
module load plink/1.9-beta7


## 0. sex inconsistency # Require at least one polymorphic X chromosome locus. 
plink --bfile keepiid --check-sex  

grep "PROBLEM" plink.sexcheck| awk '{print$1,$2}'> sex_discrepancy.txt

# 1.SNPs were filtered out on the basis of call rate (<98%), 

plink \
--bfile all \
--set-missing-var-ids @:#[b37] \
 --geno 0.02 \
 --make-bed \
 --out all_geno

#2. minor allele frequency (MAF<0.01), 
plink --bfile all_geno --maf 0.01 --make-bed --out all_geno_maf

# 3. genome-wide relatedness (PI_HAT>0.20), 
# plink --bfile  all_geno_maf  --genome --min 0.2 --out pihat_min0.2
# plink --bfile all_geno_maf  --remove 0.2_pihat.txt --make-bed --out  all_geno_maf_pi_hat

plink --bfile all_geno_maf  --genome --out all_geno_maf_relatedness 

## check relatedness file:
awk '$10 > 0.20 {print $1, $2,$3, $4, $10}' all_geno_maf_relatedness.genome > high_relatedness_file

sort -k$(awk '{print NF; exit}' high_relatedness_file) high_relatedness_file > sort_highrelate.txt
# no highly related pairs in this Guatamalan study
module purge
module load    gcc  openblas
module load r
Rscript --no-save Relatedness.R


#4. genome heterozygosity rate (mean heterozygosity± 6Std)
plink --bfile all_geno_maf --het --out R_check


## The following code generates a list of individuals who deviate more than 6 standard deviations from the heterozygosity rate mean.
module load gcc/11.3.0  openblas/0.3.20
module load r/4.2.1
Rscript --no-save heterozygosity_outliers_list_6sd.R
Rscript --no-save check_heterozygosity_rate.R
# Output of the command above: fail-het-qc.txt
# Adapt this file to make it compatible for PLINK, by removing all quotation marks from the file and selecting only the first two columns.
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt

# Remove heterozygosity rate outliers.

plink --bfile all_geno_maf --remove het_fail_ind.txt --make-bed --out all_geno_maf_hetero6sd

#5 deviation from Hardy-Weinberg equilibrium in controls (P<10-5). 
plink \
 --bfile all_geno_maf \
 --hwe 0.00001 \
 --make-bed \
 --out all_geno_maf_hwe

#6 Samples with call rate < 95% were also removed.  1 people removed due to missing genotype data,--mind).
plink \
 --bfile all_geno_maf_hwe \
 --mind 0.05 \
 --make-bed \
 --out  all_geno_maf_hwe_mind

# Total genotyping rate is 0.999144.
# 589564 variants and 87 people pass filters and QC.
# Among remaining phenotypes, 38 are cases and 49 are controls.

######################  TripleliftOver from hg 19 to hg38 coordinates ################################# 
perl /project/desmith_488/tanxinli/tools/triple-liftOver-main/tripleliftover_v133.pl \
 --bim /scratch1/tanxinli/GUAT/GRECERE_GWAS_data/combined/all_geno_maf_hwe_mind.bim \
 --base hg19 \
 --target hg38 \
 --outprefix triple

#Then use the variant name in the output to exclude any variant that is not in the triple-liftOver output 
# (these variants either failed the lift over or were lifted to a different chromosome). 
# It can be executed with a command in Linux followed by a PLINK command:

perl /project/desmith_488/tanxinli/tools/triple-liftOver-main/tripleliftover_v133.pl \
 --bim /scratch1/tanxinli/GUAT/GRECERE_GWAS_data/combined/all_geno_maf_hwe_mind.bim \
 --base hg19 \
 --target hg38 \
 --outprefix test

awk '{if (NR>1) print $1}' test.invr.txt > variant.name.in.liftOver.output

plink --bfile /scratch1/tanxinli/GUAT/GRECERE_GWAS_data/combined/all_geno_maf_hwe_mind \
      --extract variant.name.in.liftOver.output \
      --make-bed \
      --out all_geno_maf_hwe_mind.lifted

awk '{if (NR>1) print $1 "\t" $3}' test.invr.txt > variant.name.hg38.position

plink --bfile all_geno_maf_hwe_mind.lifted \
      --update-map variant.name.hg38.position \
      --make-bed --out all_geno_maf_hwe_mind.lifted.hg38

awk '{if (NR>1 && $5 == "inverted") print $1}' test.invr.txt > inverted.SNV.name
plink --bfile all_geno_maf_hwe_mind.lifted.hg38 --flip inverted.SNV.name --make-bed \
     --out all_geno_maf_hwe_mind.lifted.hg38.final

#589195 variants and 87 people pass filters and QC.


########################################################################################################
### If there allele mismatches errors from the imputation server: use the following command: 
### REF/ALT switches : fix the ref/alt allele switch by comparing the bim file to human ref .fa file
# 1. Need to first run the fix_ref_alt.switch.R script to generate file_contains_ref.alt_switch_SNP file, and then update PLINK
# 2. Correct allele switch based on .fa file, (including XY, NOT FOR THIS STUDY)


plink    --bfile all_geno_maf_hwe_mind.lifted.hg38.final \
         --a2-allele file_contains_ref.alt_switch_SNP 2 1 \
         --make-bed --out all_geno_maf_hwe_mind.lifted.hg38_ref

############################################################
# Prepare for imputation
Path to the file: 
top_path='/scratch1/tanxinli/GUAT/GRECERE_GWAS_data/combined/topmed'
name='all_geno_maf_hwe_mind.lifted.hg38_ref'

#### For auto chr1-22 
module load bcftools
module load htslib

### Convert into chr1-chr22 vcf.gz  
for chr in {1..22}
do
plink --bfile $name \
 --chr $chr --recode vcf-iid --keep-allele-order --output-chr chrM \
 --out ${top_path}/${name}_hg38_chr$chr

 bcftools sort  ${top_path}/${name}_hg38_chr"$chr".vcf | bgzip -c > ${top_path}/${name}_hg38_chr"$chr".vcf.gz
 rm ${top_path}/${name}_hg38_chr"$chr".vcf
done   

### remove ambiguous SNPs (I think there is no need to remove those variants if we know the strand)
for chr in {1..22}
do
 bcftools view -e 'ALT="N"' ${top_path}/${name}_hg38_chr"$chr".vcf.gz -Ou | \
 bcftools view -e 'ALT="A" && REF="T" || ALT="T" && REF="A" || ALT="C" && REF="G" || ALT="G" && REF="C"' \
  -Oz -o ${top_path}/${name}_hg38_amb_chr"$chr".vcf.gz
done 

### upload for imputation 



