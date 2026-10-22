#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=/scratch1/tanxinli/output/IndelsVariantRecalibrator_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/IndelsVariantRecalibrator_%A_%a.errcd
#SBATCH --job-name=IndelsVariantRecalibrator
#SBATCH -c 20
#SBATCH --ntasks=1
#SBATCH --mem=206G
#SBATCH --partition=epyc-64

#export PATH="/project/desmith_488/kerenxu/tools/gatk-4.1.8.1/:$PATH"


module load gcc/11.3.0
module load openjdk
module load gatk/4.2.6.1
module load samtools
module load htslib

# references
mills_resource_vcf=/project/desmith_488/tanxinli/refs/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
axiomPoly_resource_vcf=/project/desmith_488/tanxinli/refs/GATK_resource_bundle/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
dbsnp_resource_vcf=/project/desmith_488/tanxinli/refs/GATK_resource_bundle/Homo_sapiens_assembly38.dbsnp138.vcf

# input
sites_only_variant_filtered_vcf=smoking.38cohort.SitesOnlyGatherVcf.vcf.gz
max_gaussians=4

# output
recalibration_filename=smoking.38cohort.indels.recal
tranches_filename=smoking.38cohort.indels.tranches

cd /project/desmith_128/recordstudy/batch3/germline/HaplotypeCaller.out

gatk --java-options "-Xms180G -Djava.io.tmpdir=/scratch1/tanxinli/tmp/ " \
    VariantRecalibrator \
    -V $sites_only_variant_filtered_vcf \
    --trust-all-polymorphic \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
    -mode INDEL \
    --max-gaussians $max_gaussians \
    -resource:mills,known=false,training=true,truth=true,prior=12 $mills_resource_vcf \
    -resource:axiomPoly,known=false,training=true,truth=false,prior=10 $axiomPoly_resource_vcf \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2 $dbsnp_resource_vcf \
    -O $recalibration_filename \
    --tranches-file $tranches_filename \
    --tmp-dir /scratch1/tanxinli/tmp


