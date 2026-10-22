#severus output 

setwd("/project/desmith_128/recordstudy/cedars/SV/882949T/somatic_SVs")
.libPaths("/project/desmith_488/tanxinli/R_lib")
library(vcfR)
library(data.table)
library(dplyr)
library(tidyverse)
#library(VariantAnnotation)



# Load the VCF file
vcf <- read.vcfR("severus_somatic.vcf")

# Extract fixed fields (e.g., CHROM, POS, REF, ALT)
variants <- as.data.frame(vcf@fix)
head(variants)
fwrite(variants, "/project/desmith_128/recordstudy/cedars/SV/882949T/somatic_SVs/variants_fix.csv")

# Extract genotype data
#However, the element argument is typically used to specify one field at a time. 


gt <- extract.gt(vcf, element = "GT")  # Genotype 
gq <- extract.gt(vcf, element = "GQ")  # Genotype quality
vaf <- extract.gt(vcf, element = "VAF")  # Variant allele frequency
haf <- extract.gt(vcf, element = "hAF")  # Haplotype allele frequency
dr <- extract.gt(vcf, element = "DR")  # Reference depth
dv <- extract.gt(vcf, element = "DV")  # Variant depth

gt <- data.frame(
  GT = as.character(gt),
  GQ = as.numeric(gq),
  VAF = as.numeric(vaf),
  hAF = as.numeric(haf),
  DR = as.numeric(dr),
  DV = as.numeric(dv),
  stringsAsFactors = FALSE
)

colnames(gt) <- c("GT", "GQ", "VAF", "hAF", "DR","DV")

head(gt)


combined_df <- cbind(variants, gt)
fwrite(combined_df,"882949T_severus_somatic_SV.csv")



# Plot variant quality
variants$QUAL <- as.numeric(variants$QUAL)
# Plot the QUAL field
plot(variants$QUAL, main = "Variant Quality", xlab = "Variant Index", ylab = "Quality")

vcfR::plot(vcf@fix$QUAL)





#### read five samples:
sample_names <- c("882949T", "900129T", "906484T", "908423T")

# Initialize empty lists to store data
fix_data_list <- list()
gt_data_list <- list()
merged_data_list <- list()
# Loop through each sample, process its VCF file, and store the data
for (sample_name in sample_names) {
  
  # Construct the file path
  file_path <- paste0("/project/desmith_128/recordstudy/cedars/SV/", sample_name, "/somatic_SVs/severus_somatic.vcf")
  
  # Read the VCF file
  vcf <- read.vcfR(file_path)
  
  # Extract the fixed fields (e.g., CHROM, POS, REF, ALT)
  fix_data <- as.data.frame(vcf@fix)
  fix_data$Sample <- sample_name  # Add sample name as a column
  
  # Extract the genotype (gt) data
  
  gt <- extract.gt(vcf, element = "GT")  # Genotype 
  gq <- extract.gt(vcf, element = "GQ")  # Genotype quality
  vaf <- extract.gt(vcf, element = "VAF")  # Variant allele frequency
 # haf <- extract.gt(vcf, element = "hAF")  # Haplotype allele frequency
  dr <- extract.gt(vcf, element = "DR")  # Reference depth
  dv <- extract.gt(vcf, element = "DV")  # Variant depth
  gt_data <- data.frame(
    GT = as.character(gt),
    GQ = as.numeric(gq),
    VAF = as.numeric(vaf),
    hAF = as.character(haf),
    DR = as.numeric(dr),
    DV = as.numeric(dv),
    stringsAsFactors = FALSE
  )
  
  colnames(gt_data) <- c("GT", "GQ", "VAF", "hAF", "DR","DV")
  gt_data$Sample <- sample_name
  
  # Store in the lists
  merged_data <- cbind(fix_data, gt_data)
  merged_data_list[[sample_name]] <- merged_data  
  
}
# Combine all samples into a single data frame
final_combined_data <- do.call(rbind, merged_data_list)
fwrite(final_combined_data, "/project/desmith_128/recordstudy/cedars/SV/combined.csv")


