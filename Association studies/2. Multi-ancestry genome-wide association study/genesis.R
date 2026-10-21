################################################################################
# GENESIS GWAS Workflow for Somatic Deletion Analysis
# Author: Tanxin Liu
# 
# This script performs GWAS analysis for 4 deletion-related outcomes:
# 1. Total deletions (all types)
# 2. Off-target RAG deletions (binary and count)
# 3. Total RAG deletions
# 4. Proportion of off-target RAG deletions
################################################################################
This vignette provides a description of how to use the GENESIS package to run genetic association tests on array (SNP) data. 
GENESIS uses mixed models for genetic association testing, as PC-AiR PCs can be used as fixed effect covariates to adjust for population stratification, 
and a kinship matrix (or genetic relationship matrix) estimated from PC-Relate can be used to account for phenotype correlation due to genetic similarity among samples.
# ==============================================================================
# SECTION 0: SETUP AND DEPENDENCIES
# ==============================================================================
.libPaths("/home1/tanxinli/R_lib")

gd_new.covariate_cohort <- fread("/scratch1/tanxinli/project3_vdj/RAG_Recombination/large_data/gd_new.covariate.rag.csv" )


suppressPackageStartupMessages({
  library(GENESIS)
  library(GWASTools)
  library(SNPRelate)
  library(SeqArray)
  library(SeqVarTools)
  library(MASS)
  library(AER)
  library(dplyr)
  library(ggplot2)
  library(tidyverse)
  library(data.table)
  library(qqman)
})

# ==============================================================================
# SECTION 1: FILE PATHS AND PARAMETERS
# ==============================================================================
# Set your working directory
setwd("/scratch1/tanxinli/mp2prt_gwas/filtered/")



# gd_new.covariate_cohort

# an R matrix of SNP genotype data
# a GDS file
# PLINK files
# Example R code for creating a GenotypeData object is presented below. Much more detail can be found in the GWASTools package reference manual.
# https://bioconductor.org/packages/3.21/bioc/vignettes/GENESIS/inst/doc/assoc_test.html#preparing-scan-annotation-data

/scratch1/tanxinli/mp2prt_gwas/filtered

#### Reading in Genotype Data:
setwd("/scratch1/tanxinli/mp2prt_gwas/filtered")
# seqVCF2GDS("MP2PRT.impute.vcf.gz", "MP2PRT.impute.gds", storage.option = "LZMA_RA")
# seqSummary("MP2PRT.impute.gds")


### prepare my phenotype model
# https://cran.r-project.org/web/packages/RNOmni/vignettes/RNOmni.html
# https://cran.r-project.org/web/packages/bestNormalize/vignettes/bestNormalize.html#the-box-cox-transformation
# Rank Normal Omnibus Association Test
# The rank based inverse normal transformation (INT)
# install.packages("RNOmni")
library(RNOmni)
library(bestNormalize)
covar <- gd_new.covariate_cohort %>% select(c("sample", "Age","Sex", "n.nonig.least.rag.del"  ,"prop.nonig.least.rag.del" ))
covar$n.nonig.least.rag.del_INT2 <- RankNorm(covar$n.nonig.least.rag.del)

# Rank-based inverse normal transformation
rankNorm_obj <- orderNorm(covar$n.nonig.least.rag.del)
covar$n.nonig.least.rag.del_INT <- predict(rankNorm_obj)

df <- data.frame(
  Original = covar$n.nonig.least.rag.del,
  RankNorm = covar$n.nonig.least.rag.del_INT2
) %>%
  pivot_longer(cols = everything(),
               names_to = "Data",
               values_to = "Measurement")

ggplot(df, aes(x = Measurement, fill = Data)) +
  geom_density(alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("Original" = "coral", "RankNorm" = "royalblue")) +
  labs(
    title = expression("Phenotype Before and After INT"),
    x = "Measurement",
    y = "Density"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )
### Plot the density:

### Linear mixed model GWAS 
# read in GDS data

# Once the PLINK files have been converted to a GDS file, then a GenotypeData object can be created as described above.
# gdsfile <- "/scratch1/tanxinli/mp2prt_gwas/filtered/MP2PRT.impute.gds"
# The issue is that seqVCF2GDS() from the SeqArray package creates a SeqArray GDS file, which has a different structure than SNPRelate GDS files. The node is called "variant.id" in SeqArray, not "snp.id".
showfile.gds(closeall = FALSE)
showfile.gds(closeall = TRUE)
seqClose(gds)
library(GENESIS)
library(GWASTools)

common_samples <- intersect(sample_ids, covar$sample) ##1491


# Filter phenotype data
covar<- covar[sample %in% common_samples]
covar <- covar[order(match(sample, sample_ids))]



# ==============================================================================
# STEP 4: Calculate Principal Components (for population structure)
# ==============================================================================
# Filter for LD pruning (use high quality variants)
seqSetFilter(gds, 
             sample.id = common_samples,
             variant.sel = NULL)

# LD pruning
set.seed(123)
snpset <- snpgdsLDpruning(gds,
                          method = "corr",
                          ld.threshold = 0.2,
                          maf = 0.05,
                          missing.rate = 0.05,
                          autosome.only = TRUE,
                          num.thread = 4)

snpset_pruned <- unlist(snpset, use.names = FALSE)
cat(sprintf("Variants for PCA: %d\n", length(snpset_pruned)))
# Variants for PCA: 105438
# Calculate PCA
pca <- snpgdsPCA(gds,
                 snp.id = snpset_pruned,
                 sample.id = common_samples,
                 num.thread = 4,
                 eigen.cnt = 10)

# ### 
# Principal Component Analysis (PCA) on genotypes:
#   Calculating allele counts/frequencies ...
# # of selected variants: 105,438
# # of samples: 1,491
# # of SNVs: 105,438
# using 4 threads
# # of principal components: 10


# Create PC dataframe
pcs <- data.frame(
  IID = pca$sample.id,
  PC1 = pca$eigenvect[,1],
  PC2 = pca$eigenvect[,2],
  PC3 = pca$eigenvect[,3],
  PC4 = pca$eigenvect[,4],
  PC5 = pca$eigenvect[,5],
  PC6 = pca$eigenvect[,6],
  PC7 = pca$eigenvect[,7],
  PC8 = pca$eigenvect[,8],
  PC9 = pca$eigenvect[,9],
  PC10 = pca$eigenvect[,10]
)

# Plot PCs
png("pca_plot.png", width=800, height=600, res=100)
plot(pca$eigenvect[,1], pca$eigenvect[,2],
     xlab = sprintf("PC1 (%.2f%%)", pca$varprop[1]*100),
     ylab = sprintf("PC2 (%.2f%%)", pca$varprop[2]*100),
     main = "Principal Component Analysis",
     pch = 19, col = rgb(0.2, 0.4, 0.8, 0.6))  # Blue color
dev.off()

covar$IID <- covar$sample
covar$sample <- NULL
# Merge phenotype with PCs
pheno <- merge(covar, pcs, by = "IID")

# ==============================================================================
# STEP 4: Calculate Genetic Relationship Matrix (GRM)
# ==============================================================================

cat("Calculating Genetic Relationship Matrix...\n")

# Use LD-pruned SNPs for GRM
#### Genetic Relationship Matrix (GRM, GCTA):
grm <- snpgdsGRM(gds,
                 method = "GCTA",
                 snp.id = snpset_pruned,
                 sample.id = common_samples,
                 num.thread = 4)

# Create covariance matrix
kinship <- grm$grm
rownames(kinship) <- colnames(kinship) <- grm$sample.id

# Convert to format for GENESIS
kinship_mat <- as.matrix(kinship)

# ==============================================================================
# STEP 5: Prepare annotation data for GENESIS
# ==============================================================================

if(any(duplicated(pheno$IID))) {
  cat("WARNING: Duplicate IDs found in phenotype file. Keeping first occurrence.\n")
  dup_ids <- pheno$IID[duplicated(pheno$IID)]
  cat(sprintf("Duplicate IDs: %s\n", paste(unique(dup_ids), collapse=", ")))
  pheno <- pheno[!duplicated(pheno$IID), ]
}

# Recompute common samples after removing duplicates
common_samples <- intersect(sample_ids, pheno$IID)
# Ensure phenotype data is ordered correctly
pheno_ordered <- pheno[match(common_samples, covar$IID)]

## recode M/F for sex

pheno_ordered <- pheno_ordered %>%
  mutate(sex1 = case_when(
    Sex == "Male" ~ "M",
    Sex == "Female" ~ "F",
    TRUE ~ NA_character_  # handle missing or other values
  ))

pheno_ordered$sample.id <- pheno_ordered$IID
pheno_ordered <- pheno_ordered[!duplicated(pheno_ordered$sample.id), ]
scanAnnot <- ScanAnnotationDataFrame(
  data.frame(
    scanID = unique(pheno_ordered$sample.id),
    phenotype = pheno_ordered$n.nonig.least.rag.del_INT2,
    sex = pheno_ordered$sex1,  # MODIFY based on your covariates
    age = pheno_ordered$Age,  # MODIFY based on your covariates
    PC1 = pheno_ordered$PC1,
    PC2 = pheno_ordered$PC2,
    PC3 = pheno_ordered$PC3,
    PC4 = pheno_ordered$PC4,
    PC5 = pheno_ordered$PC5,
    PC6 = pheno_ordered$PC6,
    PC7 = pheno_ordered$PC7,
    PC8 = pheno_ordered$PC8,
    PC9 = pheno_ordered$PC9,
    PC10 = pheno_ordered$PC10,
    stringsAsFactors = FALSE
  )
)

# ==============================================================================
# STEP 6: Fit null model
# ==============================================================================


# Define model formula
# MODIFY covariates based on your data
nullmod <- fitNullModel(
  scanAnnot,
  outcome = "phenotype",
  covars = c( "PC1", "PC2", "PC3", "PC4", "PC5", 
             "PC6", "PC7", "PC8", "PC9", "PC10"),
  cov.mat = kinship_mat,
  family = gaussian,  # Use binomial() for binary traits
  verbose = TRUE
)

# Save null model
saveRDS(nullmod, "null_model.rds")

# ==============================================================================
# STEP 7: Run GWAS
# ==============================================================================

# Reset filters for full GWAS
seqResetFilter(gds)
# of selected samples: 1,505
# of selected variants: 8,202,031
seqSetFilter(gds, sample.id = common_samples)
# of selected samples: 1,491

# Apply QC filters
seqSetFilterCond(gds,
                 maf = 0.01,      # MAF > 1%
                 mac = 10,         # MAC >= 10
                 missing.rate = 0.1)  # Missing rate < 10%
# 8,174,068

# Create SeqVarData object
seqData <- SeqVarData(gds, sampleData = scanAnnot)

# Run association test
# For large datasets, this can take a while
assoc <- assocTestSingle(
  seqData,
  null.model = nullmod,
  test = "Score",  # Score test is fastest; can use "Wald" for effect sizes
  verbose = TRUE
)

# Save results
fwrite(assoc, "gwas_results.txt", sep = "\t", quote = FALSE)

# ==============================================================================
# STEP 8: Visualizations and Summary
# ==============================================================================

cat("Creating visualizations...\n")

# Load results if needed
results <- assoc

# Add -log10 p-values
results$log10p <- -log10(results$Score.pval)

# Calculate lambda (genomic inflation factor)
chisq <- qchisq(1 - results$Score.pval, 1)
lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
cat(sprintf("\nGenomic inflation factor (lambda): %.4f\n", lambda))

# Manhattan plot
library(qqman)

manhattan_data <- data.frame(
  CHR = results$chr,
  BP = results$pos,
  SNP = results$variant.id,
  P = results$Score.pval
)

# Remove missing values
manhattan_data <- manhattan_data[complete.cases(manhattan_data), ]

png("manhattan_plot.png", width = 1400, height = 700, res = 100)
manhattan(manhattan_data,
          chr = "CHR",
          bp = "BP",
          snp = "SNP",
          p = "P",
          main = "GWAS Manhattan Plot",
          suggestiveline = -log10(1e-5),
          genomewideline = -log10(5e-8),
          col = c("blue4", "orange3"))
dev.off()

# QQ plot
png("qq_plot.png", width = 700, height = 700, res = 100)
qq(manhattan_data$P, main = "QQ Plot")
text(x = 3, y = max(-log10(manhattan_data$P), na.rm = TRUE) * 0.9,
     labels = sprintf("λ = %.4f", lambda), cex = 1.5)
dev.off()

# Identify significant variants
sig_threshold <- 5e-8
sig_variants <- results[results$Score.pval < sig_threshold, ]

if(nrow(sig_variants) > 0) {
  cat(sprintf("\nGenome-wide significant variants (P < 5e-8): %d\n", nrow(sig_variants)))
  sig_variants_sorted <- sig_variants[order(sig_variants$Score.pval), ]
  print(head(sig_variants_sorted[, c("variant.id", "chr", "pos", "ref", "alt", 
                                     "n.obs", "freq", "Score.Stat", "Score.pval")], 20))
  fwrite(sig_variants_sorted, "significant_variants.txt", sep = "\t", quote = FALSE)
} else {
  cat("\nNo genome-wide significant variants found.\n")
}

# Suggestive variants
sug_threshold <- 1e-5
sug_variants <- results[results$Score.pval < sug_threshold, ]

if(nrow(sug_variants) > 0) {
  cat(sprintf("Suggestive variants (P < 1e-5): %d\n", nrow(sug_variants)))
  fwrite(sug_variants, "suggestive_variants.txt", sep = "\t", quote = FALSE)
}

# Close GDS file
seqClose(gds)

cat("\n=== GWAS Analysis Complete ===\n")
cat("Output files:\n")
cat("  - gwas_results.txt\n")
cat("  - manhattan_plot.png\n")
cat("  - qq_plot.png\n")
cat("  - pca_plot.png\n")
cat("  - null_model.rds\n")
if(nrow(sig_variants) > 0) cat("  - significant_variants.txt\n")
if(nrow(sug_variants) > 0) cat("  - suggestive_variants.txt\n")

