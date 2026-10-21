# Path to the MIGI_1_2_MXB_Post_QC dataset
# /project2/desmith_128/datasets/For_Tanxin/MIGI_1_2_MXB_Post_QC
.libPaths("/home1/tanxinli/R_lib")
suppressMessages({
  library(ggplot2); library(dplyr); library(data.table); library(ggrepel)
})
covar <- read.table("/project2/desmith_128/datasets/For_Tanxin/MIGI_1_2_MXB_Post_QC/covar_sex_pca_fixed.txt", header = TRUE)
covar <- covar %>% select(-c("PHENO"))
covar <- covar %>% rename("#FID"="FID")
fwrite(covar, "/project2/desmith_128/datasets/For_Tanxin/MIGI_1_2_MXB_Post_QC/covar_20PCs_sex.txt")
covar2 <- covar %>% select(-c("PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20"))
fwrite(covar2, "/project2/desmith_128/datasets/For_Tanxin/MIGI_1_2_MXB_Post_QC/covar_10PCs_sex.txt")

covar2 <- fread( "/project2/desmith_128/datasets/For_Tanxin/MIGI_1_2_MXB_Post_QC/covar_10PCs_sex.txt")

### generate CHLA_CALSEC_MEC 10 PCs and sec covariate file

.libPaths("/home1/tanxinli/R_lib")
library(data.table)
library(dplyr)

covar <- fread("/project2/desmith_128/new_tx/chla_calsec_mec/2vcf/filter/merge/covar_20pca_chla_calsec_mec.txt")
covar2 <- covar %>% select(-c("PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20"))
fwrite(covar2, "/project2/desmith_128/new_tx/chla_calsec_mec/2vcf/filter/merge/covar_10pca_chla_calsec_mec.txt")




