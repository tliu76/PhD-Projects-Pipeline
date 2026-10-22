#Get installed list of pacakges
installedPackages=  installed.packages()[,"Package"]

#Install CRAN R packages
list.of.packages=c("devtools", "BiocManager","dplyr","stringr","Seurat",
                   "ggplot2","ggrepel","cowplot","umap",
                   "shiny","shinyjs","shinydashboard")

new.packages = list.of.packages[!(list.of.packages %in% installedPackages)]
if(length(new.packages)){
  install.packages(new.packages)
}

#Install Bioconductor R packages
list.of.bioconductor.packages = c("DESeq2", "SingleR","SummarizedExperiment")
new.bioconductor.packages = list.of.bioconductor.packages[!(list.of.bioconductor.packages %in% installedPackages)]

if (length(new.bioconductor.packages)) {
  BiocManager::install(new.bioconductor.packages)
}

#Install PhenoGrapth
if(!"Rphenograph" %in% installedPackages){
  devtools::install_github("JinmiaoChenLab/Rphenograph")
}

#Install MD-ALL
devtools::install_github("gu-lab20/MD-ALL")

library(MDALL)
library(cli)
install.packages("cli")
install.packages("Matrix")
remotes::install_version("Matrix", version = "1.6-5")
.libPaths()


Sys.getenv("LD_LIBRARY_PATH")

remotes::install_version("Matrix", version = "1.6-5")

.libPaths()


###############################################################################

library(MDALL)
library(data.table)
setwd("/scratch1/tanxinli/rnaseq_out/")
run_shiny_MDALL()

###comand line:#################################################################
#https://github.com/gu-lab20/MD-ALL

library(MDALL)
df_count=read_input("./HTSeq/200101_T.DUX4patched.HTSeq",delimiter = "\t",header = F)
names(df_count)
df_count <- df_count[grepl("^ENSG", df_count$feature), ]

df_vst=get_vst_values(obj_in = obj_234_HTSeq,df_count = df_count)



df_feateure_exp=get_geneExpression(df_vst = df_vst,genes = c("CDX2","CRLF2","NUTM1"))




###############################################################################