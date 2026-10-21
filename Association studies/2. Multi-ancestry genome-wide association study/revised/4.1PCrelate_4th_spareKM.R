#!/usr/bin/env Rscript


# ===============================================================
# 4.run_pcrelate_full.R
# Inputs:
#   - MP2PRT.filtered.merged.gds
#   - ldpruned_allchr.RDS
#   - pcair_PCs_top10.csv
#   - pcair_unrels.txt
# Outputs:
#   - pcrelate_obj.RDS
#   - pcrelate_kinship_mat_sparse_4th.RDS
# ===============================================================


.libPaths(c("/home1/tanxinli/R_lib", .libPaths()))
suppressPackageStartupMessages({
  library(SeqArray)
  library(SNPRelate)
  library(GENESIS)
  library(SeqVarTools)
  library(Biobase)
  library(data.table)
  library(BiocParallel)
})




message("==== START: PC-Relate (full genome, no chr split) ====")


## -----------------------------------------------------------
## 1. 路径
## -----------------------------------------------------------

## -----------------------------------------------------------
## 2. 打开 GDS & 读取样本
## -----------------------------------------------------------
load("/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0210/kingMat.RData")
load("/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0210/pruned.RData")


out_dir<- "/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0210/"

gds_merged <- "/project2/desmith_128/tanxin/scratch1_new//mp2prt_gwas/filtered/split/MP2PRT.filtered.merged.gds"
library(SeqArray)
gds <- seqOpen(gds_merged)
sample.id <- seqGetData(gds, "sample.id")
### 1491 samples ###

df <- fread( "/project2/desmith_128/tanxin/scratch1_new/project3_vdj/RAG_Recombination/Aim1-3/regression_df.csv")
table(df$new_subtype)

keep.samples = df$sample[df$new_subtype == "ETV6::RUNX1"]

seqSetFilter(gds, sample.id = keep.samples)
## of selected samples: 554

library(SeqVarTools)
seqData <- SeqVarData(gds)
seqSetFilter(gds, sample.id = keep.samples)


## -----------------------------------------------------------
## 3. 读取 PC-AiR PCs & unrels，并按 GDS 顺序对齐
## -----------------------------------------------------------
### 3.2 PC-AiR ###
### The next step is PC-AiR, in which we select a set of unrelated samples that is maximally informative about all ancestries in the sample. We use this unrelated set for Principal Component Analysis (PCA), then project the relatives onto the PCs.
load("/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0210/sampset.RData")

names(sampset)

sapply(sampset, length)


# project values for relatives

load("/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0210/pcs.RData")


### 3.3 PC-Relate ###

## The first 2 PCs separate populations, so we use them to compute kinship estimates adjusting for ancestry. The PC-Relate function expects a SeqVarIterator object.
seqSetFilter(gds, sample.id = keep.samples) #Applies a filter to the GDS object; Restricts analyses to a subset of samples (or variants, if you filter those); The data are not deleted — they’re just hidden from view

### seqResetFilter(gds, verbose=FALSE) remove filters

library(SeqVarTools)
seqData <- SeqVarData(gds)

seqSetFilter(seqData, sample.id=keep.samples, variant.id=pruned)
# of selected samples: 554
# of selected variants: 369,356

iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)

pcrel <- pcrelate(iterator, pcs=pcs[,1:2], training.set=sampset$unrels, 
                  sample.include=keep.samples)
names(pcrel)

save(pcrel, file="/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0210/pcrelate_kinship.RData")



## -----------------------------------------------------------
## 8. 转换为 kinship matrix（给 fitNullModel 用）
## -----------------------------------------------------------
message("Converting PC-Relate output to kinship matrix (pcrelateToMatrix)...")

kinMat1 <- pcrelateToMatrix(pcrel)

kin_out <- file.path(out_dir, "pcrelate_kinship_mat.RDS")
saveRDS(kinMat1, kin_out)


## -----------------------------------------------------------
## 9. 收尾
## -----------------------------------------------------------

seqClose(gds)

