
.libPaths("/home1/tanxinli/R_lib")
suppressPackageStartupMessages({
  library(SeqArray)
  library(SNPRelate)
  library(GENESIS)
  library(Biobase)
  library(data.table)
})
##https://uw-gac.github.io/SISG_2019/pc-relate.html#king
##file:///Users/tanxinliu/Downloads/03_pop_structure_relatedness.html
# ---------------------------------------------------------------
# 1. Paths
# ---------------------------------------------------------------


gds_merged <- "/scratch1/tanxinli/mp2prt_gwas/filtered/split/MP2PRT.filtered.merged.gds"
pruned_dir <- "/scratch1/tanxinli/mp2prt_gwas/filtered/result/pcair/pruned_by_chr"
out_dir    <- "/scratch1/tanxinli/mp2prt_gwas/filtered/result/pcrelate_test10"

# ---------------------------------------------------------------
# 2. Merge LD-pruned variant lists
# ---------------------------------------------------------------
message("Merging LD-pruned SNPs from chr1–22...")

pruned_all <- c()
for (chr in 1:22) {
  f <- file.path(pruned_dir, paste0("ldpruned_chr", chr, ".RDS"))
  if (!file.exists(f)) next
  v <- readRDS(f)
  pruned_all <- c(pruned_all, v)
}

pruned <- unique(pruned_all)


gds <- seqOpen(gds_merged)
sample.id <- seqGetData(gds, "sample.id")

# KING
king <- snpgdsIBDKING(gds, sample.id=sample.id, snp.id=pruned)

# extract the kinship estimates
kingMat <- king$kinship
colnames(kingMat) <- rownames(kingMat) <- king$sample.id
dim(kingMat)


saveRDS(kingMat, file.path(out_dir, "kingMat.RDS"))

# ---------------------------------------------------------------
# 4. PC-AiR
# ---------------------------------------------------------------
library(GENESIS)
pca <- pcair(gds,
            kinobj = kingMat,
            kin.thresh=2^(-9/2),
            divobj = kingMat,
            div.thresh=-2^(-9/2), snp.include=pruned)
names(pca)

length(pca$unrels)
head(pca$unrels)
length(pca$rels)

# extract the top 10 PCs and make a data.frame
pcs <- data.frame(pca$vectors[,1:10])
colnames(pcs) <- paste0('PC', 1:10)
pcs$sample.id <- pca$sample.id
dim(pcs)

save(pcs, file="/scratch1/tanxinli/mp2prt_gwas/filtered/result/pcrelate_test10/pcs.RData")

pc.df <- as.data.frame(pca$vectors)
names(pc.df) <- paste0("PC", 1:ncol(pca$vectors))
pc.df$sample.id <- row.names(pca$vectors)
pc.df <- left_join(pc.df, pData(annot), by="sample.id")

#######################plot #########
library(Biobase)
annotfile <- "sample_annotation_1KG.rda"
repo_path <- "https://github.com/UW-GAC/GENESIS/raw/devel/data"
if (!file.exists(annotfile)) download.file(file.path(repo_path, annotfile), annotfile)
annot <- get(load(annotfile))
pc.df <- merge(pcs, pData(annot)[,c("sample.id", "Population")], by = "sample.id")



annot <- AnnotatedDataFrame(annot)
pData(annot)
pc.df <- merge(pcs, pData(annot)[,c("sample.id", "Population")], by = "sample.id")

library(GGally)
library(RColorBrewer)
pop.cols <- setNames(brewer.pal(12, "Paired"),
                 c("ACB", "ASW", "CEU", "GBR", "CHB", "JPT", "CLM", "MXL", "LWK", "YRI", "GIH", "PUR"))
ggparcoord(pc.df, columns=2:11, groupColumn="Population", scale="uniminmax") +
    scale_color_manual(values=pop.cols) +
    xlab("PC") + ylab("")
    
#######################PC relate #########
seqResetFilter(gds, verbose=FALSE)
library(SeqVarTools)
seqData <- SeqVarData(gds)

# filter the GDS object to our LD-pruned variants
seqSetFilter(seqData, variant.id=pruned)

iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)

pcrel <- pcrelate(iterator,  pcs=pca$vectors[,1:2],   training.set=pca$unrels)

names(pcrel)
# self-kinship estimates
dim(pcrel$kinSelf)
head(pcrel$kinSelf)
# save output
save(pcrel, file="/scratch1/tanxinli/mp2prt_gwas/filtered/result/pcrelate_test10/pcrelate.RData")

##### pc-relate file#####
library(ggplot2)
ggplot(pcrel$kinBtwn, aes(k0, kin)) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
    geom_hex(bins = 100) +
    geom_abline(intercept = 0.25, slope = -0.25) +
    ylab("kinship estimate") +
    theme_bw()
    
ggsave("/scratch1/tanxinli/mp2prt_gwas/filtered/result/pcrelate_test10/kinship_pc-relate2.png", 
       width = 6, height = 5, dpi = 300)
       
############# 6. Null model #############
covar <- readRDS("/scratch1/tanxinli/mp2prt_gwas/filtered/split/covar_full.RDS")
metadata <- data.frame(
  labelDescription = c(
    "sample.id", "Age", "Sex",  "n.all.del"  ,        "n.total.least.rag.del"  ,    "n.total.both.rag.del"  ,    
  "n.nonig.least.rag.del"    ,  "n.ig.least.rag.del" ,       
 "prop.total.least.rag.del"  , "prop.total.both.rag.del"   ,
 "prop.nonig.least.rag.del"  , "prop.ig.least.rag.del"   ,   "Inverse normal transformed RAG-del count"),
  row.names = names(covar)
)

annot <- AnnotatedDataFrame(covar, metadata)

sampleData(seqData) <- annot

save(seqData, file = "/scratch1/tanxinli/mp2prt_gwas/filtered/result/pcrelate_test10/phenotype_annotation.RData")
load("/scratch1/tanxinli/mp2prt_gwas/filtered/result/pcrelate_test10/phenotype_annotation.RData")
saveRDS(annot, file = "/scratch1/tanxinli/mp2prt_gwas/filtered/result/pcrelate_test10/annot.rds")
annot <- readRDS("/scratch1/tanxinli/mp2prt_gwas/filtered/result/pcrelate_test10/annot.rds")
grm <- pcrelateToMatrix(pcrel, scaleKin=2)


nullmod <- fitNullModel(seqData, 
                        outcome = "n.nonig.least.rag.del", 
                        covars=c("SEX", paste0("PC", 1:10)),
                        cov.mat=grm, norm.option="all",
                        family = "gaussian",
                        verbose=FALSE)
save(nullmod, file="/scratch1/tanxinli/mp2prt_gwas/filtered/result/pcrelate_test10/null_model_n.nonig.least.rag.del_INT2.RData")




                        
       