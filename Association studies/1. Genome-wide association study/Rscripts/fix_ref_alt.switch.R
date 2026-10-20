#!/mnt/linuxtools/R/R-4.0.0/bin/R

# pzhang13@dell740-statgen> pwd
# /mnt/research/active/S_Beitelshees_CanagliflozinResponse_GDA_040723_1/DataCleaning_ivy/Imputation
# pzhang13@dell740-statgen> Rscript Imputation_summary_plots.R \
# -i /mnt/research/active/S_Beitelshees_CanagliflozinResponse_GDA_040723_1/DataCleaning_ivy/Imputation/Results \
# -s 1 -e 23  >print.log.txt
# script to run ExomeDepth
# how to run:
# runExomeDepth.r -b bamlist.txt -o /output/path/ -v > runExomeDepthOutput.txt

start_time = Sys.time()
.libPaths("/home1/tanxinli/R_lib")

install.packages("optparse")

require(optparse) # http://cran.r-project.org/web/packages/optparse/optparse.pdf
options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings
options(bitmapType="cairo")


option_list = list(
  make_option(c("-b", "--bim"), action="store", default='', 
              type='character', help="PLINK bim file"),
  # make_option(c("-s", "--startchr"), action="store", default=1L, 
  #              type='integer', help="bed file"),
  # make_option(c("-e", "--endchr"), action="store", default=23L, 
  #              type='integer', help="sample name")
  make_option(c("-o", "--outdir"), action="store", default='./',
              type='character', help="output directory [default %default]")
)
opt = parse_args(OptionParser(option_list=option_list))
print(opt,file=stderr())

plink_bim = opt$bim
out_dir = opt$outdir

####function to summarize and plot only the non-monomorphic site
library(Rsamtools)
library(BSgenome)
path <- "/scratch1/tanxinli/GUAT/GRECERE_GWAS_data/combined/triple"
# outPrefix <- "beitelshees"
bim <- read.table("/scratch1/tanxinli/GUAT/GRECERE_GWAS_data/combined/triple/all_geno_maf_hwe_mind.lifted.hg38.final.bim", header=F, stringsAsFactors=F)

#setnames(bim,c("chr","snp","gpos","pos","A1","A2"))
colnames(bim) <- c("chr", "snp", "gpos", "pos", "A1", "A2")
bim$chr <- ifelse(bim$chr == 23|bim$chr == 25, "chrX", paste0("chr", bim$chr))
# bim[, chr:= ifelse(chr == 23,"chrX",paste0("chr",chr))]

# head(bim)
#   chr                   rsID gen    pos ref alt
# 1   1            1:49554-G-A   0  49554   G   A
# 2   1           JHU_1.534237   0 598858   C   A
# 3   1              rs9701055   0 630053   T   C
uniqueChr <- unique(bim$chr)
# i <- 21
ref <- NULL
for (i in 1: length(uniqueChr)) {
  # for (i in 1:23) {
  chrom <- uniqueChr[i]
  message(chrom)
  dat <- bim[bim$chr == chrom,]
  faFileName <- paste0("/project/desmith_488/tanxinli/gwasref/GRCh38.",chrom,".fa")
  fasta_file <- FaFile(file=faFileName)
  gr1 <- GRanges(dat$chr,IRanges(start=as.numeric(dat$pos), end=as.numeric(dat$pos)))
  # gr1
  # GRanges object with 22555 ranges and 0 metadata columns:
  #        seqnames    ranges strand
  #           <Rle> <IRanges>  <Rle>
  #   [1]       21   8537605      *
  # -------
  #   seqinfo: 1 sequence from an unspecified genome; no seqlengths
  refbase <- getSeq(fasta_file, gr1)
  refbase <- as.data.frame(refbase)$x
  dat$REF_hg38 <- refbase
  if (i ==1) {
    ref <- dat
    # write.table(dat, paste0(path,outPrefix,".refAlleles"), quote=F,  na = "",sep="\t", col.names=T, row.names=F )
  } else{
    ref <- rbind(ref, dat)
    # write.table(dat, paste0(path,outPrefix,".refAlleles"), append=T, quote=F,  na = "",sep="\t", col.names=F, row.names=F )
  }
}
write.table(ref,"/scratch1/tanxinli/GUAT/GRECERE_GWAS_data/combined/triple/chr1_22_refAlleles", quote=F,  na = "",sep="\t", col.names=T, row.names=F )
#write.table(ref, paste0(out_dir,"chr1_23_refAlleles"), quote=F,  na = "",sep="\t", col.names=T, row.names=F )
# alleles <- read.table("beitelshees.refAlleles", header = T)
# head(alleles)
# chr                    snp gpos    pos A1 A2 REF_hg38
# 1 chr1            1:49554-G-A    0  49554  G  A        A
# 2 chr1           JHU_1.534237    0 598858  C  A        A
# 3 chr1              rs9701055    0 630053  T  C        C
ref_a1 <- ref[ref$REF_hg38 == ref$A1, ] #142868 variants
switch_a1 <- ref_a1[, c("snp", "REF_hg38")]
write.table(switch_a1, "/scratch1/tanxinli/GUAT/GRECERE_GWAS_data/combined/triple/file_contains_ref.alt_switch_SNP", quote=FALSE, row.names=FALSE, col.names = FALSE)

#write.table(switch_a1, paste0(out_dir,"file_contains_ref.alt_switch_SNP"), quote=FALSE, row.names=FALSE, col.names = FALSE)

message("\n file_contains_ref.alt_switch_SNP is generated!")
