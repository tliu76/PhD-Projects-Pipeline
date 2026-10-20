library(RNAseqCNV)
#load("/ref_genomes/RNAseqCNV/human/GRch38/df_for_RNAseqCNV_hg38.rdata")

#file_config="/scratch/zuhu/project/GordanaRaca/ALL/out/COH002937_D1/TRANSCRIPTOME/RNAseqCNV/config"
#file_meta="/scratch/zuhu/project/GordanaRaca/ALL/out/COH002937_D1/TRANSCRIPTOME/RNAseqCNV/metadata"

args <- commandArgs(trailingOnly = TRUE)
print(args)

file_config=args[1]
file_meta=args[2]

#RNAseqCNV_wrapper(config = file_config, metadata = file_meta, snv_format = "vcf",referData = refData, keptSNP = keepSNP,mafRange=c(0.05,0.85))
  
RNAseqCNV_wrapper(config = file_config, metadata = file_meta, snv_format = "vcf",genome_version="hg38",mafRange=c(0.01,0.85))