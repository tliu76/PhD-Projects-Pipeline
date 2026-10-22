library(VennDiagram)
library(data.table)
library(tidyverse)
library(pals)
library(maftools)


setwd("/scratch2/tanxinli/wgs.smoking.out/survivor.out/AnnotSV_out_two_caller")

# read file path  
all_paths <- list.files(pattern = "*.annotated.tsv",full.names = TRUE)

# read file content
all_content <-
  all_paths %>%
  lapply(fread,
         header = TRUE, colClasses="character")

# read file name
all_filenames <- all_paths %>%
  as.list()

# combine file content list and file name list
all_lists <- mapply(c, all_content, all_filenames, SIMPLIFY = FALSE)

# unlist all lists and change column name
all_result <- rbindlist(all_lists, fill = T)

# check SV counts across chroms for each sample use MODE =FULL 
# MODE = FULL have the full records even those without gene anno
all_result %>% filter(Annotation_mode == "full") %>%separate(V1, c("sample_name", "junk"), sep = ".merged") %>% separate(sample_name, c("junk", "sample_name"), sep = "./") %>% count(sample_name, SV_chrom) %>% ggplot(aes(x = factor(SV_chrom, c("1", "2", "3", "4", "5", "6", "7",  "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")), y = n, group = sample_name,color = sample_name)) + geom_line() + theme(legend.position = "none", axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) + scale_color_manual(values=as.vector(tol.rainbow(38))) + guides(color = guide_legend(title = "")) + labs(x = "chromosome") + scale_y_continuous(expand = c(0, 0))

# check SV_length

# mandatory columns: Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Classification, Variant_Type and Tumor_Sample_Barcode.

df_SV_MAF <- all_result %>% separate(V1, c("sample_name", "junk"), sep = ".merged") %>% separate(sample_name, c("junk", "sample_name"), sep = "./") %>% mutate(Tumor_Sample_Barcode = paste0("T", sample_name)) %>% filter(Annotation_mode == "split") %>% select(c("Tumor_Sample_Barcode","SV_chrom","SV_start","SV_end", "SV_type", "REF", "ALT", "Gene_name") )%>% mutate(Variant_Type = "SV") %>% dplyr::rename("Hugo_Symbol" = "Gene_name", "Chromosome" = "SV_chrom", "Start_Position" = "SV_start", "End_Position" = "SV_end", "Reference_Allele" = "REF", "Tumor_Seq_Allele2" = "ALT", "Variant_Classification" = "SV_type")

df_SV_MAF %>% fwrite("smoking_maf_SV.maf", sep = "\t")

SV_MAF = read.maf(maf = "./smoking_maf_SV.maf", removeDuplicatedVariants = F, vc_nonSyn=c("DEL", "DUP", "TRA", "INV", "INS"))

SV_MAF_DEL = read.maf(maf = "./smoking_maf_SV.maf", removeDuplicatedVariants = F, vc_nonSyn=c("DEL"))
SV_MAF_DUP = read.maf(maf = "./smoking_maf_SV.maf", removeDuplicatedVariants = F, vc_nonSyn=c("DUP"))
SV_MAF_TRA = read.maf(maf = "./smoking_maf_SV.maf", removeDuplicatedVariants = F, vc_nonSyn=c("TRA"))
SV_MAF_INV = read.maf(maf = "./smoking_maf_SV.maf", removeDuplicatedVariants = F, vc_nonSyn=c("INV"))
getGeneSummary(SV_MAF)

# check subtypes 
# oncoplot - important genes

genes_important <- fread("/project/desmith_488/tanxinli/code/gene/genes.subtypes/Studd.master_driver_gene_list.txt") %>% pull(Gene)

oncoplot(maf = SV_MAF, showTumorSampleBarcodes = T, fontSize = 0.5, genes = genes_important, gene_mar = 6, barcode_mar = 6, legend_height = 2, legendFontSize = 0.8, titleFontSize = 1)

oncoplot(maf = SV_MAF_DEL, showTumorSampleBarcodes = T, fontSize = 0.5, genes = genes_important, gene_mar = 6, barcode_mar = 6, legend_height = 2, legendFontSize = 0.8, titleFontSize = 1)
oncoplot(maf = SV_MAF_DUP, showTumorSampleBarcodes = T, fontSize = 0.5, genes = genes_important, gene_mar = 6, barcode_mar = 6, legend_height = 2, legendFontSize = 0.8, titleFontSize = 1)
oncoplot(maf = SV_MAF_TRA, showTumorSampleBarcodes = T, fontSize = 0.5, genes = genes_important, gene_mar = 6, barcode_mar = 6, legend_height = 2, legendFontSize = 0.8, titleFontSize = 1)
oncoplot(maf = SV_MAF_INV, showTumorSampleBarcodes = T, fontSize = 0.5, genes = genes_important, gene_mar = 6, barcode_mar = 6, legend_height = 2, legendFontSize = 0.8, titleFontSize = 1)

# identify subtypes from Gu 
Gu.master_gene_list <- fread("/project/desmith_488/tanxinli/code/gene/genes.subtypes/Gu.master_gene_list.txt") 
Gu.master_gene_list.gene1 <- Gu.master_gene_list  %>% pull(Gene1) %>% unique()
Gu.master_gene_list.gene2 <- Gu.master_gene_list  %>% pull(Gene2) %>% unique()

## check subtypes
Gu.master_gene_list.AnnotSV_ID <- all_result %>% filter(Annotation_mode == "split") %>%separate(V1, c("sample_name", "junk"), sep = ".merged") %>% separate(sample_name, c("junk", "sample_name"), sep = "./") %>% select(-"junk") %>% filter(Gene_name %in% c(Gu.master_gene_list.gene1, Gu.master_gene_list.gene2)) 


# # Gu.master_gene_list.AnnotSV_ID %>% select(
# #   c("sample_name", "AnnotSV_ID","SV_chrom", "SV_start", "SV_end", "SV_length",
# #     "SV_type", "QUAL", "INFO", "Gene_name", "Tx", "Tx_start", "Tx_end", "Overlapped_tx_length",
# #     "Overlapped_CDS_length", "Overlapped_CDS_percent", "Frameshift",
# #     "Exon_count", "Location", "Location2", "Dist_nearest_SS", "Nearest_SS_type",
# #     "Intersect_start", "Intersect_end", "RE_gene", "P_gain_phen",
# #     "P_gain_hpo", "P_gain_source", "P_gain_coord", "P_loss_phen",
# #     "P_loss_hpo", "P_loss_source", "P_loss_coord", "P_ins_phen",
# #     "P_ins_hpo", "P_ins_source", "P_ins_coord", "P_snvindel_nb",
# #     "P_snvindel_phen", "B_gain_source", "B_gain_coord", "B_loss_source",
# #     "B_loss_coord", "B_ins_source", "B_ins_coord", "B_inv_source",
# #     "B_inv_coord", "TAD_coordinate", "ENCODE_experiment", "GC_content_left",
# #     "GC_content_right", "Repeat_coord_left", "Repeat_type_left",
# #     "Repeat_coord_right", "Repeat_type_right", "Gap_left", "Gap_right",
# #     "SegDup_left", "SegDup_right", "ENCODE_blacklist_left", "ENCODE_blacklist_characteristics_left",
# #     "ENCODE_blacklist_right", "ENCODE_blacklist_characteristics_right",
# #     "ACMG", "HI", "TS", "DDD_HI_percent", "DDD_status", "DDD_mode",
# #     "DDD_consequence", "DDD_disease", "DDD_pmid", "ExAC_delZ", "ExAC_dupZ",
# #     "ExAC_cnvZ", "ExAC_synZ", "ExAC_misZ", "OMIM_ID", "OMIM_phenotype",
# #     "OMIM_inheritance", "OMIM_morbid", "OMIM_morbid_candidate", "LOEUF_bin",
# #     "GnomAD_pLI", "ExAC_pLI", "Exomiser_gene_pheno_score", "Human_pheno_evidence",
# #     "Mouse_pheno_evidence", "Fish_pheno_evidence", "AnnotSV_ranking_score",
# #     "AnnotSV_ranking_criteria", "ACMG_class"
# #   )
# # ) %>% fwrite("Gu.master_gene_list.AnnotSV_ID.split.csv")


 Gu.master_gene_list.AnnotSV_ID %>% select(
   c("sample_name", "AnnotSV_ID","SV_chrom", "SV_start", "SV_end", "SV_length",
     "SV_type", "QUAL", "INFO", "Gene_name", "Tx", "Tx_start", "Tx_end", "Overlapped_tx_length",
    "Overlapped_CDS_length", "Overlapped_CDS_percent", "Frameshift",
    "Exon_count", "Location", "Location2", "Dist_nearest_SS", "Nearest_SS_type",
     "Intersect_start", "Intersect_end", "RE_gene", "P_gain_phen",
     "P_gain_hpo", "P_gain_source", "P_gain_coord", "P_loss_phen",
     "P_loss_hpo", "P_loss_source", "P_loss_coord", "P_ins_phen",
     "P_ins_hpo", "P_ins_source", "P_ins_coord", "P_snvindel_nb",
     "P_snvindel_phen", "B_gain_source", "B_gain_coord", "B_loss_source",
     "B_loss_coord", "B_ins_source", "B_ins_coord", "B_inv_source",
     "B_inv_coord", "TAD_coordinate", "ENCODE_experiment", "GC_content_left",
    "GC_content_right", "Repeat_coord_left", "Repeat_type_left",
     "Repeat_coord_right", "Repeat_type_right", "Gap_left", "Gap_right",
     "SegDup_left", "SegDup_right", "ENCODE_blacklist_left", "ENCODE_blacklist_characteristics_left",
     "ENCODE_blacklist_right", "ENCODE_blacklist_characteristics_right",
     "ACMG", "HI", "TS","DDD_HI_percent","ExAC_delZ","ExAC_dupZ","ExAC_cnvZ","ExAC_synZ",
    "ExAC_misZ","GenCC_disease","GenCC_moi","GenCC_classification","GenCC_pmid","OMIM_ID","OMIM_phenotype",
    "OMIM_inheritance","OMIM_morbid","OMIM_morbid_candidate","LOEUF_bin","GnomAD_pLI","ExAC_pLI",
    "AnnotSV_ranking_score","AnnotSV_ranking_criteria","ACMG_class"
   )
 ) %>% fwrite("Gu.master_gene_list.AnnotSV_ID.split.csv")


write.csv(Gu.master_gene_list.AnnotSV_ID, "/scratch2/tanxinli/wgs.smoking.out/survivor.out/AnnotSV_out_two_caller/Gu.master_gene_list.AnnotSV_ID",row.names=FALSE)


Gu.master_gene_list.AnnotSV_ID.vector <- Gu.master_gene_list.AnnotSV_ID  %>% distinct(AnnotSV_ID) %>% pull()

Gu.master_gene_list.AnnotSV_ID.full <- all_result %>% filter(Annotation_mode == "full") %>% separate(V1, c("sample_name", "junk"), sep = ".merged") %>% separate(sample_name, c("junk", "sample_name"), sep = "./") %>% select(-"junk") %>% filter(AnnotSV_ID %in% Gu.master_gene_list.AnnotSV_ID.vector) 

Gu.master_gene_list.AnnotSV_ID.full %>% select(
  c("sample_name", "AnnotSV_ID","SV_chrom", "SV_start", "SV_end", "SV_length", 
    "SV_type", "QUAL", "INFO", "Gene_name", "Tx", "Tx_start", "Tx_end", "Overlapped_tx_length", 
    "Overlapped_CDS_length", "Overlapped_CDS_percent", "Frameshift", 
    "Exon_count", "Location", "Location2", "Dist_nearest_SS", "Nearest_SS_type", 
    "Intersect_start", "Intersect_end", "RE_gene", "P_gain_phen", 
    "P_gain_hpo", "P_gain_source", "P_gain_coord", "P_loss_phen", 
    "P_loss_hpo", "P_loss_source", "P_loss_coord", "P_ins_phen", 
    "P_ins_hpo", "P_ins_source", "P_ins_coord", "P_snvindel_nb", 
    "P_snvindel_phen", "B_gain_source", "B_gain_coord", "B_loss_source", 
    "B_loss_coord", "B_ins_source", "B_ins_coord", "B_inv_source", 
    "B_inv_coord", "TAD_coordinate", "ENCODE_experiment", "GC_content_left", 
    "GC_content_right", "Repeat_coord_left", "Repeat_type_left", 
    "Repeat_coord_right", "Repeat_type_right", "Gap_left", "Gap_right", 
    "SegDup_left", "SegDup_right", "ENCODE_blacklist_left", "ENCODE_blacklist_characteristics_left", 
    "ENCODE_blacklist_right", "ENCODE_blacklist_characteristics_right", 
    "ACMG", "HI", "TS", "DDD_HI_percent","ExAC_delZ","ExAC_dupZ","ExAC_cnvZ","ExAC_synZ",
    "ExAC_misZ","GenCC_disease","GenCC_moi","GenCC_classification","GenCC_pmid","OMIM_ID","OMIM_phenotype",
    "OMIM_inheritance","OMIM_morbid","OMIM_morbid_candidate","LOEUF_bin","GnomAD_pLI","ExAC_pLI",
    "AnnotSV_ranking_score","AnnotSV_ranking_criteria","ACMG_class"
  )
) %>% fwrite("Gu.master_gene_list.AnnotSV_ID.full.csv")



# Gu.master_gene_list.AnnotSV_ID.full %>% select(
#   c("sample_name", "AnnotSV_ID","SV_chrom", "SV_start", "SV_end", "SV_length", 
#     "SV_type", "QUAL", "INFO", "Gene_name", "Tx", "Tx_start", "Tx_end", "Overlapped_tx_length", 
#     "Overlapped_CDS_length", "Overlapped_CDS_percent", "Frameshift", 
#     "Exon_count", "Location", "Location2", "Dist_nearest_SS", "Nearest_SS_type", 
#     "Intersect_start", "Intersect_end", "RE_gene", "P_gain_phen", 
#     "P_gain_hpo", "P_gain_source", "P_gain_coord", "P_loss_phen", 
#     "P_loss_hpo", "P_loss_source", "P_loss_coord", "P_ins_phen", 
#     "P_ins_hpo", "P_ins_source", "P_ins_coord", "P_snvindel_nb", 
#     "P_snvindel_phen", "B_gain_source", "B_gain_coord", "B_loss_source", 
#     "B_loss_coord", "B_ins_source", "B_ins_coord", "B_inv_source", 
#     "B_inv_coord", "TAD_coordinate", "ENCODE_experiment", "GC_content_left", 
#     "GC_content_right", "Repeat_coord_left", "Repeat_type_left", 
#     "Repeat_coord_right", "Repeat_type_right", "Gap_left", "Gap_right", 
#     "SegDup_left", "SegDup_right", "ENCODE_blacklist_left", "ENCODE_blacklist_characteristics_left", 
#     "ENCODE_blacklist_right", "ENCODE_blacklist_characteristics_right", 
#     "ACMG", "HI", "TS", "DDD_HI_percent", "DDD_status", "DDD_mode", 
#     "DDD_consequence", "DDD_disease", "DDD_pmid", "ExAC_delZ", "ExAC_dupZ", 
#     "ExAC_cnvZ", "ExAC_synZ", "ExAC_misZ", "OMIM_ID", "OMIM_phenotype", 
#     "OMIM_inheritance", "OMIM_morbid", "OMIM_morbid_candidate", "LOEUF_bin", 
#     "GnomAD_pLI", "ExAC_pLI", "Exomiser_gene_pheno_score", "Human_pheno_evidence", 
#     "Mouse_pheno_evidence", "Fish_pheno_evidence", "AnnotSV_ranking_score", 
#     "AnnotSV_ranking_criteria", "ACMG_class"
#   )
# ) %>% fwrite("Gu.master_gene_list.AnnotSV_ID.full.csv")


# identify subtypes from Roberts 
Roberts.master_gene_list <- fread("/project/desmith_488/tanxinli/code/gene/genes.subtypes/Roberts.master_gene_list.txt")

Roberts.master_gene_list.gene1 <- Roberts.master_gene_list %>% pull(Gene1) %>% unique()
Roberts.master_gene_list.gene2 <- Roberts.master_gene_list %>% pull(Gene2) %>% unique()

Roberts.master_gene_list.AnnotSV_ID <- all_result %>% filter(Annotation_mode == "split") %>% separate(V1, c("sample_name", "junk"), sep = ".merged") %>% separate(sample_name, c("junk", "sample_name"), sep = "./") %>% select(-"junk") %>% filter(Gene_name %in% c(Roberts.master_gene_list.gene1, Roberts.master_gene_list.gene2)) 
Roberts.master_gene_list.AnnotSV_ID %>% select(
  c("sample_name", "AnnotSV_ID","SV_chrom", "SV_start", "SV_end", "SV_length", 
    "SV_type", "QUAL", "INFO", "Gene_name", "Tx", "Tx_start", "Tx_end", "Overlapped_tx_length", 
    "Overlapped_CDS_length", "Overlapped_CDS_percent", "Frameshift", 
    "Exon_count", "Location", "Location2", "Dist_nearest_SS", "Nearest_SS_type", 
    "Intersect_start", "Intersect_end", "RE_gene", "P_gain_phen", 
    "P_gain_hpo", "P_gain_source", "P_gain_coord", "P_loss_phen", 
    "P_loss_hpo", "P_loss_source", "P_loss_coord", "P_ins_phen", 
    "P_ins_hpo", "P_ins_source", "P_ins_coord", "P_snvindel_nb", 
    "P_snvindel_phen", "B_gain_source", "B_gain_coord", "B_loss_source", 
    "B_loss_coord", "B_ins_source", "B_ins_coord", "B_inv_source", 
    "B_inv_coord", "TAD_coordinate", "ENCODE_experiment", "GC_content_left", 
    "GC_content_right", "Repeat_coord_left", "Repeat_type_left", 
    "Repeat_coord_right", "Repeat_type_right", "Gap_left", "Gap_right", 
    "SegDup_left", "SegDup_right", "ENCODE_blacklist_left", "ENCODE_blacklist_characteristics_left", 
    "ENCODE_blacklist_right", "ENCODE_blacklist_characteristics_right", 
    "ACMG", "HI", "TS", "DDD_HI_percent","ExAC_delZ","ExAC_dupZ","ExAC_cnvZ","ExAC_synZ",
    "ExAC_misZ","GenCC_disease","GenCC_moi","GenCC_classification","GenCC_pmid","OMIM_ID","OMIM_phenotype",
    "OMIM_inheritance","OMIM_morbid","OMIM_morbid_candidate","LOEUF_bin","GnomAD_pLI","ExAC_pLI",
    "AnnotSV_ranking_score","AnnotSV_ranking_criteria","ACMG_class"
  )
) %>% fwrite("Roberts.master_gene_list.AnnotSV_ID.split.csv")

# Roberts.master_gene_list.AnnotSV_ID %>% select(
#   c("sample_name", "AnnotSV_ID","SV_chrom", "SV_start", "SV_end", "SV_length", 
#     "SV_type", "QUAL", "INFO", "Gene_name", "Tx", "Tx_start", "Tx_end", "Overlapped_tx_length", 
#     "Overlapped_CDS_length", "Overlapped_CDS_percent", "Frameshift", 
#     "Exon_count", "Location", "Location2", "Dist_nearest_SS", "Nearest_SS_type", 
#     "Intersect_start", "Intersect_end", "RE_gene", "P_gain_phen", 
#     "P_gain_hpo", "P_gain_source", "P_gain_coord", "P_loss_phen", 
#     "P_loss_hpo", "P_loss_source", "P_loss_coord", "P_ins_phen", 
#     "P_ins_hpo", "P_ins_source", "P_ins_coord", "P_snvindel_nb", 
#     "P_snvindel_phen", "B_gain_source", "B_gain_coord", "B_loss_source", 
#     "B_loss_coord", "B_ins_source", "B_ins_coord", "B_inv_source", 
#     "B_inv_coord", "TAD_coordinate", "ENCODE_experiment", "GC_content_left", 
#     "GC_content_right", "Repeat_coord_left", "Repeat_type_left", 
#     "Repeat_coord_right", "Repeat_type_right", "Gap_left", "Gap_right", 
#     "SegDup_left", "SegDup_right", "ENCODE_blacklist_left", "ENCODE_blacklist_characteristics_left", 
#     "ENCODE_blacklist_right", "ENCODE_blacklist_characteristics_right", 
#     "ACMG", "HI", "TS", "DDD_HI_percent", "DDD_status", "DDD_mode", 
#     "DDD_consequence", "DDD_disease", "DDD_pmid", "ExAC_delZ", "ExAC_dupZ", 
#     "ExAC_cnvZ", "ExAC_synZ", "ExAC_misZ", "OMIM_ID", "OMIM_phenotype", 
#     "OMIM_inheritance", "OMIM_morbid", "OMIM_morbid_candidate", "LOEUF_bin", 
#     "GnomAD_pLI", "ExAC_pLI", "Exomiser_gene_pheno_score", "Human_pheno_evidence", 
#     "Mouse_pheno_evidence", "Fish_pheno_evidence", "AnnotSV_ranking_score", 
#     "AnnotSV_ranking_criteria", "ACMG_class"
#   )
# ) %>% fwrite("Roberts.master_gene_list.AnnotSV_ID.split.csv")

Roberts.master_gene_list.AnnotSV_ID.vector <- Roberts.master_gene_list.AnnotSV_ID  %>% distinct(AnnotSV_ID) %>% pull()

Roberts.master_gene_list.AnnotSV_ID.full <- all_result %>% filter(Annotation_mode == "full") %>% separate(V1, c("sample_name", "junk"), sep = ".merged") %>% separate(sample_name, c("junk", "sample_name"), sep = "./") %>% select(-"junk") %>% filter(AnnotSV_ID %in% Roberts.master_gene_list.AnnotSV_ID.vector) 

Roberts.master_gene_list.AnnotSV_ID.full %>% select(
  c("sample_name", "AnnotSV_ID","SV_chrom", "SV_start", "SV_end", "SV_length", 
    "SV_type", "QUAL", "INFO", "Gene_name", "Tx", "Tx_start", "Tx_end", "Overlapped_tx_length", 
    "Overlapped_CDS_length", "Overlapped_CDS_percent", "Frameshift", 
    "Exon_count", "Location", "Location2", "Dist_nearest_SS", "Nearest_SS_type", 
    "Intersect_start", "Intersect_end", "RE_gene", "P_gain_phen", 
    "P_gain_hpo", "P_gain_source", "P_gain_coord", "P_loss_phen", 
    "P_loss_hpo", "P_loss_source", "P_loss_coord", "P_ins_phen", 
    "P_ins_hpo", "P_ins_source", "P_ins_coord", "P_snvindel_nb", 
    "P_snvindel_phen", "B_gain_source", "B_gain_coord", "B_loss_source", 
    "B_loss_coord", "B_ins_source", "B_ins_coord", "B_inv_source", 
    "B_inv_coord", "TAD_coordinate", "ENCODE_experiment", "GC_content_left", 
    "GC_content_right", "Repeat_coord_left", "Repeat_type_left", 
    "Repeat_coord_right", "Repeat_type_right", "Gap_left", "Gap_right", 
    "SegDup_left", "SegDup_right", "ENCODE_blacklist_left", "ENCODE_blacklist_characteristics_left", 
    "ENCODE_blacklist_right", "ENCODE_blacklist_characteristics_right", 
    "ACMG", "HI", "TS",  "DDD_HI_percent","ExAC_delZ","ExAC_dupZ","ExAC_cnvZ","ExAC_synZ",
    "ExAC_misZ","GenCC_disease","GenCC_moi","GenCC_classification","GenCC_pmid","OMIM_ID","OMIM_phenotype",
    "OMIM_inheritance","OMIM_morbid","OMIM_morbid_candidate","LOEUF_bin","GnomAD_pLI","ExAC_pLI",
    "AnnotSV_ranking_score","AnnotSV_ranking_criteria","ACMG_class"
  )
) %>% fwrite("Roberts.master_gene_list.AnnotSV_ID.full.csv")

# merge driver genes, gu genes and Roberts genes togethor  
ALL_driver_genes <- fread("/project/desmith_488/tanxinli/code/gene/genes.subtypes/Studd.master_driver_gene_list.txt") %>% pull(Gene)
Gu.master_gene_list <- fread("/project/desmith_488/tanxinli/code/gene/genes.subtypes/Gu.master_gene_list.txt") 
Gu.master_gene_list.gene1 <- Gu.master_gene_list  %>% pull(Gene1) %>% unique()
Gu.master_gene_list.gene2 <- Gu.master_gene_list  %>% pull(Gene2) %>% unique()
Roberts.master_gene_list <- fread("/project/desmith_488/tanxinli/code/gene/genes.subtypes/Roberts.master_gene_list.txt")
Roberts.master_gene_list.gene1 <- Roberts.master_gene_list %>% pull(Gene1) %>% unique()
Roberts.master_gene_list.gene2 <- Roberts.master_gene_list %>% pull(Gene2) %>% unique()

ALL_master_genes <- c(ALL_driver_genes, Gu.master_gene_list.gene1, Gu.master_gene_list.gene2, Roberts.master_gene_list.gene1, Roberts.master_gene_list.gene2) %>% unique()

# data.table(candidate_gene = ALL_master_genes) %>% filter(candidate_gene!="X") %>% fwrite("ALL_master_genes.txt")

ALL.master_gene_list.AnnotSV_ID <- all_result %>% filter(Annotation_mode == "split") %>% separate(V1, c("sample_name", "junk"), sep = ".merged") %>% separate(sample_name, c("junk", "sample_name"), sep = "./") %>% select(-"junk") %>% filter(Gene_name %in% ALL_master_genes) 

ALL.master_gene_list.AnnotSV_ID %>% select(
  c("sample_name", "AnnotSV_ID","ID","SV_chrom", "SV_start", "SV_end", "SV_length", 
    "SV_type", "QUAL", "INFO", "Gene_name", "Tx", "Tx_start", "Tx_end", "Overlapped_tx_length", 
    "Overlapped_CDS_length", "Overlapped_CDS_percent", "Frameshift", 
    "Exon_count", "Location", "Location2", "Dist_nearest_SS", "Nearest_SS_type", 
    "Intersect_start", "Intersect_end", "RE_gene", "P_gain_phen", 
    "P_gain_hpo", "P_gain_source", "P_gain_coord", "P_loss_phen", 
    "P_loss_hpo", "P_loss_source", "P_loss_coord", "P_ins_phen", 
    "P_ins_hpo", "P_ins_source", "P_ins_coord", "P_snvindel_nb", 
    "P_snvindel_phen", "B_gain_source", "B_gain_coord", "B_loss_source", 
    "B_loss_coord", "B_ins_source", "B_ins_coord", "B_inv_source", 
    "B_inv_coord", "TAD_coordinate", "ENCODE_experiment", "GC_content_left", 
    "GC_content_right", "Repeat_coord_left", "Repeat_type_left", 
    "Repeat_coord_right", "Repeat_type_right", "Gap_left", "Gap_right", 
    "SegDup_left", "SegDup_right", "ENCODE_blacklist_left", "ENCODE_blacklist_characteristics_left", 
    "ENCODE_blacklist_right", "ENCODE_blacklist_characteristics_right", 
    "ACMG", "HI", "TS",  "DDD_HI_percent","ExAC_delZ","ExAC_dupZ","ExAC_cnvZ","ExAC_synZ",
    "ExAC_misZ","GenCC_disease","GenCC_moi","GenCC_classification","GenCC_pmid","OMIM_ID","OMIM_phenotype",
    "OMIM_inheritance","OMIM_morbid","OMIM_morbid_candidate","LOEUF_bin","GnomAD_pLI","ExAC_pLI",
    "AnnotSV_ranking_score","AnnotSV_ranking_criteria","ACMG_class"
  )
) %>% fwrite("ALL.master_gene_list.AnnotSV_ID.split.csv")

ALL.master_gene_list.AnnotSV_ID.sample_name_AnnotSV_ID <- ALL.master_gene_list.AnnotSV_ID %>% select(sample_name, AnnotSV_ID) %>% distinct() 

all_result.full <- all_result %>% filter(Annotation_mode == "full") %>% separate(V1, c("sample_name", "junk"), sep = ".merged") %>% separate(sample_name, c("junk", "sample_name"), sep = "./") %>% select(-"junk") 

ALL.master_gene_list.AnnotSV_ID.full <- ALL.master_gene_list.AnnotSV_ID.sample_name_AnnotSV_ID %>% left_join(all_result.full, by = c("sample_name", "AnnotSV_ID")) 

ALL.master_gene_list.AnnotSV_ID.full %>% filter(AnnotSV_ID %in% c("13_48304047_48307278_DEL_1", "9_21971206_21974676_DEL_1","X_1216185_1536603_DEL_1")) %>% fwrite("ALL.master_gene_list.AnnotSV_ID.full.duplicates.csv")

ALL.master_gene_list.AnnotSV_ID.full %>% select(
  c("sample_name", "AnnotSV_ID", "ID", "SV_chrom", "SV_start", "SV_end", "SV_length", 
    "SV_type", "QUAL", "INFO", "Gene_name", "Tx", "Tx_start", "Tx_end", "Overlapped_tx_length", 
    "Overlapped_CDS_length", "Overlapped_CDS_percent", "Frameshift", 
    "Exon_count", "Location", "Location2", "Dist_nearest_SS", "Nearest_SS_type", 
    "Intersect_start", "Intersect_end", "RE_gene", "P_gain_phen", 
    "P_gain_hpo", "P_gain_source", "P_gain_coord", "P_loss_phen", 
    "P_loss_hpo", "P_loss_source", "P_loss_coord", "P_ins_phen", 
    "P_ins_hpo", "P_ins_source", "P_ins_coord", "P_snvindel_nb", 
    "P_snvindel_phen", "B_gain_source", "B_gain_coord", "B_loss_source", 
    "B_loss_coord", "B_ins_source", "B_ins_coord", "B_inv_source", 
    "B_inv_coord", "TAD_coordinate", "ENCODE_experiment", "GC_content_left", 
    "GC_content_right", "Repeat_coord_left", "Repeat_type_left", 
    "Repeat_coord_right", "Repeat_type_right", "Gap_left", "Gap_right", 
    "SegDup_left", "SegDup_right", "ENCODE_blacklist_left", "ENCODE_blacklist_characteristics_left", 
    "ENCODE_blacklist_right", "ENCODE_blacklist_characteristics_right", 
    "ACMG", "HI", "TS", "DDD_HI_percent","ExAC_delZ","ExAC_dupZ","ExAC_cnvZ","ExAC_synZ",
    "ExAC_misZ","GenCC_disease","GenCC_moi","GenCC_classification","GenCC_pmid","OMIM_ID","OMIM_phenotype",
    "OMIM_inheritance","OMIM_morbid","OMIM_morbid_candidate","LOEUF_bin","GnomAD_pLI","ExAC_pLI",
    "AnnotSV_ranking_score","AnnotSV_ranking_criteria","ACMG_class"
  )
) %>% fwrite("ALL.master_gene_list.AnnotSV_ID.full.txt")


# identify all the translocations
all_result %>% filter(Annotation_mode == "split") %>% separate(V1, c("sample_name", "junk"), sep = ".merged") %>% separate(sample_name, c("junk", "sample_name"), sep = "./") %>% select(-"junk") %>% filter(SV_type == "TRA") %>% fwrite("TRA.AnnotSV_ID.csv")


# identify the 8 most commonly deleted genes from MLPA
## check subtypes
MLPA_gene_list.AnnotSV_ID <- all_result %>% filter(Annotation_mode == "split") %>% separate(V1, c("sample_name", "junk"), sep = ".merged") %>% separate(sample_name, c("junk", "sample_name"), sep = "./") %>% select(-"junk") %>% filter(Gene_name %in% c("CDKN2A", "IKZF1", "PAX5", "ETV6", "BTG1", "RB1", "EBF1", "CRLF2", "CSF2RA", "IL3RA")) %>% filter(SV_type == "DEL") %>% select(
  c("sample_name", "ID" ,"AnnotSV_ID", "SV_chrom", "SV_start", "SV_end", "SV_length", 
    "SV_type", "QUAL", "INFO", "Gene_name", "Tx", "Tx_start", "Tx_end", "Overlapped_tx_length", 
    "Overlapped_CDS_length", "Overlapped_CDS_percent", "Frameshift", 
    "Exon_count", "Location", "Location2", "Dist_nearest_SS", "Nearest_SS_type", 
    "Intersect_start", "Intersect_end", "RE_gene", "P_gain_phen", 
    "P_gain_hpo", "P_gain_source", "P_gain_coord", "P_loss_phen", 
    "P_loss_hpo", "P_loss_source", "P_loss_coord", "P_ins_phen", 
    "P_ins_hpo", "P_ins_source", "P_ins_coord", "P_snvindel_nb", 
    "P_snvindel_phen", "B_gain_source", "B_gain_coord", "B_loss_source", 
    "B_loss_coord", "B_ins_source", "B_ins_coord", "B_inv_source", 
    "B_inv_coord", "TAD_coordinate", "ENCODE_experiment", "GC_content_left", 
    "GC_content_right", "Repeat_coord_left", "Repeat_type_left", 
    "Repeat_coord_right", "Repeat_type_right", "Gap_left", "Gap_right", 
    "SegDup_left", "SegDup_right", "ENCODE_blacklist_left", "ENCODE_blacklist_characteristics_left", 
    "ENCODE_blacklist_right", "ENCODE_blacklist_characteristics_right", 
    "ACMG", "HI", "TS","DDD_HI_percent","ExAC_delZ","ExAC_dupZ","ExAC_cnvZ","ExAC_synZ",
    "ExAC_misZ","GenCC_disease","GenCC_moi","GenCC_classification","GenCC_pmid","OMIM_ID","OMIM_phenotype",
    "OMIM_inheritance","OMIM_morbid","OMIM_morbid_candidate","LOEUF_bin","GnomAD_pLI","ExAC_pLI",
    "AnnotSV_ranking_score","AnnotSV_ranking_criteria","ACMG_class"
  )
) %>% separate(INFO, c("junk", "SUPP_VEC"), sep = "SUPP_VEC=", remove = F) %>% select(-"junk") %>% separate(SUPP_VEC, c("SUPP_VEC", "junk"), sep = ";SVLEN=") %>% select(-"junk") %>% mutate(SUPP_Caller = case_when(
  SUPP_VEC == "011" ~ "lumpy, manta",
  SUPP_VEC == "101" ~ "delly, manta",
  SUPP_VEC == "110" ~ "delly, lumpy",
  SUPP_VEC == "111" ~ "delly, lumpy, manta"
)) %>% mutate(SV_length = as.numeric(SV_length)) %>% 
  mutate(SV_length_category = 
                case_when(
                  abs(SV_length) > 10000000  ~ "above 10mb",
                  abs(SV_length) > 1000000 & abs(SV_length) < 10000000 ~ "1mb~10mb",
                  abs(SV_length) < 50 & SV_type %in% c("DEL", "DUP") ~ "below 50b",
                  TRUE ~ "PASSED length",
                ) 
) 


####0214#######
MLPA_gene_list.AnnotSV_ID <- all_result %>% filter(Annotation_mode == "split") %>% separate(V1, c("sample_name", "junk"), sep = ".merged") %>% separate(sample_name, c("junk", "sample_name"), sep = "./") %>% select(-"junk") %>% filter(Gene_name %in% c("CDKN2A", "IKZF1", "PAX5", "ETV6", "BTG1", "RB1", "EBF1", "CRLF2", "CSF2RA", "IL3RA")) %>% filter(SV_type == "DEL") %>% select(
  c("sample_name", "ID" ,"AnnotSV_ID", "SV_chrom", "SV_start", "SV_end", "SV_length", 
    "SV_type", "QUAL", "INFO", "Gene_name", "Tx", "Tx_start", "Tx_end", "Overlapped_tx_length", 
    "Overlapped_CDS_length", "Overlapped_CDS_percent", "Frameshift", 
    "Exon_count", "Location", "Location2", "Dist_nearest_SS", "Nearest_SS_type", 
    "Intersect_start", "Intersect_end", "RE_gene", "P_gain_phen", 
    "P_gain_hpo", "P_gain_source", "P_gain_coord", "P_loss_phen", 
    "P_loss_hpo", "P_loss_source", "P_loss_coord", "P_ins_phen", 
    "P_ins_hpo", "P_ins_source", "P_ins_coord", "P_snvindel_nb", 
    "P_snvindel_phen", "B_gain_source", "B_gain_coord", "B_loss_source", 
    "B_loss_coord", "B_ins_source", "B_ins_coord", "B_inv_source", 
    "B_inv_coord", "TAD_coordinate", "ENCODE_experiment", "GC_content_left", 
    "GC_content_right", "Repeat_coord_left", "Repeat_type_left", 
    "Repeat_coord_right", "Repeat_type_right", "Gap_left", "Gap_right", 
    "SegDup_left", "SegDup_right", "ENCODE_blacklist_left", "ENCODE_blacklist_characteristics_left", 
    "ENCODE_blacklist_right", "ENCODE_blacklist_characteristics_right", 
    "ACMG", "HI", "TS","DDD_HI_percent","ExAC_delZ","ExAC_dupZ","ExAC_cnvZ","ExAC_synZ",
    "ExAC_misZ","GenCC_disease","GenCC_moi","GenCC_classification","GenCC_pmid","OMIM_ID","OMIM_phenotype",
    "OMIM_inheritance","OMIM_morbid","OMIM_morbid_candidate","LOEUF_bin","GnomAD_pLI","ExAC_pLI",
    "AnnotSV_ranking_score","AnnotSV_ranking_criteria","ACMG_class"
  )
) %>% separate(INFO, c("junk", "SUPP_VEC"), sep = "SUPP_VEC=", remove = F) %>% select(-"junk") %>% separate(SUPP_VEC, c("SUPP_VEC", "junk"), sep = ";SVLEN=") %>% select(-"junk") %>% mutate(SUPP_Caller = case_when(
  SUPP_VEC == "011" ~ "lumpy, manta",
  SUPP_VEC == "101" ~ "delly, manta",
  SUPP_VEC == "110" ~ "delly, lumpy",
  SUPP_VEC == "111" ~ "delly, lumpy, manta"
)) %>% mutate(SV_length = as.numeric(SV_length)) %>% 
  mutate(SV_length_category = 
           case_when(
             abs(SV_length) > 10000000  ~ "above 10mb",
             abs(SV_length) > 1000000 & abs(SV_length) < 10000000 ~ "1mb~10mb",
             abs(SV_length) < 50 & SV_type %in% c("DEL", "DUP") ~ "below 50b",
             TRUE ~ "PASSED length",
           ) 
  ) 

all_content.delly.master <- fread("/scratch2/tanxinli/wgs.smoking.out/survivor.out/all_content.delly.master.txt") %>% separate(Indiv, c("junk", "sample_name"), "TM_") %>% select(-"junk") %>% mutate(sample_name = as.numeric(sample_name)) %>% separate(CHROM, c("junk", "SV_chrom"), "chr") %>% select(-"junk") %>% mutate(VAF_delly = case_when(
  IMPRECISE == "TRUE" ~ gt_DV/(gt_DR + gt_DV),
  IMPRECISE == "FALSE" ~ gt_RV/(gt_RR + gt_RV)
)) 

all_content.manta.master <- fread("/scratch2/tanxinli/wgs.smoking.out/survivor.out/all_content.manta.master.txt") %>% separate(Indiv, c("junk", "sample_name"), "TM_") %>% select(-"junk") %>% mutate(sample_name = as.numeric(sample_name)) %>% separate(CHROM, c("junk", "SV_chrom"), "chr") %>% select(-"junk") %>% mutate(gt_SR = ifelse(gt_SR=="", "0,0", gt_SR)) %>% separate(gt_PR, c("gt_PR_ref", "gt_PR_alt"), sep = ",") %>% separate(gt_SR, c("gt_SR_ref", "gt_SR_alt"), sep = ",") %>% mutate_at(c("gt_PR_ref", "gt_PR_alt", "gt_SR_ref", "gt_SR_alt"), as.numeric) %>% mutate(VAF_manta = (gt_PR_alt + gt_SR_alt)/(gt_PR_ref + gt_SR_ref + gt_PR_alt + gt_SR_alt))

all_content.manta.master <- all_content.manta.master %>% select(sample_name, ID, VAF_manta)
all_content.delly.master <- all_content.delly.master %>% select(sample_name, ID, VAF_delly)
all_content.lumpy.master <- all_content.lumpy.master %>% select(sample_name, ID, VAF_lumpy)

#don't have MLPA_gene_list
MLPA_gene_list.AnnotSV_ID$sample_name <- as.numeric(MLPA_gene_list.AnnotSV_ID$sample_name)

MLPA_gene_list.AnnotSV_ID.VAF_manta <- MLPA_gene_list.AnnotSV_ID %>% inner_join(all_content.manta.master, c("sample_name", "ID"))  # 11 samples
MLPA_gene_list.AnnotSV_ID.VAF_lumpy <- MLPA_gene_list.AnnotSV_ID %>% inner_join(all_content.lumpy.master, c("sample_name", "ID"))  # 63 samples
# MLPA_gene_list.AnnotSV_ID.VAF_manta_lumpy <- rbind(
#   MLPA_gene_list.AnnotSV_ID.VAF_manta,
#   MLPA_gene_list.AnnotSV_ID.VAF_lumpy,
#   fill = T
# )

MLPA_gene_list.AnnotSV_ID.VAF_manta_lumpy <- rbind(
  MLPA_gene_list.AnnotSV_ID.VAF_manta,
  MLPA_gene_list.AnnotSV_ID.VAF_lumpy, fill=T
)
MLPA_gene_list.AnnotSV_ID.VAF_manta_lumpy <- MLPA_gene_list.AnnotSV_ID.VAF_manta_lumpy %>% arrange(sample_name, SV_chrom, SV_start, SV_end)

# print out all the SVs on commonly deleted genes 
fwrite(MLPA_gene_list.AnnotSV_ID.VAF_manta_lumpy, "MLPA_gene_list.AnnotSV_ID.csv")


# *********
ALL_smoking_wgs.sampleid_MLPAdata <- fread("ALL_smoking_wgs.sampleid_MLPAdata.csv") %>% mutate(subjectId = as.character(subjectId)) %>% mutate(subjectId = ifelse( subjectId == "36077", "036077", subjectId))

MLPA_gene_list.AnnotSV_ID.summary <- MLPA_gene_list.AnnotSV_ID %>% select(sample_name, Gene_name) %>% mutate(Gene_name = fct_relevel(Gene_name, c("CDKN2A", "IKZF1", "PAX5", "ETV6", "BTG1", "RB1", "EBF1", "CRLF2", "CSF2RA", "IL3RA"))) %>% arrange(Gene_name) %>% group_by(sample_name) %>% distinct(Gene_name) %>% summarise(Gene_name = paste(Gene_name, collapse =",")) %>% rename("subjectId" = "sample_name") %>% rename("SV_analysis_del" = "Gene_name")

ALL_smoking_wgs.sampleid_MLPAdata %>% left_join(MLPA_gene_list.AnnotSV_ID.summary, by = "subjectId") %>% fwrite("ALL_smoking_wgs.sampleid.MLPAdata.SV_summary.csv")



# report all the SVs
ALL.master.AnnotSV_ID <- all_result %>% filter(Annotation_mode == "full") %>% separate(V1, c("sample_name", "junk"), sep = ".merged") %>% separate(sample_name, c("junk", "sample_name"), sep = "./") %>% select(-"junk")

ALL.master.AnnotSV_ID %>% select(
  c("sample_name", "AnnotSV_ID","ID","SV_chrom", "SV_start", "SV_end", "SV_length", 
    "SV_type", "QUAL", "INFO", "Gene_name", "Tx", "Tx_start", "Tx_end", "Overlapped_tx_length", 
    "Overlapped_CDS_length", "Overlapped_CDS_percent", "Frameshift", 
    "Exon_count", "Location", "Location2", "Dist_nearest_SS", "Nearest_SS_type", 
    "Intersect_start", "Intersect_end", "RE_gene", "P_gain_phen", 
    "P_gain_hpo", "P_gain_source", "P_gain_coord", "P_loss_phen", 
    "P_loss_hpo", "P_loss_source", "P_loss_coord", "P_ins_phen", 
    "P_ins_hpo", "P_ins_source", "P_ins_coord", "P_snvindel_nb", 
    "P_snvindel_phen", "B_gain_source", "B_gain_coord", "B_loss_source", 
    "B_loss_coord", "B_ins_source", "B_ins_coord", "B_inv_source", 
    "B_inv_coord", "TAD_coordinate", "ENCODE_experiment", "GC_content_left", 
    "GC_content_right", "Repeat_coord_left", "Repeat_type_left", 
    "Repeat_coord_right", "Repeat_type_right", "Gap_left", "Gap_right", 
    "SegDup_left", "SegDup_right", "ENCODE_blacklist_left", "ENCODE_blacklist_characteristics_left", 
    "ENCODE_blacklist_right", "ENCODE_blacklist_characteristics_right", 
    "ACMG", "HI", "TS", "DDD_HI_percent", "DDD_status", "DDD_mode", 
    "DDD_consequence", "DDD_disease", "DDD_pmid", "ExAC_delZ", "ExAC_dupZ", 
    "ExAC_cnvZ", "ExAC_synZ", "ExAC_misZ", "OMIM_ID", "OMIM_phenotype", 
    "OMIM_inheritance", "OMIM_morbid", "OMIM_morbid_candidate", "LOEUF_bin", 
    "GnomAD_pLI", "ExAC_pLI", "Exomiser_gene_pheno_score", "Human_pheno_evidence", 
    "Mouse_pheno_evidence", "Fish_pheno_evidence", "AnnotSV_ranking_score", 
    "AnnotSV_ranking_criteria", "ACMG_class"
  )
) %>% fwrite("ALL.master.AnnotSV_ID.full.txt")


setwd("./2021_WGS_ALL_smoking/sv")
fread("ALL.master.AnnotSV_ID.full.anno.txt") %>% select(
  c("sample_name", "AnnotSV_ID", "ID", "SV_chrom", "SV_start", 
    "SV_end", "SV_length", "SV_type", "QUAL", 
    "Gene_name", "RE_gene",  "SUPP_Caller", "SV_length_category", 
    "VAF_manta", "VAF_lumpy", "VAF_delly", "VAF")
) %>% fwrite("ALL.master.AnnotSV_ID.full.anno.cleaned.csv")

fread("ALL.master_gene_list.AnnotSV_ID.full.anno.txt") %>% select(
  c("sample_name", "AnnotSV_ID", "ID", "SV_chrom", "SV_start", 
    "SV_end", "SV_length", "SV_type", "QUAL", 
    "Gene_name", "RE_gene",  "SUPP_Caller", "SV_length_category", 
    "VAF_manta", "VAF_lumpy", "VAF")
) %>% fwrite("ALL.master_gene_list.AnnotSV_ID.full.anno.cleaned.csv")


