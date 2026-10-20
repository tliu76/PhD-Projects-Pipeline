library(data.table)
library(tidyverse)

# sex data 
smoking.sex <- fread("sample_sex_cnvkit.txt", header = T) 

# subtype data
smoking.subtype <- fread("20220316.summary_smoking.sigprofiler.subtype.csv") %>% separate(sampleid, c("junk", "sampleid"), sep = "T") %>% select(-"junk") %>% mutate(sampleid = as.numeric(sampleid))

# METH score
smoking.meth <- fread("mclust.results.summary.txt") %>% rename("sampleid" = "sample_id") %>% rename("sclust.cluster_number" = "cluster_number", "sclust.MATH_score" = "MATH_score") %>% separate(sampleid, c("junk", "sampleid"), sep = "TM_") %>% select(-"junk") %>% mutate(sampleid = as.numeric(sampleid))

# signature data 
smoking.sigprofiler.sbs96.4denov <- fread("./BCFtools.filter.tumoronly.autosomal.vcf/results/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities_refit.txt") %>% separate(Samples, c("junk", "sampleid"), sep = "sample") %>% select(-"junk") %>% mutate(sampleid = as.numeric(sampleid)) %>% rename(setNames(names(.[,2:dim(.)[2]]), paste0('sigprofiler.sbs96.4denovo.', names(.[,2:dim(.)[2]]))))

smoking.sigprofiler.sbs96.5denov <- fread("./BCFtools.filter.tumoronly.autosomal.vcf/results/SBS96/output_5sig/Decompose_Solution/Activities/Decompose_Solution_Activities.txt") %>% separate(Samples, c("junk", "sampleid"), sep = "sample") %>% select(-"junk") %>% mutate(sampleid = as.numeric(sampleid)) %>% rename(setNames(names(.[,2:dim(.)[2]]), paste0('sigprofiler.sbs96.5denovo.', names(.[,2:dim(.)[2]]))))

smoking.sigprofiler.dbs78 <- fread("./BCFtools.filter.tumoronly.autosomal.vcf/results/DBS78/Suggested_Solution/COSMIC_DBS78_Decomposed_Solution/Activities/COSMIC_DBS78_Activities_refit.txt") %>% separate(Samples, c("junk", "sampleid"), sep = "sample") %>% select(-"junk") %>% mutate(sampleid = as.numeric(sampleid)) %>% rename(setNames(names(.[,2:dim(.)[2]]), paste0('sigprofiler.dbs78.', names(.[,2:dim(.)[2]]))))

smoking.sigprofiler.id83 <- fread("./BCFtools.filter.tumoronly.autosomal.vcf/results/ID83/Suggested_Solution/COSMIC_ID83_Decomposed_Solution/Activities/COSMIC_ID83_Activities_refit.txt") %>% separate(Samples, c("junk", "sampleid"), sep = "sample") %>% select(-"junk") %>% mutate(sampleid = as.numeric(sampleid)) %>% rename(setNames(names(.[,2:dim(.)[2]]), paste0('sigprofiler.id83.', names(.[,2:dim(.)[2]]))))


smoking.pmsignature <- fread("./pmsignature.out/MembershipValuek5.csv") %>% separate(V1, c("junk", "sampleid"), sep = "sample") %>% select(-"junk") %>% mutate(sampleid = as.numeric(sampleid)) %>% rename("pmsignature.SBS45" = "signature_1", "pmsignature.SBS40" = "signature_2", "pmsignature.SBS57" = "signature_3", "pmsignature.SBS1" = "signature_4", "pmsignature.SBS40.again" = "signature_5")

smoking.deconstructSigsCosmicV2 <- fread("./deconstructSigs.out/deconvolute_df_output_summary.txt") %>% separate(sampleid, c("junk", "sampleid"), sep = "sample") %>% select(-"junk") %>% mutate(sampleid = as.numeric(sampleid)) %>% rename(setNames(names(.[,2:dim(.)[2]]), paste0('deconstructSigsCosmicV2.', names(.[,2:dim(.)[2]]))))

# remove nonexisting signatures from smoking.deconstructSigsCosmicV2
drop_columns <- sapply(smoking.deconstructSigsCosmicV2, function(x) all(x == 0)) %>% .[. == "TRUE"] %>% names(.)
smoking.deconstructSigsCosmicV2 <- smoking.deconstructSigsCosmicV2 %>% select(-drop_columns)

# kataegis fom signal
smoking.kataegis.signal <- fread("signal_kataegis_regions.txt") %>% unite(regions, c("start_pos", "end_pos"), sep = "-") %>% unite(regions, c("chr", "regions"), sep = ":") %>% group_by(sample_id) %>% arrange(regions) %>% summarize(regions = paste(regions, collapse = ",")) %>% rename("kataegis.signal.regions" = "regions") %>% separate(sample_id, c("junk", "sampleid"), sep = "TM_") %>% select(-"junk") %>% mutate(sampleid = as.numeric(sampleid)) %>% mutate(kataegis.signal = "Yes")

# kataegis from maftools
## read file path
all_paths <- list.files(path = "vcf2maf.maftools.out",pattern = "*_Kataegis.tsv",full.names = TRUE)
## read file content
all_content <- all_paths %>% lapply(fread, header = TRUE)
# unlist all lists and change column name
all_result <- rbindlist(all_content, fill = T)
smoking.kataegis.maftools <- all_result %>% select(Chromosome, Start_Position, End_Position, Tumor_Sample_Barcode) %>% unite(regions, c("Start_Position", "End_Position"), sep = "-") %>% unite(regions, c("Chromosome", "regions"), sep = ":") %>% group_by(Tumor_Sample_Barcode) %>% arrange(regions) %>% summarize(regions = paste(regions, collapse = ",")) %>% rename("sampleid" = "Tumor_Sample_Barcode", "kataegis.maftools.regions" = "regions")  %>% separate(sampleid, c("junk", "sampleid"), sep = "TM_") %>% select(-"junk") %>% mutate(sampleid = as.numeric(sampleid)) %>% mutate(kataegis.maftools = "Yes")

# APOBEC enrichment score > 2
smoking.apobec_enrich <- fread("./smoking.APOBEC_scores.txt") %>% select(Tumor_Sample_Barcode, APOBEC_Enrichment, fraction_APOBEC_mutations, APOBEC_Enriched) %>% dplyr::rename("APOBEC_Enrichment_score" = "APOBEC_Enrichment") %>% dplyr::rename("sampleid" = "Tumor_Sample_Barcode") %>% separate(sampleid, c("junk", "sampleid"), sep = "TM_") %>% select(-"junk") %>% mutate(sampleid = as.numeric(sampleid)) 

# Oncogenic Signaling Pathways
smoking_maf <- maftools::read.maf(maf = "./vcf2maf.maftools.out/smoking.38.merged.maf", isTCGA = F)
smoking.nonsyn <- smoking_maf@variant.classification.summary %>% dplyr::rename(setNames(names(.[,2:dim(.)[2]]), paste0('nonsyn.', names(.[,2:dim(.)[2]])))) %>% dplyr::rename("sampleid" = "Tumor_Sample_Barcode") %>% separate(sampleid, c("junk", "sampleid"), sep = "TM_") %>% select(-"junk") %>% mutate(sampleid = as.numeric(sampleid)) 

df_OncogenicPathways_RTK_RAS <- data.frame(sampleid=c("114694", "314630", "712286", "737126", "871719", "855078", "125806", "474264", "211349", "439723", "680110", "855606","217737","947007","882476","942742","705771"), OncogenicPathways_RTK_RAS = rep(c("Yes"), 17)) %>% mutate(sampleid = as.numeric(sampleid)) 
df_OncogenicPathways_NOTCH <- data.frame(sampleid=c("114694","749755","321666","129204","375924"), OncogenicPathways_NOTCH = rep(c("Yes"), 5)) %>% mutate(sampleid = as.numeric(sampleid)) 
df_OncogenicPathways_WNT <- data.frame(sampleid=c("211349","217737","374219"), OncogenicPathways_WNT = rep(c("Yes"), 3)) %>% mutate(sampleid = as.numeric(sampleid)) 
df_OncogenicPathways_Hippo <- data.frame(sampleid=c("855606","253161"), OncogenicPathways_Hippo = rep(c("Yes"), 2)) %>% mutate(sampleid = as.numeric(sampleid)) 
df_OncogenicPathways_MYC <- data.frame(sampleid=c("712286"), OncogenicPathways_MYC = rep(c("Yes"), 1)) %>% mutate(sampleid = as.numeric(sampleid)) 
df_OncogenicPathways_Cell_Cycle <- data.frame(sampleid=c("808946", "464694"), OncogenicPathways_Cell_Cycle = rep(c("Yes"), 2)) %>% mutate(sampleid = as.numeric(sampleid)) 

smoking.OncogenicPathways <- smoking.nonsyn %>% left_join(df_OncogenicPathways_RTK_RAS, by = "sampleid") %>% left_join(df_OncogenicPathways_NOTCH, by = "sampleid") %>% left_join(df_OncogenicPathways_WNT, by = "sampleid") %>% left_join(df_OncogenicPathways_Hippo, by = "sampleid") %>% left_join(df_OncogenicPathways_MYC, by = "sampleid") %>% left_join(df_OncogenicPathways_Cell_Cycle, by = "sampleid") 

smoking.OncogenicPathways[is.na(smoking.OncogenicPathways)] <- "No"

# number of SVs

## read file path
all_paths <- list.files(path = "./sv/AnnotSV_out_two_caller", pattern = "*.annotated.tsv",full.names = TRUE)
## read file content
all_content <- all_paths %>% lapply(fread, header = TRUE, colClasses="character")
## read file name
all_filenames <- all_paths %>% as.list()
## combine file content list and file name list
all_lists <- mapply(c, all_content, all_filenames, SIMPLIFY = FALSE)
## unlist all lists and change column name
all_result <- rbindlist(all_lists, fill = T)
## check SV counts for each sample use MODE =FULL 
# MODE = FULL have the full records even those without gene anno
smoking.SV <- all_result %>% filter(Annotation_mode == "full") %>%separate(V1, c("sample_name", "junk"), sep = ".merged") %>% separate(sample_name, c("junk", "sampleid"), sep = "./sv/AnnotSV_out_two_caller/") %>% count(sampleid, SV_type) %>% mutate(sampleid = as.numeric(sampleid)) %>% pivot_wider(names_from = SV_type, values_from = n) %>% group_by(sampleid) %>% mutate(total=sum(DEL,DUP,INV,TRA, na.rm = T)) %>% ungroup() %>% dplyr::rename(setNames(names(.[,2:dim(.)[2]]), paste0('svtype.', names(.[,2:dim(.)[2]]))))
## replace all the na values with 0
smoking.SV[is.na(smoking.SV)] <- 0

# number of SVs on candidate genes 
ALL_master_genes <- fread("./sv/AnnotSV_out_two_caller/ALL.master_gene_list.AnnotSV_ID.split.csv")

smoking.SV.candidate_gene <- ALL_master_genes %>% mutate(sampleid = as.numeric(sample_name)) %>% count(sampleid, SV_type) %>% pivot_wider(names_from = SV_type, values_from = n) %>% group_by(sampleid) %>% mutate(total=sum(DEL,DUP,INV,TRA, na.rm = T)) %>% ungroup() %>% dplyr::rename(setNames(names(.[,2:dim(.)[2]]), paste0('svtype.candidate_gene.', names(.[,2:dim(.)[2]]))))
## replace all the na values with 0
smoking.SV.candidate_gene[is.na(smoking.SV.candidate_gene)] <- 0

# tellength
smoking.38cohort.T.N.ratio <- fread("./telomere.length/smoking.38cohort.T.N.ratio.txt") %>% mutate(sampleid = Sample) %>% select(-"Sample") %>% mutate(sampleid = as.numeric(sampleid)) 

# convert mutations counts to proportions 
smoking.sigprofiler.sbs96.4denov$row_sum <- rowSums(smoking.sigprofiler.sbs96.4denov[,2:10])
smoking.sigprofiler.sbs96.5denov$row_sum <- rowSums(smoking.sigprofiler.sbs96.5denov[,2:10])
smoking.sigprofiler.dbs78$row_sum <- rowSums(smoking.sigprofiler.dbs78[,2])
smoking.sigprofiler.id83$row_sum <- rowSums(smoking.sigprofiler.id83[,2:3])

smoking.sigprofiler.sbs96.4denov <- smoking.sigprofiler.sbs96.4denov %>% mutate(across(c(2:10), ~ .x / row_sum)) %>% dplyr::rename("SBS_total" = "row_sum") %>% relocate(SBS_total)
smoking.sigprofiler.sbs96.5denov <- smoking.sigprofiler.sbs96.5denov %>% mutate(across(c(2:10), ~ .x / row_sum)) %>% select(-"row_sum")
smoking.sigprofiler.dbs78 <- smoking.sigprofiler.dbs78 %>% mutate(across(c(2), ~ .x / row_sum)) %>% dplyr::rename("DBS_total" = "row_sum") %>% relocate(DBS_total)
smoking.sigprofiler.id83 <- smoking.sigprofiler.id83 %>% mutate(across(c(2:3), ~ .x / row_sum)) %>% dplyr::rename("ID_total" = "row_sum") %>% relocate(ID_total)

# signal.out results 
signal.organ <- fread("./signal.out.contribution.txt") %>% mutate(proportion = Contribution/total_base) %>% select(sample_id, Signature, proportion) %>% distinct() %>% pivot_wider(names_from = Signature, values_from = proportion) %>% rename(setNames(names(.[,2:dim(.)[2]]), paste0('signal.organ.', names(.[,2:dim(.)[2]])))) %>% separate(sample_id, c("junk", "sampleid"), sep = "TM_") %>% select(-"junk") %>% mutate(sampleid = as.numeric(sampleid))

## replace all the na values with 0
signal.organ[is.na(signal.organ)] <- 0

signal.reference <- fread("./signal.out.reference_contribution.txt") %>% mutate(proportion = Contribution/total_base) %>% select(sample_id, Signature, proportion) %>% distinct() %>% pivot_wider(names_from = Signature, values_from = proportion) %>% rename(setNames(names(.[,2:dim(.)[2]]), paste0('signal.reference.', names(.[,2:dim(.)[2]])))) %>% separate(sample_id, c("junk", "sampleid"), sep = "TM_") %>% select(-"junk") %>% mutate(sampleid = as.numeric(sampleid))

signal.reference[is.na(signal.reference)] <- 0

# fitms results 
fitms.reference <- fread(".fitms.out.reference_contribution.txt") %>% mutate(proportion = Contribution/total_base) %>% select(sample_id, Signature, proportion) %>% distinct() %>% pivot_wider(names_from = Signature, values_from = proportion) %>% rename(setNames(names(.[,2:dim(.)[2]]), paste0('fitms.reference.', names(.[,2:dim(.)[2]])))) %>% rename("sampleid" = "sample_id") %>% mutate(sampleid = as.numeric(sampleid))

fitms.reference[is.na(fitms.reference)] <- 0

# exposure group 
Tobacco.exposure.group <- fread("smoking.Tobacco.exposure.group.csv")

# MLPA deletion data 
MLPA_deletion <- fread("smoking.MLPA_deletion.csv") %>% rename("sampleid" = "subjectId") %>% rename("MLPA_deletion_total" = "TOTAL")

# merge datasets 
smoking.sex %>% left_join(smoking.subtype, by = c("sampleid"))  %>% left_join(smoking.sigprofiler.sbs96.4denov, by = c("sampleid")) %>% left_join(smoking.sigprofiler.sbs96.5denov, by = c("sampleid")) %>% left_join(smoking.sigprofiler.dbs78, by = c("sampleid")) %>% left_join(smoking.sigprofiler.id83, by = c("sampleid")) %>% left_join(smoking.pmsignature, by = c("sampleid")) %>% left_join(smoking.deconstructSigsCosmicV2, by = c("sampleid"))%>% left_join(smoking.kataegis.signal, by = c("sampleid"))  %>%  left_join(smoking.kataegis.maftools, by = c("sampleid")) %>% left_join(smoking.apobec_enrich, by = c("sampleid")) %>% left_join(smoking.OncogenicPathways, by = c("sampleid")) %>% left_join(smoking.SV, by = c("sampleid")) %>% left_join(smoking.SV.candidate_gene, by = c("sampleid")) %>% left_join(smoking.meth, by = c("sampleid")) %>% left_join(smoking.38cohort.T.N.ratio , by = c("sampleid")) %>% left_join(signal.organ, by = c("sampleid")) %>% left_join(signal.reference, by = c("sampleid")) %>% left_join(fitms.reference, by = c("sampleid")) %>% left_join(Tobacco.exposure.group, by = c("sampleid")) %>% left_join(MLPA_deletion, by = c("sampleid")) %>% fwrite("smoking.38cohort.master.txt")


smoking.sex %>% left_join(smoking.subtype, by = c("sampleid"))  %>% left_join(smoking.sigprofiler.sbs96.4denov, by = c("sampleid")) %>% left_join(smoking.sigprofiler.sbs96.5denov, by = c("sampleid")) %>% left_join(smoking.sigprofiler.dbs78, by = c("sampleid")) %>% left_join(smoking.sigprofiler.id83, by = c("sampleid"))%>% left_join(smoking.pmsignature, by = c("sampleid")) %>% left_join(smoking.deconstructSigsCosmicV2, by = c("sampleid"))%>% left_join(smoking.kataegis.signal, by = c("sampleid"))  %>%  left_join(smoking.kataegis.maftools, by = c("sampleid")) %>% left_join(smoking.apobec_enrich, by = c("sampleid")) %>% left_join(smoking.OncogenicPathways, by = c("sampleid")) %>% left_join(smoking.SV, by = c("sampleid")) %>% left_join(smoking.SV.candidate_gene, by = c("sampleid")) %>% left_join(smoking.meth, by = c("sampleid")) %>% left_join(smoking.38cohort.T.N.ratio , by = c("sampleid")) %>% left_join(signal.organ, by = c("sampleid")) %>% left_join(signal.reference, by = c("sampleid")) %>% left_join(fitms.reference, by = c("sampleid")) %>% left_join(Tobacco.exposure.group, by = c("sampleid")) %>% left_join(MLPA_deletion, by = c("sampleid")) %>% fwrite("smoking.38cohort.master.csv")

fread("smoking.38cohort.master.txt") %>% select(SBS_total) %>% summarize(min(SBS_total), max(SBS_total), median(SBS_total))
