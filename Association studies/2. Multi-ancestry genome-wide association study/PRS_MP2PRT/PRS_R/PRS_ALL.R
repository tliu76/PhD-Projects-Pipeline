### Compare ALL PRS with Nathan's PRS:
.libPaths("/home1/tanxinli/R_lib")
library(data.table)
prs <- fread("/scratch1/tanxinli/mp2prt_gwas/PRS/prs_ALL_results/ALL_prs.sscore")
prs$PRS_z <- scale(prs$SCORE1_AVG)
fwrite(prs,"/scratch1/tanxinli/mp2prt_gwas/PRS/prs_ALL_results/ALL_prs_standardized.csv")



# preview
head(prs)
# ALLELE_CT NAMED_ALLELE_DOSAGE_SUM SCORE1_AVG                                  IID      PRS_z
# <int>                   <num>      <num>                               <char>      <num>
#   1:        36                  12.991  0.0125909 003d3cb6-b01b-4d2c-9aa7-89e23c01cc8f -1.0076292
# 2:        36                  15.029  0.0477730 0042381f-abda-412d-b7fb-1ad3377338b7  0.8896892
# 3:        36                  10.992  0.0240764 004b16da-cc31-4c90-ac1d-8d860b20f10b -0.3882332
# 4:        36                   9.052  0.0416678 00e2c2f7-f98e-491a-8731-74302deedf77  0.5604448
# 5:        36                  13.003  0.0332797 0178bd7f-9e21-4d5e-a5a5-ae72ff5e38ae  0.1080871
# 6:        36                   9.981  0.0201475 01916793-36cd-4412-a8b6-d2bf6a06a39b -0.6001130
library(data.table)
library(ggplot2)

# Read the PRS file

prs$IID <- prs$`#IID`
prs$`#IID` <- NULL
# Density plot of PRS
ggplot(prs, aes(x = PRS_z)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  theme_bw() +
  labs(
    title = "Acute lymphoblastic leukemia PRS",
    x = "PRS (PRS_z)",
    y = "Density"
  )


###Read Nathan's PRS###

library(readxl)
Transethnic_PRS_MP2PRT <- read_excel("/scratch1/tanxinli/mp2prt_gwas/PRS/Transethnic_PRS_MP2PRT.xlsx")
View(Transethnic_PRS_MP2PRT)

prs_merge <- prs %>% left_join(Transethnic_PRS_MP2PRT , by ="IID")
fwrite(prs_merge,"/scratch1/tanxinli/mp2prt_gwas/PRS/prs_ALL_results/ALL_prs_merge_standardized.csv")

library(dplyr)
library(tidyr)
library(ggplot2)

df_long <- prs_merge %>%
  select(PRS_z, standardized_PRS_MP2PRT_only) %>%
  pivot_longer(
    cols = everything(),
    names_to = "PRS",
    values_to = "score"
  )


ggplot(df_long, aes(x = score, color = PRS, fill = PRS)) +
  geom_density(alpha = 0.3) +
  geom_vline(
    data = df_long %>% group_by(PRS) %>% summarize(mu = mean(score, na.rm = TRUE)),
    aes(xintercept = mu, color = PRS),
    linetype = "dashed"
  ) +
  theme_bw() +
  labs(title = "PRS Density with Mean Lines")


ggplot(df_long, aes(x = score, fill = PRS)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~PRS, scales = "free") +
  theme_bw()



ggplot(df_long, aes(x = score, color = PRS)) +
  geom_density(size = 1.2) +
  theme_bw() +
  labs(
    title = "Overlapped Density Plot of Two PRS",
    x = "PRS Value",
    y = "Density"
  )


#### Correlation between two ALL PRS###
cor_test <- cor.test(
  prs_merge$PRS_z,
  prs_merge$standardized_PRS_MP2PRT_only,
  method = "pearson"
)

r_val <- cor_test$estimate
p_val <- cor_test$p.value


library(ggplot2)

ggplot(prs_merge , aes(x = PRS_z, y = standardized_PRS_MP2PRT_only)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red") +
  theme_bw() +
  labs(
    title = "Correlation Between Two PRS",
    x = "PRS: SCORE1_AVG",
    y = "PRS: standardized_PRS_MP2PRT_only"
  )+
  annotate(
    "text",
    x = 0, y = Inf,
    label = paste0(
      "r = ", round(r_val, 3),
      "\nP = ", format.pval(p_val, digits = 3)
    ),
    hjust = 1.1, vjust = 1.1,
    size = 4
  )

library(ggpubr)

ggplot(prs_merge, aes(x = PRS_z, y = standardized_PRS_MP2PRT_only)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  theme_bw() +
  labs(
    title = "Correlation Between Two PRS",
    x = "PRS: PRS_z",
    y = "PRS: standardized_PRS_MP2PRT_only"
  )

####################################################################################################################
#### Blood cell trait PRS ###
library(data.table)
library(dplyr)
library(stringr)

# PRS directory
prs_dir <- "/scratch1/tanxinli/mp2prt_gwas/PRS/prs_results/"

# List all .sscore files
files <- list.files(
  path = prs_dir,
  pattern = "\\.sscore$",
  full.names = TRUE
)

prs_list <- lapply(files, function(f) {
  tmp <- fread(f)
  
  # Extract PRS name (file name without .sscore)
  prs_name <- f %>%
    basename() %>%
    str_remove("\\.sscore$")
  
  # Ensure columns exist
  if (!"#IID" %in% names(tmp)) stop("Missing #IID in: ", f)
  if (!"SCORE1_AVG" %in% names(tmp)) stop("Missing SCORE1_AVG in: ", f)
  
  # Extract ID + PRS value
  tmp <- tmp %>%
    select(`#IID`, SCORE1_AVG) %>%
    rename(!!prs_name := SCORE1_AVG)
  
  return(tmp)
})




# Merge all PRS tables by #IID
prs_all <- Reduce(function(x, y) merge(x, y, by = "#IID", all = TRUE), prs_list)

head(prs_all)
### Created 9 standardized PRS ###
prs_cols <- c(
  "Baso.z_prs",
  "Eosin.z_prs",
  "LMR_prs",
  "Lymph.z_prs",
  "Mono.z_prs",
  "Neutro.z_prs",
  "NLR_prs",
  "Platelet.z_prs",
  "PLR_prs"
)

prs_all$Baso.z_prs_z <- scale(prs_all$Baso.z_prs)
prs_all$Eosin.z_prs_z <- scale(prs_all$Eosin.z_prs)
prs_all$LMR_prs_z <- scale(prs_all$LMR_prs)
prs_all$Lymph.z_prs_z <- scale(prs_all$Lymph.z_prs)
prs_all$Mono.z_prs_z <- scale(prs_all$Mono.z_prs)
prs_all$Neutro.z_prs_z <- scale(prs_all$Neutro.z_prs)
prs_all$NLR_prs_z <- scale(prs_all$NLR_prs)
prs_all$Platelet.z_prs_z <- scale(prs_all$Platelet.z_prs)
prs_all$PLR_prs_z <- scale(prs_all$PLR_prs)



library(ggplot2)

ggplot(prs_all, aes(x = Lymph.z_prs_z)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  theme_classic() +
  labs(title = "Distribution of Lymphocyte PRS")


### Read Nathan's PRS first table ###
library(readxl)
lymphocyte_and_birthweight_MP2PRT_PRS <- read_excel("/scratch1/tanxinli/mp2prt_gwas/PRS/lymphocyte and birthweight MP2PRT PRS.xlsx") 
lymphocyte_prs <- lymphocyte_and_birthweight_MP2PRT_PRS[,c(1,3:6,8:11)]

names(lymphocyte_prs)
# 1] "IID"                                      "lymph_standardized_PRS"                   "lmr_standardized_PRS"                    
# [4] "nlr_standardized_PRS"                     "plr_standardized_PRS"                     "lymph_standardized_with_SMILES_ctrls_PRS"
# [7] "lmr_standardized_with_SMILES_ctrls_PRS"   "nlr_standardized_with_SMILES_ctrls_PRS"   "plr_standardized_with_SMILES_ctrls_PRS"

names(prs_all)
prs_all_standardized <- prs_all[,c(1,11:19)]
# [1] "#IID"             "Baso.z_prs"       "Eosin.z_prs"      "LMR_prs"          "Lymph.z_prs"      "Mono.z_prs"       "Neutro.z_prs"    
# [8] "NLR_prs"          "Platelet.z_prs"   "PLR_prs"          "Eosin.z_prs_z"    "Baso.z_prs_z"     "PLR_prs_z"        "Platelet.z_prs_z"
# [15] "NLR_prs_z"        "Neutro.z_prs_z"   "Mono.z_prs_z"     "Lymph.z_prs_z"    "LMR_prs_z"   

prs_all_standardized$IID <- prs_all_standardized$`#IID`
prs_all_standardized$`#IID` <- NULL

prs_blood_merge <- lymphocyte_prs %>% left_join(prs_all_standardized , by ="IID")
names(prs_blood_merge)
fwrite(prs_blood_merge,"/scratch1/tanxinli/mp2prt_gwas/PRS/prs_ALL_results/blood_prs_merge_standardized.csv")


### 画图 ###
prs_map <- tibble::tribble(
  ~group,  ~prs1,              ~prs2,
  "Lymph", "Lymph.z_prs_z",    "lymph_standardized_PRS",
  "LMR",   "LMR_prs_z",        "lmr_standardized_PRS",
  "NLR",   "NLR_prs_z",        "nlr_standardized_PRS",
  "PLR",   "PLR_prs_z",        "plr_standardized_PRS"
)


plot_density <- function(prs_blood_merge, prs1_col, prs2_col, group_name) {
  ggplot(prs_blood_merge) +
    geom_density(aes(x = .data[[prs1_col]], color = "My PRS"), size = 1) +
    geom_density(aes(x = .data[[prs2_col]], color = "Standardized PRS"), size = 1) +
    theme_bw() +
    labs(
      title = paste0(group_name, " PRS Density"),
      x = "PRS value",
      y = "Density",
      color = ""
    )
}
density_plots <- purrr::pmap(
  list(prs_map$prs1, prs_map$prs2, prs_map$group),
  ~ plot_density(prs_blood_merge, ..1, ..2, ..3)
)
library(patchwork)
density_panel <- wrap_plots(density_plots, ncol = 2)
density_plots[[1]]   # Lymph
density_plots[[2]]   # LMR
density_plots[[3]]   # NLR
density_plots[[4]]   # PLR

###Correlation plot###
plot_correlation <- function(prs_blood_merge, prs1_col, prs2_col, group_name) {
  
  # 计算 Pearson correlation
  ct <- cor.test(
    prs_blood_merge[[prs1_col]],
    prs_blood_merge[[prs2_col]],
    method = "pearson"
  )
  
  r_val <- round(ct$estimate, 3)
  #p_val <- format.pval(ct$p.value, digits = 3)
  
  ggplot(
    prs_blood_merge,
    aes(x = .data[[prs1_col]], y = .data[[prs2_col]])
  ) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    annotate(
      "text",
      x = -2, y = Inf,
      label = paste0("r = ", r_val),
      hjust = 1.1, vjust = 1.1,
      size = 4
    ) +
    theme_bw() +
    labs(
      title = paste0(group_name, " PRS Correlation"),
      x = "My PRS",
      y = "Standardized PRS"
    )
}

corr_plots <- purrr::pmap(
  list(prs_map$prs1, prs_map$prs2, prs_map$group),
  ~ plot_correlation(prs_blood_merge, ..1, ..2, ..3)
)
library(patchwork)

corr_panel <- wrap_plots(corr_plots, ncol = 2)
ggsave(
  "PRS_correlation_panel.pdf",
  corr_panel,
  width = 10,
  height = 8
)

####Step 4（强烈推荐）：Density + Correlation 合成一个大 panel
full_panel <- wrap_plots(
  density_plots,
  corr_plots,
  ncol = 2
) +
  plot_annotation(
    tag_levels = "A",
    title = "Comparison of two Standardized Polygenic Risk Scores"
  )
ggsave(
  "PRS_density_correlation_panel.pdf",
  full_panel,
  width = 12,
  height = 14
)


