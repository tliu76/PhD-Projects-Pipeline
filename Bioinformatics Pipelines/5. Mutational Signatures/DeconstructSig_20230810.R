# deconstructSigs 
# https://github.com/raerose01/deconstructSigs
# https://jakeconway.github.io/docs/deconstructSigs/  

library(deconstructSigs)
library(data.table)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
library(pals)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("deconstructSigs", quietly = TRUE))
  install.packages("deconstructSigs",)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
install.packages("deconstructSigs", repos = "https://cran.rstudio.com/")
setwd("/Users/tanxinliu/Library/CloudStorage/OneDrive-UniversityofSouthernCalifornia/wgs/mutationalsig")
sample.mut.ref = fread("smoking_tumor_autosomal_MPF.txt", header = T)

sample.mut.ref %>% count(Indiv) %>% pull(Indiv) %>% dput() -> sampleid



# Convert to deconstructSigs input
names(sample.mut.ref)
sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref, 
                                sample.id = "Indiv", 
                                chr = "CHROM", 
                                pos = "POS", 
                                ref = "REF", 
                                alt = "ALT",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)
#> Creating a generic function for 'nchar' from package 'base' in package 'S4Vectors'

# signature analysis 
# https://www.rdocumentation.org/packages/deconstructSigs/versions/1.8.0/topics/whichSignatures
plot_function <- function(i){
  plot_dataset <- whichSignatures(tumor.ref = sigs.input, 
                                  signatures.ref = signatures.cosmic, 
                                  sample.id = i,
                                  contexts.needed = TRUE)
  plot_dataset 
}

sampleid %>% map(
  plot_function
) -> plot_function_output


# make pie graphs
for (i in 1:35){
  filename <-paste0(rownames(plot_function_output[[i]]$diff),"_pie.png")
  png(file = filename, bg = "white")
  makePie(plot_function_output[[i]])
  dev.off()
}


# make signatures graphs
for (i in 1:35){
  filename <-paste0(rownames(plot_function_output[[i]]$diff),"_plotSignatures.png")
  png(file = filename, bg = "white")
  deconstructSigs::plotSignatures(plot_function_output[[i]]) 
  dev.off()
}

# make dataframe
data_function <- function(i){
  deconvolute_df <- plot_function_output[[i]]$weights
  deconvolute_df$unknwon <- plot_function_output[[i]]$unknown 
  deconvolute_df$SSE <- sum(plot_function_output[[i]]$diff^2)
  deconvolute_df
}

c(1:35) %>% map(
  data_function
) -> deconvolute_df_output

deconvolute_df_output_summary <- deconvolute_df_output %>% bind_rows() 
deconvolute_df_output_summary <- cbind(rownames(deconvolute_df_output_summary),deconvolute_df_output_summary)
rownames(deconvolute_df_output_summary) <- NULL
colnames(deconvolute_df_output_summary) <- c(names(deconvolute_df_output_summary)) #to not write all the column names

colnames(deconvolute_df_output_summary)[1] <- "sampleid" 
names(deconvolute_df_output_summary)
#recode sample id 1:35
char_vector <- sapply(1:35, function(x) paste0(x))

id<- as.array(char_vector)
deconvolute_df_output_summary$id <- id

library(data.table)

setwd("/Users/tanxinliu/Library/CloudStorage/OneDrive-UniversityofSouthernCalifornia/wgs/mutationalsig")
deconvolute_df_output_summary <- fread("deconvolute_df_output_summary.txt" )

# deconvolute_df_output_summary %>% select(c("unknwon",  "Signature.30",  "Signature.26","Signature.25",  "Signature.18","Signature.16",  "Signature.12", "Signature.9",  "Signature.8", "Signature.3", "Signature.7","Signature.13", "Signature.2", "Signature.5","Signature.1")) %>% rownames_to_column("sampleid") %>% fwrite("/dir/deconstructSigs/cosmic_default.out/deconvolute_df_output_summary.txt")


deconvolute_df_output_summary %>% fwrite("deconvolute_df_output_summary.txt")

# plot out signature proportions  
#ggplot_df <- deconvolute_df_output_summary %>% rownames_to_column("sampleid") %>% pivot_longer(!c("sampleid", "SSE"), names_to = "signatures", values_to = "proportion") %>% filter(proportion!=0)
ggplot_df <- deconvolute_df_output_summary  %>% pivot_longer(!c("sampleid","id", "SSE"), names_to = "signatures", values_to = "proportion") %>% filter(proportion!=0)

ggplot_df %>% pull(signatures) %>% unique() %>% dput()

ggplot_df <- ggplot_df %>%
  mutate(sampleid = str_sub(sampleid, 7))

library(RColorBrewer)
nb.cols <- 13
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)


##percent:
png(file = "deconstructsig2024.png", bg = "white", width = 700, height = 500, units = "px")
ggplot(ggplot_df, aes(fill=factor(signatures, levels=c("Signature.1","Signature.5","Signature.2" ,"Signature.13","Signature.7","Signature.18",  "Signature.8","Signature.3",  
                                                       "Signature.9", "Signature.12","Signature.16",  "Signature.25", 
                                                       "unknwon")), y= proportion, x = reorder(id,as.numeric(id)) ))+ 
  geom_bar(position="fill", stat="identity")+ theme_light() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+scale_fill_manual(values=mycolors) + guides(fill = guide_legend(title = "Signature"))+labs(x = "\nSample id", title = "DeconstructSigs signature deconvolution")+
  labs(y = "Proportion")

dev.off()







# previous: overall percent
png(file = "stacked_bar.png", bg = "white", width = 700, height = 500, units = "px")

ggplot(ggplot_df, aes(fill=factor(signatures, levels=c("Signature.1", "Signature.3", "Signature.5", "Signature.8", 
                                                       "unknwon", "Signature.9", "Signature.18", "Signature.12", "Signature.25", 
                                                       "Signature.16", "Signature.13", "Signature.7", "Signature.2")), y= proportion, x=reorder(sampleid, sampleid))) + 
  geom_bar(position="fill", stat="identity")+ theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_fill_manual(values=as.vector(cols25(35))) + guides(fill = guide_legend(title = "Signature"))+labs(x = "sample id", title = "deconstructSigs signature deconvolution") 
dev.off()

###by smoking group
setwd("/Users/tanxinliu/Library/CloudStorage/OneDrive-UniversityofSouthernCalifornia/wgs")
#including smoking variable
library(readxl)
library(dplyr)
library(tidyr)
library(tidyverse)
smk <- read_excel("/Users/tanxinliu/Library/CloudStorage/OneDrive-UniversityofSouthernCalifornia/wgs/TRDRP_smoking_variable_4May23.xlsx", 
                  col_types = c("numeric", "text"))
smk <- smk %>% mutate( smok = 
                         case_when(
                           Tobacco_exposure_group =="High" ~ 1,
                           Tobacco_exposure_group=="Low" ~ 0
                         )  ) 



  new <- merge(x = ggplot_df, y = smk, by = c("sampleid"), all = TRUE)
new <- apply_labels(new, smok=c("High tobacco exposure"=1, "Low tobacco exposure"=0))
  
  png(file = "smokgroup_stacked_bar2024.png", bg = "white", width = 900, height = 500, units = "px")
  
  ggplot(new, aes(fill=factor(signatures, levels=c("Signature.1","Signature.5","Signature.2" ,"Signature.13","Signature.7","Signature.18",  "Signature.8","Signature.3",  
                                                   "Signature.9", "Signature.12","Signature.16",  "Signature.25", 
                                                   "unknwon")), y= proportion, x=reorder(id, as.numeric(id)))) + 
    geom_bar(position="fill", stat="identity")+ theme_light()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),strip.text = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14))+
  scale_fill_manual(values=mycolors)+ 
  guides(fill = guide_legend(title = "Signature"))+
  labs(x = "\nSample id", title = "DeconstructSigs signature deconvolution by tobacco exposure group")+
    labs(y = "Proportion\n")+
  facet_wrap(~smok,scales = "free")

dev.off()
  
  
  
# create a master spreadsheet to print out 
deconvolute_df_output_summary %>% select(c("Signature.1", "Signature.3", "Signature.5", "Signature.8", 
                                           "unknwon", "Signature.9", "Signature.18", "Signature.12", "Signature.25", 
                                           "Signature.16", "Signature.13", "Signature.7", "Signature.2")) %>% rownames_to_column("sampleid")
deconvolute_df_output_summary %>% fwrite("deconvolute_df_output_summary")
# plot out error rate  
png(file = "SSE.png", bg = "white", width = 800, height = 600, units = "px")
deconvolute_df_output_summary %>% rownames_to_column("sampleid") %>% 
  ggplot(aes(y= SSE, x=reorder(sampleid, sampleid), label = round(SSE, 4))) + 
  geom_point(color = "red")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +labs(x = "sample id", title = "deconstructSigs SSE") + geom_text(size = 0.1) + ggrepel::geom_text_repel()
dev.off()

