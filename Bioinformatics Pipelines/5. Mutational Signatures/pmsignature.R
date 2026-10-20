# https://github.com/friend1ws/pmsignature

library(tidyverse)
library(data.table)
library(pmsignature)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(devtools)
devtools::install_github("friend1ws/pmsignature")
library(pmsignature)


inputFile<-paste0("smoking_tumor_autosomal_MPF.txt")

G <- readMPFile(inputFile, numBases = 5, trDir = TRUE, 
                bs_genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                txdb_transcript = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene)

set.seed(2021)
Param2 <- getPMSignature(G, K = 2, numInit=100)
Param3 <- getPMSignature(G, K = 3, numInit=100)
Param4 <- getPMSignature(G, K = 4, numInit=100)
Param5 <- getPMSignature(G, K = 5, numInit=100)
Param6 <- getPMSignature(G, K = 6, numInit=100)
Param7 <- getPMSignature(G, K = 7, numInit=100)
Param8 <- getPMSignature(G, K = 8, numInit=100)


set.seed(2021)
bootParam2 <- bootPMSignature(G, Param2, bootNum = 100)
bootParam3 <- bootPMSignature(G, Param3, bootNum = 100)
bootParam4 <- bootPMSignature(G, Param4, bootNum = 100)
bootParam5 <- bootPMSignature(G, Param5, bootNum = 100)
bootParam6 <- bootPMSignature(G, Param6, bootNum = 100)
bootParam7 <- bootPMSignature(G, Param7, bootNum = 100)
bootParam8 <- bootPMSignature(G, Param8, bootNum = 100)


# improved correlation matrix
library(corrplot)

cor2 <- round(cor(getMembershipValue(Param2)),
              digits = 2 # rounded to 2 decimals
) 
cor3 <- round(cor(getMembershipValue(Param3)),
              digits = 2 # rounded to 2 decimals
) 
cor4 <- round(cor(getMembershipValue(Param4)),
              digits = 2 # rounded to 2 decimals
) 
cor5 <- round(cor(getMembershipValue(Param5)),
              digits = 2 # rounded to 2 decimals
) 
cor6 <- round(cor(getMembershipValue(Param6)),
              digits = 2 # rounded to 2 decimals
) 

cor7 <- round(cor(getMembershipValue(Param7)),
              digits = 2 # rounded to 2 decimals
) 
cor8 <- round(cor(getMembershipValue(Param8)),
              digits = 2 # rounded to 2 decimals
) 


# create dataset for loglikelihood, bootstrap error and correlation
df_loglikelihood <- data.frame(K=c(2,3,4,5,6,7,8), loglikelihood = c(Param2@loglikelihood,Param3@loglikelihood,Param4@loglikelihood,Param5@loglikelihood,Param6@loglikelihood, Param7@loglikelihood, Param8@loglikelihood), bootstraperror =  c(mean(bootParam2[[1]]), mean(bootParam3[[1]]), mean(bootParam4[[1]]), mean(bootParam5[[1]]), mean(bootParam6[[1]]), mean(bootParam7[[1]]),mean(bootParam8[[1]])), max.cor = c(max(cor2[cor2!=max(cor2)]), max(cor3[cor3!=max(cor3)]), max(cor4[cor4!=max(cor4)]), max(cor5[cor5!=max(cor5)]), max(cor6[cor6!=max(cor6)]), max(cor7[cor7!=max(cor7)]), max(cor8[cor8!=max(cor8)])))


# df_loglikelihood <- data.frame(K=c(2,3,4,5,6,7,8), loglikelihood = c(Param2@loglikelihood,Param3@loglikelihood,Param4@loglikelihood,Param5@loglikelihood,Param6@loglikelihood, Param7@loglikelihood, Param8@loglikelihood), bootstraperror =  c(sum(bootParam2[[1]]^2), sum(bootParam3[[1]]^2), sum(bootParam4[[1]]^2), sum(bootParam5[[1]]^2), sum(bootParam6[[1]]^2), sum(bootParam7[[1]]^2),sum(bootParam8[[1]]^2)), max.cor = c(max(cor2[cor2!=max(cor2)]), max(cor3[cor3!=max(cor3)]), max(cor4[cor4!=max(cor4)]), max(cor5[cor5!=max(cor5)]), max(cor6[cor6!=max(cor6)]), max(cor7[cor7!=max(cor7)]), max(cor8[cor8!=max(cor8)])))


# df_loglikelihood %>% fwrite("df_loglikelihood_paramter.txt")

# create plots to decide the best number of signatures
df_loglikelihood %>% ggplot(aes(x = K, y = loglikelihood)) + geom_line()
ggsave("loglikelihood.png", width = 6,   height = 4)

df_loglikelihood %>% ggplot(aes(x = K, y = bootstraperror)) + geom_line()
ggsave("bootstraperror.png", width = 6,   height = 4)

df_loglikelihood %>% ggplot(aes(x = K, y = max.cor)) + geom_line()
ggsave("max.cor.png", width = 6,   height = 4)


corrplot(cor(getMembershipValue(Param2)),
         method = "number",
         type = "upper" # show only upper side
)

corrplot(cor(getMembershipValue(Param3)),
         method = "number",
         type = "upper" # show only upper side
)
corrplot(cor(getMembershipValue(Param4)),
         method = "number",
         type = "upper" # show only upper side
)
corrplot(cor(getMembershipValue(Param5)),
         method = "number",
         type = "upper" # show only upper side
)
corrplot(cor(getMembershipValue(Param6)),
         method = "number",
         type = "upper" # show only upper side
)

corrplot(cor(getMembershipValue(Param7)),
         method = "number",
         type = "upper" # show only upper side
)
corrplot(cor(getMembershipValue(Param8)),
         method = "number",
         type = "upper" # show only upper side
)



# visualize signature

visPMSignature(Param4, 1)
ggsave("k4.1.png", width = 6, height = 4)
visPMSignature(Param4, 2)
ggsave("k4.2.png", width = 6, height = 4)
visPMSignature(Param4, 3)
ggsave("k4.3.png", width = 6, height = 4)
visPMSignature(Param4, 4)
ggsave("k4.4.png", width = 6, height = 4)

visPMSignature(Param2, 1)
ggsave("k2.1.png", width = 6, height = 4)
visPMSignature(Param2, 2)
ggsave("k2.2.png", width = 6, height = 4)


visPMSignature(Param3, 1)
ggsave("k3.1.png", width = 6, height = 4)
visPMSignature(Param3, 2)
ggsave("k3.2.png", width = 6, height = 4)
visPMSignature(Param3, 3)
ggsave("k3.3.png", width = 6, height = 4)

visPMSignature(Param5, 1)
ggsave("k5.1.png", width = 6, height = 4)
visPMSignature(Param5, 2)
ggsave("k5.2.png", width = 6, height = 4)
visPMSignature(Param5, 3)
ggsave("k5.3.png", width = 6, height = 4)
visPMSignature(Param5, 4)
ggsave("k5.4.png", width = 6, height = 4)
visPMSignature(Param5, 5)
ggsave("k5.5.png", width = 6, height = 4)

visPMSignature(Param6, 1)
ggsave("k6.1.png", width = 6, height = 4)
visPMSignature(Param6, 2)
ggsave("k6.2.png", width = 6, height = 4)
visPMSignature(Param6, 3)
ggsave("k6.3.png", width = 6, height = 4)
visPMSignature(Param6, 4)
ggsave("k6.4.png", width = 6, height = 4)
visPMSignature(Param6, 5)
ggsave("k6.5.png", width = 6, height = 4)
visPMSignature(Param6, 6)
ggsave("k6.6.png", width = 6, height = 4)


visMembership(G, Param4)
ggsave("Membershipk4.png", width = 6, height = 4)
visMembership(G, Param3)
ggsave("Membershipk3.png", width = 6, height = 4)
visMembership(G, Param2)
ggsave("Membershipk2.png", width = 6, height = 4)
visMembership(G, Param5)
ggsave("Membershipk5.png", width = 6, height = 4)
visMembership(G, Param6)
ggsave("Membershipk6.png", width = 6, height = 4)



getSignatureValue(Param4, 1) %>% fwrite("Param4.sig1.csv")
getSignatureValue(Param4, 2) %>% fwrite("Param4.sig2.csv")
getSignatureValue(Param4, 3) %>% fwrite("Param4.sig3.csv")
getSignatureValue(Param4, 4) %>% fwrite("Param4.sig4.csv")

getSignatureValue(Param5, 1) %>% fwrite("Param5.sig1.csv")
getSignatureValue(Param5, 2) %>% fwrite("Param5.sig2.csv")
getSignatureValue(Param5, 3) %>% fwrite("Param5.sig3.csv")
getSignatureValue(Param5, 4) %>% fwrite("Param5.sig4.csv")
getSignatureValue(Param5, 5) %>% fwrite("Param5.sig5.csv")

getSignatureValue(Param3, 1) %>% fwrite("Param3.sig1.csv")
getSignatureValue(Param3, 2) %>% fwrite("Param3.sig2.csv")
getSignatureValue(Param3, 3) %>% fwrite("Param3.sig3.csv")

getSignatureValue(Param2, 1) %>% fwrite("Param2.sig1.csv")
getSignatureValue(Param2, 2) %>% fwrite("Param2.sig2.csv")

getSignatureValue(Param6, 1) %>% fwrite("Param6.sig1.csv")
getSignatureValue(Param6, 2) %>% fwrite("Param6.sig2.csv")
getSignatureValue(Param6, 3) %>% fwrite("Param6.sig3.csv")
getSignatureValue(Param6, 4) %>% fwrite("Param6.sig4.csv")
getSignatureValue(Param6, 5) %>% fwrite("Param6.sig5.csv")
getSignatureValue(Param6, 6) %>% fwrite("Param6.sig6.csv")

getMembershipValue(Param2) %>% fwrite("MembershipValuek2.csv",row.names=T)
getMembershipValue(Param3) %>% fwrite("MembershipValuek3.csv",row.names=T)
getMembershipValue(Param4) %>% fwrite("MembershipValuek4.csv",row.names=T)
getMembershipValue(Param5) %>% fwrite("MembershipValuek5.csv",row.names=T)
getMembershipValue(Param6) %>% fwrite("MembershipValuek6.csv",row.names=T)

# check correlation of different individuals
getMembershipValue(Param4) %>% ggplot(aes(x = signature_1,  y = signature_2)) + geom_point()





