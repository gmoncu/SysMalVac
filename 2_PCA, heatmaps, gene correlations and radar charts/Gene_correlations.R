##################################################### Gene correlation analysis #####################################################

# This is the code to perform the gene correlation analysis with cumulative parasitemia and prepatency of the CPS study shown in the manuscript.

rm(list=ls())

# packages
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape2)


## Load log2 FC CSP/DMSO PfRBC/uRBC data of CPS immunized subjects
FC_CPS<- read.csv(file = "data/CPS_FC_microarray_data.csv")
summary(FC_CPS$stim)
summary(FC_CPS$timepoint)

# Select only PfRBC/uRBC post-immunization
FC_CPS_C1_iRBC <- FC_CPS[FC_CPS$stim=="iRBC/uRBC" & FC_CPS$timepoint=="CI",]
dim(FC_CPS_C1_iRBC) # 23 individuals


############################################################################################
####################### Correlations using GO terms ########################################
############################################################################################


############### Correlation with cumulative parasitemia ###############

## Read data
df <- read.csv("data/CPS_studydata.csv")
df_parasit_prepat<- read.csv2("data/CPS_parasit_prepatency.csv")

# modify and merge databases
donor <- as.character(df$donor)
df$donor <- sapply(strsplit(donor, split='_', fixed=TRUE), function(x) (x[1]))

df <- df[, c(3,11)]
df <- df[!duplicated(df), ]

df_parasit <- df_parasit_prepat[, 1:5]

df_parasit_merge <- merge(df_parasit, df, by="donor", x.all=TRUE)

# merge with gene data

FC_CPS_C1_iRBC$ind_prot <- as.character(FC_CPS_C1_iRBC$ind_prot)

s1 <- sapply(strsplit(FC_CPS_C1_iRBC$ind_prot, split='_', fixed=TRUE), function(x) (x[2]))
FC_CPS_C1_iRBC$malaria <- as.factor(s1)


FC_CPS_C1_iRBC$indn <- sapply(strsplit(FC_CPS_C1_iRBC$ind_prot, split='_', fixed=TRUE), function(x) (x[1]))

tm1_m <- merge(FC_CPS_C1_iRBC, df_parasit_merge, by="indn", x.all=TRUE)

head(tm1_m[, 9780:9795])
head(tm1_m[, 1:3])


# Spearman correlations
c1_all <- sapply(2:9786, function(x) cor(tm1_m[,x], tm1_m$cumulative, method = "spearman"))
c2_all <- sapply(2:9786, function(x) cor.test(tm1_m[,x], tm1_m$cumulative, method = "spearman")$p.value)
c3_all <- p.adjust(c2_all, method= "BH")
corr_C1_iRBC_all <- data.frame( GeneSymbol = colnames(tm1_m[,2:9786]),  r = unlist(c1_all), pval = unlist(c2_all), padj=unlist(c3_all), row.names = NULL)
dim(corr_C1_iRBC_all)
head(corr_C1_iRBC_all)

# add GO terms
anotation <- read.csv2(file ="data/RMAAFFY_complete_annotation.csv",  header = TRUE)
anotation <- anotation[,c(4,21)]
anotation <- anotation[!duplicated(anotation), ]
dim(anotation)

corr_C1_iRBC_all_m <- merge(corr_C1_iRBC_all, anotation, all.x=TRUE,by="GeneSymbol")
dim(corr_C1_iRBC_all_m)
View(corr_C1_iRBC_all_m)
corr_C1_iRBC_all_m_sign <- corr_C1_iRBC_all_m[corr_C1_iRBC_all_m$pval<=0.05,]
dim(corr_C1_iRBC_all_m_sign)

# save dataframes
write.table(corr_C1_iRBC_all_m, "corr_C1_iRBC_all_CumulativeParasitemia.csv", row.names = FALSE, dec = ".", sep = ";")
write.table(corr_C1_iRBC_all_m_sign, "corr_C1_iRBC_all_significant_CumulativeParasitemia.csv", row.names = FALSE, dec = ".", sep = ";")


# plot gene correlations
p <- ggplot(corr_C1_iRBC_all_m_sign, aes(x=reorder(GeneSymbol, -r), r)) + geom_bar(stat="identity") +
  ggtitle("Spearman cor of gene PfRBC/uRBC log2(FC) vs cumulative parasitemia at post-immunization (CI)") +
  xlab("Genes") + ylab("Rho") +
  theme(axis.title =element_text(size=15), axis.text.y = element_text( size=5), axis.text.x = element_text(size=10)) + 
  coord_flip()
p

pdf("Spearman correlations CI iRBC all genes.pdf", width = 9, height = 20)
p
dev.off()


# Select genes involved in immune response (based on GO terms)
corr_C1_iRBC_all_m_sign <- corr_C1_iRBC_all_m_sign %>%
  mutate(immune.term = factor(
    ifelse(str_detect(GO.Biological.Process.Term, 
                      'immune|B cell|T cell|interferon|cell proliferation|cell cycle|leukocyte|cytokine|antigen|interleukin|immunoglobulin|T-helper|chemotaxis|tumour necrosis factor|defense response|NF-KappaB|cytotoxicity|memory|antigenic|humoral|lymphocyte|natural killer|T cell'),
           'yes', 'no')))
    
summary(corr_C1_iRBC_all_m_sign$immune.term)
    
corr_C1_iRBC_all_m_sign_immune <- corr_C1_iRBC_all_m_sign[corr_C1_iRBC_all_m_sign$immune.term=="yes",]

write.csv2(corr_C1_iRBC_all_m_sign_immune , "corr_C1_iRBC_allsign_immuneterms.csv")


### the saved file "corr_C1_iRBC_allsign_immuneterms.csv" was externally modified to get separate GO terms in different columns: "corr_C1_iRBC_allsign_immuneterms_sep.csv"

corr_immuneterms_sep <- read.csv2( "corr_C1_iRBC_allsign_immuneterms_sep.csv", dec = ".")

corr_immuneterms_sep_long <- melt(corr_immuneterms_sep, id.vars = 1:4, na.rm=FALSE)
dim(corr_immuneterms_sep_long)
corr_immuneterms_sep_long$value <- str_trim(corr_immuneterms_sep_long$value)
terms.summary <- as.data.frame(table(corr_immuneterms_sep_long$value))

colnames(terms.summary)[1] <- "terms"
terms.summary$terms <- as.character(terms.summary$terms)

terms.summary <- terms.summary %>%
  mutate(immune.term = factor(
    ifelse(str_detect(terms, 
                      'immune|B cell|T cell|interferon|cell proliferation|cell cycle|leukocyte|cytokine|antigen|interleukin|immunoglobulin|T-helper|chemotaxis|tumour necrosis factor|defense response|NF-KappaB|cytotoxicity|memory|antigenic|humoral|lymphocyte|natural killer|T cell'),
           'yes', 'no')))

summary(terms.summary$immune.term)

write.csv2(terms.summary , "corr_C1_iRBC_allsign_immuneterms_summary.csv")


terms.summary.yes <- terms.summary[terms.summary$immune.term=="yes",]

# Plot for manuscript figure showing frequencies of GO terms related to immune response
p <- ggplot(terms.summary.yes, aes(x=terms, Freq)) + geom_bar(stat="identity") +
  ggtitle("GO terms realted to immune response PfRBC/uRBC log2(FC) vs cumulative parasitemia at post-immunization (CI) ") +
  xlab("GO Biological processes terms") + ylab("Frequency") +
  theme(plot.title = element_text(size=20), axis.title =element_text(size=20), axis.text.y = element_text( size=19), axis.text.x = element_text(size=15, angle = 90, hjust = 1))
p

pdf("GO terms immuno corr cumulative parasitemia CI iRBC.pdf", width = 18, height = 18)
p
dev.off()



############### Correlation with prepatency  ###############


df_prepat <- df_parasit_prepat[, c(1,7:8)]

# merge with df databases
df_prepat_merge <- merge(df_prepat, df, by="donor", x.all=TRUE)
dim(df_prepat_merge)

# merge with gene data
dim(tm1)
tm1_m <- merge(tm1, df_prepat_merge, by="indn", x.all=TRUE)
dim(tm1_m)

# Spearman correlations with TS
c1_all <- sapply(2:9786, function(x) cor(tm1_m[,x], tm1_m$TS, use="complete.obs", method = "spearman"))
c2_all <- sapply(2:9786, function(x) cor.test(tm1_m[,x], tm1_m$TS, use="complete.obs", method = "spearman")$p.value)
c3_all <- p.adjust(c2_all, method= "BH")
corr_C1_iRBC_all_TS <- data.frame( GeneSymbol = colnames(tm1_m[,2:9786]),  r = unlist(c1_all), pval = unlist(c2_all), padj=unlist(c3_all), row.names = NULL)
dim(corr_C1_iRBC_all_TS)
head(corr_C1_iRBC_all_TS)

# Spearman correlations with PCR
c1_all <- sapply(2:9786, function(x) cor(tm1_m[,x], tm1_m$PCR, use="complete.obs", method = "spearman"))
c2_all <- sapply(2:9786, function(x) cor.test(tm1_m[,x], tm1_m$PCR, use="complete.obs", method = "spearman")$p.value)
c3_all <- p.adjust(c2_all, method= "BH")
corr_C1_iRBC_all_PCR <- data.frame( GeneSymbol = colnames(tm1_m[,2:9786]),  r = unlist(c1_all), pval = unlist(c2_all), padj=unlist(c3_all), row.names = NULL)
dim(corr_C1_iRBC_all_PCR)
head(corr_C1_iRBC_all_PCR)

plot(tm1_m$UBC, tm1_m$PCR)
plot(tm1_m$IGSF8, tm1_m$PCR)
plot(tm1_m$PSD4, tm1_m$PCR)


#########  TS ######### 

corr_C1_iRBC_all_TS_m <- merge(corr_C1_iRBC_all_TS, anotation, by="GeneSymbol")
dim(corr_C1_iRBC_all_TS_m)
corr_C1_iRBC_all_TS_m_sign <- corr_C1_iRBC_all_TS_m[corr_C1_iRBC_all_TS_m$pval<=0.05,]
dim(corr_C1_iRBC_all_TS_m_sign)

write.table(corr_C1_iRBC_all_TS_m, "corr_C1_iRBC_all_TS.csv", row.names = FALSE, dec = ".", sep = ";")
write.table(corr_C1_iRBC_all_TS_m_sign, "corr_C1_iRBC_all_TS_significant.csv", row.names = FALSE, dec = ".", sep = ";")

# plot gene correlations
p <- ggplot(corr_C1_iRBC_all_TS_m_sign, aes(x=reorder(GeneSymbol, -r), r)) + geom_bar(stat="identity") +
  ggtitle("Spearman cor of gene PfRBC/uRBC log2(FC) vs TS prepatency at post-immunization (CI) ") +
  xlab("Genes") + ylab("Rho") +
  theme(axis.title =element_text(size=15), axis.text.y = element_text( size=5), axis.text.x = element_text(size=10)) + 
  coord_flip()
p

pdf("Spearman correlations TS vs iRBC all genes.pdf", width = 9, height = 25)
p
dev.off()

# Select genes with GO terms related to immune response
corr_C1_iRBC_all_TS_m_sign <- corr_C1_iRBC_all_TS_m_sign %>%
  mutate(immune.term = factor(
    ifelse(str_detect(GO.Biological.Process.Term, 
                      'immune|B cell|T cell|interferon|cell proliferation|cell cycle|leukocyte|cytokine|antigen|interleukin|immunoglobulin|T-helper|chemotaxis|tumour necrosis factor|defense response|NF-KappaB|cytotoxicity|memory|antigenic|humoral|lymphocyte|natural killer|T cell'),
           'yes', 'no')))

summary(corr_C1_iRBC_all_TS_m_sign$immune.term)

corr_C1_iRBC_all_TS_m_sign_immune <- corr_C1_iRBC_all_TS_m_sign[corr_C1_iRBC_all_TS_m_sign$immune.term=="yes",]

write.csv2(corr_C1_iRBC_all_TS_m_sign_immune , "corr_C1_iRBC_allsign_TS_immuneterms.csv")


### the saved file "corr_C1_iRBC_allsign_TS_immuneterms.csv" was externally modified to get separate GO terms in different columns: "corr_C1_iRBC_allsign_TS_immuneterms_sep.csv"

corr_immuneterms_TS_sep <- read.csv2( "corr_C1_iRBC_allsign_TS_immuneterms_sep.csv", dec = ".")

corr_immuneterms_TS_sep_long <- melt(corr_immuneterms_TS_sep, id.vars = 1:4, na.rm=FALSE)
dim(corr_immuneterms_TS_sep_long)
corr_immuneterms_TS_sep_long$value <- str_trim(corr_immuneterms_TS_sep_long$value)
terms.summary <- as.data.frame(table(corr_immuneterms_TS_sep_long$value))

colnames(terms.summary)[1] <- "terms"
terms.summary$terms <- as.character(terms.summary$terms)

terms.summary <- terms.summary %>%
  mutate(immune.term = factor(
    ifelse(str_detect(terms, 
                      'immune|B cell|T cell|interferon|cell proliferation|cell cycle|leukocyte|cytokine|antigen|interleukin|immunoglobulin|T-helper|chemotaxis|tumour necrosis factor|defense response|NF-KappaB|cytotoxicity|memory|antigenic|humoral|lymphocyte|natural killer|T cell'),
           'yes', 'no')))

summary(terms.summary$immune.term)

write.csv2(terms.summary , "corr_C1_iRBC_allsign_TS_immuneterms_summary.csv")


terms.summary.yes <- terms.summary[terms.summary$immune.term=="yes",]

# plot frquency of GO terms related to immune response
p <- ggplot(terms.summary.yes, aes(x=terms, Freq)) + geom_bar(stat="identity") +
  ggtitle("Frequency of genes in GO terms realted to immune response PfRBC/uRBC log2(FC) vs TS prepatency at post-immunization (CI) ") +
  xlab("GO Biological processes terms") + ylab("Frequency") +
  theme(plot.title = element_text(size=20), axis.title =element_text(size=20), axis.text.y = element_text( size=19), axis.text.x = element_text(size=15, angle = 90, hjust = 1))
p

pdf("GO terms immuno corr TS prepatency CI iRBC.pdf", width = 25, height = 18)
p
dev.off()

############ 


######  PCR ###### 

corr_C1_iRBC_all_PCR_m <- merge(corr_C1_iRBC_all_PCR, anotation, by="GeneSymbol")
dim(corr_C1_iRBC_all_PCR_m)
corr_C1_iRBC_all_PCR_m_sign <- corr_C1_iRBC_all_PCR_m[corr_C1_iRBC_all_PCR_m$pval<=0.05,]
dim(corr_C1_iRBC_all_PCR_m_sign)

write.table(corr_C1_iRBC_all_PCR_m, "corr_C1_iRBC_all_PCR.csv", row.names = FALSE, dec = ".", sep = ";")
write.table(corr_C1_iRBC_all_PCR_m_sign, "corr_C1_iRBC_all_PCR_significant.csv", row.names = FALSE, dec = ".", sep = ";")

# plot gene correlations
p <- ggplot(corr_C1_iRBC_all_PCR_m_sign, aes(x=reorder(GeneSymbol, -r), r)) + geom_bar(stat="identity") +
  ggtitle("Spearman cor of gene PfRBC/uRBC log2(FC) vs PCR prepatency at post-immunization (CI) ") +
  xlab("Genes") + ylab("Rho") +
  theme(axis.title =element_text(size=15), axis.text.y = element_text( size=5), axis.text.x = element_text(size=10)) + 
  coord_flip()
p

pdf("Spearman correlations PCR vs iRBC all genes.pdf", width = 9, height = 35)
p
dev.off()


# select genes with GO terms related with immune response
corr_C1_iRBC_all_PCR_m_sign <- corr_C1_iRBC_all_PCR_m_sign %>%
  mutate(immune.term = factor(
    ifelse(str_detect(GO.Biological.Process.Term, 
                      'immune|B cell|T cell|interferon|cell proliferation|cell cycle|leukocyte|cytokine|antigen|interleukin|immunoglobulin|T-helper|chemotaxis|tumour necrosis factor|defense response|NF-KappaB|cytotoxicity|memory|antigenic|humoral|lymphocyte|natural killer|T cell'),
           'yes', 'no')))

summary(corr_C1_iRBC_all_PCR_m_sign$immune.term)

corr_C1_iRBC_all_PCR_m_sign_immune <- corr_C1_iRBC_all_PCR_m_sign[corr_C1_iRBC_all_PCR_m_sign$immune.term=="yes",]

write.csv2(corr_C1_iRBC_all_PCR_m_sign_immune , "corr_C1_iRBC_allsign_PCR_immuneterms.csv")

### the saved file "corr_C1_iRBC_allsign_PCR_immuneterms.csv" was externally modified to get separate GO terms in different columns: "corr_C1_iRBC_allsign_PCR_immuneterms_sep.csv"

corr_immuneterms_PCR_sep <- read.csv2( "corr_C1_iRBC_allsign_PCR_immuneterms_sep.csv", dec = ".")

corr_immuneterms_PCR_sep_long <- melt(corr_immuneterms_PCR_sep, id.vars = 1:4, na.rm=FALSE)
dim(corr_immuneterms_PCR_sep_long)
corr_immuneterms_PCR_sep_long$value <- str_trim(corr_immuneterms_PCR_sep_long$value)
terms.summary <- as.data.frame(table(corr_immuneterms_PCR_sep_long$value))

colnames(terms.summary)[1] <- "terms"
terms.summary$terms <- as.character(terms.summary$terms)

terms.summary <- terms.summary %>%
  mutate(immune.term = factor(
    ifelse(str_detect(terms, 
                      'immune|B cell|T cell|interferon|cell proliferation|cell cycle|leukocyte|cytokine|antigen|interleukin|immunoglobulin|T-helper|chemotaxis|tumour necrosis factor|defense response|NF-KappaB|cytotoxicity|memory|antigenic|humoral|lymphocyte|natural killer|T cell'),
           'yes', 'no')))

summary(terms.summary$immune.term)

write.csv2(terms.summary , "corr_C1_iRBC_allsign_PCR_immuneterms_summary.csv")


terms.summary.yes <- terms.summary[terms.summary$immune.term=="yes",]


# plot frequency of GO terms related to immune response
p <- ggplot(terms.summary.yes, aes(x=terms, Freq)) + geom_bar(stat="identity") +
  ggtitle("Frequency of genes in GO terms realted to immune response PfRBC/uRBC log2(FC) vs PCR prepatency at post-immunization (CI) ") +
  xlab("GO Biological processes terms") + ylab("Frequency") +
  theme(plot.title = element_text(size=20), axis.title =element_text(size=20), axis.text.y = element_text( size=19), axis.text.x = element_text(size=15, angle = 90, hjust = 1))
p

pdf("GO terms immuno corr PCR prepatency CI iRBC.pdf", width = 25, height = 18)
p
dev.off()





############################################################################################
####################### Correlations using BTM  ############################################
############################################################################################

############### Correlation with cumulative parasitemia ###############

# BTM info
anotation <- read.csv2(file ="data/btm_annotation_table.csv",  header = TRUE)
anotation <- anotation[,c(1:3,5,8)] 

anotation$Module.member.genes <- as.character(anotation$Module.member.genes)

corr_C1_iRBC_all$GeneSymbol <- as.character(corr_C1_iRBC_all$GeneSymbol)
corr_C1_iRBC_all$BTM <- NA

for (i in 1:length(corr_C1_iRBC_all$GeneSymbol)) {
  
  listbtm <-list()
  
  for (j in 1:length(anotation$Module.member.genes)) {
    
    if (str_detect(anotation$Module.member.genes[j],  corr_C1_iRBC_all$GeneSymbol[i])) {
      listbtm[[length(listbtm)+1]] = as.character(anotation$Composite.name[j])
    } else {
      listbtm[[length(listbtm)+1]] = "NA"}
    
    corr_C1_iRBC_all$BTM[i] <- paste0(unlist(listbtm),collapse="/")
    
  }
}

corr_C1_iRBC_all$BTM <- gsub("NA/", "", corr_C1_iRBC_all$BTM)
corr_C1_iRBC_all$BTM <- gsub("/NA", "", corr_C1_iRBC_all$BTM)

corr_C1_iRBC_all_sign <- corr_C1_iRBC_all[corr_C1_iRBC_all$pval<=0.05,]


write.table(corr_C1_iRBC_all, "corr_C1_iRBC_all_CumulativeParasitemia_btm.csv", row.names = FALSE, dec = ".", sep = ";")
write.table(corr_C1_iRBC_all_sign, "corr_C1_iRBC_all_significant_CumulativeParasitemia_btm.csv", row.names = FALSE, dec = ".", sep = ";")

### the saved file "corr_C1_iRBC_all_significant_CumulativeParasitemia_btm.csv" was externally modified to get separate GO terms in different columns: "corr_C1_iRBC_all_significant_CumulativeParasitemia_btm_sep.csv""

corr_btm_sep <- read.csv2( "corr_C1_iRBC_all_significant_CumulativeParasitemia_btm_sep.csv", dec = ".")


corr_btm_sep_long <- melt(corr_btm_sep, id.vars = 1:4, na.rm=FALSE) 
dim(corr_btm_sep_long)
btm.summary <- as.data.frame(table(corr_btm_sep_long$value))

colnames(btm.summary)[1] <- "BTM"

btm.summary <- btm.summary[-1,]

write.csv2(btm.summary , "corr_C1_iRBC_allsign_CumulativeParasitemia_BTM_summary.csv")


p <- ggplot(btm.summary, aes(x=BTM, Freq)) + geom_bar(stat="identity") +
  ggtitle("BTM of correlations PfRBC/uRBC log2(FC) vs cumulative parasitemia at post-immunization (CI) ") +
  xlab("BTM") + ylab("Frequency") +
  theme(plot.title = element_text(size=20), axis.title =element_text(size=20), axis.text.y = element_text( size=19), axis.text.x = element_text(size=15, angle = 90, hjust = 1))
p

pdf("BTM corr cumulative parasitemia CI iRBC.pdf", width = 20, height = 18)
p
dev.off()


############### Correlation with prepatency  ###############

### TS ###
corr_C1_iRBC_all_TS$GeneSymbol <- as.character(corr_C1_iRBC_all_TS$GeneSymbol)

corr_C1_iRBC_all_TS$BTM <- NA

for (i in 1:length(corr_C1_iRBC_all_TS$GeneSymbol)) {
  
  listbtm <-list()
  
  for (j in 1:length(anotation$Module.member.genes)) {
    
    if (str_detect(anotation$Module.member.genes[j],  corr_C1_iRBC_all_TS$GeneSymbol[i])) {
      listbtm[[length(listbtm)+1]] = as.character(anotation$Composite.name[j])
    } else {
      listbtm[[length(listbtm)+1]] = "NA"}
    
    corr_C1_iRBC_all_TS$BTM[i] <- paste0(unlist(listbtm),collapse="/")
    
  }
}

corr_C1_iRBC_all_TS$BTM <- gsub("NA/", "", corr_C1_iRBC_all_TS$BTM)
corr_C1_iRBC_all_TS$BTM <- gsub("/NA", "", corr_C1_iRBC_all_TS$BTM)

corr_C1_iRBC_all_TS_sign <- corr_C1_iRBC_all_TS[corr_C1_iRBC_all_TS$pval<=0.05,]


write.table(corr_C1_iRBC_all_TS, "corr_C1_iRBC_all_TS_btm.csv", row.names = FALSE, dec = ".", sep = ";")
write.table(corr_C1_iRBC_all_TS_sign, "corr_C1_iRBC_all_significant_TS_btm.csv", row.names = FALSE, dec = ".", sep = ";")

### the saved file "corr_C1_iRBC_all_significant_TS_btm.csv" was externally modified to get separate GO terms in different columns: "corr_C1_iRBC_all_significant_TS_btm_sep.csv""

corr_btm_sep <- read.csv2( "corr_C1_iRBC_all_significant_TS_btm_sep.csv", dec = ".")


corr_btm_sep_long <- melt(corr_btm_sep, id.vars = 1:4, na.rm=FALSE) 
dim(corr_btm_sep_long)
btm.summary <- as.data.frame(table(corr_btm_sep_long$value))

colnames(btm.summary)[1] <- "BTM"

btm.summary <- btm.summary[-1,]

write.table(btm.summary, "corr_C1_iRBC_allsign_TS_BTM_summary.csv", row.names = FALSE, dec = ".", sep = ";")

p <- ggplot(btm.summary, aes(x=BTM, Freq)) + geom_bar(stat="identity") +
  ggtitle("BTM of correlations PfRBC/uRBC log2(FC) vs TS prepatency at post-immunization (CI) ") +
  xlab("BTM") + ylab("Frequency") +
  theme(plot.title = element_text(size=20), axis.title =element_text(size=20), axis.text.y = element_text( size=19), axis.text.x = element_text(size=15, angle = 90, hjust = 1))
p

pdf("BTM corr TS prepatency CI iRBC.pdf", width = 18, height = 18)
p
dev.off()



### PCR ###
corr_C1_iRBC_all_PCR$GeneSymbol <- as.character(corr_C1_iRBC_all_PCR$GeneSymbol)
corr_C1_iRBC_all_PCR$BTM <- NA

for (i in 1:length(corr_C1_iRBC_all_PCR$GeneSymbol)) {
  
  listbtm <-list()
  
  for (j in 1:length(anotation$Module.member.genes)) {
    
    if (str_detect(anotation$Module.member.genes[j],  corr_C1_iRBC_all_PCR$GeneSymbol[i])) {
      listbtm[[length(listbtm)+1]] = as.character(anotation$Composite.name[j])
    } else {
      listbtm[[length(listbtm)+1]] = "NA"}
    
    corr_C1_iRBC_all_PCR$BTM[i] <- paste0(unlist(listbtm),collapse="/")
    
  }
}

corr_C1_iRBC_all_PCR$BTM <- gsub("NA/", "", corr_C1_iRBC_all_PCR$BTM)
corr_C1_iRBC_all_PCR$BTM <- gsub("/NA", "", corr_C1_iRBC_all_PCR$BTM)

corr_C1_iRBC_all_PCR_sign <- corr_C1_iRBC_all_PCR[corr_C1_iRBC_all_PCR$pval<=0.05,]


write.table(corr_C1_iRBC_all_PCR, "corr_C1_iRBC_all_PCR_btm.csv", row.names = FALSE, dec = ".", sep = ";")
write.table(corr_C1_iRBC_all_PCR_sign, "corr_C1_iRBC_all_significant_PCR_btm.csv", row.names = FALSE, dec = ".", sep = ";")


### the saved file "corr_C1_iRBC_all_significant_PCR_btm.csv" was externally modified to get separate GO terms in different columns: "corr_C1_iRBC_all_significant_PCR_btm_sep.csv""

corr_btm_sep <- read.csv2( "corr_C1_iRBC_all_significant_PCR_btm_sep.csv", dec = ".")

corr_btm_sep_long <- melt(corr_btm_sep, id.vars = 1:4, na.rm=FALSE) 
dim(corr_btm_sep_long)
btm.summary <- as.data.frame(table(corr_btm_sep_long$value))

colnames(btm.summary)[1] <- "BTM"

btm.summary <- btm.summary[-1,]

write.table(btm.summary, "corr_C1_iRBC_allsign_PCR_BTM_summary.csv", row.names = FALSE, dec = ".", sep = ";")

p <- ggplot(btm.summary, aes(x=BTM, Freq)) + geom_bar(stat="identity") +
  ggtitle("BTM of correlations PfRBC/uRBC log2(FC) vs PCR prepatency at post-immunization (CI) ") +
  xlab("BTM") + ylab("Frequency") +
  theme(plot.title = element_text(size=20), axis.title =element_text(size=20), axis.text.y = element_text( size=19), axis.text.x = element_text(size=15, angle = 90, hjust = 1))
p

pdf("BTM corr PCR prepatency CI iRBC.pdf", width = 25, height = 18)
p
dev.off()




####### Merge anotations GO term and BTM for correlations ########

# read files correlations with Go terms
cp <- read.csv2( "corr_C1_iRBC_CumulativeParasitemia.csv", dec = ".")
TS <- read.csv2( "corr_C1_iRBC_TS.csv", dec = ".")
PCR <- read.csv2( "corr_C1_iRBC_PCR.csv", dec = ".")
dim(cp)
dim(TS)
dim(PCR)

cp <- cp[,c(1,5)]
TS <- TS[,c(1,5)]
PCR <- PCR[,c(1,5)]

# read files correlations with BTM
cp_2 <- read.csv2( "corr_C1_iRBC_all_CumulativeParasitemia_btm.csv", dec = ".")
TS_2 <- read.csv2( "corr_C1_iRBC_all_TS_btm.csv", dec = ".")
PCR_2 <- read.csv2( "corr_C1_iRBC_all_PCR_btm.csv", dec = ".")
dim(cp_2)
dim(TS_2)
dim(PCR_2)

# merge with corr dataframes with BTM info
cp_m <- merge(cp, cp_2, all=TRUE, by="GeneSymbol")
TS_m <- merge(TS, TS_2, all=TRUE, by="GeneSymbol")
PCR_m <- merge(PCR, PCR_2, all=TRUE, by="GeneSymbol")
dim(cp_m)
dim(TS_m)
dim(PCR_m)

cp_m <- cp_m[,c(1,3,4,5,2, 6)]
TS_m <- TS_m[,c(1,3,4,5,2, 6)]
PCR_m <- PCR_m[,c(1,3,4,5,2, 6)]

# Save dataframes for Manuscript Dataset S2
write.table(cp_m, "corr_C1_iRBC_CumulativeParasitemia_btm&GO.csv", row.names = FALSE, dec = ".", sep = ";")
write.table(TS_m, "corr_C1_iRBC_TS_btm&GO.csv", row.names = FALSE, dec = ".", sep = ";")
write.table(PCR_m, "corr_C1_iRBC_PCR_btm&GO.csv", row.names = FALSE, dec = ".", sep = ";")




#### Plots of correlations with genes with p-adjus <0.1 highlighted
# Not shown in the manuscript

cp_m_sign <- cp_m[cp_m$pval<=0.05,]
TS_m_sign <- TS_m[TS_m$pval<=0.05,]
PCR_m_sign <- PCR_m[PCR_m$pval<=0.05,]
dim(cp_m_sign)
dim(TS_m_sign)
dim(PCR_m_sign)

p <- ggplot(cp_m_sign, aes(x=reorder(GeneSymbol, -r), r)) + geom_bar(stat="identity" , aes(fill = padj<0.1), position = 'dodge' ) +
  scale_fill_manual(values = c("darkgrey", "dodgerblue4")) +
  ggtitle("Spearman cor of gene PfRBC/uRBC log2(FC) vs cumulative parasitemia at post-immunization (CI) ") +
  xlab("Genes") + ylab("Rho") +
  theme(axis.title =element_text(size=15), axis.text.y = element_text( size=5), axis.text.x = element_text(size=10)) + 
  coord_flip()
p

pdf("Spearman correlations CI iRBC Cumumulative Parasitemia signif.pdf", width = 9, height = 27)
p
dev.off()


p <- ggplot(TS_m_sign, aes(x=reorder(GeneSymbol, -r), r)) + geom_bar(stat="identity" , aes(fill= padj<0.1), position = 'dodge' ) +
  scale_fill_manual(values = c("darkgrey", "dodgerblue4")) +
  ggtitle("Spearman cor of gene PfRBC/uRBC log2(FC) vs TS at post-immunization (CI) ") +
  xlab("Genes") + ylab("Rho") +
  theme(axis.title =element_text(size=15), axis.text.y = element_text( size=5), axis.text.x = element_text(size=10)) + 
  coord_flip()
p

pdf("Spearman correlations CI iRBC TS signif.pdf", width = 9, height = 27)
p
dev.off()


p <- ggplot(PCR_m_sign, aes(x=reorder(GeneSymbol, -r), r)) + geom_bar(stat="identity" , aes(fill = padj<0.1), position = 'dodge' ) +
  scale_fill_manual(values = c("darkgrey", "dodgerblue4")) +
  ggtitle("Spearman cor of gene PfRBC/uRBC log2(FC) vs PCR at post-immunization (CI) ") +
  xlab("Genes") + ylab("Rho") +
  theme(axis.title =element_text(size=15), axis.text.y = element_text( size=5), axis.text.x = element_text(size=10)) + 
  coord_flip()
p

pdf("Spearman correlations CI iRBC PCR signif.pdf", width = 9, height = 40)
p
dev.off()

