#################################### Heatmaps #################################################################

# this is the code to perform the Heatmaps shown in the manuscript figures

#packages
rm(list=ls())
library(dplyr)
library(RColorBrewer)
library(heatmap3)


##################################################### CPS #####################################################

# Load FC CSP/DMSO and PfRBC/uRBC microarray data
all_CPS <- read.csv("data/CPS_FC_microarray_data.csv")


# read study data
demog <- read.csv("data/CPS_studydata.csv")

head(all_CPS[,9785:9788])
summary(all_CPS[,9785:9788])

all_CPS$ind_prot <- as.character(all_CPS$ind_prot)

all_CPS$indn <- sapply(strsplit(all_CPS$ind_prot, split='_', fixed=TRUE), function(x) (x[1]))
all_CPS$protected <- sapply(strsplit(all_CPS$ind_prot, split='_', fixed=TRUE), function(x) (x[2]))
all_CPS$indn
all_CPS$protected

head(all_CPS[,9785:9790])

sample_plate_indn <- demog[,c(2,4,8,11)]
sample_plate_indn <- sample_plate_indn[!duplicated(sample_plate_indn), ]

# merge
dim(all_CPS)
all_CPS<- merge(all_CPS, sample_plate_indn, by="indn")
dim(all_CPS)
summary(all_CPS$Sex)

table(all_CPS$group, all_CPS$timepoint)
table(all_CPS$Sex, all_CPS$group)

m <- all_CPS[,2:9786]
m1<- sapply(m, as.numeric)
m1 <- as.matrix(m1)
m2 <- t(m1)


# color for variables 

all_CPS$protected <- as.factor(all_CPS$protected)
all_CPS$timepoint <- as.factor(all_CPS$timepoint)
all_CPS$stim <- as.factor(all_CPS$stim)
all_CPS$indn <- as.factor(all_CPS$indn)
summary(all_CPS$indn)

Dose <- brewer.pal(3, "Set2")[1:3][all_CPS$group]
Stimulation <- brewer.pal(3, "Set2")[1:2][all_CPS$stim]
Protection <- brewer.pal(3, "Set2")[1:2][as.factor(all_CPS$protected)]
Plate <- brewer.pal(3, "Set2")[1:2][all_CPS$plate]
Timepoint <- brewer.pal(3, "Set2")[1:2][all_CPS$timepoint]
Sex <- brewer.pal(3, "Set2")[1:2][all_CPS$Sex]
Subject <- colorRampPalette(brewer.pal(8, "Set2"))(24)[all_CPS$indn]

RowSideColors<-cbind(Plate, Stimulation,Timepoint, Dose, Protection, Sex, Subject)

# heatmaps

heatmap3(m2,ColSideColors=RowSideColors, method="average", scale="row" , cexRow= 5, margins=c(2,6), 
         labRow="", labCol="", balanceColor=TRUE, showRowDendro=TRUE,showColDendro=TRUE)

heatmap3(m2,ColSideColors=RowSideColors, method="average", scale="column" , cexRow= 5, margins=c(2,6), 
         labRow="", labCol="", balanceColor=TRUE, showRowDendro=TRUE,showColDendro=TRUE)

heatmap3(m2,ColSideColors=RowSideColors, method="average", scale="none" , cexRow= 5, margins=c(2,6), 
         labRow="", labCol="", balanceColor=TRUE, showRowDendro=TRUE,showColDendro=TRUE)



#################################################### RTS,S #####################################################

rm(list=ls())
library(dplyr)
library(RColorBrewer)
library(heatmap3)
library(tidyverse)

## Read RMA normalized microarray data (8 plates)

df_rtss1<- read.delim("data/SysMalVac001_fulltxtRMA_plate1.TXT", check.names=FALSE)
df_rtss2<- read.delim("data/SysMalVac001_fulltxtRMA_plate2.TXT", check.names=FALSE)
df_rtss3<- read.delim("data/SysMalVac001_fulltxtRMA_plate3.TXT", check.names=FALSE)
df_rtss4<- read.delim("data/SysMalVac001_fulltxtRMA_plate4.TXT", check.names=FALSE)
df_rtss5<- read.delim("data/SysMalVac001_fulltxtRMA_plate5.TXT", check.names=FALSE)
df_rtss6<- read.delim("data/SysMalVac001_fulltxtRMA_plate6.TXT", check.names=FALSE)
df_rtss7<- read.delim("data/SysMalVac001_fulltxtRMA_plate7.TXT", check.names=FALSE)
df_rtss8<- read.delim("data/SysMalVac001_fulltxtRMA_plate8.TXT", check.names=FALSE)

colnames(df_rtss1)

coldel <- c( "Gene Description", "mRNA Accession" , "mRNA  Source",
             "mRna - Description", "UniGene ID",                   
             "Chromosome" ,  "GO Biological Process ID",     
             "GO Biological Process Term","GO Cellular Component ID",      
             "GO Cellular Component Term",  "GO Molecular Function ID",      
             "GO Molecular Function Term", "Pathway Source",                
             "Pathway Name")

df_rtss1 <- select(df_rtss1, -coldel)
df_rtss2 <- select(df_rtss2, -coldel)
df_rtss3 <- select(df_rtss3, -coldel)
df_rtss4 <- select(df_rtss4, -coldel)
df_rtss5 <- select(df_rtss5, -coldel)
df_rtss6 <- select(df_rtss6, -coldel)
df_rtss7 <- select(df_rtss7, -coldel)
df_rtss8 <- select(df_rtss8, -coldel)

97+95+95+94+95+96+97+25-8

# merge all plates
all_rtss <- Reduce(function(x,y) merge(x,y,by="Probe Set ID",all=TRUE) ,list(df_rtss1,df_rtss2,df_rtss3,df_rtss4,df_rtss5,df_rtss6,df_rtss7,df_rtss8))

dim(all_rtss)
colnames(all_rtss)

s <- colnames(all_rtss)
s1 <- sapply(strsplit(s, split='.', fixed=TRUE), function(x) (x[1]))
s1
s2 <- sapply(strsplit(s1, split='_', fixed=TRUE), function(x) (x[1]))
s2
s3 <- sapply(strsplit(s2, split='rep', fixed=TRUE), function(x) (x[1]))
s3

length(unique(s3))
length(s3)

names(all_rtss) <- s3
colnames(all_rtss)


## remove samples that not passed  QC
drops <- c("13SE3831", "13SE3906", "13SE3914", "13SE3926", "13SE3992", "13SE3993", 
           "13SE4022", "13SE4038", "13SE4046", "13SE4054", "13SE4077", "13SE4085", "13SE4093",
           "13SE4185", "13SE4188", "13SE4189", "13SE4196", "13SE4205", "13SE4211", "13SE4431",
           "13SE4432", "13SE4433", "13SE4434", "13SE4498", "13SE4499", "13SE4500", "13SE4502",
           "13SE4503", "13SE4504", "13SE4508", "13SE4510")

df_rtss <- all_rtss[ , !(names(all_rtss) %in% drops)]
dim(all_rtss)
dim(df_rtss)

colnames(df_rtss[448]) # drop sample repeated in plate 7 13SE4300
df_rtss <- df_rtss[ , -448]
colnames(df_rtss[448])
dim(df_rtss)    
str_detect(colnames(df_rtss),"13SE4300" )

# filter genes 
probes<- read.csv2("data/PGK_RTSS_ListFilterGenes.csv", check.names=FALSE)
probes$`Probe set ID` <- as.factor(probes$`Probe set ID`)
names(df_rtss)[1]<-"ID"
dim(df_rtss)

s1 <- probes$`Probe set ID`
df_rtss <- df_rtss %>%
  filter(ID %in% s1)
dim(df_rtss)


t_df_rtss <- t(df_rtss)
t_df_rtss <- as.data.frame(t_df_rtss)
t_df_rtss$sample <- row.names(t_df_rtss)

# drop row of probe set ID
t_df_rtss <-t_df_rtss[-1,]
head(t_df_rtss[1:3, 1:5])
dim(t_df_rtss)

# change sample name with .
str_detect(t_df_rtss$sample,"13SE4300.1" )
t_df_rtss[548,10568] <- "13SE4300"
str_detect(t_df_rtss$sample,"13SE4300.1")

# read clinical data

demog <- read.csv("data/RTSS_studydata.csv", check.names=FALSE)
dim(demog)
summary(demog)

demog <- demog[,-c(2,4:7,15,16)]
names(demog)[2]<-"sample"
demog$vaccine[demog$group_nb==1] <- "RTS,S"
demog$vaccine[demog$group_nb==2] <- "RTS,S"
demog$vaccine[demog$group_nb==3] <- "Comparator"

# merge 
m <- merge(t_df_rtss, demog ,by="sample", all.x=TRUE)
dim(m)
m.test <- merge(t_df_rtss, demog ,by="sample", all.x=FALSE)
dim(m.test)

head(t_df_rtss[,1:15])
head(m[,1:15])

m[,2:10568] <- sapply(m[,2:10568], as.character)
m[,2:10568] <- sapply(m[,2:10568], as.numeric)
head(m[,1:15])


# Add info of plate
p1 <- colnames(df_rtss1)
p2 <- colnames(df_rtss2)
p3 <- colnames(df_rtss3)
p4 <- colnames(df_rtss4)
p5 <- colnames(df_rtss5)
p6 <- colnames(df_rtss6)
p7 <- colnames(df_rtss7)


p1 <- sapply(strsplit(p1, split='.', fixed=TRUE), function(x) (x[1]))
p2 <- sapply(strsplit(p2, split='.', fixed=TRUE), function(x) (x[1]))
p2 <-sapply(strsplit(p2, split='rep', fixed=TRUE), function(x) (x[1]))
p3 <- sapply(strsplit(p3, split='.', fixed=TRUE), function(x) (x[1]))
p3 <-sapply(strsplit(p3, split='_', fixed=TRUE), function(x) (x[1]))
p4 <- sapply(strsplit(p4, split='.', fixed=TRUE), function(x) (x[1]))
p5 <- sapply(strsplit(p5, split='.', fixed=TRUE), function(x) (x[1]))
p6 <- sapply(strsplit(p6, split='.', fixed=TRUE), function(x) (x[1]))
p7 <- sapply(strsplit(p7, split='.', fixed=TRUE), function(x) (x[1]))

m$sample <- as.character(m$sample)

m <- m %>%
  mutate(plate = factor(
    ifelse(str_detect(sample, paste0(unlist(p1),collapse="|")),
           '1', ifelse(str_detect(sample, paste0(unlist(p2),collapse="|")),
                       '2', ifelse(str_detect(sample, paste0(unlist(p3),collapse="|")),
                                   '3', ifelse(str_detect(sample, paste0(unlist(p4),collapse="|")),
                                               '4', ifelse(str_detect(sample, paste0(unlist(p5),collapse="|")),
                                                           '5', ifelse(str_detect(sample, paste0(unlist(p6),collapse="|")),         
                                                                       '6', ifelse(str_detect(sample, paste0(unlist(p7),collapse="|")),
                                                                                   '7',  '8')))))))))


summary(m$plate)


head(m[,10568:10579])


# generate dataframe of FC CSP/DMSO and PfRBC/uRBC

## m to wide by stim m[,2:10568]
colnames(m[,10569:10579]) # "Stimulation" "site" "pid" "agec""visit" "casecon"  "group_nb""treatmnt" "vaccine" "sex" "plate"
summary(m$Stimulation) 
#CSP DMSO iRBC uRBC 
#184  172  151  147 

m_CSPvsDMSO <- m[m$Stimulation=="CSP" | m$Stimulation=="DMSO",]
m_PfRBCvsuRBC <- m[m$Stimulation=="iRBC" | m$Stimulation=="uRBC",]
356+298
dim(m)

m_CSPvsDMSO_w <-  reshape(m_CSPvsDMSO, timevar = "Stimulation", idvar = c("pid","site", "agec", "visit", "casecon", "group_nb","treatmnt", "vaccine", "sex", "plate"), direction="wide")
m_PfRBCvsuRBC_w <-  reshape(m_PfRBCvsuRBC, timevar = "Stimulation", idvar = c("pid","site", "agec", "visit", "casecon", "group_nb","treatmnt", "vaccine", "sex", "plate"), direction="wide")

test <- m_CSPvsDMSO_w[is.na(m_CSPvsDMSO_w$sample.CSP),]
dim(test) # 9
test <- m_CSPvsDMSO_w[is.na(m_CSPvsDMSO_w$sample.DMSO),]
dim(test) # 21
test <- m_PfRBCvsuRBC_w[is.na(m_PfRBCvsuRBC_w$sample.iRBC),]
dim(test) # 3
test <- m_PfRBCvsuRBC_w[is.na(m_PfRBCvsuRBC_w$sample.uRBC),]
dim(test) # 7

m_csp <- m_CSPvsDMSO_w[,1:10578]
m_dmso <- m_CSPvsDMSO_w[,c(1:10,10579:21146)]
head(m_csp[,10575:10578])
head(m_dmso[,10575:10578])
head(m_dmso[,1:12])
head(m_csp[,1:12])

# m_csp <- m_csp[,12:10578] gene expression data (log2)
# m_dmso<- m_dmso[,12:10578] gene expression data (log2)

ratio_CSP <- cbind(m_csp[,1:10], m_csp[,12:10578]-m_dmso[,12:10578]) # we obtain log2 FC iRBC/uRBC
head(ratio_CSP[, 1:13])


m_irbc <- m_PfRBCvsuRBC_w[,1:10578]
m_urbc <- m_PfRBCvsuRBC_w[,c(1:10,10579:21146)]
head(m_irbc[,10575:10578])
head(m_urbc[,1:12])
head(m_irbc[,1:12])

# m_irbc <- m_irbc[,12:10578] gene expression data (log2)
# m_urbc<- m_urbc[,12:10579] gene expression data (log2)

ratio_PfRBC<- cbind(m_irbc[,1:10], m_irbc[,12:10578]-m_urbc[,12:10578]) # we obtain log2 FC iRBC/uRBC

head(m_irbc[, 1:13])
head(m_urbc[, 1:13])
head(ratio_PfRBC[, 1:13])



# change names of variables to do rbind
s <- colnames(ratio_CSP)
s1 <- sapply(strsplit(s, split='.', fixed=TRUE), function(x) (x[1]))
s1
names(ratio_CSP) <- s1
colnames(ratio_CSP)

s <- colnames(ratio_PfRBC)
s1 <- sapply(strsplit(s, split='.', fixed=TRUE), function(x) (x[1]))
s1
names(ratio_PfRBC) <- s1
colnames(ratio_PfRBC)

# add stim variable
ratio_CSP$stim <- "CSP/DMSO"
ratio_PfRBC$stim <- "PfRBC/uRBC"
dim(ratio_CSP)
dim(ratio_PfRBC)

# rbind ratio dataframes
ratio_all <- rbind(ratio_CSP, ratio_PfRBC)
dim(ratio_all)
dim(ratio_CSP) + dim(ratio_PfRBC)

# Remove NA
ratio_all_final <- ratio_all[complete.cases(ratio_all[ ,]),]
dim(ratio_all_final) #347 is ratio_all (40 removed as expected)
head(ratio_all_final[,1:13])

m1 <- ratio_all_final[,11:10577]
head(m1[,1:3])


m2 <- t(m1)


#Heatmap

# color for variables 
ratio_all_final$vaccine <- as.factor(ratio_all_final$vaccine) 
ratio_all_final$stim <- as.factor(ratio_all_final$stim)
summary(ratio_all_final$plate)

Vaccine <- brewer.pal(3, "Set2")[1:2][ratio_all_final$vaccine]
Stimulation <- brewer.pal(3, "Set2")[1:2][ratio_all_final$stim]
Age <- brewer.pal(3, "Set2")[1:2][ratio_all_final$agec]
Plate <- brewer.pal(8, "Set2")[1:8][ratio_all_final$plate]
Timepoint <- brewer.pal(3, "Set2")[1:2][ratio_all_final$visit]
Sex <- brewer.pal(3, "Set2")[1:2][ratio_all_final$sex]
Site <-  brewer.pal(3, "Set2")[1:3][ratio_all_final$site]
Protection <- brewer.pal(3, "Set2")[1:2][ratio_all_final$casecon]

RowSideColors<-cbind(Plate, Stimulation, Vaccine, Timepoint, Age, Site, Sex, Protection)


heatmap3(m2,ColSideColors=RowSideColors, balanceColor=TRUE, method="average", scale="row" , cexRow= 5, margins=c(2,6), 
         labRow="", labCol="", showRowDendro=TRUE,showColDendro=TRUE)

heatmap3(m2,ColSideColors=RowSideColors, method="average", scale="none" , cexRow= 5, margins=c(2,6), 
         labRow="", labCol="", balanceColor=TRUE, showRowDendro=TRUE,showColDendro=TRUE)

heatmap3(m2,ColSideColors=RowSideColors, method="average", scale="none" , cexRow= 5, margins=c(2,6), 
         labRow="", labCol="", balanceColor=TRUE, showRowDendro=FALSE,showColDendro=TRUE)

