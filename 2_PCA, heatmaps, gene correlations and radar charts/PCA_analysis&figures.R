#################################### PCA #################################################################

# this is the code to perform PCA analysis and figures shown in the manuscript

rm(list=ls())

# packages
library(factoextra)
library(FactoMineR)
library(corrplot)
library(dplyr)
library(gridExtra)
library(grid)
library(ggpubr)
library(gdata)
library(reshape2)
library(ggplot2)
library(tidyverse)


########################################
################# CPS ##################
########################################

## Read RMA normalized microarray data (2 plates)

df_cps1 <- read.delim("data/SysMalVac_CPS_plate1_RMAnormalized.TXT", check.names=FALSE)
df_cps2 <- read.delim("data/SysMalVac_CPS_plate2_RMAnormalized.TXT", check.names=FALSE)

# merge plates

colnames(df_cps1)
colnames(df_cps2)

coldel <- c( "Gene Description", "mRNA Accession" , "mRNA  Source",
             "mRna - Description", "UniGene ID",                   
             "Chromosome" ,  "GO Biological Process ID",     
             "GO Biological Process Term","GO Cellular Component ID",      
             "GO Cellular Component Term",  "GO Molecular Function ID",      
             "GO Molecular Function Term", "Pathway Source",                
             "Pathway Name")

df_cps1 <- select(df_cps1, -coldel)
df_cps2 <- select(df_cps2, -coldel)
df_cps <- merge(df_cps1, df_cps2, by="Probe Set ID",all=TRUE)

s <- colnames(df_cps)
s1 <- sapply(strsplit(s, split='.', fixed=TRUE), function(x) (x[1]))
s1

names(df_cps) <- s1
colnames(df_cps)

# remove samples that did not pass QC
dim(df_cps)
df_cps <- within(df_cps, rm("13SE1071", "13SE2260"))
dim(df_cps)

# Only genes after filter (PGK)
probes<- read.csv2("data/PGK_CPS_ListFilterGenes.csv", check.names=FALSE)
probes$`Probe set ID` <- as.factor(probes$`Probe set ID`)

names(df_cps)[1]<-"ID"

s1 <- probes$`Probe set ID`
df_cps_f <- df_cps %>%
filter(ID %in% s1)


t_df_cps <- t(df_cps_f)
t_df_cps <- as.data.frame(t_df_cps)
t_df_cps$sample <- row.names(t_df_cps)

# drop row of probe set ID
t_df_cps <-t_df_cps[-1,]
head(t_df_cps[1:3, 1:5])
dim(t_df_cps)


# Load study data
demog <- read.csv("data/CPS_studydata.csv", check.names=FALSE)
dim(demog)
summary(demog)

# merge 
m <- merge(t_df_cps, demog ,by="sample", all.x=TRUE)
m.test <- merge(t_df_cps, demog ,by="sample", all.x=FALSE)
dim(m.test)
dim(m)

head(t_df_cps[1:3, 1:3])

m[,2:10794] <- sapply(m[,2:10794], as.character)
m[,2:10794] <- sapply(m[,2:10794], as.numeric)

m1 <- m[,2:10794]

# PCA

pca <- PCA(m1, graph = FALSE, scale.unit= TRUE) 
print(pca)

# PCA results
## get eigenvalues to check amount of variation retained by each principal component
eig.val <- get_eigenvalue(pca) 
print(eig.val) #the variance.percent indicates the % of the variation that is explained by each eigenvalue

## get the results for the variables
var <- get_pca_var(pca) 

## plot the contribution of each variable

### as a correlation plot (useful with few variables)
g <- corrplot(var$contrib, is.corr=FALSE) 

### with lots of variables it is easier to visualize it with a barplot of the top 10 (for example) contributing variables
pdf("CPS_Contributions _barplot_crude.pdf")
fviz_contrib(pca, choice = "var", axes = 1, top = 10) # contributions of variables to PC1
fviz_contrib(pca, choice = "var", axes = 2, top = 10) # contributions of variables to PC2
fviz_contrib(pca, choice = "var", axes = 3, top = 10) # contributions of variables to PC3
fviz_contrib(pca, choice = "var", axes = 1:2, top = 10)
fviz_contrib(pca, choice = "var", axes = c(1,3), top = 10)
fviz_contrib(pca, choice = "var", axes = 2:3, top = 10)
dev.off()

#The red dashed line on the graph above indicates the expected average contribution. 
#If the contribution of the variables were uniform, the expected value would be 
#1/length(variables) = 1/10 = 10%. For a given component, a variable with a contribution 
#larger than this cutoff could be considered as important in contributing to the component.

contrib <- data.frame(g) #make it a dataframe

pdf("CPS_contributions_table.pdf", height=26, width=20) #create a table
grid.table(contrib[1:14, ])
dev.off()

# PCA graphs
## plot the percentage of each component that explains the variance

pdf("CPS_Scree plot.pdf")
fviz_eig(pca, addlabels = TRUE, ylim = c(0, 100)) 
dev.off()

## variable correlation plot

fviz_pca_var(pca, col.var = "black")


#get coordinates of dimensions
coord <- data.frame(pca$ind$coord)


# PCA plots components 1-3 by each study variable of interest

p1 <- fviz_pca_ind(pca, axes = c(1, 2),#axes indicate which components to represent
                   geom.ind = "point", # show points only (but not "text")
                   palette = "npg",
                   col.ind = m$stimulus,
                   addEllipses = TRUE, 
                   legend.title = "Stim")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$stimulus,
                   addEllipses = TRUE, 
                   legend.title = "Stim")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$stimulus,
                   addEllipses = TRUE, 
                   legend.title = "Stim")

# arrange the plots the way you want the to appear in the document
pdf("CPS_PCA_stim.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)#save as pdf
dev.off()


p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$timepoint,
                   addEllipses = TRUE, 
                   legend.title = "timepoint")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$timepoint,
                   addEllipses = TRUE, 
                   legend.title = "timepoint")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$timepoint,
                   addEllipses = TRUE, 
                   legend.title = "timepoint")


pdf("CPS_PCA_timepoint.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()


p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$protected,
                   addEllipses = TRUE, 
                   legend.title = "Protection")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$protected,
                   addEllipses = TRUE, 
                   legend.title = "Protection")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$protected,
                   addEllipses = TRUE, 
                   legend.title = "Protection")

pdf("CPS_PCA_protected.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()

table(m$group, m$timepoint)

p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   col.ind = m$donor,
                   addEllipses = TRUE, 
                   legend.title = "donor")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   col.ind = m$donor,
                   addEllipses = TRUE, 
                   legend.title = "donor")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   col.ind = m$donor,
                   addEllipses = TRUE, 
                   legend.title = "donor")

pdf("CPS_PCA_donor.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()


p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   col.ind = m$group,
                   addEllipses = TRUE, 
                   legend.title = "group")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   col.ind = m$group,
                   addEllipses = TRUE, 
                   legend.title = "group")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   col.ind = m$group,
                   addEllipses = TRUE, 
                   legend.title = "group")

pdf("CPS_PCA_group dose.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()

m$plate <- as.factor(m$plate)

p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$plate,
                   addEllipses = TRUE, 
                   legend.title = "Plate")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$plate,
                   addEllipses = TRUE, 
                   legend.title = "Plate")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$plate,
                   addEllipses = TRUE, 
                   legend.title = "Plate")

pdf("CPS_PCA_plate.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()


p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$Sex,
                   addEllipses = TRUE, 
                   legend.title = "Sex")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$Sex,
                   addEllipses = TRUE, 
                   legend.title = "Sex")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$Sex,
                   addEllipses = TRUE, 
                   legend.title = "Sex")

pdf("CPS_PCA_sex.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()



#############################################################
########## CPS Fold change CSP/DMSO and PfRBC/uRBC ##########
#############################################################

# Load FC CSP/DMSO and PfRBC/uRBC microarray data
all_CPS <- read.csv("data/CPS_FC_microarray_data.csv")

# read study data
demog <- read.csv("CPS_studydata.csv")
summary(demog)

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


# PCA

pca <- PCA(m1, graph = FALSE, scale.unit= TRUE) 
print(pca)

# PCA results
## get eigenvalues to check amount of variation retained by each principal component
eig.val <- get_eigenvalue(pca) 
print(eig.val) #the variance.percent indicates the % of the variation that is explained by each eigenvalue

## get the results for the variables
var <- get_pca_var(pca) 

## plot the contribution of each variable

### as a correlation plot (useful with few variables)
g <- corrplot(var$contrib, is.corr=FALSE) 

### with lots of variables it is easier to visualize it with a barplot of the top 10 (for example) contributing variables
pdf("RTSS_Contributions _barplot_crude.pdf")
fviz_contrib(pca, choice = "var", axes = 1, top = 10) # contributions of variables to PC1
fviz_contrib(pca, choice = "var", axes = 2, top = 10) # contributions of variables to PC2
fviz_contrib(pca, choice = "var", axes = 3, top = 10) # contributions of variables to PC3
fviz_contrib(pca, choice = "var", axes = 1:2, top = 10)
fviz_contrib(pca, choice = "var", axes = c(1,3), top = 10)
fviz_contrib(pca, choice = "var", axes = 2:3, top = 10)
dev.off()

#The red dashed line on the graph above indicates the expected average contribution. 
#If the contribution of the variables were uniform, the expected value would be 
#1/length(variables) = 1/10 = 10%. For a given component, a variable with a contribution 
#larger than this cutoff could be considered as important in contributing to the component.

contrib <- data.frame(g) #make it a dataframe

pdf("CPS_FC_contributions_table.pdf", height=26, width=20) #create a table
grid.table(contrib[1:14, ])
dev.off()

# PCA graphs
## plot the percentage of each component that explains the variance

pdf("CPS_FC_Scree plot.pdf")
fviz_eig(pca, addlabels = TRUE, ylim = c(0, 100)) 
dev.off()

## variable correlation plot

fviz_pca_var(pca, col.var = "black")

## plot by compartment

#get coordinates of dimensions
coord <- data.frame(pca$ind$coord)


p1 <- fviz_pca_ind(pca, axes = c(1, 2),#axes indicate which components to represent
                   geom.ind = "point", # show points only (but not "text")
                   palette = "npg",
                   col.ind = all_CPS$stim,
                   addEllipses = TRUE, 
                   legend.title = "Stimulation")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = all_CPS$stim,
                   addEllipses = TRUE, 
                   legend.title = "Stimulation")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = all_CPS$stim,
                   addEllipses = TRUE, 
                   legend.title = "Stimulation")

# arrange the plots the way you want the to appear in the document
pdf("CPS_FC_PCA_stimulation.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)#save as pdf
dev.off()



p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = all_CPS$protected,
                   addEllipses = TRUE, 
                   legend.title = "Protection")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = all_CPS$protected,
                   addEllipses = TRUE, 
                   legend.title = "Protection")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = all_CPS$protected,
                   addEllipses = TRUE, 
                   legend.title = "Protection")

pdf("CPS_FC_PCA_protection.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()

p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = all_CPS$timepoint,
                   addEllipses = TRUE, 
                   legend.title = "Timepoint")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = all_CPS$timepoint,
                   addEllipses = TRUE, 
                   legend.title = "Timepoint")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = all_CPS$timepoint,
                   addEllipses = TRUE, 
                   legend.title = "Timepoint")

pdf("CPS_FC_PCA_timepoint.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()

p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = all_CPS$indn,
                   addEllipses = TRUE, 
                   legend.title = "Subject")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = all_CPS$indn,
                   addEllipses = TRUE, 
                   legend.title = "Subject")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = all_CPS$indn,
                   addEllipses = TRUE, 
                   legend.title = "Subject")

pdf("CPS_FC_PCA_subject.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()




p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point",
                   palette = "npg",
                   col.ind = as.factor(all_CPS$plate),
                   addEllipses = TRUE, 
                   legend.title = "Plate")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = as.factor(all_CPS$plate),
                   addEllipses = TRUE, 
                   legend.title = "Plate")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = as.factor(all_CPS$plate),
                   addEllipses = TRUE, 
                   legend.title = "Plate")


pdf("CPS_FC_PCA_plate_ok.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()



p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = as.factor(all_CPS$group),
                   addEllipses = TRUE, 
                   legend.title = "Group")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = as.factor(all_CPS$group),
                   addEllipses = TRUE, 
                   legend.title = "Group")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = as.factor(all_CPS$group),
                   addEllipses = TRUE, 
                   legend.title = "Group")


pdf("CPS_FC_PCA_group dose.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()


p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = all_CPS$Sex,
                   addEllipses = TRUE, 
                   legend.title = "Sex")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = all_CPS$Sex,
                   addEllipses = TRUE, 
                   legend.title = "Sex")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = all_CPS$Sex,
                   addEllipses = TRUE, 
                   legend.title = "Sex")


pdf("CPS_FC_PCA_sex.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()





###########################################
################# RTS,S ###################
###########################################


rm(list=ls())

# packages
library(devtools)
library(factoextra)
library(FactoMineR)
library(corrplot)
library(dplyr)
library(gridExtra)
library(grid)
library(ggpubr)
library(gdata)
library(reshape2)
library(ggplot2)
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

m1 <- m[,2:10568]



# PCA

pca <- PCA(m1, graph = FALSE, scale.unit= TRUE) 
print(pca)

# PCA results
## get eigenvalues to check amount of variation retained by each principal component
eig.val <- get_eigenvalue(pca) 
print(eig.val) #the variance.percent indicates the % of the variation that is explained by each eigenvalue

## get the results for the variables
var <- get_pca_var(pca) 

## plot the contribution of each variable

### as a correlation plot (useful with few variables)
g <- corrplot(var$contrib, is.corr=FALSE) 

### with lots of variables it is easier to visualize it with a barplot of the top 10 (for example) contributing variables
pdf("RTSS_Contributions _barplot_crude.pdf")
fviz_contrib(pca, choice = "var", axes = 1, top = 10) # contributions of variables to PC1
fviz_contrib(pca, choice = "var", axes = 2, top = 10) # contributions of variables to PC2
fviz_contrib(pca, choice = "var", axes = 3, top = 10) # contributions of variables to PC3
fviz_contrib(pca, choice = "var", axes = 1:2, top = 10)
fviz_contrib(pca, choice = "var", axes = c(1,3), top = 10)
fviz_contrib(pca, choice = "var", axes = 2:3, top = 10)
dev.off()

#The red dashed line on the graph above indicates the expected average contribution. 
#If the contribution of the variables were uniform, the expected value would be 
#1/length(variables) = 1/10 = 10%. For a given component, a variable with a contribution 
#larger than this cutoff could be considered as important in contributing to the component.

contrib <- data.frame(g) #make it a dataframe

pdf("RTSS_contributions_table.pdf", height=26, width=20) #create a table
grid.table(contrib[1:14, ])
dev.off()

# PCA graphs
## plot the percentage of each component that explains the variance

pdf("RTSS_Scree plot.pdf")
fviz_eig(pca, addlabels = TRUE, ylim = c(0, 100)) 
dev.off()

## variable correlation plot

fviz_pca_var(pca, col.var = "black")

## plot by compartment

#get coordinates of dimensions
coord <- data.frame(pca$ind$coord)



p1 <- fviz_pca_ind(pca, axes = c(1, 2),#axes indicate which components to represent
                   geom.ind = "point", # show points only (but not "text")
                   palette = "npg",
                   col.ind = m$Stimulation,
                   addEllipses = TRUE, 
                   legend.title = "Stim")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$Stimulation,
                   addEllipses = TRUE, 
                   legend.title = "Stim")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$Stimulation,
                   addEllipses = TRUE, 
                   legend.title = "Stim")

# arrange the plots the way you want the to appear in the document
pdf("RTSS_PCA_stim.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)#save as pdf
dev.off()


p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$agec,
                   addEllipses = TRUE, 
                   legend.title = "Age group")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$agec,
                   addEllipses = TRUE, 
                   legend.title = "Age group")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$agec,
                   addEllipses = TRUE, 
                   legend.title = "Age group")

pdf("RTSS_PCA_age.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()



p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$casecon,
                   addEllipses = TRUE, 
                   legend.title = "Protection")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$casecon,
                   addEllipses = TRUE, 
                   legend.title = "Protection")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$casecon,
                   addEllipses = TRUE, 
                   legend.title = "Protection")


pdf("RTSS_PCA_protection.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()



p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$site,
                   addEllipses = TRUE, 
                   legend.title = "Site")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$site,
                   addEllipses = TRUE, 
                   legend.title = "Site")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$site,
                   addEllipses = TRUE, 
                   legend.title = "Site")

pdf("RTSS_PCA_site.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()



p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = as.factor(m$vaccine),
                   addEllipses = TRUE, 
                   legend.title = "Vaccine")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = as.factor(m$vaccine),
                   addEllipses = TRUE, 
                   legend.title = "Vaccine")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = as.factor(m$vaccine),
                   addEllipses = TRUE, 
                   legend.title = "Vaccine")

pdf("RTSS_PCA_vaccine.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()



p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = as.factor(m$visit),
                   addEllipses = TRUE, 
                   legend.title = "Timepoint")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = as.factor(m$visit),
                   addEllipses = TRUE, 
                   legend.title = "Timepoint")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = as.factor(m$visit),
                   addEllipses = TRUE, 
                   legend.title = "Timepoint")

pdf("RTSS_PCA_timepoint.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()



p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$plate,
                   addEllipses = TRUE, 
                   legend.title = "Plate")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$plate,
                   addEllipses = TRUE, 
                   legend.title = "Plate")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$plate,
                   addEllipses = TRUE, 
                   legend.title = "Plate")

pdf("RTSS_PCA_plate.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()     
            

p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$sex,
                   addEllipses = TRUE, 
                   legend.title = "Sex")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$sex,
                   addEllipses = TRUE, 
                   legend.title = "Sex")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = m$sex,
                   addEllipses = TRUE, 
                   legend.title = "Sex")

pdf("RTSS_PCA_sex.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()     



    

###################################   
########## RTSS FC stims ########## 
###################################  

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
  
# PCA

pca <- PCA(m1, graph = FALSE, scale.unit= TRUE) 
print(pca)

# PCA results
## get eigenvalues to check amount of variation retained by each principal component
eig.val <- get_eigenvalue(pca) 
print(eig.val) #the variance.percent indicates the % of the variation that is explained by each eigenvalue

## get the results for the variables
var <- get_pca_var(pca) 

## plot the contribution of each variable

### as a correlation plot (useful with few variables)
g <- corrplot(var$contrib, is.corr=FALSE) 

### with lots of variables it is easier to visualize it with a barplot of the top 10 (for example) contributing variables
pdf("RTSS_FC_Contributions _barplot_crude.pdf")
fviz_contrib(pca, choice = "var", axes = 1, top = 10) # contributions of variables to PC1
fviz_contrib(pca, choice = "var", axes = 2, top = 10) # contributions of variables to PC2
fviz_contrib(pca, choice = "var", axes = 3, top = 10) # contributions of variables to PC3
fviz_contrib(pca, choice = "var", axes = 1:2, top = 10)
fviz_contrib(pca, choice = "var", axes = c(1,3), top = 10)
fviz_contrib(pca, choice = "var", axes = 2:3, top = 10)
dev.off()

#The red dashed line on the graph above indicates the expected average contribution. 
#If the contribution of the variables were uniform, the expected value would be 
#1/length(variables) = 1/10 = 10%. For a given component, a variable with a contribution 
#larger than this cutoff could be considered as important in contributing to the component.

contrib <- data.frame(g) #make it a dataframe

pdf("RTSS_FC_contributions_table.pdf", height=26, width=20) #create a table
grid.table(contrib[1:14, ])
dev.off()

# PCA graphs
## plot the percentage of each component that explains the variance

pdf("RTSS_FC_Scree plot.pdf")
fviz_eig(pca, addlabels = TRUE, ylim = c(0, 100)) 
dev.off()

## variable correlation plot

fviz_pca_var(pca, col.var = "black")

## plot by compartment

#get coordinates of dimensions
coord <- data.frame(pca$ind$coord)



p1 <- fviz_pca_ind(pca, axes = c(1, 2),#axes indicate which components to represent
                   geom.ind = "point", # show points only (but not "text")
                   palette = "npg",
                   col.ind = ratio_all_final$stim,
                   addEllipses = TRUE, 
                   legend.title = "Stimulation")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = ratio_all_final$stim,
                   addEllipses = TRUE, 
                   legend.title = "Stimulation")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = ratio_all_final$stim,
                   addEllipses = TRUE, 
                   legend.title = "Stimulation")

# arrange the plots the way you want the to appear in the document
pdf("RTSS_FC_PCA_stimulation.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)#save as pdf
dev.off()

p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = ratio_all_final$visit,
                   addEllipses = TRUE, 
                   legend.title = "Timepoint")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = ratio_all_final$visit,
                   addEllipses = TRUE, 
                   legend.title = "Timepoint")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = ratio_all_final$visit,
                   addEllipses = TRUE, 
                   legend.title = "Timepoint")

pdf("RTSS_FC_PCA_timepoint.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()


p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = ratio_all_final$casecon,
                   addEllipses = TRUE, 
                   legend.title = "Protection")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = ratio_all_final$casecon,
                   addEllipses = TRUE, 
                   legend.title = "Protection")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = ratio_all_final$casecon,
                   addEllipses = TRUE, 
                   legend.title = "Protection")

pdf("RTSS_FC_PCA_protection.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()

p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = ratio_all_final$vac,
                   addEllipses = TRUE, 
                   legend.title = "Vaccine")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = ratio_all_final$vac,
                   addEllipses = TRUE, 
                   legend.title = "Vaccine")


p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = ratio_all_final$vac,
                   addEllipses = TRUE, 
                   legend.title = "Vaccine")


pdf("RTSS_FC_PCA_vaccine.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()




summary(ratio_all_final$plate)

p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = ratio_all_final$plate,
                   addEllipses = TRUE, 
                   legend.title = "Plate")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = ratio_all_final$plate,
                   addEllipses = TRUE, 
                   legend.title = "Plate")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = ratio_all_final$plate,
                   addEllipses = TRUE, 
                   legend.title = "Plate")

pdf("RTSS_FC_PCA_plate.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()


p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = ratio_all_final$agec,
                   addEllipses = TRUE, 
                   legend.title = "Age cohort")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = ratio_all_final$agec,
                   addEllipses = TRUE, 
                   legend.title = "Age cohort")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = ratio_all_final$agec,
                   addEllipses = TRUE, 
                   legend.title = "Age cohort")

pdf("RTSS_FC_PCA_Agecohort.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)#save as pdf
dev.off()



p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = as.factor(ratio_all_final$vaccine),
                   addEllipses = TRUE, 
                   legend.title = "Vaccine group")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = as.factor(ratio_all_final$vaccine),
                   addEllipses = TRUE, 
                   legend.title = "Vaccine group")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = as.factor(ratio_all_final$vaccine),
                   addEllipses = TRUE, 
                   legend.title = "Vaccine group")

pdf("RTSS_FC_PCA_vaccinegroup.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()

p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = as.factor(ratio_all_final$site),
                   addEllipses = TRUE, 
                   legend.title = "Site")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = as.factor(ratio_all_final$site),
                   addEllipses = TRUE, 
                   legend.title = "Site")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point", 
                   palette = "npg",
                   col.ind = as.factor(ratio_all_final$site),
                   addEllipses = TRUE, 
                   legend.title = "Site")

pdf("RTSS_FC_PCA_site.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()


p1 <- fviz_pca_ind(pca, axes = c(1, 2),
                   geom.ind = "point",
                   palette = "npg",
                   col.ind = ratio_all_final$sex,
                   addEllipses = TRUE, 
                   legend.title = "Sex")

p2 <- fviz_pca_ind(pca, axes = c(1, 3),
                   geom.ind = "point",
                   palette = "npg",
                   col.ind =ratio_all_final$sex,
                   addEllipses = TRUE, 
                   legend.title = "Sex")

p3 <- fviz_pca_ind(pca, axes = c(2, 3),
                   geom.ind = "point",
                   palette = "npg",
                   col.ind = ratio_all_final$sex,
                   addEllipses = TRUE, 
                   legend.title = "Sex")

pdf("RTSS_FC_PCA_sex.pdf", width = 15, height = 4)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()
