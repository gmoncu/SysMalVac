############################################################################################
################## Analysis of cell phenotyping by flow cytometry ########################
############################################################################################


# This is the code to perform descriptive plots and analysis for the manuscript 

rm(list=ls())
library(gdata)
library(Hmisc)
library(plyr)
library(ggplot2)
library(reshape2)
library(compareGroups)
library(RColorBrewer)
library(ggsci)


# load data
pheno.all <- read.csv2("data/pheno_all.csv", dec=",")


# exclude by low viability or low total cell counts
pheno.all$exclusion <- with(pheno.all, ifelse(viable<50 | total<60,1,0))
with(pheno.all, table(exclusion, useNA="always"))
summary(pheno.all$exclusion)


############################ STATISTICAL ANALYSIS #########################################

########################
######### CPS ##########
########################

# Comparison of CPS timepoints (i7 is pre-immunization, c1 is post-immunization)

## prepare dataframe
pheno.cps <- pheno.all[pheno.all$cohort=="CPS",]
pheno.cps <- pheno.cps[,c(1:2,5:14, 17,20)]
pheno.cps <- pheno.cps[pheno.cps$dose!="3x8",]
pheno.cps.wide <- reshape(pheno.cps, idvar = "pid", timevar = "time.point", direction = "wide")

## wilcoxon test
c <- sapply(2:11, function(x) wilcox.test(pheno.cps.wide[,x],pheno.cps.wide[,x+12],paired=TRUE)[[3]])
r <- data.frame(cell.sub.i7 = colnames(pheno.cps.wide[,2:11]), cell.sub.c1 = colnames(pheno.cps.wide[,14:23]), pval = round(unlist(c),3), padj=round(p.adjust(unlist(c), method = "holm"),3), row.names = NULL)
head(r)
write.csv(r, "CPS_phenotyping_timepoint.csv", row.names=FALSE)



# Comparison of CPS protected and non-protected ("casecon" variable)

summary(pheno.all$dose)
pheno.all <- subset(pheno.all, dose!="3x8") # remove validation individuals
summary(pheno.all$dose)

#post-immunization
c <- compareGroups(casecon ~ cd14p + cd14m + cd19b + cd4t +cd8t + gdt + nkt + cd56d + cd56h + hladr, data = pheno.all, 
                   subset = cohort=="CPS" & time.point=="c1" , method = 2) # method 2 is Wilcoxon rank sum test
t <- createTable(c)
p <- print(t)
pvals <- getResults(t, "p.overall")
padj <- p.adjust(pvals, method = "holm")
padj <- as.data.frame(padj)
padj <- rbind(" ",padj)
export2csv(t, file="CPS_phenotyping_casecon_c1.csv", sep = ";")
tmp <- read.csv2("CPS_phenotyping_casecon_c1.csv", dec=".")
tmp <- cbind(tmp, padj)
write.csv2(tmp, "CPS_phenotyping_casecon_c1.csv")

# pre-immunization
c <- compareGroups(casecon ~ cd14p + cd14m + cd19b + cd4t +cd8t + gdt + nkt + cd56d + cd56h + hladr, data = pheno.all, 
                   subset = cohort=="CPS" & time.point=="i7" , method = 2)
t <- createTable(c)
p <- print(t)
pvals <- getResults(t, "p.overall")
padj <- p.adjust(pvals, method = "holm")
padj <- as.data.frame(padj)
padj <- rbind(" ",padj)
export2csv(t, file="CPS_phenotyping_casecon_i7.csv", sep = ";")
tmp <- read.csv2("CPS_phenotyping_casecon_i7.csv", dec=".")
tmp <- cbind(tmp, padj)
write.csv2(tmp, "CPS_phenotyping_casecon_i7.csv")



#########################
######### RTS,S #########
#########################

# load data
pheno.all <- read.csv2("pheno_all.csv", dec=",")


# exclude by low viability or low total cell counts
pheno.all$exclusion <- with(pheno.all, ifelse(viable<50 | total<60,1,0))
with(pheno.all, table(exclusion, useNA="always"))
summary(pheno.all$exclusion)

# Comparison of protected and non-protected RTS,S-vaccinees one month post-vaccination (M3)

c <- compareGroups(casecon ~ cd14p + cd14m + cd19b + cd4t +cd8t + gdt + nkt + cd56d + cd56h + hladr, data = pheno.all, 
              subset = vaccine=="rtss" & time.point=="m3", method = 2)
t <- createTable(c)
p <- print(t)
pvals <- getResults(t, "p.overall")
padj <- p.adjust(pvals, method = "holm")
padj <- as.data.frame(padj)
padj <- rbind(" ",padj)
export2csv(t, file="RTSS_phenotyping_casecon_m3.csv", sep = ";")
tmp <- read.csv2("RTSS_phenotyping_casecon_m3.csv", dec=".")
tmp <- cbind(tmp, padj)
write.csv2(tmp, "RTSS_phenotyping_casecon_m3.csv")


# Comparison of protected and non-protected RTS,S-vaccinees at baseline (M0)
c <- compareGroups(casecon ~ cd14p + cd14m + cd19b + cd4t +cd8t + gdt + nkt + cd56d + cd56h + hladr, data = pheno.all, 
                   subset = vaccine=="rtss" & time.point=="m0", method = 2)
t <- createTable(c)
p <- print(t)
pvals <- getResults(t, "p.overall")
padj <- p.adjust(pvals, method = "holm")
padj <- as.data.frame(padj)
padj <- rbind(" ",padj)
export2csv(t, file="RTSS_phenotyping_casecon_m0.csv", sep = ";")
tmp <- read.csv2("RTSS_phenotyping_casecon_m0.csv", dec=".")
tmp <- cbind(tmp, padj)
write.csv2(tmp, "RTSS_phenotyping_casecon_m0.csv")

# Comparison of protected and non-protected comparator-vaccinees at baseline (M0)
c <- compareGroups(casecon ~ cd14p + cd14m + cd19b + cd4t +cd8t + gdt + nkt + cd56d + cd56h + hladr, data = pheno.all, 
                   subset = vaccine=="comparator" & time.point=="m0", method = 2)
t <- createTable(c)
p <- print(t)
pvals <- getResults(t, "p.overall")
padj <- p.adjust(pvals, method = "holm")
padj <- as.data.frame(padj)
padj <- rbind(" ",padj)
export2csv(t, file="RTSS_phenotyping_casecon_m0_comparator.csv", sep = ";")
tmp <- read.csv2("RTSS_phenotyping_casecon_m0_comparator.csv", dec=".")
tmp <- cbind(tmp, padj)
write.csv2(tmp, "RTSS_phenotyping_casecon_m0_comparator.csv")

# Comparison of protected and non-protected comparator-vaccinees one month post-vaccination (M3)
c <- compareGroups(casecon ~ cd14p + cd14m + cd19b + cd4t +cd8t + gdt + nkt + cd56d + cd56h + hladr, data = pheno.all, 
                   subset = vaccine=="comparator" & time.point=="m3", method = 2)
t <- createTable(c)
p <- print(t)
pvals <- getResults(t, "p.overall")
padj <- p.adjust(pvals, method = "holm")
padj <- as.data.frame(padj)
padj <- rbind(" ",padj)
export2csv(t, file="RTSS_phenotyping_casecon_m3_comparator.csv", sep = ";")
tmp <- read.csv2("RTSS_phenotyping_casecon_m3_comparator.csv", dec=".")
tmp <- cbind(tmp, padj)
write.csv2(tmp, "RTSS_phenotyping_casecon_m3_comparator.csv")


# Comparison of RTS,S and comparators one month post-vaccination (M3)
c <- compareGroups(vaccine ~ cd14p + cd14m + cd19b + cd4t +cd8t + gdt + nkt + cd56d + cd56h + hladr, data = pheno.all, 
                   subset = cohort=="RTSS" & time.point=="m3", method = 2)
t <- createTable(c)
p <- print(t)
pvals <- getResults(t, "p.overall")
padj <- p.adjust(pvals, method = "holm")
padj <- as.data.frame(padj)
padj <- rbind(" ",padj)
export2csv(t, file= "RTSS_phenotyping_vaccine_m3.csv", sep = ";")
tmp <- read.csv2("RTSS_phenotyping_vaccine_m3.csv", dec=".")
tmp <- cbind(tmp, padj)
write.csv2(tmp, "RTSS_phenotyping_vaccine_m3.csv")


# Comparison of RTS,S and comparators one month post-vaccination (M0)
c <- compareGroups(vaccine ~ cd14p + cd14m + cd19b + cd4t +cd8t + gdt + nkt + cd56d + cd56h + hladr, data = pheno.all, 
                   subset = cohort=="RTSS" & time.point=="m0", method = 2)
t <- createTable(c)
p <- print(t)
pvals <- getResults(t, "p.overall")
padj <- p.adjust(pvals, method = "holm")
padj <- as.data.frame(padj)
padj <- rbind(" ",padj)
export2csv(t, file= "RTSS_phenotyping_vaccine_m0.csv", sep = ";")
tmp <- read.csv2("RTSS_phenotyping_vaccine_m0.csv", dec=".")
tmp <- cbind(tmp, padj)
write.csv2(tmp, "RTSS_phenotyping_vaccine_m0.csv")


### Immunogenicity and Case-Control analysis stratified by AGE COHORT (Infants 6-12 weeks old and Children 5-17 month old)

c <- compareGroups(vaccine ~ cd14p + cd14m + cd19b + cd4t +cd8t + gdt + nkt + cd56d + cd56h + hladr, data = pheno.all, 
                   subset = cohort=="RTSS" & time.point=="m3" & agec=="6-12w", method = 2)
t <- createTable(c)
p <- print(t)
pvals <- getResults(t, "p.overall")
padj <- p.adjust(pvals, method = "holm")
padj <- as.data.frame(padj)
padj <- rbind(" ",padj)
export2csv(t, file= "RTSS_phenotyping_vaccine_m3_Infants.csv", sep = ";")
tmp <- read.csv2("RTSS_phenotyping_vaccine_m3_Infants.csv", dec=".")
tmp <- cbind(tmp, padj)
write.csv2(tmp, "RTSS_phenotyping_vaccine_m3_Infants.csv")


c <- compareGroups(vaccine ~ cd14p + cd14m + cd19b + cd4t +cd8t + gdt + nkt + cd56d + cd56h + hladr, data = pheno.all, 
                   subset = cohort=="RTSS" & time.point=="m3" & agec=="5-17m", method = 2)
t <- createTable(c)
p <- print(t)
pvals <- getResults(t, "p.overall")
padj <- p.adjust(pvals, method = "holm")
padj <- as.data.frame(padj)
padj <- rbind(" ",padj)
export2csv(t, file= "RTSS_phenotyping_vaccine_m3_Children.csv", sep = ";")
tmp <- read.csv2("RTSS_phenotyping_vaccine_m3_Children.csv", dec=".")
tmp <- cbind(tmp, padj)
write.csv2(tmp, "RTSS_phenotyping_vaccine_m3_Children.csv")



c <- compareGroups(casecon ~ cd14p + cd14m + cd19b + cd4t +cd8t + gdt + nkt + cd56d + cd56h + hladr, data = pheno.all, 
                   subset = vaccine=="rtss" & time.point=="m3" & agec=="6-12w", method = 2)
t <- createTable(c)
p <- print(t)
pvals <- getResults(t, "p.overall")
padj <- p.adjust(pvals, method = "holm")
padj <- as.data.frame(padj)
padj <- rbind(" ",padj)
export2csv(t, file="RTSS_phenotyping_casecon_m3_RTSS_Infants.csv", sep = ";")
tmp <- read.csv2("RTSS_phenotyping_casecon_m3_RTSS_Infants.csv", dec=".")
tmp <- cbind(tmp, padj)
write.csv2(tmp, "RTSS_phenotyping_casecon_m3_RTSS_Infants.csv")

c <- compareGroups(casecon ~ cd14p + cd14m + cd19b + cd4t +cd8t + gdt + nkt + cd56d + cd56h + hladr, data = pheno.all, 
                   subset = vaccine=="rtss" & time.point=="m3"& agec=="5-17m", method = 2)
t <- createTable(c)
p <- print(t)
pvals <- getResults(t, "p.overall")
padj <- p.adjust(pvals, method = "holm")
padj <- as.data.frame(padj)
padj <- rbind(" ",padj)
export2csv(t, file="RTSS_phenotyping_casecon_m3_RTSS_Children.csv", sep = ";")
tmp <- read.csv2("RTSS_phenotyping_casecon_m3_RTSS_Children.csv", dec=".")
tmp <- cbind(tmp, padj)
write.csv2(tmp, "RTSS_phenotyping_casecon_m3_RTSS_Children.csv")


c <- compareGroups(casecon ~ cd14p + cd14m + cd19b + cd4t +cd8t + gdt + nkt + cd56d + cd56h + hladr, data = pheno.all, 
                   subset = vaccine=="comparator" & time.point=="m3" & agec=="6-12w", method = 2)
t <- createTable(c)
p <- print(t)
pvals <- getResults(t, "p.overall")
padj <- p.adjust(pvals, method = "holm")
padj <- as.data.frame(padj)
padj <- rbind(" ",padj)
export2csv(t, file="RTSS_phenotyping_casecon_m3_comparator_Infants.csv", sep = ";")
tmp <- read.csv2("RTSS_phenotyping_casecon_m3_comparator_Infants.csv", dec=".")
tmp <- cbind(tmp, padj)
write.csv2(tmp, "RTSS_phenotyping_casecon_m3_comparator_Infants.csv")

c <- compareGroups(casecon ~ cd14p + cd14m + cd19b + cd4t +cd8t + gdt + nkt + cd56d + cd56h + hladr, data = pheno.all, 
                   subset = vaccine=="comparator" & time.point=="m3" & agec=="5-17m", method = 2)
t <- createTable(c)
p <- print(t)
pvals <- getResults(t, "p.overall")
padj <- p.adjust(pvals, method = "holm")
padj <- as.data.frame(padj)
padj <- rbind(" ",padj)
export2csv(t, file="RTSS_phenotyping_casecon_m3_comparator_Children.csv", sep = ";")
tmp <- read.csv2("RTSS_phenotyping_casecon_m3_comparator_Children.csv", dec=".")
tmp <- cbind(tmp, padj)
write.csv2(tmp, "RTSS_phenotyping_casecon_m3_comparator_Children.csv")



############################ BOX PLOTS #########################################

# reshape dataframe to long
pheno.all <- pheno.all[,c(1:3,5:14,16:20)]
pheno.all.long <- melt(pheno.all, id.vars = c(1:3,14:18))
colnames(pheno.all.long)[9] <- "cell.subset"
colnames(pheno.all.long)[10] <- "frequency"
pheno.all.long$cohort <- as.factor(pheno.all.long$cohort)


########################
######### CPS ##########
########################

pheno.CPS <- pheno.all.long[pheno.all.long$cohort=="CPS",]
pheno.CPS <- pheno.CPS[pheno.CPS$dose!="3x8",]

pheno.CPS$time.point <- as.factor(pheno.CPS$time.point) 
pheno.CPS$time.point <- relevel(pheno.CPS$time.point, "i7")

# by timepoint
ggplot(pheno.CPS , aes(cell.subset, frequency, color=time.point)) +
  geom_boxplot(outlier.colour=NA, 
               position = position_dodge(width=0.9))+ 
  geom_point(aes(fill=time.point),position=position_jitterdodge(dodge.width=0.9),
             size = 0.5)+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.text = element_text(face="plain", size=9,lineheight=5.0),
        strip.background = element_rect(fill="white", colour="white",
                                        size=1))+
  scale_color_npg()+
  scale_y_log10()  + 
  ylab("% of viable PBMC") +
  xlab("Cell subsets") +
  labs(title = "")
ggsave("CPS_timepoint.pdf", width = 7, height = 5, units="in", dpi=600) 



# by protected and non-protected at post-immunization

pheno.CPS.c1 <- pheno.CPS[pheno.CPS$time.point=="c1",]

ggplot(pheno.CPS.c1 , aes(cell.subset, frequency, color=casecon)) +
  geom_boxplot(outlier.colour=NA, 
               position = position_dodge(width=0.9))+ 
  geom_point(aes(fill=casecon),position=position_jitterdodge(dodge.width=0.9),
             size = 0.5)+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.text = element_text(face="plain", size=9,lineheight=5.0),
        strip.background = element_rect(fill="white", colour="white",
                                        size=1))+
  scale_color_npg()+
  scale_y_log10()  + 
  ylab("% of viable PBMC") +
  xlab("Cell subsets") +
  labs(title = "")
ggsave("CPS_casecon_c1.pdf", width = 7, height = 5, units="in", dpi=600) 


# by protected and non-protected at baseline

pheno.CPS.i7 <- pheno.CPS[pheno.CPS$time.point=="i7",]

ggplot(pheno.CPS.i7 , aes(cell.subset, frequency, color=casecon)) +
  geom_boxplot(outlier.colour=NA, 
               position = position_dodge(width=0.9))+ 
  geom_point(aes(fill=casecon),position=position_jitterdodge(dodge.width=0.9),
             size = 0.5)+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.text = element_text(face="plain", size=9,lineheight=5.0),
        strip.background = element_rect(fill="white", colour="white",
                                        size=1))+
  scale_color_npg()+
  scale_y_log10()  + 
  ylab("% of viable PBMC") +
  xlab("Cell subsets") +
  labs(title = "")
ggsave("CPS_casecon_i7.pdf", width = 7, height = 5, units="in", dpi=600) 


#########################
######### RTS,S #########
#########################


pheno.all.long$agec <- as.factor(pheno.all.long$agec)
pheno.all.long$agec <-factor(pheno.all.long$agec, levels = c("6-12w","5-17m"))

# RTS,S vs comparators post-vaccination (M3)

pheno.rtss.vaccine.m3 <- pheno.all.long[pheno.all.long$cohort=="RTSS"  & pheno.all.long$time.point=="m3",]

ggplot(pheno.rtss.vaccine.m3 , aes(cell.subset, frequency, color=vaccine)) +
  geom_boxplot(outlier.colour=NA, 
               position = position_dodge(width=0.9))+ 
  geom_point(aes(fill=vaccine),position=position_jitterdodge(dodge.width=0.9),
             size = 0.5)+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.text = element_text(face="plain", size=9,lineheight=5.0),
        strip.background = element_rect(fill="white", colour="white",
                                        size=1))+
  scale_color_npg()+
  scale_y_log10()  + 
  ylab("% of viable PBMC") +
  xlab("Cell subsets") +
  labs(title = "")
ggsave("RTSS_vaccine_m3.pdf", width = 7, height = 5, units="in", dpi=600) 


# RTS,S vs comparators baseline (M3)
pheno.rtss.vaccine.m0 <- pheno.all.long[pheno.all.long$cohort=="RTSS"  & pheno.all.long$time.point=="m0",]

ggplot(pheno.rtss.vaccine.m0 , aes(cell.subset, frequency, color=vaccine)) +
  geom_boxplot(outlier.colour=NA, 
               position = position_dodge(width=0.9))+ 
  geom_point(aes(fill=vaccine),position=position_jitterdodge(dodge.width=0.9),
             size = 0.5)+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.text = element_text(face="plain", size=9,lineheight=5.0),
        strip.background = element_rect(fill="white", colour="white",
                                        size=1))+
  scale_color_npg()+
  scale_y_log10()  + 
  ylab("% of viable PBMC") +
  xlab("Cell subsets") +
  labs(title = "")
ggsave("RTSS_vaccine_m0.pdf", width = 7, height = 5, units="in", dpi=600) 


# case-control comparison within RTS,S vaccinees (post-vaccination)

pheno.rtss.m3 <- pheno.all.long[pheno.all.long$cohort=="RTSS" & pheno.all.long$vaccine=="rtss" & pheno.all.long$time.point=="m3",]

ggplot(pheno.rtss.m3, aes(cell.subset, frequency, color=casecon)) +
      geom_boxplot(outlier.colour=NA, 
             position = position_dodge(width=0.9))+ 
      geom_point(aes(fill=casecon),position=position_jitterdodge(dodge.width=0.9),
             size = 0.5)+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.text = element_text(face="plain", size=9,lineheight=5.0),
        strip.background = element_rect(fill="white", colour="white",
                                        size=1))+
  scale_color_npg()+
 scale_y_log10()  + 
  ylab("% of viable PBMC") +
  xlab("Cell subsets") +
  labs(title = "")
ggsave("RTSS_casecon_rtss_m3.pdf", width = 7, height = 5, units="in", dpi=600) 

# case-control comparison within RTS,S vaccinees (baseline)
pheno.rtss.m0 <- pheno.all.long[pheno.all.long$cohort=="RTSS" & pheno.all.long$vaccine=="rtss" & pheno.all.long$time.point=="m0",]

ggplot(pheno.rtss.m0, aes(cell.subset, frequency, color=casecon)) +
  geom_boxplot(outlier.colour=NA, 
               position = position_dodge(width=0.9))+ 
  geom_point(aes(fill=casecon),position=position_jitterdodge(dodge.width=0.9),
             size = 0.5)+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.text = element_text(face="plain", size=9,lineheight=5.0),
        strip.background = element_rect(fill="white", colour="white",
                                        size=1))+
  scale_color_npg()+
  scale_y_log10()  + 
  ylab("% of viable PBMC") +
  xlab("Cell subsets") +
  labs(title = "")
ggsave("RTSS_casecon_rtss_m0.pdf", width = 7, height = 5, units="in", dpi=600) 


# case-control comparison within comparators (post-vaccination)

pheno.comparator.m3 <- pheno.all.long[pheno.all.long$cohort=="RTSS" & pheno.all.long$vaccine=="comparator" & pheno.all.long$time.point=="m3",]


ggplot(pheno.comparator.m3, aes(cell.subset, frequency, color=casecon)) +
  geom_boxplot(outlier.colour=NA, 
               position = position_dodge(width=0.9))+ 
  geom_point(aes(fill=casecon),position=position_jitterdodge(dodge.width=0.9),
             size = 0.5)+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.text = element_text(face="plain", size=9,lineheight=5.0),
        strip.background = element_rect(fill="white", colour="white",
                                        size=1))+
  scale_color_npg()+
  scale_y_log10()  + 
  ylab("% of viable PBMC") +
  xlab("Cell subsets") +
  labs(title = "")
ggsave("RTSS_casecon_comparator_m3.pdf", width = 7, height = 5, units="in", dpi=600) 

# case-control comparison within comparators (baseline)

pheno.comparator.m0 <- pheno.all.long[pheno.all.long$cohort=="RTSS" & pheno.all.long$vaccine=="comparator" & pheno.all.long$time.point=="m0",]

ggplot(pheno.comparator.m0, aes(cell.subset, frequency, color=casecon)) +
  geom_boxplot(outlier.colour=NA, 
               position = position_dodge(width=0.9))+ 
  geom_point(aes(fill=casecon),position=position_jitterdodge(dodge.width=0.9),
             size = 0.5)+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.text = element_text(face="plain", size=9,lineheight=5.0),
        strip.background = element_rect(fill="white", colour="white",
                                        size=1))+
  scale_color_npg()+
  scale_y_log10()  + 
  ylab("% of viable PBMC") +
  xlab("Cell subsets") +
  labs(title = "")
ggsave("RTSS_casecon_comparator_m0.pdf", width = 7, height = 5, units="in", dpi=600) 



# Analysis stratified by AGE COHORT (Infants 6-12 weeks old and Children 5-17 months old)

## RTS,S vs Comparator comparison STRATIFIED BY AGE COHORT

### infants
pheno.rtss.vaccine.m3.infants <- pheno.all.long[pheno.all.long$cohort=="RTSS"  & pheno.all.long$time.point=="m3" & pheno.all.long$agec=="6-12w",]

ggplot(pheno.rtss.vaccine.m3.infants , aes(cell.subset, frequency, color=vaccine)) +
  geom_boxplot(outlier.colour=NA, 
               position = position_dodge(width=0.9))+ 
  geom_point(aes(fill=vaccine),position=position_jitterdodge(dodge.width=0.9),
             size = 0.5)+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.text = element_text(face="plain", size=9,lineheight=5.0),
        strip.background = element_rect(fill="white", colour="white",
                                        size=1))+
  scale_color_npg()+
  scale_y_log10()  + 
  ylab("% of viable PBMC") +
  xlab("Cell subsets") +
  labs(title = "")
ggsave("RTSS_vaccine_m3_infants.pdf", width = 7, height = 5, units="in", dpi=600) 


### children

pheno.rtss.vaccine.m3.children <- pheno.all.long[pheno.all.long$cohort=="RTSS"  & pheno.all.long$time.point=="m3" & pheno.all.long$agec=="5-17m",]

ggplot(pheno.rtss.vaccine.m3.children , aes(cell.subset, frequency, color=vaccine)) +
  geom_boxplot(outlier.colour=NA, 
               position = position_dodge(width=0.9))+ 
  geom_point(aes(fill=vaccine),position=position_jitterdodge(dodge.width=0.9),
             size = 0.5)+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.text = element_text(face="plain", size=9,lineheight=5.0),
        strip.background = element_rect(fill="white", colour="white",
                                        size=1))+
  scale_color_npg()+
  scale_y_log10()  + 
  ylab("% of viable PBMC") +
  xlab("Cell subsets") +
  labs(title = "")
ggsave("RTSS_vaccine_m3_children.pdf", width = 7, height = 5, units="in", dpi=600) 


# case-control comparison within RTS,S vaccinees STRATIFIED BY AGE COHORT

### infants
pheno.rtss.m3.infants <- pheno.all.long[pheno.all.long$cohort=="RTSS" & pheno.all.long$vaccine=="rtss" & pheno.all.long$time.point=="m3" & pheno.all.long$agec=="6-12w",]
dim(pheno.rtss.m3.infants)
dim(pheno.rtss.m3)

ggplot(pheno.rtss.m3.infants, aes(cell.subset, frequency, color=casecon)) +
  geom_boxplot(outlier.colour=NA, 
               position = position_dodge(width=0.9))+ 
  geom_point(aes(fill=casecon),position=position_jitterdodge(dodge.width=0.9),
             size = 0.5)+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.text = element_text(face="plain", size=9,lineheight=5.0),
        strip.background = element_rect(fill="white", colour="white",
                                        size=1))+
  scale_color_npg()+
  scale_y_log10(limits=c(0.01, 100))+ 
  ylab("% of viable PBMC") +
  xlab("Cell subsets") +
  labs(title = "")
ggsave("RTSS_casecon_rtss_m3_Infants.pdf", width = 7, height = 5, units="in", dpi=600) 

### children
pheno.rtss.m3.children <- pheno.all.long[pheno.all.long$cohort=="RTSS" & pheno.all.long$vaccine=="rtss" & pheno.all.long$time.point=="m3" & pheno.all.long$agec=="5-17m",]
dim(pheno.rtss.m3.children)

ggplot(pheno.rtss.m3.children, aes(cell.subset, frequency, color=casecon)) +
  geom_boxplot(outlier.colour=NA, 
               position = position_dodge(width=0.9))+ 
  geom_point(aes(fill=casecon),position=position_jitterdodge(dodge.width=0.9),
             size = 0.5)+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.text = element_text(face="plain", size=9,lineheight=5.0),
        strip.background = element_rect(fill="white", colour="white",
                                        size=1))+
  scale_color_npg()+
  scale_y_log10(limits=c(0.01, 100))  + 
  ylab("% of viable PBMC") +
  xlab("Cell subsets") +
  labs(title = "")
ggsave("RTSS_casecon_rtss_m3_Children.pdf", width = 7, height = 5, units="in", dpi=600) 



# case-control comparison within comparators STRATIFIED BY AGE COHORT

### infants

pheno.comparator.m3.infants <- pheno.all.long[pheno.all.long$cohort=="RTSS" & pheno.all.long$vaccine=="comparator" & pheno.all.long$time.point=="m3" & pheno.all.long$agec=="6-12w",]


ggplot(pheno.comparator.m3.infants, aes(cell.subset, frequency, color=casecon)) +
  geom_boxplot(outlier.colour=NA, 
               position = position_dodge(width=0.9))+ 
  geom_point(aes(fill=casecon),position=position_jitterdodge(dodge.width=0.9),
             size = 0.5)+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.text = element_text(face="plain", size=9,lineheight=5.0),
        strip.background = element_rect(fill="white", colour="white",
                                        size=1))+
  scale_color_npg()+
  scale_y_log10(limits=c(0.01, 100))  + 
  ylab("% of viable PBMC") +
  xlab("Cell subsets") +
  labs(title = "")
ggsave("RTSS_casecon_comparator_m3_Infants.pdf", width = 7, height = 5, units="in", dpi=600) 


### children

pheno.comparator.m3.children <- pheno.all.long[pheno.all.long$cohort=="RTSS" & pheno.all.long$vaccine=="comparator" & pheno.all.long$time.point=="m3" & pheno.all.long$agec=="5-17m",]

ggplot(pheno.comparator.m3.children, aes(cell.subset, frequency, color=casecon)) +
  geom_boxplot(outlier.colour=NA, 
               position = position_dodge(width=0.9))+ 
  geom_point(aes(fill=casecon),position=position_jitterdodge(dodge.width=0.9),
             size = 0.5)+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.text = element_text(face="plain", size=9,lineheight=5.0),
        strip.background = element_rect(fill="white", colour="white",
                                        size=1))+
  scale_color_npg()+
  scale_y_log10(limits=c(0.01, 100))  + 
  ylab("% of viable PBMC") +
  xlab("Cell subsets") +
  labs(title = "")
ggsave("RTSS_casecon_comparator_m3_Children.pdf", width = 7, height = 5, units="in", dpi=600) 
