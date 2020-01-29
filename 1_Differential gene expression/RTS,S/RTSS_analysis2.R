##############################################################################################################
####################   RTS,S analysis (DGE in CSP vs DMSO and PfRBC vs uRBC) #################################
##############################################################################################################

# This is code to obtain the differential gene expression of antigen stimulations over they background controls 
# separately in 
                #  protected and non-protected RTS,S vaccinated subjects
                # protected and non-protected Comparator vaccinated subjects
                # RTS,S vaccinated subjects and in the Comparator vaccinated subjects

# The results obtained were used to perform the radar plots shown in the manuscript.


# Load functions

sapply(list.files(path="/BIO/scripts/active/Affymetrix/Filtrado",pattern="[.]R$", full.names=TRUE), source)

sapply(list.files(path="/BIO/scripts/active/Affymetrix/ControlCalidad",pattern="[.]R$", full.names=TRUE), source)

sapply(list.files(path="/BIO/scripts/active/Affymetrix/Plot",pattern="[.]R$", full.names=TRUE), source)

sapply(list.files(path="/BIO/scripts/active/Affymetrix/statistics",pattern="[.]R$", full.names=TRUE), source)


#Packages

pkgs<-c("ggplot2","GGally","CCA","oligo","biomaRt","topGO","multtest","Rgraphivz","hugene21sttranscriptcluster.db","annotate",
"GO.db""statmod","oligoData","gplots","arrayQualityMetrics","hugene21stcdf","genefilter","affy","limma","marray","amap")

lapply(pkgs, require, character.only=T)	 

options(warn=-1)


# Load data from RTS,S main analysis (RTSS_analysis1)

load("Results/SysMalVac_December.RData")

system("mkdir -p Results/genelist1_March16")
system("mkdir -p Results/genelist2_March16")
system("mkdir -p Results/genelist3_March16")
system("mkdir -p Results/genelist4_March16")
system("mkdir -p Results/genelist5_March16")
system("mkdir -p Results/genelist6_March16")
system("mkdir -p Results/genelist7_March16")
system("mkdir -p Results/genelist8_March16")
system("mkdir -p Results/genelist9_March16")
system("mkdir -p Results/genelist10_March16")
system("mkdir -p Results/genelist11_March16")
system("mkdir -p Results/genelist12_March16")

# contrasts

contrast.matrix <- makeContrasts(
  
  genelist1=((non_malariaCSPM31-non_malariaDMSOM31)), # CSP vs DMSO in protected RTS,S vaccinees 
  
  genelist2=((malariaCSPM31-malariaDMSOM31)), # CSP vs DMSO in non-protected RTS,S vaccinees 
  
  genelist5=((malariaCSPM32-malariaDMSOM32)), # CSP vs DMSO in non-protected Comparator vaccinees 
  
  genelist6=((non_malariaCSPM32-non_malariaDMSOM32)), # CSP vs DMSO in protected Comparator vaccinees 
    
  genelist7=((malariaiRBCM32-malariauRBCM32)), # PfRBC vs uRBC in non-protected Comparator vaccinees 
  
  genelist8=((non_malariaiRBCM32-non_malariauRBCM32)), # PfRBC vs uRBC in protected Comparator vaccinees 
  
  genelist9=((non_malariaiRBCM31-non_malariauRBCM31)), # PfRBC vs uRBC in protected RTS,S vaccinees 
  
  genelist10=((malariaiRBCM31-malariauRBCM31)), # PfRBC vs uRBC in non-protected RTS,S vaccinees 
  
  levels=design)
  
fit2<-eBayes(contrasts.fit(fit,contrast.matrix))



save(contrast.matrix,fit,design,myPData2,fit2,exprslistafinal_04,listafinal_04,file="Results/SysMalVac_March16.RData")

############################################################## Genelist1: CSP vs DMSO in protected RTS,S vaccinees ###########################################################################################################

load("Results/SysMalVac_March16.RData")

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist1",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist1_March16_bonferroni"

write.table(genes_bonferroni[,c(2,5,6)],file=paste(paste("Results/genelist1_March16",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist1",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist1_March16_fdr"

write.table(genes_fdr[,c(2,5,6)],file=paste(paste("Results/genelist1_March16",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist1_March16_total"

write.table(genes_total,file=paste(paste("Results/genelist1_March16",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist1_March16/genelist1_March16.RData")



############################################################## Genelist2: CSP vs DMSO in non-protected vaccinees ###########################################################################################################

load("Results/SysMalVac_March16.RData")

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist2",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist2_March16_bonferroni"

write.table(genes_bonferroni[,c(2,5,6)],file=paste(paste("Results/genelist2_March16",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist2",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist2_March16_fdr"

write.table(genes_fdr[,c(2,5,6)],file=paste(paste("Results/genelist2_March16",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist2_March16_total"

write.table(genes_total,file=paste(paste("Results/genelist2_March16",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist2_March16/genelist2_March16.RData")



################################################ Genelist5: CSP vs DMSO in non-protected comparator vaccinees ###########################################################################################################

load("Results/SysMalVac_March16.RData")

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist5",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist5_March16_bonferroni"

write.table(genes_bonferroni[,c(2,5,6)],file=paste(paste("Results/genelist5_March16",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist5",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist5_March16_fdr"

write.table(genes_fdr[,c(2,5,6)],file=paste(paste("Results/genelist5_March16",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist5_March16_total"

write.table(genes_total,file=paste(paste("Results/genelist5_March16",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist5_March16/genelist5_March16.RData")


############################################################## Genelist6: CSP vs DMSO in protected comparator vaccinees ###########################################################################################################

load("Results/SysMalVac_March16.RData")

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist6",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist6_March16_bonferroni"

write.table(genes_bonferroni[,c(2,5,6)],file=paste(paste("Results/genelist6_March16",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist6",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist6_March16_fdr"

write.table(genes_fdr[,c(2,5,6)],file=paste(paste("Results/genelist6_March16",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist6_March16_total"

write.table(genes_total,file=paste(paste("Results/genelist6_March16",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist6_March16/genelist6_March16.RData")



############################################################## Genelist7: iRBC vs uRBC non-protected comparator vaccinees ###########################################################################################################

load("Results/SysMalVac_March16.RData")

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist7",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist7_March16_bonferroni"

write.table(genes_bonferroni[,c(2,5,6)],file=paste(paste("Results/genelist7_March16",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist7",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist7_March16_fdr"

write.table(genes_fdr[,c(2,5,6)],file=paste(paste("Results/genelist7_March16",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist7_March16_total"

write.table(genes_total,file=paste(paste("Results/genelist7_March16",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist7_March16/genelist7_March16.RData")



############################################################## Genelist8: PfRBC vs uRBC in protected Comparator vaccinees ###########################################################################################################

load("Results/SysMalVac_March16.RData")

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist8",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist8_March16_bonferroni"

write.table(genes_bonferroni[,c(2,5,6)],file=paste(paste("Results/genelist8_March16",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist8",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist8_March16_fdr"

write.table(genes_fdr[,c(2,5,6)],file=paste(paste("Results/genelist8_March16",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist8_March16_total"

write.table(genes_total,file=paste(paste("Results/genelist8_March16",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist8_March16/genelist8_March16.RData")


############################################################## Genelist9: PfRBC vs uRBC in protected RTS,S vaccinees ###########################################################################################################

load("Results/SysMalVac_March16.RData")

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist9",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist9_March16_bonferroni"

write.table(genes_bonferroni[,c(2,5,6)],file=paste(paste("Results/genelist9_March16",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist9",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist9_March16_fdr"

write.table(genes_fdr[,c(2,5,6)],file=paste(paste("Results/genelist9_March16",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist9_March16_total"

write.table(genes_total,file=paste(paste("Results/genelist9_March16",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist9_March16/genelist9_March16.RData")



############################################################## Genelist10: PfRBC vs uRBC in non-protected RTS,S vaccinees  ###########################################################################################################

load("Results/SysMalVac_March16.RData")

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist10",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist10_March16_bonferroni"

write.table(genes_bonferroni[,c(2,5,6)],file=paste(paste("Results/genelist10_March16",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist10",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist10_March16_fdr"

write.table(genes_fdr[,c(2,5,6)],file=paste(paste("Results/genelist10_March16",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist10_March16_total"

write.table(genes_total,file=paste(paste("Results/genelist10_March16",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist10_March16/genelist10_March16.RData")






# Continuation
# More contrasts

load("Results/SysMalVac_December_2comp.RData")

contrast.matrix <- makeContrasts(
  
  genelist3=((CSPM31-DMSOM31)), # CSP vs uRBC in all RTS,S vaccinees
  
  genelist4=((CSPM32-DMSOM32)), # CSP vs uRBC in all comparator vaccinees
  
  genelist11=((iRBCM31-uRBCM31)), # PfRBC vs uRBC in all RTS,S vaccinees
  
  genelist12=((iRBCM32-uRBCM32)), # PfRBC vs uRBC in all comparator RTS,S vaccinees
  
  levels=design)
  
fit2<-eBayes(contrasts.fit(fit,contrast.matrix))

save(contrast.matrix,fit,design,myPData2,fit2,exprslistafinal_04,listafinal_04,file="Results/SysMalVac_110316_2comp.RData")



############################################################## Genelist3:  CSP vs uRBC in all RTS,S vaccinees ###########################################################################################################

load("Results/SysMalVac_110316_2comp.RData")


adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist3",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist3_March16_bonferroni"

write.table(genes_bonferroni[,c(1,4,5)],file=paste(paste("Results/genelist3_March16",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist3",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist3_March16_fdr"

write.table(genes_fdr[,c(1,4,5)],file=paste(paste("Results/genelist3_March16",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,4,5)],genes_bonferroni[,c(5)])

colnames(genes_total)<-c("logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist3_March16_total"

write.table(genes_total,file=paste(paste("Results/genelist3_March16",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist3_March16/genelist3_March16.RData")



############################################################## Genelist4:  CSP vs uRBC in all comparator vaccinees ###########################################################################################################

load("Results/SysMalVac_110316_2comp.RData")

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist4",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist4_March16_bonferroni"

write.table(genes_bonferroni[,c(1,4,5)],file=paste(paste("Results/genelist4_March16",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist4",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist4_March16_fdr"

write.table(genes_fdr[,c(1,4,5)],file=paste(paste("Results/genelist4_March16",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,4,5)],genes_bonferroni[,c(5)])

colnames(genes_total)<-c("logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist4_March16_total"

write.table(genes_total,file=paste(paste("Results/genelist4_March16",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist4_March16/genelist4_March16.RData")



############################################################## Genelist11:  PfRBC vs uRBC in all RTS,S vaccinees ###########################################################################################################

load("Results/SysMalVac_110316_2comp.RData")

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist11",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist11_March16_bonferroni"

write.table(genes_bonferroni[,c(1,4,5)],file=paste(paste("Results/genelist11_March16",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist11",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist11_March16_fdr"

write.table(genes_fdr[,c(1,4,5)],file=paste(paste("Results/genelist11_March16",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,4,5)],genes_bonferroni[,c(5)])

colnames(genes_total)<-c("logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist11_March16_total"

write.table(genes_total,file=paste(paste("Results/genelist11_March16",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist11_March16/genelist11_March16.RData")


############################################################## Genelist12:  PfRBC vs uRBC in all comparator vaccinees ###########################################################################################################

load("Results/SysMalVac_110316_2comp.RData")

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist12",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist12_March16_bonferroni"

write.table(genes_bonferroni[,c(1,4,5)],file=paste(paste("Results/genelist12_March16",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist12",sort.by="P",number=dim(exprslistafinal_04)[1])

fileout<-"genelist12_March16_fdr"

write.table(genes_fdr[,c(1,4,5)],file=paste(paste("Results/genelist12_March16",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,4,5)],genes_bonferroni[,c(5)])

colnames(genes_total)<-c("logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist12_March16_total"

write.table(genes_total,file=paste(paste("Results/genelist12_March16",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist12_March16/genelist12_March16.RData")

###################################################################################################################################################################################################################
###################################################################################################################################################################################################################
###################################################################################################################################################################################################################
