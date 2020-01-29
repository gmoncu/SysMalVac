############################################################################################
#################################### CPS analysis ##########################################
############################################################################################

# This is code to perform the differential gene analysis of the CPS study reported in the manuscript
# Datasets of differential gene expression for the analysis of CPS immunogenicity and protection provided in the manuscript are based on this analysis
# We used the ranked genes obtained here for the gene set enrichment analyses based on differential gene expression 

# This is code to obtain also the differential gene expression of antigen stimulations over they background controls 
# separately in 
                # all subjects pre-immunization
                # all subjects post-immunization
                # protected subjects pre-immunization
                # non-protected subjects post-immunization
# The results obtained were used to perform the radar plots shown in the manuscript.




# Load R functions

sapply(list.files(path="/BIO/scripts/active/Affymetrix/Filtrado",pattern="[.]R$", full.names=TRUE), source)

sapply(list.files(path="/BIO/scripts/active/Affymetrix/ControlCalidad",pattern="[.]R$", full.names=TRUE), source)

sapply(list.files(path="/BIO/scripts/active/Affymetrix/Plot",pattern="[.]R$", full.names=TRUE), source)

sapply(list.files(path="/BIO/scripts/active/Affymetrix/statistics",pattern="[.]R$", full.names=TRUE), source)


#Packages 

pkgs<-c("ggplot2","GCally","CCA","affy","limma","annotate","hugene21sttranscriptcluster.db","biomaRt","GO.db","topGO","multtest","Rgraphviz","marray","amap","gplots")

lapply(pkgs, require, character.only=T)	  

options(warn=-1)

path=getwd()

names <- list.celfiles(path=path)

for (i in c(1:length(names))){
    j<-0
    while (substr(names[i],j,j)!=".")
      {
      j=j+1
      }
    names[i]<-substr(names[i],1,(j-1))
  }

#QC

myPData<-read.table("phenoDATA.txt", sep="",header=TRUE,as.is=TRUE)

rownames(myPData)<-myPData$sample


QualityControl(path=getwd(),method="RMA",control="Both",intgroup="stimulus")

# Remove 13SE1071  

QualityControl(path=getwd(),method="RMA",control="Both",intgroup="stimulus")

# Remove 13SE2260 

QualityControl(path=getwd(),method="RMA",control="Both",intgroup="stimulus")


listafinal<-read.table("listaTrabajoRma_lessStringent.txt") # table obtained in CPS_analysis1 

#10841 genes 190 samples

colnames(listafinal)<-names



########################################################################################################################
##################### Analysis of differential gene expression in CSP vs DMSO and PfRBC vs uRBC ######################
########################################################################################################################

### Analysis of differential gene expression of CSP vs DMSO and PfRBC vs uRBC in all subjects at baseline and post-immunization

system("mkdir -p Results")
system("mkdir -p Results/genelist1_new")
system("mkdir -p Results/genelist3_new")
system("mkdir -p Results/genelist3_new")
system("mkdir -p Results/genelist4_new")


myPData<-read.table("phenoDATA.txt", sep="",header=TRUE,as.is=TRUE)

stimulus<-factor(myPData$stimulus)

protected<-factor(myPData$protected)

group<-factor(myPData$group)

timepoint<-factor(myPData$timepoint)

donor<-factor(myPData$donor)

experimento <- factor(paste(stimulus,timepoint,sep="."))

design<-model.matrix(~0+experimento)

colnames(design) <- c(levels(experimento))

#Duplicate correction
corfit <- duplicateCorrelation(listafinal,design,block=donor)

#corfit$consensus=0.3559615

fit<- lmFit(listafinal,design,block=donor,correlation=corfit$consensus)


# contrasts

contrast.matrix <- makeContrasts(
  
  genelist1=(CS_peptide_pool.I_7 - DMSO_control.I_7), # CSP vs DMSO pre-immunization
  
  genelist2=(CS_peptide_pool.C_1 - DMSO_control.C_1), # CSP vs DMSO post-immunization
  
  genelist3=(PfRBC.I_7 - uRBC.I_7), # PfRBC vs uRBC pre-immunization
  
  genelist4=(PfRBC.C_1 - uRBC.C_1), # PfRBC vs uRBC post-immunization
  
  levels=design)
  
 
fit2<-eBayes(contrasts.fit(fit,contrast.matrix))


################################################################### Genelist1 ######################################################################################################

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist1",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist1_new_bonferroni"

write.table(genes_bonferroni[,c(1,2,5,6)],file=paste(paste("Results/genelist1_new",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist1",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist1_new_fdr"

write.table(genes_fdr[,c(1,2,5,6)],file=paste(paste("Results/genelist1_new",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("ID","logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist1_new_total"

write.table(genes_total,file=paste(paste("Results/genelist1_new",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist1_new/genelist1_new.RData")



################################################################### Genelist2 ######################################################################################################


adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist2",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist2_new_bonferroni"

write.table(genes_bonferroni[,c(1,2,5,6)],file=paste(paste("Results/genelist2_new",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist2",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist2_new_fdr"

write.table(genes_fdr[,c(1,2,5,6)],file=paste(paste("Results/genelist2_new",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("ID","logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist2_new_total"

write.table(genes_total,file=paste(paste("Results/genelist2_new",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist2_new/genelist2_new.RData")


################################################################### Genelist3 ######################################################################################################

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist3",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist3_new_bonferroni"

write.table(genes_bonferroni[,c(1,2,5,6)],file=paste(paste("Results/genelist3_new",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist3",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist3_new_fdr"

write.table(genes_fdr[,c(1,2,5,6)],file=paste(paste("Results/genelist3_new",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("ID","logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist3_new_total"

write.table(genes_total,file=paste(paste("Results/genelist3_new",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist3_new/genelist3_new.RData")

################################################################### Genelist4 ######################################################################################################

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist4",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist4_new_bonferroni"

write.table(genes_bonferroni[,c(1,2,5,6)],file=paste(paste("Results/genelist4_new",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist4",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist4_new_fdr"

write.table(genes_fdr[,c(1,2,5,6)],file=paste(paste("Results/genelist4_new",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("ID","logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist4_new_total"

write.table(genes_total,file=paste(paste("Results/genelist4_new",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist4_new/genelist4_new.RData")



### Analysis of differential gene expression of CSP vs DMSO and PfRBC vs uRBC in protected and in non-protected subjects at both timepoints


system("mkdir -p Results/genelist5_new")
system("mkdir -p Results/genelist6_new")
system("mkdir -p Results/genelist7_new")
system("mkdir -p Results/genelist8_new")
system("mkdir -p Results/genelist9_new")
system("mkdir -p Results/genelist10_new")
system("mkdir -p Results/genelist11_new")
system("mkdir -p Results/genelist12_new")



experimento <- factor(paste(stimulus,timepoint,protected,sep="."))

design<-model.matrix(~0+experimento)

colnames(design) <- c(levels(experimento))

#Duplicate correction
corfit <- duplicateCorrelation(listafinal,design,block=donor)

#corfit$consensus=0.3492553

fit<- lmFit(listafinal,design,block=donor,correlation=corfit$consensus)

# contrasts for the comparisons of CSP vs DMSO and PfRBC vs uRBC in protected and non-protected subjects before and after immunization

contrast.matrix <- makeContrasts(
  
  genelist5=(CS_peptide_pool.C_1.yes - DMSO_control.C_1.yes), # CSP vs DMSO in protected CPS subjects post-immunization
  
  genelist6=(CS_peptide_pool.C_1.no - DMSO_control.C_1.no), # CSP vs DMSO  in non-protected CPS subjects post-immunization
  
  genelist7=(PfRBC.C_1.yes - uRBC.C_1.yes), # PfRBC vs uRBC in protected CPS subjects post-immunization
  
  genelist8=(PfRBC.C_1.no - uRBC.C_1.no), # PfRBC vs uRBC in non-protected CPS subjects post-immunization
  
  genelist9=(CS_peptide_pool.I_7.yes - DMSO_control.I_7.yes), # CSP vs DMSO in protected  subjects pre-immunization
  
  genelist10=(CS_peptide_pool.I_7.no - DMSO_control.I_7.no), # CSP vs DMSO  in non-protected subjects pre-stimmunization
  
  genelist11=(PfRBC.I_7.yes - uRBC.I_7.yes), # PfRBC vs uRBC in protected subjects pre-stimmunization
  
  genelist12=(PfRBC.I_7.no -uRBC.I_7.no ), # PfRBC vs uRBC in non-protected subjects pre-stimmunization
  
  levels=design)
  
 
fit2<-eBayes(contrasts.fit(fit,contrast.matrix))


################################################################### Genelist5 ######################################################################################################


adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist5",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist5_new_bonferroni"

write.table(genes_bonferroni[,c(1,2,5,6)],file=paste(paste("Results/genelist5_new",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist5",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist5_new_fdr"

write.table(genes_fdr[,c(1,2,5,6)],file=paste(paste("Results/genelist5_new",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("ID","logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist5_new_total"

write.table(genes_total,file=paste(paste("Results/genelist5_new",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist5_new/genelist5_new.RData")

################################################################### Genelist6 ######################################################################################################


adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist6",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist6_new_bonferroni"

write.table(genes_bonferroni[,c(1,2,5,6)],file=paste(paste("Results/genelist6_new",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist6",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist6_new_fdr"

write.table(genes_fdr[,c(1,2,5,6)],file=paste(paste("Results/genelist6_new",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("ID","logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist6_new_total"

write.table(genes_total,file=paste(paste("Results/genelist6_new",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist6_new/genelist6_new.RData")

################################################################### Genelist7 ######################################################################################################

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist7",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist7_new_bonferroni"

write.table(genes_bonferroni[,c(1,2,5,6)],file=paste(paste("Results/genelist7_new",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist7",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist7_new_fdr"

write.table(genes_fdr[,c(1,2,5,6)],file=paste(paste("Results/genelist7_new",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("ID","logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist7_new_total"

write.table(genes_total,file=paste(paste("Results/genelist7_new",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist7_new/genelist7_new.RData")


################################################################### Genelist8 ######################################################################################################

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist8",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist8_new_bonferroni"

write.table(genes_bonferroni[,c(1,2,5,6)],file=paste(paste("Results/genelist8_new",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist8",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist8_new_fdr"

write.table(genes_fdr[,c(1,2,5,6)],file=paste(paste("Results/genelist8_new",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("ID","logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist8_new_total"

write.table(genes_total,file=paste(paste("Results/genelist8_new",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist8_new/genelist8_new.RData")



################################################################### Genelist9 ######################################################################################################


adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist9",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist9_new_bonferroni"

write.table(genes_bonferroni[,c(1,2,5,6)],file=paste(paste("Results/genelist9_new",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist9",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist9_new_fdr"

write.table(genes_fdr[,c(1,2,5,6)],file=paste(paste("Results/genelist9_new",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("ID","logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist9_new_total"

write.table(genes_total,file=paste(paste("Results/genelist9_new",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist9_new/genelist9_new.RData")


################################################################### Genelist10 ######################################################################################################


adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist10",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist10_new_bonferroni"

write.table(genes_bonferroni[,c(1,2,5,6)],file=paste(paste("Results/genelist10_new",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist10",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist10_new_fdr"

write.table(genes_fdr[,c(1,2,5,6)],file=paste(paste("Results/genelist10_new",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("ID","logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist10_new_total"

write.table(genes_total,file=paste(paste("Results/genelist10_new",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist10_new/genelist10_new.RData")

################################################################### Genelist11 ######################################################################################################

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist11",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist11_new_bonferroni"

write.table(genes_bonferroni[,c(1,2,5,6)],file=paste(paste("Results/genelist11_new",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist11",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist11_new_fdr"

write.table(genes_fdr[,c(1,2,5,6)],file=paste(paste("Results/genelist11_new",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("ID","logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist11_new_total"

write.table(genes_total,file=paste(paste("Results/genelist11_new",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist11_new/genelist11_new.RData")


################################################################### Genelist12 ######################################################################################################


adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist12",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist12_new_bonferroni"

write.table(genes_bonferroni[,c(1,2,5,6)],file=paste(paste("Results/genelist12_new",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist12",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist12_new_fdr"

write.table(genes_fdr[,c(1,2,5,6)],file=paste(paste("Results/genelist12_new",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("ID","logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist12_new_total"

write.table(genes_total,file=paste(paste("Results/genelist12_new",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist12_new.RData")



########################################################################################################################################################################################
##################################################### CPS immunogenicity ##############################################################################################
########################################################################################################################################################################################

system("mkdir -p Results/genelist13_new")
system("mkdir -p Results/genelist14_new")


experimento <- factor(paste(stimulus,timepoint,sep="."))

design<-model.matrix(~0+experimento)

colnames(design) <- c(levels(experimento))

#Corrección de duplicados
corfit <- duplicateCorrelation(listafinal,design,block=donor)

#corfit$consensus=0.3559615

fit<- lmFit(listafinal,design,block=donor,correlation=corfit$consensus)

# contrasts

contrast.matrix <- makeContrasts(
  
  genelist13=(CS_peptide_pool.C_1 - CS_peptide_pool.I_7)-(DMSO_control.C_1 - DMSO_control.I_7), # Post- vs pre-immunization for CSP stimulations
  
  genelist14=(PfRBC.C_1 - PfRBC.I_7)-(uRBC.C_1 - uRBC.I_7), # Post- vs pre-immunization for PfRBC stimulations
  
  levels=design)


fit2<-eBayes(contrasts.fit(fit,contrast.matrix))


################################################################### Genelist13 ######################################################################################################


adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist13",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist13_new_bonferroni"

write.table(genes_bonferroni[,c(1,2,5,6)],file=paste(paste("Results/genelist13_new",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist13",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist13_new_fdr"

write.table(genes_fdr[,c(1,2,5,6)],file=paste(paste("Results/genelist13_new",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("ID","logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist13_new_total"

write.table(genes_total,file=paste(paste("Results/genelist13_new",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist13_new/genelist13_new.RData")


################################################################### Genelist14 ######################################################################################################

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist14",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist14_new_bonferroni"

write.table(genes_bonferroni[,c(1,2,5,6)],file=paste(paste("Results/genelist14_new",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist14",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist14_new_fdr"

write.table(genes_fdr[,c(1,2,5,6)],file=paste(paste("Results/genelist14_new",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("ID","logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist14_new_total"

write.table(genes_total,file=paste(paste("Results/genelist14_new",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist14_new/genelist14_new.RData")


########################################################################################################################################################################################
##################################################### CPS protection ##############################################################################################
########################################################################################################################################################################################


system("mkdir -p Results/genelist15_new")
system("mkdir -p Results/genelist16_new")
system("mkdir -p Results/genelist17_new")
system("mkdir -p Results/genelist18_new")

experimento <- factor(paste(stimulus,timepoint,protected,sep="."))

design<-model.matrix(~0+experimento)

colnames(design) <- c(levels(experimento))

#Corrección de duplicados
corfit <- duplicateCorrelation(listafinal,design,block=donor)

#corfit$consensus=0.3492553

fit<- lmFit(listafinal,design,block=donor,correlation=corfit$consensus)

# contrasts

contrast.matrix <- makeContrasts(
  
  genelist15=(CS_peptide_pool.C_1.yes - CS_peptide_pool.C_1.no)-(DMSO_control.C_1.yes - DMSO_control.C_1.no), # CPS protection post-immunization for CSP stimulation
  
  genelist16=(PfRBC.C_1.yes - PfRBC.C_1.no)-(uRBC.C_1.yes - uRBC.C_1.no), # CPS protection post-immunization for PfRBC stimulation
  
  genelist17=(CS_peptide_pool.I_7.yes - CS_peptide_pool.I_7.no)-(DMSO_control.I_7.yes - DMSO_control.I_7.no), # CPS protection pre-immunization for CSP stimulation
  
  genelist18=(PfRBC.I_7.yes - PfRBC.I_7.no)-(uRBC.I_7.yes - uRBC.I_7.no), # CPS protection pre-immunization for PfRBC stimulation
  
  levels=design)


fit2<-eBayes(contrasts.fit(fit,contrast.matrix))




################################################################### Genelist15 ######################################################################################################


adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist15",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist15_new_bonferroni"

write.table(genes_bonferroni[,c(1,2,5,6)],file=paste(paste("Results/genelist15_new",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist15",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist15_new_fdr"

write.table(genes_fdr[,c(1,2,5,6)],file=paste(paste("Results/genelist15_new",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("ID","logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist15_new_total"

write.table(genes_total,file=paste(paste("Results/genelist15_new",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist15_new/genelist15_new.RData")


################################################################### Genelist16 ######################################################################################################


adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist16",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist16_new_bonferroni"

write.table(genes_bonferroni[,c(1,2,5,6)],file=paste(paste("Results/genelist16_new",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist16",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist16_new_fdr"

write.table(genes_fdr[,c(1,2,5,6)],file=paste(paste("Results/genelist16_new",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("ID","logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist16_new_total"

write.table(genes_total,file=paste(paste("Results/genelist16_new",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist16_new/genelist16_new.RData")



################################################################### Genelist17 ######################################################################################################


adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist17",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist17_new_bonferroni"

write.table(genes_bonferroni[,c(1,2,5,6)],file=paste(paste("Results/genelist17_new",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist17",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist17_new_fdr"

write.table(genes_fdr[,c(1,2,5,6)],file=paste(paste("Results/genelist17_new",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("ID","logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist17_new_total"

write.table(genes_total,file=paste(paste("Results/genelist17_new",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist17_new/genelist17_new.RData")



################################################################### Genelist18 ######################################################################################################


adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="genelist18",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist18_new_bonferroni"

write.table(genes_bonferroni[,c(1,2,5,6)],file=paste(paste("Results/genelist18_new",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="genelist18",sort.by="P",number=dim(listafinal)[1])

fileout<-"genelist18_new_fdr"

write.table(genes_fdr[,c(1,2,5,6)],file=paste(paste("Results/genelist18_new",fileout,sep="/"),".txt",sep=""))

genes_total<-cbind(genes_fdr[,c(1,2,5,6)],genes_bonferroni[,c(6)])

colnames(genes_total)<-c("ID","logFC","Pvalue","FDR_Pvalue","Bonferroni_Pvalue")

fileout<-"genelist18_new_total"

write.table(genes_total,file=paste(paste("Results/genelist18_new",fileout,sep="/"),".txt",sep=""))

save(list=ls(all=TRUE),file="Results/genelist18_new/genelist18_new.RData")





###############################################################################################################################################################################
######################################################################### end #################################################################################################
###############################################################################################################################################################################
