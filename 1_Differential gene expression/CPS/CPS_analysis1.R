############################################################################################
#################################### CPS analysis ##########################################
############################################################################################

# This is code to perform the microarray data analysis of the CPS study
# it includes normalization using RMA method, filtering of background noise and invariant genes 

# Load functions
sapply(list.files(path="utils/Filter",pattern="[.]R$", full.names=TRUE), source)

sapply(list.files(path="utils/QualityControl",pattern="[.]R$", full.names=TRUE), source)

sapply(list.files(path="utils/Plot",pattern="[.]R$", full.names=TRUE), source)

sapply(list.files(path="utils/statistic",pattern="[.]R$", full.names=TRUE), source)

# Packages
install.packages("qtl_1.30-4.tar.gz")
install.packages("hugene21stcdf_1.0.0.tar.gz")

pkgs<-c("ggplot2","oligo","statmod","oligoData","gplots","arrayQualityMetrics","hugene21stcdf","genefilter","affy","limma")

lapply(pkgs, require, character.only=T)	  

# Load data
raw.data=ReadAffy(verbose=TRUE, filenames=list.celfiles())
dat <- read.celfiles(list.celfiles())
eset <- rma(dat)

# normalization
rma.table=exprs(eset)
dim(rma.table)
#[1] 53617   192
write.table(rma.table, file="fullRMAdata_192.txt", sep="\t",row.names=FALSE)

#[1] 53617   192
rma.table[1:5,1:5]
#         13SE1059.CEL 13SE1060.CEL 13SE1063.CEL 13SE1064.CEL 13SE1065.CEL
#16650001     1.626163     2.140718     2.340795     1.125010     3.015840
#16650003     3.287487     2.748844     2.533877     3.949070     2.618608
#16650005     2.730029     4.848114     4.314844     2.542559     5.008024
#16650007     4.201780     4.510551     5.031684     3.682314     3.573254
#16650009     1.472716     1.121580     1.253914     1.515542     1.326307


write.table(file="nombres_cels.txt", cel.names, sep="\n")



arrayQualityMetrics(expressionset=data.rma.norm,outdir="QC_rma",split.plots=50)
#13SE1071 is an outlier. Remove from directory and rerun

dat <- read.celfiles(list.celfiles())
cels = list.files("/media/data/Sysmalvac_192/cels", pattern = "CEL")
sample_names<-gsub(".CEL", "", cels)
length(sample_names)
#[1] 191

sampleNames(dat)<-sample_names

# Load phenoData (CPS study variables and subjects characteristics)
myPData<-read.table("CPS/phenoData.txt", sep="\t",header=TRUE,row.names=1,as.is=TRUE) 
dim(myPData)
#[1] 192   6

#remove sample outlier 13SE1071
myPData <- myPData[!rownames(myPData)=="13SE1071",]

dim(myPData)
#[1] 191   6
head(myPData)
#           donor timepoint group protected     stimulus stimulus_time
#13SE1059 1169_18       I_7     3        no DMSO_control           24h
#13SE1060 1169_18       C_1     3        no DMSO_control           24h
#13SE1063 1317_38       I_7     1       yes DMSO_control           24h
#13SE1064 1317_38       C_1     1       yes DMSO_control           24h
#13SE1065 1296_50       I_7     3       yes DMSO_control           24h
#13SE1066 1296_50       C_1     3       yes DMSO_control           24h

# Group is dose group, timepoint is I_7 (pre-immunization) and C_1 (post-immunization)
 pData(dat)<-myPData

eset <- rma(dat)

rma.table=exprs(eset)
dim(rma.table)
#[1] 53617   191

arrayQualityMetrics(expressionset=eset,outdir="QC_rma_191",split.plots=50,intgroup="stimulus")
#The directory 'QC_rma_191' has been created.

#remove sample in position 189 from loaded files. Keep same phenoData

dat_190<-dat[,-189]

#Normalize
eset_190 <- rma(dat_190)
rma.table_190=exprs(eset_190)

write.table(rma.table_190, file="fullRMAdata_190.txt", sep="\t",row.names=FALSE)

#Filter 

database <- db(pd.hugene.2.1.st)
IDscontrol<- dbGetQuery(database, "select fsetid from featureSet where type in ('2','3','4','5','6','7','9');")[,1]
IDs<-as.character(IDscontrol)
rma_190.filtrado = featureFilter(eset_190,
					     require.entrez=FALSE, 
					     require.GOBP=FALSE, 
					     require.GOCC=FALSE, 
					     require.GOMF=FALSE, 
					     remove.dupEntrez=FALSE,
					     feature.exclude=IDs)



myControls<-read.table("/media/data/Sysmalvac_192/unnanotated_controls.txt", sep="\t",header=FALSE,as.is=TRUE)
dim(myControls)
#[1] 5849    1
rma_190.filtrado = featureFilter(eset_190,
      require.entrez=FALSE, 
      require.GOBP=FALSE, 
      require.GOCC=FALSE, 
      require.GOMF=FALSE, 
      remove.dupEntrez=FALSE,
      feature.exclude=myControls[,1])


#Filter intensity log2(100) and IQR. Begin IQR>0.5
f1=pOverA(1, log2(100))
f2=function(x) (IQR(x)>0.5)
ff=filterfun(f1,f2)
pasan_filtros=genefilter(rma_190.filtrado,ff)

# Modify IQR

f2=function(x) (IQR(x)>0.2)
ff=filterfun(f1,f2)
pasan_filtros=genefilter(rma_190.filtrado,ff)
listaTrabajoRma=rma_190.filtrado[pasan_filtros,]

#normalized expression table

listaTrabajoRma.table=exprs(listaTrabajoRma)

#plotCluster
system("mkdir Results")
plotCluster(listaTrabajoRma.table,path=getwd(),as.factor(phenoData(listaTrabajoRma)$protected),as.factor(phenoData(listaTrabajoRma)$stimulus),class.legend=c("CS_peptide_pool","DMSO_control","PfRBC","uRBC","no","yes"),fileout="cluster.listaTrabajo",pos.legend="topleft",horiz=TRUE)

#PCA
plotPCA_prcomp(eset=listaTrabajoRma.table,path=getwd(),fileout="PCA.listaTrabajo",pos.legend="topleft",nPCs=3)


# test other filters

save(rma_190.filtrado,file="SysMalVac_190_cambioFiltros.RData")


#intensity filter log2(100) and IQR. comienzo con IQR>0.2
f1=pOverA(0.25, log2(100))
f2=function(x) (IQR(x)>0.2)
ff=filterfun(f1,f2)
pasan_filtros=genefilter(rma_190.filtrado,ff)
sum(pasan_filtros[pasan_filtros=="TRUE"])
#[1] 10841
listaTrabajoRma=rma_190.filtrado[pasan_filtros,]
rma.table=exprs(listaTrabajoRma)
write.table(rma.table, file=" listaTrabajoRma_lessStringent.txt ", sep="\t",row.names=FALSE)


# Differential gene expression for CSP vs DMSO and PfRBC vs uRBC post-immunization for group dose 1 
# this analysis was performed to check if a more relaxed filter would lead to detection of some genes detected in set up experiments with samples of subjects immunized with this same dose.
contrast.sinI7.paired <- makeCo

#limma
tratamientoANDgroup <- factor(paste(listaTrabajoRma@phenoData$stimulus,listaTrabajoRma@phenoData$timepoint,listaTrabajoRma@phenoData$group,sep="."))
design.sinI7.paired <- model.matrix(~0+tratamientoANDgroup)
colnames(design.sinI7.paired) <- levels(tratamientoANDgroup)
corfit <- duplicateCorrelation(listaTrabajoRma,design.sinI7.paired,block=listaTrabajoRma@phenoData$donor)
fit.sinI7.paired <- lmFit(listaTrabajoRma,design.sinI7.paired,block=listaTrabajoRma@phenoData$donor,correlation=corfit$consensus)

#contrast for CSP vs DMSO and PfRBC vs uRBC post-immunization for group dose 1

contrast.sinI7.paired <- makeContrasts(
  CS_peptide_poolvsDMSO_control_sinI7_group1.paired = CS_peptide_pool.C_1.1-DMSO_control.C_1.1,
  PfRBCvsuRBC_sinI7_group1.paired = PfRBC.C_1.1-uRBC.C_1.1,
  levels=design.sinI7.paired)

fit2.sinI7.paired <- contrasts.fit(fit.sinI7.paired, contrast.sinI7.paired)
fit2.sinI7.paired <- eBayes(fit2.sinI7.paired)
xxx<-topTable(fit2.sinI7.paired, coef="CS_peptide_poolvsDMSO_control_sinI7_group1.paired", adjust = "BH",sort.by="P",number=6431);signiCS_peptide_poolsDMSO_control_sinI7_group1.paired<-xxx[xxx$P.Value<0.05,]

xxx<-topTable(fit2.sinI7.paired, coef="PfRBCvsuRBC_sinI7_group1.paired", adjust = "BH",sort.by="P",number=6431);signiPfRBCvsuRBC_sinI7_group1.paired<-xxx[xxx$P.Value<0.05,]

write.table(signiPfRBCvsuRBC_sinI7_group1.paired, file="signiPfRBCvsuRBC_C1alone_group1.paired_mod.txt", sep="\t",row.names=FALSE)
write.table(signiCS_peptide_poolsDMSO_control_sinI7_group1.paired, file="signiCS_peptide_poolvsDMSO_control_C1alone_group1.paired_mod.txt", sep="\t",row.names=FALSE)

# We will continue the analysis with this filter
save.image("SysMalVac_190_lessStringentFiltering.RData")


