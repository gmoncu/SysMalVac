############################################################################################
#################################### CPS analysis ##########################################
############################################################################################


###########################################################################
#######  Analysis of CPS immunogenicity stratified by dose group  ##########
###########################################################################

# Results of this analysis were not used in the manuscript

# This is code to perform  differential gene analysis of post-CPS immunization vs pre-immunization stratified for dose group


load("SysMalVac_190_lessStringentFiltering.RData")

pkgs<-c("affy","limma")

lapply(pkgs, require, character.only=T)	  


stimulus.f<-factor(listaTrabajoRma@phenoData$stimulus)
timepoint.f<-factor(listaTrabajoRma@phenoData$timepoint)
group.f<-factor(listaTrabajoRma@phenoData$group)
protected.f<-factor(listaTrabajoRma@phenoData$protected)
colnames(final.design)
#[1] "stimulus.fCS_peptide_pool" "stimulus.fDMSO_control"   
#[3] "stimulus.fPfRBC"           "stimulus.fuRBC"           
#[5] "timepoint.fI_7"            "group.f2"                 
#[7] "group.f3"                  "protected.fyes"           



#full design

stimulusANDtimepointANDgroupANDprotected <- factor(paste(listaTrabajoRma@phenoData$stimulus,listaTrabajoRma@phenoData$timepoint,listaTrabajoRma@phenoData$group,listaTrabajoRma@phenoData$protected,sep="."))

full.design <- model.matrix(~0+stimulusANDtimepointANDgroupANDprotected)

colnames(full.design) <- levels(stimulusANDtimepointANDgroupANDprotected)

corfit <- duplicateCorrelation(listaTrabajoRma,full.design,block=listaTrabajoRma@phenoData$donor)

#control for donor. Inter-subject correlation is input in the linear model.
fit.full.design<- lmFit(listaTrabajoRma,full.design,block=listaTrabajoRma@phenoData$donor,correlation=corfit$consensus)

#contrasts

contrast.full.design<-makeContrasts(
  #GeneList1.1
  C1vsI7inCSPeptide_DMSOremoved_3x15only = ( CS_peptide_pool.C_1.1.yes - CS_peptide_pool.I_7.1.yes ) - (DMSO_control.C_1.1.yes - DMSO_control.I_7.1.yes ),
  #GeneList1.2
  C1vsI7inCSPeptide_DMSOremoved_3x10only = ( CS_peptide_pool.C_1.2.yes - CS_peptide_pool.I_7.2.yes ) - (DMSO_control.C_1.2.yes - DMSO_control.I_7.2.yes ),
  #GeneList1.3
  C1vsI7inCSPeptide_DMSOremoved_3x5only = ( CS_peptide_pool.C_1.3.yes - CS_peptide_pool.I_7.3.yes ) - (DMSO_control.C_1.3.yes - DMSO_control.I_7.3.yes ),
  #GeneList2.1
  C1vsI7inPfRBC_uRBCremoved_3x15only = ( PfRBC.C_1.1.yes - PfRBC.I_7.1.yes ) - (uRBC.C_1.1.yes - uRBC.I_7.1.yes ),
  #GeneList2.2
  C1vsI7inPfRBC_uRBCremoved_3x10only = ( PfRBC.C_1.2.yes  - PfRBC.I_7.2.yes ) - (uRBC.C_1.2.yes - uRBC.I_7.2.yes ),
  #GeneList2.3
  C1vsI7inPfRBC_uRBCremoved_3x5only = ( PfRBC.C_1.3.yes - PfRBC.I_7.3.yes ) - (uRBC.C_1.3.yes - uRBC.I_7.3.yes ),

  levels=full.design
)


fit2.full.design <- contrasts.fit(fit.full.design, contrast.full.design)
fit2.full.design <- eBayes(fit2.full.design)

#GeneList1.1
topTable(fit2.full.design, coef="C1vsI7inCSPeptide_DMSOremoved_3x15only", adjust = "BH")
xxx<-topTable(fit2.full.design, coef="C1vsI7inCSPeptide_DMSOremoved_3x15only", adjust = "BH",sort.by="P",number=10841);signiC1vsI7inCSPeptide_DMSOremoved_3x15only<-xxx[xxx$P.Value<0.05,]
dim(signiC1vsI7inCSPeptide_DMSOremoved_3x15only)  
write.table(signiC1vsI7inCSPeptide_DMSOremoved_3x15only, file="signiC1vsI7inCSPeptide_DMSOremoved_3x15only.txt", sep="\t",row.names=FALSE)

#GeneList1.2
topTable(fit2.full.design, coef="C1vsI7inCSPeptide_DMSOremoved_3x10only", adjust = "BH")
xxx<-topTable(fit2.full.design, coef="C1vsI7inCSPeptide_DMSOremoved_3x10only", adjust = "BH",sort.by="P",number=10841);signiC1vsI7inCSPeptide_DMSOremoved_3x10only<-xxx[xxx$P.Value<0.05,]
dim(signiC1vsI7inCSPeptide_DMSOremoved_3x10only)  
write.table(signiC1vsI7inCSPeptide_DMSOremoved_3x10only, file="signiC1vsI7inCSPeptide_DMSOremoved_3x10only.txt", sep="\t",row.names=FALSE)


#GeneList1.3
topTable(fit2.full.design, coef="C1vsI7inCSPeptide_DMSOremoved_3x5only", adjust = "BH")
xxx<-topTable(fit2.full.design, coef="C1vsI7inCSPeptide_DMSOremoved_3x5only", adjust = "BH",sort.by="P",number=10841);signiC1vsI7inCSPeptide_DMSOremoved_3x5only<-xxx[xxx$P.Value<0.05,]
dim(signiC1vsI7inCSPeptide_DMSOremoved_3x5only)  
write.table(signiC1vsI7inCSPeptide_DMSOremoved_3x5only, file="signiC1vsI7inCSPeptide_DMSOremoved_3x5only.txt", sep="\t",row.names=FALSE)



#GeneList2.1
topTable(fit2.full.design, coef="C1vsI7inPfRBC_uRBCremoved_3x15only", adjust = "BH")
xxx<-topTable(fit2.full.design, coef="C1vsI7inPfRBC_uRBCremoved_3x15only", adjust = "BH",sort.by="P",number=10841);signiC1vsI7inPfRBC_uRBCremoved_3x15only<-xxx[xxx$P.Value<0.05,]
dim(signiC1vsI7inPfRBC_uRBCremoved_3x15only)  
write.table(signiC1vsI7inPfRBC_uRBCremoved_3x15only, file="signiC1vsI7inPfRBC_uRBCremoved_3x15only.txt", sep="\t",row.names=FALSE)



#GeneList2.2
topTable(fit2.full.design, coef="C1vsI7inPfRBC_uRBCremoved_3x10only", adjust = "BH")
xxx<-topTable(fit2.full.design, coef="C1vsI7inPfRBC_uRBCremoved_3x10only", adjust = "BH",sort.by="P",number=10841);signiC1vsI7inPfRBC_uRBCremoved_3x10only<-xxx[xxx$P.Value<0.05,]
dim(signiC1vsI7inPfRBC_uRBCremoved_3x10only)  
write.table(signiC1vsI7inPfRBC_uRBCremoved_3x10only, file="signiC1vsI7inPfRBC_uRBCremoved_3x10only.txt", sep="\t",row.names=FALSE)


#GeneList2.3
topTable(fit2.full.design, coef="C1vsI7inPfRBC_uRBCremoved_3x5only", adjust = "BH")
xxx<-topTable(fit2.full.design, coef="C1vsI7inPfRBC_uRBCremoved_3x5only", adjust = "BH",sort.by="P",number=10841);signiC1vsI7inPfRBC_uRBCremoved_3x5only<-xxx[xxx$P.Value<0.05,]
dim(signiC1vsI7inPfRBC_uRBCremoved_3x5only)  
write.table(signiC1vsI7inPfRBC_uRBCremoved_3x5only, file="signiC1vsI7inPfRBC_uRBCremoved_3x5only.txt", sep="\t",row.names=FALSE)



save.image("SysMalVac_190_lessStringentFiltering_finalAnalysis.RData")
