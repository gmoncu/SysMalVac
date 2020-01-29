##############################################################################################################
#################################### RTS,S analysis (MAIN ANALYSIS) ##########################################
##############################################################################################################

# This is code to perform the main microarray data analysis of the RTS,S study
# It includes normalization using RMA method, filtering of background noise and invariant genes and differential gene expression
# Datasets of differential gene expression provided in the manuscript are based on this analysis
# We used the ranked genes obtained here for the gene set enrichment analyses based on differential gene expression 


# load functions

sapply(list.files(path="utils/Filter",pattern="[.]R$", full.names=TRUE), source)

sapply(list.files(path="utils/QualityControl",pattern="[.]R$", full.names=TRUE), source)

sapply(list.files(path="utils/Plot",pattern="[.]R$", full.names=TRUE), source)

sapply(list.files(path="utils/statistic",pattern="[.]R$", full.names=TRUE), source)


# Packages 

pkgs<-c("ggplot2","GGally","CCA","oligo","biomaRt","topGO","multtest","Rgraphivz","hugene21sttranscriptcluster.db","annotate",
"GO.db""statmod","oligoData","gplots","arrayQualityMetrics","hugene21stcdf","genefilter","affy","limma","marray","amap")

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

QualityControl(path=getwd(),method="RMA",control="Both",intgroup="casecon")

#Filter

#filter control probes and low intensity

exprsdataCEL<-intensityFilter(path=getwd(),kvalue=2)

#18801 genes

dfataCEL<-exprs(exprsdataCEL)

# intensity filter log2(100) and IQR. 

exprsdataNorm<-normFilter(path=getwd(),
		      method="RMA",
		      fdataCEL=exprsdataCEL)
		      		      

exprslistafinal_05<-noninformativeFilter(data=exprsdataNorm,
				  typefilter="IQR",
				  cutoff=0.5)
#7257

exprslistafinal_04<-noninformativeFilter(data=exprsdataNorm,
				  typefilter="IQR",
				  cutoff=0.4)
#10567

exprslistafinal_03<-noninformativeFilter(data=exprsdataNorm,
				  typefilter="IQR",
				  cutoff=0.3)
#14840


# Use the filter IQR<0.4
listafinal_04<-exprs(exprslistafinal_04)	

namesmatrix <- colnames(listafinal_04)

for (i in c(1:length(namesmatrix))){
    j<-0
    while (substr(namesmatrix[i],j,j)!=".")
      {
      j=j+1
      }
    namesmatrix[i]<-substr(namesmatrix[i],1,(j-1))
  }

colnames(listafinal_04)<-namesmatrix


myPData<-read.table("phenoDATA.txt", sep="",header=TRUE,as.is=TRUE)

sampleout<- c("13SE3831", "13SE3906", "13SE3914", "13SE3926","13SE3992","13SE3993","13SE4022","13SE4038", "13SE4046", "13SE4054", "13SE4077", "13SE4085",
	      "13SE4093", "13SE4185", "13SE4188", "13SE4189","13SE4196", "13SE4205","13SE4211","13SE4431", "13SE4432", "13SE4433", "13SE4434", "13SE4498",
	      "13SE4499", "13SE4500", "13SE4502", "13SE4503","13SE4504", "13SE4508","13SE4510")


for (i in c(1: length(sampleout))){

  myPData<-myPData[-which(myPData$PGK_code==sampleout[i]),]

}


myPData2<-myPData[,c(4,8,9,10)]

colnames(myPData2)<-c("Stimulation","timepoint","casecon","groupnb")

rownames(myPData2)<-myPData$PGK_code

myPData2$groupnb[which(myPData2$groupnb=="2")]<-"1" # 1 is RTS,S
myPData2$groupnb[which(myPData2$groupnb=="3")]<-"2" # 2 is comparator
myPData2$casecon[which(myPData2$casecon=="non-malaria")]<-"non_malaria" 

pData(exprslistafinal_04)<-myPData2


###########################################################################
####### Analysis of protection in RTS,S and comparator vaccinees ##########
###########################################################################

#limma

experimento <- factor(paste(exprslistafinal_04@phenoData$casecon,
			    exprslistafinal_04@phenoData$Stimulation,
			    exprslistafinal_04@phenoData$timepoint,
			    exprslistafinal_04@phenoData$groupnb,sep=""))

design <- model.matrix(~ 0 + experimento)

colnames(design) <- c(levels(experimento))

fit<-lmFit(exprslistafinal_04,design)

#GO annotation 

GO2geneID<-as.list(hugene21sttranscriptclusterGO2ALLPROBES)

geneID2GO<-inverseList(GO2geneID)

topDiffGenes<-function (allScore) { 
    return(allScore < 0.05)
}

save(list=ls(all=TRUE),file="GO_hugene21.RData")


system("mkdir -p Results/malaria_CSP_RTS") # Protection in RTS,S vaccinees at M3 (post-vaccination) for CSP stimulation
system("mkdir -p Results/malaria_CSP_CONTROL") # Protection in Comparator vaccinees at M3 (post-vaccination) for CSP stimulation
system("mkdir -p Results/malaria_iRBC_RTS") # Protection in RTS,S vaccinees at M3 (post-vaccination) for PfRBC stimulation
system("mkdir -p Results/malaria_iRBC_CONTROL") # Protection in Comparator vaccinees at M3 (post-vaccination) for PfRBC stimulation

# Contrasts
    ### Note that these contrasts result in the comparison of malaria vs. non-malaria (non-protected vs. non-protected) and we reported in the manuscript the opposite sense: protected vs. non-protected. 
    ### We changed to results so they would be consistent with the CPS analysis.

contrast.matrix <- makeContrasts(
  
  CSP_RTSvsDMSO_RTS=((malariaCSPM31-malariaDMSOM31)-(non_malariaCSPM31-non_malariaDMSOM31)), # Protection in RTS,S vaccinees at M3 (post-vaccination) for CSP stimulation
  
  CSP_CONTROLvsDMSO_CONTROL=((malariaCSPM32-malariaDMSOM32)-(non_malariaCSPM32-non_malariaDMSOM32)), # Protection in Comparator vaccinees at M3 (post-vaccination) for CSP stimulation
  
  iRBC_RTSvsuRBC_RTS=((malariaiRBCM31-malariauRBCM31)-(non_malariaiRBCM31-non_malariauRBCM31)), # Protection in RTS,S vaccinees at M3 (post-vaccination) for PfRBC stimulation
  
  iRBC_CONTROLvsuRBC_CONTROL=((malariaiRBCM32-malariauRBCM32)-(non_malariaiRBCM32-non_malariauRBCM32)),  # Protection in Comparator vaccinees at M3 (post-vaccination) for PfRBC stimulation
  
  levels=design)
  
fit2<-eBayes(contrasts.fit(fit,contrast.matrix))


save(contrast.matrix,fit,design,myPData2,fit2,exprslistafinal_04,listafinal_04,file="Results/SysMalVac_December.RData")

########################################### CSP_RTSvsDMSO_RTS (malaria_CSP_RTS) ####################################################################################

load("Results/SysMalVac_December.RData")

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="CSP_RTSvsDMSO_RTS",sort.by="P",number=dim(exprslistafinal_04)[1])

CSP_RTSvsDMSO_RTS_bonferroni_pval005<-genes_bonferroni[genes_bonferroni$P.Value<0.05,]

CSP_RTSvsDMSO_RTS_bonferroni_adjpval005<-genes_bonferroni[genes_bonferroni$adj.P.Val<0.05,]

fileout<-"CSP_RTSvsDMSO_RTS_bonferroni_adjpval005"

write.table(CSP_RTSvsDMSO_RTS_bonferroni_adjpval005[,c(1,2,5,6)],file=paste(paste("Results/malaria_CSP_RTS",fileout,sep="/"),".txt",sep=""))


fileout<-"CSP_RTSvsDMSO_RTS_bonferroni_pval005"

write.table(CSP_RTSvsDMSO_RTS_bonferroni_pval005[,c(1,2,4,5,6)],file=paste(paste("Results/malaria_CSP_RTS",fileout,sep="/"),".txt",sep=""))


adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="CSP_RTSvsDMSO_RTS",sort.by="P",number=dim(exprslistafinal_04)[1])

CSP_RTSvsDMSO_RTS_fdr_pval005<-genes_fdr[genes_fdr$P.Value<0.05,]

CSP_RTSvsDMSO_RTS_fdr_adjpval005<-genes_fdr[genes_fdr$adj.P.Val<0.05,]

fileout<-"CSP_RTSvsDMSO_RTS_fdr_adjpval005"

write.table(CSP_RTSvsDMSO_RTS_fdr_adjpval005[,c(1,2,5,6)],file=paste(paste("Results/malaria_CSP_RTS",fileout,sep="/"),".txt",sep=""))


fileout<-"CSP_RTSvsDMSO_RTS_fdr_pval005"

write.table(CSP_RTSvsDMSO_RTS_fdr_pval005[,c(1,2,4,5,6)],file=paste(paste("Results/malaria_CSP_RTS",fileout,sep="/"),".txt",sep=""))

genes<-genes_fdr

rownames(genes_fdr)<-genes_fdr$ID

namesDMSOnon<-which(myPData2$Stimulation=="DMSO" & myPData2$casecon=="non_malaria" & myPData2$groupnb=="1" & myPData2$timepoint=="M3")

namesDMSOmalaria<-which(myPData2$Stimulation=="DMSO" & myPData2$casecon=="malaria" & myPData2$groupnb=="1" & myPData2$timepoint=="M3")

namesCSPnon<-which(myPData2$Stimulation=="CSP" & myPData2$casecon=="non_malaria" & myPData2$groupnb=="1" & myPData2$timepoint=="M3")

namesCSPmalaria<-which(myPData2$Stimulation=="CSP" & myPData2$casecon=="malaria" & myPData2$groupnb=="1" & myPData2$timepoint=="M3")

names<- colnames(listafinal_04)[c(namesCSPnon,namesCSPmalaria,namesDMSOnon,namesDMSOmalaria)]

rownames(CSP_RTSvsDMSO_RTS_fdr_pval005)<-CSP_RTSvsDMSO_RTS_fdr_pval005$ID

png(paste(paste(paste(paste("Results/malaria_CSP_RTS","heatmap_CSP_RTSvsDMSO_RTS",sep="/"),adjust,sep="_"),pvalue,sep="_"),".png",sep=""),width = 1200, height = 1200)

  heatmap(listafinal_04[rownames(CSP_RTSvsDMSO_RTS_fdr_pval005),names],distfun=function(c)
  {Dist(c,method="correlation")},hclustfun=function(d)
  {hclust(d,method="average")},Colv=cov(listafinal_04[rownames(CSP_RTSvsDMSO_RTS_fdr_pval005),names]),Rowv=cov(listafinal_04[rownames(CSP_RTSvsDMSO_RTS_fdr_adjpval005),names]))

dev.off()

save(list=ls(all=TRUE),file="Results/malaria_CSP_RTS/CSP_RTSvsDMSO_RTS.RData")

#Volcano Plot

load("Results/malaria_CSP_RTS/CSP_RTSvsDMSO_RTS.RData")

gene_list<-genes_fdr

gene_list$threshold = as.factor( gene_list$P.Value < 0.05)

gene_list$threshold_fdr = as.factor( gene_list$adj.P.Val < 0.05)

gene_list$colour[which(gene_list$threshold=="FALSE")]<-"black"

gene_list$colour[which(gene_list$threshold=="TRUE" & gene_list$threshold_fdr=="FALSE")]<-"grey"

gene_list$colour[which(gene_list$threshold=="TRUE" & gene_list$threshold_fdr=="TRUE")]<-"red"

png("Results/malaria_CSP_RTS/VolcanoPlot_CSP_RTSvsDMSO_RTS.png")
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(P.Value))) +
    geom_point(alpha=1,pch=20, size=3,col=gene_list$colour) +
    theme(legend.position = "none") +
    xlim(c(-10, 10)) + ylim(c(0, 15)) +
    xlab("log2 fold change") + ylab("log p-value")
    main="CSP_RTSvsDMSO_RTS"
g
dev.off()

#
load("GO_hugene21.RData")

geneList<-genes$P.Value

names(geneList)<-rownames(genes)

GOdataBP <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=topDiffGenes ,annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOdataMF <- new("topGOdata", ontology = "MF", allGenes = geneList,geneSel=topDiffGenes ,annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOdataCC <- new("topGOdata", ontology = "CC", allGenes = geneList,geneSel=topDiffGenes ,annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Enrichment Analysis

resultFisherBP<-runTest(GOdataBP,algorithm="classic",statistic="fisher")

resultFisherMF<-runTest(GOdataMF,algorithm="classic",statistic="fisher")

resultFisherCC<-runTest(GOdataCC,algorithm="classic",statistic="fisher")

allResBP<-GenTable(GOdataBP,
		 pValue_Fisher=resultFisherBP,
		 orderBy="pValueFisher",
		 ranksOf="pValueFisher",
		 numChar=1000,
		 topNodes=sum(resultFisherBP@score<0.05))
		 
allResMF<-GenTable(GOdataMF,
		 pValue_Fisher=resultFisherMF,
		 orderBy="pValueFisher",
		 ranksOf="pValueFisher",
		 numChar=1000,
		 topNodes=sum(resultFisherMF@score<0.05))
		 
allResCC<-GenTable(GOdataCC,
		 pValue_Fisher=resultFisherCC,
		 orderBy="pValueFisher",
		 ranksOf="pValueFisher",
		 numChar=1000,
		 topNodes=sum(resultFisherCC@score<0.05))

allResBP$p.adj_fdr<-round(p.adjust(allResBP$pValue_Fisher, method = "fdr", n =length(GOdataBP@graph@nodes)),4)

allResMF$p.adj_fdr<-round(p.adjust(allResMF$pValue_Fisher, method = "fdr", n =length(GOdataMF@graph@nodes)),4)

allResCC$p.adj_fdr<-round(p.adjust(allResCC$pValue_Fisher, method = "fdr", n =length(GOdataCC@graph@nodes)),4)

allResBP$p.adj_bonferroni<-round(p.adjust(allResBP$pValue_Fisher, method = "bonferroni", n = length(GOdataBP@graph@nodes)),4)

allResMF$p.adj_bonferroni<-round(p.adjust(allResMF$pValue_Fisher, method = "bonferroni", n = length(GOdataMF@graph@nodes)),4)

allResCC$p.adj_bonferroni<-round(p.adjust(allResCC$pValue_Fisher, method = "bonferroni", n = length(GOdataCC@graph@nodes)),4)


goIDBP<-allResBP[,"GO.ID"]

goIDMF<-allResMF[,"GO.ID"]

goIDCC<-allResCC[,"GO.ID"]



gtBP<-printGenes(GOdataBP,whichTerms=goIDBP,chip="hugene21sttranscriptcluster.db",geneCutOff=max(allResBP$Significant),pvalCutOff=0.05)

gtMF<-printGenes(GOdataMF,whichTerms=goIDMF,chip="hugene21sttranscriptcluster.db",geneCutOff=max(allResMF$Significant),pvalCutOff=0.05)

gtCC<-printGenes(GOdataCC,whichTerms=goIDCC,chip="hugene21sttranscriptcluster.db",geneCutOff=max(allResCC$Significant),pvalCutOff=0.05)


for (i in c(1:length(goIDBP))){
  
  allResBP$gene_SYMBOL[[i]]<-paste(as.vector(gtBP[[i]]$Symbol.id), collapse=",")
  
  }

for (i in c(1:length(goIDMF))){
  
  allResMF$gene_SYMBOL[[i]]<-paste(as.vector(gtMF[[i]]$Symbol.id), collapse=",")
  
  }
  
for (i in c(1:length(goIDCC))){
  
  allResCC$gene_SYMBOL[[i]]<-paste(as.vector(gtCC[[i]]$Symbol.id), collapse=",")
  
  }
  
  
write.csv(allResCC,file="Results/malaria_CSP_RTS/GO_CC_CSP_RTSvsDMSO_RTS.csv")
  
write.csv(allResMF,file="Results/malaria_CSP_RTS/GO_MF_CSP_RTSvsDMSO_RTS.csv")
  
write.csv(allResBP,file="Results/malaria_CSP_RTS/GO_BP_CSP_RTSvsDMSO_RTS.csv")

save(list=ls(all=TRUE),file="Results/malaria_CSP_RTS/GO_CSP_RTSvsDMSO_RTS.RData")



######################################## CSP_CONTROLvsDMSO_CONTROL (malaria_CSP_CONTROL) ###########################################################################################################

load("Results/SysMalVac_December.RData")

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="CSP_CONTROLvsDMSO_CONTROL",sort.by="P",number=dim(exprslistafinal_04)[1])

CSP_CONTROLvsDMSO_CONTROL_bonferroni_pval005<-genes_bonferroni[genes_bonferroni$P.Value<0.05,]

CSP_CONTROLvsDMSO_CONTROL_bonferroni_adjpval005<-genes_bonferroni[genes_bonferroni$adj.P.Val<0.05,]

fileout<-"CSP_CONTROLvsDMSO_CONTROL_bonferroni_adjpval005"

write.table(CSP_CONTROLvsDMSO_CONTROL_bonferroni_adjpval005[,c(1,2,5,6)],file=paste(paste("Results/malaria_CSP_CONTROL",fileout,sep="/"),".txt",sep=""))

fileout<-"CSP_CONTROLvsDMSO_CONTROL_bonferroni_pval005"

write.table(CSP_CONTROLvsDMSO_CONTROL_bonferroni_pval005[,c(1,2,4,5,6)],file=paste(paste("Results/malaria_CSP_CONTROL",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="CSP_CONTROLvsDMSO_CONTROL",sort.by="P",number=dim(exprslistafinal_04)[1])

CSP_CONTROLvsDMSO_CONTROL_fdr_pval005<-genes_fdr[genes_fdr$P.Value<0.05,]

CSP_CONTROLvsDMSO_CONTROL_fdr_adjpval005<-genes_fdr[genes_fdr$adj.P.Val<0.05,]

fileout<-"CSP_CONTROLvsDMSO_CONTROL_fdr_adjpval005"

write.table(CSP_CONTROLvsDMSO_CONTROL_fdr_adjpval005[,c(1,2,5,6)],file=paste(paste("Results/malaria_CSP_CONTROL",fileout,sep="/"),".txt",sep=""))


fileout<-"CSP_CONTROLvsDMSO_CONTROL_fdr_pval005"

write.table(CSP_CONTROLvsDMSO_CONTROL_fdr_pval005[,c(1,2,4,5,6)],file=paste(paste("Results/malaria_CSP_CONTROL",fileout,sep="/"),".txt",sep=""))

genes<-genes_fdr

rownames(genes_fdr)<-genes_fdr$ID

namesDMSOnon<-which(myPData2$Stimulation=="DMSO" & myPData2$casecon=="non_malaria" & myPData2$groupnb=="2" & myPData2$timepoint=="M3")

namesDMSOmalaria<-which(myPData2$Stimulation=="DMSO" & myPData2$casecon=="malaria" & myPData2$groupnb=="2" & myPData2$timepoint=="M3")

namesCSPnon<-which(myPData2$Stimulation=="CSP" & myPData2$casecon=="non_malaria" & myPData2$groupnb=="2" & myPData2$timepoint=="M3")

namesCSPmalaria<-which(myPData2$Stimulation=="CSP" & myPData2$casecon=="malaria" & myPData2$groupnb=="2" & myPData2$timepoint=="M3")

names<- colnames(listafinal_04)[c(namesCSPnon,namesCSPmalaria,namesDMSOnon,namesDMSOmalaria)]

rownames(CSP_CONTROLvsDMSO_CONTROL_fdr_pval005)<-CSP_CONTROLvsDMSO_CONTROL_fdr_pval005$ID

png(paste(paste(paste(paste("Results/malaria_CSP_CONTROL","heatmap_CSP_CONTROLvsDMSO_CONTROL",sep="/"),adjust,sep="_"),pvalue,sep="_"),".png",sep=""),width = 1200, height = 1200)

  heatmap(listafinal_04[rownames(CSP_CONTROLvsDMSO_CONTROL_fdr_pval005),names],distfun=function(c)
  {Dist(c,method="correlation")},hclustfun=function(d)
  {hclust(d,method="average")},Colv=cov(listafinal_04[rownames(CSP_CONTROLvsDMSO_CONTROL_fdr_pval005),names]),Rowv=cov(listafinal_04[rownames(CSP_CONTROLvsDMSO_CONTROL_fdr_adjpval005),names]))

dev.off()

save(list=ls(all=TRUE),file="Results/malaria_CSP_CONTROL/CSP_CONTROLvsDMSO_CONTROL.RData")

#Volcano Plot

load("Results/malaria_CSP_CONTROL/CSP_CONTROLvsDMSO_CONTROL.RData")

gene_list<-genes_fdr

gene_list$threshold = as.factor( gene_list$P.Value < 0.05)

gene_list$threshold_fdr = as.factor( gene_list$adj.P.Val < 0.05)

gene_list$colour[which(gene_list$threshold=="FALSE")]<-"black"

gene_list$colour[which(gene_list$threshold=="TRUE" & gene_list$threshold_fdr=="FALSE")]<-"grey"

gene_list$colour[which(gene_list$threshold=="TRUE" & gene_list$threshold_fdr=="TRUE")]<-"red"

png("Results/malaria_CSP_CONTROL/VolcanoPlot_CSP_CONTROLvsDMSO_CONTROL.png")
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(P.Value))) +
    geom_point(alpha=1,pch=20, size=3,col=gene_list$colour) +
    theme(legend.position = "none") +
    xlim(c(-10, 10)) + ylim(c(0, 15)) +
    xlab("log2 fold change") + ylab("log p-value")
    main="CSP_CONTROLvsDMSO_CONTROL"
g
dev.off()

# GO annotation 

load("GO_hugene21.RData")

geneList<-genes$P.Value

names(geneList)<-rownames(genes)

GOdataBP <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=topDiffGenes ,annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOdataMF <- new("topGOdata", ontology = "MF", allGenes = geneList,geneSel=topDiffGenes ,annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOdataCC <- new("topGOdata", ontology = "CC", allGenes = geneList,geneSel=topDiffGenes ,annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Enrichment Analysis

resultFisherBP<-runTest(GOdataBP,algorithm="classic",statistic="fisher")

resultFisherMF<-runTest(GOdataMF,algorithm="classic",statistic="fisher")

resultFisherCC<-runTest(GOdataCC,algorithm="classic",statistic="fisher")

allResBP<-GenTable(GOdataBP,
		 pValue_Fisher=resultFisherBP,
		 orderBy="pValueFisher",
		 ranksOf="pValueFisher",
		 numChar=1000,
		 topNodes=sum(resultFisherBP@score<0.05))
		 
allResMF<-GenTable(GOdataMF,
		 pValue_Fisher=resultFisherMF,
		 orderBy="pValueFisher",
		 ranksOf="pValueFisher",
		 numChar=1000,
		 topNodes=sum(resultFisherMF@score<0.05))
		 
allResCC<-GenTable(GOdataCC,
		 pValue_Fisher=resultFisherCC,
		 orderBy="pValueFisher",
		 ranksOf="pValueFisher",
		 numChar=1000,
		 topNodes=sum(resultFisherCC@score<0.05))

allResBP$p.adj_fdr<-round(p.adjust(allResBP$pValue_Fisher, method = "fdr", n =length(GOdataBP@graph@nodes)),4)

allResMF$p.adj_fdr<-round(p.adjust(allResMF$pValue_Fisher, method = "fdr", n =length(GOdataMF@graph@nodes)),4)

allResCC$p.adj_fdr<-round(p.adjust(allResCC$pValue_Fisher, method = "fdr", n =length(GOdataCC@graph@nodes)),4)

allResBP$p.adj_bonferroni<-round(p.adjust(allResBP$pValue_Fisher, method = "bonferroni", n = length(GOdataBP@graph@nodes)),4)

allResMF$p.adj_bonferroni<-round(p.adjust(allResMF$pValue_Fisher, method = "bonferroni", n = length(GOdataMF@graph@nodes)),4)

allResCC$p.adj_bonferroni<-round(p.adjust(allResCC$pValue_Fisher, method = "bonferroni", n = length(GOdataCC@graph@nodes)),4)


goIDBP<-allResBP[,"GO.ID"]

goIDMF<-allResMF[,"GO.ID"]

goIDCC<-allResCC[,"GO.ID"]



gtBP<-printGenes(GOdataBP,whichTerms=goIDBP,chip="hugene21sttranscriptcluster.db",geneCutOff=max(allResBP$Significant),pvalCutOff=0.05)

gtMF<-printGenes(GOdataMF,whichTerms=goIDMF,chip="hugene21sttranscriptcluster.db",geneCutOff=max(allResMF$Significant),pvalCutOff=0.05)

gtCC<-printGenes(GOdataCC,whichTerms=goIDCC,chip="hugene21sttranscriptcluster.db",geneCutOff=max(allResCC$Significant),pvalCutOff=0.05)


for (i in c(1:length(goIDBP))){
  
  allResBP$gene_SYMBOL[[i]]<-paste(as.vector(gtBP[[i]]$Symbol.id), collapse=",")
  
  }

for (i in c(1:length(goIDMF))){
  
  allResMF$gene_SYMBOL[[i]]<-paste(as.vector(gtMF[[i]]$Symbol.id), collapse=",")
  
  }
  
for (i in c(1:length(goIDCC))){
  
  allResCC$gene_SYMBOL[[i]]<-paste(as.vector(gtCC[[i]]$Symbol.id), collapse=",")
  
  }
  
  
write.csv(allResCC,file="Results/malaria_CSP_CONTROL/GO_CC_CSP_CONTROLvsDMSO_CONTROL.csv")
  
write.csv(allResMF,file="Results/malaria_CSP_CONTROL/GO_MF_CSP_CONTROLvsDMSO_CONTROL.csv")
  
write.csv(allResBP,file="Results/malaria_CSP_CONTROL/GO_BP_CSP_CONTROLvsDMSO_CONTROL.csv")

save(list=ls(all=TRUE),file="Results/malaria_CSP_CONTROL/GO_CSP_CONTROLvsDMSO_CONTROL.RData")


################################################## iRBC_RTSvsuRBC_RTS (malaria_iRBC_RTS) #########################################################################

load("Results/SysMalVac_December.RData")

adjust="bonferroni";pvalue=0.05

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="iRBC_RTSvsuRBC_RTS",sort.by="P",number=dim(exprslistafinal_04)[1])

iRBC_RTSvsuRBC_RTS_bonferroni_pval005<-genes_bonferroni[genes_bonferroni$P.Value<0.05,]

iRBC_RTSvsuRBC_RTS_bonferroni_adjpval005<-genes_bonferroni[genes_bonferroni$adj.P.Val<0.05,]

fileout<-"iRBC_RTSvsuRBC_RTS_bonferroni_adjpval005"

write.table(iRBC_RTSvsuRBC_RTS_bonferroni_adjpval005[,c(1,2,5,6)],file=paste(paste("Results/malaria_iRBC_RTS",fileout,sep="/"),".txt",sep=""))


fileout<-"iRBC_RTSvsuRBC_RTS_bonferroni_pval005"

write.table(iRBC_RTSvsuRBC_RTS_bonferroni_pval005[,c(1,2,4,5,6)],file=paste(paste("Results/malaria_iRBC_RTS",fileout,sep="/"),".txt",sep=""))


adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="iRBC_RTSvsuRBC_RTS",sort.by="P",number=dim(exprslistafinal_04)[1])

iRBC_RTSvsuRBC_RTS_fdr_pval005<-genes_fdr[genes_fdr$P.Value<0.05,]

iRBC_RTSvsuRBC_RTS_fdr_adjpval005<-genes_fdr[genes_fdr$adj.P.Val<0.05,]

fileout<-"iRBC_RTSvsuRBC_RTS_fdr_adjpval005"

write.table(iRBC_RTSvsuRBC_RTS_fdr_adjpval005[,c(1,2,5,6)],file=paste(paste("Results/malaria_iRBC_RTS",fileout,sep="/"),".txt",sep=""))


fileout<-"iRBC_RTSvsuRBC_RTS_fdr_pval005"

write.table(iRBC_RTSvsuRBC_RTS_fdr_pval005[,c(1,2,4,5,6)],file=paste(paste("Results/malaria_iRBC_RTS",fileout,sep="/"),".txt",sep=""))

genes<-genes_fdr

rownames(genes_fdr)<-genes_fdr$ID

namesuRBCnon<-which(myPData2$Stimulation=="uRBC" & myPData2$casecon=="non_malaria" & myPData2$groupnb=="1" & myPData2$timepoint=="M3")

namesuRBCmalaria<-which(myPData2$Stimulation=="uRBC" & myPData2$casecon=="malaria" & myPData2$groupnb=="1" & myPData2$timepoint=="M3")

namesiRBCnon<-which(myPData2$Stimulation=="iRBC" & myPData2$casecon=="non_malaria" & myPData2$groupnb=="1" & myPData2$timepoint=="M3")

namesiRBCmalaria<-which(myPData2$Stimulation=="iRBC" & myPData2$casecon=="malaria" & myPData2$groupnb=="1" & myPData2$timepoint=="M3")

names<- colnames(listafinal_04)[c(namesiRBCnon,namesiRBCmalaria,namesuRBCnon,namesuRBCmalaria)]

rownames(iRBC_RTSvsuRBC_RTS_fdr_pval005)<-iRBC_RTSvsuRBC_RTS_fdr_pval005$ID

png(paste(paste(paste(paste("Results/malaria_iRBC_RTS","heatmap_iRBC_RTSvsuRBC_RTS",sep="/"),adjust,sep="_"),pvalue,sep="_"),".png",sep=""),width = 1200, height = 1200)

  heatmap(listafinal_04[rownames(iRBC_RTSvsuRBC_RTS_fdr_pval005),names],distfun=function(c)
  {Dist(c,method="correlation")},hclustfun=function(d)
  {hclust(d,method="average")},Colv=cov(listafinal_04[rownames(iRBC_RTSvsuRBC_RTS_fdr_pval005),names]),Rowv=cov(listafinal_04[rownames(iRBC_RTSvsuRBC_RTS_fdr_adjpval005),names]))

dev.off()

save(list=ls(all=TRUE),file="Results/malaria_iRBC_RTS/iRBC_RTSvsuRBC_RTS.RData")

#Volcano Plot

load("Results/malaria_iRBC_RTS/iRBC_RTSvsuRBC_RTS.RData")

gene_list<-genes_fdr

gene_list$threshold = as.factor( gene_list$P.Value < 0.05)

gene_list$threshold_fdr = as.factor( gene_list$adj.P.Val < 0.05)

gene_list$colour[which(gene_list$threshold=="FALSE")]<-"black"

gene_list$colour[which(gene_list$threshold=="TRUE" & gene_list$threshold_fdr=="FALSE")]<-"grey"

gene_list$colour[which(gene_list$threshold=="TRUE" & gene_list$threshold_fdr=="TRUE")]<-"red"

png("Results/malaria_iRBC_RTS/VolcanoPlot_iRBC_RTSvsuRBC_RTS.png")
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(P.Value))) +
    geom_point(alpha=1,pch=20, size=3,col=gene_list$colour) +
    theme(legend.position = "none") +
    xlim(c(-10, 10)) + ylim(c(0, 15)) +
    xlab("log2 fold change") + ylab("log p-value")
    main="iRBC_RTSvsuRBC_RTS"
g
dev.off()

# GO annotation 

load("GO_hugene21.RData")

geneList<-genes$P.Value

names(geneList)<-rownames(genes)

GOdataBP <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=topDiffGenes ,annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOdataMF <- new("topGOdata", ontology = "MF", allGenes = geneList,geneSel=topDiffGenes ,annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOdataCC <- new("topGOdata", ontology = "CC", allGenes = geneList,geneSel=topDiffGenes ,annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Enrichment Analysis

resultFisherBP<-runTest(GOdataBP,algorithm="classic",statistic="fisher")

resultFisherMF<-runTest(GOdataMF,algorithm="classic",statistic="fisher")

resultFisherCC<-runTest(GOdataCC,algorithm="classic",statistic="fisher")

allResBP<-GenTable(GOdataBP,
		 pValue_Fisher=resultFisherBP,
		 orderBy="pValueFisher",
		 ranksOf="pValueFisher",
		 numChar=1000,
		 topNodes=sum(resultFisherBP@score<0.05))
		 
allResMF<-GenTable(GOdataMF,
		 pValue_Fisher=resultFisherMF,
		 orderBy="pValueFisher",
		 ranksOf="pValueFisher",
		 numChar=1000,
		 topNodes=sum(resultFisherMF@score<0.05))
		 
allResCC<-GenTable(GOdataCC,
		 pValue_Fisher=resultFisherCC,
		 orderBy="pValueFisher",
		 ranksOf="pValueFisher",
		 numChar=1000,
		 topNodes=sum(resultFisherCC@score<0.05))

allResBP$p.adj_fdr<-round(p.adjust(allResBP$pValue_Fisher, method = "fdr", n =length(GOdataBP@graph@nodes)),4)

allResMF$p.adj_fdr<-round(p.adjust(allResMF$pValue_Fisher, method = "fdr", n =length(GOdataMF@graph@nodes)),4)

allResCC$p.adj_fdr<-round(p.adjust(allResCC$pValue_Fisher, method = "fdr", n =length(GOdataCC@graph@nodes)),4)

allResBP$p.adj_bonferroni<-round(p.adjust(allResBP$pValue_Fisher, method = "bonferroni", n = length(GOdataBP@graph@nodes)),4)

allResMF$p.adj_bonferroni<-round(p.adjust(allResMF$pValue_Fisher, method = "bonferroni", n = length(GOdataMF@graph@nodes)),4)

allResCC$p.adj_bonferroni<-round(p.adjust(allResCC$pValue_Fisher, method = "bonferroni", n = length(GOdataCC@graph@nodes)),4)


goIDBP<-allResBP[,"GO.ID"]

goIDMF<-allResMF[,"GO.ID"]

goIDCC<-allResCC[,"GO.ID"]



gtBP<-printGenes(GOdataBP,whichTerms=goIDBP,chip="hugene21sttranscriptcluster.db",geneCutOff=max(allResBP$Significant),pvalCutOff=0.05)

gtMF<-printGenes(GOdataMF,whichTerms=goIDMF,chip="hugene21sttranscriptcluster.db",geneCutOff=max(allResMF$Significant),pvalCutOff=0.05)

gtCC<-printGenes(GOdataCC,whichTerms=goIDCC,chip="hugene21sttranscriptcluster.db",geneCutOff=max(allResCC$Significant),pvalCutOff=0.05)


for (i in c(1:length(goIDBP))){
  
  allResBP$gene_SYMBOL[[i]]<-paste(as.vector(gtBP[[i]]$Symbol.id), collapse=",")
  
  }

for (i in c(1:length(goIDMF))){
  
  allResMF$gene_SYMBOL[[i]]<-paste(as.vector(gtMF[[i]]$Symbol.id), collapse=",")
  
  }
  
for (i in c(1:length(goIDCC))){
  
  allResCC$gene_SYMBOL[[i]]<-paste(as.vector(gtCC[[i]]$Symbol.id), collapse=",")
  
  }
  
  
write.csv(allResCC,file="Results/malaria_iRBC_RTS/GO_CC_iRBC_RTSvsuRBC_RTS.csv")
  
write.csv(allResMF,file="Results/malaria_iRBC_RTS/GO_MF_iRBC_RTSvsuRBC_RTS.csv")
  
write.csv(allResBP,file="Results/malaria_iRBC_RTS/GO_BP_iRBC_RTSvsuRBC_RTS.csv")

save(list=ls(all=TRUE),file="Results/malaria_iRBC_RTS/GO_iRBC_RTSvsuRBC_RTS.RData")


################################### iRBC_CONTROLvsuRBC_CONTROL (malaria_iRBC_CONTROL) ##########################################################################################

load("Results/SysMalVac_December.RData")

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="iRBC_CONTROLvsuRBC_CONTROL",sort.by="P",number=dim(exprslistafinal_04)[1])

iRBC_CONTROLvsuRBC_CONTROL_bonferroni_pval005<-genes_bonferroni[genes_bonferroni$P.Value<0.05,]

iRBC_CONTROLvsuRBC_CONTROL_bonferroni_adjpval005<-genes_bonferroni[genes_bonferroni$adj.P.Val<0.05,]

fileout<-"iRBC_CONTROLvsuRBC_CONTROL_bonferroni_adjpval005"

write.table(iRBC_CONTROLvsuRBC_CONTROL_bonferroni_adjpval005[,c(1,2,5,6)],file=paste(paste("Results/malaria_iRBC_CONTROL",fileout,sep="/"),".txt",sep=""))

fileout<-"iRBC_CONTROLvsuRBC_CONTROL_bonferroni_pval005"

write.table(iRBC_CONTROLvsuRBC_CONTROL_bonferroni_pval005[,c(1,2,4,5,6)],file=paste(paste("Results/malaria_iRBC_CONTROL",fileout,sep="/"),".txt",sep=""))

adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="iRBC_CONTROLvsuRBC_CONTROL",sort.by="P",number=dim(exprslistafinal_04)[1])

iRBC_CONTROLvsuRBC_CONTROL_fdr_pval005<-genes_fdr[genes_fdr$P.Value<0.05,]

iRBC_CONTROLvsuRBC_CONTROL_fdr_adjpval005<-genes_fdr[genes_fdr$adj.P.Val<0.05,]

fileout<-"iRBC_CONTROLvsuRBC_CONTROL_fdr_adjpval005"

write.table(iRBC_CONTROLvsuRBC_CONTROL_fdr_adjpval005[,c(1,2,5,6)],file=paste(paste("Results/malaria_iRBC_CONTROL",fileout,sep="/"),".txt",sep=""))


fileout<-"iRBC_CONTROLvsuRBC_CONTROL_fdr_pval005"

write.table(iRBC_CONTROLvsuRBC_CONTROL_fdr_pval005[,c(1,2,4,5,6)],file=paste(paste("Results/malaria_iRBC_CONTROL",fileout,sep="/"),".txt",sep=""))

genes<-genes_fdr

rownames(genes_fdr)<-genes_fdr$ID

namesuRBCnon<-which(myPData2$Stimulation=="uRBC" & myPData2$casecon=="non_malaria" & myPData2$groupnb=="2" & myPData2$timepoint=="M3")

namesuRBCmalaria<-which(myPData2$Stimulation=="uRBC" & myPData2$casecon=="malaria" & myPData2$groupnb=="2" & myPData2$timepoint=="M3")

namesiRBCnon<-which(myPData2$Stimulation=="iRBC" & myPData2$casecon=="non_malaria" & myPData2$groupnb=="2" & myPData2$timepoint=="M3")

namesiRBCmalaria<-which(myPData2$Stimulation=="iRBC" & myPData2$casecon=="malaria" & myPData2$groupnb=="2" & myPData2$timepoint=="M3")

names<- colnames(listafinal_04)[c(namesiRBCnon,namesiRBCmalaria,namesuRBCnon,namesuRBCmalaria)]

rownames(iRBC_CONTROLvsuRBC_CONTROL_fdr_pval005)<-iRBC_CONTROLvsuRBC_CONTROL_fdr_pval005$ID

png(paste(paste(paste(paste("Results/malaria_iRBC_CONTROL","heatmap_iRBC_CONTROLvsuRBC_CONTROL",sep="/"),adjust,sep="_"),pvalue,sep="_"),".png",sep=""),width = 1200, height = 1200)

  heatmap(listafinal_04[rownames(iRBC_CONTROLvsuRBC_CONTROL_fdr_pval005),names],distfun=function(c)
  {Dist(c,method="correlation")},hclustfun=function(d)
  {hclust(d,method="average")},Colv=cov(listafinal_04[rownames(iRBC_CONTROLvsuRBC_CONTROL_fdr_pval005),names]),Rowv=cov(listafinal_04[rownames(iRBC_CONTROLvsuRBC_CONTROL_fdr_adjpval005),names]))

dev.off()

save(list=ls(all=TRUE),file="Results/malaria_iRBC_CONTROL/iRBC_CONTROLvsuRBC_CONTROL.RData")

#Volcano Plot

load("Results/malaria_iRBC_CONTROL/iRBC_CONTROLvsuRBC_CONTROL.RData")

gene_list<-genes_fdr

gene_list$threshold = as.factor( gene_list$P.Value < 0.05)

gene_list$threshold_fdr = as.factor( gene_list$adj.P.Val < 0.05)

gene_list$colour[which(gene_list$threshold=="FALSE")]<-"black"

gene_list$colour[which(gene_list$threshold=="TRUE" & gene_list$threshold_fdr=="FALSE")]<-"grey"

gene_list$colour[which(gene_list$threshold=="TRUE" & gene_list$threshold_fdr=="TRUE")]<-"red"

png("Results/malaria_iRBC_CONTROL/VolcanoPlot_iRBC_CONTROLvsuRBC_CONTROL.png")
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(P.Value))) +
    geom_point(alpha=1,pch=20, size=3,col=gene_list$colour) +
    theme(legend.position = "none") +
    xlim(c(-10, 10)) + ylim(c(0, 15)) +
    xlab("log2 fold change") + ylab("log p-value")
    main="iRBC_CONTROLvsuRBC_CONTROL"
g
dev.off()

# GO annotation 

load("GO_hugene21.RData")

geneList<-genes$P.Value

names(geneList)<-rownames(genes)

GOdataBP <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=topDiffGenes ,annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOdataMF <- new("topGOdata", ontology = "MF", allGenes = geneList,geneSel=topDiffGenes ,annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOdataCC <- new("topGOdata", ontology = "CC", allGenes = geneList,geneSel=topDiffGenes ,annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Enrichment Analysis

resultFisherBP<-runTest(GOdataBP,algorithm="classic",statistic="fisher")

resultFisherMF<-runTest(GOdataMF,algorithm="classic",statistic="fisher")

resultFisherCC<-runTest(GOdataCC,algorithm="classic",statistic="fisher")

allResBP<-GenTable(GOdataBP,
		 pValue_Fisher=resultFisherBP,
		 orderBy="pValueFisher",
		 ranksOf="pValueFisher",
		 numChar=1000,
		 topNodes=sum(resultFisherBP@score<0.05))
		 
allResMF<-GenTable(GOdataMF,
		 pValue_Fisher=resultFisherMF,
		 orderBy="pValueFisher",
		 ranksOf="pValueFisher",
		 numChar=1000,
		 topNodes=sum(resultFisherMF@score<0.05))
		 
allResCC<-GenTable(GOdataCC,
		 pValue_Fisher=resultFisherCC,
		 orderBy="pValueFisher",
		 ranksOf="pValueFisher",
		 numChar=1000,
		 topNodes=sum(resultFisherCC@score<0.05))

allResBP$p.adj_fdr<-round(p.adjust(allResBP$pValue_Fisher, method = "fdr", n =length(GOdataBP@graph@nodes)),4)

allResMF$p.adj_fdr<-round(p.adjust(allResMF$pValue_Fisher, method = "fdr", n =length(GOdataMF@graph@nodes)),4)

allResCC$p.adj_fdr<-round(p.adjust(allResCC$pValue_Fisher, method = "fdr", n =length(GOdataCC@graph@nodes)),4)

allResBP$p.adj_bonferroni<-round(p.adjust(allResBP$pValue_Fisher, method = "bonferroni", n = length(GOdataBP@graph@nodes)),4)

allResMF$p.adj_bonferroni<-round(p.adjust(allResMF$pValue_Fisher, method = "bonferroni", n = length(GOdataMF@graph@nodes)),4)

allResCC$p.adj_bonferroni<-round(p.adjust(allResCC$pValue_Fisher, method = "bonferroni", n = length(GOdataCC@graph@nodes)),4)


goIDBP<-allResBP[,"GO.ID"]

goIDMF<-allResMF[,"GO.ID"]

goIDCC<-allResCC[,"GO.ID"]



gtBP<-printGenes(GOdataBP,whichTerms=goIDBP,chip="hugene21sttranscriptcluster.db",geneCutOff=max(allResBP$Significant),pvalCutOff=0.05)

gtMF<-printGenes(GOdataMF,whichTerms=goIDMF,chip="hugene21sttranscriptcluster.db",geneCutOff=max(allResMF$Significant),pvalCutOff=0.05)

gtCC<-printGenes(GOdataCC,whichTerms=goIDCC,chip="hugene21sttranscriptcluster.db",geneCutOff=max(allResCC$Significant),pvalCutOff=0.05)


for (i in c(1:length(goIDBP))){
  
  allResBP$gene_SYMBOL[[i]]<-paste(as.vector(gtBP[[i]]$Symbol.id), collapse=",")
  
  }

for (i in c(1:length(goIDMF))){
  
  allResMF$gene_SYMBOL[[i]]<-paste(as.vector(gtMF[[i]]$Symbol.id), collapse=",")
  
  }
  
for (i in c(1:length(goIDCC))){
  
  allResCC$gene_SYMBOL[[i]]<-paste(as.vector(gtCC[[i]]$Symbol.id), collapse=",")
  
  }
  
  
write.csv(allResCC,file="Results/malaria_iRBC_CONTROL/GO_CC_iRBC_CONTROLvsuRBC_CONTROL.csv")
  
write.csv(allResMF,file="Results/malaria_iRBC_CONTROL/GO_MF_iRBC_CONTROLvsuRBC_CONTROL.csv")
  
write.csv(allResBP,file="Results/malaria_iRBC_CONTROL/GO_BP_iRBC_CONTROLvsuRBC_CONTROL.csv")

save(list=ls(all=TRUE),file="Results/malaria_iRBC_CONTROL/GO_iRBC_CONTROLvsuRBC_CONTROL.RData")





############################################################################################
####### Analysis of RTS,S immunogenicity (RTS,S vs Comparators post-vaccination) ###########
############################################################################################

load("Results/SysMalVac_December.RData")

system("mkdir -p Results/CSP_RTS_CONTROL")  #RTS,S vs Comparators post-vaccination (M3) for CSP stimulation

system("mkdir -p Results/iRBC_RTS_CONTROL") #RTS,S vs Comparators post-vaccination (M3) for PfRBC stimulation

# limma
experimento <- factor(paste(exprslistafinal_04@phenoData$Stimulation,
			    exprslistafinal_04@phenoData$timepoint,
			    exprslistafinal_04@phenoData$groupnb,sep=""))

design <- model.matrix(~ 0 + experimento)

colnames(design) <- c(levels(experimento))

fit<-lmFit(exprslistafinal_04,design)


# contrasts
contrast.matrix <- makeContrasts(
  
  CSP_RTS_CONTROL=((CSPM31-DMSOM31)-(CSPM32-DMSOM32)),
  
  iRBC_RTS_CONTROL=((iRBCM31-uRBCM31)-(iRBCM32-uRBCM32)),
  
  levels=design)
  
fit2<-eBayes(contrasts.fit(fit,contrast.matrix))

save(contrast.matrix,fit,design,myPData2,fit2,exprslistafinal_04,listafinal_04,file="Results/SysMalVac_December_2comp.RData")


############################################################################## CSP_RTS_CONTROL ###########################################################################################################


load("Results/SysMalVac_December_2comp.RData")

adjust="bonferroni";pvalue=0.05

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="CSP_RTS_CONTROL",sort.by="P",number=dim(exprslistafinal_04)[1])

CSP_RTS_CONTROL_bonferroni_pval005<-genes_bonferroni[genes_bonferroni$P.Value<0.05,]

CSP_RTS_CONTROL_bonferroni_adjpval005<-genes_bonferroni[genes_bonferroni$adj.P.Val<0.05,]

fileout<-"CSP_RTS_CONTROL_bonferroni_adjpval005"

write.table(CSP_RTS_CONTROL_bonferroni_adjpval005[,c(1,2,5,6)],file=paste(paste("Results/CSP_RTS_CONTROL",fileout,sep="/"),".txt",sep=""))


fileout<-"CSP_RTS_CONTROL_bonferroni_pval005"

write.table(CSP_RTS_CONTROL_bonferroni_pval005[,c(1,2,4,5,6)],file=paste(paste("Results/CSP_RTS_CONTROL",fileout,sep="/"),".txt",sep=""))


adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="CSP_RTS_CONTROL",sort.by="P",number=dim(exprslistafinal_04)[1])

CSP_RTS_CONTROL_fdr_pval005<-genes_fdr[genes_fdr$P.Value<0.05,]

CSP_RTS_CONTROL_fdr_adjpval005<-genes_fdr[genes_fdr$adj.P.Val<0.05,]

fileout<-"CSP_RTS_CONTROL_fdr_adjpval005"

write.table(CSP_RTS_CONTROL_fdr_adjpval005[,c(1,2,5,6)],file=paste(paste("Results/CSP_RTS_CONTROL",fileout,sep="/"),".txt",sep=""))


fileout<-"CSP_RTS_CONTROL_fdr_pval005"

write.table(CSP_RTS_CONTROL_fdr_pval005[,c(1,2,4,5,6)],file=paste(paste("Results/CSP_RTS_CONTROL",fileout,sep="/"),".txt",sep=""))

genes<-genes_fdr

namesDMSO1<-which(myPData2$Stimulation=="DMSO" & myPData2$groupnb=="1" & myPData2$timepoint=="M3")

namesDMSO2<-which(myPData2$Stimulation=="DMSO" & myPData2$groupnb=="2" & myPData2$timepoint=="M3")

namesCSP1<-which(myPData2$Stimulation=="CSP" & myPData2$groupnb=="1" & myPData2$timepoint=="M3")

namesCSP2<-which(myPData2$Stimulation=="CSP" & myPData2$groupnb=="2" & myPData2$timepoint=="M3")

names<- colnames(listafinal_04)[c(namesCSP1,namesCSP2,namesDMSO1,namesDMSO2)]

png(paste(paste(paste(paste("Results/CSP_RTS_CONTROL","heatmap_CSP_RTS_CONTROL",sep="/"),adjust,sep="_"),pvalue,sep="_"),".png",sep=""),width = 1200, height = 1200)

  heatmap(listafinal_04[rownames(CSP_RTS_CONTROL_fdr_pval005),names],distfun=function(c)
  {Dist(c,method="correlation")},hclustfun=function(d)
  {hclust(d,method="average")},Colv=cov(listafinal_04[rownames(CSP_RTS_CONTROL_fdr_pval005),names]),Rowv=cov(listafinal_04[rownames(CSP_RTS_CONTROL_fdr_adjpval005),names]))

dev.off()

save(list=ls(all=TRUE),file="Results/CSP_RTS_CONTROL/CSP_RTS_CONTROL.RData")

#Volcano Plot

load("Results/CSP_RTS_CONTROL/CSP_RTS_CONTROL.RData")

gene_list<-genes_fdr

gene_list$threshold = as.factor( gene_list$P.Value < 0.05)

gene_list$threshold_fdr = as.factor( gene_list$adj.P.Val < 0.05)

gene_list$colour[which(gene_list$threshold=="FALSE")]<-"black"

gene_list$colour[which(gene_list$threshold=="TRUE" & gene_list$threshold_fdr=="FALSE")]<-"grey"

gene_list$colour[which(gene_list$threshold=="TRUE" & gene_list$threshold_fdr=="TRUE")]<-"red"

png("Results/CSP_RTS_CONTROL/VolcanoPlot_CSP_RTS_CONTROL.png")
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(P.Value))) +
    geom_point(alpha=1,pch=20, size=3,col=gene_list$colour) +
    theme(legend.position = "none") +
    xlim(c(-10, 10)) + ylim(c(0, 15)) +
    xlab("log2 fold change") + ylab("log p-value")
    main="CSP_RTS_CONTROL"
g
dev.off()

# GO annotation 

load("GO_hugene21.RData")

geneList<-genes$P.Value

names(geneList)<-rownames(genes)

GOdataBP <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=topDiffGenes ,annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOdataMF <- new("topGOdata", ontology = "MF", allGenes = geneList,geneSel=topDiffGenes ,annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOdataCC <- new("topGOdata", ontology = "CC", allGenes = geneList,geneSel=topDiffGenes ,annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Enrichment Analysis

resultFisherBP<-runTest(GOdataBP,algorithm="classic",statistic="fisher")

resultFisherMF<-runTest(GOdataMF,algorithm="classic",statistic="fisher")

resultFisherCC<-runTest(GOdataCC,algorithm="classic",statistic="fisher")

allResBP<-GenTable(GOdataBP,
		 pValue_Fisher=resultFisherBP,
		 orderBy="pValueFisher",
		 ranksOf="pValueFisher",
		 numChar=1000,
		 topNodes=sum(resultFisherBP@score<0.05))
		 
allResMF<-GenTable(GOdataMF,
		 pValue_Fisher=resultFisherMF,
		 orderBy="pValueFisher",
		 ranksOf="pValueFisher",
		 numChar=1000,
		 topNodes=sum(resultFisherMF@score<0.05))
		 
allResCC<-GenTable(GOdataCC,
		 pValue_Fisher=resultFisherCC,
		 orderBy="pValueFisher",
		 ranksOf="pValueFisher",
		 numChar=1000,
		 topNodes=sum(resultFisherCC@score<0.05))

allResBP$p.adj_fdr<-round(p.adjust(allResBP$pValue_Fisher, method = "fdr", n =length(GOdataBP@graph@nodes)),4)

allResMF$p.adj_fdr<-round(p.adjust(allResMF$pValue_Fisher, method = "fdr", n =length(GOdataMF@graph@nodes)),4)

allResCC$p.adj_fdr<-round(p.adjust(allResCC$pValue_Fisher, method = "fdr", n =length(GOdataCC@graph@nodes)),4)

allResBP$p.adj_bonferroni<-round(p.adjust(allResBP$pValue_Fisher, method = "bonferroni", n = length(GOdataBP@graph@nodes)),4)

allResMF$p.adj_bonferroni<-round(p.adjust(allResMF$pValue_Fisher, method = "bonferroni", n = length(GOdataMF@graph@nodes)),4)

allResCC$p.adj_bonferroni<-round(p.adjust(allResCC$pValue_Fisher, method = "bonferroni", n = length(GOdataCC@graph@nodes)),4)


goIDBP<-allResBP[,"GO.ID"]

goIDMF<-allResMF[,"GO.ID"]

goIDCC<-allResCC[,"GO.ID"]



gtBP<-printGenes(GOdataBP,whichTerms=goIDBP,chip="hugene21sttranscriptcluster.db",geneCutOff=max(allResBP$Significant),pvalCutOff=0.05)

gtMF<-printGenes(GOdataMF,whichTerms=goIDMF,chip="hugene21sttranscriptcluster.db",geneCutOff=max(allResMF$Significant),pvalCutOff=0.05)

gtCC<-printGenes(GOdataCC,whichTerms=goIDCC,chip="hugene21sttranscriptcluster.db",geneCutOff=max(allResCC$Significant),pvalCutOff=0.05)


for (i in c(1:length(goIDBP))){
  
  allResBP$gene_SYMBOL[[i]]<-paste(as.vector(gtBP[[i]]$Symbol.id), collapse=",")
  
  }

for (i in c(1:length(goIDMF))){
  
  allResMF$gene_SYMBOL[[i]]<-paste(as.vector(gtMF[[i]]$Symbol.id), collapse=",")
  
  }
  
for (i in c(1:length(goIDCC))){
  
  allResCC$gene_SYMBOL[[i]]<-paste(as.vector(gtCC[[i]]$Symbol.id), collapse=",")
  
  }
  
  
write.csv(allResCC,file="Results/CSP_RTS_CONTROL/GO_CC_CSP_RTS_CONTROL.csv")
  
write.csv(allResMF,file="Results/CSP_RTS_CONTROL/GO_MF_CSP_RTS_CONTROL.csv")
  
write.csv(allResBP,file="Results/CSP_RTS_CONTROL/GO_BP_CSP_RTS_CONTROL.csv")

save(list=ls(all=TRUE),file="Results/CSP_RTS_CONTROL/GO_CSP_RTS_CONTROL.RData")



######################################################  iRBC_RTS_CONTROL ##################################################################

load("Results/SysMalVac_December_2comp.RData")

adjust="bonferroni";pvalue=0.05

adjust="bonferroni";pvalue=0.05

genes_bonferroni<-topTable(fit2,adjust.method=adjust,coef="iRBC_RTS_CONTROL",sort.by="P",number=dim(exprslistafinal_04)[1])

iRBC_RTS_CONTROL_bonferroni_pval005<-genes_bonferroni[genes_bonferroni$P.Value<0.05,]

iRBC_RTS_CONTROL_bonferroni_adjpval005<-genes_bonferroni[genes_bonferroni$adj.P.Val<0.05,]

fileout<-"iRBC_RTS_CONTROL_bonferroni_adjpval005"

write.table(iRBC_RTS_CONTROL_bonferroni_adjpval005[,c(1,2,5,6)],file=paste(paste("Results/iRBC_RTS_CONTROL",fileout,sep="/"),".txt",sep=""))


fileout<-"iRBC_RTS_CONTROL_bonferroni_pval005"

write.table(iRBC_RTS_CONTROL_bonferroni_pval005[,c(1,2,4,5,6)],file=paste(paste("Results/iRBC_RTS_CONTROL",fileout,sep="/"),".txt",sep=""))


adjust="fdr";pvalue=0.05

genes_fdr<-topTable(fit2,adjust.method=adjust,coef="iRBC_RTS_CONTROL",sort.by="P",number=dim(exprslistafinal_04)[1])

iRBC_RTS_CONTROL_fdr_pval005<-genes_fdr[genes_fdr$P.Value<0.05,]

iRBC_RTS_CONTROL_fdr_adjpval005<-genes_fdr[genes_fdr$adj.P.Val<0.05,]

fileout<-"iRBC_RTS_CONTROL_fdr_adjpval005"

write.table(iRBC_RTS_CONTROL_fdr_adjpval005[,c(1,2,5,6)],file=paste(paste("Results/iRBC_RTS_CONTROL",fileout,sep="/"),".txt",sep=""))


fileout<-"iRBC_RTS_CONTROL_fdr_pval005"

write.table(iRBC_RTS_CONTROL_fdr_pval005[,c(1,2,4,5,6)],file=paste(paste("Results/iRBC_RTS_CONTROL",fileout,sep="/"),".txt",sep=""))

genes<-genes_fdr

namesuRBC1<-which(myPData2$Stimulation=="uRBC" & myPData2$groupnb=="1" & myPData2$timepoint=="M3")

namesuRBC2<-which(myPData2$Stimulation=="uRBC" & myPData2$groupnb=="2" & myPData2$timepoint=="M3")

namesiRBC1<-which(myPData2$Stimulation=="iRBC" & myPData2$groupnb=="1" & myPData2$timepoint=="M3")

namesiRBC2<-which(myPData2$Stimulation=="iRBC" & myPData2$groupnb=="2" & myPData2$timepoint=="M3")

names<- colnames(listafinal_04)[c(namesiRBC1,namesiRBC2,namesuRBC1,namesuRBC2)]

png(paste(paste(paste(paste("Results/iRBC_RTS_CONTROL","heatmap_iRBC_RTS_CONTROL",sep="/"),adjust,sep="_"),pvalue,sep="_"),".png",sep=""),width = 1200, height = 1200)

  heatmap(listafinal_04[rownames(iRBC_RTS_CONTROL_fdr_pval005),names],distfun=function(c)
  {Dist(c,method="correlation")},hclustfun=function(d)
  {hclust(d,method="average")},Colv=cov(listafinal_04[rownames(iRBC_RTS_CONTROL_fdr_pval005),names]),Rowv=cov(listafinal_04[rownames(iRBC_RTS_CONTROL_fdr_adjpval005),names]))

dev.off()

save(list=ls(all=TRUE),file="Results/iRBC_RTS_CONTROL/iRBC_RTS_CONTROL.RData")

#Volcano Plot

load("Results/iRBC_RTS_CONTROL/iRBC_RTS_CONTROL.RData")

gene_list<-genes_fdr

gene_list$threshold = as.factor( gene_list$P.Value < 0.05)

gene_list$threshold_fdr = as.factor( gene_list$adj.P.Val < 0.05)

gene_list$colour[which(gene_list$threshold=="FALSE")]<-"black"

gene_list$colour[which(gene_list$threshold=="TRUE" & gene_list$threshold_fdr=="FALSE")]<-"grey"

gene_list$colour[which(gene_list$threshold=="TRUE" & gene_list$threshold_fdr=="TRUE")]<-"red"

png("Results/iRBC_RTS_CONTROL/VolcanoPlot_iRBC_RTS_CONTROL.png")
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(P.Value))) +
    geom_point(alpha=1,pch=20, size=3,col=gene_list$colour) +
    theme(legend.position = "none") +
    xlim(c(-10, 10)) + ylim(c(0, 15)) +
    xlab("log2 fold change") + ylab("log p-value")
    main="iRBC_RTS_CONTROL"
g
dev.off()

#GO annotation 

load("GO_hugene21.RData")

geneList<-genes$P.Value

names(geneList)<-rownames(genes)

GOdataBP <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=topDiffGenes ,annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOdataMF <- new("topGOdata", ontology = "MF", allGenes = geneList,geneSel=topDiffGenes ,annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOdataCC <- new("topGOdata", ontology = "CC", allGenes = geneList,geneSel=topDiffGenes ,annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Enrichment Analysis

resultFisherBP<-runTest(GOdataBP,algorithm="classic",statistic="fisher")

resultFisherMF<-runTest(GOdataMF,algorithm="classic",statistic="fisher")

resultFisherCC<-runTest(GOdataCC,algorithm="classic",statistic="fisher")

allResBP<-GenTable(GOdataBP,
		 pValue_Fisher=resultFisherBP,
		 orderBy="pValueFisher",
		 ranksOf="pValueFisher",
		 numChar=1000,
		 topNodes=sum(resultFisherBP@score<0.05))
		 
allResMF<-GenTable(GOdataMF,
		 pValue_Fisher=resultFisherMF,
		 orderBy="pValueFisher",
		 ranksOf="pValueFisher",
		 numChar=1000,
		 topNodes=sum(resultFisherMF@score<0.05))
		 
allResCC<-GenTable(GOdataCC,
		 pValue_Fisher=resultFisherCC,
		 orderBy="pValueFisher",
		 ranksOf="pValueFisher",
		 numChar=1000,
		 topNodes=sum(resultFisherCC@score<0.05))

allResBP$p.adj_fdr<-round(p.adjust(allResBP$pValue_Fisher, method = "fdr", n =length(GOdataBP@graph@nodes)),4)

allResMF$p.adj_fdr<-round(p.adjust(allResMF$pValue_Fisher, method = "fdr", n =length(GOdataMF@graph@nodes)),4)

allResCC$p.adj_fdr<-round(p.adjust(allResCC$pValue_Fisher, method = "fdr", n =length(GOdataCC@graph@nodes)),4)

allResBP$p.adj_bonferroni<-round(p.adjust(allResBP$pValue_Fisher, method = "bonferroni", n = length(GOdataBP@graph@nodes)),4)

allResMF$p.adj_bonferroni<-round(p.adjust(allResMF$pValue_Fisher, method = "bonferroni", n = length(GOdataMF@graph@nodes)),4)

allResCC$p.adj_bonferroni<-round(p.adjust(allResCC$pValue_Fisher, method = "bonferroni", n = length(GOdataCC@graph@nodes)),4)


goIDBP<-allResBP[,"GO.ID"]

goIDMF<-allResMF[,"GO.ID"]

goIDCC<-allResCC[,"GO.ID"]



gtBP<-printGenes(GOdataBP,whichTerms=goIDBP,chip="hugene21sttranscriptcluster.db",geneCutOff=max(allResBP$Significant),pvalCutOff=0.05)

gtMF<-printGenes(GOdataMF,whichTerms=goIDMF,chip="hugene21sttranscriptcluster.db",geneCutOff=max(allResMF$Significant),pvalCutOff=0.05)

gtCC<-printGenes(GOdataCC,whichTerms=goIDCC,chip="hugene21sttranscriptcluster.db",geneCutOff=max(allResCC$Significant),pvalCutOff=0.05)


for (i in c(1:length(goIDBP))){
  
  allResBP$gene_SYMBOL[[i]]<-paste(as.vector(gtBP[[i]]$Symbol.id), collapse=",")
  
  }

for (i in c(1:length(goIDMF))){
  
  allResMF$gene_SYMBOL[[i]]<-paste(as.vector(gtMF[[i]]$Symbol.id), collapse=",")
  
  }
  
for (i in c(1:length(goIDCC))){
  
  allResCC$gene_SYMBOL[[i]]<-paste(as.vector(gtCC[[i]]$Symbol.id), collapse=",")
  
  }
  
  
write.csv(allResCC,file="Results/iRBC_RTS_CONTROL/GO_CC_iRBC_RTS_CONTROL.csv")
  
write.csv(allResMF,file="Results/iRBC_RTS_CONTROL/GO_MF_iRBC_RTS_CONTROL.csv")
  
write.csv(allResBP,file="Results/iRBC_RTS_CONTROL/GO_BP_iRBC_RTS_CONTROL.csv")

save(list=ls(all=TRUE),file="Results/iRBC_RTS_CONTROL/GO_iRBC_RTS_CONTROL.RData")

###########################################################################################################################################################################################
##########################################################################################################################################################################################
