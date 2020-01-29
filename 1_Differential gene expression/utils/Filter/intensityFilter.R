######################################################################################
##########################Intensity Filter RAW data without Controls##################
######################################################################################

intensityFilter<-function(path=path,kvalue,Avalue){


options(warn=-1)

pkgs<-c("affy","genefilter","GenomicRanges","oligo","oligoClasses")

sapply(pkgs,version,FUN=function(pkgs,version){
		if (is.na(packageDescription(pkgs))) {
		source("http://bioconductor.org/biocLite.R")
		biocLite(pkgs)}
		})

lapply(pkgs, require, character.only=T)	 

pkgsDB<-c("pd.hugene.2.1.st","pd.ragene.1.0.st.v1","pd.mogene.2.1.st")

 sapply(pkgsDB,version,FUN=function(pkgsDB,version){
 		if (is.na(packageDescription(pkgsDB))) {
 		source("http://bioconductor.org/biocLite.R")
 		biocLite(pkgsDB)}
 		})
		
lapply(pkgsDB, require, character.only=T)

#Read CEL files
database<-IDposcontrol<-NULL

dataCEL<-ReadAffy(celfile.path=path)

if (slot(dataCEL,"annotation")=="hugene21st"){

require("hugene21stcdf")
database <- db(pd.hugene.2.1.st)

load("posControl_hugene21st.RData")

} else if(slot(dataCEL,"annotation")=="mogene21st"){.
  
database <- db(pd.mogene.2.1.st)


} else if(slot(dataCEL,"annotation")=="ragene10stv1"){

database <- db(pd.ragene.1.0.st.v1)

} 

if (is.null(database)){

      print ("No annotation for these Platform")
      
  } else {
  
	  if (!is.na(packageDescription(paste(slot(dataCEL,"annotation"),"cdf",sep="")))){
  
	    #Remove probe control
      
	    IDscontrol<- dbGetQuery(database, "select fsetid from featureSet where type in ('2','3','4','5','6','7','9');")[,1]

	    #dbGetquery no elimina controles positivos para pd.hugene.2.1.st.  Tenemos que eliminar manulmente los controles positivos posControl_hugene21st.RData
	    
	       
	    load("posControl_hugene21st.RData")
	    
	    IDscontrol<- c(IDscontrol,IDposcontrol)
	    
	    IDs<-as.character(IDscontrol)

	    exprsdataCEL<-computeExprSet(dataCEL, 
				    pmcorrect.method="pmonly",
				    summary.method="avgdiff")
				  
				    
	    filtexprsdataCEL = featureFilter(exprsdataCEL,
					     require.entrez=FALSE, 
					     require.GOBP=FALSE, 
					     require.GOCC=FALSE, 
					     require.GOMF=FALSE, 
					     remove.dupEntrez=FALSE,
					     feature.exclude=IDs)     

	     #Low Intensity Cutoff Filter.Percentage Cutoffs 100% 

	    filter<-kOverA(k=kvalue,A=Avalue)
  
	    ffilter<-filterfun(filter)

	    goodIntensity<-genefilter(filtexprsdataCEL,ffilter)

	    out<-which(goodIntensity=="FALSE")

	    filtIntCondataCEL<-filtexprsdataCEL[-out,]

	    return(filtIntCondataCEL)
	    
	    } else {
	    
	    print (paste(" Install annotation database",paste(slot(dataCEL,"annotation"),"cdf",sep=""),sep=" "))
	    
	    }
	    
  }
}