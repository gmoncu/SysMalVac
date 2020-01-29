normFilter<-function(path,method,fdataCEL){

options(warn=-1)
pkgs<-c("affy","gcrma","vsn","genefilter")

sapply(pkgs,version,FUN=function(pkgs,version){
		if (is.na(packageDescription(pkgs))) {
		source("http://bioconductor.org/biocLite.R")
		biocLite(pkgs)}
		})

lapply(pkgs, require, character.only=T)	  
	  	  
#Reading CEL files
names <- list.celfiles(path=path)
dataCEL<-ReadAffy(celfile.path=path)

#Normalisation 
dataNorm<-switch(method,
			"vsn"=normVSN(dataCEL=dataCEL),
			"RMA"=normRMA(names=names,path=path),
			"GCRMA"=normGCRMA(names=names,path=path),
			"None"=print ("No normalisation"))

fdataNorm<-dataNorm[rownames(exprs(fdataCEL)),]
		
exprsdataNorm<-fdataNorm

return(exprsdataNorm)
}



