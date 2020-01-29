######################################################################################
###############################Quality Control########################################
######################################################################################

QualityControl<-function(path,method,control,intgroup){

options(warn=-1)
 
pkgs<-c("affy","gcrma","vsn","arrayQualityMetrics","hgu133plus2probe","oligo")

version<-c("1.36.1","2.30.0","3.26.0","3.14.0","2.11.0","1.22.0")

for (i in c(1:length(pkgs)) ) {   

      if (is.na(packageDescription(pkgs[i])) || packageVersion(pkgs[i])<version[i]) {

      source("http://bioconductor.org/biocLite.R")
      
      biocLite(pkgs[i])}

      }

lapply(pkgs, require, character.only=T)	  
	  
	  
#Reading CEL files
names <- list.celfiles(path=path)

dataCEL <- read.celfiles(list.celfiles())


print ("Reading Done 1/3")
#Normalisation 
dataNORM<-switch(method,
			"vsn"=normVSN(dataCEL=dataCEL),
			"RMA"=normRMA(names=names,path=path),
			"GCRMA"=normGCRMA(names=names,path=path),
            "None"=print ("No normalisation"))

print ("Normalisation Done 2/3")

if (is.null(intgroup)==FALSE){
  
  if (file.exists("phenoDATA.txt")){
  
  #myPData<-read.table("phenoDATA.txt", sep="\t",header=TRUE,row.names=1,as.is=TRUE)
  
  myPData<-read.table("phenoDATA.txt",header=TRUE)
  #Proyecto HypOrth
  #myPData<-myPData[order(myPData$Customer_code),]
  
  pData(dataCEL)<-myPData
  
  pData(dataNORM)<-myPData
  
  #Quality control
 
  print ("All Done 3/3")
  switch(control,
	      "Pre"= preQualityControl(data=dataCEL,intgroup),
	      "Post"= postQualityControl(data=dataNORM,intgroup),
	      "Both"= bothQualityControl(dataCEL=dataCEL,dataNORM=dataNORM,intgroup),
	      "None"= print ("No quality Control"),
	      )
  }else{
  print("Warning no file phenoDATA.txt")
    }
  }else{
  switch(control,
	      "Pre"= preQualityControl(data=dataCEL,intgroup),
	      "Post"= postQualityControl(data=dataNORM,intgroup),
	      "Both"= bothQualityControl(dataCEL=dataCEL,dataNORM=dataNORM,intgroup),
	      "None"= print ("No quality Control"),
	      )
  print ("All Done 3/3")
  }
 } 