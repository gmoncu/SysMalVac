fixed.effect<-function(eset,design,adjust,fileout,pvalue=0.05,dataNorm,path=getwd()){ 

  options(warn=-1)
  
  pkgs<-c("limma","marray","amap","affy")

  sapply(pkgs,version,FUN=function(pkgs,version){
		if (is.na(packageDescription(pkgs))) {
		source("http://bioconductor.org/biocLite.R")
		biocLite(pkgs)}
		})

  lapply(pkgs, require, character.only=T)	 



   contrasts<-readline("Contrast:") 
   contrasts<-as.factor(unlist(strsplit(tolower(contrasts), ",")))
   ok<-0
   for (i in c(1:length(contrasts))){
      if (summary(strsplit(as.character(contrasts[i]),"-")[[1]][1]==colnames(design))[2]!=length(colnames(design)) & 
  	   summary(strsplit(as.character(contrasts[i]),"-")[[1]][2]==colnames(design))[2]!=length(colnames(design)) ){
   ok<-ok+1
   	}
    }

   
   if (ok==length(contrasts)){
    
      names <- list.celfiles(path=path)
 
      for (i in c(1:length(names))){
	j<-0
	  while (substr(names[i],j,j)!=".")
	  {
	  j=j+1
	  }
	  names[i]<-substr(names[i],1,(j-1))
      }
    
      
      fit<-lmFit(eset,design)
      contrast.matrix<-makeContrasts(contrasts=contrasts,levels=design)
      fit2<-eBayes(contrasts.fit(fit,contrast.matrix))
      genes<-topTable(fit2,adjust.method=adjust,number=dim(eset)[1],p.value=pvalue)
      genes[,"P.Value"]<-round(genes[,"P.Value"],digits=3)
      write.table(genes,file=paste(paste("Resultados",fileout,sep="/"),".txt",sep=""))
      rgb<-maPalette(low="green",high="red",mid="black")
      png(paste(paste(paste(paste("Resultados",as.character(contrasts),sep="/"),adjust,sep="_"),pvalue,sep="_"),".png",sep=""),width = 1200, height = 1200)
      heatmap(dataNorm[genes$ID,],col=rgb,distfun=function(c)
	{Dist(c,method="correlation")},hclustfun=function(d)
	{hclust(d,method="average")},labCol=names,Colv=cov(dataNorm[genes$ID,]),Rowv=cov(dataNorm[genes$ID,]))
	
      dev.off()
    
    save(list=ls(all=TRUE),file=paste(paste("Resultados",fileout,sep="/"),".RData",sep=""))
    
  } else {
    print ("Contrasts and matrix design names are not equal")
  }

} 