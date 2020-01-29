F.statistic<-function(eset,Factor1,Factor2=NULL,names.class,fileout){
    options(warn=-1)
  
    pkgs<-c("limma")

    sapply(pkgs,version,FUN=function(pkgs,version){
		if (is.na(packageDescription(pkgs))) {
		source("http://bioconductor.org/biocLite.R")
		biocLite(pkgs)}
		})

    lapply(pkgs, require, character.only=T)	 

    F<-NULL
    if (!is.null(Factor2)){
   
    X<-model.matrix(~Factor1+Factor2+Factor1*Factor2)
    fit<-lmFit(lista.final[,c(1:length(Factor1))],X)
    
    for ( i in c(1:attr(X,"assign")[length(attr(X,"assign"))])){
	cont.ia<-diag(ncol(X))[,attr(X,"assign")==i]
	fit2<-eBayes(contrasts.fit(fit,cont.ia))
	F<-c(F,mean(fit2$F))
      }
    } else {
     
    X<-model.matrix(~Factor1)
    fit<-lmFit(lista.final[,c(1:length(Factor1))],X)
    cont.ia<-diag(ncol(X))[,attr(X,"assign")==1]
    fit2<-eBayes(contrasts.fit(fit,cont.ia))
    F<-mean(fit2$F)
    }
    
    F<-c(F,1)
    png(paste("Resultados",paste(fileout,".png",sep=""),sep="/"))
     
	bplt<-barplot(F,names.arg=names.class,col=c(2,3,4,5),border=1,angle=60,ylab="Median F ratio",main="Sources of Variation")
	abline(h=1,col=1)
	text(y= F+0.05,x=bplt, labels=as.numeric(round(F,digits=3)), xpd=TRUE)
	
    dev.off()
        
}