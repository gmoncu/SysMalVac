  plotPCA_prcomp<-function(eset,path,fileout,pos.legend,nPCs=nPCs,Xmax,Xmin,Ymax,Ymin){
 
 x <- readline("Number of functional categories in this gene expression study:")
 x <- as.numeric(unlist(strsplit(x, ",")))
 
  if (x>2) {

      print ("The maximum number of categories for PCA plot are 2")
  
  } else {
  
      y <- readline("Levels points:")
      class.factor1 <- as.factor(unlist(strsplit(y, ",")))
      z <- readline("Levels elipse:")
      class.factor2 <- as.factor(unlist(strsplit(z, ",")))
      p<-readline("Names Legend:")
      class.legend <- as.character(unlist(strsplit(p, ",")))
  
  
  options(warn=-1)

  pkgs<-c("pcaMethods","pls","stats","graphics","rgl","affy","FactorMineR")

  sapply(pkgs,version,FUN=function(pkgs,version){
		if (is.na(packageDescription(pkgs))) {
		source("http://bioconductor.org/biocLite.R")
		biocLite(pkgs)}
		})

  lapply(pkgs, require, character.only=T)	  

  names <- list.celfiles(path=path)
 
  for (i in c(1:length(names))){
    j<-0
    while (substr(names[i],j,j)!=".")
      {
      j=j+1
      }
    names[i]<-substr(names[i],1,(j-1))
  }

  dataPCA<-prcomp(t(eset),cor=TRUE)
  variance<-sprintf("%.2f",((dataPCA$sd)^2/sum(dataPCA$sd^2)*100))
  rownames(dataPCA$x)<-names

  #labelColors= topo.colors( length(levels(factor(class.factor1))))
  
  labelColors=class.factor1
  
  png(paste("Resultados",paste(paste("pca",fileout,sep=""),".png",sep=""),sep="/"))

  par(mfrow=c(nPCs*0.5*(nPCs-1),1))
  
  for (ii in  c(1:nPCs)){
  
  jj<-ii+1
     
    while(jj<=nPCs){
    
    if (is.na(class.factor2[1])==FALSE){
	boxColors=rainbow(length(levels(factor(class.factor2))))
	xmax<-xmin<-ymax<-ymin<-NULL
	
	for (k in c(1:length(levels(factor(class.factor2))))){
	  
	  y<-dataPCA$x[names[which(class.factor2==k)],ii]
	  x<-dataPCA$x[names[which(class.factor2==k)],jj]
	  
	  fit<-fit.ellipse (x,y)
	  elipse<-get.ellipse(fit,n=360)
      
		    }
          
      
       y<-dataPCA$x[names[which(class.factor2==1)],ii]
       x<-dataPCA$x[names[which(class.factor2==1)],jj]
	fit<-fit.ellipse (x,y)
	  elipse<-get.ellipse(fit,n=360)
      
      plot(elipse,type="l",main="Principal Component Analysis",xlim=c(Xmin,Xmax),ylim=c(Ymin,Ymax),ylab=paste(paste(paste("pca",ii,sep= ""),variance[ii],sep="(")," %)",sep=""),xlab=paste(paste(paste("pca",jj,sep= ""),variance[jj],sep=" ("),"%)",sep=""),col=1)
     
      
      
  legend(pos.legend,class.legend,
	    col=c(topo.colors( length(levels(factor(class.factor1)))),c(1:length(levels(factor(class.factor2))))),
	    lty=c(rep(0,length(levels(factor(class.factor1)))),rep(1,length(levels(factor(class.factor2))))),
	    lwd=c(rep(2,length(levels(factor(class.factor1)))),rep(2,length(levels(factor(class.factor2))))),
	    pch=c(rep(10,length(levels(factor(class.factor1)))),rep(-1,length(levels(factor(class.factor2))))),
	    horiz=FALSE,cex=1,merge=FALSE)
      
      
      for (i in c(2:length(levels(factor(class.factor2))))){
  
	y<-dataPCA$x[names[which(class.factor2==i)],ii]
	x<-dataPCA$x[names[which(class.factor2==i)],jj]	
	fit<-fit.ellipse (x[-12],y[-12])
	elipse<-get.ellipse(fit,n=360)
	points(elipse,type="l",col=i)  
  
      }
      
      for (i in c(1: length(levels(factor(class.factor1))))){

	y<-dataPCA$x[names[which(class.factor1==i)],ii]
	x<-dataPCA$x[names[which(class.factor1==i)],jj]
	points(x,y,col=topo.colors(i)[i],lwd=2,pch=10)
	
	}
      
      jj<-jj+1
    
    } else { 
    
      y<-dataPCA$x[names[which(class.factor1==1)],ii]
      x<-dataPCA$x[names[which(class.factor1==1)],jj]
      
      plot(x,y,col=topo.colors(1),xlim=c(Xmin,Xmax),ylim=c(Ymin,Ymax),lwd=2,pch=10,main="Principal Component Analysis",ylab=paste(paste(paste("pca",ii,sep= ""),variance[ii],sep="(")," %)",sep=""),xlab=paste(paste(paste("pca",jj,sep= ""),variance[jj],sep=" ("),"%)",sep=""))
	
     for (i in c(2: length(levels(factor(class.factor1))))){

	y<-dataPCA$x[names[which(class.factor1==i)],ii]
	x<-dataPCA$x[names[which(class.factor1==i)],jj]
	points(x,y,col=topo.colors(i)[i],lwd=2,pch=10)
	
	}
	jj<-jj+1
      
      }
    }
 }

 dev.off()

plot3d(dataPCA$x[1:dim(dataPCA$x)[2],],size=1.5,type="s", col=labelColors,main="PCA",
       xlab=paste(paste("pca 1 (",variance[1],sep=""),"%)",sep=""),
       ylab=paste(paste("pca 2 (",variance[2],sep=""),"%)",sep=""),
       zlab=paste(paste("pca 3 (",variance[3],sep=""),"%)",sep=""))
       
}	

}