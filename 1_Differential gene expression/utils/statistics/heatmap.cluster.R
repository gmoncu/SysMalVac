heatmap.cluster<-function(data,eset,filein,fileout,path,class.factor1,class.factor2,pos.legend,class.legend,sample){
 
options(warn=-1)

pkgs<-c("affy","marray","amap","gplots")


sapply(pkgs,version,FUN=function(pkgs,version){
		if (is.na(packageDescription(pkgs))) {
		source("http://bioconductor.org/biocLite.R")
		biocLite(pkgs)}
		})

lapply(pkgs, require, character.only=T)	  

y<-read.table(paste("Resultados",filein,sep="/"))
names <- list.celfiles(path=path)

for (i in c(1:length(names))){
    j<-0
    while (substr(names[i],j,j)!=".")
      {
      j=j+1
      }
    names[i]<-substr(names[i],1,(j-1))
  }

rgb<-maPalette(low="green",high="red",mid="black")
temp<-data[as.character(y$ID),paste(sample,".CEL",sep="")]

if (!is.null(class.factor2)){
  class.factor1<-as.factor(class.factor1)
  class.factor2<-as.factor(class.factor2)
  
  color<-rainbow(length(levels(class.factor2))*length(levels(class.factor1)))
  
  class.color<-vector(mode="character",length=length(class.factor1))
  legend.names<-vector(mode="character",length=length(color))

  for (i in c(1:length(levels(class.factor1)))){
    for (j in c(1:length(levels(class.factor2)))){
     indice<-which(class.factor2==levels(class.factor2)[j]) 
     class.color[indice[which(class.factor1[indice]==levels(class.factor1)[i])]]<-color[length(levels(class.factor2))*(i-1)+j]
     legend.names[length(levels(class.factor2))*(i-1)+j]<-paste(class.legend[i],class.legend[length(levels(class.factor1))+j],sep="_")
    }
  }

  }else{
  class.factor1<-as.factor(class.factor1)
  color<-rainbow(length(levels(class.factor1)))
  class.color<-vector(mode="character",length=length(class.factor1))
  legend.names<-vector(mode="character",length=length(color))
  
  for (i in c(1:length(levels(class.factor1)))){
    indice<-which(class.factor1==levels(class.factor1)[i]) 
    class.color[indice]<-color[i]
    legend.names[i]<-class.legend[i]
    }
  
  }

png(paste(paste("Resultados",fileout,sep="/"),".png",sep=""),width = 1000, height = 1000)

row.names(temp) <- rep("", nrow(temp))

heatmap.2(temp,distfun=function(c)
  {Dist(c,method="correlation")},hclustfun=function(d)
  {hclust(d,method="average")},key=TRUE,col=rgb,labCol=sample,trace="none",scale="col",ColSideColors=class.color)

legend(pos.legend,legend.names,
	    col=color,lty=rep(0,length(color)),lwd=rep(2,length(color)),pch=rep(15,length(color)),horiz=FALSE,cex=1,merge=FALSE)
	    	    
dev.off()	 




}

