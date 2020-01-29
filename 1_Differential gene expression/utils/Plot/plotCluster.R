plotCluster<-function(eset,path,class.factor1,class.factor2,class.factor3,class.legend,fileout,pos.legend="topleft",horiz=TRUE){

options(warn=-1)
pkgs<-c("dendroextras","affy","hyperSpec")

sapply(pkgs,version,FUN=function(pkgs,version){
		if (is.na(packageDescription(pkgs))) {
		source("http://bioconductor.org/biocLite.R")
		biocLite(pkgs)}
		})

lapply(pkgs, require, character.only=T)	  
d<-1-cor(eset,method="pearson")
hc1<-hclust(as.dist(d),method="average")

names <- list.celfiles(path=path)
 
for (i in c(1:length(names))){
  j<-0
  while (substr(names[i],j,j)!=".")
    {
    j=j+1
    }
   names[i]<-substr(names[i],1,(j-1))
 }

hc1$labels<-names
hcd<-as.dendrogram(hc1)

local({
        colLab <<- function(n) {
            if(is.leaf(n)) {
              a <- attributes(n)
              i <<- i+1
              attr(n, "nodePar") <- c(a$nodePar, list(lab.col = labelColors[i], lab.font=3,lab.cex=1, col=boxColors[i], cex=1, pch=16 ))
            }
            n
        }
        labelColors= topo.colors ( length(levels(factor(class.factor1))))[class.factor1[hc1$order]]
        boxColors= rainbow (length(levels(factor(class.factor2))))[class.factor2[hc1$order]]
        i <- 0
       })
       
clusDendro=dendrapply(hcd,colLab)
png(paste("Resultados",paste(fileout,".png",sep=""),sep="/"),width = 800, height = 900)
par (xpd = TRUE, mar = c (10, 4, 4, 3))
plot(clusDendro,axes=FALSE,sub=paste(paste("Average linkage G=",dim(eset)[1],sep=""),"genes",sep=" "),horiz=horiz,xlab="",main="Hierarchical clustering") 
if (is.null(class.factor3)==FALSE){
mark.dendrogram(hc1,groups=class.factor3,pos.marker=-0.05 ,pos.text=-0.055,height=0.002)
}

legend(pos.legend,class.legend,
      col=c( topo.colors( length(levels(factor(class.factor1)))),rainbow (length(levels(factor(class.factor2))))),
      lty=c(rep(1,length(levels(factor(class.factor1)))),rep(0,length(levels(factor(class.factor2))))),
      lwd=c(rep(3,length(levels(factor(class.factor1)))),rep(5,length(levels(factor(class.factor2))))),
      pch=c(rep(-1,length(levels(factor(class.factor1)))),rep(20,length(levels(factor(class.factor2))))),
      horiz=FALSE,cex=1,merge=TRUE)
dev.off()
}
