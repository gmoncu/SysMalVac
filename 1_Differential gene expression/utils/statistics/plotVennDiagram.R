plotVennDiagram<-function(filein1,filein2,filein3,names.labels,color,colint=colors()[c(137)]){

options(warn=-1)

pkgs<-c("venneuler")

sapply(pkgs,version,FUN=function(pkgs,version){
		if (is.na(packageDescription(pkgs))) {
		source("http://bioconductor.org/biocLite.R")
		biocLite(pkgs)}
		})

lapply(pkgs, require, character.only=T)	  

GENES1<-read.table(paste("Resultados",filein1,sep="/"))
GENES2<-read.table(paste("Resultados",filein2,sep="/"))
if (!is.null(filein3)){
GENES3<-read.table(paste("Resultados",filein3,sep="/"))


pos<-pos2<-pos3<-pos123<-NULL
for (i in c(1:dim(GENES1))){
    pos<-c(pos,which(GENES1[i,1]==GENES2[,1]))
    pos123<-c(pos123,which(GENES1[which(GENES1[i,1]==GENES2[,1]),1]==GENES3[,1]))
    pos3<-c(pos,which(GENES1[i,1]==GENES3[,1]))
  }
  
for (i in c(1:dim(GENES2))){
    pos2<-c(pos,which(GENES2[i,1]==GENES3[,1]))
    
  }
  
 v<-venneuler(c("table1"=dim(GENES1)[1],"table2"=dim(GENES2)[1],"table3"=dim(GENES3)[1],"table1&table2"=length(pos)),
		"table2&table3"=length(pos2),"table1&table3"=length(pos3),"table1&table2&table3"=length(pos123),
		names=c("table1","table2","table3"))
		
v$labels<-c(as.character(paste (names.labels[1],dim(GENES1)[1],sep=" ")),as.character(paste (names.labels[2],dim(GENES2)[1],as.character(paste (names.labels[2],dim(GENES3)[1],sep=" ")),sep=" ")))

png("Resultados/VennDiagram.png")
plot.VennDiagram(v,color,alpha=0.3,border=TRUE,main="Venn Diagram")
legend("bottomright",legend=names.labels, fill=c(color,colint)) 
dev.off()

}else{

  pos<-NULL
for (i in c(1:dim(GENES1))){
    pos<-c(pos,which(GENES1[i,1]==GENES2[,1]))
  }

 v<-venneuler(c("table1"=dim(GENES1)[1],"table2"=dim(GENES2)[1],"table1&table2"=length(pos)),
		names=c("table1","table2"))
 v$labels<-c("+/+","-/-")

png("../../../../miRNA4.1/CEL/Resultados/1.CellStimulus/VennDiagrammiRNAHCT116.png")
 
plot.VennDiagram(v,color,alpha=0.3,border=TRUE,main="Venn Diagram miRNA HCT116")

		
v$labels<-c(as.character(paste (names.labels[1],dim(GENES1)[1],sep=" ")),
	    as.character(paste (names.labels[2],dim(GENES2)[1],sep=" ")),
	    as.character(paste(paste(names.labels[1],names.labels[2],sep=" ∩ "),length(pos),sep=" ")))

legend("bottomleft",legend=v$labels, fill=c("rosybrown1","gray48","rosybrown4"))

dev.off()
	
}
  
  
}