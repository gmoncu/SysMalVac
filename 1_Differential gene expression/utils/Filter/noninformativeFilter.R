noninformativeFilter<-function(data,typefilter,cutoff){

options(warn=-1)
pkgs<-c("genefilter")

sapply(pkgs,version,FUN=function(pkgs,version){
		if (is.na(packageDescription(pkgs))) {
		source("http://bioconductor.org/biocLite.R")
		biocLite(pkgs)}
		})

lapply(pkgs, require, character.only=T)	  


filter<-switch(typefilter,
			"IQR"=function(x)(IQR(x)>cutoff),
			"sd"=function(x)(sd(x)>cutoff),
			"variance"=function(x)(var(x)>cutoff),
			"None"=print ("No filter"))
			
ffun<-filterfun(filter)
list.filter<-genefilter(data,ffun)			

list.eset=data[list.filter,]

return(list.eset)
} 