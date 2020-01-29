model.effect<-function(){

  options(warn=-1)
  
  pkgs<-c("limma")

  version<-c("3.14.4")


  sapply(pkgs,version,FUN=function(pkgs,version){
		if (is.na(packageDescription(pkgs)) || packageVersion(pkgs)<version) {
		source("http://bioconductor.org/biocLite.R")
		biocLite(pkgs)}
		})

  lapply(pkgs, require, character.only=T)	 
  
  x<-readline("Have you got random effects? (Yes/No)") 
  x<-as.factor(unlist(strsplit(x, ",")))
 
   
 if (tolower(x)=="no") {
    
    y<-readline("Number of factors:")
    y<-as.numeric(unlist(strsplit(y, ",")))
    if (y<=3){ 
	  z<-readline("Factor 1:")
	  z<-as.factor(unlist(strsplit(z, ",")))
	  
	  if (y==1){
	      
		design<-model.matrix(~-1+z)
		r<-readline("Column names:")
		r<-as.factor(unlist(strsplit(r, ",")))
		colnames(design)<-tolower(r)
		
		
		return(design)
		
		
	      } else if (y==2){
		
		w<-readline("Factor 2:")
		w<-as.factor(unlist(strsplit(w, ",")))
		
		#r<-readline("Column names:")
		#r<-as.factor(unlist(strsplit(r, ",")))
		
		yy<-readline("Which factors do you want to compare? (1/2/both)")
		if (is.na(as.numeric(yy))=="TRUE"){
		
			r<-readline("Column names:")
			r<-as.factor(unlist(strsplit(r, ",")))
		
			k<-0
			zw<-rr<-vector()
		
		
			for (i in c(1:length(levels(z)))){
		  
			    for (j in c(1:length(levels(w)))){
		    
			      k<-k+1
		    
			      zw[which(z==levels(z)[i])[ which(w[which(z==levels(z)[i])]==levels(w)[j]) ]]<- k
		    
			      if (is.na(summary(zw==k)["TRUE"])=="FALSE"){
			      rr<-cbind(rr,paste(r[i],r[length(levels(z))+j],sep=""))
			      
			      }
			    }
			  }
		
			design<-model.matrix(~-1+factor(zw))
			colnames(design)<-tolower(as.vector(rr))
	        
		} else {
		
		    if (as.numeric(yy)=="1"){
		
			design<-model.matrix(~-1+z)
			r<-readline("Column names:")
			r<-as.factor(unlist(strsplit(r, ",")))
			colnames(design)<-tolower(as.vector(r))
			
		    } else { 
		    
			r<-readline("Column names:")
			r<-as.factor(unlist(strsplit(r, ",")))
			design<-model.matrix(~-1+factor(w))
			colnames(design)<-tolower(as.vector(r))
		    }
		 }
		 
	         return(design)
	         
	       } else {
	        
		  w<-readline("Factor 2:")
		  w<-as.factor(unlist(strsplit(w, ",")))
		 
		  q<-readline("Factor 3:")
		  q<-as.factor(unlist(strsplit(q, ",")))
	       
		  yy<-readline("Which factors do you want to compare? (1/2/3/both)")
		 
		  if (is.na(as.numeric(yy))=="TRUE"){
		     
		      r<-readline("Column names:")
		      r<-as.factor(unlist(strsplit(r, ",")))
	        
		      k<-0
		      zwq<-indiceA<-indiceB<-vector()
		      rr<-NULL
		      for (i in c(1:length(levels(z)))){
		
			for (j in c(1:length(levels(w)))){
		  
			    for (p in c(1:length(levels(q)))){
		    
			      k<-k+1
		    
			      indiceA<-which(z==levels(z)[i])
			      indiceB<-indiceA[which(w[indiceA]==levels(w)[j])]
			      indiceC<-indiceB[which(w[indiceB]==levels(w)[p])]
			      zwq[indiceC]<-k
		      
			      if (is.na(summary(zwq==k)["TRUE"])=="FALSE"){
		      
			      rr<-cbind(rr,paste(paste(r[i],r[length(levels(z))+j],sep=""),r[length(levels(z))+length(levels(w))+p],sep=""))
			      }
		      
			    }
			  }
			}
		      design<-model.matrix(~-1+factor(zwq))
		      colnames(design)<-tolower(as.vector(rr))
		      
		      } else {
			
			if (as.numeric(yy)=="1"){
			    design<-model.matrix(~-1+factor(z))  
			 }
			 if (as.numeric(yy)=="2"){
		
			    design<-model.matrix(~-1+factor(w))
			 }
			 
			 if (as.numeric(yy)=="3"){ 
			    design<-model.matrix(~-1+factor(q))
	
			 }
			 r<-readline("Column names:")
			 r<-as.factor(unlist(strsplit(r, ",")))
			 colnames(design)<-tolower(as.vector(r))
			 }
		      
	        return(design)
		
	       }
    } else {
    
	  print ("The maximum number of factors for this fixed model are 3")
    }
  } else {
  
    print ("Random Effects")
    
     }
 }