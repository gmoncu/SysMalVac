 postQualityControl<-function(data,intgroup){
if (is.null(intgroup)){
  arrayQualityMetrics(expressionset=data,outdir=paste(getwd(),"QC_post",sep="/"),force=TRUE,do.logtransform=TRUE)
  }else{
  arrayQualityMetrics(expressionset=data,outdir=paste(getwd(),"QC_post",sep="/"),force=TRUE,do.logtransform=TRUE, intgroup=intgroup)
  }
pathOUT<-paste(paste(getwd(),"QC_post",sep="/"),"index.html",sep="/")
if (sessionInfo()$platform=="x86_64-pc-linux-gnu (64-bit)"){
  system(paste('gnome-open',pathOUT,sep=" "))
  }else{
  print ("Open html manually")
  }
}