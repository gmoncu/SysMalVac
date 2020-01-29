bothQualityControl<-function(dataCEL,dataNORM,intgroup){

    preQualityControl(data=dataCEL,intgroup)
    postQualityControl(data=dataNORM,intgroup)

}