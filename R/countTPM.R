#' countTPM
#'
#'Analyse RNA-seq data
#' @param DATA raw data
#' @param clean whether to clean the data
#' @param x what method would you like to use to clean your data
#'
#' @example
#' countTPM(rawdata,clean,x)
countTPM<-function(rawdata,clean,x){
  a<-ncol(tocount_TPM)
  tpm1<-(tocount_TPM[,6:a]/tocount_TPM[,5])/colSums(tocount_TPM[,6:a]/tocount_TPM[,5])*1e6
  tpm1
  if(clean=="MEDIAN"){
    tpm1<-as.matrix(tpm1)
    TPM_medianed<-tpm1[rowMedians(tpm1)>x,]
  }else if(clean=="SIMPLE"){
    TPM_simple<-tpm1[rowSums(tpm1)!=x,]}
  else{TPM=tpm1}
}
