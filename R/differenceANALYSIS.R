#' differenceANALYSIS
#'
#'Analyse RNA-seq data
#' @param rawdata raw data
#' @param group experiment group
#' @param omitNA omitNAs
#' @example
#' differenceANALYSIS(rawdata,group,omitNA)
differenceANALYSIS<-function(rawdata,group,omitNA){
  a<-ncol(tocount_TPM)
  library(DESeq2)
  COUNTDATA<-(rawdata[,c(6:a)])
  COUNTDATA_1<-as.matrix(COUNTDATA)
  COUNTDATA<-COUNTDATA_1[rowMedians(COUNTDATA_1)>10,]
  as.matrix(COUNTDATA)
  colnames(COUNTDATA)==COLDATA$id#判断是否一致，-/.
  colnames(COUNTDATA)<-COLDATA$id
  dds<-DESeqDataSetFromMatrix(countData = COUNTDATA, colData = group,design = ~dex)
  dds<-DESeq(dds)
  res1<-results(dds)
  head(res1)
  res2<-data.frame(res1)
  if(omitNA==T){
    RESULT<-na.omit(res2)#去掉缺失值
  }else{RESULT=res2}
}
