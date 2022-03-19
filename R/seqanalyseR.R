#' seqanalyseR : A package for RNA-seq
#'
#' A package to provide RNA-seq downstream analysis
#'
#' @section Functions :
#' 这个R包共分为八个函数，涵盖归一、差异分析、可视化、富集分析,下面是函数的案例
#' nyq2022<-read.table（file = "nyq2022.txt",header = T,row.names = 1
#' group<-read.csv(file = "group.csv",header = T)
#' @section 1.转换ID
#' aa<-IDchange( geneIDtoCHANGE = row.names(nyq2022),
#' CHANGEfromTYPE = "ENSEMBL", CHANGEtoTYPE = "ENTREZID",
#' SPECIES = "RAT",
#' omitNA = TRUE,
#' MERGEtoDATA = nyq2022,
#' drop = TRUE )
#' @section 2.differenceALALYSIS :调用DESeq做差异分析
#' @section 做富集分析之前要转换ID成ENTRIZID
#' @section 3.IDchange :可以改变gene ID的类型，可以给列表，也可以直接把RAWDATA的gene ID列放进来，如果有需要可以设置merge到rawdata和去除未匹配的值
#' @section 4.GOenrich :GO富集
#' @section 5.KEGGenrich :KEGG富集
#' @section 上述富集分析默认会在工作目录下直接输出可视化结果，不用指定输出对象
#' @section 6.plotVOCANOL :画火山图
#' @section 7.plotHEATMAP1 :画热图，比较符合审美
#' @section 8.plotHEATMAP2 :画热图，更加符合审美
#' @section 9.DOWNstreamANALYSIS : 走完整个下游分析
#'
#'
#' @encoding UTF-8
#' @docType pacakge
#' @name seqanalyseR

NULL
