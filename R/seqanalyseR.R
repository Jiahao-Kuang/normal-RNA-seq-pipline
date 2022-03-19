#' seqanalyseR : A package for RNA-seq
#'
#' A package to provide RNA-seq downstream analysis
#'
#' @section Functions :
#' RNA-seq 下游分析
#' 用help函数查看下列每一个函数的帮助
#' @section 1.countTPM :算TPM值，输入RAWDATA
#' @section 2.differenceALALYSIS :调用DESeq做差异分析
#' @section 做富集分析之前要转换ID成ENTRIZID
#' @section 3.IDchange :可以改变gene ID的类型，可以给列表，也可以直接把RAWDATA的gene ID列放进来，如果有需要可以设置merge到rawdata和去除未匹配的值
#' @section 4.GOenrich :GO富集
#' @section 5.KEGGenrich :KEGG富集
#' @section 上述富集分析默认会在工作目录下直接输出可视化结果，不用指定输出对象
#' @section 6.plotVOCANOL :画火山图
#' @section 7.plotHEATMAP1 :画热图，比较符合审美
#' @section 8.plotHEATMAP2 :画热图，更加符合审美
#' @section 9.DOWNstreamANALYSIS 禁忌的函数，一点不对都会报错，测试过程中画出的图都很诡异
#' 但是直接用它可以走完整个下游流程，不需要指定输出，运气好可以直接看到图片在工作目录生成
#'
#' @encoding UTF-8
#' @docType pacakge
#' @name seqanalyseR

NULL
