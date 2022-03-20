#' plotDOT2
#'
#'Analyse RNA-seq data
#' @param x GO or KEGG enrich analized data
#' @example
#' plotDOT2(x)
plotDOT2<-function(x){
  library(ggplot2)
  library(ggpubr)
  matrix<-matrix[order(matrix$GeneCounts,decreasing = T),]
  matrix$row<-c(1:25)
  ggplot(matrix, aes(x=GeneCounts, y=reorder(GO_Terms,GeneCounts), fill=row)) +
    geom_dotplot(binaxis='y', stackdir='center')+
    labs(title="GO Terms",x="Number of Genes", y = "")->p1
  p1+
    #主题
    theme_bw()+
    theme(
      panel.border = element_rect(linetype = 1, color="black",size=2),
      axis.ticks = element_line(size = 1,colour = 'black'),#坐标轴格子
      axis.title = element_text(size = 12, color = 'black',face = "bold"),#坐标轴标题字体
      axis.text.x = element_text(size = 10,color = 'black',face = "bold"),#坐标轴字体
      axis.text.y = element_text(size = 10,color = a,face = "bold"),
      plot.title = element_text(size = 15 , hjust=0.5,face = "bold"),#图表标题
      legend.position = "none"
    )->p2
  p2+
    scale_fill_gradientn(values = seq(0,1,0.001),colours = c("#EDEF5C","#255668"))->p3
  return(p3)
}