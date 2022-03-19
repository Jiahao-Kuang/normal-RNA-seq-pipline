#' plotHEATMAP2
#'
#'Analyse RNA-seq data
#' @param difference_analysised_DATA difference analysised data with differenceANALYSE
#' @param HEATMAP_COL_GROUP colgroup
#' @example
#' plotHEATMAP1(difference_analysised_DATA,HEATMAP_COL_GROUP)
plotHEATMAP2<-function(difference_analysised_DATA,HEATMAP_COL_GROUP){
  library(dplyr)
  DF=difference_analysised_DATA
  DF_up<-DF[DF$group=="Up",]
  DF_down<-DF[DF$group=="Down",]
  up_genename<-DF_up[,1]
  down_genename<-DF_down[,1]
  different_genename<-c(up_genename,down_genename)
  HEATMAP_MATRIX<-tocount_TPM[different_genename,]
  HEATMAP_MATRIX<-HEATMAP_MATRIX[,6:a]
  HEATMAP_MATRIX=t(scale(t(HEATMAP_MATRIX)))
  HEATMAP_COL_GROUP=HEATMAP_COL_GROUP
  colnames(HEATMAP_COL_GROUP)=c('group')
  annotation_color<-list(group=c(control='dodgerblue',treat='firebrick'))
  colnames(HEATMAP_MATRIX)<-rownames(HEATMAP_COL_GROUP)
  colnames(HEATMAP_COL_GROUP)<-c("annotation")
  col <- list(annotation=c(control = "#2fa1dd",treat = "#f87669"))
  library(aplot)
  library(tidyr)
  library(ggtree)
  p <- scale(HEATMAP_MATRIX) %>% data.frame()
  b<-ncol(p)
  #绘制行聚类树
  phr <- hclust(dist(p)) %>% ggtree(layout="rectangular", branch.length="none",color='white')
  #绘制列聚类树
  phc <- hclust(dist(t(p))) %>% ggtree() + layout_dendrogram()
  p$mtxars <- rownames(p)
  #宽表转长表
  p1 <- gather(p, 1:b, key="condition", value='pvalue')
  #绘制热图
  pp <- ggplot(p1,aes(condition,mtxars,fill=pvalue)) + geom_tile()+
    theme_tree() +
    scale_fill_viridis_c() +
    scale_y_discrete(position = POSITION)+xlab(NULL) + ylab(NULL)+
    theme(panel.background = element_blank(),#空白背景
          axis.title = element_text(size = 10, color = 'black'),#坐标轴标题字体
          axis.text.y = element_text(color = 'white'),#y轴字体
          axis.text.x = element_text(angle = 45,hjust = 0.8,vjust = 0.8 ,colour = 'black',size = 10),
          legend.position = 'right'
    )
  #利用aplot包将聚类树与热图拼接
  pp %>% insert_left(phr, width=.1) %>%
    insert_top(phc, height=.1)->p_ggplot2_heatmap
  ggsave(file= "ggplot2heatmap.png", plot = p_ggplot2_heatmap ,dpi = DPI,width=WIDTH, height=HEIGHT)
}
