#' plotHEATMAP1
#'
#'Analyse RNA-seq data
#' @param difference_analysised_DATA difference analysised data with differenceANALYSE
#' @param HEATMAP_COL_GROUP colgroup
#' @example
#' plotHEATMAP1(difference_analysised_DATA,HEATMAP_COL_GROUP)
plotHEATMAP1<-function(difference_analysised_DATA,HEATMAP_COL_GROUP){
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
  ggheatmap::ggheatmap(HEATMAP_MATRIX,cluster_rows = T,cluster_cols = T,
                       color = colorRampPalette(c("#2fa1dd", "white", "#f87669"))(100),
                       annotation_cols = HEATMAP_COL_GROUP,
                       annotation_color = col,
                       text_show_rows = F,
                       show_cluster_rows = F,
                       scale = "row")%>%
    ggheatmap_theme(1,theme =list(
      theme(axis.text.x = element_text(angle = 45, color = 'black', hjust = 0.8,vjust = 0.8, size = 10)))
    )->p_ggheatmap
  ggsave(file= "goodlooking_heatmap.png", plot = p_ggheatmap,dpi = DPI,width=WIDTH, height=HEIGHT)
}
