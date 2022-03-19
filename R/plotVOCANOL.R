#' plotVOCANOL
#'
#'Analyse RNA-seq data
#' @param differenceDATA analised data with differenceANALYSE
#' @param LOG2FC logFC cutoff
#' @param PVALUE Pvalue cutoff
#' @param TAQ add taqs in your interested genes
#' @param TITLE the title of the plot
#' @example
#' plotVOCANOL(differenceDATA,LOG2FC,PVALUE,TAQ,TITLE)
plotVOCANOL<-function(differenceDATA,LOG2FC,PVALUE,TAQ,TITLE){
  DF_1<-cbind(gene_name=row.names(differenceDATA),differenceDATA)
  row.names(DF_1)=NULL
  DF_1$group<-ifelse(DF_1$log2FoldChange>=LOG2FC& DF_1$pvalue<=PVALUE,"Up",
                     ifelse(DF_1$log2FoldChange<=-LOG2FC&DF_1$pvalue<=PVALUE,"Down","Not sig"))
  DF_1<-na.omit(DF_1)
  DF_1_sig<-DF_1[-grep("Not sig",DF_1$group),]
  library(ggplot2)
  ggplot(data = DF_1,aes(x=log2FoldChange,y=-log10(pvalue)))+
    geom_point() ->p
  p+
    geom_point(aes(color=group))+
    scale_color_manual(values = c("#2fa1dd", "grey", "#f87669"))+
    scale_x_continuous(limits = c(-6,6), breaks = c(-6,-4,-2,0,2,4,6))+
    labs(x='logFC',y='-log10(pvalue)',title = TITLE)+
    scale_x_continuous(limits = c(-6,6), breaks = c(-6,-4,-2,0,2,4,6))+
    geom_vline(xintercept = c(LOG2FC,-LOG2FC), linetype = 5,size=1)+
    geom_hline(yintercept=PVALUE, linetype=5,size=1) +
    theme_bw()+
    theme(
      panel.border = element_rect(linetype = 1, color="black",size=2),
      axis.ticks = element_line(size = 1,colour = 'black'),
      axis.title = element_text(size = 15, color = 'black',face = 'bold'),
      axis.text = element_text(size = 10,color = 'black',face = 'bold'),
      plot.title = element_text(size = 15 , hjust=0.5,face = 'bold'),
      legend.title = element_text(size = 10,color = 'black',face = 'bold'),
      legend.text = element_text(size = 10,color = 'black',face = 'bold')
    ) ->p1
  if(mode(TAQ)=="character"){
    TAQ<-as.data.frame(TAQ)
    colnames(TAQ)<-"gene_name"
    TAQ<-merge(DF_1,TAQ,by="gene_name")
    library(ggrepel)
    p1+
      geom_label_repel(data=TAQ,aes(x=log2FoldChange, y=-log10(pvalue),label=gene_name))->p2
    p2
  }else{
    p1->p2
    p2
  }
  ggsave(file= "vocanol.png", plot = p2,dpi = 800,width=6, height=6)
}
