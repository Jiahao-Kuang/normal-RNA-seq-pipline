#' plotBAR
#'
#'Analyse RNA-seq data
#' @param x GO or KEGG enrich analized data
#' @param enrichTYPE GO or KEGG enrich type
#' @example
#' plotBAR(x,enrichTYPE)
plotBAR<-function(x,enrichTYPE){
  type_list<-c("GO","KEGG")
  if (!enrichTYPE%in%species_list) {
    ERRORinSPECIES<-paste("input should be one of",paste(type_list,collapse = ","))
    stop('"enrichTYPE"', ERRORinSPECIES)
  }
  if (enrichTYPE=="GO") {
    GO_2<-x
    GO_3<-GO_2[,c(4,7,11)]
    GO_4<--log10(GO_3$pvalue)
    GO_5<-GO_3[,c(1,3)]
    GO_6<-cbind(GO_5,GO_4)
    colnames(GO_6)<-c('GO_terms','Gene_counts','-log10(Pvalue)')
    GO_7<-GO_6[order(GO_6$`-log10(Pvalue)`,decreasing = TRUE),]
    GO_8<-GO_7[c(1:20),]#取最显著的20个通路
    ggplot(GO_8)+
      geom_col(aes(x=reorder(GO_terms,`-log10(Pvalue)`),y=Gene_counts),
               color='black',#柱边框颜色
               width= 0.6,#柱宽
               fill='#d5a478')+
      labs(x=NULL,y='Gene count', #自定义x、y轴、标题内容
           title = 'Enriched GO Terms (expanded gene families)')+
      theme_test(base_size = 15)+ #主题基本大小
      theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 8),
            axis.text = element_text(color = 'black',face = 'bold'),
            #plot.margin = margin(1,0.5,0.5,2.5,'cm'),
            panel.border = element_rect(size = 1),
            axis.title = element_text(face = 'bold'),
            plot.title = element_text(face = 'bold',
                                      size=13,hjust = 0.5))->p_A
    p_B <- p_A+
      geom_text(aes(x=reorder(GO_terms,`-log10(Pvalue)`),
                    y=Gene_counts,
                    label=Gene_counts),
                vjust=-0.5,size=3.5,fontface='bold')
    p_C <- p_B+
      scale_y_continuous(sec.axis = sec_axis(~./42,
                                             name = '-log10(Pvalue)',
      ))+
      geom_line(aes(x= GO_terms,
                    y=`-log10(Pvalue)`*42,
                    group=1),
                linetype=3,cex=1)+
      geom_point(aes(x= GO_terms,
                     y=`-log10(Pvalue)`*42),
                 color='#589c47',size=3.5)
    return(p_C)
  }
  if (enrichTYPE=="KEGG") {
  KEGG_2<-x
  KEGG_3<-KEGG_2[,c(3,6,10)]
  KEGG_4<--log10(KEGG_3$pvalue)
  KEGG_5<-KEGG_3[,c(1,3)]
  KEGG_6<-cbind(KEGG_5,KEGG_4)
  colnames(KEGG_6)<-c('KEGG_terms','Gene_counts','-log10(Pvalue)')
  KEGG_7<-KEGG_6[order(KEGG_6$`-log10(Pvalue)`,decreasing = TRUE),]
  KEGG_8<-KEGG_7[c(1:20),]#取最显著的十个通路
  ggplot(KEGG_8)+
    geom_col(aes(x= reorder(KEGG_terms,`-log10(Pvalue)`) ,y=Gene_counts),
             color='black',#柱边框颜色
             width= 0.6,#柱宽
             fill='#d5a478')+
    labs(x=NULL,y='Gene count', #自定义x、y轴、标题内容
         title = 'Enriched KEGG Terms (expanded gene families)')+
    theme_test(base_size = 15)+ #主题基本大小
    theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 8),
          axis.text = element_text(color = 'black',face = 'bold'),
          #plot.margin = margin(1,0.5,0.5,2.5,'cm'),
          panel.border = element_rect(size = 1),
          axis.title = element_text(face = 'bold'),
          plot.title = element_text(face = 'bold',
                                    size=13,hjust = 0.5))->p_D
  p_E <- p_D+
    geom_text(aes(x= reorder(KEGG_terms,`-log10(Pvalue)`),
                  y=Gene_counts,
                  label=Gene_counts),
              vjust=-0.5,size=3.5,fontface='bold')
  p_F <- p_E+
    scale_y_continuous(sec.axis = sec_axis(~./42,
                                           name = '-log10(Pvalue)',
    ))+
    geom_line(aes(x= KEGG_terms,
                  y=`-log10(Pvalue)`*42,
                  group=1),
              linetype=3,cex=1)+
    geom_point(aes(x= KEGG_terms,
                   y=`-log10(Pvalue)`*42),
               color='#589c47',size=3.5)
  }
}
