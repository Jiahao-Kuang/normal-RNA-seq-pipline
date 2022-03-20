#' DOWNSTREAManalysis
#'
#'Analyse RNA-seq data
#' @param DATA raw data
#' @param GROUP experiment group
#' @param GENE_TYPE gene type annotated in raw data
#' @param SPECIES species
#' @examples
#' read.csv
#' DOWNSTREAManalysis(DATA = "DATA",
#'                    GROUP = "group",
#'                    GENE_TYPE = "ENSEMBL",
#'                    SPECIES = "RAT")
DOWNSTREAManalysis<-function(DATA, GROUP, GENE_TYPE, SPECIES){
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(DOSE)
  library(ggheatmap)
  library("clusterProfiler")
  library("org.Hs.eg.db")
  library("org.Mm.eg.db")
  library("org.Rn.eg.db")
  #输入的基因命名类型（ENSEMBL/SYMBOL）
  #若counts文件的基因名为ENSEMBL_id，则需要id转换，选错后面运行不了
  #导出图片的设定
  DPI = 800
  WIDTH = 6
  HEIGHT = 6
  #火山图设定
  #火山图阈值调整
  LOG2FC <- 1
  PVALUE <- 0.05
  #火山图阈值线调整
  XLINE<-c(-1,1)
  YLINE<-c(0.05)
  #热图聚类显示的方向
  POSITION<-'right'
  #设定GO富集分析的数据库
  species_list<-c("HUMAN","MOUSE","RAT")
  if (!SPECIES%in%species_list) {
    ERRORinSPECIES<-paste("input should be one of",paste(species_list,collapse = ","))
    stop('"SPECIES"', ERRORinSPECIES)
  }
  if (SPECIES=="HUMAN") {
    OrgDb="org.Hs.eg.db"
    KEGG_database="hsa"
  }
  if (SPECIES=="MOUSE") {
    OrgDb="org.Mn.eg.db"
    KEGG_database="mmu"
  }
  if (SPECIES=="RAT") {
    OrgDb="org.Rn.eg.db"
    KEGG_database="rno"
  }
  #一、归一化
  a<-ncol(tocount_TPM)
  tpm1<-(tocount_TPM[,6:a]/tocount_TPM[,5])/colSums(tocount_TPM[,6:a]/tocount_TPM[,5])*1e6
  tpm1<-as.matrix(tpm1)
  tpm1[rowMedians(tpm1)>10,]
  TPM<-tpm1[rowSums(tpm1)!=0,]
  #三、差异分析
  COUNTDATA<-(tocount_TPM[,c(6:a)])
  COUNTDATA_1<-as.matrix(COUNTDATA)
  COUNTDATA<-COUNTDATA_1[rowMedians(COUNTDATA_1)>10,]
  as.matrix(COUNTDATA)
  COLDATA<-GROUP
  colnames(COUNTDATA)==COLDATA$id#判断是否一致，-/.
  colnames(COUNTDATA)<-COLDATA$id
  dds<-DESeqDataSetFromMatrix(countData = COUNTDATA, colData = COLDATA,design = ~dex)
  dds<-DESeq(dds)
  res1<-results(dds)
  head(res1)
  res2<-data.frame(res1)
  res2 %>% #将数据集传递给后面的函数mutate,mutate的作用为给数据集增加新的列
    mutate(group = case_when(
      log2FoldChange >= 1 & padj <= 0.05 ~ "UP" ,
      log2FoldChange <=-1 & padj <= 0.05 ~ "DOWN" ,
      TRUE ~ "NOT_CHANGE"
    )) ->RESULT
  RESULT<-na.omit(res2)#去掉缺失值
  RESULT<-na.omit(RESULT)#去掉缺失值
  write.csv(RESULT,file = "different_express.csv",quote = F)
  diff_RESULT_output<-as.data.frame(RESULT)
  diff_RESULT_output$ensembl_gene_id<-gsub("\\.\\d*","",row.names(RESULT))
  OUTPUT_DIFF<-merge(my_result,diff_RESULT_output,by= intersect(names(my_result),names(diff_RESULT_output)))
  OUTPUT_DIFF$ensembl_gene_id
  TPM$ensembl_gene_id<-rownames(TPM)
  output<-merge(OUTPUT_DIFF,TPM,by='ensembl_gene_id')
  write.csv(output,file = "diff_express_gene.csv",row.names = F)
  #五、绘制火山图
  #绘制火山图前
  DF<-cbind(gene_name=row.names(res2),res2)
  row.names(DF)=NULL
  DF$group<-ifelse(DF$log2FoldChange>=LOG2FC& DF$pvalue<=PVALUE,"Up",
                   ifelse(DF$log2FoldChange<=-LOG2FC&DF$pvalue<=PVALUE,"Down","Not sig"))
  DF<-na.omit(DF)#去掉缺失值
  DF$dramatically_sig <- (-log10(DF$pvalue)) #挑选极为显著的基因
  DRAMATICALLY_sig<-DF[DF$dramatically_sig>=10,]
  DF_sig<-DF[-grep("Not sig",DF$group),]
  #火山图
  ggplot(data = DF,aes(x=log2FoldChange,y=-log10(pvalue)))+
    geom_point() ->p
  #火山图美化
  p+
    geom_point(aes(color=group))+
    scale_color_manual(values = c("#2fa1dd", "grey", "#f87669"))+
    #坐标轴
    scale_x_continuous(limits = c(-6,6), breaks = c(-6,-4,-2,0,2,4,6))+
    scale_x_continuous(limits = c(-6,6), breaks = c(-6,-4,-2,0,2,4,6))+
    #阈值线
    geom_vline(xintercept = XLINE, linetype = 5,size=1)+
    geom_hline(yintercept=-log10(YLINE), linetype=5,size=1) +
    #主题
    theme_bw()+
    theme(
      panel.border = element_rect(linetype = 1, color="black",size=2),
      axis.ticks = element_line(size = 1,colour = 'black'),#坐标轴格子
      axis.title = element_text(size = 15, color = 'black'),#坐标轴标题字体
      axis.text = element_text(size = 10,color = 'black')#坐标轴字体
    ) ->p1
  ggsave(file= "vocanol.png", plot = p1,dpi = DPI,width=WIDTH, height=HEIGHT)
  #六、绘制热图
  #热图准备
  #挑出上调和下调的基因名做成矩阵
  DF_up<-DF[DF$group=="Up",]
  DF_down<-DF[DF$group=="Down",]
  up_genename<-DF_up[,1]
  down_genename<-DF_down[,1]
  different_genename<-c(up_genename,down_genename)
  HEATMAP_MATRIX<-tocount_TPM[different_genename,]
  HEATMAP_MATRIX<-HEATMAP_MATRIX[,6:a]
  HEATMAP_MATRIX=t(scale(t(HEATMAP_MATRIX)))
  HEATMAP_COL_GROUP<-read.csv("分组文件.CSV",header = T,row.names = 1)
  colnames(HEATMAP_COL_GROUP)=c('group')
  #确定分组颜色
  annotation_color<-list(group=c(control='dodgerblue',treat='firebrick'))
  #行列名统一
  colnames(HEATMAP_MATRIX)<-rownames(HEATMAP_COL_GROUP)
  #开始绘制热图
  #热图(pheatmap)
  pheatmap(HEATMAP_MATRIX,
           annotation_col = HEATMAP_COL_GROUP,
           color = colorRampPalette(c("#003366", "white"))(50),
           cutree_rows = 2,
           cluster_rows = F,
           show_rownames = F,
           annotation_colors = annotation_color
  ) ->p_pheatmap
  ggsave(file= "heatmap.png", plot = p_pheatmap ,dpi = DPI,width=WIDTH, height=HEIGHT)
  #六、富集分析
  #转换id
  GENE_TO_ENRICH <- bitr(DF$gene_name, fromType = GENE_TYPE,
                         toType = c("ENTREZID"),
                         OrgDb = OrgDb)
  #GO富集
  GO<-enrichGO(GENE_TO_ENRICH$ENTREZID,
               OrgDb = OrgDb,
               keyType = 'ENTREZID',
               ont = 'ALL',
               pvalueCutoff = 0.05,#设定p值阈值
               qvalueCutoff = 0.05,#设定q值阈值
               readable = TRUE
  )
  #GO可视化
  barplot(GO, showCategory = 10)->p_GO_barplot#点图
  write.csv(GO,file = 'GO.csv')
  barplot(GO,split='ONTOLOGY')+facet_grid(ONTOLOGY~.,scales = 'free')->p_GO_barplot_distributed#柱状图(分组)
  dotplot(GO, showCategory = 10)->p_GO_dotplot#点图
  dotplot(GO,split='ONTOLOGY')+facet_grid(ONTOLOGY~.,scales = 'free')->p_GO_dotplot_distributed#点图（分组）
  enrichplot::cnetplot(GO, circular= FALSE, colorEdge= TRUE)-> p_GO_pathway#基因-相关通路网络图
  enrichplot::cnetplot(GO, circular= TRUE, colorEdge= TRUE)-> p_GO_circled_pathway#基因-相关通路网络图
  #GO可视化（2）
  GO_2<-read.csv('GO.csv',header = T)
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
  #KEGG富集
  KEGG<-enrichKEGG(GENE_TO_ENRICH$ENTREZID,
                   organism = KEGG_database,
                   pvalueCutoff = 0.05,#设定p值阈值
                   qvalueCutoff = 0.05,#设定q值阈值
  )
  #转换ID成可读ID
  KEGG_readable<-setReadable(KEGG,
                             OrgDb = OrgDb,
                             keyType = "ENTREZID")
  write.csv(KEGG_readable,'KEGG.csv')
  #KEGG可视化
  barplot(KEGG, showCategory = 10)->p_KEGG_barplot#柱状图
  dotplot(KEGG, showCategory = 10)->p_KEGG_dotplot#点图
  enrichplot::cnetplot(KEGG, circular= FALSE, colorEdge= TRUE)-> p_KEGG_pathway#基因-相关通路网络图
  enrichplot::cnetplot(KEGG, circular= TRUE, colorEdge= TRUE)-> p_KEGG_circled_pathway#基因-相关通路网络图
  #KEGG(可视化2)
  KEGG_2<-read.csv('KEGG.csv',header = T)
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
  ggsave(file= "goodlooking_GO.png", plot = p_C,dpi = DPI,width=WIDTH, height=HEIGHT)
  ggsave(file= "goodlooking_KEGG.png", plot = p_F,dpi = DPI,width=WIDTH, height=HEIGHT)
  #导出GO、KEGG富集分析结果
  ggsave(file= "GO_barplot.png", plot = p_GO_barplot,dpi = DPI,width=WIDTH, height=HEIGHT)
  ggsave(file= "GO_barplot_distributed.png", plot = p_GO_barplot_distributed,dpi = DPI,width=8, height=12)
  ggsave(file= "GO_dotplot.png", plot = p_GO_dotplot,dpi = DPI,width=WIDTH, height=HEIGHT)
  ggsave(file= "GO_dotplot_distributed.png", plot = p_GO_dotplot_distributed,dpi = DPI,width=8, height=15)
  ggsave(file= "GO_pathway.png", plot = p_GO_pathway,dpi = DPI,width=12, height=12)
  ggsave(file= "GO_circled_pathway.png", plot = p_GO_circled_pathway,dpi = DPI,width=24, height=24)
  ggsave(file= "KEGG_barplot.png", plot = p_KEGG_barplot,dpi = DPI,width=WIDTH, height=HEIGHT)
  ggsave(file= "KEGG_dotplot.png", plot = p_KEGG_dotplot,dpi = DPI,width=WIDTH, height=HEIGHT)
  ggsave(file= "KEGG_pathway.png", plot = p_KEGG_pathway,dpi = DPI,width=12, height=12)
  ggsave(file= "KEGG_pathway_circled.png", plot = p_KEGG_circled_pathway,dpi = DPI,width=24, height=24)
}
