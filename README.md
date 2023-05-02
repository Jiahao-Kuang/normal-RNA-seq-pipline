# seqanalyseR
# 示例文档
在data中添加了示例文件，count.txt为生数据，group.csv为分组文件  
data<-read.table(file = "count.txt",header = T)   
group<-read.csv(file = "group.csv",header = T)  
group2<-read.csv(file = "group.csv",header = T,row.names = 1)   
# 下游分析全过程函数
DOWNSTREAManalysis(DATA = data,GROUP = group,SPECIES = "RAT",GENE_TYPE = "ENSEMBL")  

# 转换ID
aa<-IDchange(
  geneIDtoCHANGE = row.names(data),  
  CHANGEfromTYPE = "ENSEMBL", 
  CHANGEtoTYPE = "ENTREZID",  
  SPECIES = "RAT",  
  omitNA = TRUE,  
  MERGEtoDATA = data,  
  drop = TRUE 
) 
# 差异分析
bb<-differenceANALYSIS1(rawdata = data,group = group,omitNA = T)  
# GO富集（KEGG同理）
GO<-GOenrich(aa$ENTREZID, SPECIES = "RAT")  
# 火山图
plotVOCANOL(bb,    
            LOG2FC = 0.05,   
            PVALUE =0.05, 
            TAQ = FALSE,#(可以提供感兴趣基因的列表)     
            TITAL = "火山图")
# 热图
plotHEATMAP1(bb, group2) #分组文件2（行名为样品名，第一列为分组）)
# 点图
f<- plotDOT1(GO)
f<- plotDOT2(GO)
# 条状图
g<- plotBAR(KEGG)
