# seqanalyseR
这个R包共分为八个函数，涵盖归一、差异分析、可视化、富集分析,下面是函数的案例  
setwd("C:/Users/colet/Desktop/nyq2022")  



group<-read.csv(file = "group.csv",header = T)  
# 下游分析全过程
DOWNSTREAManalysis(DATA = nyq2022,GROUP = group,SPECIES = "RAT",GENE_TYPE = "ENSEMBL")  

# 转换ID
aa<-IDchange(
  geneIDtoCHANGE = row.names(nyq2022),  
  CHANGEfromTYPE = "ENSEMBL", 
  CHANGEtoTYPE = "ENTREZID",  
  SPECIES = "RAT",  
  omitNA = TRUE,  
  MERGEtoDATA = nyq2022,  
  drop = TRUE 
) 
# 差异分析
aa<-differenceANALYSIS1(rawdata = nyq2022,group = group,omitNA = T)  
# GO富集（KEGG同理）
GO<-GOenrich(aa$ENTREZID, SPECIES = "RAT")  
# 火山图
plotVOCANOL(差异分析后数据,    
            LOG2FC阈值,   
            PVALUE阈值, 
            标记感兴趣的基因,     
            标题)
# 热图
plotHEATMAP1(差异分析后的数据, 分组文件2（行名为样品名，第一列为分组）)
