#' plotDOT1
#'
#'Analyse RNA-seq data
#' @param x GO or KEGG enrich analized data
#' @example
#' plotDOT1(x)
plotDOT1<-function(x){
  library(clusterProfiler)
  dotplot(GO,split='ONTOLOGY')+facet_grid(ONTOLOGY~.,scales = 'free')
}
