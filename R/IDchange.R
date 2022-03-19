#' IDchange
#'
#'Analyse RNA-seq data
#' @param DATA raw data
#' @param geneIDtoCHANGE gene list, may from rawdata
#' @param CHANGEfromTYPE type of gene
#' @param CHANGEtoTYPE type of gene
#' @param SPECIES specie
#' @param omitNA omitNA
#' @param MERGEtoDATA do you want to merge chanegd geneid to the rawdata
#' @example
#' IDchange(geneIDtoCHANGE, CHANGEfromTYPE, CHANGEtoTYPE, SPECIES,omitNA = TRUE, MERGEtoDATA = FALSE,drop = TRUE)
IDchange<-function (geneIDtoCHANGE, CHANGEfromTYPE, CHANGEtoTYPE, SPECIES,omitNA = TRUE, MERGEtoDATA = FALSE,drop = TRUE)
{ library(dplyr)
  species_list<-c("HUMAN","MOUSE","RAT")
  if (!SPECIES%in%species_list) {
    ERRORinSPECIES<-paste("input should be one of",paste(species_list,collapse = ","))
    stop('"SPECIES"', ERRORinSPECIES)
  }
  if (SPECIES=="HUMAN") {
    library("org.Hs.eg.db")
    OrgDb="org.Hs.eg.db"
  }
  if (SPECIES=="MOUSE") {
    library("org.Mn.eg.db")
    OrgDb="org.Mn.eg.db"
  }
  if (SPECIES=="RAT") {
    library("org.Rn.eg.db")
    OrgDb="org.Rn.eg.db"
  }
  db <- GOSemSim:::load_OrgDb(OrgDb)
  idTypes<-keytypes(db)
  if (!CHANGEfromTYPE %in% idTypes) {
    stop("'CHANGEfromTYPE' ", ERRORinTYPE)
  }
  if (!all(CHANGEtoTYPE %in% idTypes)) {
    stop("'CHANGEtoTYPE' ", ERRORinTYPE)
  }
  ERRORinTYPE <- paste0("should be one of ", paste(idTypes, collapse = ", "),
                        ".")
  geneID<-geneIDtoCHANGE
  geneID %<>% as.character %>% unique
  db <- GOSemSim:::load_OrgDb(OrgDb)
  res <- suppressWarnings(AnnotationDbi::select(db, keys = geneID,
                                                keytype = CHANGEfromTYPE, columns = c(CHANGEfromTYPE, CHANGEtoTYPE)))
  faied_to_change <- which(is.na(res[, 2]))
  if (length(faied_to_change)) {
    n <- res[faied_to_change, 1] %>% unique %>% length
    if (n) {
      warning(paste0(n/length(geneID) * 100,
                     "%"), " of input genes are failed to changeIDs...")
    }
    if (drop) {
      res <- res[-faied_to_change, ]
    }
  }
  if(omitNA==T){
    res<-na.omit(res)
  }
  RES<-res
  if(mode(MERGEtoDATA)=="list"){
    colnames(RES)[1]<-colnames(MERGEtoDATA)[1]
    RES<-merge(RES,MERGEtoDATA,by=colnames(MERGEtoDATA)[1])
    RES<-RES[,-1]
  }
  return(RES)
}
