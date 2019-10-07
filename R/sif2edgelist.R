#' Clean Pathway Commons SIF file
#'
#' Remove some ubiquitious IDs from Pathway Commons SIF file, and transform to edge-list.
#'
#' @param sif SIF file from Pathway Commons
#' @param rm.ids vector of analyte IDs to remove. Default \code{CHEBI:15377} corresponds to water.
#' @return Edge-list matrix with two columns.
#' @export

sif2edgelist <- function(sif, rm.ids="CHEBI:15377"){
  stopifnot(ncol(sif) >= 3)
  #yields 2 conjugate acids that should be switched
  #rm rows with CHEBI:15377=water
  rm.rows <- which(apply(sif, MARGIN=1, FUN=function(v){
    length(grep(paste(paste0("^", rm.ids, "$"), collapse="|"), v)) > 0
  }))
  if (length(rm.rows) > 0) sif <- sif[-rm.rows,,drop=FALSE]
  #some interactions are A->B and duplicated as B->A
  sif <- as.matrix(sif[,c(1,3)],drop=FALSE)
  sif <- t(apply(sif, MARGIN=1, FUN=sort))
  #don't use which() in case none are duplicated
  sif <- sif[!duplicated(sif),,drop=FALSE]
  return(sif)
}
