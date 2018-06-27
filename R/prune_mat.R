#' Prune matrix
#'
#' Prune (subset) matrix before plotting as heatmap.
#'
#' @param object Matrix-like object.
#' @param labrows Labels corresponding to rows, e.g. gene symbols.
#' @param only.labrows Logical, if \code{TRUE} only include rows with a non-missing symbol.
#' @param unique.rows Logical, indicating if to remove duplicated row labels to make rows unique.
#' @param ntop Scalar number of rows to include.
#' @param na.lab Character vector of labels in \code{lab.col} to treat as missing, in addition to \code{NA}.
#' @param verbose Logical indicating whether to print messages to console.
#' @return Pruned matrix.

prune_mat <- function(object, labrows=NULL, only.labrows=FALSE, unique.rows=FALSE, ntop=NULL, na.lab=c('---', ''), verbose=TRUE){
  if (!is.matrix(object)) object <- data.matrix(object)
  if (!is.null(labrows)){
    stopifnot(length(labrows)==nrow(object), names(labrows)==rownames(object))
    na.sym.ind <- which(is.na(labrows) | labrows %in% na.lab)
    if (length(na.sym.ind) > 0){
      rownames(object)[-na.sym.ind] <- labrows[-na.sym.ind]
    } else {
      rownames(object) <- labrows
    }

    if (only.labrows){
      if (verbose) message('Removing ', length(na.sym.ind), ' rows without gene labrows.')
      object <- object[-na.sym.ind,]
    }
    if (unique.rows){
      if (verbose) message('Removing ', sum(duplicated(rownames(object))), ' rows with duplicated names.')
      object <- object[!duplicated(rownames(object)),]
    }
  }#end if !is.null(sym)

  if (!is.null(ntop)){
    if (verbose) message('Selecting top ', ntop, ' rows.')
    if (ntop <= nrow(object)){
      object <- object[1:ntop,]
    } else {
      if (verbose) message("After processing, object has only ", nrow(object), " rows, so cannot subset to ", ntop, " rows.")
    }
  }
  return(object)
}
