#' Prune matrix
#'
#' Prune (subset) rows of matrix for plotting as heatmap.
#'
#' @inheritParams ezheat
#' @return Vector of pruned matrix rownames; vectors names are from rownames(object).

prune_mat <- function(object, labrows=rownames(object), only.labrows=FALSE, unique.rows=FALSE, ntop=NULL,
                      na.lab=c("---", ""), verbose=TRUE){
  if (!is.matrix(object)) object <- as.matrix(object)
  stopifnot(length(labrows)==nrow(object), names(labrows)==rownames(object))

  ret <- stats::setNames(labrows, nm=rownames(object))

  if (only.labrows){
    na.sym.ind <- which(is.na(labrows) | labrows %in% na.lab)
    if (length(na.sym.ind) > 0){
      if (verbose) message("Removing ", length(na.sym.ind), " rows without gene labrows.")
      ret <- ret[-na.sym.ind]
      object <- object[-na.sym.ind,]
    }
  }
  if (unique.rows){
    dup <- duplicated(ret)
    if (verbose) message("Removing ", sum(dup), " rows with duplicated names.")
    object <- object[!dup]
    ret <- ret[!dup]
  }
  if (!is.null(ntop)){
    if (verbose) message("Selecting top ", ntop, " rows.")
    if (ntop <= length(ret)){
      ret <- ret[1:ntop]
    } else {
      warning("After processing, object has only ", length(ret), " rows, so cannot subset to ", ntop, " rows.")
    }
  }
  return(ret)
}
