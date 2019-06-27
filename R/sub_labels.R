#' Replace names with non-missing labels
#'
#' Replace names, esp. rownames, with non-missing labels.
#'
#' @param nms Vector of names, e.g. rownames.
#' @param labels Vector of labels, e.g. \code{labrows}.
#' @inheritParams ezheat
#' @return Vector of replaced rownames; vector's names are rownames

# currently unused
sub_labels <- function(nms, labels, na.lab=c("---", "")){
  ret <- stats::setNames(nms, nm=nms)
  na.sym.ind <- which(is.na(labels) | labels %in% na.lab)
  if (length(na.sym.ind) > 0){
    ret[-na.sym.ind] <- labels[-na.sym.ind]
  } else {
    ret <- labels
  }
  ret
}
