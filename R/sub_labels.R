#' Replace names with non-missing label.v
#'
#' Replace names, esp. rownames, with non-missing label.v. \code{nms} and \code{label.v} must be of same length,
#' and \code{nms} must not have \code{NA}s.
#'
#' @param nms Vector of names, e.g. rownames.
#' @param label.v Vector of label.v, e.g. \code{labrows}.
#' @inheritParams ezheat
#' @return Vector of replaced rownames; vector's names are rownames

# currently unused
sub_labels <- function(nms, label.v, na.lab=c("---", "")){
  stopifnot(length(nms)==length(label.v), !is.na(nms), !(nms %in% na.lab))

  ret <- stats::setNames(nms, nm=nms)
  na.sym.ind <- which(is.na(label.v) | label.v %in% na.lab)
  if (length(na.sym.ind) > 0){
    ret[-na.sym.ind] <- label.v[-na.sym.ind]
  } else {
    ret <- label.v
  }
  ret
}
