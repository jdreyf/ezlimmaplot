#' Extract prefix from names
#'
#' Extract prefixes from names, presumably colnames.
#'
#' @param nm Character vector of names to search.
#' @param suffix Character string to search for at end of \code{cnames}.
#' @param sep Separator in \code{cnames} preceeding \code{suffix}. If \code{NA} or absent from \code{nm}, it is ignored.
#' @details If \code{sep} is ignored, elements of \code{nm} must match \code{suffix} exactly.
#' @return A character vector of prefixes

extract_prefix <- function(nm, suffix, sep='.'){
  if (!is.na(sep)){
    sep.cols <- grep(pattern=paste0('\\', sep), x=nm)
    if (length(sep.cols)==0) sep <- NA
  }

  #recheck is.na(sep), since may have assigned it to NA above
  if (!is.na(sep)){
    patt <- paste0('\\', sep, '(', suffix, ')$')
    ind <- grep(pattern=patt, x=nm)
    prefix.v <- sub(patt, '', nm[ind])
  } else {
    patt <- paste0("^", suffix, '$')
    ind <- grep(pattern=patt, x=nm)
    #could return string of multiple ""
    prefix.v <- sub(patt, '', nm[ind])
  }
  if (length(prefix.v)==0) stop("No prefix found.")
  return(prefix.v)
}
