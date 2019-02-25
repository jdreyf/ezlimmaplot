#' Extract prefix from names
#'
#' Extract prefixes from names, presumably colnames. For example, the prefix of \code{"x.p"} would be \code{"x"}.
#'
#' @param nm Character vector of names to search.
#' @param suffix Character string to search for at end of \code{nm}.
#' @param sep Separator in \code{nm} preceeding \code{suffix}. One of \code{".", "_", NA}, If \code{NA}, elements
#' of \code{nm} must match \code{suffix} exactly.
#' @return A character vector of prefixes. If there is no prefix and no sep, \code{NA} is returned.

extract_prefix <- function(nm, suffix="p", sep="."){
  stopifnot(sep %in% c(".", "_") || is.na(sep))

  if (length(grep(paste0("\\", sep), nm)) == 0){
    sep <- NA
    message("sep not found, so assigned NA.")
  }

  if (!is.na(sep)){
    patt <- paste0("\\", sep, "(", suffix, ")$")
    ind <- grep(pattern=patt, x=nm)
    prefix.v <- sub(patt, "", nm[ind])
  } else {
    patt <- paste0("^", suffix, "$")
    ind <- grep(pattern=patt, x=nm)
    #could return string of NAs
    #NA instead of "" distinguishes "p" from ".p"
    prefix.v <- sub(patt, NA, nm[ind])
  }
  if (length(prefix.v)==0) stop("Suffix ", suffix, " not found.")
  return(prefix.v)
}
