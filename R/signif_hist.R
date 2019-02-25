#' Plot histograms of significance (p-value & FDR) columns
#'
#' Plot histograms of significance (p-value & FDR) columns.
#'
#' @param p.suffix Suffix for p-value columns. P-value column names cannot be duplicated.
#' @param fdr.suffix Suffix for FDR columns. Set to \code{NA} if no FDR columns. FDR column names cannot be duplicated.
#' @param sep Separator for column names before \code{p} or \code{FDR}, passed to \code{\link[ezlimma]{extract_prefix}}.
#' If not found, it is assumed to be \code{NA}.
#' @param pi0 Logical indicating if proportion of null hypotheses should be calculated per p-value histogram. If
#' \code{TRUE}, \code{\link[limma]{propTrueNull}} with \code{method="convest"} is used.
#' @inheritParams ezheat
#' @inheritParams ezvenn
#' @details Some p-value columns must be identifiable using \code{p.suffix}. If \code{!is.na(fdr.suffix)}, FDR
#' colnames must have same prefix.
#' @return Invisibly, a subset of \code{tab} with only columns that contain significances.
#' @export

#assume each comparison has a p-value & q-value column unless fdr.suffix=NA
#could allow for no prefix, ie colnames(tab)=c("p", "FDR")
signif_hist <- function(tab, p.suffix="p", fdr.suffix="FDR", sep=".", pi0 = FALSE, name="signif_hist", plot=TRUE){
  stopifnot(nrow(tab) > 0, ncol(tab) > 0, !is.null(colnames(tab)))

  prefix.v <- extract_prefix(colnames(tab), suffix=p.suffix, sep=sep)
  if (any(duplicated(prefix.v))) stop("p-value column names are duplicated.")
  if (is.na(prefix.v[1])){
    p.cols <- match(p.suffix, colnames(tab))
  } else {
    p.cols <- match(paste0(prefix.v, sep, p.suffix), colnames(tab))
  }

  if (!is.na(fdr.suffix)){
    if (is.na(prefix.v[1])){
      fdr.cols <- match(fdr.suffix, colnames(tab))
    } else {
      fdr.cols <- match(paste0(prefix.v, sep, fdr.suffix), colnames(tab))
    }
    if (any(is.na(fdr.cols))) stop("!is.na(fdr.suffix) but FDR columns not found.")
    tab.ss <- tab[,sort(c(p.cols, fdr.cols))]
  } else {
    tab.ss <- tab[,p.cols]
  }

  #set name=NA to turn off pdf
  if (plot){
    if (!is.na(name)){ grDevices::pdf(paste0(name, ".pdf")) }
    graphics::par(mfrow=c(2,2))
    for (ind.tmp in 1:length(p.cols)){
      prefix <- prefix.v[ind.tmp]
      p.col <- p.cols[ind.tmp]
      stopifnot(length(p.col)==1)
      subtitle <- NULL
      if(pi0) {
        if (!requireNamespace("limma", quietly = TRUE)){
          stop("Package 'limma' needed to estimate pi0. Please install it.", call. = FALSE)
        }
        prop.null <- limma::propTrueNull(tab[,p.col], method = "convest")
        subtitle <- paste("Proportion of True Null = ", signif(prop.null, 3))
      }
      graphics::hist(tab[,p.col], xlab="P-value", main=prefix, sub = subtitle)

      if (!is.na(fdr.suffix)){
        fdr.col <- fdr.cols[ind.tmp]
        stopifnot(length(fdr.col)==1)
        graphics::hist(tab[,fdr.col], xlab="FDR", main=prefix)
      }
    }
    if (!is.na(name)){ grDevices::dev.off() }
  }
  return(invisible(tab.ss))
}
