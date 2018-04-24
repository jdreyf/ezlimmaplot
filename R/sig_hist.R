#' Plot histograms of significance (p-value & FDR) columns
#'
#' Plot histograms of significance (p-value & FDR) columns.
#'
#' @param tab Table of output from \code{ezlimma}.
#' @param p.ext Extension for p-value columns.
#' @param fdr.ext Extension for FDR columns.
#' @param sep Separator for column names before \code{p} or \code{FDR}.
#' @param pi0 Logical indicating if proportion of null hypotheses be calculated per p-value histogram.
#' @param name Name of PDF file to write. Set to \code{NA} to suppress writing to PDF.
#' @export

#assume each comparison has a p-value & q-value column unless fdr.ext=NA
sig_hist <- function(tab, p.ext='pval|p', fdr.ext='q|FDR', sep='.', pi0 = FALSE, name='sig_hist'){
  p.cols <- grep(paste0('\\', sep, '(', p.ext, ')$'), colnames(tab))
  prefix.v <- sub(paste0('\\', sep, '(', p.ext, ')$'), '', colnames(tab)[p.cols])

  #set name=NA to turn off pdf
  if (!is.na(name)){ pdf(paste0(name, '.pdf')) }
  par(mfrow=c(2,2))
  for (contr in prefix.v){
    p.col <- which(colnames(tab) %in% paste(contr, unlist(strsplit(p.ext, split='|', fixed=TRUE)), sep=sep))
    stopifnot(length(p.col)==1)
    subtitle <- NULL
    if(pi0) {
      prop.null <- limma::propTrueNull(tab[,p.col], method = 'convest')
      subtitle <- paste('Proportion of True Null = ', signif(prop.null, 3))
    }
    hist(tab[,p.col], xlab='P-value', main=contr, sub = subtitle)
    if (!is.na(fdr.ext)){
      fdr.col <- which(colnames(tab) %in% paste(contr, unlist(strsplit(fdr.ext, split='|', fixed=TRUE)), sep=sep))
      stopifnot(length(fdr.col)==1)
      hist(tab[,fdr.col], xlab='FDR', main=contr)
    }
  }
  if (!is.na(name)){ dev.off() }
  invisible(TRUE)
}
