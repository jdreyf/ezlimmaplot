#' Multiple Volcano plots
#'
#' Multiple Volcano plots from a table output from \code{ezlimma}.
#'
#' @param tab Table of output from \code{ezlimma}.
#' @param lab.col Column with labels, such as gene symbol, annotating features.
#' @param ntop.sig Number of top significant features to annotate.
#' @param ntop.lfc Number of top logFC features to annotate.
#' @param name Name of PDF file to write to. Set to \code{NA} to suppress writing to file.
#' @param add.rnames Additional rownames of features to annotate. These must be in \code{rownames(tab)}.
#' @param up.color Color for points that are upregulated (\code{logFC>0}).
#' @param down.color Color for points that are downregulated (\code{logFC<0}).
#' @param same.scale Logical indicating if different volcano plots should have the same x-limits and y-limits.
#' @param type.sig Type of significance y-axis should use, either "p" or "FDR".
#' @param cut.color Color of points that meet both \code{cut.lfc} and \code{cut.sig}. If \code{NULL}, cutoffs are ignored.
#' @param cut.lfc Points need to have \code{|logFC| >= cut.lfc} to have \code{cut.color}.
#' @param cut.sig Points need to have significance \code{<= cut.sig} to have \code{cut.color}. Significance type is of
#' \code{type.sig}.
#' @param sep Separator string between contrast names and suffix such as \code{logFC}.
#' @param na.lab Character vector of labels in \code{lab.col} to treat as missing, in addition to \code{NA}.
#' @details Uses colnames(tab) that have suffix \code{logFC} to infer comparisons.
#' @return List of ggplot objects from \code{\link{ezvolcano}}, invisibly.
#' @export
#' @import grDevices

multi_volcano <- function(tab, lab.col='Gene.Symbol', ntop.sig=10, ntop.lfc=0, name='volcanoes', add.rnames=NULL,
                          up.color='black', down.color='black', same.scale=FALSE, type.sig=c('p', 'FDR'), cut.color=NULL,
                          cut.lfc=1, cut.sig=0.05, sep='.', na.lab=c('---', '')){

  lfc.cols <- grep(paste0('\\', sep, 'logFC$'), colnames(tab))
  if (length(lfc.cols)==0) stop("No logFC columns detected.")
  if (type.sig=="p"){
    sig.cols <- grep(paste0('\\', sep, 'p'), colnames(tab))
  } else {
    sig.cols <- grep(paste0('\\', sep, 'FDR'), colnames(tab))
  }
  if (length(sig.cols)==0) stop("No significance columns with suffix", type.sig, "detected.")
  #use logFC cols instead of p cols to get contr.names, in case tab also has cor cols
  contr.names <- sub(paste0('\\', sep, '(logFC)$'), '', colnames(tab)[lfc.cols])

  if(same.scale){
    x.bound <- max(abs(tab[,lfc.cols]))
    y.bound <- max(-log10(tab[,sig.cols]))
  } else {
    x.bound <- y.bound <- NULL
  }

  if (!is.na(name)) grDevices::pdf(paste0(name, '.pdf'))
  ret.lst <- list()
  for (contr in contr.names){
    ret.lst[[contr]] <- ezvolcano(tab=tab, lab.col=lab.col, ntop.sig=ntop.sig, ntop.lfc=ntop.lfc, comparison=contr,
                                  name=NA, add.rnames=add.rnames, up.color=up.color, down.color=down.color, x.bound=x.bound,
                                  y.bound=y.bound, type.sig=type.sig, cut.color=cut.color, cut.lfc=cut.lfc, cut.sig=cut.sig,
                                  sep=sep, na.lab=na.lab)
  }
  if (!is.na(name)) grDevices::dev.off()
  return(invisible(ret.lst))
}
