#' Multiple volcano plots
#'
#' Multiple volcano plots, where one or more comparisons are inferred from columns of \code{tab}
#' e.g. \code{AvsB.logFC}, and each comparison is plotted with \code{ezvolcano}.
#'
#' @param same.scale Logical, should axes' scales be the same for the different plots?
#' @inheritParams ezheat
#' @inheritParams ezvenn
#' @inheritParams ezvolcano
#' @details Uses \code{colnames(tab)} that have suffix \code{logFC} to infer comparisons.
#' @return Invisibly, a list of ggplot objects from \code{\link{ezvolcano}}.
#' @seealso ezvolcano
#' @export

multi_volcano <- function(tab, lab.col='Gene.Symbol', ntop.sig=0, ntop.lfc=0, alpha=0.4, name="volcanoes", ann.rnames=NULL,
                          base.size=14, up.ann.color="black", down.ann.color="black", same.scale=FALSE, type.sig=c("p", "FDR"),
                          up.cut.color="black", down.cut.color="black", cut.lfc=1, cut.sig=0.05, lines.sig=ifelse(type.sig[1]=="p", yes = 0.05, no=NA),
                          raster=FALSE, sep=".", na.lab=c("---", ""), seed=0, plot=TRUE){
  type.sig <- match.arg(type.sig)
  lfc.cols <- grep(paste0("\\", sep, "logFC$"), colnames(tab))
  if (length(lfc.cols)==0) stop("No logFC columns detected.")
  if (type.sig=="p"){
    sig.cols <- grep(paste0("\\", sep, "p"), colnames(tab))
  } else {
    sig.cols <- grep(paste0("\\", sep, "FDR"), colnames(tab))
  }
  if (length(sig.cols)==0) stop("No significance columns with suffix ", type.sig, " detected.")
  #use logFC cols instead of p cols to get contr.names, in case tab also has cor cols, whose p cols are not wanted
  contr.names <- sub(paste0("\\", sep, "(logFC)$"), "", colnames(tab)[lfc.cols])

  if(same.scale){
    x.bound <- max(abs(tab[,lfc.cols]))
    y.bound <- max(-log10(tab[,sig.cols]))
  } else {
    x.bound <- y.bound <- NULL
  }

  if (!is.na(name)) {
    grDevices::pdf(paste0(name, ".pdf"))
    on.exit(grDevices::dev.off())
  }
  ret.lst <- list()
  for (contr in contr.names){
    ret.lst[[contr]] <- ezvolcano(tab=tab, lab.col=lab.col, ntop.sig=ntop.sig, ntop.lfc=ntop.lfc, alpha=alpha, comparison=contr,
                                  name=NA, ann.rnames=ann.rnames, base.size=base.size, x.bound=x.bound, y.bound=y.bound, type.sig=type.sig,
                                  up.cut.color=up.cut.color, down.cut.color = down.cut.color,
                                  cut.lfc=cut.lfc, cut.sig=cut.sig, up.ann.color=up.ann.color, down.ann.color=down.ann.color,
                                  lines.sig=lines.sig, raster=raster, sep=sep, na.lab=na.lab, seed=seed, plot=plot)
  }
  return(invisible(ret.lst))
}
