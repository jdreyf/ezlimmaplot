#' Venn diagrams
#'
#' Venn diagrams using output from \code{ezlimma} package.
#'
#' @param tab Table of output from \code{ezlimma}.
#' @param prefix.v A vector of prefixes that prefix \code{.p}, \code{.FDR}, or \code{.logFC} in \code{colnames(tab)}.
#' These will be the columns that are used in the Venn diagram. If \code{NULL} these are inferred from
#' \code{colnames(tab)} that end with \code{.p}.
#' @param p.cutoff p-value cutoff to threshold features in \code{tab}. p-values < \code{p.cutoff} are significant.
#' @param fdr.cutoff FDR cutoff to threshold features in \code{tab}. FDRs < \code{fdr.cutoff} are significant.
#' @param logfc.cutoff log fold-change cutoff to threshold features in \code{tab}, when \code{tab} came from
#' \code{\link[ezlimma]{limma_contrasts}}.
#' @param circle.names Names of circles corresponding to \code{prefix.v}.
#' @param main Main title of plot.
#' @param name Name of PDF that gets "_venn" appended, then written. Set to \code{NA} to suppress writing to file.
#' @param cex A numerical value giving the amount by which plotting text and symbols should be magnified relative to
#' the default. See \code{\link[graphics]{par}}.
#' @param plot If the Venn diagram should be plotted.
#' @return Binary matrix indicating which features (rows) of \code{tab} were significant with specified cutoffs.
#' @details One of \code{fdr.cutoff} or \code{p.cutoff} must be given. If both are given, only \code{fdr.cutoff} is used.
#' \code{logfc.cutoff} if given is used in addition to these.
#' @return A matrix with elements: -1, 0, 1, invisibly. 0 indicates no significant change; -1 indicates down; and 1
#' indicates up if corresponding \code{logFC} columns are found, otherwise it indicates significance.
#' @export

ezvenn <- function(tab, prefix.v=NULL, p.cutoff = NULL, fdr.cutoff = NULL, logfc.cutoff = NULL, circle.names = prefix.v,
                   main = '', name = NA, cex = c(1, 1, 1), plot = TRUE){
  if (!requireNamespace("limma", quietly = TRUE)){
    stop("Package 'limma' needed for this function to work. Please install it.", call. = FALSE)
  }
  if (is.null(fdr.cutoff) & is.null(p.cutoff)){
    stop("One of 'p.cutoff' or 'fdr.cutoff' must be given.")
  }

  if (is.null(prefix.v)){
    p.cols <- grep(paste0('\\.p$'), colnames(tab))
    prefix.v <- sub(paste0('\\.p$'), '', colnames(tab)[p.cols])
  }

  if (!is.null(fdr.cutoff)){
    fdr.col <- paste0(prefix.v, '.FDR')
    stopifnot(fdr.col %in% colnames(tab))
    tab.sig <- tab[, fdr.col]
    tab.sig <- (tab.sig < fdr.cutoff)
  } else {
    p.col <- paste0(prefix.v, '.p')
    stopifnot(p.col %in% colnames(tab))
    tab.sig <- tab[, p.col]
    tab.sig <- (tab.sig < p.cutoff)
  }

  if (!is.na(name) & plot){
    name <- paste0(name, "_venn.pdf")
    grDevices::pdf(name)
  }

  logfc.col <- paste0(prefix.v, '.logFC')
  if (all(logfc.col %in% colnames(tab))){
    tab.logfc <- tab[,logfc.col]
    if(!is.null(logfc.cutoff)){
      tab.sig <- tab.sig * (abs(tab.logfc) > logfc.cutoff)
    }

    tab.sig <- tab.sig * sign(tab.logfc)
    tab.sig[is.na(tab.sig)] <- 0

    if (plot) limma::vennDiagram(tab.sig, include = c('up', 'down'), names = circle.names,
                         circle.col = rainbow(length(prefix.v)), counts.col = c('red', 'blue'),
                         main = main, cex = cex)
  } else {
    tab.sig[is.na(tab.sig)] <- 0
    if (plot) limma::vennDiagram(tab.sig, include = 'both', names = circle.names,
                         circle.col = rainbow(length(prefix.v)), counts.col = 'blue',
                         main = main, cex = cex)
  }

  if (!is.na(name) & plot) grDevices::dev.off()

  colnames(tab.sig) <- gsub('\\.(logFC|p|FDR)$', '.sig', colnames(tab.sig))
  tab.sig <- tab.sig[order(-abs(rowSums(tab.sig))), ]
  return(invisible(tab.sig))
}
