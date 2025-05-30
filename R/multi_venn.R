#' Multiple Venn diagrams
#'
#' Plot multiple Venn diagrams using output from \code{ezlimma} package. Each is a page in a PDF, and a single table
#' with all results is invisibly returned.
#'
#' @param prefix.lst A list of \code{prefix.v} vectors passed to \code{\link{ezvenn}}.
#' @param circle.names.lst List of length \code{length(prefix.lst)}, where each element is a \code{circle.names} vector.
#' Passed to \code{\link{ezvenn}}.
#' @param main.v Vector of names of length \code{length(prefix.lst)}.
#' @inheritParams ezheat
#' @inheritParams ezvenn
#' @details One of \code{fdr.cutoff} or \code{p.cutoff} must be given. If both are given, only \code{fdr.cutoff} is used.
#' \code{logfc.cutoff} if given is used in addition to these.
#'
#' \code{"_venn.pdf"} is appended to \code{name}.
#'
#' If \code{main.v==""} and \code{names(prefix.lst)} is not \code{NULL}, \code{main.v} is assigned
#' \code{names(prefix.lst)}.
#'
#' @return Invisibly, a matrix with elements \code{-1, 0, 1} indicating which features (rows) of \code{tab} were significant
#' with specified cutoffs. 0 indicates no significant change; -1 indicates down; and 1 indicates up if corresponding
#' \code{logFC} columns are found, otherwise 1 indicates significance.
#' @export

multi_venn <- function(tab, prefix.lst, p.cutoff = NULL, fdr.cutoff = NULL, logfc.cutoff = NULL, circle.names.lst=prefix.lst,
                   main.v = "", name = "multi", cex = c(1, 1, 1)){
  if (is.null(circle.names.lst)) circle.names.lst <- prefix.lst
  stopifnot(is.list(prefix.lst), length(prefix.lst) > 0, length(prefix.lst)==length(circle.names.lst))

  if (!is.null(names(prefix.lst)) && main.v == "") main.v <- names(prefix.lst)

  if (!is.na(name)) {
    grDevices::pdf(paste0(name, "_venn.pdf"))
    on.exit(grDevices::dev.off())
  }
  tab.lst <- list()
  for (ind in 1:length(prefix.lst)){
    tab.lst[[ind]] <- ezvenn(tab=tab, prefix.v=prefix.lst[[ind]], p.cutoff = p.cutoff, fdr.cutoff = fdr.cutoff,
                             logfc.cutoff = logfc.cutoff, circle.names = circle.names.lst[[ind]], main = main.v[ind],
                             name = NA, cex = cex)
    if (ind == 1){
      tab.sig <- tab.lst[[ind]]
    } else {
      new.cols <- setdiff(colnames(tab.lst[[ind]]), colnames(tab.sig))
      if (length(new.cols) >= 1){
        tab.sig <- cbind(tab.sig, tab.lst[[ind]][rownames(tab.sig), new.cols, drop=FALSE])
      }
    }#end else
  }#end for
  tab.sig <- tab.sig[order(rowSums(abs(tab.sig)), decreasing = TRUE), ]
  return(invisible(tab.sig))
}
