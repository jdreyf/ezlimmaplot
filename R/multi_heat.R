#' Plot one heatmap per comparison
#'
#' Plot one heatmap per comparison with \code{ezheat}, where features in \code{object} are reordered per comparison
#' using \code{tab}.
#'
#' @inheritParams ezheat
#' @inheritParams ezvenn
#' @details \code{rownames(tab)} and \code{rownames(object)} should overlap, and some \code{colnames(tab)}
#' should end in \code{.p}, so they can be identified.
#' @export

multi_heat <- function(tab, object, labrows=NULL, pheno.df=NULL, main="Log2 Expression", name="heats",
                   sc="ctr", clip=NA, color.v=NULL, unique.rows=FALSE, only.labrows=FALSE, ntop=50, stat.tab = NULL,
                   cutoff = 0.05, labcols=NULL, reorder_rows=FALSE, reorder_cols=FALSE, fontsize_row=10, fontsize_col=10,
                   na.lab=c("---", ""), plot=TRUE, width=NA, height=NA, verbose=FALSE){

  p.cols <- grep(paste0("\\.p$"), colnames(tab), value=TRUE)
  contr.names <- sub(paste0("\\.(p)$"), "", p.cols)

  stopifnot(length(intersect(rownames(tab), rownames(object))) > 1, length(p.cols) > 0)

  if (!is.na(name)) grDevices::pdf(paste0(name, ".pdf"))
  ret.lst <- list()
  for (contr in contr.names){
    main.tmp <- paste(main, contr)
    p.col <- paste0(contr, ".p")
    rows.tmp <- rownames(tab)[order(tab[,p.col])]
    object.tmp <- object[rows.tmp,]

    ret.lst[[contr]] <- ezheat(object=object.tmp, labrows=labrows, pheno.df=pheno.df, main=main.tmp, sc=sc, clip=clip,
                               color.v=color.v, unique.rows=unique.rows, only.labrows=only.labrows, ntop=ntop,
                               stat.tab = stat.tab, cutoff = cutoff, labcols=labcols, reorder_rows=reorder_rows,
                               reorder_cols=reorder_cols, fontsize_row=fontsize_row, fontsize_col=fontsize_col,
                               na.lab=na.lab, plot=plot, width=width, height=height, verbose=verbose, name=NA)
  }
  if (!is.na(name)) grDevices::dev.off()
  return(invisible(ret.lst))
}
