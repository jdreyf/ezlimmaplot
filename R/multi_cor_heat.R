#' Plot one heatmap per phenotype correlation
#'
#' Plot one heatmap per phenotype correlation with \code{ezheat}, where features in \code{object}, e.g. genes, are reordered per comparison
#' using \code{tab}. This function is based on \code{multi_heat}, even though that function should probably be called \code{multi_contasts_heat}.
#'
#' @param p.thresh Only features with p-value at or below this will be plotted.
#' @param fdr.thresh Only features with FDR at or below this will be plotted.
#' @inheritParams ezheat
#' @inheritParams ezvenn
#' @inheritParams multi_heat
#' @inheritParams ezlimma::multi_cor
#' @details \code{rownames(tab)} and \code{rownames(object)} should overlap, \code{labrows} should correspond to \code{object}
#' and some \code{colnames(tab)} should end in \code{.p}, so they can be identified. If \code{fdr.thresh < 1}, then the
#' \code{colnames(tab)} that end in \code{.p} should be matched by \code{colnames(tab)} that end in \code{.FDR} instead of \code{.p}.
#'
#' To prevent this function from being called with an unnamed \code{labrows} that corresponds to \code{tab} instead of \code{object},
#' which is incorrect, if \code{labrows} is not \code{names(object)} (the default) then it must be named.
#' @export

multi_cor_heat <- function(tab, object, pheno.tab=NULL, labrows=rownames(object), labcols=colnames(object),
                       main="Log2 Expression", name="heats", sc="ctr", clip=NA, color.v=NULL,
                       unique.rows=FALSE, only.labrows=FALSE, ntop=50, stat.tab = NULL, p.thresh = 1, fdr.thresh=1,
                       cutoff = 0.05, reorder_rows=TRUE, reorder_cols=FALSE, fontsize_row=10, fontsize_col=10,
                       na.lab=c("---", ""), plot=TRUE, width=NA, height=NA, verbose=FALSE){

  if (length(labrows)==1) labrows <- rep(x=labrows, nrow(object))
  stopifnot(length(labrows)==nrow(object), names(labrows)==rownames(object))
  if (any(labrows != rownames(object))) stopifnot(!is.null(names(labrows)))
  if (all(labrows == rownames(object))) names(labrows) <- rownames(object)

  p.cols <- grep(paste0("\\.p$"), colnames(tab), value=TRUE)
  cor.names <- sub(paste0("\\.(p)$"), "", p.cols)

  stopifnot(length(intersect(rownames(tab), rownames(object))) > 1, length(p.cols) > 0, cor.names %in% colnames(pheno.tab))
  rows.int <- intersect(rownames(object), rownames(tab))
  tab <- tab[rows.int,, drop=FALSE]
  object <- object[rows.int,, drop=FALSE]
  labrows <- labrows[rows.int]

  if (!is.na(name)) {
    grDevices::pdf(paste0(name, ".pdf"))
    on.exit(grDevices::dev.off())
  }
  ret.lst <- list()
  for (ph.nm in cor.names){
    reorder_rows_tmp <- reorder_rows
    p.col <- paste0(ph.nm, ".p")
    tab.tmp <- tab |> dplyr::arrange(!!rlang::sym(p.col)) |>
      dplyr::filter(!!rlang::sym(p.col) <= p.thresh)
    if (fdr.thresh < 1){
      fdr.col <- paste0(ph.nm, ".FDR")
      stopifnot(fdr.col %in% colnames(tab))
      tab.tmp <- tab.tmp |> dplyr::filter(!!rlang::sym(fdr.col) <= fdr.thresh)
    }
    rows.tmp <- tab.tmp |> rownames()
    if (length(rows.tmp) > 0){
      pheno.df <- pheno.tab |> dplyr::select(!!rlang::sym(ph.nm)) |>
        dplyr::filter(!is.na(!!rlang::sym(ph.nm))) |>
        dplyr::arrange(!!rlang::sym(ph.nm))
      object.tmp <- object[rows.tmp, rownames(pheno.df), drop=FALSE]
      labrows.tmp <- labrows[rows.tmp]
      labcols.tmp <- labcols[match(colnames(object.tmp), colnames(object))]
      ntop.tmp <- min(ntop, length(rows.tmp))
      if (length(rows.tmp) == 1) reorder_rows_tmp <- FALSE

      ret.lst[[ph.nm]] <- ezheat(object=object.tmp, labrows=labrows.tmp, pheno.df=pheno.df, main=main, sc=sc, clip=clip,
                                 color.v=color.v, unique.rows=unique.rows, only.labrows=only.labrows, ntop=ntop.tmp,
                                 stat.tab = stat.tab, cutoff = cutoff, labcols = labcols.tmp, reorder_rows=reorder_rows_tmp,
                                 reorder_cols=reorder_cols, fontsize_row=fontsize_row, fontsize_col=fontsize_col,
                                 na.lab=na.lab, plot=FALSE, width=width, height=height, verbose=verbose, name=NA)
      if (plot) plot(ret.lst[[ph.nm]])
    }
  }
  return(invisible(ret.lst))
}
