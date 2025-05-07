#' Plot one heatmap per comparison
#'
#' Plot one heatmap per comparison with \code{ezheat}, where features in \code{object} are reordered per comparison
#' using \code{tab}.
#'
#' @param only.contr.cols Logical; only include samples involved in each comparison? Otherwise shows all samples in all comparisons to
#' information content.
#' @param grp.var Character; name of group variable in \code{pheno.df} to select samples if \code{only.contr.cols} is true.
#' @inheritParams ezheat
#' @inheritParams ezvenn
#' @details \code{rownames(tab)} and \code{rownames(object)} should overlap, \code{labrows} should correspond to \code{object}
#' and some \code{colnames(tab)} should end in \code{.p}, so they can be identified. If \code{only.contr.cols}, comparisons are parsed from
#' p-value columns, whose names should include \code{vs}.
#'
#' To prevent this function from being called with an unnamed \code{labrows} that corresponds to \code{tab} instead of \code{object},
#' which is incorrect, if \code{labrows} is not \code{names(object)} (the default) then it must be named.
#' @export

multi_heat <- function(tab, object, pheno.df=NULL, labrows=rownames(object), labcols=colnames(object),
                       main="Log2 Expression", name="heats", sc="ctr", clip=NA, color.v=NULL,
                       unique.rows=FALSE, only.labrows=FALSE, only.contr.cols=FALSE, grp.var="grp", ntop=50, stat.tab = NULL,
                       cutoff = 0.05, reorder_rows=TRUE, reorder_cols=FALSE, fontsize_row=10, fontsize_col=10,
                       na.lab=c("---", ""), plot=TRUE, width=7, height=7, verbose=FALSE){
  if (length(labrows)==1) labrows <- rep(x=labrows, nrow(object))
  stopifnot(length(labrows)==nrow(object), names(labrows)==rownames(object))
  if (any(labrows != rownames(object))) stopifnot(!is.null(names(labrows)))
  if (all(labrows == rownames(object))) names(labrows) <- rownames(object)

  p.cols <- grep(paste0("\\.p$"), colnames(tab), value=TRUE)
  contr.names <- sub(paste0("\\.(p)$"), "", p.cols)

  stopifnot(length(intersect(rownames(tab), rownames(object))) > 1, length(p.cols) > 0)
  rows.int <- intersect(rownames(object), rownames(tab))
  tab <- tab[rows.int,, drop=FALSE]
  object <- object[rows.int,, drop=FALSE]
  labrows <- labrows[rows.int]

  ret.lst <- list()
  for (contr in contr.names){
    main.tmp <- paste(main, contr)
    p.col <- paste0(contr, ".p")
    rows.tmp <- rownames(tab)[order(tab[,p.col])]
    if (only.contr.cols){
      cols.tmp <- rownames(pheno.df)[pheno.df[[grp.var]] %in% unlist(strsplit(contr, split="(_|)vs(_|)"))]
      object.tmp <- object[rows.tmp, cols.tmp]
      pheno.df <- pheno.df[cols.tmp,, drop=FALSE]
      # intersect appears to keep order of first input
      labcols <- intersect(cols.tmp, labcols)
    } else {
      object.tmp <- object[rows.tmp,]
    }

    labrows.tmp <- labrows[rows.tmp]

    ret.lst[[contr]] <- ezheat(object=object.tmp, labrows=labrows.tmp, pheno.df=pheno.df, main=main.tmp, sc=sc, clip=clip,
                               color.v=color.v, unique.rows=unique.rows, only.labrows=only.labrows, ntop=ntop,
                               stat.tab = stat.tab, cutoff = cutoff, labcols=labcols, reorder_rows=reorder_rows,
                               reorder_cols=reorder_cols, fontsize_row=fontsize_row, fontsize_col=fontsize_col,
                               na.lab=na.lab, plot=FALSE, verbose=verbose, name=NA)
  }
  if (plot){
    if (!is.na(name)) {
      grDevices::pdf(paste0(name, ".pdf"), width = width, height = height)
      on.exit(grDevices::dev.off())
    }
    for (contr in contr.names){
      grid::grid.newpage()
      grid::grid.draw(ret.lst[[contr]]$gtable)
    }
  }
  return(invisible(ret.lst))
}
