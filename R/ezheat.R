#' Plot heatmap
#'
#' Draw heatmap using \code{pheatmap} package.
#'
#' @param object Matrix-like object with samples as columns and features as rows. Cannot have \code{NA}s.
#' @param symbols Labels corresponding to rows, e.g. gene symbols.
#' @param pheno.df Phenotype dataframe.
#' @param main Title of plot appended to \code{data.type}.
#' @param data.type Description of data values in \code{object}; this is included in main title.
#' @param color.v Color palette for heatmap. If \code{NULL}, it's set to
#' \code{colorRampPalette(rev(brewer.pal(n=9, name='RdYlBu')))(50)}.
#' @param stat.tab Matrix-like object with statistics, such as p-values, per column. If given, its dimensions should
#' match \code{object} and it also cannot have \code{NA}s.
#' @param cutoff Cutoff such that elements with \code{stats.tab < cutoff} show asterisk if \code{stats.tab} and
#' \code{cutoff} are not \code{NULL}.
#' @param name Name of PDF to plot. Set to \code{NA} to suppress writing to PDF.
#' @param reorder_rows Logical indicating if rows should be reordered.
#' @param reorder_cols Logical indicating if columns should be reordered.
#' @param fontsize_row Font size for row labels.
#' @param fontsize_col Font size for column labels.
#' @param only.symbols Logical, if \code{TRUE} only include rows where \code{symbols} isn't missing, as per \code{na.lab}.
#' @param unique.rows Logical, indicating if to remove duplicated row labels to make rows unique.
#' @param ntop Scalar number of rows to include.
#' @param sc Character string. Should rows be centered ('ctr'), z-scored ('z'), or neither ('none').
#' @param main Main title for plot.
#' @param na.lab Character vector of labels in \code{lab.col} to treat as missing, in addition to \code{NA}.
#' @param clip Values with magnitude > \code{clip} are reset to value \code{clip}. If given, must be > 0.
#' @details If the data after scaling and clipping (if they are used) has positive and negative values, the key is made
#' symmetric about zero.
#' @export
#' @import grDevices
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal

#can avoid dendro (& clustering) by eg cluster_rows=FALSE
ezheat <- function(object, symbols=NULL, pheno.df=NULL, main='Gene Expression', data.type='Log2', name='topgenes_heat.pdf',
                   color.v=NULL, unique.rows=FALSE, only.symbols=FALSE, ntop=NULL, stat.tab = NULL, cutoff = 0.05,
                   reorder_rows=FALSE, reorder_cols=FALSE, sc='ctr', clip=NA, fontsize_row=10, fontsize_col=10){

  stopifnot(sum(is.na(object)) == 0, sc %in% c('ctr', 'z', 'none'), is.na(clip)|(length(clip)==1 && clip > 0))
  if (!is.matrix(object)) object <- data.matrix(object)
  mat <- prune_mat(object, symbols = symbols, only.symbols = only.symbols, unique.rows = unique.rows, ntop = ntop)

  if (is.null(color.v)){
    if (!requireNamespace("RColorBrewer", quietly = TRUE)){
      stop("Package 'RColorBrewer' needed for this function to work. Please install it.", call. = FALSE)
    }
    color.v <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n=9, name='RdYlBu')))(50)
  }

  asterisk <- matrix(data = '', nrow = nrow(mat), ncol = ncol(mat))
  if (!is.null(stat.tab) & !is.null(cutoff)){
    stat.tab <- prune_mat(stat.tab, symbols = symbols, only.symbols = only.symbols, unique.rows = unique.rows,
                          ntop = ntop, verbose = FALSE)
    stopifnot(dim(stat.tab) == dim(mat), sum(is.na(stat.tab)) == 0)
    asterisk[stat.tab < cutoff] <- '*'
  } else {
    asterisk <- FALSE
  }

  # scale mat
  if (sc=='ctr'){
    #center but don't scale, so can see logFC
    mat <- t(scale(x=t(mat), scale=FALSE))
    data.type <- paste('Centered', data.type)
  } else if (sc=='z') {
    mat <- t(scale(x=t(mat), center=TRUE, scale=TRUE))
    data.type <- 'Z-score'
  }

  #reorder for heat w/o dendro
  if (reorder_rows){
    hc <- hclust(dist(mat, method = "euclidean"), method = "complete")
    mat <- mat[tree_row$order, , drop = FALSE]
  }
  if (reorder_cols){
    hc <- hclust(dist(t(mat), method = "euclidean"), method = "complete")
    mat <- mat[,tree_row$order, drop = FALSE]
    pheno.df <- pheno.df[tree_row$order, drop = FALSE]
  }

  # clip
  if (!is.na(clip)){
    mat[mat < -clip] <- -clip
    mat[mat > clip] <- clip
    data.type <- paste('Clipped', data.type)
  }

  #make key symmetric if has both pos & neg values
  if (any(mat > 0) & any(mat < 0)){
    len.colors <- length(color.v)
    ma <- max(abs(mat))
    breaks <- seq(from = -ma, to = ma, length.out = len.colors+1)
  } else {
    breaks <- NA
  }

  main <- paste(data.type, main)

  # params after name sent to grid::grid.text for asterisks, but vjust doesn't work
  pheatmap::pheatmap(mat, col=color.v, breaks = breaks, annotation_col = pheno.df, main=main, cluster_rows=FALSE,
                     cluster_cols=FALSE, fontsize_row=fontsize_row, fontsize_col=fontsize_col, display_numbers = asterisk,
                     filename=name)
}
