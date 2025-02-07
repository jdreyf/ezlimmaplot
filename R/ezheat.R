#' Plot heatmap
#'
#' Draw heatmap using \pkg{pheatmap} package.
#'
#' @param object Matrix-like object with features (e.g. genes) as rows and samples as columns.
#' @param pheno.df Data frame with rows as samples and columns as phenotypes.
#' @param labrows Labels for rows, e.g. gene symbols. This can be of length 1 (in which case it is recycled),
#' of length \code{nrow(object)}.
#' @param labcols Labels for columns. This can be of length 1 (in which case it is recycled),
#' of length \code{ncol(object)}.
#' @param main Main title of plot.
#' @param name Name of file to create. Set to \code{NA} to plot to screen instead of to file.
#' @param sc Row scaling. Should rows be centered ('ctr'), z-scored ('z'), or neither ('none').
#' @param clip Values with magnitude > \code{clip} are reset to value \code{clip}. If given, must be > 0.
#' @param color.v Color palette for heatmap. If \code{NULL}, it's set to
#' \code{colorRampPalette(rev(brewer.pal(n=9, name='RdYlBu')))(50)}.
#' @param unique.rows Logical; remove duplicated row labels to make rows unique?
#' @param only.labrows Logical; only include rows where \code{labrows} aren't missing
#' (missing is defined by \code{na.lab})?
#' @param ntop Scalar number of rows to include.
#' @param stat.tab Matrix-like object with statistics, such as p-values, per column, used to add "*". If given, its
#' dimensions should match \code{object} and it cannot have \code{NA}s.
#' @param cutoff Cutoff such that elements with \code{stats.tab < cutoff} show asterisk if \code{stats.tab} and
#' \code{cutoff} are not \code{NULL}.
#' @param reorder_rows Logical; should rows be reordered with hierarchical clustering?
#' @param reorder_cols Logical; should columns be reordered with hierarchical clustering?
#' @param annotation_colors list for specifying annotation_row and annotation_col track colors manually.
#' It is possible to define the colors for only some of the features. See \code{\link[pheatmap]{pheatmap}} examples.
#' @param gaps_row vector of row indices after which to insert gaps.
#' @param gaps_col vector of column indices after which to insert gaps.
#' @param fontsize_row Font size for row labels.
#' @param fontsize_col Font size for column labels.
#' @param na.lab Character vector of labels in \code{lab.col} to treat as missing, in addition to \code{NA}.
#' @param plot Logical; should plot be generated?
#' @param width Manual option for determining the output file width in inches.
#' @param height Manual option for determining the output file height in inches.
#' @param verbose Logical; print pruning messages to console?
#' @inheritParams pheatmap::pheatmap
#' @return A list with element \code{mat} of the matrix of values plotted, and if \code{plot=TRUE} element
#' \code{gtable}, containing the \code{gtable} object returned by \code{\link[pheatmap]{pheatmap}}.
#' @details \code{main} is modified to describe data transformations.
#'
#' If the data after scaling and clipping (if they are used) has positive and negative values, the key is made
#' symmetric about zero.
#'
#' \code{object} cannot have \code{NA}s.
#' @export

# can avoid dendro (& clustering) by eg cluster_rows=FALSE
# keep object & labrows separate, to use annotate_row w/o labrows
# labrows works like labcols for consistency
ezheat <- function(object, pheno.df=NULL, labrows=rownames(object), labcols=colnames(object),
                   main="Log2 Expression", name="topgenes_heat", sc="ctr", clip=NA,
                   color.v=NULL, unique.rows=FALSE, only.labrows=FALSE, ntop=NULL, stat.tab = NULL,
                   cutoff = 0.05, reorder_rows=FALSE, reorder_cols=FALSE, gaps_row = NULL, gaps_col = NULL,
                   annotation_row = NA, annotation_colors = NA, angle_col=c("270", "0", "45", "90", "315"),
                   fontsize_row=10, fontsize_col=10, na.lab=c("---", ""), plot=TRUE, width=7, height=7, verbose=FALSE){
  angle_col <- match.arg(angle_col)
  if (!is.matrix(object)) object <- data.matrix(object)
  stopifnot(sum(is.na(object)) == 0, sc %in% c("ctr", "z", "none"), is.na(clip) | (length(clip)==1 && clip > 0),
            length(labrows) %in% c(1, nrow(object)), length(labcols) %in% c(1, ncol(object)),
            is.null(stat.tab) || !any(is.na(stat.tab)))
  if (length(labrows)==1) labrows <- rep(x=labrows, nrow(object))
  if (length(labcols)==1) labcols <- rep(x=labcols, ncol(object))

  if (!is.null(pheno.df)){
    stopifnot(ncol(object)==nrow(pheno.df), colnames(object)==rownames(pheno.df))
  }
  # ok if labrows NULL
  labrows <- prune_mat(object, labrows = labrows, only.labrows = only.labrows, unique.rows = unique.rows, ntop = ntop,
                   na.lab=na.lab, verbose=verbose)
  mat <- object[names(labrows),, drop=FALSE]

  if (is.null(color.v)){
    color.v <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n=9, name="RdYlBu")))(50)
  }

  # scale mat
  if (sc=="ctr"){
    #center but do not scale, so can see logFC
    mat <- t(scale(x=t(mat), scale=FALSE))
    main <- paste("Centered", main)
  } else if (sc=="z") {
    mat <- t(scale(x=t(mat), center=TRUE, scale=TRUE))
    main <- paste("Z-scored", main)
  }

  #reorder for heat w/o dendro
  if (reorder_rows){
    hc <- stats::hclust(stats::dist(mat, method = "euclidean"), method = "ward.D2")
    mat <- mat[hc$order,, drop = FALSE]
    labrows <- labrows[hc$order]
  }

  if (reorder_cols){
    hc <- stats::hclust(stats::dist(t(mat), method = "euclidean"), method = "ward.D2")
    mat <- mat[, hc$order, drop = FALSE]
    pheno.df <- pheno.df[hc$order,, drop = FALSE]
    if (!is.null(labcols)){
      labcols <- labcols[hc$order]
    }
  }

  asterisk <- matrix(data = "", nrow = nrow(mat), ncol = ncol(mat))
  if (!is.null(stat.tab) & !is.null(cutoff)){
    stat.tab <- stat.tab[rownames(mat),, drop=FALSE]
    stopifnot(dim(stat.tab) == dim(mat), !is.na(stat.tab))
    asterisk[stat.tab < cutoff] <- "*"
  } else {
    asterisk <- FALSE
  }

  # clip
  if (!is.na(clip)){
    mat[mat < -clip] <- -clip
    mat[mat > clip] <- clip
    main <- paste("Clipped", main)
  }

  #make key symmetric if has both pos & neg values
  if (any(mat > 0) & any(mat < 0)){
    len.colors <- length(color.v)
    ma <- max(abs(mat))
    breaks <- seq(from = -ma, to = ma, length.out = len.colors+1)
  } else {
    breaks <- NA
  }

  # params after name sent to grid::grid.text for asterisks, but vjust doesn't work
  ph <- pheatmap::pheatmap(mat, col=color.v, breaks = breaks, annotation_col = pheno.df, main=main,
                           cluster_rows=FALSE, cluster_cols=FALSE, gaps_col=gaps_col, gaps_row=gaps_row,
                           annotation_row=annotation_row, annotation_colors=annotation_colors,
                           labels_row = labrows, labels_col = labcols, angle_col=angle_col,
                           fontsize_row=fontsize_row, fontsize_col=fontsize_col,
                           display_numbers = asterisk, filename=NA, silent=TRUE)

  if (plot){
    # fname <- ifelse(is.na(name), NA, paste0(name, ".pdf"))
    if (!is.na(name)) {  # Only open PDF if a filename is provided
      grDevices::pdf(file = paste0(name, ".pdf"), width = width, height = height) # Open PDF device *manually*
      on.exit(grDevices::dev.off())  # Ensure device is closed even if errors occur
    }

    grid::grid.newpage()
    grid::grid.draw(ph$gtable)
  }
  ret <- list(mat=mat, gtable=ph$gtable)
  return(invisible(ret))
}
