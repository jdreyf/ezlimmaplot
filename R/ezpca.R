#' PCA plot of top two components
#'
#' PCA plot of top two principal components using \code{ggplot2}.
#'
#' @param alpha Transparency, passed to \code{\link[ggplot2]{geom_point}}.
#' @param all.size Passed to \code{\link[ggplot2]{geom_point}} \code{size} parameter to give size for all points without
#' appearing in legend. \code{ggplot2} default is size=2.
#' @param facet A formula with columns in \code{pheno.df} to facet by.
#' @param title Title text; suppressed if it is \code{NULL}.
#' @param subtitle Subtitle text; suppressed if it is \code{NULL} or \code{title} is \code{NULL}. If you'd like a
#' \code{subtitle} but no \code{title}, set \code{title = ""}.
#' @param rm.leg.title Logical indicating if legend title should be removed.
#' @param labels Logical, should sample labels be added next to points?
#' @param manual.color Vector passed to \code{\link[ggplot2:scale_manual]{scale_colour_manual}}.
#' @param manual.shape Vector passed to \code{\link[ggplot2:scale_manual]{scale_shape_manual}}.
#' @inheritParams ezheat
#' @param ... Passed to \code{\link[ggplot2:aes_]{aes_string}}.
#' @details PCA is calculated with \code{\link[stats]{prcomp}}. \code{object} must have colnames, and if \code{pheno.df}
#' is given, it is checked that \code{colnames(object)==rownames(pheno.df)}.
#' @return Invisibly, first two principal components appended to \code{pheno.df}.
#' @export

ezpca <- function(object, pheno.df=NULL, name="pca", alpha=1, all.size=NULL, facet=NULL, title=NULL, subtitle=NULL,
                  rm.leg.title=FALSE, labels=FALSE, manual.color = NULL, manual.shape = NULL, plot=TRUE, ...){
  stopifnot(limma::isNumeric(object), nrow(object[rowSums(is.na(object))==0,]) > 0, ncol(object) > 0,
            !is.null(colnames(object)), is.logical(plot))

  pca <- stats::prcomp(t(object[rowSums(is.na(object))==0,]), scale. = FALSE)
  pve <- signif(summary(pca)$importance["Proportion of Variance", 1:2]*100, 2)

  if (!is.null(pheno.df)){
    stopifnot(ncol(object)==nrow(pheno.df), colnames(object)==rownames(pheno.df))
    dat <- data.frame(pca$x[rownames(pheno.df), 1:2], pheno.df, check.names = FALSE)
    dat$row_names <- rownames(pheno.df)
  } else{
    dat <- data.frame(pca$x[colnames(object), 1:2], check.names = FALSE)
    dat$row_names <- colnames(object)
  }

  dots <- list(...)
  if (is.null(names(dots))){
    n <- 0
  } else {
    chars <- vector("list", 2*length(dots))
    for (i in seq_along(dots)){
      chars[[2*i]] <- dots[[i]]
      # want length of longest element, but dots[[i]] can be an expression not in colnames(dat)
      # then don"t parse & ignore length
      chars[[2*i-1]] <- ifelse(dots[[i]] %in% colnames(dat), as.character(dat[, dots[[i]] ]), "")
    }
    n <- max(nchar(unlist(chars)), na.rm = TRUE)
  }

  width <- 7 + n/12
  if (!is.na(name)){
    grDevices::pdf(paste0(name, ".pdf"), width = width, height = 7)
    on.exit(grDevices::dev.off())
  }

  # need to set alpha/all.size in geom_point, else it appears in legend
  qp <- ggplot2::ggplot(dat, mapping=ggplot2::aes_string(x="PC1", y="PC2", ...)) + ggplot2::theme_bw()
  if (!is.null(all.size)){
    qp <- qp + ggplot2::geom_point(size=all.size, alpha=alpha)
  } else {
    qp <- qp + ggplot2::geom_point(alpha=alpha)
  }
  if (!is.null(facet)){ qp <- qp + ggplot2::facet_grid(facet) }
  qp <- qp + ggplot2::xlab(paste0("PC1 (", pve[1], "%)")) + ggplot2::ylab(paste0("PC2 (", pve[2], "%)", sep=""))
  if (rm.leg.title){ qp <- qp + ggplot2::theme(legend.title=ggplot2::element_blank()) }
  if (!is.null(title)){ qp <- qp + ggplot2::ggtitle(label=title, subtitle=subtitle) }

  if (labels){
    qp <- qp + ggplot2::geom_text(data=dat, mapping=ggplot2::aes_string(x="PC1", y="PC2", label="row_names"),
                                  size=2, vjust=-.7)
  }

  if(!is.null(manual.color)) qp <- qp + ggplot2::scale_colour_manual(values = manual.color)
  if(!is.null(manual.shape)) qp <- qp + ggplot2::scale_shape_manual(values = manual.shape)

  if (plot) graphics::plot(qp)
  return(invisible(dat))
}
