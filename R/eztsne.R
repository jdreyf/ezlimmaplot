#' tSNE plot of first two dimensions
#'
#' tSNE plot of first two dimensions using \code{ggplot2}.
#'
#' @param seed Random seed.
#' @inheritParams ezpca
#' @inheritParams Rtsne::Rtsne
#' @param ... Passed to \code{\link[ggplot2:aes_]{aes_string}}.
#' @details \code{object} must have colnames, and if \code{pheno.df}
#' is given, it is checked that \code{colnames(object)==rownames(pheno.df)}.
#' @return Invisibly, first two dimensionss appended to \code{pheno.df}.
#' @export

eztsne<- function(object, pheno.df, name='tsne', check_duplicates=FALSE, pca=TRUE,
                  perplexity=min(30, round(ncol(object)/5)), theta=0.5, seed = 123, initial_dims=50, max_iter=1000,
                  alpha=1, all.size=NULL, facet=NULL, title = NULL, subtitle = NULL, rm.leg.title=FALSE, labels=FALSE,
                  manual.color = NULL, manual.shape = NULL, ...){

  stopifnot(ncol(object)==nrow(pheno.df), colnames(object)==rownames(pheno.df))

  set.seed(seed)
  tsne1=Rtsne::Rtsne(t(as.matrix(object)), check_duplicates=check_duplicates, pca=pca, perplexity=perplexity,
                     theta=theta, dims=2, initial_dims=initial_dims, max_iter=max_iter)
  tsne1 <- as.data.frame(tsne1$Y)
  colnames(tsne1) <- c("tSNE1", "tSNE2")
  dat <- data.frame(pheno.df, tsne1)

  dots <- list(...)
  if(is.null(names(dots))){
    n <- 0
  }else{
    chars <- vector("list", 2*length(dots))
    for(i in seq_along(dots)){
      chars[[2*i]] <- dots[[i]]
      chars[[2*i-1]] <- as.character(dat[, dots[[i]]])
    }
    n <- max(nchar(unlist(chars)))
  }

  width <- 7 + n / 12
  if (!is.na(name)){ pdf(paste0(name, ".pdf"), width = width, height = 7) }

  #need to set alpha/all.size in geom_point, else it appears in legend
  qp <- ggplot2::ggplot(dat, mapping=ggplot2::aes_string(x='tSNE1', y='tSNE2', ...)) + ggplot2::theme_bw()

  if (!is.null(all.size)){
    qp <- qp + ggplot2::geom_point(size=all.size, alpha=alpha)
  } else {
    qp <- qp + ggplot2::geom_point(alpha=alpha)
  }

  if (!is.null(facet)){ qp <- qp + ggplot2::facet_grid(facet) }

  if (rm.leg.title){ qp <- qp + ggplot2::theme(legend.title=ggplot2::element_blank()) }

  if (!is.null(title)) { qp <- qp + ggplot2::ggtitle(label = title, subtitle = subtitle) }

  if (labels){
    dat2 <- dat
    dat2$row_names <- rownames(pheno.df)
    qp <- qp + ggplot2::geom_text(data = dat2, mapping=ggplot2::aes_string(label='row_names'), size=2, vjust=-.7)
  }

  if(!is.null(manual.color)) qp <- qp + ggplot2::scale_colour_manual(values = manual.color)
  if(!is.null(manual.shape)) qp <- qp + ggplot2::scale_shape_manual(values = manual.shape)

  graphics::plot(qp)
  if (!is.na(name)){ dev.off() }

  return(invisible(dat))

}
