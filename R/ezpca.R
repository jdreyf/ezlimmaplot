#' PCA plot of top two components
#'
#' PCA plot of top two principal components using \code{ggplot2}.
#'
#' @param object Matrix-like object with samples as columns
#' @param pheno Dataframe of phenotypes
#' @param name Name of PNG to write out. Set to \code{NA} to suppress writing to file.
#' @param scale. Logical indicating whether to scale in \code{\link[stats]{prcomp}}.
#' @param alpha Transparency, passed to \code{\link[ggplot2]{geom_point}}.
#' @param all.size Passed to \code{\link[ggplot2]{geom_point}} \code{size} parameter to give size for all points without appearing in legend. \code{ggplot2} default is size=2.
#' @param facet Variable to facet by.
#' @param rm.leg.title Logical indicating if legend title should be removed.
#' @param labels Logical indicating if sample labels should be added next to points.
#' @param manual.color Values passed to \code{\link[ggplot2]{scale_colour_manual}}.
#' @param manual.shape Values passed to \code{\link[ggplot2]{scale_shape_manual}}.
#' @param ... Passed to \code{\link[ggplot2]{ggplot}} \code{aes_string} parameter.
#' @export
#' @import ggplot2
#' @import stats

ezpca <- function(object, pheno, name='pca', scale.=FALSE, alpha=1, all.size=NULL, facet=NULL, rm.leg.title=FALSE,
                  labels=FALSE, manual.color = NULL, manual.shape = NULL, ...){
  if (!requireNamespace("ggplot2", quietly = TRUE)){
    stop("Package 'ggplot2' needed for this function to work. Please install it.", call. = FALSE)
  }
  stopifnot(ncol(object)==nrow(pheno), colnames(object)==rownames(pheno))

  pca <- stats::prcomp(t(object[rowSums(is.na(object))==0,]), scale.=scale.)
  pve <- signif(summary(pca)$importance['Proportion of Variance', 1:2]*100, 2)

  dat <- data.frame(pca$x[rownames(pheno), 2:1], pheno)

  #need to set alpha/all.size in geom_point, else it appears in legend
  qp <- ggplot(dat, aes_string(y='PC1', x='PC2', ...)) + theme_bw() + theme(panel.grid=element_line(color='black'))
  if (!is.null(all.size)){ qp <- qp + geom_point(size=all.size, alpha=alpha) } else { qp <- qp + geom_point(alpha=alpha) }
  if (!is.null(facet)){ qp <- qp + facet_grid(facet) }
  qp <- qp + ylab(paste0('PC1 (', pve[1], '%)')) + xlab(paste0('PC2 (', pve[2], '%)', sep=''))
  if (rm.leg.title){ qp <- qp + theme(legend.title=element_blank()) }
  if (labels){
    dat2 <- dat
    dat2$row_names <- rownames(pheno)
    qp <- qp + geom_text(data=dat2, mapping=aes_string(y='PC1', x='PC2', label='row_names'), size=2, vjust=-.7)
  }

  if(!is.null(manual.color)) qp <- qp + scale_colour_manual(values = manual.color)
  if(!is.null(manual.shape)) qp <- qp + scale_shape_manual(values = manual.shape)

  if (!is.na(name)){ ggsave(filename=paste0(name, '.png'), plot=qp) } else { print(qp) }

  return(invisible(dat))
}
