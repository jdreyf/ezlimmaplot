#' PCA plots of top two components colored by each covariate
#'
#' PCA plot of top two principal components using \code{ggplot2} whose shape is defined by group and color in each plot
#' by iterating over each covariate.
#'
#' @param object Matrix-like object with samples as columns.
#' @param pheno.df Dataframe with rows as samples and columns as phenotypes.
#' @param name Name of PDF to write out.
#' @param grp.var Column of \code{pheno.df} to include as group variable in PCA plots.
#' @param covars Columns of \code{pheno.df} to include as covariates in PCA plots.
#' @param alpha Transparency, passed to \code{\link[ggplot2]{geom_point}}.
#' @param all.size Passed to \code{\link[ggplot2]{geom_point}} \code{size} parameter to give size for all points without
#' appearing in legend. \code{ggplot2} default is size=2.
#' @param facet A formula with columns in \code{pheno.df} to facet by.
#' @param rm.leg.title Logical indicating if legend title should be removed.
#' @param labels Logical indicating if sample labels should be added next to points.
#' @param manual.color Vector passed to \code{\link[ggplot2]{scale_colour_manual}}.
#' @param manual.shape Vector passed to \code{\link[ggplot2]{scale_shape_manual}}.
#' @details PCA is calculated with \code{\link[stats]{prcomp}}.
#' @return Invisibly, a list of the first two principal components appended to \code{pheno.df} for each covariate.
#' @export

multi_covar_pca <- function(object, pheno.df, name='covar_pca', grp.var='grp', covars=setdiff(colnames(pheno.df), grp.var),
                            alpha=1, all.size=NULL, facet=NULL, rm.leg.title=FALSE, labels=FALSE, manual.color = NULL,
                            manual.shape = NULL){
  stopifnot(c(grp.var, covars) %in% colnames(pheno.df))

  pca.lst <- list()

  pdf(paste0(name, ".pdf"))
  for (cvr.ind in 1:length(covars)){
    pca.lst[[ covars[cvr.ind] ]] <- ezpca(object=object, pheno.df=pheno.df, shape=grp.var, color=covars[cvr.ind], name=NA,
                                  alpha=alpha, all.size=all.size, facet=facet, rm.leg.title=rm.leg.title, labels=labels,
                                  manual.color=manual.color, manual.shape=manual.shape)
  }
  dev.off()

  return(invisible(pca.lst))
}
