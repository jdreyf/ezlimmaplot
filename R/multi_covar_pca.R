#' PCA plots of top two components colored by each covariate
#'
#' PCA plot of top two principal components using \code{ggplot2} whose shape is defined by group and color in each plot
#' by iterating over each covariate.
#'
#' @param pheno.df Dataframe with rows as samples and columns as phenotypes.
#' @param grp.var Column of \code{pheno.df} to include as group variable in PCA plots.
#' @param covars Columns of \code{pheno.df} to include as covariates in PCA plots.
#' @inheritParams ezheat
#' @inheritParams ezpca
#' @details PCA is calculated with \code{\link[stats]{prcomp}}.
#' @return Invisibly, a list of the first two principal components appended to \code{pheno.df} for each covariate.
#' @export

multi_covar_pca <- function(object, pheno.df, name='covar_pca', grp.var='grp', covars=setdiff(colnames(pheno.df), grp.var),
                            alpha=1, all.size=NULL, facet=NULL, rm.leg.title=FALSE, labels=FALSE, manual.color = NULL,
                            manual.shape = NULL, plot=TRUE){
  stopifnot(c(grp.var, covars) %in% colnames(pheno.df))
  pca.lst <- list()
  grDevices::pdf(paste0(name, ".pdf"))
  for (cvr.ind in 1:length(covars)){
    pca.lst[[ covars[cvr.ind] ]] <- ezpca(object=object, pheno.df=pheno.df, shape=grp.var, color=covars[cvr.ind], name=NA,
                                  alpha=alpha, all.size=all.size, facet=facet, rm.leg.title=rm.leg.title, labels=labels,
                                  manual.color=manual.color, manual.shape=manual.shape, plot=plot)
  }
  grDevices::dev.off()
  return(invisible(pca.lst))
}
