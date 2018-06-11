#' Make boxplot or dotplots for a feature per group
#'
#' Make boxplot or dotplots for a feature per group using \code{ggplot2}.
#'
#' @param object A matrix-like data object containing log-ratios or log-expression values for a series of arrays, with
#' rows corresponding to genes and columns to samples.
#' @param grp Vector of phenotype groups of the samples, which represent valid variable names in R. Should be same
#' length as \code{ncol(object)}. If the vector is named, names should match \code{colnames(object)}.
#' @param name Name of PDF that gets written. Set to \code{NA} to suppress writing to file. "_dotplots" or "_boxplots"
#' is appended to plot filename, based on \code{type}.
#' @param main.v Character vector of main titles for plot.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param type Type of plot. Either \code{"dot"} or \code{"box"}.
#' @param color Passed to \code{\link[ggplot2]{scale_fill_manual}} and \code{\link[ggplot2]{scale_color_manual}}.
#' @param x.angle Passed to \code{theme(axis.text.x = element_text(angle))}.
#' @param add.se Add standard error of the mean.
#' @param dotsize Passed to \code{\link[ggplot2]{geom_dotplot}} \code{dotsize}.
#' @param bins Used to calculate binwidth, which is passed to \code{\link[ggplot2]{geom_dotplot}} \code{binwidth}.
#' @export

plot_by_grp <- function(object, grp, name='top_genes', main.v='', xlab = 'Group',  ylab='Log2 Expression', type='dot',
                        color = NULL, x.angle = 0, add.se = FALSE, dotsize = 0.7, bins = 30){
  if (!requireNamespace("ggplot2", quietly = TRUE)){
    stop("Package ggplot2 needed for this function to work. Please install it.", call. = FALSE)
  }

  if (is.vector(object)){ object <- t(as.matrix(object)) }
  if (all(main.v == '') & !is.null(rownames(object))){ main.v <- rownames(object) }
  grp <- factor(grp)

  stopifnot(ncol(object) == length(grp), nrow(object) == length(main.v), colnames(object) == names(grp), type %in% c('dot', 'box'))

  if (!is.na(name)){ grDevices::pdf(paste0(name, '_', type, 'plots.pdf')) }

  for (i in 1:nrow(object)){
    object2p <- data.frame(Exprs = object[i, ], Group = grp)
    binwidth <- (max(object2p$Exprs, na.rm = TRUE) - min(object2p$Exprs, na.rm = TRUE))/bins
    ggp <- ggplot2::ggplot(data = object2p, mapping = ggplot2::aes(x = Group, y = Exprs)) + ggplot2::theme_bw()
    ggp <- ggp + ggplot2::labs(title = main.v[i], x = xlab, y = ylab) + ggplot2::theme(legend.position = "none")
    if (type == 'dot'){
      ggp <- ggp + ggplot2::geom_dotplot(mapping = ggplot2::aes(fill = Group), binaxis='y', stackdir='center',
                                         binwidth = binwidth, dotsize = dotsize)
      if(add.se) { ggp <- ggp + ggplot2::stat_summary(fun.data = mean_se, geom = "crossbar", width = 0.3) }
    } else {
      ggp <- ggp + ggplot2::geom_boxplot(mapping = ggplot2::aes(fill = Group))
    }
    if(!is.null(color)) { ggp <- ggp + ggplot2::scale_fill_manual(values = color) + ggplot2::scale_color_manual(values = color) }
    if(x.angle != 0){ ggp <- ggp + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x.angle, hjust = 1)) }
    graphics::plot(ggp)
  }
  if (!is.na(name)){ grDevices::dev.off() }
  return(invisible(ggp))
}
