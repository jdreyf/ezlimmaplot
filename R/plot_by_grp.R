#' Make boxplot or dotplots for a feature per group
#'
#' Make boxplot or dotplots for a feature per group using \pkg{ggplot2}.
#'
#' @param grp Vector of phenotype groups of the samples, which represent valid variable names in R. Should be same
#' length as \code{ncol(object)}. If the vector is named, names should match \code{colnames(object)}.
#' @param main.v Character vector of main titles for plot.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param type Type of plot. Either \code{"dot"} or \code{"box"}.
#' @param manual.color Vector passed to both \code{\link[ggplot2:scale_manual]{scale_fill_manual}} and
#' \code{\link[ggplot2:scale_manual]{scale_color_manual}}.
#' @param x.angle Angle of x-axis text, passed to \code{theme(axis.text.x = element_text(angle))}.
#' @param add.se Logical indicating if to add standard error of the mean to dotplot.
#' @param dotsize Passed to \code{\link[ggplot2]{geom_dotplot}} \code{dotsize}.
#' @param bins Used to calculate binwidth, which is passed to \code{\link[ggplot2]{geom_dotplot}} \code{binwidth}.
#' @inheritParams ezheat
#' @details Based on \code{type}, \code{"_dotplots"} or \code{"_boxplots"} is appended to \code{name}.
#'
#' If \code{object} is a data frame, it is coerced to a \code{matrix} using \code{data.matrix}.
#' @return Invisibly, a \code{ggplot} object from the last row of \code{object} that was plotted.
#' @export

plot_by_grp <- function(object, grp, name="topgenes", main.v="", xlab = "Group",  ylab="Log2 Expression", type="dot",
                        manual.color = NULL, x.angle = 0, add.se = FALSE, dotsize = 1, bins = 30){
  if (is.vector(object)){ object <- t(as.matrix(object)) }
  if (is.data.frame(object)){ object <- data.matrix(object) }
  if (all(main.v == "") & !is.null(rownames(object))){ main.v <- rownames(object) }
  grp <- factor(grp)

  stopifnot(nrow(object) > 0, ncol(object) > 0, ncol(object) == length(grp), nrow(object) == length(main.v),
            colnames(object) == names(grp), type %in% c("dot", "box"))

  if (!is.na(name)){ grDevices::pdf(paste0(name, "_", type, "plots.pdf")) }

  for (i in 1:nrow(object)){
    object2p <- data.frame(Exprs = object[i, ], Group = grp)
    binwidth <- (max(object2p$Exprs, na.rm = TRUE) - min(object2p$Exprs, na.rm = TRUE))/bins
    ggp <- ggplot2::ggplot(data = object2p, mapping = ggplot2::aes(x = Group, y = Exprs)) + ggplot2::theme_bw()
    ggp <- ggp + ggplot2::labs(title = main.v[i], x = xlab, y = ylab) + ggplot2::theme(legend.position = "none")
    if (type == "dot"){
      ggp <- ggp + ggplot2::geom_dotplot(mapping = ggplot2::aes(fill = Group), binaxis="y", stackdir="center",
                                         binwidth = binwidth, dotsize = dotsize)
      if(add.se) { ggp <- ggp + ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "crossbar", width = 0.3) }
    } else {
      ggp <- ggp + ggplot2::geom_boxplot(mapping = ggplot2::aes(fill = Group))
    }
    if(!is.null(manual.color)) { ggp <- ggp + ggplot2::scale_fill_manual(values = manual.color) +
      ggplot2::scale_color_manual(values = manual.color) }
    if(x.angle != 0){ ggp <- ggp + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x.angle, hjust = 1)) }
    graphics::plot(ggp)
  }
  if (!is.na(name)){ grDevices::dev.off() }
  return(invisible(ggp))
}
