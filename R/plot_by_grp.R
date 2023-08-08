#' Make box, dot, violin or bar plots for a feature per group
#'
#' Make box, dot, violin or bar plotsfor a feature per group using \pkg{ggplot2}.
#'
#' @param grp Vector of phenotype groups of the samples, which represent valid variable names in R. Should be same
#' length as \code{ncol(object)}. If the vector is named, names should match \code{colnames(object)}.
#' @param main.v Character vector of main titles for plot.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param type Type of plot. Either \code{"dot"}, \code{"box"}, \code{"violin"}, or \code{"bar"}.
#' @param manual.color Vector passed to both \code{\link[ggplot2:scale_manual]{scale_fill_manual}} and
#' \code{\link[ggplot2:scale_manual]{scale_color_manual}}.
#' @param x.angle Angle of x-axis text, passed to \code{theme(axis.text.x=element_text(angle))}.
#' @param add.se Logical indicating if to add standard error of the mean to dot plot.
#' @param dotsize Passed to \code{\link[ggplot2]{geom_dotplot}} \code{dotsize} for dot plot.
#' @param bins Used to calculate binwidth, which is passed to \code{\link[ggplot2]{geom_dotplot}} \code{binwidth} for dot plot.
#' @param violin.alpha Colour transparency in [0,1] for violin in violin plot.
#' @param dot.alpha olour transparency in [0,1] for dots in violin plot.
#' @param add.dot Logical indicating if to add individual data points as dots to bar plot.
#' @param bar.width Width of the bars in bar plot.
#' @param errorbar.width Width of the errorbars in bar plot.
#' @inheritParams ezheat
#' @details Based on \code{type}, \code{"_dotplots"}, \code{"_boxplots"}, \code{"_violinplots"}, or \code{"_barplots"}is appended to \code{name}.
#'
#' If \code{object} is a data frame, it is coerced to a \code{matrix} using \code{data.matrix}.
#' @return Invisibly, a \code{ggplot} object from the last row of \code{object} that was plotted.
#' @export

plot_by_grp <- function(object, grp, name="topgenes", width=7, height=7, main.v="", xlab="Group",  ylab="Log2 Expression", type="dot",
                        manual.color=NULL, x.angle=0, add.se=FALSE, dotsize=0.7, bins=30, violin.alpha=0.3, dot.alpha=0.3, add.dot=FALSE,
                        bar.width=0.7, errorbar.width=0.3) {
  if (is.vector(object)) { object <- t(as.matrix(object)) }
  if (is.data.frame(object)) { object <- data.matrix(object) }
  if (all(main.v=="") & !is.null(rownames(object))) { main.v <- rownames(object) }
  grp <- factor(grp)

  stopifnot(nrow(object) > 0, ncol(object) > 0, ncol(object)==length(grp), nrow(object)==length(main.v),
            colnames(object)==names(grp), type %in% c("dot", "box", "violin", "bar"))

  if (!is.na(name)) {
    grDevices::pdf(paste0(name, "_", type, "plots.pdf"), width=width, height=height)
    on.exit(grDevices::dev.off())
  }

  for (i in 1:nrow(object)) {
    object2p <- data.frame(Exprs=object[i, ], Group=grp)
    object2p <- object2p[stats::complete.cases(object2p), ]

    if (type != "bar") {
      ggp <- ggplot2::ggplot(data=object2p, mapping=ggplot2::aes(x=Group, y=Exprs))

      if (type=="dot") {
        binwidth <- (max(object2p$Exprs) - min(object2p$Exprs))/bins
        ggp <- ggp + ggplot2::geom_dotplot(mapping=ggplot2::aes(fill=Group), binaxis="y", stackdir="center",
                                           binwidth=binwidth, dotsize=dotsize)
        if (add.se) { ggp <- ggp + ggplot2::stat_summary(fun.data=ggplot2::mean_se, geom="crossbar", width=0.3) }
      } else if (type=="box") {
        ggp <- ggp + ggplot2::geom_boxplot(mapping=ggplot2::aes(fill=Group))
      } else {
        ggp <- ggp + ggplot2::geom_violin(mapping=ggplot2::aes(fill=Group), trim=FALSE, alpha=violin.alpha)
        ggp <- ggp + ggplot2::geom_jitter(mapping=ggplot2::aes(fill=Group), shape=21, position=ggplot2::position_jitter(0.2), alpha=dot.alpha)
      }
    } else {
      # object2p.avg <- plyr::ddply(object2p, "Group", plyr::summarise, N=length(Exprs), Mean=mean(Exprs), SE=stats::sd(Exprs)/sqrt(N))
      object2p.avg <- object2p |> dplyr::group_by(Group) |>
        dplyr::summarize(N = dplyr::n(), Mean = mean(Exprs), SE=stats::sd(Exprs)/sqrt(N))
      ggp <- ggplot2::ggplot(data=object2p.avg, mapping=ggplot2::aes(x=Group, y=Mean))
      ggp <- ggp + ggplot2::geom_col(mapping=ggplot2::aes(fill=Group), position=ggplot2::position_dodge(), width=bar.width, color="black")
      ggp <- ggp + ggplot2::geom_errorbar(mapping=ggplot2::aes(ymin=Mean, ymax=Mean+SE), position=ggplot2::position_dodge(), width=errorbar.width)
      if (add.dot) {ggp <- ggp + ggplot2::geom_jitter(mapping=ggplot2::aes(x=Group, y=Exprs), data=object2p,
                                                      shape=21, fill="black", size=3, position=ggplot2::position_jitter(0.2)) }
    }
    ggp <- ggp + ggplot2::theme_bw() + ggplot2::labs(title=main.v[i], x=xlab, y=ylab) + ggplot2::theme(legend.position="none")

    if(!is.null(manual.color)) { ggp <- ggp + ggplot2::scale_fill_manual(values=manual.color) +
      ggplot2::scale_color_manual(values=manual.color) }
    if(x.angle !=0) { ggp <- ggp + ggplot2::theme(axis.text.x=ggplot2::element_text(angle=x.angle, hjust=1)) }
    graphics::plot(ggp)
  }
  return(invisible(ggp))
}
