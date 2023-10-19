#' Make dot plots for the significance and DE proportion of pathway analysis
#'
#' Make dot plots like clusterProfiler showing the significance and proportion of genes in each pathway
#' with p-value below 0.05.
#'
#' @param tab Table of output from \code{ezlimma::roast_*}.
#' @param prefix.v A vector of prefixes that prefix \code{.p}, \code{.FDR}, or \code{.logFC} in \code{colnames(tab)}.
#' These will be the comparisons shown in the plot. If \code{NULL} these are inferred from \code{colnames(tab)} that end with \code{.p}.
#' @param type.sig Either "p" or "FDR"; type of significance to show.
#' @param cut.sig Numeric in [0, 1]. Pathways need to have significance of type \code{type.sig < cut.sig} in a comparison to be shown on the dot plot.
#' @param ntop Integer number of top pathways to show.
#' @param name Name of file to create. Set to \code{NA} to plot to screen instead of to file; otherwise, "_dotplot.pdf" is appended to the name.
#' @param mixed Character string. Should mixed statistics be included, should they be the only statistics included, or should they be excluded?
#' @param caption Logical; should the caption explaining the color bar title be included?
#' @param colorbar.title Character title of color bar.
#' @inheritParams ezheat
#' @inheritParams ezvenn
#' @inheritParams barplot_pwys
#' @return Invisibly, the ggplot object.
#' @export

# 3 columns per contrast: up, down, mixed
# color is prop; size is -log10 of p or q
# https://github.com/YuLab-SMU/enrichplot/blob/0a01deaa5901b8af37d78c561e5f8200b4748b59/R/dotplot.R
# prefix.v=NULL; name = NA; type.sig="p"; cut.sig=0.05; ntop = 20; pwys_nm_size = 100; width = 8; height = 8; mixed="include"; caption=TRUE
dotplot_pwys <- function(tab, prefix.v=NULL, name = NA, type.sig=c("p", "FDR"), cut.sig=0.05, ntop = 20, pwys_nm_size = 100, width = 8, height = 8,
                         mixed = c("include", "exclude", "only"), caption=TRUE, colorbar.title = "Prop P<5%"){
  type.sig <- match.arg(type.sig)
  mixed <- match.arg(mixed)
  stopifnot(mixed == "exclude" | any(grepl("Mixed", colnames(tab))), cut.sig > 0, cut.sig <= 1)
  tab <- as.data.frame(tab)
  rownames(tab) <- substr(rownames(tab), 1, pwys_nm_size)
  # extract prefix
  if (is.null(prefix.v)){
    p.colnms <- ezlimma:::grep_cols(tab=tab, p.cols="p")
    p.colnms <- p.colnms[-grep("Mixed", p.colnms)]
    prefix.v <- gsub(pattern = "(\\.|^)p$", replacement = "", x=p.colnms)
  }
  stopifnot(sapply(tab[, paste0(prefix.v, ".Direction"), drop=FALSE], FUN=function(vv) all(vv %in% c("Up", "Down") )))

  if (!is.na(name)){
    name <- paste0(name, "_dotplot.pdf")
    grDevices::pdf(name, width=width, height=height)
    on.exit(grDevices::dev.off())
  }

  dir.v <- switch(mixed, include = c("Up", "Down", "Mixed"), exclude = c("Up", "Down"), only = "Mixed")
  col.labs <- paste(rep(prefix.v, each=length(dir.v)), dir.v, sep=".")

  ds <- pivot_roast_longer(tab=tab, prefix.v = prefix.v, direction.v = dir.v)

  # filter out non-significant comparison x pathway rows
  # could do this in for loop separately for Mixed and not Mixed, but here I only need to do it once
  # then subset to 1st ntop rows: this does not cause error or warning even if ntop > nrow(ds)
  ds <- ds |> dplyr::arrange(!!rlang::sym(type.sig)) |>
    dplyr::filter(!!rlang::sym(type.sig) < cut.sig) |>
    dplyr::filter(Pwy %in% (ds |> dplyr::slice(1:ntop) |> dplyr::pull(Pwy)))
  if (nrow(ds) == 0) stop("No pathways had ", type.sig, " < ", cut.sig, ".")

  # i probably need to modify the font and/or use str_wrap() for long pwy nms
  # enrichplot::dotplot uses theme_dose(font.size) w/ default font size 12
  # angle axis labels: https://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2
  # order legends: https://stackoverflow.com/questions/11393123/controlling-ggplot2-legend-display-order, but `color=guide_legend(order=1)` treats it as factor :-/
  ggp <- ggplot2::ggplot(data = ds, mapping=ggplot2::aes(x=factor(Comparison, levels = col.labs, ordered = TRUE),
                                                         y=factor(Pwy, levels = rev(unique(Pwy)), ordered = TRUE),
                                                         size = -log10(!!rlang::sym(type.sig)), color=Prop_p05)) +
    ggplot2::geom_point() + ggplot2::xlab(NULL) + ggplot2::ylab(NULL) + ggplot2::labs(color = colorbar.title) +
    ggplot2::theme_bw() + ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 45)) +
    ggplot2::guides(size = ggplot2::guide_legend(order = 2))
  if (caption){
    capt <- paste(colorbar.title, "= proportion of genes in pwy with P < 0.05 in given direction or for *Mixed* in either direction")
    ggp <- ggp + ggplot2::labs(caption = capt)
  }
  graphics::plot(ggp)
  return(invisible(ggp))
}
