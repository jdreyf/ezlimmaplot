#' Make dot plots for the significance and direction of pathway analysis
#'
#' Make dot plots like clusterProfiler showing the significance and proportion of genes in each pathway
#' with p-value below 5%.
#'
#' @param tab Table of output from \code{ezlimma::roast_*}.
#' @param prefix.v A vector of prefixes that prefix \code{.p}, \code{.FDR}, or \code{.logFC} in \code{colnames(tab)}.
#' These will be the comparisons shown in the plot. If \code{NULL} these are inferred from \code{colnames(tab)} that end with \code{.p}.
#' @param type.sig Either "p" or "FDR"; type of significance to show.
#' @param cut.sig Numeric in [0, 1]. Pathways need to have significance of type \code{type.sig < cut.sig} in a comparison to be shown on the dot plot.
#' @param ntop Integer number of top pathways to show.
#' @inheritParams ezheat
#' @inheritParams ezvenn
#' @inheritParams barplot_pwys
#' @return Invisibly, the ggplot object.
#' @importFrom rlang !!
#' @export

# 3 columns per contrast: up, down, mixed
# color is prop; size is -log10 of p or q
# https://github.com/YuLab-SMU/enrichplot/blob/0a01deaa5901b8af37d78c561e5f8200b4748b59/R/dotplot.R
dotplot_pwys <- function(tab, prefix.v=NULL, name = NA, type.sig=c("p", "FDR"), cut.sig=0.05, ntop = 20, pwys_nm_size = 100, width = 8, height = 8){
  type.sig <- match.arg(type.sig)
  if (is.null(prefix.v)){
    p.colnms <- ezlimma:::grep_cols(tab=tab, p.cols="p")
    p.colnms <- p.colnms[-grep("Mixed", p.colnms)]
    prefix.v <- gsub(pattern = "(\\.|^)p$", replacement = "", x=p.colnms)
    stopifnot(sapply(tab[, grep("Direction$", colnames(tab)), drop=FALSE], FUN=function(vv) all(vv %in% c("Up", "Down"))))
  } else {
    stopifnot(sapply(tab[, paste0(prefix.v, ".Direction"), drop=FALSE], FUN=function(vv) all(vv %in% c("Up", "Down") )))
  }

  if (!is.na(name)){
    name <- paste0(name, "_dotplots.pdf")
    grDevices::pdf(name, width, height)
    on.exit(grDevices::dev.off())
  }

  rownames(tab) <- pwy.nms <- substr(rownames(tab), 1, pwys_nm_size)

  dir.v <- c("Up", "Down", "Mixed")
  col.labs <- paste(rep(prefix.v, each=3), dir.v, sep=".")
  ds <- data.frame()

  for (prefix in prefix.v){
    for (dir.tmp in dir.v){
      col.lab.nm <- paste(prefix, dir.tmp, sep=".")
      if (dir.tmp != "Mixed"){
        prop.suffix <- paste0("Prop", dir.tmp, "P05")
        dir.col.tmp <- paste0(prefix, ".Direction")
        rows.tmp <- tab[, dir.col.tmp] == dir.tmp
        if (!any(rows.tmp)) next
        tab.tmp <- tab[rows.tmp,]
        prop.v <- tab.tmp[, paste(prefix, prop.suffix, sep=".")]
      } else {
        tab.tmp <- tab
        prop.suffix.v <- paste0("Prop", c("Up", "Down"), "P05")
        prop.v <- rowSums(tab.tmp[, paste(prefix, prop.suffix.v, sep=".")])
        prefix <- paste(prefix, "Mixed", sep=".")
      }
      ds <- dplyr::bind_rows(ds, dplyr::bind_cols(Comparison = col.lab.nm, Pwy = rownames(tab.tmp), `Prop P<5%` = prop.v,
                              p=tab.tmp[, paste0(prefix, ".p")], FDR=tab.tmp[, paste0(prefix, ".FDR")]))
    }
  }

  # filter out non-significant comparison x pathway rows
  # could do this in for loop separately for Mixed and not Mixed, but here I only need to do it once
  # then subset to 1st ntop rows: this does not cause error or warning even if ntop > nrow(ds)
  ds <- ds |> dplyr::filter(!!rlang::sym(type.sig) < cut.sig) |>
    dplyr::slice(1:ntop)
  if (nrow(ds) == 0) stop("No pathways had ", type.sig, " < ", cut.sig, ".")

  # i probably need to modify the font for long pwy nms
  # enrichplot::dotplot uses theme_dose(font.size) w/ default font size 12
  # angle axis labels: https://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2
  ggp <- ggplot2::ggplot(data = ds, mapping=ggplot2::aes(x=factor(Comparison, levels = col.labs, ordered = TRUE), y=factor(Pwy),
                                                         size = -log10(!!rlang::sym(type.sig)), color=`Prop P<5%`)) +
    ggplot2::geom_point() + ggplot2::xlab(NULL) + ggplot2::ylab(NULL) +
    ggplot2::labs(caption = "Prop P<5% = proportion of genes in pwy with P < 0.05 in given direction or for *Mixed* in either direction") +
    ggplot2::theme_bw() + ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 45))
  graphics::plot(ggp)
  invisible(ggp)
}
