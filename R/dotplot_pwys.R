#' Make dot plots for the significance and direction of pathway analysis
#'
#' Make dot plots like clusterProfiler showing the significance and proportion of genes in each pathway
#' with p-value below 5%.
#'
#' @inheritParams ezheat
#' @inheritParams ezvenn
#' @inheritParams barplot_pwys
#' @return Invisibly, the ggplot object
#' @export


# 3 columns per contrast: up, down, mixed
# color is prop; size is -log10 of p or q
# https://github.com/YuLab-SMU/enrichplot/blob/0a01deaa5901b8af37d78c561e5f8200b4748b59/R/dotplot.R
dotplot_pwys <- function(tab, prefix.v=NULL, name = NA, width = 8, height = 8, ntop = 20, pwys_nm_size = 100){
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

  npwys <- min(ntop, nrow(tab))
  # add prop for mixed
  tab <- tab |> dplyr::slice(1:npwys)
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
      ds <- dplyr::bind_rows(ds, dplyr::bind_cols(Comparison = col.lab.nm, Pwy = rownames(tab.tmp), PropP05 = prop.v,
                              p=tab.tmp[, paste0(prefix, ".p")], FDR=tab.tmp[, paste0(prefix, ".FDR")]))
    }
  }

  # i probably need to modify the font for long pwy nms
  # enrichplot::dotplot uses theme_dose(font.size) w/ default font size 12
  # angle axis labels: https://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2
  ggp <- ggplot2::ggplot(data = ds, mapping=ggplot2::aes(x=factor(Comparison, levels = col.labs, ordered = TRUE), y=factor(Pwy), size = -log10(p), color=PropP05)) +
    ggplot2::geom_point() + ggplot2::xlab(NULL) + ggplot2::ylab(NULL) +
    ggplot2::labs(caption = "PropP05 = proportion of genes in pwy in given direction with P < 0.05 or for *Mixed* in either direction") +
    ggplot2::theme_bw() + ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 45))
  graphics::plot(ggp)
  invisible(ggp)
}
