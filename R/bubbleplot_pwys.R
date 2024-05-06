#' Make bubble plots for the significance, DE proportion, and counts of pathway analysis
#'
#' Make bubble plots like clusterProfiler showing the significance, proportion of genes with p-value below 0.05, and counts per comparison
#' per pathway. These are called "bubble plots" instead of "dot plots"  because they represent the counts as an additional dimension of
#' numerical (i.e. non-factor) data.
#'
#' @param name Name of file to create. Set to \code{NA} to plot to screen instead of to file; otherwise, "_bubbleplots.pdf" is appended to the name.
#' @param ntop Integer; number of top pathways to show.
#' @param caption Logical; should the caption explaining the x-axis be included?
#' @inheritParams dotplot_pwys
#' @inheritParams ezheat
#' @inheritParams ezvenn
#' @inheritParams barplot_pwys
#' @return Invisibly, a list of ggplot objects.
#' @export

# color = p; size = counts; x = DE percent
# one plot per comparison x direction
bubbleplot_pwys <- function(tab, prefix.v=NULL, name = NA, type.sig=c("p", "FDR"), cut.sig=0.05, ntop = 20, pwys_nm_size = 100, width = 8, height = 8,
                            caption=TRUE, mixed = c("include", "exclude", "only")){
  type.sig <- match.arg(type.sig)
  mixed <- match.arg(mixed)

  tab <- as.data.frame(tab)
  rownames(tab) <- gsub("_", " ", substr(rownames(tab), 1, pwys_nm_size))
  pwy.max.nchar <- max(nchar(rownames(tab)))
  stopifnot(mixed == "exclude" | any(grepl("Mixed", colnames(tab))), cut.sig > 0, cut.sig <= 1)
  # extract prefix
  if (is.null(prefix.v)){
    p.colnms <- ezlimma:::grep_cols(tab=tab, p.cols="p")
    p.colnms <- p.colnms[!grepl("Mixed", p.colnms)]
    prefix.v <- gsub(pattern = "(\\.|^)p$", replacement = "", x=p.colnms)
  }
  stopifnot(sapply(tab[, paste0(prefix.v, ".Direction"), drop=FALSE], FUN=function(vv) all(vv %in% c("Up", "Down") )))

  if (!is.na(name)){
    name <- paste0(name, "_bubbleplots.pdf")
    grDevices::pdf(name, width=width, height=height)
    on.exit(grDevices::dev.off())
  }

  dir.v <- switch(mixed, include = c("Up", "Down", "Mixed"), exclude = c("Up", "Down"), only = "Mixed")

  ds <- pivot_roast_longer(tab=tab, prefix.v = prefix.v, direction.v = dir.v) |>
    dplyr::filter(!!rlang::sym(type.sig) < cut.sig)
  if (nrow(ds) == 0) stop(paste("All", type.sig, ">=", cut.sig))
  cmprs.v <- ds |> dplyr::pull(Comparison) |> unique()

  ggp.lst <- list()
  for (cmpr in cmprs.v){
    # select top pwys by p but order pwys in plot by DE ratio
    ds.tmp <- ds |> dplyr::filter(Comparison == cmpr) |>
      dplyr::arrange(p) |>
      dplyr::slice(1:ntop) |>
      dplyr::arrange(-Prop_p05)
    pwy.max.nchar.tmp <- ds.tmp |> dplyr::pull(Pwy) |> nchar() |> max()

    ggp <- ggplot2::ggplot(data = ds.tmp, mapping=ggplot2::aes(x=100*Prop_p05, y=factor(Pwy, levels = rev(unique(Pwy)), ordered = TRUE),
                                                           size = Count, color = -log10(!!rlang::sym(type.sig)))) +
      ggplot2::geom_point() + ggplot2::xlab("DE percent") + ggplot2::ylab(NULL) + ggplot2::ggtitle(cmpr) +
      ggplot2::scale_size(range = c(9, 14)) + ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = 3)) +
      ggplot2::scale_y_discrete(labels = scales::label_wrap(width = ceiling(pwy.max.nchar.tmp/2) + 10)) +
      ggplot2::guides(color = ggplot2::guide_colorbar(order = 1)) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(color = "black", size = 12))
    if (caption){
      capt <- paste("DE percent = percent of genes in pwy with P < 0.05 in given direction or for *Mixed* in either direction")
      ggp <- ggp + ggplot2::labs(caption = capt)
    }
    graphics::plot(ggp)
    ggp.lst[[cmpr]] <- ggp
  }

  return(invisible(ggp.lst))
}
