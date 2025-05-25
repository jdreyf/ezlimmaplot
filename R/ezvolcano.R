#' Volcano plot
#'
#' Volcano plot in ggplot2 using output from \code{ezlimma} package, e.g. \code{limma_contrasts}.
#' The logFC and significance columns are supplied by the user, or inferred from the supplied \code{comparison},
#' e.g. if \code{comparison="AvsB"}, \code{ezvolcano} expects column names \code{AvsB.logFC} and,
#' depending on \code{type.sig}, \code{AvsB.p} or \code{AvsB.FDR}.
#'
#' @param lfc.col Column name or index of tab with logFC. Some features should be > 0 and others < 0.
#' @param sig.col Column name or index of tab with p-values or FDRs.
#' @param lab.col Column name or index of tab with labels, such as gene symbol, annotating features. If \code{NULL},
#' no points are annotated.
#' @param ntop.sig Number of top significant features to annotate.
#' @param ntop.lfc Number of top logFC features to annotate.
#' @param comparison Name of contrast to plot. If given, it's assumed that \code{lfc.col=paste0(comparison, '.logFC')}
#' and \code{sig.col=paste0(comparison, '.p') or paste0(comparison, '.FDR')}, and these are over-ridden.
#' @param alpha Transparency for non-annotated points, passed to \code{\link[ggplot2]{geom_point}}.
#' @param ann.rnames Character vector of additional rownames of \code{tab} to annotate; elements must be in \code{rownames(tab)}.
#' @param up.ann.color Color for annotated points that are upregulated (\code{logFC>0}).
#' @param down.ann.color Color for annotated points that are downregulated (\code{logFC<0}).
#' @param shape Shape of non-annotated dots. Set to 1 for empty dots.
#' @param x.bound x-axis limits are set to \code{c(-x.bound, x.bound)}. If \code{NULL, x.bound=max(abs(tab[,lfc.col]))}.
#' @param y.bound y-axis limits are set to \code{c(0, y.bound)}. If \code{NULL, y.bound=max(tab[,'nlg10sig'])}.
#' @param type.sig Type of significance y-axis should use, either "p" or "FDR".
#' @param cut.color Color of points that meet both \code{cut.lfc} and \code{cut.sig} but are not labeled.
#' @param cut.lfc Points need to have \code{|logFC| >= cut.lfc} to have \code{cut.color}.
#' @param cut.sig Points need to have significance \code{tab[,sig.col] <= cut.sig} to have \code{cut.color}.
#' @param lines.sig Numeric vector of values of \code{sig.type} at which to draw lines. For example, if
#' \code{type.sig="p"}, you may want to set \code{lines.sig = 0.05}, which will draw a line at \code{y = -log10(0.05)}.
#' @param axis.text.size Numeric; font size of axis text.
#' @param text.repel.size Numeric; font size of labels, which are repelled from each other.
#' @param raster Rasterize points using \code{ggrastr} so plot is lighter.
#' @param sep Separator string between contrast names and suffix such as \code{logFC}.
#' @param na.lab Character vector of labels in \code{lab.col} to treat as missing, in addition to \code{NA}.
#' @param seed Numeric seed for reproducibility since \code{ggrepel} uses a random algorithm.
#' @inheritParams ezheat
#' @inheritParams ezvenn
#' @details If \code{ntop.sig>0} or \code{ntop.lfc>0}, then \code{lab.col} must be in \code{colnames(tab)}.
#' For annotated points, \code{up.ann.color} & \code{down.ann.color} dominate \code{cut.color}.
#' @return Invisibly, a \code{ggplot} object.
#' @seealso multi_volcano
#' @export

ezvolcano <- function(tab, lfc.col=NA, sig.col=NA, lab.col='Gene.Symbol', ntop.sig=0, ntop.lfc=0, comparison=NULL, alpha=0.4,
                      name='volcano', ann.rnames=NULL, up.ann.color='black', down.ann.color='black', shape = 19,
                      x.bound=NULL, y.bound=NULL, type.sig=c('p', 'FDR'), cut.color="black", cut.lfc=1, cut.sig=0.05, lines.sig=NA,
                      axis.text.size = 12, text.repel.size=3,
                      raster = FALSE, sep='.', na.lab=c('---', ''), seed = 0, plot=TRUE){
  set.seed(seed = seed) # ggrepel is random
  type.sig <- match.arg(type.sig)
  # can't annot if no lab.col
  # if null, need to halt and return FALSE
  have.labs <- !(is.null(lab.col) || !(lab.col %in% colnames(tab)))
  if (!have.labs & (ntop.lfc > 0 | ntop.sig > 0 | !is.null(ann.rnames))){
    stop("ntop.lfc > 0 | ntop.sig > 0 | !is.null(ann.rnames) but lab.col is null or not in the column names of tab.")
  }

  if (type.sig=="p"){
    y.lab <- expression("-"*log[10]~p*"-"*value)
  } else {
    y.lab <- expression("-"*log[10]~FDR)
  }

  # infer columns
  if (!is.null(comparison)){
    lfc.col <- paste0(comparison, sep, 'logFC')
    if (type.sig=="p"){
      sig.col <- paste0(comparison, sep, 'p')
    } else {
      sig.col <- paste0(comparison, sep, 'FDR')
    }
    if (!is.na(name)) name <- paste(comparison, name, sep='_')
  }

  stopifnot(ntop.sig==as.integer(ntop.sig), lfc.col %in% colnames(tab), sig.col %in% colnames(tab),
            length(x.bound)<=1, length(y.bound)<=1, all(is.na(lines.sig)) || (is.numeric(lines.sig) && length(lines.sig)<=5), is.logical(plot))

  # There is no need for prefixing !!! or !! or {{ with rlang::. These operators are not function calls, they are specially interpreted by rlang in data-masked arguments.
  tab <- data.frame(tab, nlg10sig = -log10(tab[,sig.col]), check.names = FALSE) |>
    dplyr::mutate(rnm = rownames(tab), drctn.color = ifelse(!!rlang::sym(lfc.col) > 0, up.ann.color, down.ann.color), annot=FALSE,
                  na.lab = dplyr::case_when(!have.labs ~ TRUE,
                                            is.na(!!rlang::sym(lab.col)) | !!rlang::sym(lab.col) %in% na.lab ~ TRUE,
                                            .default = FALSE))

  rnms.top.lfc <- rnms.top.sig <- NULL
  if (ntop.lfc > 0) rnms.top.lfc <- tab |> dplyr::arrange(-abs(!!rlang::sym(lfc.col))) |> dplyr::slice(1:ntop.lfc) |> dplyr::pull(rnm)
  if (ntop.sig > 0) rnms.top.sig <- tab |> dplyr::arrange(!!rlang::sym(sig.col)) |> dplyr::slice(1:ntop.sig) |> dplyr::pull(rnm)

  # ?geom_text_repel has example to hide many labels but repel from all points where the labels they don't want are ""
  # default shape=19; point size = 1.5
  tab <- tab |>
    dplyr::mutate(map.grp = dplyr::case_when(rnm %in% c(ann.rnames, rnms.top.lfc, rnms.top.sig) ~ "annot",
                                             abs(!!rlang::sym(lfc.col)) > cut.lfc & !!rlang::sym(sig.col) <= cut.sig ~ "cut",
                                             .default = "rest"),
                  color.point = dplyr::case_match(map.grp, "annot" ~ drctn.color, "cut" ~ cut.color, .default = "black"),
                  label.point = dplyr::case_when(map.grp == "annot" & !na.lab ~ !!rlang::sym(lab.col), .default = ""),
                  alpha.point = dplyr::case_when(label.point == "" ~ alpha, .default = 1),
                  size.point = dplyr::case_when(map.grp == "annot" ~ 2, .default = 1.5),
                  shape.point = dplyr::case_when(map.grp == "annot" ~ 19, .default = shape))

  # want symmetric x-axis
  if (is.null(x.bound)) x.bound <- max(abs(tab[,lfc.col]))
  if (is.null(y.bound)) y.bound <- max(tab[,'nlg10sig'])

  # construct ggplot object
  vol <- ggplot2::ggplot(data=tab, mapping=ggplot2::aes(x=!!rlang::sym(lfc.col), y = nlg10sig, color = color.point, alpha=alpha.point, shape=shape.point)) +
    ggplot2::theme_bw() + ggplot2::theme(axis.text=ggplot2::element_text(size=axis.text.size, face="bold")) + ggplot2::xlab(expression(log[2]~fold~change)) + ggplot2::ylab(y.lab) +
    ggplot2::xlim(c(-x.bound, x.bound)) + ggplot2::ylim(c(0, y.bound))

  if (all(!is.na(lines.sig))){
    # y is already -log10(sig)
    # linetype values that we want = 2-6
    # 0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash
    lty.v <- c("dashed", "dotted", "dotdash", "longdash", "twodash")
    lty.leg <- lty.v[order(as.character(lines.sig))]
    names(lty.leg) <- sort(as.character(lines.sig))
    vol <- vol + ggplot2::geom_hline(data=data.frame(lines.sig),
                                     mapping = ggplot2::aes(yintercept = -log10(lines.sig), linetype=as.character(lines.sig)),
                                     show.legend = TRUE) + ggplot2::scale_linetype_manual(values=lty.leg) +
      ggplot2::guides(linetype=ggplot2::guide_legend(title=type.sig))
  }

  if (!is.null(comparison)){
    # rm "In|IN|in|OF|Of|of" since "vs" should be enough
    first.grp <- unlist(strsplit(x=gsub("_", "", comparison), split = "vs|VS|Vs"))[1]
    # label left & right sides
    vol <- vol + ggplot2::ggtitle(comparison) +
      ggplot2::geom_text(mapping=ggplot2::aes(x=2*x.bound/3, y = -Inf, label = paste0("Up in ", first.grp)), color="darkgrey",
                         vjust = -0.5, show.legend=FALSE) +
      ggplot2::geom_text(mapping=ggplot2::aes(x=-2*x.bound/3, y = -Inf, label = paste0("Down in ", first.grp)), color="darkgrey",
                         vjust = -0.5, show.legend=FALSE)
  }

  # finalize plot
  if (raster){
    vol <- vol + ggrastr::rasterize(ggplot2::geom_point(data=tab |> dplyr::filter(map.grp == "rest"), size=2)) +
      ggplot2::geom_point(data = tab |> dplyr::filter(map.grp != "rest"), mapping = ggplot2::aes(size = size.point), show.legend = FALSE)
  } else {
    vol <- vol + ggplot2::geom_point(mapping = ggplot2::aes(size = size.point), show.legend = FALSE)
  }

  vol <- vol + ggplot2::scale_color_identity() + ggplot2::scale_alpha_identity() + ggplot2::scale_shape_identity() + ggplot2::scale_size_identity() +
    ggrepel::geom_text_repel(mapping=ggplot2::aes(label=label.point), show.legend = FALSE, size=text.repel.size)

  if (plot){
    if (!is.na(name)) ggplot2::ggsave(filename=paste0(name, ".png"), plot=vol) else graphics::plot(vol)
  }
  return(invisible(vol))
}
