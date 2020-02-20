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
#' @param cut.color Color of points that meet both \code{cut.lfc} and \code{cut.sig}. If \code{NULL}, cutoffs are ignored.
#' @param cut.lfc Points need to have \code{|logFC| >= cut.lfc} to have \code{cut.color}.
#' @param cut.sig Points need to have significance \code{tab[,sig.col] <= cut.sig} to have \code{cut.color}.
#' @param lines.sig Numeric vector of values of \code{sig.type} at which to draw lines. For example, if
#' \code{type.sig="p"}, you may want to set \code{lines.sig = 0.05}, which will draw a line at \code{y = -log10(0.05)}.
#' @param sep Separator string between contrast names and suffix such as \code{logFC}.
#' @param na.lab Character vector of labels in \code{lab.col} to treat as missing, in addition to \code{NA}.
#' @inheritParams ezheat
#' @inheritParams ezvenn
#' @details If \code{ntop.sig>0} or \code{ntop.lfc>0}, then \code{lab.col} must be in \code{colnames(tab)}.
#' For annotated points, \code{up.ann.color} & \code{down.ann.color} dominate \code{cut.color}.
#' @return Invisibly, a \code{ggplot} object.
#' @seealso multi_volcano
#' @export

ezvolcano <- function(tab, lfc.col=NA, sig.col=NA, lab.col='Gene.Symbol', ntop.sig=0, ntop.lfc=0, comparison=NULL, alpha=0.4,
                      name='volcano', ann.rnames=NULL, up.ann.color='black', down.ann.color='black', shape = 16,
                      x.bound=NULL, y.bound=NULL, type.sig=c('p', 'FDR'), cut.color=NULL, cut.lfc=1, cut.sig=0.05, lines.sig=NA,
                      sep='.', na.lab=c('---', ''), plot=TRUE){
  # can't annot if no lab.col
  if (is.null(lab.col)){
    if (ntop.lfc > 0){
      ntop.lfc <- 0
      warning("ntop.lfc > 0 but lab.col is null.")
    }
    if (ntop.sig > 0){
      ntop.sig <- 0
      warning("ntop.sig > 0 but lab.col is null.")
    }
    if (!is.null(ann.rnames)){
      ann.rnames <- NULL
      warning("ann.rnames is not null, but lab.col is null.")
    }
  }# end if lab.col null

  type.sig <- match.arg(type.sig)
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

  stopifnot((ntop.sig==0 & ntop.lfc==0) | lab.col %in% colnames(tab), ntop.sig==as.integer(ntop.sig),
            ntop.lfc==as.integer(ntop.lfc), is.null(ann.rnames)|ann.rnames %in% rownames(tab),
            lfc.col %in% colnames(tab), sig.col %in% colnames(tab), any(tab[,lfc.col]<0), any(tab[,lfc.col]>=0),
            length(x.bound)<=1, length(y.bound)<=1, all(is.na(lines.sig)) || (is.numeric(lines.sig) && length(lines.sig)<=5),
            is.logical(plot))

  tab <- data.frame(tab, nlg10sig = -log10(tab[,sig.col]), check.names = FALSE)
  # want symmetric x-axis
  if (is.null(x.bound)) x.bound <- max(abs(tab[,lfc.col]))
  if (is.null(y.bound)) y.bound <- max(tab[,'nlg10sig'])

  ## need to order tab st annotated points are plotted last, so they are seen even if in dense areas
  ## however, not concerned about cut points, since they are in outer region
  ## https://stackoverflow.com/questions/15706281/controlling-order-of-points-in-ggplot2-in-r

  # ann.rnames to plot with symbol
  ind.annot <- NULL
  if (!is.null(ann.rnames)){
    ind.ann.rnames <- which(rownames(tab) %in% ann.rnames)
    tab <- tab[c(setdiff(1:nrow(tab), ind.ann.rnames), ind.ann.rnames),, drop=FALSE]
    ind.annot <- (nrow(tab) - length(ind.ann.rnames) + 1):nrow(tab)
  }

  # ntop indices to plot with symbol
  if (ntop.sig > 0 | ntop.lfc > 0){
    na.lab.ind <- which(is.na(tab[,lab.col]) | tab[,lab.col] %in% na.lab)
    if (ntop.lfc > 0) top.lfc.ind <- order(-abs(tab[,lfc.col]))[1:ntop.lfc] else top.lfc.ind <- NULL
    if (ntop.sig > 0) top.sig.ind <- order(tab[,sig.col])[1:ntop.sig] else top.sig.ind <- NULL
    ind.ntop <- setdiff(union(top.sig.ind, top.lfc.ind), na.lab.ind)
    ind.annot <- union(ind.annot, ind.ntop)
  }

  # construct ggplot object
  vol <- ggplot2::ggplot(data=tab, mapping=ggplot2::aes_string(x=lfc.col, y='nlg10sig')) + ggplot2::theme_bw() +
    ggplot2::theme(axis.text=ggplot2::element_text(size=12, face="bold")) + ggplot2::xlab(expression(log[2]~fold~change)) +
    ggplot2::xlim(c(-x.bound, x.bound)) + ggplot2::ylim(c(0, y.bound)) + ggplot2::ylab(y.lab)

  if (!is.null(comparison)){
    first.grp <- unlist(strsplit(x=gsub("_", "", comparison), split = "vs|VS|Vs|In|IN|in|OF|Of|of"))[1]

    vol <- vol + ggplot2::ggtitle(comparison) +
      ggplot2::geom_text(mapping=ggplot2::aes(x=2*x.bound/3, y = -Inf, label = paste0("Up in ", first.grp)), color="darkgrey",
                         vjust = -0.5, show.legend=FALSE) +
      ggplot2::geom_text(mapping=ggplot2::aes(x=-2*x.bound/3, y = -Inf, label = paste0("Down in ", first.grp)), color="darkgrey",
                         vjust = -0.5, show.legend=FALSE)
  }

  # plot annotated with color
  if (!is.null(ind.annot)){
    ind.annot.up <- ind.annot[which(tab[ind.annot, lfc.col] >= 0)]
    if (length(ind.annot.up) > 0){
      vol <- vol + ggplot2::geom_point(data=tab[ind.annot.up,], size=2, color = up.ann.color) +
        ggrepel::geom_text_repel(data=tab[ind.annot.up,], mapping=ggplot2::aes_string(x=lfc.col, y='nlg10sig', label=lab.col),
                           size=3, vjust=2, color = up.ann.color)
    }

    ind.annot.down <- ind.annot[which(tab[ind.annot, lfc.col] < 0)]
    if (length(ind.annot.down) > 0){
      vol <- vol + ggplot2::geom_point(data=tab[ind.annot.down,], size=2, color = down.ann.color) +
        ggrepel::geom_text_repel(data=tab[ind.annot.down,], mapping=ggplot2::aes_string(x=lfc.col, y='nlg10sig', label=lab.col),
                           size=3, vjust=2, color = down.ann.color)
    }
  }

  # cut points
  if (!is.null(cut.color)){
    ind.cut <- which(abs(tab[,lfc.col]) > cut.lfc & tab[,sig.col] <= cut.sig)
    ind.cut <- setdiff(ind.cut, ind.annot)
    if (length(ind.cut) > 0){
      vol <- vol + ggplot2::geom_point(data=tab[ind.cut,], alpha=alpha, size=2, color = cut.color, shape=shape)
    }
  } else {
    ind.cut <- NULL
  }

  # plot rest
  ind.rest <- setdiff(1:nrow(tab), union(ind.annot, ind.cut))
  vol <- vol + ggplot2::geom_point(data=tab[ind.rest,], alpha=alpha, size=2, shape=shape)

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

  if (plot){
    if (!is.na(name)) ggplot2::ggsave(filename=paste0(name, ".png"), plot=vol) else graphics::plot(vol)
  }
  return(invisible(vol))
}
