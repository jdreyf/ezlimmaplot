% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ezvolcano.R
\name{ezvolcano}
\alias{ezvolcano}
\title{Volcano plot}
\usage{
ezvolcano(
  tab,
  lfc.col = NA,
  sig.col = NA,
  lab.col = "Gene.Symbol",
  ntop.sig = 0,
  ntop.lfc = 0,
  comparison = NULL,
  alpha = 0.4,
  name = "volcano",
  ann.rnames = NULL,
  up.ann.color = "black",
  down.ann.color = "black",
  shape = 19,
  base.size = 14,
  x.bound = NULL,
  y.bound = NULL,
  type.sig = c("p", "FDR"),
  cut.lfc = 1,
  cut.sig = 0.05,
  up.cut.color = "black",
  down.cut.color = "black",
  lines.sig = NA,
  raster = FALSE,
  sep = ".",
  na.lab = c("---", ""),
  seed = 0,
  plot = TRUE
)
}
\arguments{
\item{tab}{Table of output from \code{ezlimma}.}

\item{lfc.col}{Column name or index of tab with logFC. Some features should be > 0 and others < 0.}

\item{sig.col}{Column name or index of tab with p-values or FDRs.}

\item{lab.col}{Column name or index of tab with labels, such as gene symbol, annotating features. If \code{NULL},
no points are annotated.}

\item{ntop.sig}{Number of top significant features to annotate.}

\item{ntop.lfc}{Number of top logFC features to annotate.}

\item{comparison}{Name of contrast to plot. If given, it's assumed that \code{lfc.col=paste0(comparison, '.logFC')}
and \code{sig.col=paste0(comparison, '.p') or paste0(comparison, '.FDR')}, and these are over-ridden.}

\item{alpha}{Transparency for non-annotated points, passed to \code{\link[ggplot2]{geom_point}}.}

\item{name}{Name of file to create. Set to \code{NA} to plot to screen instead of to file.}

\item{ann.rnames}{Character vector of additional rownames of \code{tab} to annotate; elements must be in \code{rownames(tab)}.}

\item{up.ann.color}{Color for annotated points that are upregulated (\code{logFC>0}).}

\item{down.ann.color}{Color for annotated points that are downregulated (\code{logFC<0}).}

\item{shape}{Shape of non-annotated dots. Set to 1 for empty dots.}

\item{base.size}{Numeric; font size text inherits from e.g. axis.title is 80 percent of this, and axis tick text is 64 percent.}

\item{x.bound}{x-axis limits are set to \code{c(-x.bound, x.bound)}. If \code{NULL, x.bound=max(abs(tab[,lfc.col]))}.}

\item{y.bound}{y-axis limits are set to \code{c(0, y.bound)}. If \code{NULL, y.bound=max(tab[,'nlg10sig'])}.}

\item{type.sig}{Type of significance y-axis should use, either "p" or "FDR".}

\item{cut.lfc}{Points need to have \code{|logFC| >= cut.lfc} to have \code{cut.color}.}

\item{cut.sig}{Points need to have significance \code{tab[,sig.col] <= cut.sig} to have \code{cut.color}.}

\item{up.cut.color}{Color for up-regulated points that meet both \code{cut.lfc} and \code{cut.sig} but are not labeled.}

\item{down.cut.color}{Color for down-regulated points that meet both \code{cut.lfc} and \code{cut.sig} but are not labeled.}

\item{lines.sig}{Numeric vector of values of \code{sig.type} at which to draw lines. For example, if
\code{type.sig="p"}, you may want to set \code{lines.sig = 0.05}, which will draw a line at \code{y = -log10(0.05)}.}

\item{raster}{Rasterize points using \code{ggrastr} so plot is lighter.}

\item{sep}{Separator string between contrast names and suffix such as \code{logFC}.}

\item{na.lab}{Character vector of labels in \code{lab.col} to treat as missing, in addition to \code{NA}.}

\item{seed}{Numeric seed for reproducibility since \code{ggrepel} uses a random algorithm.}

\item{plot}{Logical; should plot be generated?}
}
\value{
Invisibly, a \code{ggplot} object.
}
\description{
Volcano plot in ggplot2 using output from \code{ezlimma} package, e.g. \code{limma_contrasts}.
The logFC and significance columns are supplied by the user, or inferred from the supplied \code{comparison},
e.g. if \code{comparison="AvsB"}, \code{ezvolcano} expects column names \code{AvsB.logFC} and,
depending on \code{type.sig}, \code{AvsB.p} or \code{AvsB.FDR}.
}
\details{
If \code{ntop.sig>0} or \code{ntop.lfc>0}, then \code{lab.col} must be in \code{colnames(tab)}.
For annotated points, \code{up.ann.color} & \code{down.ann.color} dominate \code{cut.color}.
}
\seealso{
multi_volcano
}
