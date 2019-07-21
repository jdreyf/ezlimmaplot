#' Plot network diagram for a pathway
#'
#' Plot nodes most impacting a pathway as a network diagram with node color corresponding to z-score.
#' The most impactful nodes are inferred by assuming the input here is the same as
#' was used to calculate pathway significance.
#'
#' @param feat.tab Matrix-like object for all analytes with rownames as analyte names, which should
#' have some overlap with \code{G.pwy$genes}. 1st column with stats.
#' Stats should be finite unless they are \code{NA}. Stats for all nodes are used to define color scale,
#' for consistency across pathway visualizations. 1st column name used as stat descriptor in plot.
#' @param G.pwy Element of list object \code{G} used in \pkg{ezlimma} or \pkg{{Pame}}. Should have elements
#' \code{name}, \code{description}, \code{genes}. Accessible via \code{G[[pathway.name]]}.
#' @param gr graph object of class \code{igraph}.
#' @param stat.colnm Column name in \code{feat.tab} indicating scores to select top analytes for this pathway and plot.
#' @param annot.col Column index or name in \code{feat.tab} of annotations to replace \code{rownames(G.pwy)}.
#' \code{NA}s in \code{annot.col} are ignored, and these nodes retain their rowname as a label.
#' @param ntop Number of top impactful analytes to plot. Their network neighbors may also be included.
#' @param name Name of file to plot to. If \code{NULL}, a filename is created using \code{colnames(Gmat.pwy)[1]}.
#' Set to \code{NA} to plot to screen instead of to file.
#' @param colorbar.nm Title of color bar.
#' @param repel Logcal; use the repel functionality from \pkg{ggrepel} to avoid overlapping text?
#' @param plot Logical; should plot be generated?
#' @inheritParams ezlimma::roast_contrasts
#' @details Unmeasured nodes have stat of \code{NA} and are drawn gray.
#'
#' @return Invisibly, a \code{\link[tidygraph]{tbl_graph}}, a subclass of
#' \pkg{igraph} so every \pkg{igraph} method will work as expected.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @importFrom rlang :=
#' @export

# alternative?
# do not calculate redundantly to p(h)
# check if data 2-sided instead of alternative, which may be diff than used in p(h), so confusing
# cannot infer graph from ker
plot_pwy <- function(feat.tab, G.pwy, gr, stat.colnm, annot.col=NULL, ntop = 7, name = NULL, colorbar.nm=stat.colnm,
                     alternative=c("two.sided", "greater", "less"), repel=FALSE, plot = TRUE, seed = 0){
  stopifnot(ncol(feat.tab)>=1, limma::isNumeric(feat.tab[, stat.colnm]), !is.null(colnames(feat.tab)),
            !is.null(rownames(feat.tab)), all(is.na(feat.tab[, stat.colnm]) | is.finite(feat.tab[, stat.colnm])),
            length(intersect(rownames(feat.tab), G.pwy$genes)) > 0, names(G.pwy)==c("name", "description", "genes"),
            stat.colnm %in% colnames(feat.tab), annot.col %in% colnames(feat.tab),
            igraph::is_simple(gr), is.numeric(ntop), is.logical(repel), is.logical(plot), is.numeric(seed))

  pwy.nm <- G.pwy$name
  if (ntop > nrow(feat.tab)) ntop <- nrow(feat.tab)
  if (is.null(name)) name <- paste0(ezlimma::clean_filenames(pwy.nm), "_ntop", ntop)

  # expand graph to include all features, even if they are isolated
  new.v <- setdiff(G.pwy$genes, igraph::V(gr)$name)
  gr <- igraph::add_vertices(graph=gr, nv=length(new.v), name=new.v)

  pwy.nodes <- G.pwy$genes
  # look for 1st degree neighbors to plot
  pwy.neighbors <- setdiff(neighbor_nms(gr, pwy.nodes), pwy.nodes)

  # want highest positive impact
  top.nodes <- switch(alternative,
   greater = rownames(feat.tab)[order(feat.tab[, stat.colnm], decreasing = TRUE)][1:ntop],
   two.sided = rownames(feat.tab)[order(abs(feat.tab[, stat.colnm]), decreasing = TRUE)][1:ntop],
   less = rownames(feat.tab)[order(feat.tab[, stat.colnm], decreasing = FALSE)][1:ntop],
  )

  # get neighbors
  pwy.neighbors.ss <- intersect(pwy.neighbors, top.nodes)
  # get pwy nodes with large stats OR pwy nodes connected to neighbors with large stats
  pwy.nodes.ss <- union(intersect(pwy.nodes, top.nodes),
                        intersect(pwy.nodes, neighbor_nms(gr, pwy.neighbors.ss)))

  gr.pwy <- igraph::induced_subgraph(gr, vid=which(igraph::V(gr)$name %in% c(pwy.nodes.ss, pwy.neighbors.ss)))

  gg.pwy <- tidygraph::as_tbl_graph(gr.pwy) %>%
    dplyr::mutate(Pathway = c("outside", "inside")[(igraph::V(gr.pwy)$name %in% pwy.nodes)+1])  %>%
    dplyr::mutate(!!stat.colnm := feat.tab[igraph::V(gr.pwy)$name, stat.colnm])

  # sub V(gg.pwy)$name w/ feat.tab[,2]
  if (!is.null(annot.col) && any(!is.na(feat.tab[, annot.col]))){
    na.annot <- is.na(feat.tab[, annot.col])
    feat.tab[na.annot, annot.col] <- names(feat.tab[na.annot, annot.col])
    nms.int <- intersect(names(feat.tab[, annot.col]), igraph::V(gg.pwy)$name)
    if (length(nms.int) > 0){
      igraph::V(gg.pwy)$name[match(nms.int, igraph::V(gg.pwy)$name)] <- feat.tab[nms.int, annot.col]
    }
  }

  if (plot){
    if (!is.na(name)) grDevices::pdf(paste0(name, ".pdf"))
    try({
      set.seed(seed)
      # need to set font fam to avoid "font fam not found" errors
      # geom_label() draws a rectangle behind the text
      # repel repels labels from node centers?
      ggg <- ggraph::ggraph(gg.pwy, layout = "nicely") + ggraph::geom_edge_link() +
        ggraph::theme_graph(base_family = "sans") +
        ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size=4))) +
        ggraph::geom_node_point(mapping=ggplot2::aes(shape=Pathway, color = !!rlang::ensym(stat.colnm)), size=12) +
        ggraph::geom_node_text(mapping=ggplot2::aes(label=I(name)), repel = repel)

      # creates common lim using all feat.tab[, stat.colnm], so pwys have consistent colorbar
      if (min(feat.tab[, stat.colnm])<0 && max(feat.tab[, stat.colnm])>0){
        # use yellow in middle to distinguish NAs, which are grey
        ggg <- ggg + ggplot2::scale_colour_distiller(type="div", palette = "RdYlBu", direction = -1,
                      limits=c(-max(abs(feat.tab[, stat.colnm])), max(abs(feat.tab[, stat.colnm]))))
      } else {
        ggg <- ggg + ggplot2::scale_colour_distiller(type="seq", palette = "Reds", direction = 1,
                                                   limits=range(feat.tab[, stat.colnm]))
      }
      graphics::plot(ggg)
    })
    if (!is.na(name)) grDevices::dev.off()
  }
  return(invisible(gg.pwy))
}
