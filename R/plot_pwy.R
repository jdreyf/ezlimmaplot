#' Plot network diagram for a pathway
#'
#' Plot nodes most impacting a pathway as a network diagram with node color corresponding to z-score.
#' The most impactful nodes are inferred by assuming the input here is the same as
#' was used to calculate pathway significance. The \code{feat.tab} is assumed to have all nodes to consider
#' for plotting. Including nodes with \code{NA} stats allows these to be included in networks as a connector.
#' An \code{annot.col} should have preferred labels for all nodes.
#'
#' @param feat.tab Matrix-like object for all analytes with rownames as analyte names, which should
#' have some overlap with \code{G.pwy$genes}.
#' @param G.pwy Element of list object \code{G} used in \pkg{ezlimma} or \pkg{{Pame}}. Should have elements
#' \code{name}, \code{description}, \code{genes}. Accessible via \code{G[[pathway.name]]}.
#' @param gr graph object of class \code{igraph}.
#' @param stat.colnm Column name in \code{feat.tab} indicating scores to select top analytes for this pathway and plot.
#' Stats should be finite unless they are \code{NA}. Stats for all nodes are used to define color scale,
#' for consistency across pathway visualizations.
#' @param annot.col Column index or name in \code{feat.tab} of annotations to replace \code{rownames(G.pwy)}.
#' \code{NA}s in \code{annot.col} are ignored, and these nodes retain their rowname as a label.
#' @param ntop Number of top impactful analytes to plot, >=2. Their network neighbors may also be included.
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
plot_pwy <- function(feat.tab, G.pwy, gr, stat.colnm, annot.col, ntop = 7, name = NULL, colorbar.nm=stat.colnm,
                     alternative=c("two.sided", "greater", "less"), repel=FALSE, plot = TRUE, seed = 0){

  stopifnot(ncol(feat.tab)>=1, limma::isNumeric(feat.tab[, stat.colnm]), !is.null(colnames(feat.tab)),
            !is.null(rownames(feat.tab)), all(is.na(feat.tab[, stat.colnm]) | is.finite(feat.tab[, stat.colnm])),
            length(intersect(rownames(feat.tab), G.pwy$genes)) > 0, names(G.pwy)==c("name", "description", "genes"),
            stat.colnm %in% colnames(feat.tab), annot.col %in% colnames(feat.tab), any(!is.na(feat.tab[, annot.col])),
            igraph::is_simple(gr), is.numeric(ntop), ntop>=2, is.logical(repel), is.logical(plot), is.numeric(seed))

  pwy.nm <- G.pwy$name
  pwy.genes <- intersect(rownames(feat.tab)[!is.na(feat.tab[, stat.colnm])], G.pwy$genes)
  feat.pwy <- feat.tab[pwy.genes,]
  if (ntop > nrow(feat.pwy)) ntop <- nrow(feat.pwy)
  if (is.null(name)) name <- paste0(ezlimma::clean_filenames(pwy.nm), "_ntop", ntop)

  # want highest impact
  top.nodes <- switch(alternative,
   greater = rownames(feat.pwy)[order(feat.pwy[, stat.colnm], decreasing = TRUE)][1:ntop],
   two.sided = rownames(feat.pwy)[order(abs(feat.pwy[, stat.colnm]), decreasing = TRUE)][1:ntop],
   less = rownames(feat.pwy)[order(feat.pwy[, stat.colnm], decreasing = FALSE)][1:ntop])

  # pwy nodes not in gr
  missing.v <- setdiff(pwy.genes, igraph::V(gr)$name)
  gr <- igraph::add_vertices(graph=gr, nv=length(missing.v), name=missing.v)

  # add edges
  gr.pwy <- igraph::induced_subgraph(gr, vid=which(igraph::V(gr)$name %in% top.nodes))
  feat.pwy <- feat.tab[igraph::V(gr.pwy)$name,]

  top.neighbors <- neighbor_nms(gr, top.nodes)
  pwy.gene.neighbors <- intersect(pwy.genes, top.neighbors)
  # inc genes in pwy & gr but are NA; could help connect pwy
  pwy.na.genes <- intersect(rownames(feat.tab)[is.na(feat.tab[, stat.colnm])], G.pwy$genes)
  pwy.na.gene.neighbors <- intersect(pwy.na.genes, top.neighbors)

  n.comp <- igraph::components(gr.pwy)$no
  if (n.comp > 1){
    # order nodes w/ lowest eigen centrality
    ec <- igraph::eigen_centrality(gr)$vector
    pgn.o <- pwy.gene.neighbors[order(ec[pwy.gene.neighbors])]
    pngn.o <- pwy.na.gene.neighbors[order(ec[pwy.na.gene.neighbors])]
    candidate.nodes <- c(pgn.o, pngn.o)

    n.add <- min(floor(sqrt(n.comp)), length(candidate.nodes))
    # what na nodes reduce n.comp most?
    new.ncomp <- apply(as.matrix(candidate.nodes), MARGIN=1, FUN=function(xx){
      igraph::components(igraph::induced_subgraph(gr, vid=union(top.nodes, xx)))$no
    })
    names(new.ncomp) <- candidate.nodes
    # only want those that reduce n.comp
    new.ncomp <- new.ncomp[new.ncomp < n.comp]
    if (length(new.ncomp) >= 1){
      # accounts for order of analytes, which accounts for non-NA or NA
      central.v <- names(new.ncomp)[order(new.ncomp)][1:min(n.add, length(new.ncomp))]
      gr.pwy <- igraph::induced_subgraph(gr, vid=which(igraph::V(gr)$name %in% c(top.nodes, central.v)))
      feat.pwy <- feat.tab[igraph::V(gr.pwy)$name,]
    }
  }

  gg.pwy <- tidygraph::as_tbl_graph(gr.pwy) %>%
    dplyr::mutate(!!stat.colnm := feat.pwy[, stat.colnm])
  igraph::V(gg.pwy)$name[match(rownames(feat.pwy), igraph::V(gg.pwy)$name)] <- feat.pwy[, annot.col]

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
        ggraph::geom_node_point(mapping=ggplot2::aes(color = !!rlang::ensym(stat.colnm)), size=12) +
        ggraph::geom_node_text(mapping=ggplot2::aes(label=I(name)), repel = repel)

      # colorbar title
      guide <- ggplot2::guide_colourbar(title=colorbar.nm, title.theme=ggplot2::element_text(face="bold", size=16))

      # creates common lim using all feat.tab[, stat.colnm], so pwys have consistent colorbar
      if (min(feat.tab[, stat.colnm], na.rm = TRUE)<0 && max(feat.tab[, stat.colnm], na.rm = TRUE)>0){
        # use yellow in middle to distinguish NAs, which are grey
        ggg <- ggg + ggplot2::scale_colour_distiller(type="div", palette = "RdYlBu", direction = -1, guide=guide,
                      limits=c(-max(abs(feat.tab[, stat.colnm]), na.rm = TRUE),
                               max(abs(feat.tab[, stat.colnm]), na.rm = TRUE)))
      } else {
        ggg <- ggg + ggplot2::scale_colour_distiller(type="seq", palette = "Reds", direction = 1, guide=guide,
                      limits=range(feat.tab[, stat.colnm], na.rm = TRUE))
      }
      graphics::plot(ggg)
    })
    if (!is.na(name)) grDevices::dev.off()
  }
  return(invisible(gg.pwy))
}
