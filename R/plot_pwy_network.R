#' Plot network diagram for a pathway
#'
#' Plot nodes most impacting a pathway as a network diagram with node color corresponding to z-score.
#' The most impactful nodes are inferred by assuming the input here is the same as
#' was used to calculate pathway significance. The \code{feat.tab} is assumed to have all nodes to consider
#' for plotting. Including nodes with \code{NA} stats allows these to be included in networks as a connector.
#' An \code{annot.colnm} should have preferred labels for all nodes.
#'
#' @param feat.tab Matrix-like object for all analytes with row names as analyte names, which should
#' have some overlap with \code{G.pwy$genes}.
#' @param G.pwy Element of list object \code{G} used in \pkg{ezlimma}. Should have elements
#' \code{name}, \code{description}, \code{genes}. Accessible via \code{G[[pathway.name]]}.
#' @param gr graph object of class \code{igraph}, possibly a result of \code{\link{edgelist2graph}}.
#' @param stat.colnm Column name in \code{feat.tab} indicating scores to select top analytes for this pathway and plot.
#' Stats should be finite unless they are \code{NA}. Stats for all nodes are used to define color scale,
#' for consistency across pathway visualizations.
#' @param annot.colnm Column name in \code{feat.tab} of node names to possibly replace \code{rownames(feat.tab)}.
#' Set this column to be \code{rownames(feat.tab)} if these are desired node names. Must not be \code{NA}.
#' @param ntop Number of top impactful analytes to plot, >=2. Their network neighbors may also be included.
#' @param name Name of file to plot to. If \code{NULL}, a filename is created using \code{colnames(Gmat.pwy)[1]}.
#' Set to \code{NA} to plot to screen instead of to file.
#' @param colorbar.nm Title of color bar.
#' @param repel Logcal; use the repel functionality from \pkg{ggrepel} to avoid overlapping text?
#' @param plot Logical; should plot be generated?
#' @param seed Numeric seed value for reproducibility.
#' @inheritParams ezlimma::roast_contrasts
#' @details Unmeasured nodes have stat of \code{NA} and are drawn gray.
#'
#' @return Invisibly, a list with \code{\link[tidygraph]{tbl_graph}}, a subclass of
#' \pkg{igraph} so every \pkg{igraph} method will work as expected. If \code{plot=TRUE}, there
#' is a second element of the list, which is a \code{ggplot} object.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @importFrom rlang :=
#' @export

# alternative?
# do not calculate redundantly to p(h)
# check if data 2-sided instead of alternative, which may be diff than used in p(h), so confusing
# cannot infer graph from ker
plot_pwy_network <- function(feat.tab, G.pwy, gr, stat.colnm, annot.colnm, ntop = 7, name = NULL, colorbar.nm=stat.colnm,
                     alternative=c("two.sided", "greater", "less"), repel=FALSE, plot = TRUE, seed = 0){

  stopifnot(ncol(feat.tab)>=1, limma::isNumeric(feat.tab[, stat.colnm]), !is.null(colnames(feat.tab)),
            !is.null(rownames(feat.tab)), all(is.na(feat.tab[, stat.colnm]) | is.finite(feat.tab[, stat.colnm])),
            length(intersect(rownames(feat.tab), G.pwy$genes)) > 0, names(G.pwy)==c("name", "description", "genes"),
            length(intersect(rownames(feat.tab), igraph::V(gr)$name)) > 0,
            stat.colnm %in% colnames(feat.tab), annot.colnm %in% colnames(feat.tab), !is.na(feat.tab[, annot.colnm]),
            igraph::is_simple(gr), is.numeric(ntop), ntop>=2, is.logical(repel), is.logical(plot), is.numeric(seed))

  pwy.nm <- G.pwy$name
  n.nona <- sum(!is.na(feat.tab[G.pwy$genes, stat.colnm]))
  if (n.nona == 0) stop("All pathway nodes are NA", call. = FALSE) else ntop <- min(ntop, n.nona)

  # pwy nodes not in gr
  missing.v <- setdiff(G.pwy$genes, igraph::V(gr)$name)
  gr <- igraph::add_vertices(graph=gr, nv=length(missing.v), name=missing.v)
  gr.pwy <- igraph::induced_subgraph(gr, vid=which(igraph::V(gr)$name %in% G.pwy$genes))

  top.nodes <- order_stats_by_alternative(stat.v = stats::setNames(feat.tab[G.pwy$genes, stat.colnm], nm=G.pwy$genes),
                                          alternative=alternative)

  gr.pwy.top <- igraph::induced_subgraph(gr.pwy, names(top.nodes)[1:ntop])

  n.comp <- igraph::components(gr.pwy.top)$no
  if (n.comp > 1){
    # order nodes w/ lowest eigen centrality to avoid h20, atp, ...
    ec <- igraph::eigen_centrality(gr)$vector
    # NAs are ordered last
    candidate.nodes <- igraph::V(gr.pwy)$name[order(-feat.tab[igraph::V(gr.pwy)$name, stat.colnm],
                                                    ec[igraph::V(gr.pwy)$name])]
    n.add <- min(floor(sqrt(n.comp)), length(candidate.nodes))

    # want 2nd added node to reduce n.comp more, etc.
    while (n.comp > 1 && n.add > 0){
      # what na nodes reduce n.comp most?
      new.ncomp <- apply(as.matrix(candidate.nodes), MARGIN=1, FUN=function(xx){
        igraph::components(igraph::induced_subgraph(gr.pwy, vid=union(igraph::V(gr.pwy.top)$name, xx)))$no
      })
      names(new.ncomp) <- candidate.nodes
      # only want those that reduce n.comp
      new.ncomp <- new.ncomp[new.ncomp < n.comp]
      if (length(new.ncomp) >= 1){
        # accounts for order of analytes, which accounts for non-NA or NA
        # central.v <- names(new.ncomp)[order(new.ncomp)][1:min(n.add, length(new.ncomp))]
        new.top.v <- names(new.ncomp)[which.min(new.ncomp)]
        gr.pwy.top <- igraph::induced_subgraph(gr.pwy, vid=union(igraph::V(gr.pwy.top)$name, new.top.v))
      }
      n.comp <- igraph::components(gr.pwy.top)$no
      n.add <- n.add - 1
    } # end while
  }

  gg.pwy <- tidygraph::as_tbl_graph(gr.pwy.top) %>%
    dplyr::mutate(!!stat.colnm := feat.tab[igraph::V(gr.pwy.top)$name, stat.colnm]) %>%
    dplyr::mutate(!!annot.colnm := feat.tab[igraph::V(gr.pwy.top)$name, annot.colnm])

  if (plot){
    if (is.null(name)) name <- paste0(ezlimma::clean_filenames(pwy.nm), "_ntop", ntop)
    if (!is.na(name)) {
      grDevices::pdf(paste0(name, ".pdf"))
      on.exit(expr = grDevices::dev.off())
    }

    set.seed(seed)
    # need to set font fam to avoid "font fam not found" errors
    # geom_label() draws a rectangle behind the text
    # repel repels labels from node centers?
    ggg <- ggraph::ggraph(gg.pwy, layout = "nicely") + ggraph::geom_edge_link() +
      ggraph::theme_graph(base_family = "sans") +
      ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size=4))) +
      ggraph::geom_node_point(mapping=ggplot2::aes(color = !!rlang::ensym(stat.colnm)), size=12) +
      ggraph::geom_node_text(mapping=ggplot2::aes(label = !!rlang::ensym(annot.colnm)), repel = repel, size=6)

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
    return(invisible(list(gg.pwy=gg.pwy, ggplot=ggg)))
  } else {
    return(invisible(list(gg.pwy=gg.pwy)))
  }
}
