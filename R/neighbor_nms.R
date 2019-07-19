#' Get names of nodes' neighbors
#'
#' Get node names of nodes' direct neighbors
#'
#' @param nodes Node names in the graph
#' @inheritParams plot_pwy

neighbor_nms <- function(gr, nodes){
  return(igraph::V(gr)$name[unlist(igraph::ego(gr, order=1, match(nodes, igraph::V(gr)$name)))])
}
