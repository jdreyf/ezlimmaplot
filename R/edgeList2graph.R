#' Transform a data frame edge list into an \code{igraph} graph
#'
#' Transforms an edge list into an \code{igraph} graph. The edge list can come from
#' columns 1 and 3 of a SIF (Simple Interaction Format) file, or can otherwise
#' contain pairs of nodes connected in the graph.
#'
#' @param edge.lst Matrix-like object with two columns that define undirected edge-list.
#' @details This function applies \code{\link[igraph]{graph_from_edgelist}} then \code{\link[igraph]{simplify}}.
#' @return An \code{igraph} graph.
#' @export

edgelist2graph <- function(edge.lst){
  stopifnot(length(dim(edge.lst))==2, ncol(edge.lst)==2)
  net <- igraph::simplify(igraph::graph_from_edgelist(as.matrix(edge.lst), directed = FALSE))
  return(net)
}
