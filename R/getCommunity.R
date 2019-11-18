#' getCommunity
#'
#' Get community of proteins that interact with seed proteins.
#'
#' @param graph - igraph object
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none
#'
#' @export
#'
#' @examples
#' getCommunity(graph, seeds)
getCommunity <- function(graph, seeds, k = 2) {
  suppressPackageStartUpMessages({
    library(igraph)
    library(dplyr)
  })
  # Get neighborhood of every seed.
  ego_graphs <- make_ego_graph(graph, nodes = seeds)
  # Extract community nodes and create new graph.
  subg <- induced_subgraph(graph, unlist(sapply(ego_graphs, function(x) names(V(x)))))
  all_nodes <- names(V(subg))
  # Remove nodes with less than 2 connections to seed nodes.
  adjm <- as.matrix(as_adjacency_matrix(subg))
  dm <- adjm[rownames(adjm) %in% seeds, ]
  ds <- apply(dm, 2, sum)
  keep <- names(ds)[ds >= k]
  community <- c(seeds, all_nodes[all_nodes %in% keep])
  community_graph <- induced_subgraph(graph, community)
  return(community_graph)
}
