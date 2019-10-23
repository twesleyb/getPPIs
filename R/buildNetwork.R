#' buildNetwork
#'
#' Given a HitPredict data.frame (PPIs) and a vector of entrez IDs, build a
#' protein-protein interaction graph.
#'
#' @param none
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
#' buildNetwork()
buildNetwork <- function(hitpredict, mygenes, taxid = 10090, save = TRUE) {
  # Fixme: Only works for mouse!!!
  # Imports.
  suppressPackageStartupMessages({
    require(getPPIs)
    require(dplyr)
    require(org.Mm.eg.db)
    require(igraph)
  })
  # Load proteins of interest.
  if (is.character(mygenes)) {
    # if character, then read file into R.
    myprots <- data.table::fread(mygenes)$Entrez
    idy <- grep("entrez", tolower(colnames(myprots)))
    mygenes <- unlist(myprots[, ..idy])
    names(mygenes) <- NULL
  } else {
    # Use genes passed by user.
    myprots <- mygenes
  }

  # Get interactions among genes of interest.
  ppis <- hitpredict %>% filter(osEntrezA %in% mygenes & osEntrezB %in% mygenes)
  # keep relevant columns.
  sif <- ppis %>% dplyr::select(
    osEntrezA, osEntrezB, Source_database,
    Interaction_detection_methods, Methods,
    EntrezA, EntrezB,
    Interactor_A_Taxonomy, Interactor_B_Taxonomy,
    Publications
  )
  # Node attributes.
  entrez <- unique(c(sif$osEntrezA, sif$osEntrezB, mygenes))
  # Get gene symbols, suppress output with sink.
  symbols <- AnnotationDbi::mapIds(org.Mm.eg.db,
    keys = as.character(entrez),
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )
  # Check that all nodes (entrez) are mapped to gene symbols.
  not_mapped <- entrez[is.na(symbols)]
  if (sum(is.na(symbols)) != 0) {
    message(paste("Unable to map", length(not_mapped), "Entrez IDs to gene symbols!"))
  }
  # Remove nodes that were not successfully mapped.
  nodes <- symbols
  names(nodes) <- entrez
  nodes <- na.omit(nodes)
  # Node attribute table.
  noa <- data.table::data.table(entrez = names(nodes), symbol = nodes)
  rownames(noa) <- noa$entrez
  # Insure that any interactions among unmapped genes are removed.
  out <- sif$osEntrezA %in% not_mapped | sif$osEntrezB %in% not_mapped
  sif <- subset(sif, !out)
  # Build igraph object.
  g <- graph_from_data_frame(sif, directed = FALSE, vertices = noa)
  g <- simplify(g) # remove any redundant edges.
  # Change node names to symbols.
  g <- set.vertex.attribute(g, "name", value = vertex_attr(g, "symbol"))
  # Status report.
  nNodes <- format(length(V(g)), 1, nsmall = 1, big.mark = ",")
  nEdges <- format(length(E(g)), 1, nsmall = 1, big.mark = ",")
  message(paste(nEdges, "edges identified among", nNodes, "nodes!"))
  # Save the noa and sif files.
  if (save == TRUE) {
    data.table::fwrite(noa, "noa.csv")
    data.table::fwrite(sif, "sif.csv")
  }
  # Return igraph object and noa an sif, which can be imported into Cytoscape..
  return(list("network" = g, "noa" = noa, "sif" = sif))
}
