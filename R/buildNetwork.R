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
  # Create graph.
  sif <- ppis %>% dplyr::select(osEntrezA, osEntrezB)
  nodes <- unique(c(sif$osEntrezA, sif$osEntrezB))
  # Get gene symbols, suppress output with sink.
  symbols <- AnnotationDbi::mapIds(org.Mm.eg.db,
    keys = as.character(nodes),
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )
  # Node attribute table.
  noa <- data.table::data.table(node = nodes, symbol = symbols)
  # Check.
  not_mapped <- nodes[is.na(symbols)]
  if (sum(is.na(symbols)) != 0) {
    stop("Unable to map all Entrez IDs to gene symbols!")
  }
  # Remove from sif and noa.
  sif <- filter(sif, sif$osEntrezA %in% not_mapped | sif$osEntrezB %in% not_mapped) 
  noa <- filter(noa, nodes %in% not_mapped)
# Build igraph object.
  g <- graph_from_data_frame(sif, directed = FALSE, vertices = noa)
  g <- simplify(g) # remove any redundant edges.
  nNodes <- length(V(g))
  nEdges <- length(E(g))
  message(paste(nEdges, "edges identified among", nNodes, "nodes!"))
  # Save the noa and sif files.
  if (save == TRUE) {
    data.table::fwrite(noa, "noa.csv")
    data.table::fwrite(sif, "sif.csv")
  }
  # Return igraph object.
  return(g)
}
